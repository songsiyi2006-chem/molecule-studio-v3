"use strict";
let activeJob = null, lastJobId = null, jobGeneration = 0, pollTimer = null, pollFailures = 0, historyGeneration = 0;
let selectedConformer = 0, currentConformers = null, toastTimer = null;
const comparison = [];
const viewer = new MoleculeViewer($("molecule-canvas"), $("element-legend"));
const pageTitles = { analyze: "分析工作台", compare: "分子对比", database: "参考数据库", evidence: "模型证据" };
function toast(message) { clearTimeout(toastTimer); $("toast").textContent = message; $("toast").hidden = false; toastTimer = setTimeout(() => { $("toast").hidden = true; }, 3600); }
function selectPage(page) {
  if (!pageTitles[page]) return;
  document.querySelectorAll("[data-page]").forEach((button) => { const selected = button.dataset.page === page; button.classList.toggle("active", selected); if (selected) button.setAttribute("aria-current", "page"); else button.removeAttribute("aria-current"); });
  Object.keys(pageTitles).forEach((key) => { $(`page-${key}`).hidden = key !== page; }); $("page-name").textContent = pageTitles[page];
  if (page === "database" && !databaseLoaded) searchDatabase(0);
  if (page === "compare") renderComparison();
}
document.querySelectorAll("[data-page]").forEach((button) => button.addEventListener("click", () => selectPage(button.dataset.page)));
function bindSwitch(selector, callback) {
  const buttons = [...document.querySelectorAll(selector)];
  buttons.forEach((button, index) => { button.addEventListener("click", () => callback(button)); button.addEventListener("keydown", (event) => { let next; if (event.key === "ArrowRight") next = (index + 1) % buttons.length; if (event.key === "ArrowLeft") next = (index - 1 + buttons.length) % buttons.length; if (event.key === "Home") next = 0; if (event.key === "End") next = buttons.length - 1; if (next !== undefined) { event.preventDefault(); callback(buttons[next]); buttons[next].focus(); } }); });
}
function selectAnalysisTab(tab) { document.querySelectorAll("[data-analysis-tab]").forEach((button) => { const selected = button.dataset.analysisTab === tab; button.setAttribute("aria-selected", String(selected)); button.tabIndex = selected ? 0 : -1; $(`${button.dataset.analysisTab}-workspace`).hidden = !selected; }); }
bindSwitch("[data-analysis-tab]", (button) => selectAnalysisTab(button.dataset.analysisTab));
function setTheme(theme) { document.documentElement.dataset.theme = theme; $("theme-toggle").setAttribute("aria-label", theme === "dark" ? "切换浅色模式" : "切换深色模式"); try { localStorage.setItem("molecule-studio-theme", theme); } catch {} document.dispatchEvent(new Event("themechange")); }
try { if (localStorage.getItem("molecule-studio-theme") === "dark") setTheme("dark"); } catch {}
$("theme-toggle").addEventListener("click", () => setTheme(document.documentElement.dataset.theme === "dark" ? "light" : "dark"));
function safeSvg(text) {
  if (typeof text !== "string") return null;
  const parsed = new DOMParser().parseFromString(text, "image/svg+xml"), root = parsed.documentElement;
  if (root.localName !== "svg" || parsed.querySelector("parsererror")) return null;
  const permitted = new Set(["svg", "g", "path", "rect", "circle", "ellipse", "line", "polyline", "polygon", "text", "tspan", "defs", "clipPath", "title", "desc"]);
  for (const element of [root, ...root.querySelectorAll("*")]) {
    if (!permitted.has(element.localName)) { element.remove(); continue; }
    for (const attribute of [...element.attributes]) { const name = attribute.name.toLowerCase(); if (name.startsWith("on") || name.includes("href") || name === "id" || /url\s*\(|javascript:|@import|expression\s*\(/i.test(attribute.value)) element.removeAttribute(attribute.name); }
  }
  root.removeAttribute("width"); root.removeAttribute("height"); root.setAttribute("aria-hidden", "true"); return document.importNode(root, true);
}
function showStructure(destination, svg) { destination.replaceChildren(); const sanitized = safeSvg(svg); destination.append(sanitized || node("p", "small-empty", "结构图暂不可用")); }
function updateInputCount() { $("input-count").textContent = `${$("smiles").value.length.toLocaleString("zh-CN")} / 2,000`; }
function setSingleView(view) { $("empty-state").hidden = view !== "empty"; $("job-panel").hidden = view !== "job"; $("result-content").hidden = view !== "result"; }
function setBusy(busy) {
  singleBusy = busy;
  $("analysis-form").querySelectorAll("input,textarea,button").forEach((input) => { input.disabled = busy; });
  $("analyze-button").textContent = busy ? "分析进行中…" : "开始分析 ↗";
  $("cancel-job").disabled = busy && !activeJob;
}
function invalidateResult() { lastResult = null; lastJobId = null; currentConformers = null; $("view-3d").disabled = true; $("export-sdf").disabled = true; hideMessage("analysis-error"); $("smiles").removeAttribute("aria-invalid"); setSingleView("empty"); }
function setInput(smiles) { $("smiles").value = smiles; updateInputCount(); invalidateResult(); }
$("smiles").addEventListener("input", () => { updateInputCount(); invalidateResult(); });
$("smiles").addEventListener("keydown", (event) => { if (event.key === "Enter" && (event.ctrlKey || event.metaKey)) { event.preventDefault(); $("analysis-form").requestSubmit(); } });
document.querySelectorAll('#analysis-form input').forEach((input) => input.addEventListener("change", invalidateResult));
document.querySelectorAll("[data-example]").forEach((button) => button.addEventListener("click", () => { const example = examples[button.dataset.example]; setInput(example.smiles); $("smiles").focus(); announce(`已填入${example.name}的结构。`); }));
$("clear-button").addEventListener("click", () => { if (!singleBusy) { setInput(""); $("smiles").focus(); } });
const stageLabels = { queued: "任务已排队", validating: "正在校验分子结构", validation: "正在校验分子结构", prediction: "正在分析分子性质", predicting: "正在分析分子性质", descriptors: "正在计算分子描述符", conformers: "正在探索三维构象", embedding: "正在生成三维构象", optimizing: "正在优化构象", conformer_generation: "正在探索三维构象", completed: "分析完成", failed: "分析未完成", cancelled: "计算已取消" };
function showJobProgress(job) {
  $("job-stage").textContent = stageLabels[job.stage] || job.stage || "正在计算";
  $("job-time").textContent = finite(job.elapsed_seconds) ? `已运行 ${numeric(job.elapsed_seconds, 1)} 秒` : "等待运行时间…";
  if (finite(job.progress)) $("job-progress").value = Math.min(100, Math.max(0, job.progress)); else $("job-progress").removeAttribute("value");
  $("job-note").textContent = job.mode === "deep" ? "深入模式将采样并优化三维构象，耗时取决于分子的复杂程度。" : "标准模式计算二维结构与性质，不生成三维构象。";
}
function schedulePoll(id, generation) { clearTimeout(pollTimer); pollTimer = setTimeout(() => pollJob(id, generation), 800); }
function finishJob() { clearTimeout(pollTimer); pollTimer = null; activeJob = null; jobGeneration++; setBusy(false); loadHistory(); }
function handleJob(job, generation) {
  if (generation !== jobGeneration) return;
  activeJob = job.id; showJobProgress(job); $("cancel-job").disabled = false;
  if (job.status === "completed") {
    if (!job.result) { showMessage("analysis-error", "任务已结束，但没有返回分析结果。请从最近计算重新打开。"); setSingleView("empty"); }
    else { lastJobId = job.id; renderResult(job.result, job.smiles || $("smiles").value.trim()); }
    finishJob(); return;
  }
  if (job.status === "failed") { showMessage("analysis-error", job.error?.message || "本次分析未完成，请检查结构后重试。"); setSingleView("empty"); finishJob(); return; }
  if (job.status === "cancelled") { setSingleView("empty"); finishJob(); toast("计算已取消"); return; }
  schedulePoll(job.id, generation);
}
async function pollJob(id, generation) {
  if (generation !== jobGeneration) return;
  try { const job = await fetchJSON(`/analyses/${encodeURIComponent(id)}`, { signal: AbortSignal.timeout(15000) }); if (generation !== jobGeneration) return; pollFailures = 0; handleJob(job, generation); }
  catch (error) {
    if (generation !== jobGeneration) return;
    if (++pollFailures < 4) { $("job-note").textContent = "服务连接暂时中断，正在重新读取任务状态…"; schedulePoll(id, generation); }
    else { showMessage("analysis-error", "暂时无法读取任务状态。后台计算可能仍在继续，可刷新最近计算后重新打开。"); setSingleView("empty"); finishJob(); }
  }
}
$("analysis-form").addEventListener("submit", async (event) => {
  event.preventDefault(); if (singleBusy) return;
  const smiles = $("smiles").value.trim(), mode = document.querySelector('input[name="mode"]:checked').value, properties = [...document.querySelectorAll('input[name="properties"]:checked')].map((input) => input.value);
  hideMessage("analysis-error"); $("smiles").removeAttribute("aria-invalid");
  if (!smiles || smiles.length > 2000) { showMessage("analysis-error", smiles ? "单个 SMILES 最多 2,000 字符。" : "请先输入结构或选择一个示例分子。"); $("smiles").setAttribute("aria-invalid", "true"); $("smiles").focus(); return; }
  if (!properties.length) { showMessage("analysis-error", "请至少选择一种预测性质。"); return; }
  invalidateResult(); const generation = ++jobGeneration; activeJob = null; pollFailures = 0; setBusy(true); setSingleView("job"); $("job-stage").textContent = "正在提交任务"; $("job-time").textContent = "等待服务响应…"; $("job-progress").removeAttribute("value");
  try { const job = await fetchJSON("/analyses", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ smiles, properties, mode }), signal: AbortSignal.timeout(20000) }); handleJob(job, generation); loadHistory(); }
  catch (error) { if (generation !== jobGeneration) return; showMessage("analysis-error", networkMessage(error)); setSingleView("empty"); finishJob(); if (error.status === 422) { $("smiles").setAttribute("aria-invalid", "true"); $("smiles").focus(); } }
});
$("cancel-job").addEventListener("click", async () => {
  if (!activeJob) return; const id = activeJob, generation = jobGeneration; $("cancel-job").disabled = true; $("cancel-job").textContent = "正在取消…";
  try { const job = await fetchJSON(`/analyses/${encodeURIComponent(id)}`, { method: "DELETE", signal: AbortSignal.timeout(15000) }); if (generation === jobGeneration) handleJob(job, generation); }
  catch (error) { if (generation === jobGeneration) { toast("取消请求未确认，请重试或等待任务结果。"); $("cancel-job").disabled = false; } }
  finally { $("cancel-job").textContent = "取消计算"; }
});
function resultName(result, input) { return Object.values(examples).find((example) => example.smiles === input)?.name || (typeof result.molecule === "string" && result.molecule !== input && result.molecule !== result.canonical_smiles && result.molecule.length < 70 ? result.molecule : result.molecule?.name) || "分子分析结果"; }
function renderResult(result, input) {
  lastResult = result; const name = resultName(result, input); $("molecule-name").textContent = name; $("molecular-formula").textContent = result.formula || "分子式未提供";
  showStructure($("structure-box"), result.structure_svg); $("canonical-smiles").textContent = result.canonical_smiles || input;
  $("prediction-grid").replaceChildren(); Object.keys(targets).filter((id) => result.predictions?.[id]).forEach((id) => $("prediction-grid").append(renderPrediction(id, result.predictions[id])));
  renderDescriptors(result.descriptor || {}); renderMeasurements(result.measurements); hideMessage("prediction-warnings");
  if (Array.isArray(result.warnings) && result.warnings.length) { const list = node("ul"); result.warnings.forEach((warning) => list.append(node("li", "", stringValue(warning)))); $("prediction-warnings").append(list); $("prediction-warnings").hidden = false; }
  $("raw-result").textContent = JSON.stringify(result, null, 2); setSingleView("result"); renderConformers(result.conformers); announce(`${name}分析完成。`);
}
function switchViewer(mode) {
  const is3d = mode === "3d" && currentConformers?.items?.length; $("structure-box").hidden = Boolean(is3d); $("canvas-wrap").hidden = !is3d;
  $("view-2d").classList.toggle("active", !is3d); $("view-3d").classList.toggle("active", Boolean(is3d)); $("view-2d").setAttribute("aria-pressed", String(!is3d)); $("view-3d").setAttribute("aria-pressed", String(Boolean(is3d)));
  $("viewer-caption").textContent = is3d ? "力场优化构象 · 非量化计算" : "RDKit 二维结构"; if (is3d) viewer.render();
}
function renderConformers(data) {
  currentConformers = data || null; selectedConformer = 0; hideMessage("conformer-warning"); $("conformer-select").replaceChildren(); $("energy-chart").replaceChildren();
  const items = Array.isArray(data?.items) ? data.items : [], available = items.length > 0;
  $("view-3d").disabled = !available; $("export-sdf").disabled = !available || !data.sdf_url; $("conformer-section").hidden = !available;
  if (!available) { switchViewer("2d"); const reasons = { standard_mode: "标准分析未生成三维构象。选择“深入分析”后可探索构象与相对能量。", no_converged_conformers: "本次没有得到收敛构象，未展示三维坐标或最低能量。", too_large_or_flexible: "分子过大或过于柔性，超出本地构象探索的规模范围。性质预测结果单独保留。", too_few_heavy_atoms: "该结构的重原子太少，未进行构象探索。", too_many_atoms_with_hydrogens: "补全氢原子后的规模超过构象计算上限。", missing_mmff94_parameters: "当前分子力场缺少此结构所需的参数，未生成三维构象。", embedding_failed: "没有生成可用的初始三维构象。", optimization_failed: "构象优化未完成，未报告最低能量。", worker_failed: "构象计算未完成，性质预测结果单独保留。", calculation_error: "构象计算遇到问题，未生成可用结果。" }; showMessage("conformer-warning", reasons[data?.reason] || (data?.status === "timeout" ? "构象计算超过时间上限，性质预测结果仍单独保留。" : "本次未返回可用三维构象。")); return; }
  $("conformer-method").textContent = `${data.method || "构象优化"} · ${items.length} 个已收敛构象 · 坐标 Å`;
  items.forEach((item, index) => { const option = node("option", "", `构象 ${index + 1} · ${numeric(item.relative_energy_kcal_mol, 2)} kcal/mol`); option.value = String(index); $("conformer-select").append(option); });
  const energies = items.map((item) => item.relative_energy_kcal_mol), maximum = Math.max(.001, ...energies.filter(finite));
  items.forEach((item, index) => { const button = node("button", "energy-bar"), bar = node("i"); button.type = "button"; button.dataset.index = String(index); button.title = `构象 ${index + 1}：相对能量 ${numeric(item.relative_energy_kcal_mol, 3)} kcal/mol`; button.setAttribute("aria-label", button.title); bar.style.height = `${finite(item.relative_energy_kcal_mol) ? 4 + Math.max(0, item.relative_energy_kcal_mol) / maximum * 58 : 4}px`; button.append(node("strong", "", numeric(item.relative_energy_kcal_mol, 2)), bar, node("span", "", `C${index + 1}`)); button.addEventListener("click", () => selectConformer(index)); $("energy-chart").append(button); });
  $("conformer-note").textContent = `能量相对本分子已收敛集合的最低值，单位 kcal/mol；不是实验构象占比，也不能与其他分子的力场能量比较。${finite(data.unconverged_count) && data.unconverged_count > 0 ? `另有 ${data.unconverged_count} 个未收敛构象未展示。` : ""}`;
  if (data.warnings?.length) showMessage("conformer-warning", data.warnings.map(stringValue).join(" "));
  const ready = selectConformer(0); viewer.reset(); switchViewer(ready ? "3d" : "2d");
}
function selectConformer(index) {
  const item = currentConformers?.items?.[index]; if (!item) return; selectedConformer = index; $("conformer-select").value = String(index);
  $("energy-chart").querySelectorAll("button").forEach((button) => { const selected = Number(button.dataset.index) === index; button.classList.toggle("active", selected); button.setAttribute("aria-pressed", String(selected)); });
  try { viewer.setStructure(currentConformers.atoms, currentConformers.bonds, item.coordinates); return true; } catch (error) { showMessage("conformer-warning", error.message); $("view-3d").disabled = true; switchViewer("2d"); return false; }
}
$("conformer-select").addEventListener("change", () => selectConformer(Number($("conformer-select").value)));
$("view-2d").addEventListener("click", () => switchViewer("2d")); $("view-3d").addEventListener("click", () => switchViewer("3d"));
$("zoom-in").addEventListener("click", () => viewer.scale(1.2)); $("zoom-out").addEventListener("click", () => viewer.scale(1 / 1.2)); $("reset-view").addEventListener("click", () => viewer.reset());
$("copy-smiles").addEventListener("click", async () => { try { await navigator.clipboard.writeText($("canonical-smiles").textContent); toast("标准化 SMILES 已复制"); } catch { const range = document.createRange(); range.selectNodeContents($("canonical-smiles")); const selection = getSelection(); selection.removeAllRanges(); selection.addRange(range); toast("已选中结构表达式，请手动复制"); } });
$("export-json").addEventListener("click", () => { if (lastResult) saveFile(JSON.stringify(lastResult, null, 2), "application/json;charset=utf-8", `molecule-analysis-${timestamp()}.json`); });
$("export-sdf").addEventListener("click", async () => {
  const path = currentConformers?.sdf_url; if (!path) return;
  try { const url = new URL(path, location.origin); if (url.origin !== location.origin) throw new Error("构象下载地址不属于本地服务。"); const response = await fetch(url); if (!response.ok) throw new Error("构象文件暂不可用，请稍后重试。"); saveFile(await response.text(), "chemical/x-mdl-sdfile;charset=utf-8", `molecule-conformers-${timestamp()}.sdf`); } catch (error) { toast(error.message); }
});
async function loadHistory() {
  const generation = ++historyGeneration;
  try { const data = await fetchJSON("/analyses?limit=12", { signal: AbortSignal.timeout(12000) }); if (generation !== historyGeneration) return; $("history-list").replaceChildren(); const items = Array.isArray(data.items) ? data.items : [];
    if (!items.length) $("history-list").append(node("p", "small-empty", "还没有计算记录，开始第一次分析吧。"));
    items.forEach((job) => { const button = node("button", "history-item"); button.type = "button"; const detail = node("div"), statuses = { queued: "排队中", running: "计算中", completed: "已完成", failed: "未完成", cancelled: "已取消" }; detail.append(node("strong", "", job.molecule || job.smiles || job.input?.smiles || `计算 ${job.id.slice(0, 8)}`), node("small", "", `${job.mode === "deep" ? "深入" : "标准"} · ${statuses[job.status] || job.status}${job.created_at ? ` · ${new Date(job.created_at).toLocaleTimeString("zh-CN", { hour: "2-digit", minute: "2-digit" })}` : ""}`)); button.append(node("i"), detail, node("span", "", "↗")); button.addEventListener("click", () => openHistory(job.id)); $("history-list").append(button); });
  } catch { if (generation === historyGeneration) $("history-list").replaceChildren(node("p", "small-empty", "暂时无法读取计算记录。")); }
}
async function openHistory(id) {
  if (singleBusy) { toast("请先等待当前计算完成，或取消当前计算。"); return; }
  const generation = ++jobGeneration; selectPage("analyze"); selectAnalysisTab("single"); invalidateResult(); setBusy(true); setSingleView("job"); $("job-stage").textContent = "正在读取计算记录";
  try { const job = await fetchJSON(`/analyses/${encodeURIComponent(id)}`, { signal: AbortSignal.timeout(15000) }); if (generation !== jobGeneration) return; const input = job.smiles || job.input?.smiles || job.result?.canonical_smiles; if (input) { $("smiles").value = input; updateInputCount(); } const modeInput = document.querySelector(`input[name="mode"][value="${job.mode === "deep" ? "deep" : "standard"}"]`); modeInput.checked = true; const properties = job.properties || (job.result?.predictions ? Object.keys(job.result.predictions) : Object.keys(targets)); document.querySelectorAll('input[name="properties"]').forEach((checkbox) => { checkbox.checked = properties.includes(checkbox.value); }); handleJob(job, generation); }
  catch (error) { if (generation === jobGeneration) { showMessage("analysis-error", networkMessage(error)); setSingleView("empty"); finishJob(); } }
}
$("refresh-history").addEventListener("click", loadHistory);
$("add-comparison").addEventListener("click", () => {
  if (!lastResult) return; const existing = comparison.findIndex((entry) => entry.result.canonical_smiles === lastResult.canonical_smiles), entry = { name: $("molecule-name").textContent, result: lastResult, jobId: lastJobId };
  if (existing >= 0) { comparison[existing] = entry; toast("已更新该分子的对比结果"); }
  else if (comparison.length >= 4) { toast("最多对比 4 个分子，请先移除一项。"); return; }
  else { comparison.push(entry); toast(`已加入对比（${comparison.length}/4）`); }
  renderComparison();
});
function renderComparison() {
  $("compare-count").textContent = String(comparison.length); const destination = $("comparison-results"); destination.replaceChildren();
  if (!comparison.length) { const empty = node("section", "panel comparison-empty"), button = node("button", "secondary-button", "前往分析工作台 ↗"); button.type = "button"; button.addEventListener("click", () => selectPage("analyze")); empty.append(node("span", "", "▥"), node("h2", "", "还没有加入分子"), node("p", "", "完成分析后，点击“加入对比”。"), button); destination.append(empty); return; }
  const wrap = node("section", "panel table-wrap"), table = node("table", "comparison-table"), body = node("tbody"), header = node("tr"); header.append(node("th", "", "分子"));
  comparison.forEach((entry, index) => { const td = node("td"), title = node("div", "compare-molecule"), remove = node("button", "remove-compare", "×"); remove.type = "button"; remove.setAttribute("aria-label", `移除 ${entry.name}`); remove.addEventListener("click", () => { comparison.splice(index, 1); renderComparison(); }); title.append(node("h3", "", entry.name), remove); const picture = node("div", "compare-svg"); showStructure(picture, entry.result.structure_svg); const reopen = node("button", "text-button", "重新打开 ↗"); reopen.type = "button"; reopen.addEventListener("click", () => entry.jobId ? openHistory(entry.jobId) : null); td.append(title, picture, node("small", "", entry.result.formula || ""), reopen); header.append(td); }); body.append(header);
  Object.keys(targets).forEach((id) => { const row = node("tr"); row.append(node("th", "", targets[id].label)); comparison.forEach((entry) => { const td = node("td"), prediction = entry.result.predictions?.[id]; td.append(node("div", hasPrediction(prediction) ? "compare-value" : "compare-label", hasPrediction(prediction) ? numeric(prediction.value) : prediction?.status === "out_of_domain" ? "超出适用范围" : "未提供"), node("small", "", targets[id].unit), node("small", "", stateLabels[prediction?.status] || "暂无结果")); if (hasPrediction(prediction) && prediction.interval90?.every(finite)) td.append(node("p", "compare-interval", `90% 区间 [${numeric(prediction.interval90[0], 2)}, ${numeric(prediction.interval90[1], 2)}]`)); row.append(td); }); body.append(row); });
  descriptorDefinitions.forEach(([key, label, , unit, digits]) => { const row = node("tr"); row.append(node("th", "", `${label}${unit ? ` / ${unit}` : ""}`)); comparison.forEach((entry) => row.append(node("td", "", numeric(entry.result.descriptor?.[key], digits)))); body.append(row); });
  table.append(body); wrap.append(table); destination.append(wrap);
}
$("clear-comparison").addEventListener("click", () => { comparison.length = 0; renderComparison(); });
updateInputCount(); renderComparison(); loadHistory();
