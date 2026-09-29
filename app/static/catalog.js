"use strict";
let catalogKind = "experimental", experimentalStats = null, quantumStats = null;
function renderCatalogStats() {
  const destination = $("database-stats"); destination.replaceChildren();
  const quantum = catalogKind === "quantum", data = quantum ? quantumStats : experimentalStats;
  if (!data) { destination.append(node("div", "stat-card", "数据库信息正在读取或暂不可用")); return; }
  const rows = quantum ? [["量化计算记录", count(data.records), "独立于实验性质库，不合并计数"], ["规范结构", count(data.unique_structures), "同一结构保留原始来源编号"], ["计算方法", data.method || "未提供", "QM9 · 计算参考，非实验值"]] : [["实验库结构", count(data.molecules), "按标准化结构汇总"], ["实测记录", count(data.observations), "同一结构可有多条记录"], ["来源数据集", count(data.dataset_count ?? data.datasets?.length), "可在模型证据中查看来源"]];
  rows.forEach(([label, value, note], index) => { const card = node("div", "stat-card"), number = node("strong", "", value); if (quantum && index === 2) { number.style.fontSize = "15px"; number.style.overflowWrap = "anywhere"; } card.append(node("span", "", label), number, node("p", "", note)); destination.append(card); });
}
function setExperimentalStats(data) {
  experimentalStats = data; const select = $("database-source"), previous = select.value; select.replaceChildren(node("option", "", "全部来源")); select.firstChild.value = "";
  (data.datasets || []).forEach((source) => { const option = node("option", "", source.name || source.id); option.value = source.id; select.append(option); }); if ([...select.options].some((option) => option.value === previous)) select.value = previous; renderCatalogStats();
}
function updateCatalogHelp() {
  const quantum = catalogKind === "quantum", mode = $("database-mode").value;
  $("database-help").textContent = quantum ? "精确查询来源编号（如 gdb_1）、分子式、SMILES 或已支持的常见名，不提供模糊名称检索。" : mode === "text" ? "支持英文名称、已收录的中文常见名、分子式或 SMILES 文本。" : mode === "exact" ? "输入 SMILES 或已支持的常见名，按规范化结构精确匹配。" : "输入 SMILES 或已支持的常见名，按分子指纹相似程度检索。相似性不等于性质相同。";
}
function selectCatalog(kind) {
  catalogKind = kind; invalidateDatabaseSearch();
  document.querySelectorAll("[data-db-tab]").forEach((button) => { const selected = button.dataset.dbTab === kind; button.setAttribute("aria-selected", String(selected)); button.tabIndex = selected ? 0 : -1; });
  const quantum = kind === "quantum"; document.querySelectorAll(".db-experimental-field").forEach((field) => { field.hidden = quantum; }); $("database-form").classList.toggle("quantum-form", quantum);
  $("database-kind").textContent = quantum ? "量化计算参考 · 非实测" : "实验参考";
  $("database-query-label").textContent = quantum ? "来源编号、分子式或结构" : "名称、分子式或 SMILES";
  $("database-query").placeholder = quantum ? "例如 gdb_1、C2H6O 或 CCO" : "例如 乙醇、caffeine 或 CCO";
  $("database-note").textContent = quantum ? "QM9 是指定量化方法的计算参考，不是实验测量，也不是本工作台对输入结构的新预测。每条记录保留来源编号，跨库结构可能重叠。" : "实验记录与模型输出分别保留。同一结构的多条观测可能来自不同测量条件。";
  renderCatalogStats(); updateCatalogHelp(); searchDatabase(0);
}
bindSwitch("[data-db-tab]", (button) => selectCatalog(button.dataset.dbTab));
function updateDatabasePagination(busy = false) { $("database-pagination").hidden = databaseTotal <= 0; $("database-page-status").textContent = `第 ${Math.floor(databaseOffset / databasePageSize) + 1} / ${Math.max(1, Math.ceil(databaseTotal / databasePageSize))} 页 · 每页 20 条`; $("database-previous").disabled = busy || databaseOffset <= 0; $("database-next").disabled = busy || !databaseVisibleCount || databaseOffset + databaseVisibleCount >= databaseTotal; }
function loadCatalogStructure(smiles) { if (singleBusy) { toast("请先结束当前分析，再载入其他结构。"); return; } setInput(smiles); selectPage("analyze"); selectAnalysisTab("single"); $("smiles").focus(); toast("已载入结构，可选择分析深度后开始计算。"); }
function loadButton(item) { const button = node("button", "secondary-button compact", "载入结构"); button.type = "button"; button.disabled = !item.canonical_smiles; button.addEventListener("click", () => loadCatalogStructure(item.canonical_smiles)); return button; }
function renderDatabaseItems(data, offset, kind) {
  const destination = $("database-results"), items = Array.isArray(data.items) ? data.items : []; destination.replaceChildren(); databaseOffset = finite(data.offset) ? data.offset : offset; databaseTotal = finite(data.total) ? data.total : 0; databaseVisibleCount = items.length;
  $("database-count").textContent = `第 ${items.length ? count(databaseOffset + 1) : "0"}–${items.length ? count(databaseOffset + items.length) : "0"} / ${count(databaseTotal)} 条记录`; updateDatabasePagination();
  if (!items.length) { destination.append(node("p", "small-empty spacious", "没有匹配记录。请检查查询方式，或尝试其他结构与筛选条件。")); return; }
  const quantum = kind === "quantum", { table, body } = makeTable(quantum ? ["来源记录与结构", "量化计算参考", "操作"] : ["分子与结构", "分子式 / 分子量", "性质与来源", "操作"]);
  items.forEach((item) => {
    const row = node("tr"), molecule = node("td"), action = node("td"); molecule.append(node("span", "table-name", quantum ? item.id : item.name || `分子 ${item.id}`), node("code", "", item.canonical_smiles || "结构未提供")); if (quantum) molecule.append(node("small", "", item.formula || ""));
    row.append(molecule);
    if (quantum) {
      const td = node("td"), values = node("div", "quantum-cell"), properties = item.properties || {};
      [["dipole_moment", "偶极矩", "D"], ["polarizability", "极化率", "a₀³"], ["homo", "HOMO", "eV"], ["lumo", "LUMO", "eV"], ["gap", "能隙", "eV"]].forEach(([id, label, unit]) => values.append(node("span", "", `${label} ${numeric(properties[id], 3)} ${data.units?.[id] || unit}`)));
      td.append(values); if (item.quality_note) td.append(node("small", "coverage-caution", stringValue(item.quality_note))); row.append(td);
    } else {
      const description = node("td"), properties = node("td"), badges = node("div", "table-badges"); description.append(node("div", "", item.formula || "—"), node("small", "", finite(item.mw) ? `${numeric(item.mw, 2)} g/mol` : "分子量未提供"));
      (Array.isArray(item.targets) ? item.targets : Object.keys(item.targets || {})).forEach((id) => badges.append(node("span", "", targets[id]?.short || id))); properties.append(badges);
      if (Array.isArray(item.sources) && item.sources.length) properties.append(node("small", "", item.sources.map(stringValue).join(" · ")));
      if (finite(item.similarity)) properties.append(node("small", "", `结构相似度 ${numeric(item.similarity)}`)); row.append(description, properties);
    }
    action.append(loadButton(item)); row.append(action); body.append(row);
  }); destination.append(table);
}
async function searchDatabase(offset = 0) {
  if (databaseController) databaseController.abort("superseded");
  const kind = catalogKind, queryValue = $("database-query").value.trim(), mode = $("database-mode").value;
  if (kind === "experimental" && mode !== "text" && !queryValue) { showMessage("database-error", "精确结构与相似结构检索需要先输入一个 SMILES 或常见分子名。"); return; }
  const controller = new AbortController(); databaseController = controller; const timeout = setTimeout(() => controller.abort("timeout"), 45000);
  $("search-button").disabled = true; $("search-button").textContent = "检索中…"; hideMessage("database-error"); updateDatabasePagination(true); $("database-results").setAttribute("aria-busy", "true");
  const query = new URLSearchParams({ q: queryValue, limit: "20", offset: String(offset) });
  if (kind === "experimental") { query.set("mode", mode); if ($("database-source").value) query.set("source", $("database-source").value); if ($("database-target").value) query.set("target", $("database-target").value); }
  try { const data = await fetchJSON(`${kind === "quantum" ? "/quantum/search" : "/database/advanced"}?${query}`, { signal: controller.signal }); if (databaseController !== controller || controller.signal.aborted || kind !== catalogKind) return; renderDatabaseItems(data, offset, kind); databaseLoaded = true; announce(`检索完成：${$("database-count").textContent}`); }
  catch (error) { if (databaseController === controller && controller.signal.reason !== "superseded") { showMessage("database-error", networkMessage(error, controller)); $("database-results").replaceChildren(); $("database-count").textContent = "本次检索未完成"; databaseTotal = 0; databaseVisibleCount = 0; } }
  finally { clearTimeout(timeout); if (databaseController === controller) { databaseController = null; $("search-button").disabled = false; $("search-button").textContent = "检索 ↗"; $("database-results").setAttribute("aria-busy", "false"); updateDatabasePagination(); } }
}
function invalidateDatabaseSearch() { if (databaseController) { databaseController.abort("superseded"); databaseController = null; } databaseLoaded = false; databaseOffset = 0; databaseTotal = 0; databaseVisibleCount = 0; hideMessage("database-error"); $("database-results").replaceChildren(); $("database-results").setAttribute("aria-busy", "false"); $("database-count").textContent = "检索条件已更新，点击检索查看结果。"; $("search-button").disabled = false; $("search-button").textContent = "检索 ↗"; updateDatabasePagination(); }
$("database-form").addEventListener("submit", (event) => { event.preventDefault(); searchDatabase(0); }); $("database-query").addEventListener("input", invalidateDatabaseSearch);
["database-mode", "database-source", "database-target"].forEach((id) => $(id).addEventListener("change", () => { invalidateDatabaseSearch(); updateCatalogHelp(); if ($("database-mode").value === "text" || $("database-query").value.trim()) searchDatabase(0); }));
$("database-previous").addEventListener("click", () => searchDatabase(Math.max(0, databaseOffset - 20))); $("database-next").addEventListener("click", () => searchDatabase(databaseOffset + 20));
renderCatalogStats(); updateCatalogHelp();
