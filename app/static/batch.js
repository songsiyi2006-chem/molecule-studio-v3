"use strict";
let csvImportGeneration = 0;
function parseCSV(text) {
  const input = text.replace(/^\uFEFF/, ""), rows = []; let row = [], field = "", quoted = false, afterQuote = false;
  for (let i = 0; i < input.length; i++) {
    const char = input[i];
    if (quoted) { if (char === '"') { if (input[i + 1] === '"') { field += '"'; i++; } else { quoted = false; afterQuote = true; } } else field += char; continue; }
    if (char === '"') { if (field || afterQuote) throw new Error("CSV 引号格式不正确：引号必须从字段开头开始。"); quoted = true; continue; }
    if (char === ",") { row.push(field); field = ""; afterQuote = false; continue; }
    if (char === "\n" || char === "\r") { if (char === "\r" && input[i + 1] === "\n") i++; row.push(field); if (row.some((value) => value.length)) rows.push(row); row = []; field = ""; afterQuote = false; continue; }
    if (afterQuote) { if (char === " " || char === "\t") continue; throw new Error("CSV 结束引号后只能是逗号或换行。"); } field += char;
  }
  if (quoted) throw new Error("CSV 中有未闭合的引号。"); row.push(field); if (row.some((value) => value.length)) rows.push(row); return rows;
}
function batchInputs() { return $("batch-smiles").value.split(/\r?\n/).map((line) => line.trim()).filter(Boolean); }
function resetBatch() { batchResult = null; $("batch-json").disabled = true; $("batch-csv").disabled = true; $("batch-results").replaceChildren(node("p", "small-empty spacious", "尚无分析结果")); $("batch-summary").textContent = "每行独立校验，完成后可导出全部状态与数值。"; }
function updateBatchInput() { csvImportGeneration++; $("batch-count").textContent = `${batchInputs().length} / 50`; hideMessage("batch-error"); resetBatch(); }
$("batch-smiles").addEventListener("input", updateBatchInput);
$("batch-examples").addEventListener("click", () => { $("batch-smiles").value = Object.values(examples).map((example) => example.smiles).join("\n"); updateBatchInput(); });
$("import-csv").addEventListener("click", () => $("csv-file").click());
$("csv-file").addEventListener("change", async () => {
  const file = $("csv-file").files[0]; if (!file || batchBusy) return; const generation = ++csvImportGeneration;
  try {
    if (file.size > 2 * 1024 * 1024) throw new Error("CSV 文件超过 2 MB，请先拆分文件。");
    const text = await file.text(); if (generation !== csvImportGeneration || batchBusy) return;
    const rows = parseCSV(text); if (!rows.length) throw new Error("CSV 文件为空。");
    const index = rows[0].findIndex((value) => ["smiles", "canonical_smiles"].includes(value.trim().toLowerCase())); if (index < 0) throw new Error("CSV 第一行需包含 SMILES 或 canonical_smiles 列名。");
    const values = rows.slice(1).map((row, rowIndex) => { if (row.length !== rows[0].length) throw new Error(`CSV 第 ${rowIndex + 2} 行列数不一致。`); return (row[index] || "").trim(); }).filter(Boolean);
    if (!values.length) throw new Error("CSV 的 SMILES 列没有可用结构。"); if (values.length > 50) throw new Error("每批最多 50 个结构，请拆分 CSV 后再导入。"); if (values.some((value) => /[\r\n]/.test(value))) throw new Error("SMILES 字段包含换行，无法作为单行结构导入。其他带引号字段可包含换行。");
    $("batch-smiles").value = values.join("\n"); updateBatchInput(); toast(`已导入 ${values.length} 个结构`);
  } catch (error) { showMessage("batch-error", error.message); } finally { $("csv-file").value = ""; }
});
function renderBatch(data) {
  batchResult = data; const items = Array.isArray(data.items) ? data.items : [], success = items.filter((item) => item.result?.valid && !item.error).length;
  $("batch-summary").textContent = `共 ${items.length} 个结构，${success} 个完成分析，${items.length - success} 个未通过校验。各性质仍需独立查看适用范围。`;
  const { table, body } = makeTable(["分子", "logS", "logD₇.₄", "水合 ΔG / kcal/mol"], "batch-table");
  items.forEach((item, index) => { const row = node("tr"), input = node("td"); input.append(node("strong", "", `#${index + 1}`), node("code", "", item.input || item.result?.canonical_smiles || "—")); row.append(input);
    if (item.error || !item.result) { row.className = "row-error"; const td = node("td", "", stringValue(item.error) || "此行分析未完成"); td.colSpan = 3; row.append(td); }
    else Object.keys(targets).forEach((id) => { const prediction = item.result.predictions?.[id], td = node("td", "", hasPrediction(prediction) ? numeric(prediction.value) : prediction?.status === "out_of_domain" ? "超出范围" : "未提供"); td.append(node("small", "", stateLabels[prediction?.status] || "暂无结果")); row.append(td); }); body.append(row);
  }); $("batch-results").replaceChildren(table); $("batch-json").disabled = false; $("batch-csv").disabled = false; announce("批量分析完成。");
}
$("batch-form").addEventListener("submit", async (event) => {
  event.preventDefault(); if (batchBusy) return; const smiles = batchInputs(); hideMessage("batch-error");
  if (!smiles.length || smiles.length > 50) { showMessage("batch-error", !smiles.length ? "请每行输入一个结构，或导入 CSV。" : "每批最多 50 个分子，请减少输入行数。"); return; }
  resetBatch(); batchBusy = true; ["batch-button", "batch-smiles", "batch-examples", "csv-file"].forEach((id) => { $(id).disabled = true; }); $("batch-button").textContent = "正在批量分析…"; $("batch-summary").textContent = `正在分析 ${smiles.length} 个分子…`;
  try { const data = await fetchJSON("/predict/batch", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ smiles }), signal: AbortSignal.timeout(180000) }); renderBatch(data); }
  catch (error) { showMessage("batch-error", networkMessage(error)); $("batch-summary").textContent = "批量分析未完成，请查看提示后重试。"; }
  finally { batchBusy = false; ["batch-button", "batch-smiles", "batch-examples", "csv-file"].forEach((id) => { $(id).disabled = false; }); $("batch-button").textContent = "批量分析 ↗"; }
});
function csvCell(value) { let text = value === null || value === undefined ? "" : String(value); if (typeof value !== "number" && /^[\s]*[=+\-@\t\r]/.test(text)) text = `'${text}`; return `"${text.replaceAll('"', '""')}"`; }
$("batch-json").addEventListener("click", () => { if (batchResult) saveFile(JSON.stringify(batchResult, null, 2), "application/json;charset=utf-8", `molecule-batch-${timestamp()}.json`); });
$("batch-csv").addEventListener("click", () => {
  if (!batchResult) return; const headers = ["index", "input", "canonical_smiles", "error", ...Object.keys(targets).flatMap((id) => [`${id}_value`, `${id}_unit`, `${id}_status`, `${id}_interval90_low`, `${id}_interval90_high`])];
  const rows = (batchResult.items || []).map((item, index) => [index + 1, item.input || "", item.result?.canonical_smiles || "", item.error ? stringValue(item.error) : "", ...Object.keys(targets).flatMap((id) => { const prediction = item.result?.predictions?.[id], show = hasPrediction(prediction); return [show ? prediction.value : null, prediction?.unit || targets[id].unit, prediction?.status || "unavailable", show && finite(prediction.interval90?.[0]) ? prediction.interval90[0] : null, show && finite(prediction.interval90?.[1]) ? prediction.interval90[1] : null]; })]);
  saveFile("\uFEFF" + [headers, ...rows].map((row) => row.map(csvCell).join(",")).join("\r\n"), "text/csv;charset=utf-8", `molecule-batch-${timestamp()}.csv`);
});
