"use strict";

const $ = (id) => document.getElementById(id);
const targets = {
  logS: { label: "水溶解度", symbol: "logS", unit: "log₁₀(S / [mol/L])", short: "logS" },
  logD74: { label: "分配系数 · pH 7.4", symbol: "logD₇.₄", unit: "log₁₀(D₇.₄)", short: "logD₇.₄" },
  hydration_free_energy: { label: "水合自由能", symbol: "ΔG hydration", unit: "kcal/mol", short: "水合 ΔG" }
};
const examples = {
  ethanol: { name: "乙醇", smiles: "CCO" },
  aspirin: { name: "阿司匹林", smiles: "CC(=O)Oc1ccccc1C(=O)O" },
  caffeine: { name: "咖啡因", smiles: "Cn1c(=O)c2c(ncn2C)n(C)c1=O" }
};
const stateLabels = { ok: "模型估计", caution: "谨慎参考", out_of_domain: "超出适用范围", unavailable: "模型不可用" };
let lastResult = null;
let batchResult = null;
let singleController = null;
let singleBusy = false;
let batchBusy = false;
let databaseLoaded = false;
let databaseController = null;
const databasePageSize = 20;
let databaseOffset = 0;
let databaseTotal = 0;
let databaseVisibleCount = 0;

function node(tag, className, text) { const element = document.createElement(tag); if (className) element.className = className; if (text !== undefined && text !== null) element.textContent = String(text); return element; }
function finite(value) { return typeof value === "number" && Number.isFinite(value); }
function numeric(value, places = 3) { return finite(value) ? value.toFixed(places) : "未提供"; }
function count(value) { return finite(value) ? value.toLocaleString("zh-CN") : "未提供"; }
function percentage(value) { return finite(value) ? `${(value * 100).toFixed(1)}%` : "未提供"; }
function announce(text) { $("live-status").textContent = text; }
function stringValue(value) { if (typeof value === "string") return value; if (value && typeof value === "object") return value.message || value.reason || value.code || JSON.stringify(value); return String(value ?? ""); }
function safeLink(url, text, className = "source-link") {
  try { const parsed = new URL(url); if (!["https:", "http:"].includes(parsed.protocol)) return null; const link = node("a", className, text); link.href = parsed.href; link.target = "_blank"; link.rel = "noopener noreferrer"; return link; } catch { return null; }
}
function showMessage(id, message) { const element = $(id); element.textContent = message; element.hidden = false; }
function hideMessage(id) { $(id).hidden = true; $(id).replaceChildren(); }
class ApiError extends Error { constructor(message, status, code) { super(message); this.status = status; this.code = code; } }
async function fetchJSON(url, options = {}) {
  const response = await fetch(url, { ...options, headers: { Accept: "application/json", ...(options.headers || {}) } });
  let data;
  try { data = await response.json(); } catch { throw new ApiError("服务返回的数据无法读取，请检查后端是否正常运行。", response.status); }
  if (!response.ok || (data.error && !data.id)) throw new ApiError(data.error?.message || (typeof data.detail === "string" ? data.detail : "请求未完成，请检查输入或稍后重试。"), response.status, data.error?.code);
  return data;
}
function networkMessage(error, controller) {
  if (error?.name === "TimeoutError") return "服务响应超时，请检查服务状态后重试。";
  if (controller?.signal.aborted) return "请求等待时间过长，请检查后端服务后重试。";
  return error instanceof TypeError ? "无法连接服务，请确认本地后端已启动。" : error.message || "请求未完成，请稍后重试。";
}
function saveFile(content, mime, filename) {
  const url = URL.createObjectURL(new Blob([content], { type: mime })), link = node("a");
  link.href = url; link.download = filename; document.body.append(link); link.click(); link.remove(); setTimeout(() => URL.revokeObjectURL(url), 1000);
}
function timestamp() { return new Date().toISOString().slice(0, 19).replace(/[:T]/g, "-"); }
function makeTable(headers, className = "") { const table = node("table", className), head = node("thead"), row = node("tr"), body = node("tbody"); headers.forEach((label) => { const th = node("th", "", label); th.scope = "col"; row.append(th); }); head.append(row); table.append(head, body); return { table, body }; }

const descriptorDefinitions = [
  ["MW", "分子量", "MW", "g/mol", 2], ["LogP", "中性分配系数", "LogP", "", 2], ["TPSA", "极性表面积", "TPSA", "Å²", 2],
  ["HBD", "氢键供体", "HBD", "", 0], ["HBA", "氢键受体", "HBA", "", 0], ["RotatableBonds", "可旋转键", "RB", "", 0],
  ["RingCount", "环数量", "Rings", "", 0], ["HeavyAtoms", "重原子数", "HA", "", 0], ["AromaticFraction", "芳香原子占比", "AP", "", 3]
];
function renderDescriptors(descriptor) {
  $("descriptor-grid").replaceChildren();
  descriptorDefinitions.forEach(([key, label, abbreviation, unit, digits]) => {
    const group = node("div"), dt = node("dt", "", label), abbr = node("abbr", "", abbreviation), dd = node("dd", "", finite(descriptor[key]) ? numeric(descriptor[key], digits) : "—");
    abbr.title = key; dt.append(abbr); if (unit) dd.append(node("small", "", unit)); group.append(dt, dd); $("descriptor-grid").append(group);
  });
}
function hasPrediction(prediction) { return prediction && ["ok", "caution"].includes(prediction.status) && finite(prediction.value); }
function renderPrediction(id, prediction) {
  const definition = targets[id] || { label: prediction.label || id, symbol: id, unit: prediction.unit || "" };
  const status = Object.hasOwn(stateLabels, prediction.status) ? prediction.status : "unavailable";
  const card = node("article", `prediction-card ${status}`); card.dataset.target = id;
  card.append(node("h5", "prediction-title", definition.label), node("p", "target-symbol", definition.symbol), node("span", "status-chip", stateLabels[status]));
  const canShow = hasPrediction(prediction), valueLabel = canShow ? numeric(prediction.value) : status === "out_of_domain" ? "不提供可靠预测" : "暂无预测";
  card.append(node("div", `prediction-value${canShow ? "" : " no-value"}`, valueLabel), node("div", "prediction-unit", definition.unit));
  const interval = node("div", "prediction-interval"); interval.append(node("span", "", "90% 校准区间"));
  const intervalValues = prediction.interval90;
  interval.append(document.createTextNode(canShow && Array.isArray(intervalValues) && intervalValues.length === 2 && intervalValues.every(finite) ? `[${numeric(intervalValues[0], 2)}, ${numeric(intervalValues[1], 2)}]` : "未提供")); card.append(interval);
  const evidence = node("div", "prediction-evidence");
  if (prediction.training_match === true) evidence.append(node("div", "", "训练集已见连接结构"));
  else if (prediction.training_match === false) evidence.append(node("div", "", "训练集未见连接结构"));
  if (finite(prediction.nearest_similarity)) evidence.append(node("div", "", `最近结构相似度 ${numeric(prediction.nearest_similarity, 3)}`));
  if (prediction.model_name) evidence.append(node("div", "", prediction.model_name));
  const messages = [...(Array.isArray(prediction.domain_violations) ? prediction.domain_violations : []), ...(Array.isArray(prediction.warnings) ? prediction.warnings : [])].map(stringValue).filter(Boolean);
  if (prediction.training_match === true && !messages.some((message) => /训练|training/i.test(message))) messages.push("此结构已见于训练集，不能用于评估模型对全新结构的泛化能力。");
  if (messages.length) {
    const details = node("details"), summary = node("summary", "", `适用范围与提示（${new Set(messages).size}）`), list = node("ul");
    [...new Set(messages)].forEach((message) => list.append(node("li", "", message))); details.append(summary, list); evidence.append(details);
  }
  if (Array.isArray(prediction.nearest_neighbors) && prediction.nearest_neighbors.length) {
    const details = node("details"), summary = node("summary", "", "查看相似训练结构"), list = node("ul");
    prediction.nearest_neighbors.slice(0, 3).forEach((neighbor) => {
      const li = node("li", "", `${neighbor.name || neighbor.canonical_smiles || "参考结构"} · ${numeric(neighbor.similarity, 3)}`); list.append(li);
    }); details.append(summary, list); evidence.append(details);
  }
  const members = prediction.ensemble?.member_values;
  if (canShow && Array.isArray(members) && members.length) {
    const details = node("details"), summary = node("summary", "", `模型成员（${members.length}）`); details.append(summary);
    members.forEach((member) => { const row = node("div", "member-row"); row.append(node("span", "", member.name), node("b", "", numeric(member.value))); details.append(row); });
    if (finite(prediction.ensemble.spread)) details.append(node("p", "", `成员差异 ${numeric(prediction.ensemble.spread)}（加权标准差）`));
    details.append(node("p", "", members.length === 1 ? "仅一个模型成员，差异为零不表示预测没有误差。" : "成员差异不是经过校准的置信区间。")); evidence.append(details);
  }
  card.append(evidence); return card;
}
function renderMeasurements(measurements) {
  const destination = $("measurements"); destination.replaceChildren();
  const records = Array.isArray(measurements) ? measurements : [];
  $("measurement-count").textContent = `${records.length} 条结构匹配记录`;
  if (!records.length) { destination.append(node("p", "empty-measurement", "当前参考数据库未找到该结构的实测记录。未命中不代表该分子没有已发表实验数据。")); return; }
  records.forEach((record) => {
    const item = node("article", "measurement-item"), top = node("div", "measurement-top"), title = node("strong"), value = node("span", "", `${numeric(record.value)} ${record.unit || targets[record.target]?.unit || ""}`);
    title.append(node("span", "measurement-observed", "实测参考"), document.createTextNode(targets[record.target]?.short || record.target || "性质")); top.append(title, value);
    const source = node("p", "measurement-source", `${record.dataset || "数据集"}${record.source_id ? ` · ${record.source_id}` : ""}${record.quality ? ` · ${stringValue(record.quality)}` : ""}`);
    const link = safeLink(record.source_url, "查看来源 ↗"); if (link) { source.append(document.createTextNode(" · "), link); }
    item.append(top, source); destination.append(item);
  });
}
