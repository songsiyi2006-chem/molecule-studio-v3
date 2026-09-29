"use strict";
function detailRow(dl, label, value) { const row = node("div"); row.append(node("dt", "", label), node("dd", "", value ?? "未提供")); dl.append(row); }
function renderEvidenceModels(data) {
  const destination = $("model-cards"); destination.replaceChildren(); const models = Array.isArray(data.models) ? data.models : [];
  if (!models.length) { const empty = node("section", "panel"); empty.append(node("p", "small-empty spacious", "暂无可读取的模型评估。")); destination.append(empty); return; }
  models.forEach((model) => {
    const id = model.id || model.target, unavailable = model.available === false || model.status === "unavailable", card = node("article", "panel model-card");
    card.append(node("span", "overline", unavailable ? "MODEL UNAVAILABLE" : "INTERNAL HOLDOUT EVALUATION"), node("h3", "", targets[id]?.label || model.label || id), node("p", "", `${model.model_name || "模型名称未提供"}\n${model.version || data.version || "版本未提供"}`));
    if (unavailable) { card.append(node("p", "model-note", "该模型未加载，暂时不提供预测或评估数值。")); destination.append(card); return; }
    const test = model.metrics?.test || {}, calibration = model.calibration || {}, split = model.split_counts || {}, metric = node("div", "model-metric"); metric.append(node("strong", "", finite(test.rmse) ? numeric(test.rmse) : "—"), node("span", "", "测试集 RMSE")); card.append(metric);
    const dl = node("dl"); detailRow(dl, "单位", model.unit || targets[id]?.unit); detailRow(dl, "测试 MAE", numeric(test.mae)); detailRow(dl, "测试 R²", numeric(test.r2)); detailRow(dl, "均值基线 RMSE", numeric(model.metrics?.baseline?.rmse)); detailRow(dl, "选中方案验证 RMSE", numeric(model.metrics?.validation_selected?.rmse)); card.append(dl);
    const cal = node("div", "calibration-box"), calRows = node("dl"); cal.append(node("strong", "", "区间覆盖与校准")); detailRow(calRows, "名义覆盖率", percentage(calibration.nominal_coverage)); detailRow(calRows, "测试实际覆盖率", percentage(calibration.test_coverage)); detailRow(calRows, "区间半径", finite(calibration.radius) ? `± ${numeric(calibration.radius)}` : "未提供"); detailRow(calRows, "平均区间宽度", numeric(calibration.mean_width)); cal.append(calRows);
    if (finite(calibration.test_coverage) && finite(calibration.nominal_coverage) && calibration.test_coverage < calibration.nominal_coverage) cal.append(node("p", "coverage-caution", "实际覆盖率低于名义目标，不能将该区间理解为已达到 90% 可靠性。")); card.append(cal);
    const names = model.ensemble?.member_names, weights = model.ensemble?.weights;
    if (Array.isArray(names) && names.length) { const selected = node("details"); selected.append(node("summary", "", `选中的模型成员（${names.length}）`)); names.forEach((name, index) => { const row = node("div", "member-row"); row.append(node("span", "", name), node("b", "", finite(weights?.[index]) ? `${numeric(weights[index] * 100, 1)}%` : "权重未提供")); selected.append(row); }); selected.append(node("p", "model-note", "成员与组合权重由验证数据选择；成员差异不是校准不确定度，也不是性能提升的证明。")); card.append(selected); }
    const candidates = model.candidate_validation || model.candidates;
    if (Array.isArray(candidates) && candidates.length) {
      const details = node("details"), { table, body } = makeTable(["候选模型", "验证 RMSE"], "model-candidates"); details.append(node("summary", "", `候选评估（${candidates.length} 项）`));
      candidates.forEach((candidate) => { const row = node("tr"); row.append(node("td", "", candidate.model_name || candidate.name || "候选模型"), node("td", "", numeric(candidate.rmse ?? candidate.metrics?.rmse))); body.append(row); }); details.append(table); card.append(details);
    }
    const partitions = node("details"), parts = node("dl"); partitions.append(node("summary", "", "数据划分与选模边界")); detailRow(parts, "训练分区", count(split.train)); detailRow(parts, "选模验证分区", count(split.validation)); detailRow(parts, "最终拟合样本", count(model.fit_count)); detailRow(parts, "独立校准分区", count(split.calibration)); detailRow(parts, "测试分区", count(split.test)); detailRow(parts, "结构身份重叠", count(model.identity_overlap)); detailRow(parts, "骨架重叠", count(model.scaffold_overlap)); if (finite(model.quarantined_from_fit)) detailRow(parts, "隔离记录", count(model.quarantined_from_fit)); partitions.append(parts);
    partitions.append(node("p", "model-note", "模型与组合方案在验证集上选择，之后重拟合训练与验证数据。校准和测试标签不参与选择。内部骨架留出仍不能替代外部验证。")); card.append(partitions);
    card.append(node("p", "model-note", "以上测试指标包含全部留出样本，包括应用可能拒报的范围外样本。不同性质的误差不能横向比较。")); destination.append(card);
  });
  renderVersionComparisons(models);
}
function renderVersionComparisons(models) {
  const destination = $("comparison-evidence"); destination.replaceChildren(); const available = models.filter((model) => model.v2_same_test_comparison);
  if (!available.length) return;
  const panel = node("section", "panel comparison-panel"), heading = node("div", "panel-heading"); heading.append(node("h2", "", "V2 → V3：在相同测试数据上比较"), node("span", "overline", "DESCRIPTIVE COMPARISON")); panel.append(heading);
  const wrapper = node("div", "table-wrap"), { table, body } = makeTable(["目标性质", "V2 RMSE", "V3 RMSE", "变化", "测试样本"]);
  available.forEach((model) => { const comparison = model.v2_same_test_comparison, change = comparison.rmse_change, row = node("tr"); row.append(node("td", "", targets[model.id || model.target]?.label || model.label), node("td", "", numeric(comparison.v2_metrics?.rmse)), node("td", "", numeric(comparison.v3_metrics?.rmse)), node("td", "", finite(change) ? `${change > 0 ? "+" : ""}${numeric(change)} · ${change < 0 ? "下降" : change > 0 ? "增加" : "不变"}` : "未提供"), node("td", "", count(comparison.test_count))); body.append(row); }); wrapper.append(table); panel.append(wrapper);
  panel.append(node("p", "evidence-note", "误差下降仅表示这组测试数据上的结果更好。V2 的测试结果此前已被查看，因此本次比较是描述性比较，不是新的前瞻性外部验证，也不保证所有分子都改善。")); destination.append(panel);
}
function renderSources() {
  const destination = $("sources-list"); destination.replaceChildren(); const sources = Array.isArray(experimentalStats?.datasets) ? experimentalStats.datasets : [];
  sources.forEach((source) => { const item = node("div", "source-item"), detail = node("div"); detail.append(node("strong", "", source.name || source.id), node("p", "", `${finite(source.accepted_count) ? `收录 ${count(source.accepted_count)} 条` : "记录数未提供"}${finite(source.raw_count) ? ` · 原始 ${count(source.raw_count)} 条` : ""}${source.license_note ? ` · ${source.license_note}` : ""}`)); item.append(detail); const link = safeLink(source.source_url, "数据来源 ↗"); if (link) item.append(link); destination.append(item); });
  if (quantumStats) { const item = node("div", "source-item"), detail = node("div"); detail.append(node("strong", "", "QM9 · 量化计算参考"), node("p", "", `${count(quantumStats.records)} 条计算记录 · ${quantumStats.method || "计算方法未提供"} · 与实验库独立展示`), node("p", "", quantumStats.note || "量化计算参考不是实验测量。")); item.append(detail); const link = safeLink(quantumStats.source_url, "计算数据来源 ↗"); if (link) item.append(link); destination.append(item); }
  if (!sources.length && !quantumStats) destination.append(node("p", "small-empty", "暂时无法读取数据来源。"));
}
async function loadServiceInfo() {
  const results = await Promise.allSettled([fetchJSON("/health", { signal: AbortSignal.timeout(20000) }), fetchJSON("/model-info", { signal: AbortSignal.timeout(20000) }), fetchJSON("/database/stats", { signal: AbortSignal.timeout(20000) }), fetchJSON("/quantum/stats", { signal: AbortSignal.timeout(20000) })]);
  const [health, models, experimental, quantum] = results;
  if (health.status === "fulfilled") { const data = health.value, ready = ["ok", "healthy"].includes(data.status) && data.model_loaded !== false && data.database_ready !== false; $("connection").className = `connection ${ready ? "online" : "offline"}`; $("connection").querySelector("span").textContent = ready ? "本地服务已就绪" : "部分服务未就绪"; } else { $("connection").className = "connection offline"; $("connection").querySelector("span").textContent = "服务未连接"; }
  if (models.status === "fulfilled") renderEvidenceModels(models.value); else { const panel = node("section", "panel"); panel.append(node("p", "small-empty spacious", "暂时无法读取当前模型信息。")); $("model-cards").replaceChildren(panel); }
  if (experimental.status === "fulfilled") setExperimentalStats(experimental.value);
  if (quantum.status === "fulfilled") { quantumStats = quantum.value; renderCatalogStats(); }
  renderSources();
}
loadServiceInfo();
