# Changelog

## 1.3.0 — 2026-07-14

- 将 Seurat/full-feature RDS 设为可用时的默认计算主干，复用 Scanpy 的迭代策略而不是参数或旧标签。
- 明确区分论文候选分类、分析 parent pool 与最终发布大类，禁止 pool 名直接成为 cell type。
- 基于 2023–2026 年羊卵巢高质量单细胞/多组学研究更新大类命名、证据权重与浅层亚群策略。
- 羊卵巢默认使用 `Stromal/mesenchymal`，并对 Theca、Smooth muscle、Pericyte/mural 和 Mesenchymal progenitor-like 设置独立证据门。
- 禁止将 `Theca/follicular wall` 和 `Stromal/perivascular` 作为方便的最终兜底大类。
- 新增发布 taxonomy 审计器，并把审计哈希接入 completion gate、用户确认与自动控制器。
- 最终报告必须分开生物学大类与界面/QC/技术/待审状态统计。
- 加入输入快照注册、R-first 工作流和真实羊卵巢 forward test 推导出的通用回归测试规则，不包含样本数据或答案映射。

## 1.2.0 — 2026-07-14

- 将迭代控制器和多路线待定义细胞处理设为完成门的核心要求。
- 加入平衡锚点重聚类、局部界面审查、完整 QC 池重聚类和校准 atlas 多通道救回。
- 明确 atlas moderate-or-higher/high 是 held-out precision 校准目标，而非固定 raw-score 阈值。
- 加入 strict、inclusive、display 三视图及不可变 cell/pool/route 状态追踪。
- 加强卵母细胞等稀有谱系的 marker、anti-marker、空间局灶和重聚类上下文门。
- 最终报告要求同时输出大类与亚群的 canonical/data-specific 树状 dotplot。
- 加入最终用户确认门，确认前不生成耗时的最终发布资产。
- 完成 BANKSY 独立 forward test 和发布包脱敏审计。
