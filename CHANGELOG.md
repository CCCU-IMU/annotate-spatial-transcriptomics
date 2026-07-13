# Changelog

## 1.2.1 — 2026-07-14

- 明确 StereoPy `cellbin_PPed` 转换 RDS 是原始计数载体，而不是已完成 SCT 的分析对象。
- 加入可复用的 Seurat 全组织 SCT 前处理脚本，冻结同批次 QC、SCT、PCA、邻接图和候选分辨率参数并输出 manifest。
- 同步 Seurat 池子重聚类的 SCT v2/`glmGamPoi`、抽样上限、内存策略、cosine Annoy 和 UMAP 参数。
- 保留最终分辨率与生物学注释的样本/池子自适应选择；小簇仅标记审查，不自动重分配。

## 1.2.0 — 2026-07-14

- 将迭代控制器和多路线待定义细胞处理设为完成门的核心要求。
- 加入平衡锚点重聚类、局部界面审查、完整 QC 池重聚类和校准 atlas 多通道救回。
- 明确 atlas moderate-or-higher/high 是 held-out precision 校准目标，而非固定 raw-score 阈值。
- 加入 strict、inclusive、display 三视图及不可变 cell/pool/route 状态追踪。
- 加强卵母细胞等稀有谱系的 marker、anti-marker、空间局灶和重聚类上下文门。
- 最终报告要求同时输出大类与亚群的 canonical/data-specific 树状 dotplot。
- 加入最终用户确认门，确认前不生成耗时的最终发布资产。
- 完成 BANKSY 独立 forward test 和发布包脱敏审计。
