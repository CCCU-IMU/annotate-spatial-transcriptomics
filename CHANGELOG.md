# Changelog

## 1.4.1 — 2026-07-14

- 修复 Seurat 转换对象中 `Spatial@data` 与原始 `counts` 完全相同却被用于 Wilcoxon DEG 的证据层缺口。
- 新增独立、全基因、带输入/analysis-set/输出哈希的 LogNormalize 验证对象；SCT 聚类对象保持不变且验证对象明确禁止用于聚类。
- 初始 cluster DEG 和最终 label DEG 对 `Spatial` assay 强制校验验证 manifest、对象路径、analysis-set 哈希、层身份及 `data != counts`，不满足时 fail closed。
- 分辨率稳定性同时输出全观测 ARI/AMI 和仅用于宏观排序的 `n>=100` restricted 指标；完整保留微簇迁移表和小簇复核清单。
- 明确微簇不因大小被删除、重命名或跳过 DEG、空间图与稀有谱系/技术复核。
- 新增强制调度作业名 `SAMPLE__Pnn_STAGE[_SCOPE]__Ann`、阶段字典、生成/校验器及 registry/report 字段，使多样本并行进度可直接从调度页面识别。

## 1.4.0 — 2026-07-14

- 羊/绵羊/Ovis + 卵巢上下文现在由机器可读 resolver 自动选择 Seurat R-first；确认 StereoPy `cellbin_PPed` 转换后自动执行固定 SCT/PCA/邻接参数合同。
- 修复单 assay Seurat 检查被 `jsonlite` 自动解包为字符串时未触发固定参数合同的问题，并加入真实单 assay 检查格式回归测试。
- 固定 cellbin runner 现在拒绝未记录的参数漂移并要求可验证 SHA256，不再允许无 `digest` 时写出 `NA` 哈希。
- 将 GSE233801 明确为无配对 count-level 参考时的成年羊体细胞主公共 atlas，并限制其 Oocyte/Theca/Epithelial 自动救回范围。
- 加入 2025 Science 人鼠卵巢、2026 Advanced Science 羊–人图谱和 2025 AJOG 专家综述的跨物种边界/反过度注释规则。
- Atlas 最终阈值必须来自与目标不重叠的 query-like 当前 query held-out anchors；外部参考自分类标记为 diagnostic-only，旧 `medium_high` 路线禁止写回。
- 新增主 Agent 审计/用户沟通 + 每样本一个完整子 Agent 的 cohort 控制板、worker ownership 和状态验证脚本。
- 新增真实自动化回归测试，覆盖羊卵巢自动路由、固定参数、GSE233801 优先级、dotplot 禁止转移、held-out 来源、legacy bypass 与多样本隔离。
- 新增配对单细胞参考策略：统一大类语义但不强制空间数据复刻参考标签数量或细分层级。
- 新增可审计 source-label crosswalk、transfer ceiling、验证脚本和羊卵巢常见参考标签别名表。
- 明确证据优先级为当前 query 全基因锚点与形态、配对同阶段单细胞、同物种公共 atlas、发育期/跨组织/跨物种参考。
- 只有 marker dotplot 时仅更新候选 marker/anti-marker 与命名；具备 count-level 对象后才允许深度匹配和校准映射。
- 配对参考默认只能辅助 broad-only 回归，不能绕过 Oocyte、Theca、Smooth muscle、Pericyte/mural 或其他稀有谱系门，也不能替代最终空间全量 broad DEG/dotplot。

## 1.3.0 — 2026-07-14

- 将 Seurat/full-feature RDS 设为可用时的默认计算主干，复用 Scanpy 的迭代策略而不是参数或旧标签。
- 明确区分论文候选分类、分析 parent pool 与最终发布大类，禁止 pool 名直接成为 cell type。
- 基于 2023–2026 年羊卵巢高质量单细胞/多组学研究更新大类命名、证据权重与浅层亚群策略。
- 羊卵巢默认使用 `Stromal/mesenchymal`，并对 Theca、Smooth muscle、Pericyte/mural 和 Mesenchymal progenitor-like 设置独立证据门。
- 禁止将 `Theca/follicular wall` 和 `Stromal/perivascular` 作为方便的最终兜底大类。
- 新增发布 taxonomy 审计器，并把审计哈希接入 completion gate、用户确认与自动控制器。
- 最终报告必须分开生物学大类与界面/QC/技术/待审状态统计。
- 加入输入快照注册、R-first 工作流和真实羊卵巢 forward test 推导出的通用回归测试规则，不包含样本数据或答案映射。

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
