# Changelog

## 1.6.1 — 2026-07-15

- Removed the last active QC-reclustering prerequisites; Atlas/internal-anchor/marker/spatial consensus now audits only the exact frozen terminal residual-QC partition.
- Added evidence-content schemas and validators for cohort outcomes, direct returns and per-label support, plus an exact cell-membership closure audit bound to completion, master quality and release.
- Replaced pool resolution ranking with question-aware cohort ranking; broad-purity audits permit one compartment and can close successfully as `homogeneous_parent_confirmed` after the complete grid.
- Bound fixed sheep-ovary StereoPy cellbin grids to a hash-verified workflow profile, explicit strategy preset and multi-artifact input provenance instead of species/path/feature hints.
- Made full-feature validation depend on the loaded biological profile and require exact broad/fine support coverage.
- Standardized new-project confidence to `low`, `moderate`, `high`, retained the separate Atlas enum, and added a v1.6.0 migration utility.
- Added structured controller gap codes, business-state exit handling and content-hash dependency manifests for expensive assets.
- Isolated the retired pool/Route controller under `legacy/` and added homogeneous, subtype, state-only, mixed-lineage, closure, deterministic, leakage-safe benchmark and ablation regression tests.

## 1.6.0 — 2026-07-15

- Replaced persistent biological review pools with one default direct-lineage architecture for every new project.
- Added immutable broad-class/targeted cohort and direct-return registries plus a fail-closed workflow validator.
- Removed mandatory QC reclustering before terminal Atlas rescue; moderate-or-higher Atlas returns remain broad-only and low confidence remains QC.
- Added direct cross-lineage return without intermediate target pools or automatic repeat reclustering.
- Reworked Skill entry, R-first guide, state contract, iteration planner, completion gate, master review and report around the same architecture.
- Preserved old pool/branch registries only as migration provenance.
- Declared clean-runner regression-test dependencies, made them mandatory in both PR validation and Release packaging, and removed duplicate PR-branch `push` checks while retaining validation on `main`.

## 1.5.1 — 2026-07-15

- 新增逐样本主 Agent 注释质量审批硬门：仅在池重聚类、适用的 Route A–D、救回写回、单一最终注释及 completion gate 全部完成后触发，禁止在仅完成大类注释时提前审批。
- 审批职责保持轻量：完成门负责机械闭环，主 Agent 只复核大类合理性、marker/anti-marker 与空间形态、文献候选/易混淆谱系安全性，以及是否达到内置脱敏羊卵巢 R-first 参考流程的证据质量；允许通过并保留 concerns。
- 主审批、轻量确认页、用户确认及最终发布均使用哈希绑定；任一 ledger、路线、池、完成门或参考版本变化都会使下游审批失效。
- 多样本状态机增加 `READY_FOR_MASTER_QUALITY_GATE` 与 `MASTER_QUALITY_APPROVED`，每个样本均须独立审批后才可冻结并进入 cohort 用户确认。
- 新增显式 `sheep_ovary_same_batch_rfirst` 策略预设：复用脱敏成功流程的阶段顺序、池轴、证据门、救回顺序和报告合同，但禁止复制参考样本的分辨率、簇标签、membership、比例或亚群目录。
- 预设采用开放式谱系发现而非封闭标签表：新增机器可读的 14 项候选边界目录，覆盖卵泡/生殖系、甾体生成、基质–间充质–收缩/壁细胞、血管/淋巴、免疫、上皮/间皮及神经胶质/神经内分泌家族；逐样本阴性结论有效，目录外未解释多基因程序必须创建额外候选记录，标准池列表明确非穷尽。
- 取消新流程中的通用 rare-cell 路线。样本特异候选使用普通生物学池验证；仅 Oocyte 因 ZP/邻近污染保留按需触发的专门安全路线，旧 `strict_rare_cell_review` 记录只作兼容读取。

## 1.5.0 — 2026-07-15

- 发布结果收敛为一套最终注释：大类至少达到 calibrated moderate-or-higher，亚群仅保留 high-confidence 且可作为 fine anchor 的结果；旧 `strict/inclusive/display` 只保留为历史迁移元数据，不再进入 HTML 或发布资产。
- 报告元数据明确绑定 `primary_broad_label` / `primary_subtype_label`；修复布尔资格列被 Pandas 自动推断后与字符串比较造成最终标签意外置空的问题，并禁止生成 `Broad only: ...` 伪亚群。
- 新增确认前轻量证据 HTML 硬门：完成门通过后先展示大类/亚群支持理由、高区分度大类空间投影和 canonical broad marker dotplot；确认请求与最终发布均绑定该报告、support registry 和底层快照哈希，完整 DEG/树图/空间资产仍只在确认后生成。
- 最终 HTML 采用已验证 Scanpy 报告的信息架构：顶部整体空间注释、可展开大类/亚群树与节点跳转、大类及条件性亚群 canonical/data-specific 树状 dotplot、按细胞类型分组的空间基因图和底部中文详细流程。
- 羊卵巢 whole-tissue、所有生物学池和 QC 池统一执行 `0.1,0.2,0.3,0.4,0.6` 正常候选网格；`0.01/0.02/0.05` 及不完整网格 fail closed，图异常时先修复邻接图而不是继续降分辨率。
- 将 Atlas 范围收紧为“完整 QC-holdout 锚点重聚类后仍为 QC 的冻结残余 membership”。RCTD medium/low 先进入 QC holdout，不能直接调用 Atlas；Atlas route 必须链接前序 QC-anchor route、outcome hash 与完全一致的 residual-QC membership hash。
- GSE233801 继续作为无可用配对 count-level 参考时的成年羊体细胞主公共 Atlas，但只服务 residual-QC broad-only 救回；禁止对全对象、已定义 broad/fine 或普通生物学池常规分类。
- 大类 DEG、canonical/data-specific dotplot、UMAP 与空间统计纳入所有正式写回该大类的直接与 broad-only 救回观测；fine DEG/图只使用高置信真实 fine labels，避免锚点子集偏倚或伪造亚群。
- 新增生成作业预检、事故登记/零开放事故门、profile-role 检查、固定分辨率检查、单一最终注释构建及 v1.4→v1.5 迁移工具，并加入相应发布回归测试。
- 明确主 Agent 负责跨样本统筹、状态展示、质量审计与唯一用户确认；每个样本由一个完整流程子 Agent 独立执行，不能降级为 cluster 重命名或只审计。

## 1.4.4 — 2026-07-15

- 加入羊免疫球蛋白 LOC 常数区别名，避免在空间 cellbin 中因 `PTPRC` 稀疏而系统性漏掉抗体分泌/B-lineage 群。
- 新增受保护的替代浆细胞门：稳定多免疫球蛋白位点程序还必须同时具备 `JCHAIN`、`POU2AF1`、`TENT5C` 或 `MZB1` 等独立 B/浆细胞调控证据；单个 Ig、JCHAIN、CD74 或 MHC 仍不能定类。
- 该路线默认只支持 broad `Immune`，plasma/antibody-secreting 仅作状态标签；除非独立 fine-label 完成门通过，否则不发布浆细胞亚型。
- 修复 pool runner 在 source/state composition 汇总时把循环变量误解析为 data.table closure 的错误；已完成聚类可做后处理修复，无需重复 SCT/Leiden。

## 1.4.3 — 2026-07-15

- 将已完成的羊卵巢 Seurat/R-first 卵母细胞路线脱敏为通用的两层候选策略：完整多模块起始候选池负责召回并进入 query-only 重聚类，严格 marker/anti-program 种子和空间焦点只提供身份支持。
- 修复严格空间焦点可能被误用为唯一重聚类 membership 的过严风险；孤立但通过起始门的皮质/小卵母细胞候选继续保留在完整候选池中。
- 候选池重聚类后必须以非 ZP 身份/母源胞质程序、体细胞反程序和空间对象形态共同判定；颗粒或基质主导的邻近群回流相应体细胞池并保留 ambient/adjacent 状态标签。
- `calibrate_rare_cell_candidates.py` 与 `screen_spatial_foci.py` 现在同时输出完整候选重聚类 membership 和严格种子/焦点证据 membership，并明确 cellbin 数量不等同于生物学卵母细胞数量。
- 新增发布回归测试，验证完整候选池不会被严格种子缩减，并保持 zona-only、位置-only 和未经过体细胞分流的 Oocyte 判定 fail closed。
- 修复稳定性比较器把 Seurat/pool runner 的 TSV cluster membership 当作 CSV 读取的问题；CSV/TSV/压缩 TSV 现在按显式后缀选择分隔符，已完成的聚类不再因轻量后处理误报失败。

## 1.4.2 — 2026-07-14

- 修复 Seurat 5.5 候选 Leiden 网格在默认 sequential `future` plan 下逐分辨率串行、却申请大量空闲 CPU 的资源错配。
- 全组织和 pool runner 新增 `--resolution-workers`/`--resolution-future-plan`，在 Linux 上以多个 `future` worker 并行不同 resolution，并将同一 worker 数传给 `uwot` UMAP；manifest 记录请求/实际 worker、backend 和线程数。
- 新增强制 CPU–parallelism 合同：单线程 Leiden/直接 `FindAllMarkers` 不得占用多核；可拆为独立 resolution/cluster 作业，并用依赖汇总；完成后审计 CPU time / wall time。
- 明确皮质、皮质下、外周或切片边缘位置不得作为 Oocyte 的负证据或降权项；小/原始卵母细胞可以位于皮质。位置本身仍不能替代多模块 marker、体细胞 anti-program、空间对象和严格重聚类证据。

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
