# Changelog

## 2.0.3 — 2026-07-22

- 修复完整终态账本的控制器兼容性：cell decision 可绑定 active direct return、active Atlas route 或显式 terminal residual-QC freeze；零计数 Atlas 分区可被哈希审计；精确父亚簇 challenger 与最终 direct return 均纳入 membership closure。
- biological context 统一通过一个解析器兼容 `config/biological_context.json` 与 `config/context.json`，completion、autopilot、report 和 release audit 不再因文件别名产生互相冲突的结果。
- residual-QC validator 默认自动识别 v2 的 `final_state`/`qc_reason` 列，同时保留显式列参数；避免把完整终态 ledger 错读为旧 `annotation_status`。
- 最终 HTML 在每个大类/亚群空间高亮旁直接展示经验证的定义路线、定量支撑、排他验证、空间证据、证据来源与前五个 current-data DEG；上下文别名与外部参考下载链接也可正确解析。autopilot 同时统一识别 canonical `final_report_metadata`/`final_*_DEG` 与旧发布文件名。
- 全细胞 marker 空间图固定使用 `all_analysis_set_observations` 范围并可校验预期 observation 数；报告元数据明确要求在含 pandas 的预检运行时执行。以上均复用既有证据/发布审计，不新增公开完成门。

## 2.0.2 — 2026-07-22

- 初始大类改为“基础区室 + 连续信号记忆”：不再为了补齐最终 taxonomy 强制制造 Epithelial、Vascular 或 Smooth muscle 初始簇；有信息的混合簇归入数值占优的受支持 parent，较弱连贯信号进入 `watch` 并在二次聚类中重建。
- 新增共享的 machine-derived lineage decision table 合同。初始簇和 cohort 亚簇的 winner、runner-up 与 margin 均由 validator 重算，Agent/论文解释不能覆盖数值排序。
- cohort outcome v2 要求逐亚簇两阳性 family、无未解决实质矛盾、正 margin 和 observation-level purity；父级标签只作来源信息，不能把全部亚簇自动 `parent_return`。混合亚簇进入一次 targeted question 或 QC。
- 最终 broad completeness 复核增加本类程序、最强竞争程序、细胞级纯度和空间形态的数值内容校验，阻断低典型 marker、强竞争谱系的大标签。
- 上皮召回改为只写回高纯度重聚类亚簇，禁止由少量 marker 将混合簇或空间组件整体扩张。成熟 Smooth muscle 强制 MYH11/CNN1/ACTG2 核心与 mural 排除；RGS5/PDGFRB/NOTCH3/MCAM/CSPG4 支持的血管壁细胞归入 `Vascular-associated`。
- 不新增公开完成门；`check_completion_gate.py` 是唯一公开入口，内部检查按 input integrity、biological evidence、workflow closure 与 release audit 四阶段组织。

## 2.0.1 — 2026-07-22

- 将本地验证过的独立 SCT+BANKSY 原始计数重建合同同步到发布版：BANKSY 使用 SCT Pearson residuals，导入对象的归一化、降维、聚类和历史标签不再作为当前项目计算输入；全基因 LogNormalize 仅用于 marker 验证。
- 修复零锚点 query-only cohort、压缩 membership、稀疏矩阵求和、坐标重复合并、共享 marker 一对多展开、绝对 dotplot 单列树顺序等 6 类可复现运行缺陷。
- 新增 Atlas 路由映射绑定器，要求全细胞路由输入是校准 target 与 disjoint held-out 当前 query 的精确哈希绑定并集；任意临时合并表 fail closed。
- 将分析范围策略、release sessionInfo、review 目录、非空 workflow timeline 和绝对/归一化 dotplot 双视图纳入 stock builder 与发布审计。
- annotation contract 现在把 workflow/biological profile 与候选目录冻结为项目内副本，避免共享安装目录在另一会话更新后改变在跑项目；调度预检同时校验 scheduler-visible job name 和实际 Python/R 环境依赖。HNAICC AIP 明确使用标准输入提交，completion audit 在终态后登记，避免自注册死锁。
- 保留 v2 的开放谱系架构：全候选 selected-plus-two 扫描、大类内嵌 1–5% 程序、零 census、亚群跨谱系重建、Oocyte 完整通过簇和全细胞 Atlas challenger 均继续作为完成门。

## 2.0.0 — 2026-07-22

- 新增项目级不可变 `annotation_contract.json`，在任何生物学命名前冻结输入快照、表达对象边界、workflow/biological profile、候选谱系目录、全组织初始聚类网格、query 重聚类网格、Atlas 路由和发布 taxonomy；后续证据与完成门按 SHA256 绑定。
- 明确区分“已有 BANKSY 输入的原始候选网格”和“项目内 fresh R-first/query 重聚类网格”，禁止把 Seurat 固定网格反向套到上游 BANKSY，也禁止用子集 Leiden 参数解释全组织 BANKSY 分辨率。
- 将大类发现从少数示例 marker 扩展为所有 review-required 候选的完整 `cluster × candidate × positive-family` 全基因绝对证据矩阵；中心化分数与 one-vs-rest DEG 仅作比较证据，不能单独否定组织主体或弱而相干的谱系。
- 统一为单一全细胞 Atlas router：无原标签的冻结 QC 在类别级校准中高置信、非 OOD 且本体兼容时直接 broad-only 回写；已有大类只比较并对群体性差异触发一次正交复核。校准、crosswalk、状态和 mapping 均哈希绑定。
- 发布 taxonomy 统一血管相关 broad 为 `Vascular-associated`，并强制 broad/fine 层级、历史别名拒绝、QC/low-information/candidate 等非生物学 fine 语义拒绝。
- 新增 typed residual-QC 审计和只读项目结果审计；大规模 QC 必须给出逐类原因与上游召回证据，结果目录中的未提交 ledger、旧 completion PASS 冲突、非终态任务、跨项目表达复用和旧 taxonomy 均 fail closed。
- 提供 v1.10→v2 项目迁移器与非破坏性 taxonomy 迁移器。迁移会使旧完成状态失效并要求重新计算 v2 证据，不会把历史标签 grandfather 为可信结果。
- 删除两个未接入但与现行 Atlas/taxonomy 冲突的旧 adjudicator，避免 Agent 在多个近似脚本之间选择到过时路线。
- 新增 v2 合同、完整 marker-family 矩阵、Atlas 路由、残余 QC、结果审计 schema/validator 及回归测试。

## 1.10.0 — 2026-07-21

- 修复同一 SCT+BANKSY 输入在新版中出现大规模 QC 的规则回归：默认大类必须先声明至少两个显式、互不重叠的阳性 marker family；绝对全基因检出率/伪 bulk 用于判断大类是否存在，中心化 module score 与 one-vs-rest DEG 不得单独否定组织主体连续谱。
- 新增 `validate_broad_marker_family_contract.py`。羊卵巢 Stromal、Granulosa、Vascular、Immune、Epithelial 和 Oocyte 默认候选均补齐显式 marker family，阻止“规则要求两 family、profile 却只有一个 program”的不可满足配置。
- 恢复 Atlas 的状态感知写回：无原生物学标签的冻结 QC 若校准为中/高置信、非 OOD、无本体冲突且位于 Atlas 适用范围，可直接回写 broad-only；已有标签只做一致性比较或触发群体复核。marker/空间证据保留为审计与群体挑战，不再成为重复的逐细胞门槛。
- 新增 `route_calibrated_atlas_by_primary_state.py`，并新增 `validate_residual_qc_audit.py`：残余 QC 达到 10% 或 50,000 时必须回查初始大类分辨率、selected-plus-two-higher 全目录扫描和大簇内嵌信号，不能仅以 Atlas low confidence 结案。
- 卵母细胞规则回退到“严格种子判簇、通过簇完整纳入”的语义。新增 `materialize_oocyte_cluster_membership.py`，并强化 `validate_oocyte_context_boundary.py`，禁止用 strict seed 或 spatial object ID 逐细胞缩小已通过 Oocyte 簇。
- 上皮大类拆分 keratin backbone 与 surface/adhesion support 两个 family；空间连续性按通过 cluster/component 复核，最终纳入完整通过成员，而非仅 marker-positive seeds。
- 依据 2025 *Science* 人鼠卵巢图谱关于内皮与周细胞共同构成卵巢血管系统的定义，将血管相关细胞统一发布为 `Vascular-associated`；保留 `Blood endothelial`、`Lymphatic endothelial`、`Pericyte/mural` 为可选细分，成熟 Smooth muscle 仍独立。

## 1.9.2 — 2026-07-21

- Make state validation understand P00 frozen full-object membership tables containing both `analysis_set` and `excluded_initial_qc`, with compatible `_count` and `_n` policy fields.
- Canonicalize reviewed-mapping state aliases before ledger writes and reject unknown state values, preventing descriptive controller terms from corrupting the formal state schema.

## 1.9.1 — 2026-07-21

- Preserve all v1.9 prelabel-evidence freeze fields when reviewed cluster mappings are committed.
- Add a regression test proving that evidence artifact, hash, winner, runner-up and margin survive the state-ledger write.

## 1.9.0 — 2026-07-21

- Added a fail-closed project-local derived-expression registry and validator binding every query evidence object to project/sample, raw-count hash, analysis-set hash, parent artifact and purpose; cross-project derived expression is restricted to explicit Atlas/reference channels.
- Added BANKSY whole-tissue broad-resolution selection evidence based on complete lineage recall, zero-census review, large-cluster purity, DEG/marker coherence, spatial morphology and adjacent-resolution migration; cluster count is shortlist-only.
- Added post-Atlas query-derived broad-class completeness validation for both present labels and zero-census default tissue lineages; Atlas cannot establish absence.
- Added observation-level embedded-lineage spatial component screening to retain coherent 1–5% programs hidden inside large labels.
- Separated canonical Oocyte recall membership from evidence-only spatial context and added multichannel pregranulosa safeguards.
- Added bounded graph-sensitivity adjudication in which more clusters alone cannot rescue or name a lineage.
- Added regression fixtures for cross-project expression leakage, weak-signal memory, large-label dilution, zero-census review, graph sensitivity, Oocyte context and open-world cross-lineage reconstruction.

## 1.8.0 — 2026-07-20

- Separated biological naming thresholds from signal memory: subthreshold but coherent lineage evidence is retained as `watch` rather than silently treated as absence.
- Added continuous full-catalog lineage scans at whole-tissue, broad-class and targeted-reclustering boundaries, covering the selected resolution plus the next two higher available candidates.
- Added immutable lineage boundary/signal registries and a fail-closed coverage validator for missing catalog products, unexplained programs, large-label purity audits and unresolved signals.
- Prohibited parent broad labels from narrowing later candidate searches; cohort subclusters remain open to cross-lineage reconstruction.
- Added v1.7-to-v1.8 migration that creates empty registries and intentionally blocks completion until historical evidence is backfilled.
- Added fixed-point, label-independent all-cell canonical-marker spatial panels to the pre-confirmation report contract, including explicit missing-marker records.
- Added regression tests for complete negative scans, missing lineage products and positive evidence incorrectly recorded as absent.

## 1.7.0 — 2026-07-16

- Added a mandatory label-blind prelabel broad-evidence freeze with all-candidate positive/anti programs, winner/runner-up margin and contradiction checks to prevent paper-marker anchoring.
- Replaced the new-project residual-QC-only Atlas controller with one calibrated all-cell broad mapping after terminal QC freeze; only the frozen QC subset can write labels directly.
- Added deterministic all-cell concordance routing for QC writeback/reject, defined-label agreement, weak challenge, material disagreement, ontology conflict and coherent OOD/unknown candidates.
- Required one evidence-bound orthogonal review for every material broad disagreement, material ontology conflict or coherent OOD group; Atlas-only overwrite and fine-label transfer fail closed.
- Added reusable fixed-transform/reference/index guidance and prohibited dense query-by-reference distance matrices, per-sample joint Atlas retraining and whole-object RCTD defaults.
- Added v1.6.1-to-v1.7.0 migration, new schemas/validators, state fields, completion dependencies, report surfaces and regression tests.

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
