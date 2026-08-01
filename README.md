# Annotate Spatial Transcriptomics Skill

[![Validate](https://github.com/CCCU-IMU/annotate-spatial-transcriptomics/actions/workflows/validate.yml/badge.svg)](https://github.com/CCCU-IMU/annotate-spatial-transcriptomics/actions/workflows/validate.yml)

面向 Codex Agent 的空间转录组/单细胞转录组迭代式注释 Skill。v2.5 采用“第一轮只划分 cohort、第二轮重聚类负责正式注释”的哈希绑定控制器：每个初始 cluster 都从项目自身 raw counts 独立执行 SCT/PCA/SNN/Leiden，完整扫描开放候选目录；只有二级 mixed subcluster 才进行局部 observation 拆分。全部 cohort 合并后才冻结 broad，再执行一次固定 GSE233801 Atlas、羊卵巢三终点组织学复核，最后严格串行地对每个可评估 broad 执行独立的过召回、欠召回、分子、空间与文献边界复核。同一时刻只有一个 active cell type；一个类型闭环后才进入下一个。稳定分区默认复用，只有原 source subcluster 确实不可评估时才重聚类。普通 fine/state 候选不能反向重建 broad；只有候选目录显式授权的 parent reconstruction 路线才可在完整父类证据同时通过时提出 broad proposal。项目自编 scorer 或 subset writer 只能产生 experimental 结果，不能写正式标签。

适配 Seurat RDS、AnnData/H5AD、SingleCellExperiment、BANKSY、Scanpy/Leiden、Seurat 聚类和外部 cluster table。空间数据以可靠的大类为主要终点，亚群只在证据充分时定义。

## 一键安装

在运行 Codex 的服务器或工作站执行：

```bash
curl -fsSL https://raw.githubusercontent.com/CCCU-IMU/annotate-spatial-transcriptomics/main/install.sh | bash
```

默认安装到 `${CODEX_HOME:-$HOME/.codex}/skills/annotate-spatial-transcriptomics`，并自动运行完整性检查。重新启动 Codex 或开启一个新任务后即可使用。

安装固定版本：

```bash
curl -fsSL https://raw.githubusercontent.com/CCCU-IMU/annotate-spatial-transcriptomics/main/install.sh \
  | bash -s -- --ref v2.5.0
```

克隆后本地安装（适合内网或需要审查源码的环境）：

```bash
git clone https://github.com/CCCU-IMU/annotate-spatial-transcriptomics.git
cd annotate-spatial-transcriptomics
bash install.sh
```

指定 Codex 目录：

```bash
bash install.sh --dest /path/to/.codex/skills
```

卸载只需删除安装目录：

```bash
rm -rf "${CODEX_HOME:-$HOME/.codex}/skills/annotate-spatial-transcriptomics"
```

## 最佳使用方法：怎样获得高质量结果

高质量来自“上下文完整 + 自主迭代 + 多证据闭环 + 状态不可变”，而不是给 Agent 一张 marker 表后让它一次性重命名 cluster。推荐按下面方式启动。

Seurat 用户要特别注意表达边界：SCT+BANKSY 对象可以提供第一轮空间分区，但每个第二轮 cohort 必须从同一项目的非 SCT raw-count assay 重新执行 `SCT v2/glmGamPoi → PCA → SNN → Leiden`。同一 raw-count assay 可在 cohort 内建立全基因 LogNormalize 数据层用于 DEG/marker 证据；禁止对 SCT corrected counts 再做 SCTransform，也禁止复用其他注释项目的表达对象。极小 cohort 记录为 `underpowered_not_evaluable` 并进入 Atlas/unresolved 路径，不能静默继承第一轮 provisional 标签。

已经完成 SCT+BANKSY 的羊卵巢 Seurat RDS 使用唯一 bootstrap；只给输入、项目目录和样本名即可自动审计非 SCT raw counts、坐标与常见 BANKSY resolution 列，冻结 analysis set、标准 partition grid、固定 GSE233801 Atlas 和唯一 annotation contract。自动推断含糊时才要求显式 `--cluster-mapping`，不再为每个样本临时编写 input audit 或 contract adapter：

```bash
python annotate-spatial-transcriptomics/scripts/bootstrap_sct_banksy_project.py \
  --sample SAMPLE_ID \
  --rds /absolute/path/to/SAMPLE_sct_BANKSY_preprocessed_seurat.rds \
  --project-root /absolute/path/to/SP_ANNO/SAMPLE \
  --observation-unit cellbin \
  --rscript /absolute/path/to/the/validated/Rscript
```

每个正式 controller 调用都会原子写出 `current_stage.json` 和 `next_action_manifest.json`。`REVIEW_REQUIRED`/`ITERATION_REQUIRED` 是可恢复的成功暂停，不再在调度器中伪装成软件崩溃；只有 `FAILED_RUNTIME` 才表示输入、环境、资源或代码失败。bootstrap 还会把非 SCT raw-count assay 和 input-audit SHA 写入 contract；最终归一化 dotplot 使用当前表达 assay，绝对值检测率则固定从 contract 绑定的 raw-count assay 读取。

多样本并行时，每个调度作业必须使用 `样本__P阶段_任务[_队列或目标]__A尝试号`，例如 `SAMPLE1__P10_SCT__A01`、`SAMPLE1__P40_COHORT_stromal__A02`。Skill 内的生成器会限制阶段码和长度，并把名称写入 run registry/报告；禁止继续使用 `sct_preprocess_v0` 这类无法从调度页面判断阶段的名称。

### 1. 首次消息一次给全背景

至少提供：

- 输入文件或输入目录、输出项目目录；
- 物种、组织、发育阶段/处理条件；
- 测序平台与观测单位（真实细胞、cellbin、spot 等）；
- 当前已有的聚类结果及其生成方法；
- 主要生物学问题、需要重点审查的稀有谱系；
- 可用 R/Python 环境、调度系统和计算资源；
- 参考 atlas 的优先级（如已有），但不要把参考标签当作真值。
- 与空间样本对应的单细胞对象、原始注释列、每类细胞数、DEG/marker dotplot 及样本/阶段关系（如已有）。

推荐提示词：

```text
请使用 $annotate-spatial-transcriptomics 对以下数据进行端到端、尽量少人工介入的注释：

input_root=/path/to/input
project_root=/path/to/output
species=物种
tissue=组织
stage_or_condition=阶段或处理
platform=平台
observation_unit=cell/cellbin/spot
primary_questions=主要生物学问题
priority_lineages=需要严格验证的稀有或关键谱系
runtime=可用的 R/Python 环境与调度资源

请自主发现输入、选择聚类强度、投递并监控任务、修复失败、维护状态并持续迭代。不要复制示例参数或标签。初次大类解释前先冻结 label-blind 的全候选正反证据、winner/runner-up 与矛盾；论文 marker 只能事后解释，不能缩窄候选集。Broad freeze 后对 analysis set 只做一次 contract-bound GSE233801 Broad-only 映射：原本未标注且在 Atlas 能力矩阵内的中高置信、非 OOD 观测可回填大类；已定义大类只做一致性/OOD 比较，不能被 Atlas 直接覆盖。完成羊卵巢三终点后，对每个大类一次只开启一个专项复核，关闭后才进入下一类。最终只发布一套注释：大类至少中等置信度、亚群仅高置信度。
```

如果输入是 BANKSY 参数网格，明确要求 Agent 自主比较候选结果，不要预先指定某个 resolution/k 值：

```text
输入目录包含多个 BANKSY 参数结果。请结合簇规模、邻近参数稳定性、DEG 可解释性、UMAP 和空间组织学形态选择聚类，而不是按文件名或固定参数选择。
```

如果输入属于同一批 StereoPy `cellbin_PPed` 转换 RDS，建议在首次消息中再加入：

```text
这些 Seurat RDS 是由同一批 StereoPy cellbin_PPed H5AD 转换得到的原始计数载体。
请不要沿用转换对象内导入的 StereoPy PCA/UMAP，也不要把转换 RDS 误判为已完成 SCT。
请按 Skill 的 same-batch Seurat 规范，从 Spatial counts 重新执行统一 SCT/PCA/邻接图前处理，
保存 analysis-set membership、输入/分析集哈希和完整 preprocessing manifest；
最终大类分辨率、各大类/定向重聚类分辨率及生物学标签仍应针对当前样本自适应判断。
```

若希望严格复用已验证成功样本的**流程策略**（不复制其参数答案），再加入：

```text
请启用 strategy_preset=sheep_ovary_same_batch_rfirst。
复用内置脱敏羊卵巢 R-first 参考的阶段顺序、大类重聚类、直接跨簇回归、残余 QC 救回、卵母细胞安全门、状态与报告合同；
不要复制参考样本的最终 resolution、簇号到标签映射、membership、比例或亚群目录。
将 active preset/profile/Skill 哈希写入项目配置和最终审计。
```

当 `species` 为 sheep/Ovis/ovine/羊、`tissue` 为 ovary/ovarian/卵巢且发现全特征 Seurat RDS 时，Skill 会自动选择 R-first。只有进一步确认 `Spatial` 原始计数层和 StereoPy `cellbin_PPed` 转换来源后，才会自动启用下方固定前处理合同；其他羊卵巢平台仍采用 R-first 策略，但不会盲目套用该技术参数。

羊卵巢卵母细胞采用两层候选策略：通过预声明多模块起始门的完整候选 cohort 进入专门的 query-only 重聚类；严格的非 ZP/母源 marker、低体细胞反程序和紧凑空间焦点只作为识别富集簇的种子证据。若重聚类后出现明确前颗粒/颗粒或基质群，直接跨谱系回归相应大类/亚群并保留来源，不再建立中间池。

多个样本建议在首次消息里直接给出样本表和并行数：

```text
请由主 Agent 维护所有样本进度、与我沟通关键决策并做跨样本审计；
每个样本只分配一个完整流程子 Agent，不能把子 Agent 简化成 cluster 重命名或只做审计。
最多并行 N 个样本，资源不足时分 wave 执行；每个样本必须独立通过完成门和发布审计。
```

主 Agent 是唯一用户入口；子 Agent 每个负责一个样本的输入发现、初始 cluster cohort/定向 cohort、二级 return、terminal residual QC、状态和报告全流程。并行只缩短等待时间，不降低证据门或交付内容。

### 2. 允许 Agent 运行完整迭代，不要把任务简化成 cluster 重命名

标准流程是：

```text
输入与表达层审计
  -> 第一轮全组织分区
     -> 只选择稳定结构并建立 one-initial-cluster : one-cohort 精确 membership
     -> 只记录 provisional broad / mixed / unknown 和 lineage watch，不写正式标签或 QC
  -> 每个 initial cluster 从项目自身非 SCT raw counts 独立执行第二轮
     -> SCT v2/glmGamPoi -> PCA -> query-only SNN -> Leiden 完整网格
     -> provisional 与历史标签不可见时，逐二级 subcluster 扫描完整 broad/fine/state/exploratory 目录
     -> 高纯度亚簇整体 parent/cross-lineage/missing-broad return
     -> 只能定义大类时保留 Broad-only；亚型证据不足时 fine 留空
     -> 只有具有互斥直接身份组件的真正 mixed 亚簇，或 specific-in-generic 组件，才启动局部拆分
     -> 优先检查相邻 Leiden 分区，仍不可分才使用局部 observation 组件和 exact remainder
     -> 未进入局部 subset 的成员重新评估 parent/unresolved，不能自动变为 QC
  -> 合并全部 cohort 为 analysis set 的互斥精确覆盖，并首次冻结正式 broad
  -> 在已冻结 broad parent 内物化高置信 fine；state 使用独立列
  -> 对完整 analysis set 做一次固定 Atlas Broad-only 映射
     -> 未标注成员：校准后中高置信、非 OOD、scope/ontology 兼容才直接回填 broad
     -> 已有大类且一致：关闭；低置信差异：记录
     -> 已有大类存在群体性可信差异或 coherent OOD：完整 source subcluster/cohort 复核一次
  -> 完成羊卵巢三终点后，严格串行地执行逐大类全样本专项复核
     -> 同一时刻只允许一个 active cell type，Agent 先输出“现在开始对 <cell type> 进行专项复核。”
     -> 当前 broad 形成一个唯一证据包和一个决策
        -> precision：检查该类型现有成员是否过召回
        -> recall：在全组织其余细胞中检查是否欠召回
        -> molecular：检查多基因、DEG/pseudobulk、竞争谱系与替代解释
        -> spatial/literature：检查全切片空间一致性和文献边界
        -> source subcluster / 空间组件 / group watch 仅作内部证据与 patch 边界
     -> 每个类型只做一次全样本专项复核；retain、absence 或精确 patch 均原子关闭该类型，不因后续数量/签名变化重新排队
     -> 后续其他类型复核若发现明确属于已关闭类型的细胞，可按证据包约束的精确 cell ID 写入或移出，并记录 post-closure delta，但不重开整套专项复核
     -> 零 census 的直接两家族/三基因 identity-core 信号即使空间碎片化，也必须进入 bounded source-subcluster 复核
     -> 当前类型关闭后才进入下一个，禁止一次产生或决策多个大类的正式复核包
     -> 默认复用稳定分区和一次性 raw-count/坐标/稀疏 count 缓存；每轮只为 active broad 及其证据竞争者算 pseudobulk
     -> retain/absence 的零变化决策不重写 membership、不追加 transform chain
  -> 缺失谱系、未建模程序、空间合理性与 residual QC 最终闭环
  -> 构建单一最终注释（中等及以上大类；仅高置信亚群）
  -> 完成门审计
  -> 主 Agent 注释质量审批（仅此时；对照已验证且脱敏的羊卵巢 R-first 流程质量，不要求标签一致）
  -> 生成确认前轻量 HTML（注释支持原因 + 高区分度大类空间图 + 大类典型 marker dotplot）
  -> 用户确认冻结注释
  -> 最终 DEG、双层树点图、空间图和 HTML 报告
```

核心原则：

- 大类先于亚群；空间数据不为了“树更深”而强制细分。
- 聚类分辨率在每个 initial-cluster cohort 或触发的临时定向 cohort 中从完整合同网格重新选择，不能照搬全局参数或示例参数。先保护稳定谱系/真实亚群，再避免状态性或技术性碎片；复杂度只在证据基本等价时作为 tie-breaker。
- HVG/BANKSY 特征可以用于聚类，但最终 marker、anti-marker、开放式谱系发现和 Oocyte 抗污染判定应回到全基因表达对象。
- 一个基因、一个参考标签或空间邻近都不能单独决定身份。
- ECM-rich、contractile、cortical、ambient、low-RNA 等是状态标签，不应替代生物学大类。
- 已关闭 cohort 有不可变 membership 和来源链；跨谱系回归直接写入目标标签，不创建新的长期池。

若有同项目或生物学匹配的单细胞数据，应先建立“来源标签 -> 空间候选大类”的可审计 crosswalk。证据优先级为：当前 query 全基因锚点与形态 > 匹配单细胞参考 > 同物种同组织公共 atlas > 跨阶段/跨物种参考。计数级参考只在所有大类/定向判断结束并冻结 QC 后映射一次全 analysis set；其直接写回权限仍只属于 QC 子集。

只有点图时，可用于完善 marker/anti-marker，不能声称完成细胞级映射。具备 count-level 参考时，默认使用预计算固定特征变换、低维参考表示和 ANN 索引；禁止默认构造 query×reference 稠密距离或每个样本重做联合整合。映射 ceiling 为大类，并显式输出低 margin、混合近邻和 OOD；未知类型不能被强制分到最近的已知类别。

### 2.1 不要混淆论文分类、计算 cohort 和最终细胞类型

这是获得稳定自动注释结果的关键：

- **论文分类目录**是候选谱系检查表，用于提醒 Agent 审查可能遗漏的细胞类型；不是要求当前样本必须补齐的答案表。
- **计算 cohort**是一次明确问题所需的不可变 membership，例如某个初始大类的重聚类集合或局部混合群的定向重聚类集合；它只记录计算边界、来源和哈希，不能直接成为最终细胞类型。所有新项目都不建立持久化生物学池。
- **最终发布大类**必须由当前样本的全基因 marker、anti-marker、稳定性与空间形态独立通过证据门。论文中存在但当前样本不支持的类型应记录 negative audit，而不是降低阈值强行创建。

以羊卵巢为例，近年整卵巢单细胞研究报告的大类数量并不一致：成年发情湖羊数据识别 5 类体细胞，西藏羊 111,548 细胞图谱报告 7 类，跨五个发育时间点的图谱报告 9 类。这种差异反映阶段、取样、消化和分辨率，而不是哪篇文章可作为固定 label map。

羊卵巢 profile 的候选目录来自分层证据，而不是单篇论文的标签照搬：

| 证据层 | 近年研究与用途 | 对自动注释的约束 |
|---|---|---|
| 羊整卵巢主参考 | [成年发情湖羊/GSE233801（J Anim Sci Biotechnol, 2023）](https://pubmed.ncbi.nlm.nih.gov/37964337/)、[西藏羊整卵巢图谱（Mol Biol Evol, 2024）](https://pmc.ncbi.nlm.nih.gov/articles/PMC10980521/)、[五个发育时间点羊卵巢图谱（iScience, 2025）](https://pubmed.ncbi.nlm.nih.gov/40641558/) | 无可用配对 count-level 参考时，唯一正式 Atlas 为 contract-bound 固定 `sheep_ovary_GSE233801_split_wall_v2`。它从原始 GSE233801 `res0.4` cluster 重新建立 Granulosa、Immune、Stromal/mesenchymal、Endothelial、Pericyte/mural 和 Smooth muscle 原型；Epithelial/mesothelial 与 Theca 仅有 challenge authority，Oocyte、Luteal 无映射权限。所有直接救回仍要求当前 query 独立的类别级校准。旧合并 v1 只允许断点续跑。完整对应见 [GSE233801 crosswalk](docs/GSE233801_to_skill_broad_taxonomy_crosswalk_v2.md)。 |
| 跨物种/多组织验证 | [九物种卵巢图谱（J Anim Sci Biotechnol, 2026）](https://pubmed.ncbi.nlm.nih.gov/41975518/)、羊–人 15 组织生殖与中枢图谱（Advanced Science, 2026, DOI 10.1002/advs.202517633）、人鼠卵巢衰老比较（Science, 2025, DOI 10.1126/science.adx0659） | 用于检查跨物种保守的大类边界、羊基因符号及 glia/平滑肌/壁细胞等候选；theca、pericyte、epithelial 细分具有物种差异，不能直接把参考标签写回 query。 |
| 成人卵巢方法学边界 | 成人卵巢单细胞专家综述（AJOG, 2025, DOI 10.1016/j.ajog.2024.05.046） | 强调取样/过滤造成的谱系缺失、表面上皮难捕获、单 marker 不可靠和“无限亚型”风险，支持以可靠浅层大类为终点。 |
| 谱系专项证据 | [羊巨噬细胞–颗粒细胞互作（FASEB J, 2026）](https://pubmed.ncbi.nlm.nih.gov/41801067/)、[人卵泡 theca–stroma 连续轨迹](https://pubmed.ncbi.nlm.nih.gov/36599970/)、[人卵巢空间图谱](https://pubmed.ncbi.nlm.nih.gov/38578993/) | 用于解决 macrophage、theca/stroma、血管/壁细胞和空间界面等竞争假设；专项论文不能越过当前样本的 marker、anti-marker 与空间门。 |

发布命名遵循“最浅且足够”的原则：`Stromal/mesenchymal` 是允许的诚实大类；只有 `CYP17A1/CYP11A1/STAR/HSD3B1` 甾体生成核心与 `INSL3/ANPEP/NR5A1/FDX1/FDXR/POR/CYB5A` 雄激素支持共同成立时才单列 `Theca`，`LHCGR` 仅作支持。Theca 必须先在完整分子候选空间中发现，卵泡 ROI 与距离只做事后解剖复核，不能用于缩小或扩张候选。`Luteal` 同时要求同一生物学亚簇中 `STAR/CYP11A1/HSD3B1` 类甾体生成核心与独立 corpus-luteum identity（候选如 `OXT/PTGFR/PARM1/LDLR/PRLR`）的直接多基因联合覆盖、直接黄体 discriminator、当前 query 的正向 DEG、阶段兼容性和实性黄体样结构；一般 steroidogenic 或全组织普遍表达的 OXT/PTGFR 只是功能/背景程序，不作为 broad 名称，且 Luteal 禁止使用 dominant-identity 整簇捷径。`Endothelial`、`Pericyte/mural` 与 `Smooth muscle` 是三个独立竞争大类：前者要求内皮连接骨架与独立血管支持，周细胞要求 mural identity backbone 与独立支持，平滑肌要求 `MYH11/CNN1/ACTG2/SMTN/LMOD1` 成熟收缩核心；`ACTA2/TAGLN/MYL9` 或血管邻近本身不能定类。高置信 `Lymphatic endothelial` 是 Endothelial 的内部 fine identity，公开唯一 `final_cell_type` 可直接显示该名称；旧 `Vascular-associated` 不再允许发布。`Oocyte` 以完整 canonical cluster 判定；常规二级扫描为零但全对象多模块起始候选非零时，允许一次标签不可见的 query-only targeted cohort，并纳入通过簇全部成员，仅排除客观输入 QC 或直接多家族体细胞硬矛盾。颗粒细胞仅在完整、稳定且有文献支持的功能程序通过时使用浅层亚型，否则保留 `Granulosa` 大类。

推荐把 `Granulosa`、`Stromal/mesenchymal`、`Endothelial`、`Pericyte/mural`、`Smooth muscle`、`Immune`、`Epithelial/mesothelial` 和严格门控的 `Oocyte` 作为候选审查骨架，但审查不止这些类型。Skill 内置机器可读的羊卵巢候选谱系目录，逐样本覆盖卵泡/生殖系、甾体生成、基质–间充质–收缩/壁细胞、血管/淋巴、免疫、上皮/间皮及神经胶质/神经内分泌边界；`Theca`、`Luteal` 和 `Glial/Schwann-like` 等只有在各自多通道身份程序通过后才作为大类发布，`Mesenchymal progenitor-like` 与 `Neuroendocrine-like` 默认保持 exploratory。目录不是答案表，也不穷尽当前样本：每个候选可以得到阴性结论，目录外的相干多基因程序则必须新增审查。不要使用 `Vascular-associated`、`Theca/follicular wall`、`Stromal/perivascular` 或 generic `steroidogenic` 作为方便的最终兜底大类。公开报告、空间图、census、DEG 与主 dotplot 仅使用单列 `final_cell_type`；broad/fine 只保留为内部审计 provenance。

Agent 通过 `init_open_world_lineage_audit.py` 从当前 cluster ledger 和目录生成完整审查骨架，填入各边界的全基因 DEG、anti-marker、相邻分辨率稳定性与空间证据后，再由 `validate_open_world_lineage_audit.py` 校验。目录、审查源文件和 biological profile 均参与哈希绑定；只审查少数示例、漏掉阴性结论或在完成门后修改目录都会使发布失效。

Skill 的默认架构是：第一轮全组织聚类只生成 provisional cohort 计划；每个初始 cluster 都从项目自身 raw counts 进入第二轮 SCT/PCA/SNN/Leiden，并在标签不可见状态下扫描完整候选目录；高纯度二级亚簇整体返回，真正 mixed 的亚簇才进行局部 observation 拆分；所有 cohort 合并后冻结 broad，再做 Atlas、缺失谱系和 residual QC 复核。若在 Oocyte 定向重聚类中识别出可信前颗粒细胞，直接跨谱系写回 `Granulosa`，不创建中转 cohort，也不自动再次聚类。

Skill 内置一份[脱敏羊卵巢 R-first forward-test 参考](annotate-spatial-transcriptomics/references/profiles/sheep_ovary_rfirst_case_reference.md)，完整记录从原始转换 RDS、开放式谱系发现、大类/定向 cohort、跨谱系回归、Oocyte 抗污染门到残余 QC 救回和最终 DEG/报告的策略，但不包含样本 ID、私有路径、观测 ID、cluster 答案或历史映射。它是策略回归参考，不是 label map。

### 3. 把大类重聚类、跨谱系回归和残余 QC 处理作为流程核心

不要用“未注释比例低”作为完成标准，也不要把所有弱信号细胞直接过滤掉。

- **第一轮只建 cohort。** 选择保留主要组织结构、不过度技术碎裂的全组织分辨率；仅记录 provisional broad/mixed/unknown 与 lineage watch，不产生正式 broad/fine/QC。
- **第二轮负责注释。** 每个初始 cluster 独立运行完整 `0.1,0.2,0.3,0.4,0.6` 网格并扫描全部候选。能可靠定义亚型才写 fine；亚型证据不足时保留 supported broad，不能为了目录完整而强行命名。
- **跨谱系与局部混合。** 一个清晰完整谱系直接 parent/cross-lineage return；多个候选可见或共享背景信号不等于 mixed。清晰 specific parent 只有在另一谱系也具备独立亚簇级多基因、DEG/pseudobulk、候选特异 direct discriminator，且双方形成实质性互斥身份组件时才重开。提议的拆分还必须在至少一个相邻高分辨率 Leiden 分区中复现为富集且可分离的组件；随后 observation 组件只负责界定 selected mixed subcluster 内的精确成员。通用 Stromal 背景中的少数 specific 组件也遵循同样证据链。共享 steroidogenesis/ECM/contractile 背景或高分辨率不复现的散在尾部只记 watch；强信号完全共表达且无法消歧时保留 unresolved。未进入局部子集的成员重新评估 parent/remainder，不自动进入 QC。
- **语义修复不重跑稳定分区。** 输入、cohort membership、grid、seed 与聚类脚本哈希均未变化时，复用 controller 生成的 derived partition，只重算 scorer、边界和写回。任何修复 proposal 必须用 `apply_cell_id_membership_patch.py` 按唯一 `cell_id` 连接，禁止按行序赋值。
- **全细胞 Atlas 与逐细胞类型复核。** broad 冻结后才运行全细胞 broad mapping；Atlas 只能直接救回 unlabeled、非 OOD 的中高置信成员，不能静默覆盖已有 broad 或生成 fine。随后以“一个 broad 细胞类型”为一级单位逐个复核：在同一次复核中检查现有成员的过召回、全组织其余细胞中的欠召回，以及整个切片上的空间分布。每个已注释 broad 必须分别给出成员纯度、全组织召回、分子身份和空间合理性结论，缺一项不能关闭。source subcluster、严格逐细胞空间组件和亚簇级信号缺口 watch 只作为内部证据与 patch 边界，不拆成大量用户任务；group watch 只能触发原始 counts 定位，不能整簇写回。
- **上下文只控制能否评估。** 阶段、处理或解剖背景不能充当身份分数。`not_evaluable` 候选不得形成完整性阳性、Atlas 救回、fine 或最终标签，也不能因为其阳性程序计数而阻断一个合理的零 census；若最终 membership 已含该标签则必须失败并回到来源复核。

Release-critical 的控制器默认阈值统一来自 `annotate-spatial-transcriptomics/references/controller_thresholds_v2_2.json`，包括 direct/local scoring 权重、whole-subcluster 与 local-subset 写回、resolution selector 权重及最终 residual-QC 阈值。初始化项目、构建 annotation contract、Python adjudicator 和 R scorer 均读取或绑定该注册表；项目覆盖必须显式写入 contract。Atlas classwise 校准目标和 StereoPy 技术前处理参数仍分别属于 workflow profile，不能与控制器写回阈值混为一张表。

同批次羊卵巢固定网格只在 active workflow profile、显式 strategy preset、对象/层/marker/坐标审计和 StereoPy conversion provenance 共同通过后生效。prelabel freeze、cohort outcome、direct return、Atlas concordance/review 和最终 broad/fine support 都必须通过内容与哈希合同。

Atlas 阈值必须通过当前项目的 query-depth-matched held-out anchors 校准。默认目标精度为 `moderate-or-higher >= 0.90`、`high >= 0.95`，它们是**校准目标精度**，不是通用 raw score 截断值。输出应满足：

```text
moderate_or_higher_n = high_n + moderate_only_n
```

经过独立证据复核的 `high` 和 `moderate_only` 均可回归**大类**；`low_reject` 保留为 QC。Atlas 救回的 broad-only 细胞必须设置 `fine_anchor_eligible=false`，不能反向参与精细 marker 或亚群锚点发现。

这里的 held-out anchors 必须来自当前 query 的独立高置信锚点，并与待救回 membership 不重叠。把同一个外部 atlas 随机拆成 train/held-out 只是在测 atlas 自分类，不能校准 query 救回；旧版合并 `medium_high` 阈值只允许做诊断，不能写回。无配对 count-level 羊单细胞对象时，羊卵巢默认以 GSE233801 为公共 atlas 主通道；只有配对 marker dotplot 时仍不能执行细胞级转移。

### 4. 对稀有、易污染谱系使用更严格的上下文门

卵母细胞、浆细胞、淋巴内皮等稀有类型不能因为单一强 marker 就大范围定义。应同时要求：

- 多基因正向程序；
- 排除邻近高丰度谱系的 anti-program；
- 在全基因表达层验证，而非仅 HVG；
- 与组织学相符的局灶空间结构；
- 在候选定向 cohort 重聚类后仍形成稳定、纯净的小群；
- 明确区分 cellbin/spot 数量与真实生物学细胞数量。

未通过上下文门的候选若具有明确体细胞程序，则直接回归相应大类；没有明确身份或混合严重者进入 `qc_holdout`，而不是保留一个看似“稀有”的宽松标签。

### 5. 完成注释后先由主 Agent 审批，再进入用户确认

Agent 应自主完成常规分辨率比较、大类/定向 cohort 重聚类、跨谱系回归、残余 QC 救回、作业投递、监控、日志排错、重投、最终回写和状态审计。只有这些工作全部完成且 completion gate 通过后，主 Agent 才审批注释质量；不能在仅完成大类注释时提前审批。审批只复核大类合理性、marker/anti-marker 与空间形态支持、易混淆群安全性及整体是否达到内置脱敏 R-first 成功流程的质量水平，不重复机械完成门，允许“通过但有备注”。审批通过后再生成轻量化审阅 HTML 供用户确认；用户确认后才生成耗时的最终 DEG、完整图和发布报告。

如果冻结注释之后又修改 ledger、route 或 completion gate，主 Agent 审批和用户确认都会自动失效，必须重新审批/确认。这能避免在注释仍变化时反复重画全套报告。

### 6. 用结果合同验收，而不是只看一张 UMAP

最终交付至少应包括：

- 全部样本都输出大类 DEG；只有真实高置信亚群存在时才输出独立亚群 DEG；
- 全部样本都输出大类树状 marker dotplot；只有真实高置信亚群存在时才输出亚群树状 dotplot，零亚群是合法结果；
- 两个层级分别包含 canonical marker 与当前数据特异 marker；
- dotplot 同时输出 PNG、PDF 和源 TSV；点大小和颜色按基因内部归一化展示，但源表保留绝对检测率和平均表达量；
- 大类/亚群 UMAP、全空间图、逐节点高亮网格；
- 按支持细胞类型分组的 marker 空间表达图；
- 只发布一套最终注释：大类至少达到中等置信度，亚群必须达到高置信度；
- 最终大类 DEG 和 marker dotplot 必须纳入所有 direct return 与 Atlas broad rescue 后正式属于该大类的观测；亚群 DEG 仅使用具有真实高置信 fine label 且允许作为精细证据的细胞，禁止把 broad-only 救回强行分到亚群；
- 可展开注释树、路线/阈值/结局面板、中文详细流程时间线和原始状态记录；
- cell-level ledger、cluster decision ledger、recluster-cohort/direct-return/run/route registry、session info、manifest 和 checksums；
- completion gate 与 release audit 均通过。
- 生物学大类统计与解剖界面/QC/技术/待审状态统计分开；后者不能进入生物学大类 DEG、marker dotplot 或注释树。

证据不是“有文件即可”。每个 cohort outcome 必须通过 `schemas/cohort_outcome.schema.json`，记录 query membership/hash、完整候选网格、逐 resolution cluster membership 与证据索引、相邻 ARI/迁移、full-feature marker/anti-marker、空间形态、来源/QC 组成、选择理由、被拒候选及每个亚簇的互斥结局。每个 direct return 和最终 broad/fine support 也必须通过各自 schema。空 JSON、只有 status 的 JSON、空 TSV 或过期哈希都不能通过。

`audit_annotation_membership_partition.py` 强制逐 cell 闭合：每个 initial cluster 必须且只能进入一个二级 cohort，全部 cohort query 互斥并精确覆盖 analysis set；二级 whole-subcluster return、局部 supported subset 与 local remainder 结局互斥完备；Atlas query 精确等于 broad freeze 后的 analysis set，terminal residual QC 与 ledger 写回一致；最终注释唯一覆盖 analysis set。

## 运行环境说明

一键安装只配置 Skill 本身，不会擅自修改分析环境。真实项目可能需要 Seurat、sctransform、Scanpy、anndata、BANKSY、spacexr/RCTD 等；Skill 会先检查已有环境，并优先使用项目指定的 R/Python 环境和计算集群。内存密集型表达矩阵操作应投递到调度节点，本地仅用于发现、审计和小型汇总。

### 同批次 StereoPy cellbin 转换 RDS

由 StereoPy `cellbin_PPed` 批量转换得到的 Seurat RDS 只是原始计数输入载体，并不等同于已经完成 SCT 前处理。转换步骤通常只把 `exp_matrix@raw` 写入 Seurat `Spatial` counts，保留坐标，并可能复制 StereoPy PCA/UMAP 用于溯源；这些 reduction 不能作为新的 R 聚类图。

详细规则见 [Seurat cellbin 前处理合同](annotate-spatial-transcriptomics/references/seurat-cellbin-preprocessing.md)。全组织 SCT+BANKSY 与第二轮 cohort SNN 是两套不同图构建过程，参数不得混用。

同一生产批次的**全组织 SCT+BANKSY** 固定参数：

| 阶段 | 固定参数 |
|---|---|
| 分析入口 | `nCount_Spatial >= 100 AND nFeature_Spatial >= 75` |
| 高计数处理 | 高于样本 99.9% 分位数只加 review flag，不在入口硬删除 |
| 线粒体/双细胞 | 羊基因符号无法可靠识别时不做线粒体硬过滤；入口不做 doublet 硬删除 |
| SCTransform | `Spatial -> SCT`，`vst.flavor="v2"`，`method="glmGamPoi"` |
| SCT规模 | 4,000 variable features，`ncells=min(50000, n)` |
| SCT内存策略 | `conserve.memory=TRUE`，`return.only.var.genes=TRUE` |
| BANKSY输入 | `SCT/scale.data` 的 4,000 个 Pearson-residual HVG |
| BANKSY图 | `M=0`，`k_geom=30`，`lambda=0.2`，30 PCs，Leiden `k=50` |
| 全组织候选网格 | `0.2,0.4,0.6,0.8` |
| UMAP | 30 neighbors，`min.dist=0.3`，`spread=1`，300 epochs |

每个 initial-cluster **第二轮 cohort** 固定边界：

| 阶段 | 合同 |
|---|---|
| 表达入口 | 项目自身非 SCT raw counts；禁止 SCT corrected counts 再做 SCTransform |
| 图构建 | `SCT v2/glmGamPoi -> PCA -> query-only SNN -> Leiden` |
| PCs / k | 依据当前 cohort 规模、谱复杂度和相邻分辨率稳定性自适应；不复用 BANKSY 参数 |
| Leiden候选网格 | `0.1,0.2,0.3,0.4,0.6`，`algorithm=4` |
| Leiden并行 | 默认 `resolution-workers=1`；大型 Seurat carrier 禁止按 resolution fork 复制对象。64 CPU 用于矩阵内核或并行独立 cohort，而不是同一对象的 5 个分辨率副本 |

全组织统一前处理可直接运行：

```bash
Rscript annotate-spatial-transcriptomics/scripts/run_seurat_sct_preprocess.R \
  --rds /path/to/sample.seurat.rds \
  --out /path/to/project/reclustering/seurat_sct \
  --sample SAMPLE_ID \
  --assay Spatial \
  --resolutions 0.2,0.4,0.6,0.8 \
  --future-globals-max-gb 100
```

脚本会输出：

- 保留全部输入观测的 `analysis_scope.tsv.gz`，明确区分 `analysis_set` 与 `excluded_initial_qc`；
- 每个候选分辨率的 cluster membership 和 cluster count；
- 新生成的 SCT/PCA/UMAP Seurat RDS；
- UMAP、空间坐标、输入 SHA256、分析集 SHA256、全部前处理参数和 `sessionInfo`；
- `RUN_COMPLETE.tsv`。脚本不会替 Agent 自动选择最终分辨率或写入生物学标签。

每个 initial-cluster cohort 或触发的定向 cohort（底层复用同一标准 runner）重聚类使用：

```bash
Rscript annotate-spatial-transcriptomics/scripts/run_seurat_cohort_recluster.R \
  --rds /path/to/full_feature.seurat.rds \
  --membership /path/to/frozen_query_anchor_membership.tsv \
  --out /path/to/project/reclustering/COHORT_ID \
  --cell-id-col cell_id \
  --assay Spatial \
  --resolutions 0.1,0.2,0.3,0.4,0.6 \
  --resolution-contract sheep_ovary \
  --resolution-workers 1
```

membership 至少需要唯一 `cell_id`。锚点模式还需要 `query_or_anchor` 和 `anchor_label`；anchors 与 query 共同参与 SCT/PCA，但邻接图、Leiden、UMAP 和 DEG 必须只在 query 中计算。cohort DEG 使用全基因 `Spatial` LogNormalize 数据，而不是只在 SCT variable features 中寻找 marker。

正式控制器会在重聚类前生成 cohort 资源计划：≥20,000 observations 时强制单 resolution worker，并建议 64 CPU；若调用者声明的内存低于经验峰值安全线，任务会在读取和复制大对象前阻断。

需要固定的是**同批次技术前处理**，不能固定的是**生物学决策**：

- 最终大类分辨率必须结合 DEG、marker/anti-marker、相邻分辨率稳定性、UMAP 和空间形态选择；
- 每个 initial-cluster cohort 或触发的定向 cohort 重新决定 PCs、k 和最终分辨率；羊卵巢正式候选分辨率始终为 `0.1,0.2,0.3,0.4,0.6`，小型卵母细胞/免疫 cohort 可降低 PC 和 k，但不能改用低于 0.1 的网格；
- `<100` observations 只触发小簇审查，不能自动并入最近 PCA 簇；
- 大类、亚群、合并关系和置信度不能从示例样本复制。

如果 `glmGamPoi` 或 SHA256 依赖不可用，SCT 路线会直接失败并要求修复环境，不会静默切换模型或写出不可验证 manifest。固定技术参数的任何覆盖都必须显式传入 batch exception 和理由。只有 preprocessing manifest 与当前输入/分析集哈希匹配时，已有 SCT 计算才允许复用。

调度资源必须与真实并行度一致。固定资源档为：状态/preflight/hash 1–4 CPU，逐大类 marker/DEG/pseudobulk 8 CPU（最多 16），SCT/PCA/SNN/Leiden 或大型 RDS 物化 64 CPU。`leidenbase` 的单个 Leiden 优化本身是单线程，因此 64 CPU 重任务仍保持单对象 `resolution-workers=1`；额外核心用于矩阵内核或并行独立 cohort，不能把同一大型 Seurat 对象复制 64 份。当前 Seurat 的直接 `FindAllMarkers` 也是逐簇串行；Skill 会记录实际 worker/backend/parallel unit，并用 CPU time / wall time 审计资源是否被浪费。

## 仓库结构

```text
annotate-spatial-transcriptomics/
  SKILL.md                 # Agent 入口与强制工作流
  agents/openai.yaml       # Skill 元数据
  references/              # 聚类、路由、状态、报告、测试规范
  references/profiles/     # 组织/物种背景 profile（不是标签映射）
  assets/                  # marker、锚点与多通道救回配置模板
  schemas/                 # cohort、direct return 与逐标签支持证据 schema
  scripts/                 # 状态控制、聚类、校准、作图与报告脚本
install.sh                 # 一键安装
scripts/verify_install.py  # 安装完整性检查
scripts/validate_repo.py   # 发布前脱敏、语法与结构审计
```

仓库不包含任何测试样本的表达矩阵、注释结果、状态文件或服务器私有路径。组织 profile 仅提供背景约束与验证策略，禁止当作 cluster-to-label 映射。

## 更新、验证与开发

更新到最新版：

```bash
curl -fsSL https://raw.githubusercontent.com/CCCU-IMU/annotate-spatial-transcriptomics/main/install.sh | bash
```

验证已安装版本：

```bash
python "${CODEX_HOME:-$HOME/.codex}/skills/annotate-spatial-transcriptomics/scripts/check_runtime.py"
python scripts/verify_install.py "${CODEX_HOME:-$HOME/.codex}/skills/annotate-spatial-transcriptomics"
```

发布前验证：

```bash
python -m pip install -r requirements-ci.txt
python scripts/validate_repo.py
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -v
bash -n install.sh
```

GitHub Actions 的 PR 验证与 Release 打包均使用 Python 3.11 和仓库内的 `requirements-ci.txt` 建立干净验证环境。该依赖文件包含报告元数据与发布合同测试实际导入的包；发布测试不得因缺少可选依赖而静默跳过关键合同。

对 Skill 的重要更新必须做 leakage-safe forward test：只给原始输入与 Skill，不给预期标签；评价它能否建立稳定的第一轮 cohort、完成所有第二轮重聚类与 missing-broad 重建、局部拆解真正 mixed 的二级亚簇、执行全细胞 Atlas/OOD 复核、闭合残余 QC 并生成完整报告。历史结果只可在新 membership 冻结后作为外部生物学验收，不能作为运行时输入。

## 版本

当前版本：`2.5.0`（稳定版），项目 framework schema 仍为 `2.0.0`，内部 controller/artifact protocol 保持 `2.2.0` 以兼容已冻结项目和可复用分区。v2.5 以第二轮重聚类为正式注释主体：第一轮只生成精确 cohort；每个二级 subcluster 扫描完整 broad/fine/state/exploratory 目录；只有竞争谱系共存的二级 mixed subcluster 才局部拆分。全部 cohort 合并后冻结 broad，再执行 calibrated Atlas、羊卵巢三终点、逐大类全样本 precision/recall、fine/state parent lock 和最终 QC。该稳定版同时发布独立 Endothelial、Pericyte/mural、Smooth muscle taxonomy，修复 Atlas 90%/95% 校准、ROI 证据链、逐大类复核来源绑定、局部拆分工作量审计和高 unresolved 有界停止；任何成功样本的数量、簇号或空间 ROI 都不是运行时先验。

## 许可

MIT License。详见 [LICENSE](LICENSE)。
