# Annotate Spatial Transcriptomics Skill

[![Validate](https://github.com/CCCU-IMU/annotate-spatial-transcriptomics/actions/workflows/validate.yml/badge.svg)](https://github.com/CCCU-IMU/annotate-spatial-transcriptomics/actions/workflows/validate.yml)

面向 Codex Agent 的空间转录组/单细胞转录组迭代式注释 Skill。它不会复制旧项目答案，而是自主完成：label-blind 全候选证据冻结、聚类选择、大类 cohort 纯度审查、直接回归、定向重聚类、一次全细胞 Atlas 大类一致性/OOD 审计、残余 QC 回填、状态追踪、失败恢复、用户确认和可审计报告。

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
  | bash -s -- --ref v1.10.0
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

Seurat 用户要特别注意双对象合同：SCT 对象只负责聚类；如果原始转换对象的 `Spatial@data` 与 `counts` 完全相同，Agent 必须在计算节点另建带 manifest 的全基因 LogNormalize 验证对象，供 DEG/marker 使用，不能改写 SCT 对象。候选分辨率比较同时保留全观测 ARI/AMI；`n<100` 只触发小簇复核，并且只从额外的宏观评分中暂时隔离，绝不能据此删除细胞或跳过 DEG、空间和稀有谱系审查。

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

请自主发现输入、选择聚类强度、投递并监控任务、修复失败、维护状态并持续迭代。不要复制示例参数或标签。初次大类解释前先冻结 label-blind 的全候选正反证据、winner/runner-up 与矛盾；论文 marker 只能事后解释，不能缩窄候选集。所有大类和定向判断结束并冻结 QC 后，对 analysis set 做一次 Atlas Broad-only 映射：QC 通过多通道门后直接回填大类；已定义大类只做一致性/OOD 比较，群体性差异才复核，不能被 Atlas 直接覆盖。最终只发布一套注释：大类至少中等置信度、亚群仅高置信度。
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

主 Agent 是唯一用户入口；子 Agent 每个负责一个样本的输入发现、broad/targeted cohort、direct return、terminal residual QC、状态和报告全流程。并行只缩短等待时间，不降低证据门或交付内容。

### 2. 允许 Agent 运行完整迭代，不要把任务简化成 cluster 重命名

标准流程是：

```text
输入与表达层审计
  -> 自适应选择大类聚类
  -> 在看到论文标签前冻结全候选大类证据矩阵
     -> 正 DEG、anti-DEG、程序分数、winner/runner-up、margin、矛盾与技术状态
  -> 对选定分辨率的每个簇执行开放式谱系审查并初步注释大类
     -> 中高置信：直接写入初步大类
     -> 低信息、无特征或无法解释的严重混合：QC holdout
  -> 每个初步大类独立投递一次 broad_class_recluster cohort
     -> 按 active workflow contract 完整运行候选网格
     -> 在完整网格中选择当前生物学问题综合证据最优的分辨率
     -> 若没有稳定亚群或混杂，合法结束为 homogeneous_parent_confirmed 并全部回归父大类
  -> 亚簇判断
     -> 高置信亚群：写入大类 + 亚群
     -> 只能确定父类：直接回归父大类 Broad-only
     -> 明确属于其他谱系：跨谱系直接回归目标大类/亚群并保留来源
     -> 可解释竞争混合：临时定向重聚类；必要时低优先级 RCTD
     -> 无法解释：QC holdout
  -> 冻结最终残余 QC，使用固定 Atlas 表示/索引对 analysis set 做一次 Broad-only 映射
     -> 原标签为 QC：多通道中高置信且非 OOD，直接回填大类；否则保留 QC
     -> 已有大类且一致：关闭；低置信差异：记录
     -> 已有大类存在群体性可信差异或 coherent OOD：完整 cluster/cohort 复核一次
  -> 构建单一最终注释（中等及以上大类；仅高置信亚群）
  -> 完成门审计
  -> 主 Agent 注释质量审批（仅此时；对照已验证且脱敏的羊卵巢 R-first 流程质量，不要求标签一致）
  -> 生成确认前轻量 HTML（注释支持原因 + 高区分度大类空间图 + 大类典型 marker dotplot）
  -> 用户确认冻结注释
  -> 最终 DEG、双层树点图、空间图和 HTML 报告
```

核心原则：

- 大类先于亚群；空间数据不为了“树更深”而强制细分。
- 聚类分辨率在每个大类或临时定向 cohort 中从完整合同网格重新选择，不能照搬全局参数或示例参数。先保护稳定谱系/真实亚群，再避免状态性或技术性碎片；复杂度只在证据基本等价时作为 tie-breaker。
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
| 羊整卵巢主参考 | [成年发情湖羊/GSE233801（J Anim Sci Biotechnol, 2023）](https://pubmed.ncbi.nlm.nih.gov/37964337/)、[西藏羊整卵巢图谱（Mol Biol Evol, 2024）](https://pmc.ncbi.nlm.nih.gov/articles/PMC10980521/)、[五个发育时间点羊卵巢图谱（iScience, 2025）](https://pubmed.ncbi.nlm.nih.gov/40641558/) | 无可用配对 count-level 参考时，GSE233801 是成年羊体细胞**残余 QC** Atlas 救回的主参考，主要覆盖颗粒、基质、血管/壁细胞候选和免疫；不同文章的 5、7、9 类结果只用于候选审查与 negative audit。 |
| 跨物种/多组织验证 | [九物种卵巢图谱（J Anim Sci Biotechnol, 2026）](https://pubmed.ncbi.nlm.nih.gov/41975518/)、羊–人 15 组织生殖与中枢图谱（Advanced Science, 2026, DOI 10.1002/advs.202517633）、人鼠卵巢衰老比较（Science, 2025, DOI 10.1126/science.adx0659） | 用于检查跨物种保守的大类边界、羊基因符号及 glia/平滑肌/壁细胞等候选；theca、pericyte、epithelial 细分具有物种差异，不能直接把参考标签写回 query。 |
| 成人卵巢方法学边界 | 成人卵巢单细胞专家综述（AJOG, 2025, DOI 10.1016/j.ajog.2024.05.046） | 强调取样/过滤造成的谱系缺失、表面上皮难捕获、单 marker 不可靠和“无限亚型”风险，支持以可靠浅层大类为终点。 |
| 谱系专项证据 | [羊巨噬细胞–颗粒细胞互作（FASEB J, 2026）](https://pubmed.ncbi.nlm.nih.gov/41801067/)、[人卵泡 theca–stroma 连续轨迹](https://pubmed.ncbi.nlm.nih.gov/36599970/)、[人卵巢空间图谱](https://pubmed.ncbi.nlm.nih.gov/38578993/) | 用于解决 macrophage、theca/stroma、血管/壁细胞和空间界面等竞争假设；专项论文不能越过当前样本的 marker、anti-marker 与空间门。 |

发布命名遵循“最浅且足够”的原则：`Stromal/mesenchymal` 是允许的诚实大类；只有形成完整 `STAR/CYP11A1/CYP17A1/HSD3B/NR5A1/LHCGR` 甾体生成或雄激素程序时才单列 `Theca`，仅有 `ALPL/PTCH1`、胶原或收缩信号的卵泡壁细胞直接回归其支持的基质/壁细胞大类；`Smooth muscle` 需要 `MYH11/MYL9/TAGLN/ACTA2/CNN1` 成熟收缩骨架。内皮和周细胞统一发布为 `Vascular-associated` 大类，`Blood endothelial`、`Lymphatic endothelial`、`Pericyte/mural` 作为可选高置信亚型；成熟平滑肌仍为独立大类。`Oocyte` 先以 strict seed/空间对象判定通过簇，再纳入该簇全部 canonical cohort 成员，不能用 seed 或 object ID 逐细胞缩小最终集合。颗粒细胞仅在完整、稳定且有文献支持的功能程序通过时使用浅层亚型，否则保留 `Granulosa` 大类。

推荐把 `Granulosa`、`Stromal/mesenchymal`、`Vascular-associated`、`Immune`、`Epithelial/mesothelial` 和严格门控的 `Oocyte` 作为候选审查骨架，但审查不止这些类型。Skill 内置机器可读的羊卵巢候选谱系目录，逐样本覆盖卵泡/生殖系、甾体生成、基质–间充质–收缩/壁细胞、血管/淋巴、免疫、上皮/间皮及神经胶质/神经内分泌边界；`Theca`、`Smooth muscle`、`Mesenchymal progenitor-like`、`Luteal steroidogenic` 和 `Neural/Schwann/glia` 等只有在各自证据门通过后才作为大类发布，周细胞和内皮谱系则保留在共享血管大类下的可选亚型。目录不是答案表，也不穷尽当前样本：每个候选可以得到阴性结论，目录外的相干多基因程序则必须新增审查。不要使用 `Theca/follicular wall` 或 `Stromal/perivascular` 作为方便的最终兜底大类。

Agent 通过 `init_open_world_lineage_audit.py` 从当前 cluster ledger 和目录生成完整审查骨架，填入各边界的全基因 DEG、anti-marker、相邻分辨率稳定性与空间证据后，再由 `validate_open_world_lineage_audit.py` 校验。目录、审查源文件和 biological profile 均参与哈希绑定；只审查少数示例、漏掉阴性结论或在完成门后修改目录都会使发布失效。

Skill 的默认架构只建立三类可追溯边界：每个可信初始大类一个 `broad_class_recluster` cohort；只有局部、可解释的混合/污染问题才建立一次性 `targeted_recluster` cohort；最终仍无足够身份信息的观测进入 terminal residual QC。羊卵巢 profile 在此架构上补充固定技术前处理、开放谱系候选、Oocyte 安全门和 GSE233801 参考策略，而不是替换控制器。若在 Oocyte 定向重聚类中识别出可信前颗粒细胞，直接跨谱系写回 `Granulosa`（必要时写入高置信亚群），保留来源和证据，不创建中转 cohort，也不自动进行二次重聚类。

Skill 内置一份[脱敏羊卵巢 R-first forward-test 参考](annotate-spatial-transcriptomics/references/profiles/sheep_ovary_rfirst_case_reference.md)，完整记录从原始转换 RDS、开放式谱系发现、大类/定向 cohort、跨谱系回归、Oocyte 抗污染门到残余 QC 救回和最终 DEG/报告的策略，但不包含样本 ID、私有路径、观测 ID、cluster 答案或历史映射。它是策略回归参考，不是 label map。

### 3. 把大类重聚类、跨谱系回归和残余 QC 处理作为流程核心

不要用“未注释比例低”作为完成标准，也不要把所有弱信号细胞直接过滤掉。

- **初始大类。** 在自适应选出的全组织分辨率上，逐簇完成开放式谱系审查；可信簇直接获得初始大类，低信息/不可解释簇直接进入 `qc_holdout`，此时不建立中转 cohort。
- **大类重聚类。** 每个可信初始大类建立一个不可变 `broad_class_recluster` cohort。羊卵巢同批次 StereoPy cellbin R-first 合同必须完整运行 `0.1,0.2,0.3,0.4,0.6`；其他平台读取自己的 active workflow profile。可定义真实亚群时写入高置信亚群；没有稳定亚群时以 `homogeneous_parent_confirmed` 成功结束并全部回归父大类。只有正式记录、哈希绑定的 underpowered skip 可跳过计算。
- **跨谱系与局部混合。** 某个大类重聚类中出现另一谱系的完整程序时，直接写回目标大类/高置信亚群，不建立中转 cohort，也不自动随目标大类再次重聚类。只有界面/污染信号确实可解释且需要拆分时才建立一次性定向 cohort；RCTD 仅为低优先级辅助：canonical `high` 且有独立证据可支持亚群，`moderate` 只支持大类，`low` 进入 `qc_holdout`。
- **全细胞 Atlas 审计与残余 QC 救回。** 冻结最终 `qc_holdout` 后，对完整 analysis set 做一次 Broad-only 映射。QC 只有在校准层级、当前 query marker/anti-marker、内部锚点/空间通道和 OOD 门均通过时直接回填大类。已定义标签仅参与一致性比较；默认至少 30 个观测且占 cluster 20% 的可信差异才复核完整 cluster/cohort。coherent OOD 小群体可作为 unknown 例外触发。

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

`audit_annotation_membership_partition.py` 强制逐 cell 闭合：initial broad 与 broad cohort query 相等；cohort 结局互斥完备；direct return 不重叠；v1.7 Atlas query 精确等于 analysis set，QC 子 membership 精确等于 terminal residual QC，accepted/rejected 划分 QC 并与 ledger 写回一致；最终注释唯一覆盖 analysis set。

## 运行环境说明

一键安装只配置 Skill 本身，不会擅自修改分析环境。真实项目可能需要 Seurat、sctransform、Scanpy、anndata、BANKSY、spacexr/RCTD 等；Skill 会先检查已有环境，并优先使用项目指定的 R/Python 环境和计算集群。内存密集型表达矩阵操作应投递到调度节点，本地仅用于发现、审计和小型汇总。

### 同批次 StereoPy cellbin 转换 RDS

由 StereoPy `cellbin_PPed` 批量转换得到的 Seurat RDS 只是原始计数输入载体，并不等同于已经完成 SCT 前处理。转换步骤通常只把 `exp_matrix@raw` 写入 Seurat `Spatial` counts，保留坐标，并可能复制 StereoPy PCA/UMAP 用于溯源；这些 reduction 不能作为新的 R 聚类图。

详细规则见 [Seurat cellbin 前处理合同](annotate-spatial-transcriptomics/references/seurat-cellbin-preprocessing.md)。同一生产批次如需和已验证 R 流程保持可比，固定以下计算参数：

| 阶段 | 固定参数 |
|---|---|
| 分析入口 | `nCount_Spatial >= 100 AND nFeature_Spatial >= 75` |
| 高计数处理 | 高于样本 99.9% 分位数只加 review flag，不在入口硬删除 |
| 线粒体/双细胞 | 羊基因符号无法可靠识别时不做线粒体硬过滤；入口不做 doublet 硬删除 |
| SCTransform | `Spatial -> SCT`，`vst.flavor="v2"`，`method="glmGamPoi"` |
| SCT规模 | 3,000 variable features，`ncells=min(50000, n)` |
| SCT内存策略 | `conserve.memory=TRUE`，`return.only.var.genes=TRUE` |
| PCA | 计算50个PC，邻接图使用前30个PC |
| 邻接图 | k=30，Annoy，50 trees，cosine |
| Leiden候选网格 | `0.1,0.2,0.3,0.4,0.6`，`algorithm=4` |
| Leiden并行 | 单个 Leiden 优化不支持多线程；用 5 个 `future` worker 并行计算 5 个候选 resolution |
| UMAP | 30 neighbors，`min.dist=0.3`，cosine |

全组织统一前处理可直接运行：

```bash
Rscript annotate-spatial-transcriptomics/scripts/run_seurat_sct_preprocess.R \
  --rds /path/to/sample.seurat.rds \
  --out /path/to/project/reclustering/seurat_sct \
  --sample SAMPLE_ID \
  --assay Spatial \
  --resolutions 0.1,0.2,0.3,0.4,0.6 \
  --resolution-workers 5 \
  --resolution-future-plan auto \
  --future-globals-max-gb 100
```

脚本会输出：

- 保留全部输入观测的 `analysis_scope.tsv.gz`，明确区分 `analysis_set` 与 `excluded_initial_qc`；
- 每个候选分辨率的 cluster membership 和 cluster count；
- 新生成的 SCT/PCA/UMAP Seurat RDS；
- UMAP、空间坐标、输入 SHA256、分析集 SHA256、全部前处理参数和 `sessionInfo`；
- `RUN_COMPLETE.tsv`。脚本不会替 Agent 自动选择最终分辨率或写入生物学标签。

大类或定向 cohort（底层可复用兼容 runner）重聚类使用：

```bash
Rscript annotate-spatial-transcriptomics/scripts/run_seurat_cohort_recluster.R \
  --rds /path/to/full_feature.seurat.rds \
  --membership /path/to/frozen_query_anchor_membership.tsv \
  --out /path/to/project/reclustering/COHORT_ID \
  --cell-id-col cell_id \
  --assay Spatial \
  --resolutions 0.1,0.2,0.3,0.4,0.6 \
  --resolution-contract sheep_ovary \
  --resolution-workers 5
```

membership 至少需要唯一 `cell_id`。锚点模式还需要 `query_or_anchor` 和 `anchor_label`；anchors 与 query 共同参与 SCT/PCA，但邻接图、Leiden、UMAP 和 DEG 必须只在 query 中计算。cohort DEG 使用全基因 `Spatial` LogNormalize 数据，而不是只在 SCT variable features 中寻找 marker。

需要固定的是**同批次技术前处理**，不能固定的是**生物学决策**：

- 最终大类分辨率必须结合 DEG、marker/anti-marker、相邻分辨率稳定性、UMAP 和空间形态选择；
- 每个大类或定向 cohort 重新决定 PCs、k 和最终分辨率；羊卵巢正式候选分辨率始终为 `0.1,0.2,0.3,0.4,0.6`，小型卵母细胞/免疫 cohort 可降低 PC 和 k，但不能改用低于 0.1 的网格；
- `<100` observations 只触发小簇审查，不能自动并入最近 PCA 簇；
- 大类、亚群、合并关系和置信度不能从示例样本复制。

如果 `glmGamPoi` 或 SHA256 依赖不可用，SCT 路线会直接失败并要求修复环境，不会静默切换模型或写出不可验证 manifest。固定技术参数的任何覆盖都必须显式传入 batch exception 和理由。只有 preprocessing manifest 与当前输入/分析集哈希匹配时，已有 SCT 计算才允许复用。

调度资源必须与真实并行度一致：`leidenbase` 的单个 Leiden 优化本身是单线程，Seurat 5.5 通过 `future` 并行的是不同 resolution；因此五档网格通常申请/使用 5 个 worker，而不是占用 64 核串行运行。当前 Seurat 的直接 `FindAllMarkers` 也是逐簇串行；未实现显式并行的作业只申请 1 核，或拆成每个 resolution/cluster 的独立作业后并行。Skill 会记录实际 worker/backend/parallel unit，并用 CPU time / wall time 审计资源是否被浪费。

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

对 Skill 的重要更新必须做 leakage-safe forward test：只给原始输入与 Skill，不给预期聚类或标签；评价它能否冻结无标签证据、完成 cohort/direct return、全细胞 Atlas 一致性与 OOD 复核、残余 QC 回填、失败恢复和完整报告，而不是是否复现旧标签。

## 版本

当前版本：`1.10.0`。发布包的校验和将在 GitHub Release 中提供。

## 许可

MIT License。详见 [LICENSE](LICENSE)。
