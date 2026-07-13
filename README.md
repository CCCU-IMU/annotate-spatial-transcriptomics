# Annotate Spatial Transcriptomics Skill

[![Validate](https://github.com/CCCU-IMU/annotate-spatial-transcriptomics/actions/workflows/validate.yml/badge.svg)](https://github.com/CCCU-IMU/annotate-spatial-transcriptomics/actions/workflows/validate.yml)

面向 Codex Agent 的空间转录组/单细胞转录组迭代式注释 Skill。它不是把一次成功项目的参数和标签复制到新样本，而是让 Agent 根据当前物种、组织、平台、数据复杂度和生物学问题，自主完成：输入发现、聚类选择、大类注释、多路线待定义细胞处理、按需亚群注释、状态追踪、失败恢复、用户确认以及可审计 HTML 报告。

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
  | bash -s -- --ref v1.2.0
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

### 1. 首次消息一次给全背景

至少提供：

- 输入文件或输入目录、输出项目目录；
- 物种、组织、发育阶段/处理条件；
- 测序平台与观测单位（真实细胞、cellbin、spot 等）；
- 当前已有的聚类结果及其生成方法；
- 主要生物学问题、需要重点审查的稀有谱系；
- 可用 R/Python 环境、调度系统和计算资源；
- 参考 atlas 的优先级（如已有），但不要把参考标签当作真值。

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

请自主发现输入、选择适合当前数据的聚类强度、投递并监控任务、检查日志、修复失败、维护状态文件并持续迭代。不要复制示例项目的参数或标签。空间数据优先获得可靠大类，亚群只在证据充分时定义。所有待定义细胞必须走完适用的多路线策略；在最终耗时作图前，先向我展示冻结的第一版注释统计并等待确认。
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
最终大类分辨率、池子分辨率及生物学标签仍应针对当前样本自适应判断。
```

### 2. 允许 Agent 运行完整迭代，不要把任务简化成 cluster 重命名

标准流程是：

```text
输入与表达层审计
  -> 自适应选择大类聚类
  -> 基于 DEG、marker/anti-marker、UMAP、空间形态建立高置信锚点
  -> 生成待定义池并按失败原因分流
     -> A：大而连续、信号可用的群体：平衡锚点重聚类
     -> B：局部混合/界面群：定向重聚类，适用时再用校准 RCTD/参考辅助
     -> C：低复杂度群：完整 QC 池锚点重聚类，剩余拒绝再做校准 atlas 救回
  -> 回写大类，按需在生物学一致的大池内做亚群
  -> strict / inclusive / display 三视图
  -> 完成门审计
  -> 用户确认冻结注释
  -> 最终 DEG、双层树点图、空间图和 HTML 报告
```

核心原则：

- 大类先于亚群；空间数据不为了“树更深”而强制细分。
- 聚类分辨率在每个大池中重新选择，不能照搬全局参数或示例参数。
- HVG/BANKSY 特征可以用于聚类，但最终 marker、anti-marker 和稀有细胞判定应回到全基因表达对象。
- 一个基因、一个参考标签或空间邻近都不能单独决定身份。
- ECM-rich、contractile、cortical、ambient、low-RNA 等是状态标签，不应替代生物学大类。
- 已关闭池有不可变 membership 和来源链，不能在后续上下文中被重复聚类或重复注释。

### 3. 把多路线待定义细胞处理当作流程核心

不要用“未注释比例低”作为完成标准，也不要把所有弱信号细胞直接过滤掉。

- **路线 A：锚点重聚类。** 对数量较大、连续且仍有生物学信号的待定义细胞，使用各大类的独立高置信锚点，做深度平衡、reference-only anchor / query-only graph 的重聚类；结合 cluster DEG、细胞级程序纯度、anti-program 和空间连续性决定大类回归。
- **路线 B：界面审查。** 对局部颗粒–卵泡膜、血管–基质等界面群，先定向重聚类。RCTD/参考映射优先级较低，只有极高置信且有独立 marker、anti-marker、分辨率和空间证据时才支持亚群回归；高置信可回归大类；中低置信继续进入 atlas/内部锚点审查或保留为界面。
- **路线 C：QC holdout 救回。** 低复杂度群先作为“聚类后 QC 保留池”，必须对完整 QC 池做一次锚点重聚类；仍无法解释的细胞才进入 atlas + 内部锚点 + marker/anti-marker + 观测密度空间证据的联合校准。

Atlas 阈值必须通过当前项目的 query-depth-matched held-out anchors 校准。默认目标精度为 `moderate-or-higher >= 0.90`、`high >= 0.95`，它们是**校准目标精度**，不是通用 raw score 截断值。输出应满足：

```text
moderate_or_higher_n = high_n + moderate_only_n
```

经过独立证据复核的 high 和 moderate-only 均可回归**大类**；低于 moderate 的保留为 reject/review。Atlas 救回的 broad-only 细胞必须设置 `fine_anchor_eligible=false`，不能反向参与精细 marker 或亚群锚点发现。

### 4. 对稀有、易污染谱系使用更严格的上下文门

卵母细胞、浆细胞、淋巴内皮等稀有类型不能因为单一强 marker 就大范围定义。应同时要求：

- 多基因正向程序；
- 排除邻近高丰度谱系的 anti-program；
- 在全基因表达层验证，而非仅 HVG；
- 与组织学相符的局灶空间结构；
- 在候选池重聚类后仍形成稳定、纯净的小群；
- 明确区分 cellbin/spot 数量与真实生物学细胞数量。

未通过上下文门的候选应回滚到合适的体细胞/界面池，而不是保留一个看似“稀有”的宽松标签。

### 5. 只在一个关键节点人工确认

Agent 应自主完成常规分辨率比较、作业投递、监控、日志排错、重投和状态回写。建议只保留一个强制人工节点：完成门通过后，Agent 展示冻结的第一版注释统计、未定义/QC/界面数量、稀有细胞候选及 atlas/RCTD 回归情况，用户确认后再生成耗时的最终图和报告。

如果冻结注释之后又修改 ledger、route 或 completion gate，旧确认自动失效，必须重新确认。这能避免在注释仍变化时反复重画全套报告。

### 6. 用结果合同验收，而不是只看一张 UMAP

最终交付至少应包括：

- 大类与亚群的独立 DEG；
- 大类和亚群两套树状 marker dotplot；
- 两个层级分别包含 canonical marker 与当前数据特异 marker；
- dotplot 同时输出 PNG、PDF 和源 TSV；点大小和颜色按基因内部归一化展示，但源表保留绝对检测率和平均表达量；
- 大类/亚群 UMAP、全空间图、逐节点高亮网格；
- 按支持细胞类型分组的 marker 空间表达图；
- strict、inclusive、display 三套互斥统计与空间概览；
- 可展开注释树、路线/阈值/结局面板、中文详细流程时间线和原始状态记录；
- cell-level ledger、cluster decision ledger、pool/run/route registry、session info、manifest 和 checksums；
- completion gate 与 release audit 均通过。

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
| UMAP | 30 neighbors，`min.dist=0.3`，cosine |

全组织统一前处理可直接运行：

```bash
Rscript annotate-spatial-transcriptomics/scripts/run_seurat_sct_preprocess.R \
  --rds /path/to/sample.seurat.rds \
  --out /path/to/project/reclustering/seurat_sct \
  --sample SAMPLE_ID \
  --assay Spatial \
  --resolutions 0.1,0.2,0.3,0.4,0.6 \
  --future-globals-max-gb 100
```

脚本会输出：

- 保留全部输入观测的 `analysis_scope.tsv.gz`，明确区分 `analysis_set` 与 `excluded_initial_qc`；
- 每个候选分辨率的 cluster membership 和 cluster count；
- 新生成的 SCT/PCA/UMAP Seurat RDS；
- UMAP、空间坐标、输入 SHA256、分析集 SHA256、全部前处理参数和 `sessionInfo`；
- `RUN_COMPLETE.tsv`。脚本不会替 Agent 自动选择最终分辨率或写入生物学标签。

池子或平衡锚点重聚类使用：

```bash
Rscript annotate-spatial-transcriptomics/scripts/run_seurat_pool_recluster.R \
  --rds /path/to/full_feature.seurat.rds \
  --membership /path/to/frozen_query_anchor_membership.tsv \
  --out /path/to/project/pools/POOL_ID \
  --cell-id-col cell_id \
  --assay Spatial \
  --resolutions 0.1,0.2,0.4,0.6
```

membership 至少需要唯一 `cell_id`。锚点模式还需要 `query_or_anchor` 和 `anchor_label`；anchors 与 query 共同参与 SCT/PCA，但邻接图、Leiden、UMAP 和 DEG 必须只在 query 中计算。池子 DEG 使用全基因 `Spatial` LogNormalize 数据，而不是只在 SCT variable features 中寻找 marker。

需要固定的是**同批次技术前处理**，不能固定的是**生物学决策**：

- 最终大类分辨率必须结合 DEG、marker/anti-marker、相邻分辨率稳定性、UMAP 和空间形态选择；
- 每个父池重新决定 PCs、k、候选分辨率和最终分辨率；小型卵母细胞/免疫池可降低 PC 和 k；
- `<100` observations 只触发小簇审查，不能自动并入最近 PCA 簇；
- 大类、亚群、合并关系和置信度不能从示例样本复制。

如果 `glmGamPoi` 不可用，SCT 路线会直接失败并要求修复环境，不会静默切换模型。只有 preprocessing manifest 与当前输入/分析集哈希匹配时，已有 SCT 计算才允许复用。

## 仓库结构

```text
annotate-spatial-transcriptomics/
  SKILL.md                 # Agent 入口与强制工作流
  agents/openai.yaml       # Skill 元数据
  references/              # 聚类、路由、状态、报告、测试规范
  references/profiles/     # 组织/物种背景 profile（不是标签映射）
  assets/                  # marker、锚点与多通道救回配置模板
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
python scripts/validate_repo.py
bash -n install.sh
```

对 Skill 的重要更新必须用新 Agent 做 forward test：只给它原始输入与 Skill，不给预期聚类选择和注释答案；评价其是否能自主完成多轮分流、恢复失败、控制稀有谱系和产出完整报告，而不是是否复现某个旧项目的标签。

## 版本

当前版本：`1.2.1`。发布包的校验和将在 GitHub Release 中提供。

## 许可

MIT License。详见 [LICENSE](LICENSE)。
