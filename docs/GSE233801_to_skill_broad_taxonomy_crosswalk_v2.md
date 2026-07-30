# GSE233801 与羊卵巢 Skill broad taxonomy 对应表

## 结论

旧 `sheep_ovary_GSE233801_v1` 把 c1/c3/c6/c10/c15/c16/c17 合并为 `Vascular/endothelial`，并把 c8、c12 与 c2 合并到 `Stromal/perivascular`。这不是 GSE233801 缺少 Endothelial、Pericyte/mural 或 Smooth muscle，而是旧 Atlas 索引的发布层级过粗。

split-wall v2 使用同一 GSE233801 原始参考、同一独立 `res0.4` 审阅结果重新构建原型。它不读取任何 query 标签，也不通过重命名旧原型伪造新类别。

## 论文公开的大类标签

GSE233801 论文正式发布的五类体细胞为：

| 论文标签 | Skill 对应 |
|---|---|
| Granulosa | Granulosa |
| Endothelial | Endothelial |
| Stromal | Stromal/mesenchymal |
| Perivascular | 不能整体直译；独立重建后分别对应 c12 Pericyte/mural 与 c8 Smooth muscle |
| Immune | Immune |

论文还保留了若干未命名/未知成分。它们不能因为名称为空而被丢弃，也不能按 marker 单基因直接补名；split-wall v2 只使用独立全特征审阅后足够连贯的成分。

## GSE233801 可重建的大类/成分

| GSE233801 重建身份 | 参考 cluster | 细胞数 | Skill broad 对应 | v2 权限 |
|---|---:|---:|---|---|
| Blood/microvascular/capillary/arterial/lymphatic endothelial 及 endothelial state | 1, 3, 6, 10, 15, 16, 17 | 20,955 | Endothelial | `supported`；仅 broad，不能直接写 Lymphatic fine |
| Fibroblast stroma | 2 | 6,233 | Stromal/mesenchymal | `supported` |
| Granulosa | 7 | 3,001 | Granulosa | `supported` |
| Mature smooth muscle | 8 | 2,209 | Smooth muscle | `supported` |
| T lymphoid、myeloid/APC、inflammatory myeloid | 9, 11, 19, 21 | 4,670 | Immune | `supported` |
| Pericyte | 12 | 1,684 | Pericyte/mural | `supported` |
| Ovarian epithelial | 18 | 605 | Epithelial/mesothelial | `challenge_only` |
| Steroidogenic theca-like | 22 | 139 | Theca | `challenge_only` |
| Structural theca-stromal follicular wall | 4 | 4,909 | Stromal/Theca 边界审计 | 不建直接救回原型 |
| Peripheral glial-like | 20 | 191 | Glial/Schwann-like | `challenge_only`；冻结参考对象中无原型 |
| Low-information/cycling | 5, 13, 14 | 7,492 | 无 | 排除 |

## Skill 全部 broad 的映射能力

| Skill broad | GSE233801 对应 | Atlas 用途 |
|---|---|---|
| Granulosa | c7 | 类别校准通过后，可救回未标注 broad；已有标签只做一致性挑战 |
| Oocyte | 无可靠参考簇 | `unsupported`；完全依赖 query canonical-cluster 规则 |
| Theca | c22 稀有类固醇样成分；c4 为结构性混合边界 | `challenge_only`；不得直接救回或以空间位置扩张 |
| Luteal | 无独立黄体簇 | `unsupported`；依赖阶段、query 黄体程序与空间形态 |
| Stromal/mesenchymal | c2 | 类别校准通过后可直接救回；c4 不加入直接原型 |
| Smooth muscle | c8 | 类别校准通过后可直接救回；必须与 c12 Pericyte/mural 独立竞争 |
| Pericyte/mural | c12 | 类别校准通过后可直接救回；不得回退为泛血管标签 |
| Endothelial | c1/c3/c6/c10/c15/c16/c17 | 类别校准通过后可直接救回；c15 只增强 broad Endothelial，不赋予 fine |
| Immune | c9/c11/c19/c21 | 类别校准通过后可直接救回 |
| Epithelial/mesothelial | c18 | `challenge_only`；小簇暂不能单凭参考自分类获得直接写回权限 |
| Glial/Schwann-like | c20 的证据存在但未进入冻结参考 | `challenge_only/not_evaluable`；不能直接救回 |

## 路由约束

1. `supported` 表示具备候选原型，不等于自动写标签；仍需当前 query 的独立、类别级 held-out calibration 达到目标精度，且观测未标注、非 OOD、无 ontology conflict。
2. `challenge_only` 只能提出欠召回/冲突问题，不能直接写回。
3. 已有 broad 与 Atlas 冲突时只能重开来源 subcluster/cohort 的生物学复核，不能被 Atlas 静默覆盖。
4. 旧 `Vascular/endothelial` 或 `Stromal/perivascular` 原型不能通过别名拆成新大类；split-wall v2 必须由原始 reference cluster 重新建模。
5. Atlas 永远没有 fine-label authority；GSE233801 c15 即使为淋巴内皮，也只能映射到 broad `Endothelial`。

## v2 reference held-out 诊断

以下结果只验证 GSE233801 内部可分性，不替代当前 query 的独立校准：

| 原型 | Recall | Precision | 解释 |
|---|---:|---:|---|
| Endothelial | 99.40% | 95.95% | 可作为 supported broad 原型 |
| Smooth muscle | 86.40% | 95.79% | 精度高，但部分 c8 被 mural 原型吸收；query 中仍需完整成熟收缩核心 |
| Pericyte/mural | 93.11% | 85.96% | 未阈值时混入 Smooth muscle；在 held-out 高置信子集可达到 95.61% precision，因此保留 supported 原型但强制类别级校准 |
| Stromal/mesenchymal | 99.40% | 93.60% | 高置信阈值后 precision 96.49% |
| Granulosa | 97.20% | 99.79% | supported |
| Immune | 96.80% | 99.79% | supported |
| Epithelial/mesothelial | 94.70% | 99.31% | 内部可分，但因小簇及历史过召回风险仍为 challenge-only |
| Theca | 94.12% | 94.12% | 仅 34 个 held-out、来源样本偏倚，仍为 challenge-only |
