# Annotate Spatial Transcriptomics Skill v1.3.0

本版本把默认控制中心调整为 Seurat/R-first，同时保留 Scanpy、BANKSY、SCE 和外部聚类适配器。它复用已验证项目的迭代策略，不复制其聚类参数或注释结果。

主要更新：

- 将“论文候选分类目录”“分析 parent pool”“最终发布大类”拆成三个不可互相复制的层级。
- 基于近年羊卵巢整组织单细胞与多组学研究，更新 sheep-ovary profile、marker/anti-marker、命名和证据权重。
- `Theca` 仅保留可信类固醇生成群；成熟 `Smooth muscle`、`Pericyte/mural` 与 `Mesenchymal progenitor-like` 需要分别通过证据门。
- 新增 fail-closed release taxonomy audit，防止 pool 名、`Theca/follicular wall` 等兜底名和 QC/interface 状态进入生物学大类树。
- 最终确认会绑定 taxonomy audit 哈希；ledger 或 profile 改变后必须重新审计和确认。
- 空间数据继续以可靠大类为主，计算 cluster 不自动成为亚群。

安装最新版：

```bash
curl -fsSL https://raw.githubusercontent.com/CCCU-IMU/annotate-spatial-transcriptomics/main/install.sh | bash
```

安装固定版本：

```bash
curl -fsSL https://raw.githubusercontent.com/CCCU-IMU/annotate-spatial-transcriptomics/main/install.sh | bash -s -- --ref v1.3.0
```
