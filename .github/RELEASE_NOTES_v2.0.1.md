# v2.0.1

这是 v2 合同的运行边界和发布完整性修订，不改变 `framework_version=2.0.0` 的项目 schema。

主要修复：

- 独立 SCT+BANKSY 从原始 counts 重建，BANKSY 只使用 SCT Pearson residuals；导入表达、降维、聚类和标签不能跨项目复用。
- 零锚点 cohort、TSV.GZ membership、稀疏矩阵求和、坐标重复字段和共享 marker 展开均有 stock 修复与回归约束。
- 全细胞 Atlas 输入必须是 target 与 disjoint current-query heldout rows 的精确哈希绑定并集。
- 最终报告原位切换归一化/绝对 dotplot，发布清单覆盖确认审阅目录、sessionInfo 和工作流事件。
- workflow/biological profile 与候选目录冻结为项目内合同副本；HNAICC 作业名和实际 Python/R 依赖在提交前校验，AIP 使用标准输入提交，completion audit 不再因预注册自身而死锁。

生物学控制保持 v2 语义：开放谱系信号不会因初始阈值或 parent 标签丢失，亚群重聚类允许重建遗漏大类，Atlas 对已定义标签仅作 challenger，对中高置信且原为 QC 的 observation 才可 broad-only 回写。
