# Annotate Spatial Transcriptomics Skill v1.2.0

首个正式公开版本，提供可移植的 Codex Skill、通用分析脚本、组织背景 profile、状态机规范和报告合同。

本版本的核心是多路线待定义细胞处理：大而连续的生物学池进行平衡锚点重聚类，局部混合群进行界面审查，完整 QC holdout 池先重聚类再执行校准 atlas 多通道救回。空间数据以可靠大类为主要终点，亚群只在证据充分时定义。

安装：

```bash
curl -fsSL https://raw.githubusercontent.com/CCCU-IMU/annotate-spatial-transcriptomics/main/install.sh | bash
```

固定版本：

```bash
curl -fsSL https://raw.githubusercontent.com/CCCU-IMU/annotate-spatial-transcriptomics/main/install.sh | bash -s -- --ref v1.2.0
```

下载后可使用 `CHECKSUMS.sha256` 验证 ZIP/TAR 发布包。
