# v1.4.1 - Verified Seurat DEG layers and rare-cluster-safe stability

This patch makes the Seurat evidence boundary fail closed.

- A converted raw-count container whose `Spatial@data` equals `counts` can no
  longer enter Wilcoxon DEG.
- `prepare_seurat_full_feature_validation.R` creates a separate, full-feature
  LogNormalize validation-only object with immutable analysis membership and
  SHA256 provenance. It does not mutate or replace the SCT clustering object.
- Initial-cluster and final-label DEG require and validate that manifest for a
  `Spatial` assay.
- Clustering comparison keeps all-observation ARI/AMI and adds a separately
  labeled macro-restricted score. Clusters below 100 observations remain in
  migration, DEG, spatial and rare/technical review outputs.
- Scheduler jobs use a validated `SAMPLE__Pnn_STAGE[_SCOPE]__Ann` name and
  persist that visible name in sample/cohort run registries and reports.

Existing v1.4.0 projects should switch only at a safe pre-DEG checkpoint and
record the new Skill/profile hash; do not hot-swap a frozen Skill during SCT,
PCA or graph jobs.
