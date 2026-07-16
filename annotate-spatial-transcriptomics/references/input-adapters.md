# Input adapters

## Required logical fields

Every adapter must expose unique observation IDs, a count matrix, normalized expression or a reproducible normalization route, feature names, clustering candidates, and reductions. Spatial projects must also expose x/y coordinates.

Record whether each object is full-feature, HVG-only or targeted. Run `scripts/audit_feature_scope.R` against the active tissue profile. Reduced features are valid for graph construction but cannot establish final marker absence, complete open-world lineage discovery or pass a context-specific identity decision.

## Seurat

When a readable full-feature Seurat RDS is available, prefer it as the primary computational object unless the user chooses another backbone or the object fails validation. Prefer RNA/Spatial raw counts for DEG and absolute detection. Record active assay, layer names, object dimensions, reductions and metadata columns. Do not assume SCT residuals are suitable for absolute marker detection. Read `r-first-workflow.md` before using existing Seurat clustering or cohort artifacts.

For a StereoPy `cellbin_PPed` batch-converted RDS, read `seurat-cellbin-preprocessing.md`. The converter may preserve StereoPy reductions but does not thereby create an SCT-processed Seurat object. Use the `Spatial` count layer and generate a fresh SCT/PCA/neighbour/UMAP result with `scripts/run_seurat_sct_preprocess.R`; require its manifest before calling preprocessing comparable across the batch.

## AnnData/Scanpy

Record `X`, `raw`, layer names, `obs`, `var`, `obsm` and `uns`. Determine whether `X` is counts, log-normalized or scaled. Preserve original observation names as strings.

Some StereoPy H5AD files contain nonstandard/custom encoded groups that a given AnnData version cannot read. Inspect HDF5 structure and try a compatible existing environment; never rewrite the source in place. If a corresponding Seurat/SCE RDS is readable, use `scripts/export_r_object_subset_mtx.R` on a frozen membership and `scripts/pack_mtx_h5ad.py` to create a standard, versioned adapter H5AD. Record source/output hashes, dimensions and excluded zero-count observations.

## SingleCellExperiment/BANKSY

Record assay names, `colData`, `reducedDims`, `spatialCoords`, BANKSY parameters and every clustering column. Prefer exported cluster tables when they provide explicit cell IDs and parameters; verify them against the expression object.

## External cluster tables

Require a unique cell ID and cluster column. Accept optional x/y, UMAP, method, resolution and neighborhood fields. Fail if IDs duplicate or if overlap with the expression object is incomplete without an explicit exclusion ledger.

## Existing progress

Search for state ledgers, completion sentinels, manifests, logs and previous reports. Classify them as current, superseded or historical. Never infer that a `DONE` file means annotation is complete; it may only mean clustering finished.

## Reusable cohort runners

Use `scripts/run_seurat_cohort_recluster.R` for configurable Seurat SCT or log-normalized cohort reclustering, `scripts/run_sce_cohort_recluster.R` for SCE/BANKSY-backed cohorts and `scripts/run_scanpy_cohort_recluster.py` for Scanpy. All accept immutable cohort memberships and multiple candidate resolutions; they route zero-count observations to a QC artifact and do not decide the final resolution. BANKSY results can enter through exported cluster tables, while new BANKSY subset analysis should be configured from the current SCE/spatial object rather than copied from a whole-tissue run.

For an existing SingleCellExperiment/BANKSY clustering, use `scripts/run_sce_initial_cluster_evidence.R`. The shared dotplot, annotation-map and spatial-gene scripts accept both Seurat and SCE/SummarizedExperiment objects; for SCE specify assay names when `counts` plus `normcounts`/`logcounts` are not present.
