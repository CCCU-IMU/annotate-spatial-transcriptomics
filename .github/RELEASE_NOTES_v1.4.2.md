# Annotate Spatial Transcriptomics v1.4.2

This patch release aligns scheduler CPU requests with real Seurat parallelism and clarifies the sheep-ovary Oocyte spatial gate.

## Resource-aware Seurat execution

- `run_seurat_sct_preprocess.R` and `run_seurat_pool_recluster.R` now expose `--resolution-workers` and `--resolution-future-plan`.
- One Leiden optimization remains single-threaded, while candidate resolutions run concurrently through Seurat 5.5 `future` workers.
- `uwot` UMAP receives the same recorded worker count.
- Scheduler guidance now forbids reserving many idle CPUs for single-threaded Leiden or direct `FindAllMarkers`; use one CPU or independent per-resolution/per-cluster jobs plus validated aggregation.
- Run manifests record requested/used workers, backend and UMAP threads, and job review includes CPU-time/wall-time auditing.

## Sheep-ovary Oocyte location rule

Cortical, subcortical, peripheral or section-edge location is never negative evidence against Oocyte. Small or primordial oocytes may be cortical. Location alone remains insufficient positive evidence and does not relax multi-module markers, somatic anti-programs, spatial-object grouping or strict candidate reclustering.
