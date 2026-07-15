# Seurat cellbin preprocessing contract

Use this contract when a project contains Seurat RDS files converted from the
same StereoPy `cellbin_PPed` production batch. The converted RDS is an input
container, not a preprocessed Seurat analysis object.

## Converted-input contract

The batch converter imports `exp_matrix@raw` as the Seurat `Spatial` counts
layer, preserves observation metadata and coordinates, and may copy existing
StereoPy PCA/UMAP reductions for provenance. It does **not** run
`NormalizeData`, `FindVariableFeatures`, `ScaleData` or `SCTransform`.

Therefore:

- start R clustering from the `Spatial` raw-count layer;
- do not treat an imported StereoPy PCA/UMAP as an R preprocessing result;
- keep the full `Spatial` assay for marker detection, anti-marker checks and
  rare-lineage validation;
- create a fresh `SCT` assay and fresh PCA/neighbour/UMAP reductions for each
  whole-tissue analysis;
- do not infer preprocessing from the presence of an assay or reduction name;
  require a matching preprocessing manifest.

## Frozen whole-tissue profile for this production batch

For samples intended to be directly comparable with the validated R workflow,
use the following fixed numerical preprocessing profile unless a recorded
batch-level exception is approved:

| Stage | Parameter | Value |
|---|---|---|
| entry QC | minimum counts | `100` |
| entry QC | minimum detected features | `75` |
| entry QC | rule | `nCount_Spatial >= 100 AND nFeature_Spatial >= 75` |
| high-count review | flag | above the sample's `99.9%` count quantile; flag only |
| mitochondrial QC | hard filter | none unless sheep mitochondrial symbols are reliably mapped |
| doublet handling | hard deletion | none at entry; preserve a review flag |
| SCTransform | input/new assay | `Spatial` / `SCT` |
| SCTransform | flavour/method | `v2` / `glmGamPoi` |
| SCTransform | variable features | `3000` |
| SCTransform | fitting observations | `min(50000, n_observations)` |
| SCTransform | memory/feature policy | `conserve.memory=TRUE`, `return.only.var.genes=TRUE` |
| PCA | components computed | `50` on SCT variable features |
| neighbours | components used | first `30` |
| neighbours | graph | `k=30`, Annoy, `n.trees=50`, cosine |
| clustering | algorithm | Leiden (`algorithm=4`) |
| broad grid | candidate resolutions | `0.1,0.2,0.3,0.4,0.6` |
| Leiden execution | parallel unit | one single-threaded Leiden optimization per resolution; run candidate resolutions concurrently with `future` workers |
| UMAP | geometry | `n.neighbors=30`, `min.dist=0.3`, cosine |

The reference implementation is `scripts/run_seurat_sct_preprocess.R`. It
writes a cell-scope table, preprocessing manifest, all candidate cluster
memberships and an SCT/PCA/UMAP Seurat object. It deliberately does not choose
a final resolution or assign biological labels.

The runner enforces these values. Any numerical override requires both
`--allow-batch-exception` and a substantive `--batch-exception-reason`; the
manifest records the exception. Missing SHA256 support is fatal rather than
silently emitting an unverifiable manifest.

Computational parallelism is provenance, not a biological parameter. On Linux compute nodes, set `--resolution-workers` to at most the number of candidate resolutions and request the same number of scheduler CPUs; `--resolution-future-plan auto` selects `multicore` when supported. One Leiden optimization remains single-threaded, but the five candidate resolutions run concurrently. The same worker count is exposed to `uwot` UMAP. Do not request 64 CPUs for a five-resolution grid or for a sequential `FindAllMarkers` call. When whole-grid evidence is produced by a wrapper, prefer one one-CPU job per resolution plus a dependency-gated aggregation job unless the wrapper contains verified cluster-level parallelism.

## Separate full-feature validation object

The converted `Spatial@data` layer may be missing or may be an exact copy of
raw counts. Such a layer is not valid input for Wilcoxon DEG, even though the
SCT graph and clusters remain valid. Do not normalize or otherwise mutate the
frozen SCT clustering object to repair this.

On a compute node, run `scripts/prepare_seurat_full_feature_validation.R` from
the immutable raw-count RDS plus the frozen analysis-set membership. It creates
a separate all-feature `LogNormalize(scale.factor=10000)` object, removes all
reductions, marks the object as `full_feature_deg_marker_validation_only` and
writes input, membership, analysis-set and output hashes in a validation
manifest. The object is explicitly ineligible for clustering.

Pass that object and `--validation-manifest` to
`scripts/run_initial_cluster_evidence.R` and `scripts/run_final_label_deg.R`.
For assay `Spatial`, both scripts fail closed when the manifest is absent,
points to another object, has a mismatched analysis-set hash, or when
`data == counts`. Detection percentages and marker presence continue to use
raw counts; Wilcoxon rank evidence uses the verified LogNormalize `data`
layer. This validation job can be deferred until DEG/marker evidence is
needed, which is why a running SCT/PCA job may not encounter the issue yet.

## Parameters that remain adaptive

Batch harmonisation does not mean copying biological decisions. The following
must be selected from each sample or broad/targeted cohort:

- the final broad-class resolution;
- cohort membership and anchor composition;
- cohort-specific resolution grids and the selected cohort resolution;
- the number of PCs for small cohorts when fewer than 30 are supported;
- biological labels, merges, subtypes and confidence.

Select the lowest-complexity broad resolution that preserves defensible
lineages using DEG, absolute marker/anti-marker evidence, cluster-size
distribution, adjacent-resolution stability, UMAP and spatial morphology.

Do not automatically merge or reassign a cluster only because it contains
fewer than 100 observations. Rare oocyte, immune, endothelial or epithelial
populations may be small. Record a `small_cluster_review` flag and adjudicate it
biologically.

## Broad-class and targeted-cohort profile

Each cohort reclustering begins again from the parent object's `Spatial` counts. When
anchors are present, anchors and query are jointly SCT/PCA fitted, while the
neighbour graph, Leiden clusters, UMAP and DEG are query-only.

Use these defaults:

- `SCTransform(vst.flavor="v2", method="glmGamPoi")`;
- `ncells=min(50000, n_query_plus_anchors)`;
- `conserve.memory=TRUE`, `return.only.var.genes=TRUE`;
- 1,500 SCT variable features for a cohort below 1,000 observations, otherwise
  3,000; a balanced-anchor round may use 2,500 below 5,000 joint observations
  and 3,500 otherwise when the richer anchor boundary requires it;
- at most 30 PCs for ordinary cohorts, at most 25 for shallow immune cohorts, and
  at most 15 for the Oocyte targeted cohort, constrained by cohort size and feature
  rank;
- adaptive `k=min(30, max(5, floor(sqrt(n_query))))`; cap strict oocyte
  candidates at 20;
- Annoy cosine neighbours with 50 trees and Leiden algorithm 4;
- UMAP cosine with `min.dist=0.3`; `0.2` may be used for strict small oocyte
  candidates;
- for sheep ovary, run the complete formal grid `0.1,0.2,0.3,0.4,0.6` in every broad/targeted cohort and choose the final value from current-cohort evidence. Pass `--resolution-contract sheep_ovary`; adapt PCs and k, not the formal resolution candidates. A value below 0.1 is a graph-repair trigger, not an analysis candidate. Final QC is not reclustered.

The portable `scripts/run_seurat_cohort_recluster.R` exposes these controls. Its
default SCT route fails closed if `glmGamPoi` is unavailable; silently changing
the SCT model would break batch comparability.

## Required provenance

Every run must record the input checksum, assay/count layer, exact QC rule,
SCT flavour and method, variable-feature count, SCT fitting sample size, PCA
components computed/used, neighbour metric and implementation, candidate
resolutions, seeds, software versions, analysis-set observation IDs and their
hash. Keep initial QC exclusions in the full ledger; do not relabel them as a
post-clustering biological population.
