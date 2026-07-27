# Adaptive clustering selection

## Shortlist evidence

Evaluate:

1. Number and distribution of cluster sizes.
2. Fraction in tiny clusters and singleton/microclusters.
3. ARI or contingency consistency across adjacent resolutions and neighborhood parameters.
4. Whether higher resolution adds lineages or only splits states/technical gradients.
5. Canonical marker coherence and DEG specificity.
6. UMAP separation without relying on UMAP alone.
7. Spatial coherence, fragmentation and agreement with plausible morphology.
8. Ability to identify rare but supported biology without manufacturing many uninterpretable clusters.

## Selection principle

Choose the complete-grid candidate with the best integrated evidence for the current biological question. First avoid losing stable lineages or true subtypes, then avoid state-only/technical fragmentation. Lower complexity is only a tie-breaker when evidence is otherwise equivalent.

Automated scores only rank candidates. Freeze a candidate only after marker and spatial review. Record rejected alternatives and reasons.

Always report adjacent-resolution ARI/AMI on the sampled common observations
with every cluster included. Also report a separate macro-restricted score in
which an observation contributes only when its cluster has at least 100
members in both partitions. The restricted score is only a broad-resolution
ranking aid. Export the full cluster migration table and a
`small_cluster_review` audit: size alone never deletes, relabels or suppresses
a cluster from DEG, spatial evidence or rare-lineage/technical adjudication.

Use `scripts/rank_cohort_resolutions.py --question-mode broad_purity_audit|targeted_mixture` to shortlist cohort resolutions from cluster size, adjacent ARI and profile-program interpretability. For spatial cohorts, run `scripts/summarize_spatial_cluster_morphology.py` to quantify components, noise and largest-component fraction. Interpret these values against the profile: repeated follicles or a vascular network can legitimately be noncompact, while arbitrary salt-and-pepper splitting is a warning.

## BANKSY

Neighborhood size changes spatial smoothing and cannot be treated as another resolution value. Compare resolution within each k, then compare plausible candidates across k. High k may merge small anatomical structures; low k may fragment continuous tissue.

For a fresh verified same-batch SCT+BANKSY result, the whole-tissue grid is `0.2,0.4,0.6,0.8`. For a reusable upstream BANKSY result, the candidate grid is exactly the complete grid present in the frozen input and recorded in `annotation_contract.json`. Do not replace either BANKSY grid with the query-only cohort grid `0.1,0.2,0.3,0.4,0.6`. BANKSY resolution, BANKSY neighborhood `k`, Leiden/Louvain resolution and graph-neighbor `k` are separate parameters and must never be used interchangeably in evidence or prose.

When BANKSY memberships are the requested whole-tissue starting point, select one stable cohort-partition resolution after reviewing the complete formal grid. Preserve major spatial structure and avoid technical fragmentation; do not require the first pass to separate every lineage. Candidate programs and embedded signals become provisional `watch` evidence for second-round cohorts, not first-pass labels. Cluster count, UMAP appearance and a target number of clusters are shortlist aids only.

Record weak coherent programs seen at any reviewed resolution in the watch ledger. If no first-pass resolution separates every lineage, choose the most stable anatomical partition and rely on full-catalog second-round reconstruction; never silently close the missing lineage.

## Seurat/R-first

When a full-feature Seurat RDS exists, make Seurat the default whole-tissue and cohort-reclustering backbone. Use SCTransform for graph construction when appropriate, retain the full RNA/Spatial counts for marker validation, and compare a coarse broad-pass resolution grid rather than inheriting a BANKSY or Scanpy resolution. Existing Seurat grids may be reused after membership/hash and artifact validation, but historical biological labels remain blinded.

For mesenchymal-rich tissues, candidate selection must demonstrate whether generic stroma, mesenchymal-progenitor-like, mature smooth muscle, pericyte/mural and endothelial programs are separable. Do not prefer a resolution merely because it creates all named literature classes; require stable multi-gene and spatial evidence.

## Second-round cohort reclustering

Repeat biological selection independently inside every initial-cluster cohort. For sheep ovary, keep `0.1,0.2,0.3,0.4,0.6` while adapting selected resolution, PCs and k from current-cohort evidence. Compare the selected partition with the nearest lower/higher resolutions; at a grid edge use the two nearest candidates. Do not require two higher resolutions and do not reuse whole-tissue BANKSY graph parameters. Other profiles may declare different grids.

The second round optimizes stable broad identity before subtype yield. If many small subclusters are marker-poor or unresolved while a lower neighboring resolution restores coherent parent or cross-lineage broad programs without hiding a validated minority component, treat the selected resolution as overfragmented and prefer the lower stable candidate. Do not force a named subtype from weak markers. A subcluster with a supported broad background and no independent child discriminator remains broad-only; an independently supported different lineage returns to that broad even when no fine subtype is available.
