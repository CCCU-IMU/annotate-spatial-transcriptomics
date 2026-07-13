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

Choose the lowest-complexity candidate that preserves supported broad lineages and useful rare populations. A broad-pass clustering can be deliberately coarser than subtype-pass clustering.

Automated scores only rank candidates. Freeze a candidate only after marker and spatial review. Record rejected alternatives and reasons.

Use `scripts/rank_pool_resolutions.py` to shortlist pool resolutions from cluster size, adjacent ARI and profile-program interpretability. For spatial pools, run `scripts/summarize_spatial_cluster_morphology.py` to quantify components, noise and largest-component fraction. Interpret these values against the profile: repeated follicles or a vascular network can legitimately be noncompact, while arbitrary salt-and-pepper splitting is a warning.

## BANKSY

Neighborhood size changes spatial smoothing and cannot be treated as another resolution value. Compare resolution within each k, then compare plausible candidates across k. High k may merge small anatomical structures; low k may fragment continuous tissue.

## Subset reclustering

Repeat selection inside every broad pool. Do not copy whole-tissue resolution, BANKSY k, number of PCs or neighbor settings into a subset without re-evaluation.
