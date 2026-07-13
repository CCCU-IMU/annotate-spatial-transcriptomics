# Multi-route unresolved-population controller

The controller routes by failure mode, not by the name of the source graph cluster. Preserve every source and state tag before routing.

## Route A: biological parent-pool anchor reclustering

Use for a large population with usable counts and one or more plausible lineages.

1. Freeze query membership and balanced, trustworthy reference-only anchors separately.
2. Membership must contain `query_or_anchor`, `anchor_label`, `source_key`, `state_tags`, `spatial_tags` and `qc_tags`.
3. Fit normalization and PCA jointly. Build neighbors, clusters, UMAP, DEG and cluster statistics on query observations only. Anchors never enter query DEG or final pool counts.
4. Test a pool-specific resolution grid. Review adjacent-resolution migration/ARI, marker and anti-marker completeness, source/state/QC composition, PCA distance to anchors and spatial morphology.
5. Each subcluster becomes one of: fine definition, broad-only return, canonical-pool reroute, interface, QC/technical retention or strict rare review.
6. Roll back source-, depth-, ECM-, stress-, cell-cycle- or ribosomal-driven splits.

An analysis that subsets only query observations is not anchor-assisted even if its directory or route name says so.

## Route B: local mixed interface

Use for a spatially local narrow-layer or mixed-lineage population after Route A fails to produce a stable separation.

1. Run targeted anchor reclustering first.
2. Record deconvolution applicability from query size, UMI/depth, reference types, gene overlap, interface morphology and presence of separable held-out anchors.
3. If applicable, depth-match held-out anchors to the query and compare parameter settings. Calibrate top weight/probability and top-vs-second margin to a predeclared precision target. Record coverage, precision and recall.
4. Apply tiered evidence. Extreme-confidence RCTD plus independent full-feature marker/anti-marker, stable reclustering and spatial evidence may support a fine return. High-confidence RCTD supports broad-only return with `fine_anchor_eligible=false`.
5. Send medium/low-confidence RCTD cells to calibrated atlas/internal-anchor mapping rather than closing them as interface/QC. Atlas may rescue broad-only labels; it cannot promote a fine subtype. Only cells rejected after this fallback remain anatomical interface, QC or technical review.
6. A local interface must be spatially restricted and small relative to the analysis set. A large or diffuse reject population must be reopened as a canonical broad parent pool and/or atlas fallback; the existence of an RCTD result does not make it a valid terminal interface.

`not_applicable` requires a validation artifact and concrete failed criteria; a narrative statement alone cannot close the branch.

## Route C: post-clustering low-information/QC holdout

Keep initial QC exclusion separate. For a nontrivial post-clustering holdout:

1. Freeze the complete QC query and its source components.
2. Run a full QC-pool Route-A-style anchor reclustering before atlas mapping. Use broad anchors, query-only graph/DEG and a coarse grid aimed at broad rescue, not invented fine states.
3. Preserve technical components and singleton/microclusters as QC even if they land near a large biological cluster.
4. For remaining rejects, audit reference adequacy and depth-match held-out internal anchors. Combine external reference, internal anchors, marker programs and spatial support.
5. Compute an observed-density spatial prior from all eligible anchors. A class-balanced spatial prior may be used only as a sensitivity analysis and must not replace observed tissue density.
6. Rescue calibrated broad-only labels into `inclusive`; keep conservative/high-confidence biological definitions in `strict`. Rejects remain QC. All rescues are `fine_anchor_eligible=false`.

Atlas mapping alone is not a completed Route C for a large QC pool.

## Route D: strict rare/context-gated identity

Use full-feature multi-module positive evidence, contradictory resident programs, spatial objects/morphology, adjacent/background comparison and strict-candidate reclustering. Return contamination/adjacent cells to canonical somatic pools with explicit ambient/adjacent tags.

## Route E: diagnostic supervised state review

Use only after the appropriate unsupervised route for ambient-high, LOC-dominated, baseline-low-RNA, vascular–stromal transition or graph-outlier failures. Compare full and nuisance-excluded models or depth-matched anchors. It may support broad return or explicit state retention; it must not manufacture a fine subtype.

## Required route registry outcomes

Every applicable or reviewed-inapplicable branch records query/anchor counts, membership hash, selected parameters/resolution, validation artifact, outcome counts, retained state, fine-anchor eligibility and no-repeat policy in `route_attempt_registry.tsv` and `branch_control_board.tsv`.

For Route B additionally record the extreme/high/medium-low confidence counts, independent evidence gates, atlas/internal-anchor fallback membership, fallback route ID and final reject count. Do not merge these confidence tiers in one narrative total.

## Closure

A branch closes only when its entire frozen membership is mutually exclusively accounted for and every applicable route has a validated outcome. Large broad-only, interface or QC populations cannot close merely because one clustering or mapping job completed.
