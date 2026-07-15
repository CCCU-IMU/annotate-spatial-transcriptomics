# Multi-route unresolved-population controller

The controller routes by failure mode, not by the name of the source graph cluster. Preserve every source and state tag before routing.

## Route A: biological parent-pool anchor reclustering

Use for a large population with usable counts and one or more plausible lineages.

1. Freeze query membership and balanced, trustworthy reference-only anchors separately.
2. Membership must contain `query_or_anchor`, `anchor_label`, `source_key`, `state_tags`, `spatial_tags` and `qc_tags`.
3. Fit normalization and PCA jointly. Build neighbors, clusters, UMAP, DEG and cluster statistics on query observations only. Anchors never enter query DEG or final pool counts.
4. Test the profile-declared pool grid. Sheep ovary always uses `0.1,0.2,0.3,0.4,0.6`; review adjacent-resolution migration/ARI, marker and anti-marker completeness, source/state/QC composition, PCA distance and spatial morphology.
5. Each subcluster becomes one of: fine definition, broad-only return, canonical-pool reroute, interface, QC/technical retention or a context-specific candidate needing ordinary lineage validation.
6. Roll back source-, depth-, ECM-, stress-, cell-cycle- or ribosomal-driven splits.

An analysis that subsets only query observations is not anchor-assisted even if its directory or route name says so.

## Route B: local mixed interface

Use for a spatially local narrow-layer or mixed-lineage population after Route A fails to produce a stable separation.

1. Run targeted anchor reclustering first.
2. Record deconvolution applicability from query size, UMI/depth, reference types, gene overlap, interface morphology and presence of separable held-out anchors.
3. If applicable, depth-match held-out anchors to the query and compare parameter settings. Calibrate top weight/probability and top-vs-second margin to a predeclared precision target. Record coverage, precision and recall.
4. Apply tiered evidence. Extreme-confidence RCTD plus independent full-feature marker/anti-marker, stable reclustering and spatial evidence may support a fine return. High-confidence RCTD supports broad-only return with `fine_anchor_eligible=false`.
5. Send every medium/low-confidence RCTD cell to the frozen post-clustering QC holdout. Do **not** call Atlas from Route B. These cells first join the complete QC-holdout membership used by Route C and must undergo its broad-anchor reclustering.
6. A local interface must be spatially restricted and small relative to the analysis set. A large or diffuse reject population must be reopened as a canonical broad parent pool or sent to the same frozen QC holdout; the existence of an RCTD result does not make it a valid terminal interface.

`not_applicable` requires a validation artifact and concrete failed criteria; a narrative statement alone cannot close the branch.

## Route C: post-clustering low-information/QC holdout

Keep initial QC exclusion separate. For a nontrivial post-clustering holdout:

1. Freeze the complete QC query and its source components.
2. Run a full QC-pool Route-A-style anchor reclustering before atlas mapping. Use broad anchors, query-only graph/DEG and the profile-declared formal grid. Sheep ovary uses the complete `0.1,0.2,0.3,0.4,0.6` grid here as well; never substitute an ultra-low grid.
3. Preserve technical components and singleton/microclusters as QC even if they land near a large biological cluster.
4. Freeze the cells still labelled QC after that reclustering as a new residual-QC child snapshot. Only this exact residual membership may enter Atlas/internal-anchor review. The Atlas attempt must link the preceding QC-anchor attempt, its outcome hash and the residual-QC membership hash; an already defined broad/fine cell or an ordinary biological-pool cell is out of scope.
5. For that residual child snapshot, audit reference adequacy and depth-match held-out internal anchors. Combine external reference, internal anchors, marker programs and spatial support.
6. Compute an observed-density spatial prior from all eligible anchors. A class-balanced spatial prior may be used only as a sensitivity analysis and must not replace observed tissue density.
7. Rescue calibrated moderate-or-higher broad-only labels into the single final broad annotation. Rejects remain QC. All rescues are `fine_anchor_eligible=false` and can never receive a fine label from reference transfer alone.

Atlas mapping alone is not a completed Route C for a large QC pool.

## Context-specific candidate validation (no generic rare-cell route)

Rarity alone never defines a cell type or triggers a mandatory branch. The whole-tissue open-world audit must examine unexplained programs and literature-supported Oocyte, luteal, Smooth muscle and Neural/Schwann boundaries. Supported non-Oocyte candidates form sample-specific ordinary biological pools and use Route A or targeted validation. Unsupported candidates receive negative audits.

Oocyte is the exception because zona/adjacent contamination is a known identity confounder. When an Oocyte candidate exists, use full-feature multi-module positive evidence, contradictory somatic programs, spatial objects/morphology, adjacent/background comparison and complete-candidate reclustering. Return contamination/adjacent cells to canonical somatic pools with explicit ambient/adjacent tags. The registry class is `context_specific_identity_review`; legacy `strict_rare_cell_review` records remain readable but should not be created in new projects.

## Route E: diagnostic supervised state review

Use only after the appropriate unsupervised route for ambient-high, LOC-dominated, baseline-low-RNA, vascular–stromal transition or graph-outlier failures. Compare full and nuisance-excluded models or depth-matched anchors. It may support broad return or explicit state retention; it must not manufacture a fine subtype.

## Required route registry outcomes

Every applicable or reviewed-inapplicable branch records query/anchor counts, membership artifact and hash, selected parameters/resolution, validation artifact, outcome counts, retained state, fine-anchor eligibility and no-repeat policy in `route_attempt_registry.tsv` and `branch_control_board.tsv`.

For Route B additionally record the extreme/high/medium-low confidence counts, independent evidence gates, the QC-holdout reroute membership, its successor QC-anchor route ID and final reject count. For Route C, record the full QC input snapshot separately from the residual-QC child snapshot; the Atlas query hash must equal the QC-anchor route's frozen residual-QC hash. Do not merge confidence tiers or silently broaden the Atlas query.

## Closure and confirmation order

A branch closes only when its entire frozen membership is mutually exclusively accounted for and every applicable route has a validated outcome. The project must machine-audit whole tissue, open-world lineage discovery, Route A, Route B, Route C, an Oocyte-specific route when applicable, and the single final annotation in that order. A zero-query or not-applicable phase is recorded explicitly. Main-Agent quality approval and user confirmation are requested only after this phase audit, correct profile/preset binding and zero-open-incident gate pass.
