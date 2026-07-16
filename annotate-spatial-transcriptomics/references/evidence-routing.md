# Evidence and unresolved-observation routing

This file is the active contract. Retired routing instructions live under `references/legacy/` and cannot control a new project.

## Evidence floor

Use coherent marker families, anti-marker results, full-feature DEG/detection, adjacent-resolution migration, spatial morphology and source/QC composition. A path, status flag or empty table is not biological evidence. Cohort, direct-return and per-label support artifacts must pass their content schemas.

## Broad-class reclustering cohort

Every credible initial broad class receives exactly one immutable `broad_class_recluster` cohort, except a formally documented underpowered skip. Run the complete active workflow grid and audit purity, hidden lineages and defensible subtypes.

Use `question_mode=broad_purity_audit`. Its expected minimum compartment count is one. A cohort is not penalized for remaining one stable biological population. If the full grid shows no stable mixture or true subtype, close it successfully as `homogeneous_parent_confirmed` and return every observation to the parent broad class.

Select the resolution with the best integrated evidence for the current biological question. Preserve stable lineages or true subtypes first; avoid state-only and technical fragmentation second. Use lower complexity only when candidates are otherwise equivalent.

## Subcluster outcomes and direct returns

Every selected-resolution subcluster has exactly one outcome:

- direct parent-broad return;
- high-confidence, fine-anchor-eligible fine return;
- direct cross-lineage broad/fine return;
- one `targeted_recluster` successor for a predeclared mixture question; or
- QC successor.

These memberships must be mutually exclusive and exactly partition the cohort query. A cross-lineage return does not create an intermediate cohort and does not automatically rerun the target broad cohort.

## Targeted mixtures and RCTD

Use `question_mode=targeted_mixture` only for an interpretable local mixture with two or more predeclared competing hypotheses. Its expected compartment count may equal the number of hypotheses. RCTD/reference assistance remains lower priority: canonical high plus independent current-query evidence may support fine; moderate may support broad-only; low goes to QC. RCTD never turns a reject into an anatomical interface by itself.

Reference/query overlap must be zero. Calibration uses disjoint, query-like held-out current-query anchors; external-reference self-classification is diagnostic only.

## Terminal residual QC consensus

After every broad-class and targeted cohort is terminal, recompute the residual `qc_holdout` membership from the current cell ledger and freeze that exact set once.

**Do not recluster terminal residual QC.** Do not create a QC cohort, anchor-recluster successor or other intermediate membership. Run calibrated Atlas, internal-anchor, current-query marker/anti-marker and observed-density/spatial consensus only on the frozen exact set.

Atlas output tiers are `high`, `moderate_only` and `low_reject`. High and moderate_only consensus may return broad-only with `fine_anchor_eligible=false`; low_reject remains QC. Accepted and rejected memberships must be disjoint and exactly partition the Atlas query, and every result must match the cell-ledger writeback cell for cell.

For sheep ovary without a matched count-level reference, GSE233801 remains the primary public somatic channel only at this terminal step. It cannot serve as its own reference in a held-out benchmark.

## Retained uncertainty

Retain `qc_holdout`, `qc_reject` or `technical_state` when evidence is insufficient. A large unresolved fraction triggers evidence review, not a quota and not QC reclustering. Annotation rate is never optimized at the expense of evidence.

## Contamination and ambient signal

Require lineage-core programs and compatible anti-markers. Spatial adjacency, zona signal or a single marker cannot override a coherent resident program. Report cellbin/spot counts separately from biological cell counts.
