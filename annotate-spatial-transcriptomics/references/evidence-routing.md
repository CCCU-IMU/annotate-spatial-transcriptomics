# Evidence and unresolved-observation routing

This file is the active contract. Retired routing instructions live under `references/legacy/` and cannot control a new project.

## Evidence floor

Use coherent marker families, anti-marker results, full-feature DEG/detection, adjacent-resolution migration, spatial morphology and source/QC composition. A path, status flag or empty table is not biological evidence. Cohort, direct-return and per-label support artifacts must pass their content schemas.

Before initial broad interpretation, freeze a label-blind evidence matrix covering every eligible broad lineage. Record positive and anti-DEG artifacts, multi-family program scores, winner/runner-up margin, contradictions and technical alternatives. Do not load paper labels or select a convenient marker subset until this artifact is immutable.

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

## Terminal all-cell Atlas concordance

After every broad-class and targeted cohort is terminal, recompute residual `qc_holdout` and freeze it once. Then map the complete analysis set to one calibrated broad-only Atlas.

**Do not recluster terminal residual QC.** Use a fixed reusable Atlas transform/index when available. Never compute a dense query-by-reference distance matrix or repeat joint Atlas integration per sample.

Atlas output tiers are `high`, `moderate_only` and `low_reject`, with explicit OOD and neighbor-mixture diagnostics. For an unlabeled frozen-QC observation, high/moderate_only returns broad-only unless OOD, ontology-conflicted or excluded by the profile's Atlas scope. Current-query marker/anti-marker and spatial evidence audit the return and may reopen a coherent contradictory group, but are not required again per observation after the tier has been calibrated on held-out query anchors. Low/OOD/ontology-conflicted results remain QC. For defined cells, agreement closes, weak disagreement is logged, and material calibrated disagreement reopens the entire cluster/cohort once. Atlas alone never overwrites a defined label.

Use default material triggers of at least 30 observations and 20% of a cluster unless a profile supplies stricter thresholds. A coherent OOD group receives a rare/unknown exception rather than being forced to the nearest known class. Orthogonal review compares both hypotheses using query full-feature programs, anti-markers, pseudobulk, sample consistency, morphology and technical alternatives. Bind each decision to `atlas_discrepancy_decision.schema.json`; close as retain, supersede, downgrade, QC or `unknown_candidate`.

For sheep ovary without a matched count-level reference, GSE233801 remains the primary public somatic channel only at this terminal step. It cannot serve as its own reference in a held-out benchmark.

## Retained uncertainty

Retain `qc_holdout`, `qc_reject`, `unknown_candidate` or `technical_state` when evidence is insufficient or outside reference support. A large unresolved fraction triggers evidence review, not a quota and not QC reclustering. Annotation rate is never optimized at the expense of evidence.

## Contamination and ambient signal

Require lineage-core programs and compatible anti-markers. Spatial adjacency, zona signal or a single marker cannot override a coherent resident program. Report cellbin/spot counts separately from biological cell counts.
