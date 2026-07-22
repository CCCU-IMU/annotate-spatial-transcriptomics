# Direct-lineage annotation controller

This is the default architecture for every new project. Active entities are `broad_class_recluster cohort`, `targeted_recluster cohort`, `direct return` and `terminal residual QC`. Retired controller registries remain readable only for migration.

## Ordered state machine

1. **Whole-tissue proposal and prelabel freeze.** Select a defensible broad resolution. For BANKSY, the selection question is explicitly whole-tissue broad-class recall/purity: review every formal resolution for full-catalog recall, zero-census default lineages, embedded programs inside large clusters, marker/DEG coherence, spatial morphology, technical fragmentation and adjacent-resolution migration. Cluster count is never the selection target. Before loading paper labels or a favored interpretation, freeze full-feature positive DEG, anti-DEG, technical flags and simultaneous evidence for every eligible broad lineage. Record winner, runner-up, margin and contradictions with `prelabel_broad_evidence.schema.json`.
2. **Initial broad decisions.** Review the frozen matrix against the open-world catalog. The initial partition establishes supported coarse resident compartments and does not have to expose every final lineage. First validate that every default lineage exposes at least two explicit positive marker families. Assign a moderate-or-higher `initial_broad_label` only when the machine-derived decision table recomputes a positive winner margin, absolute full-feature detection/pseudobulk supports its backbone and support family, and alternatives/anti-programs are cleared. Return an informative mixed cluster to its dominant supported parent and retain coherent subthreshold alternatives as `watch`; only genuinely low-information, featureless or irreducibly mixed observations enter `qc_holdout`.
3. **One broad cohort per class.** Freeze one immutable `broad_class_recluster` membership for every initial broad class. Run its complete active grid and choose the integrated-evidence optimum for the purity/subtype question. A single stable population is a successful `homogeneous_parent_confirmed` endpoint, never an under-split failure. A genuinely underpowered class may close `not_applicable_reviewed` with a hash-bound evidence artifact and remain broad-only.
4. **Adjudicate every subcluster.** A coherent high-confidence program may receive one shallow fine label. A subcluster that numerically confirms its parent returns directly to that broad class. Parent identity alone is not evidence. Except for the canonical Oocyte exception, whole-subcluster writeback requires a machine-derived positive-margin winner, two positive families, contradiction clearance and observation-level purity. A mixed subcluster receives one targeted question or QC rather than wholesale parent inheritance. State-only fragmentation may merge only after that parent-purity decision passes.
   At the selected resolution and the next two higher available candidates, every subcluster is rescored against the complete open-world lineage catalog plus unexplained coherent multi-gene programs. The cohort parent is not a prior and must not narrow the candidates.
5. **Direct cross-lineage return.** If a subcluster coherently supports another lineage, write it directly to that broad class or high-confidence fine label. Preserve source cohort, subcluster, membership hash, evidence and state/spatial tags. Do not create an intermediate cohort and do not automatically recluster it with the target class.
6. **Targeted cohort only when decision-relevant.** Create one immutable `targeted_recluster` cohort only for a local, interpretable mixture, contamination boundary or context-gated identity that cannot be adjudicated from the existing broad cohort. It must state the competing hypotheses and stop after answering that question.
7. **Optional RCTD/reference assistance.** Use only when a local interface has an appropriate reference and machine-readable applicability record. Canonical high confidence plus independent marker/anti-marker, resolution and spatial evidence may support a fine return; moderate confidence supports broad-only; low enters `qc_holdout`.
8. **Terminal all-cell Atlas pass.** After all broad and targeted cohorts are terminal, freeze residual `qc_holdout`, then map the complete analysis set once to a calibrated broad-only Atlas with OOD outputs. Use the same result for QC rescue and defined-label concordance.
9. **State-aware routing.** An unlabeled frozen-QC observation with calibrated moderate-or-higher non-OOD Atlas support and no ontology/profile-scope conflict returns broad-only. Calibration already represents the validated confidence gate; marker/anti-marker and spatial channels audit or challenge coherent groups rather than becoming a duplicate per-cell prerequisite. Defined labels that agree close immediately; weak differences are logged. Material broad disagreement, material ontology conflict after crosswalk inspection or coherent OOD reopens the complete cluster/cohort for one query-evidence review. Atlas alone cannot overwrite a defined label.
10. **Query-derived broad-class completeness review.** After Atlas, audit every present class for complete-membership expression/spatial validity and observation-level purity. Audit every zero-census default tissue lineage using all-cell programs, selected-plus-two-higher scans, embedded large-label components, QC/OOD, morphology and technical missingness. Atlas cannot establish absence. Residual QC at least 10% or 50,000 observations automatically reopens upstream broad-recall review and cannot be closed solely as Atlas low confidence.
11. **Single final annotation.** Every analysis-set observation has one moderate-or-higher final broad label with an optional high-confidence fine label, or an explicit retained QC/technical/unknown state.

## Cohort semantics

A cohort is a frozen computational query boundary for one question. It has no biological name, cannot become a report node and cannot collect unrelated unresolved cells over time. A direct return closes the source question. Re-entering a completed cohort requires a versioned successor and an explicit new scientific question.

## Continuous lineage-signal memory

Whole tissue, every broad-class cohort and every targeted cohort are mandatory open-world scan boundaries. Record the exact cluster universe at the selected resolution plus the next two higher available resolutions, then create one catalog scan row per cluster and candidate lineage. Also record coherent programs not represented in the catalog.

`absent`, `watch`, `candidate`, `supported` and `refuted` are evidence states, not labels. Any positive marker-family evidence makes `absent` invalid. A `watch` signal is deliberately below the naming threshold but remains active across later boundaries. It closes only after support/writeback or explicit multichannel refutation using positive families, anti-programs, cross-resolution persistence, spatial morphology and technical alternatives. Failure to reach an initial broad-label threshold is never itself a refutation. Atlas evidence cannot close the query-derived signal ledger.

Every large boundary (default: at least 10% of the analysis set or 50,000 observations) receives a purity audit even when its aggregate DEG supports the parent. Completion fails if any required boundary/catalog product is missing, any unexplained-program audit is absent, or any `watch`/`candidate`/`supported` signal is unresolved.

## Completion contract

`validate_direct_lineage_workflow.py` must verify:

- each initial broad class has a terminal broad cohort or reviewed underpowered skip;
- every cohort and direct return has immutable membership and evidence;
- broad/fine confidence and RCTD ceilings are respected;
- the global Atlas query equals the analysis set, its frozen QC submembership equals terminal residual QC, and only QC predictions can write broad labels directly;
- all material defined-label disagreements, ontology conflicts and OOD triggers have exactly one orthogonal review decision;
- a single final annotation accounts for the whole analysis set.
- every whole-tissue/broad/targeted boundary passed complete catalog-by-cluster coverage, selected-plus-two-higher cross-resolution review, unexplained-program review and explicit signal closure.

Only after this audit, state/taxonomy/completion gates and main-Agent quality approval pass may the lightweight confirmation report be shown to the user. Full DEG/figures/report wait for explicit confirmation.
