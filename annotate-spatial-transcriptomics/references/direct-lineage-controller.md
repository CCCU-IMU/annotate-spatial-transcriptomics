# Direct-lineage annotation controller

This is the default architecture for every new project. It replaces persistent biological review pools and the old Route A/B/C controller. Legacy pool registries remain readable only for migration.

## Ordered state machine

1. **Whole-tissue proposal.** Select a defensible broad resolution from the current sample. Generate full-feature DEG, marker/anti-marker summaries, UMAP and spatial evidence for every cluster.
2. **Initial broad decisions.** Review each cluster against an open-world lineage catalog. Assign a moderate-or-higher `initial_broad_label` directly when supported. Low-information, featureless or irreducibly mixed observations enter `qc_holdout`. Do not create a biological pool.
3. **One broad cohort per class.** Freeze one immutable `broad_class_recluster` membership for every initial broad class. Recluster it with a cohort-appropriate grid and choose the lowest stable biological resolution. A genuinely underpowered class may close `not_applicable_reviewed` with an evidence artifact and remain broad-only.
4. **Adjudicate every subcluster.** A coherent high-confidence program may receive one shallow fine label. A subcluster that only confirms its parent returns directly to that broad class. State-only fragmentation is merged into the parent.
5. **Direct cross-lineage return.** If a subcluster coherently supports another lineage, write it directly to that broad class or high-confidence fine label. Preserve source cohort, subcluster, membership hash, evidence and state/spatial tags. Do not create an intermediate target pool and do not automatically recluster it with the target class.
6. **Targeted cohort only when decision-relevant.** Create one immutable `targeted_recluster` cohort only for a local, interpretable mixture, contamination boundary or context-gated identity that cannot be adjudicated from the existing broad cohort. It must state the competing hypotheses and stop after answering that question.
7. **Optional RCTD/reference assistance.** Use only when a local interface has an appropriate reference and machine-readable applicability record. Extreme confidence plus independent marker/anti-marker, resolution and spatial evidence may support a fine return; high confidence supports broad-only; medium/low enters `qc_holdout`.
8. **Terminal residual QC review.** After all broad and targeted cohorts are terminal, freeze the complete residual `qc_holdout` membership once. Do not recluster it. Apply calibrated Atlas/internal-anchor/marker/spatial consensus only to that exact membership. Moderate-or-higher returns broad-only with `fine_anchor_eligible=false`; lower confidence remains QC reject.
9. **Single final annotation.** Every analysis-set observation has one moderate-or-higher final broad label with an optional high-confidence fine label, or an explicit retained QC/technical state.

## Why cohort is not a pool

A cohort is a frozen computational query boundary for one question. It has no biological name, cannot become a report node and cannot collect unrelated unresolved cells over time. A direct return closes the source question. Re-entering a completed cohort requires a versioned successor and an explicit new scientific question.

## Completion contract

`validate_direct_lineage_workflow.py` must verify:

- each initial broad class has a terminal broad cohort or reviewed underpowered skip;
- every cohort and direct return has immutable membership and evidence;
- broad/fine confidence and RCTD ceilings are respected;
- the Atlas query equals the final residual QC membership and Atlas returns are broad-only;
- a single final annotation accounts for the whole analysis set.

Only after this audit, state/taxonomy/completion gates and main-Agent quality approval pass may the lightweight confirmation report be shown to the user. Full DEG/figures/report wait for explicit confirmation.
