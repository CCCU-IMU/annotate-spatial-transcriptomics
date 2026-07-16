# Iterative annotation controller

The first clustering is a proposal generator, not the final annotation. New projects use the direct-lineage cohort architecture in `direct-lineage-controller.md`; persistent biological pools are not part of the active state machine.

## Round 0: context and input

Validate biological context, discover and hash inputs, freeze the analysis set, verify the full-feature evidence layer and resume existing completed cohorts without recomputing them.

## Round 1: whole-tissue broad annotation

Select a broad clustering adaptively. Generate cluster DEG, marker/anti-marker programs, UMAP and spatial maps. Perform open-world lineage review and write a direct `initial_broad_label` for supported clusters. Send only low-information, featureless or irreducibly mixed observations to `qc_holdout`; do not create intermediate biological containers.

## Round 2: broad-class cohorts

For every supported initial broad class, freeze one `broad_class_recluster` cohort, run a cohort-specific candidate resolution grid and adjudicate every subcluster. Outcomes are mutually exclusive: high-confidence fine label, direct parent-broad return, direct cross-lineage return, one decision-relevant targeted cohort, or QC/technical retention. A small class may remain broad-only after a recorded underpowered skip.

## Round 3: targeted questions

Use a temporary targeted cohort only for an interpretable local mixture, contamination boundary or context-gated identity. RCTD/reference assistance is optional and lower priority. Canonical high plus independent evidence may support fine; moderate supports broad-only; low enters the final QC holdout.

## Round 4: residual QC and closure

After every broad and targeted cohort is terminal, freeze the complete residual QC membership. Do not recluster it. Apply calibrated Atlas/internal-anchor/marker/spatial consensus to that exact membership; moderate-or-higher returns broad-only and low remains QC reject. Materialize one final annotation and run `validate_direct_lineage_workflow.py`.

Run `plan_next_iteration.py` after every atomic writeback. A non-empty queue blocks completion. Automated scores and reference predictions support decisions but never write labels directly.
