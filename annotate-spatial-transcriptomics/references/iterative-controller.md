# Iterative annotation controller

The first clustering is a proposal generator, not the final annotation. New projects use the direct-lineage cohort architecture in `direct-lineage-controller.md`; persistent biological pools are not part of the active state machine.

## Round 0: context and input

Validate biological context, discover and hash inputs, freeze the analysis set, verify the full-feature evidence layer and resume existing completed cohorts without recomputing them.

## Round 1: whole-tissue broad annotation

Select a broad clustering adaptively. Before paper context or label interpretation, freeze label-blind positive/anti DEG, all eligible broad-lineage programs, winner/runner-up margin, contradictions and technical flags. Only then perform open-world review and write a supported `initial_broad_label`. Send low-information, contradicted or irreducibly mixed observations to `qc_holdout`.

## Round 2: broad-class cohorts

For every supported initial broad class, freeze one `broad_class_recluster` cohort, run a cohort-specific candidate resolution grid and adjudicate every subcluster. Outcomes are mutually exclusive: high-confidence fine label, direct parent-broad return, direct cross-lineage return, one decision-relevant targeted cohort, or QC/technical retention. A small class may remain broad-only after a recorded underpowered skip.

## Round 3: targeted questions

Use a temporary targeted cohort only for an interpretable local mixture, contamination boundary or context-gated identity. RCTD/reference assistance is optional and lower priority. Canonical high plus independent evidence may support fine; moderate supports broad-only; low enters the final QC holdout.

## Round 4: residual QC and closure

After every broad and targeted cohort is terminal, freeze residual QC. Do not recluster it. Run one calibrated broad-only mapping over the complete analysis set. Directly return only multichannel-supported moderate-or-higher QC observations; compare every defined broad label with the same Atlas result. Reopen only material cluster/cohort disagreement or coherent OOD once, using independent query evidence. Materialize one final annotation and run both Atlas-concordance and direct-workflow validators.

Run `plan_next_iteration.py` after every atomic writeback. A non-empty queue blocks completion. Automated reference predictions may directly fill an empty frozen-QC broad label under the calibrated multichannel gate; they never overwrite a defined label or write fine labels directly.
