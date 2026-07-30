# Existing-sample non-destructive resume after the execution-chain repair

This procedure is for a v2.5 project that already completed expensive clustering, broad merge, Atlas routing, unresolved review or follicle-ROI repair. The objective is to enter the repaired strictly serial per-cell-type review without repeating validated computation.

## Decide the recovery boundary

Use the latest artifact that passed biological and membership validation:

1. If a canonical controller pause already contains `membership_transform_chain`, `cell_type_review_state` and `pause_reason=single_active_cell_type_review|single_cell_type_review_step_applied`, resume it directly.
2. If the project predates the generic transform chain but has exact source/result membership and canonical stage manifests, reconstruct the chain in chronological order, validate it, and create the next serial-review pause from the same frozen membership.
3. If a transition lacks exact source/result membership, exact cell IDs or a canonical evidence manifest, do not guess it. Return only to that transition's source stage; do not rerun whole-tissue or second-round clustering.
4. If clustering input, cohort membership, raw-count ancestry, grid or seed changed, its partition cache is not reusable.

## GSE233801 split-wall Atlas migration

The deprecated `sheep_ovary_GSE233801_v1` merged endothelial and mural/contractile components. It remains valid only when resuming a review whose Atlas result and membership transition were already frozen. Never reinterpret an old `Vascular/endothelial` prediction as Endothelial, or an old `Stromal/perivascular` prediction as Stromal, Pericyte/mural or Smooth muscle.

When an existing sample must adopt the split taxonomy, keep its input, both clustering rounds, local splits and frozen pre-Atlas broad membership. Re-run only the fixed-feature projection and Atlas routing with `sheep_ovary_GSE233801_split_wall_v2`, then resume the serial cell-type review. Reopen Endothelial, Pericyte/mural, Smooth muscle and any membership actually changed by the new Atlas stage; do not repeat whole-tissue or cohort clustering unless a source subcluster is independently shown to be non-evaluable.

## Direct resume from a repaired pause

Run the same `atlas_and_completeness_review` command that created the pause and add:

```bash
--resume-review-manifest /absolute/path/to/lineage_controller_manifest.json
```

The required contract, prerequisite and frozen broad arguments must bind the same project, but Atlas mapping, unresolved rescue, ROI repair and clustering are not executed again. Supply at most one new `--lineage-review-decisions` file, and it must target the active review ID shown in the pause manifest.

The controller output gives exactly one required user-facing message:

```text
现在开始对 <cell type> 进行专项复核。
```

After that type is retained or patched, resume the newly written pause to activate the next type.

## Reconstruct an older project's transform chain

Initialize from the exact frozen broad membership:

```bash
python annotate-spatial-transcriptomics/scripts/manage_membership_transform_chain.py init \
  --membership /absolute/path/to/frozen_broad_membership.tsv.gz \
  --out /absolute/path/to/recovery/T0000
```

Append every actual membership transition in chronological order. Example operation names are:

```text
atlas_unlabeled_broad_rescue
post_merge_unresolved_return
follicle_roi_assignment
follicle_roi_depublication
follicle_roi_reconciliation
source_unit_sync
cell_type_review_patch
final_release_materialization
```

For each transition:

```bash
python annotate-spatial-transcriptomics/scripts/manage_membership_transform_chain.py append \
  --chain /absolute/path/to/previous/membership_transform_chain.json \
  --operation <canonical_operation> \
  --source /absolute/path/to/source_membership.tsv.gz \
  --result /absolute/path/to/result_membership.tsv.gz \
  --evidence-manifest /absolute/path/to/canonical_stage_manifest.json \
  --out /absolute/path/to/recovery/T0001
```

For a `cell_type_review_patch`, also supply `--target-cell-type '<broad label>'`. Validate the terminal chain against the current membership:

```bash
python annotate-spatial-transcriptomics/scripts/manage_membership_transform_chain.py validate \
  --chain /absolute/path/to/terminal/membership_transform_chain.json \
  --current-membership /absolute/path/to/current_membership.tsv.gz \
  --out /absolute/path/to/recovery/validation
```

The validator rejects a changed cell universe, discontinuous source/result paths, non-exact deltas, identity changes during source synchronization, Atlas overwrites of defined labels and evidence manifests that do not bind their result.

## What is reused

- input snapshot and raw-count ancestry;
- whole-tissue partition;
- every validated second-round SCT/PCA/SNN/Leiden partition;
- local mixed-subcluster results;
- broad freeze;
- the one fixed-Atlas mapping and its decisions, except when explicitly migrating a merged v1 project to the split-wall v2 taxonomy;
- bounded unresolved and follicle-ROI results whose transform provenance validates;
- shared raw-count, coordinate and pseudobulk caches.

## What is recomputed

- the current active cell type's target-versus-competitor evidence packet;
- its current-member precision and whole-query recall questions;
- its molecular, spatial and literature-boundary conclusion;
- exact membership deltas if a patch is approved;
- signatures for cell types touched by that patch;
- final completeness/report assets after all types close.

An unaffected, already closed cell type remains closed. The recovered project must never be used as an Atlas or as expected truth for another sample.
