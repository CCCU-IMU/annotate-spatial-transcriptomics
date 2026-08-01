# v2.2 state and artifact model

## Contract boundary

`config/annotation_contract.json` binds the controller and dependency hashes, workflow and biological profiles, complete candidate catalog, input snapshot, exact `full_object`, `analysis_set` and `excluded_initial_qc`, project-local raw-count ancestry, whole-tissue and cohort grids, master seed, Atlas route and artifact roles.

Permitted artifact roles are:

- `runtime_input`: query material allowed to drive scoring or membership;
- `external_reference`: explicitly declared Atlas/reference material;
- `experimental`: non-release analysis output;
- `failed_diagnostic`: retained failure evidence that may provide only aggregate diagnostics;
- `release_output`: formal output produced by an authorized phase.

`failed_diagnostic`, historical labels, repair membership and same-batch annotations cannot enter runtime membership or blind-regression truth.

## Authority model

| Phase | Membership authority |
|---|---|
| `whole_tissue_partition` | provisional initial-cluster partition only |
| `cluster_cohort_recluster` | second-round proposals only |
| `local_mixed_subcluster_split` | local broad candidates only |
| `merge_and_freeze_broad` | formal broad freeze |
| `atlas_and_completeness_review` | unlabeled broad rescue only |
| `materialize_final_release` | final broad/fine/state/QC release |

The first phase never writes biological labels to the cell ledger. `unresolved_biological` is a biological pending state until final materialization; it is not an early synonym for QC.

## Core artifacts

- `whole_tissue_cohort_plan.tsv`: exact initial cluster membership, one cohort per cluster, provisional record and watch programs; never formal labels.
- `raw_count_ancestry.json`: project/sample, non-SCT raw-count assay, parent object and analysis membership hashes for a cohort.
- `cohort_outcome.json`: selected/neighbor resolutions, full-catalog evidence, score-blindness declaration, source initial cluster, whole-subcluster/local/unresolved outcome and optional fine/state/unmodeled proposals.
- `local_split_artifact`: trigger, exact source subcluster, independent candidate components, overlap audit and exact local remainder.
- `frozen_broad_membership.tsv.gz`: first formal biological membership, created only after all cohort outcomes merge.
- `post_atlas_broad_membership.tsv.gz`: frozen broad plus permitted unlabeled Atlas rescue; existing labels remain unchanged.
- `membership_transform_chain.json`: ordered, contiguous membership transforms after broad freeze. Every entry binds exact source/result semantic hashes, the unchanged cell universe, exact delta cell IDs and its canonical evidence manifest. Supported operations include Atlas rescue, unresolved return, ROI assignment/depublication/reconciliation, identity-neutral source synchronization, one-active-cell-type patches and final broad/fine/state/QC plus `final_cell_type` materialization.
- `cell_type_review_state.json`: serial broad-review state. It contains at most one `active_cell_type_review`; every other type is queued, atomically closed after its single full-query decision or explicitly blocked because that decision could not be validated. One state transition can validate and apply only the active type. Closed types remain legal exact transfer destinations/sources for later active reviews but are never requeued.
- `current_stage.json` and `next_action_manifest.json`: scheduler-facing controller state. `REVIEW_REQUIRED` and `ITERATION_REQUIRED` are successful controlled pauses with a bound resume token; `DONE_PENDING_USER_REVIEW` is a completed frozen candidate; `FAILED_RUNTIME` is a genuine software/input/resource failure.
- `final_membership.tsv.gz`: one final broad or typed retained state per analysis observation, optional parent-locked fine and independent state annotations.

All membership artifacts store a deterministic SHA256 over stably sorted observation IDs and semantic columns.

## Cohort outcome fields

Every default cohort binds:

- `source_initial_cluster` and exact membership hash;
- `provisional_broad_after_score_freeze`;
- `raw_count_assay` and ancestry hash;
- complete grid and selected/neighbor resolutions;
- `full_catalog_scan=true` and evidence artifact;
- outcome from `parent_return`, `cross_lineage_return`, `missing_broad_reconstruction`, `fine_candidate`, `state_annotation`, `unmodeled_lineage_candidate`, `unresolved_biological` or `underpowered_not_evaluable`;
- local split requirement/artifact;
- broad and fine release proposals.

Provisional identity is provenance. It cannot narrow candidate scoring or establish a parent return.

## Final cell fields

Final membership includes at minimum:

- `sample_id`, `cell_id` and `analysis_scope`;
- `final_broad_label`, optional `final_fine_label` and independent `state_annotations`;
- `final_state`, confidence/tier and typed unresolved/QC reason;
- source initial cluster, second-round cohort, subcluster, assignment mode and evidence artifact;
- input, membership and semantic hashes.

An excluded-initial-QC observation stays outside the biological analysis set. A fine label requires a nonempty matching frozen broad parent and high-confidence independent evidence. Atlas rescue is broad-only and never fine-anchor eligible.

## Invariants

- Every analysis observation belongs to exactly one initial cluster and one default second-round cohort.
- Every second-round scorer sees the full catalog and cannot see provisional/historical labels before score freeze.
- Local splitting is restricted to the exact triggered second-round subcluster.
- Broad freeze exactly and disjointly covers the analysis set, including unresolved biological members.
- Fine assignment cannot alter broad membership.
- Atlas cannot overwrite an existing broad label.
- The sheep-ovary Atlas identity, capability matrix and asset hashes are immutable contract inputs. Unsupported or `not_evaluable` classes cannot rescue or challenge query labels.
- Formal per-broad review is serial. A multi-type decision file, multiple active packets or batch closure has no release authority.
- Membership-count or scope-signature changes never reopen a closed cell-type specialist review. A later active review may write an exact evidence-bounded cell set into or out of a closed type and records the change as a post-closure delta.
- A controlled review pause can resume from its bound membership transform chain without repeating either clustering round, Atlas routing, unresolved rescue or follicle-ROI repair.
- Technical QC is assigned only to objective input failures or during final typed closure.
- Same input, contract and seed reproduce identical membership and semantic hashes.
- Version changes append or supersede artifacts; they do not erase failure evidence or previous release records.
