# State model

## Cell-level fields

Required fields are `sample_id`, `cell_id`, `decision_id`, `analysis_scope`, `input_snapshot_id`, `source_method`, `source_run_id`, `source_cluster`, `source_key`, `parent_decision_id`, `parent_pool_id`, `pool_snapshot_id`, `generation`, `route`, `route_attempt_id`, `state_tags`, `spatial_tags`, `qc_tags`, `candidate_lineages`, `state`, `broad_label`, `fine_label`, `confidence`, `evidence_status`, `validation_feature_scope`, `fine_anchor_eligible`, `decision_version`, `supersedes`, `closed` and `next_action`. Every current cell/bin `decision_id` must point to an active, nonsuperseded cluster/pool decision.

Final ledgers additionally retain mutually exclusive `strict_state/broad_label/fine_label`, `inclusive_state/broad_label/fine_label`, and `display_state/broad_label/fine_label`, plus `display_policy`. Initial QC exclusion, post-clustering QC holdout, interface retention and calibrated broad rescue are distinct states.

Freeze the exact `analysis_set` membership in an immutable file and register its path, SHA256, cell-ID column, full-object count, analysis-set count and excluded-initial-QC count in `provenance/analysis_scope_policy.json`. `validate_state.py` treats this policy as authoritative when present. Observations outside that membership remain in the full-object ledger, but all three view states must be `excluded_initial_qc` with empty biological labels. They cannot enter DEG, marker discovery, anchors, rescue routes or biological proportions.

Use `build_annotation_views.py --mode derive` only when view fields have not yet been adjudicated. After a multi-evidence atomic writeback, use `--mode preserve --registry-membership-path state/cell_ledger.tsv.gz` to refresh census and hashes without overwriting biological view decisions.

Rare or context-gated identities additionally require `validation_status`, `validation_artifact` and `validation_feature_scope`. Iterative branches record `iteration`, `route_run_id` and `closure_rationale`. Use `validation_feature_scope=full_feature` only when label-level positive/negative evidence was computed from an object that passed the profile feature-scope audit.

For spot/cellbin rare identities, record `spatial_object_count` and `count_interpretation=observations_not_inferred_cells` (or an equally explicit phrase). The labeled observation count and inferred object count are separate quantities; neither automatically equals histological cell count.

Rare-labeled cellbins also retain `spatial_object_id` so their object aggregation can be reconstructed from the cell ledger.

Allowed states:

- `defined_fine`
- `defined_broad_only`
- `interface_review`
- `qc_holdout`
- `technical_state`
- `pending_review`
- `excluded_initial_qc`
- `closed_and_frozen`

## Registries

- `input_snapshot_registry.tsv`: immutable input files, dimensions and hashes.
- `clustering_decision_ledger.tsv`: every candidate and selection/rejection rationale.
- `cluster_decision_ledger.tsv`: cluster evidence, decision, confidence and next route.
- `pool_registry.tsv`: membership snapshot, parent, status and closure.
- `run_registry.tsv`: scripts, parameters, environment, scheduler job and outputs.
- `pool_snapshot_registry.tsv`: immutable versioned memberships and SHA256; parent pools are never silently expanded.
- `route_attempt_registry.tsv`: Route A–E applicability, query/anchor boundary, parameters, validation and outcome counts. Applicable RCTD/reference-assisted routes additionally record `rctd_extreme_n`, `rctd_high_n`, `rctd_medium_low_n`, `rctd_fine_return_n`, `rctd_broad_return_n`, `independent_fine_evidence` and `fallback_route_attempt_id`.
- `branch_control_board.tsv`: branch generation, selected resolution, terminal/no-repeat policy and authoritative artifact.
- `workflow_event_registry.tsv`: chronological Chinese-ready input, job, failure, repair, review and writeback events.
- `annotation_view_registry.tsv`: strict/inclusive/display census and policy artifacts.
- `provenance/analysis_scope_policy.json`: immutable full-object versus analysis-set boundary and membership hash.
- `provenance/release_taxonomy_audit.json`: profile-bound audit of final biological names versus source pools and retained states, including separate biological and retained-state censuses plus cell-ledger/profile hashes.
- `cell_ledger.tsv.gz`: final cell/bin-level provenance and labels.

## Invariants

Cell IDs are strings and unique per sample. Frozen pools cannot be reused in the same decision version. Broad-only rescues have `fine_anchor_eligible=false`. A fine label requires a nonempty broad label. Every non-initial label must point to a source run and evidence status. Version changes append records or create a new ledger; they do not erase previous decisions.

Every closed biological or interface decision requires `closure_rationale` and a resolvable `validation_artifact`. Every terminal rare/context-gated identity requires `validation_status=passed` plus its validation artifact. A non-zero large post-clustering QC holdout requires a validated full QC-pool anchor-recluster attempt before depth-matched atlas/internal-anchor review; verified zero-count observations are the exception. Every interface requires a deconvolution applicability record. Scheduler failures use terminal status `failed_preserved`; a repaired execution receives a new run ID. Pool IDs are versioned and membership hashes are immutable.

RCTD is low-priority assistance. Its tier counts must partition the complete query. `rctd_fine_return_n` cannot exceed the extreme tier and requires independent marker/anti-marker, reclustering and spatial evidence. Fine labels produced by the RCTD route itself are not fine-anchor eligible. `rctd_fine_return_n + rctd_broad_return_n` cannot exceed extreme plus high confidence. Fine, broad and rerouted outcomes partition the query; extreme/high observations failing an independent/context gate are rerouted rather than forced. Every medium/low observation is rerouted to the recorded atlas/internal-anchor fallback; the interface branch cannot close until that fallback has a validated outcome.

No-repeat policy is scoped to a frozen branch generation. A cell can enter a later scientific question only through an explicit superseding version; it cannot silently re-enter an unresolved pool in the same annotation version.

Pool identity is provenance, not biology. A pool name cannot be copied directly into a biological label. Literature categories that fail query-specific gates are recorded as negative audits, not inserted into the cell ledger. Biological broad labels and retained anatomical-interface/QC/technical/pending states have separate release censuses and evidence cohorts.
