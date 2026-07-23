# State model

## v2 project contract

Every new project first freezes `config/annotation_contract.json`. It binds the project/sample, input snapshot, workflow and biological profiles, candidate catalog, query-expression ancestry, whole-tissue partition/grid, separate query-reclustering grid, broad-family evidence, authoritative Atlas route and release taxonomy. `provenance/annotation_contract_validation.json` must match its SHA256. Modifying a bound input/profile/catalog invalidates downstream evidence and completion.

Before broad naming, v2 also requires `provenance/broad_family_evidence_validation.json` over the complete current-project full-feature `cluster × candidate × positive-family` matrix. The matrix stores absolute detection/prevalence and pseudobulk-like support; centered scores and one-vs-rest DEG remain comparative evidence only.

Skill release 2.1.0 keeps the project framework schema at 2.0.0 and adds one shared machine-derived lineage decision table for both initial clusters and selected-resolution cohort subclusters. Required columns are `cluster`, `candidate_id`, `candidate_broad_lineage`, `program_score`, `positive_family_count`, `positive_families`, `anti_program_burden`, `contradiction_count` and `eligible`; cohort tables additionally include `purity_status`, `lineage_supported_fraction` and `strongest_competing_fraction`. The candidate universe must be identical for every cluster. `contradiction_count` means unresolved **material** contradiction; a coherent weaker alternative belongs in the signal ledger as `watch` and does not by itself make the dominant resident parent ineligible. Validators recompute winner, runner-up and margin from this table. Agent prose, paper markers and cohort provenance cannot override it.

For mixed-cluster returns, one `subcluster_id` may legitimately appear in multiple disjoint `subcluster_outcomes` rows. Each biological `supported_subset` row binds its exact membership and `observation_purity_evidence`; all rows for that source cluster plus targeted/QC successors must exactly partition the selected-resolution membership. Never collapse these rows into one whole-cluster return. Present-lineage completeness evidence also stores `source_writeback_audits` and the complete parent-specific `fine_candidate_audit`; see `observation-subset-writeback.md`.

## Cell-level fields

Required fields are `sample_id`, `cell_id`, `decision_id`, `analysis_scope`, `input_snapshot_id`, `source_method`, `source_run_id`, `source_cluster`, `source_key`, `parent_decision_id`, `generation`, `route`, `route_attempt_id`, `state_tags`, `spatial_tags`, `qc_tags`, `candidate_lineages`, `state`, `initial_broad_label`, `broad_label`, `fine_label`, `confidence`, `evidence_status`, `validation_feature_scope`, `recluster_cohort_id`, `assignment_mode`, `cross_lineage_target`, `fine_anchor_eligible`, `decision_version`, `supersedes`, `closed` and `next_action`. Initial cluster rows additionally bind `prelabel_evidence_artifact`, `prelabel_evidence_sha256`, `prelabel_evidence_frozen`, `prelabel_winner`, `prelabel_runner_up` and `prelabel_winning_margin`. Every current cell/bin `decision_id` points to an active, nonsuperseded decision.

Final ledgers retain one `final_state`, `final_broad_label`, `final_fine_label`, `final_broad_confidence`, `final_fine_confidence`, compatibility field `final_confidence`, `final_assignment_tier`, `final_broad_eligible` and `final_fine_eligible`. Moderate-or-higher supports final broad assignment; only high confidence plus `fine_anchor_eligible=true` supports a final fine label. Initial QC exclusion, post-clustering QC holdout, interface retention and calibrated broad rescue remain distinct states. Legacy strict/inclusive/display columns may be preserved for old projects but are not release views.

Freeze the exact `analysis_set` membership in an immutable file and register its path, SHA256, cell-ID column, full-object count, analysis-set count and excluded-initial-QC count in `provenance/analysis_scope_policy.json`. `validate_state.py` treats this policy as authoritative when present. Observations outside that membership remain in the full-object ledger, but their single final state must be `excluded_initial_qc` with empty biological labels. They cannot enter DEG, marker discovery, anchors, rescue routes or biological proportions.

Legacy projects may use `build_annotation_views.py` only for migration diagnostics. New releases use `build_final_annotation.py`; after a multi-evidence atomic writeback it creates the sole `final` registry entry without overwriting the evidence ledger.

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
- `unknown_candidate`
- `excluded_initial_qc`
- `closed_and_frozen`

## Registries

- `input_snapshot_registry.tsv`: immutable input files, dimensions and hashes.
- `derived_expression_registry.tsv`: every normalized/SCT/full-feature query-evidence object with project/sample IDs, raw-count and analysis-set hashes, parent artifact, purpose and explicit external-reference status. Cross-project derived expression is forbidden in query-evidence channels.
- `clustering_decision_ledger.tsv`: every candidate and selection/rejection rationale.
- `cluster_decision_ledger.tsv`: cluster evidence, decision, confidence and next route.
- `recluster_cohort_registry.tsv`: immutable broad-class or targeted cohort membership, source boundary, candidate grid, selected resolution, run, terminal status and outcome artifact.
- `direct_return_registry.tsv`: direct broad/fine/cross-lineage writeback with source cohort/subcluster, target, confidence, evidence, membership hash and anchor eligibility.
- `pool_registry.tsv`: legacy migration registry only; new projects do not create persistent biological pools.
- `run_registry.tsv`: scripts, parameters, environment, scheduler-visible stage-readable name/ID and outputs plus unique active `work_key`, execution fingerprint, worker owner, attempt and `supersedes_run_id`. Submitted/running jobs require a name generated by `scheduler_job_name.py`. Run `scripts/migrate_project_v1_3_to_v1_4.py` before resuming an older registry.
- `pool_snapshot_registry.tsv`: legacy migration registry only.
- `route_attempt_registry.tsv`: RCTD/Atlas assistance. New projects use one terminal `global_atlas_broad_audit` row with `source_state=analysis_set_after_terminal_qc_freeze`; its query equals the analysis set and its `qc_membership` equals terminal residual QC. It binds all-cell concordance, cluster summaries, discrepancy decisions, calibration/index manifests and OOD counts. QC accepted/rejected memberships partition frozen QC; defined-cell outcomes partition into concordant, weak challenge, material disagreement, OOD and ontology conflict. Legacy residual-QC Atlas rows remain migration-only.
- `config/annotation_contract.json`: v2 immutable boundary contract; it is the only authority for the upstream whole-tissue grid and bound candidate catalog.
- `provenance/broad_family_evidence_validation.json`: complete full-feature broad-family matrix validation and hashes.
- `provenance/residual_qc_audit_validation.json`: typed residual-QC reason census plus evidence artifacts; large unresolved fractions cannot close with a free-text rationale.
- `provenance/project_results_audit.json`: read-only release-surface audit of the committed ledger, run/incident state, taxonomy and result-directory consistency.
- `provenance/incidents/incident_registry.tsv`: every EXIT/cancel/fail-closed/runtime/resource incident, failure boundary, reusable artifacts, repair and verification. Open rows block completion.
- `branch_control_board.tsv`: legacy migration registry only; new projects use cohort outcomes and structured controller gaps.
- `annotation_support_registry.tsv`: one validated evidence summary for every released broad class and high-confidence fine label, including positive markers, anti-markers, resolution stability, spatial evidence, route ancestry and validation artifacts. It is frozen into the pre-confirmation HTML and confirmation hash.
- `lineage_signal_boundary_registry.tsv`: every whole-tissue, broad-class and targeted boundary, its selected/candidate/audited resolutions, exact cluster universe, complete candidate catalog, unexplained-program audit and large-label purity audit.
- `lineage_signal_registry.tsv`: one row per boundary/resolution/cluster/catalog lineage plus additional unexplained programs. It preserves `watch` signals across boundaries and records explicit supported or refuted closure evidence.
- `broad_class_completeness_registry.tsv`: post-Atlas, query-derived review of every present broad class and every zero-census default tissue lineage, including full-membership purity, selected-plus-two-higher, large-label embedded programs, spatial morphology, QC/OOD and technical-missingness closure.
- `workflow_event_registry.tsv`: chronological Chinese-ready input, job, failure, repair, review and writeback events.
- `annotation_view_registry.tsv`: the validated single `final` annotation census and membership artifact. Legacy view rows may remain as non-release history.
- `provenance/analysis_scope_policy.json`: immutable full-object versus analysis-set boundary and membership hash.
- `provenance/release_taxonomy_audit.json`: profile-bound audit of final biological names versus source cohorts and retained states, including separate biological and retained-state censuses plus cell-ledger/profile hashes.
- `provenance/prelabel_broad_evidence_validation.json`: exact active-cluster coverage by label-blind evidence freezes and current hashes.
- `provenance/annotation_membership_partition_audit.json`: exact closure of initial broad cohorts, subcluster outcomes, direct returns, targeted successors, all-cell Atlas query, frozen-QC accepted/rejected outcomes and the final analysis set.
- `provenance/annotation_support_validation.json`: exact equality between final broad/fine labels and schema-validated support records.
- `provenance/whole_tissue_resolution_grid_validation.json`: hash-bound active-workflow audit of the complete whole-tissue candidate grid; mandatory for a verified same-batch StereoPy cellbin contract.
- `state/reference_registry.tsv`: immutable matched/public reference records with role, species/tissue/stage match, object/annotation/marker paths and hashes, source-label column, donor/sample composition and usability status.
- `config/matched_reference_crosswalk.tsv` plus `provenance/matched_reference_crosswalk_validation.json`: source labels, candidate review axes, candidate release names, transfer ceilings, reference markers and required current-query evidence. Source labels are provenance and are never overwritten.
- `cell_ledger.tsv.gz`: final cell/bin-level provenance and labels.

## Invariants

Cell IDs are strings and unique per sample. Frozen cohorts cannot be recomputed in the same decision version. Broad-only rescues have `fine_anchor_eligible=false`. A fine label requires a nonempty catalog-bound broad parent, high fine confidence, `fine_anchor_eligible=true`, `final_fine_eligible=true`, and `state=final_state=defined_fine`; an empty fine label cannot carry `defined_fine`. Every non-initial label must point to a source run and evidence status. Fine writeback cannot change a locked broad label. Version changes append records or create a new ledger; they do not erase previous decisions.

Every closed biological or interface decision requires `closure_rationale` and a resolvable `validation_artifact`. Every terminal context-specific identity requires `validation_status=passed` plus its validation artifact. Every supported initial broad class has one terminal broad-class cohort unless an underpowered skip is recorded. Atlas receives the complete analysis set once after residual QC is frozen; direct writeback is restricted to that QC submembership, while defined labels require a closed discrepancy decision before any change. Every interface requires a deconvolution applicability record only when RCTD is attempted. Scheduler failures use terminal status `failed_preserved`; a repaired execution receives a new run ID. Cohort IDs are versioned and membership hashes are immutable.

RCTD is low-priority assistance. Its canonical `high`, `moderate` and `low` tier counts must partition the complete query. `rctd_fine_return_n` cannot exceed the high tier and requires independent marker/anti-marker, reclustering and spatial evidence. RCTD evidence alone does not make a fine anchor. `rctd_fine_return_n + rctd_broad_return_n` cannot exceed high plus moderate confidence. Fine, broad and rerouted outcomes partition the query; high/moderate observations failing an independent/context gate are rerouted rather than forced. Every low observation enters the final QC holdout and may reach calibrated Atlas only after all broad/targeted decisions are terminal.

No-repeat policy is scoped to a frozen cohort generation. A cell can enter a later scientific question only through an explicit superseding version; it cannot silently re-enter a completed cohort in the same annotation version. Direct cross-lineage return does not itself trigger another target-lineage reclustering.

Cohort identity is provenance, not biology. A cohort name cannot be copied directly into a biological label. Legacy pool identity follows the same rule. Literature categories that fail query-specific gates are recorded as negative audits, not inserted into the cell ledger. Biological broad labels and retained anatomical-interface/QC/technical/pending states have separate release censuses and evidence cohorts.

A label threshold controls naming, not memory. Positive-family evidence cannot be stored as `absent`. Every `watch`, `candidate` or `supported` lineage signal remains open until it has a supported biological outcome or hash-bound multichannel refutation. Parent/cohort identity never narrows the catalog scanned at a later boundary.

The canonical Oocyte recall cohort and an Oocyte spatial context window are different memberships. The latter is evidence-only and cannot expand the Oocyte census. It may support direct somatic returns or a candidate review.
