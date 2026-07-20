---
name: annotate-spatial-transcriptomics
description: Independently annotate spatial transcriptomics or single-cell RNA-seq with label-blind evidence-first broad decisions, adaptive clustering, broad/targeted reclustering, direct cross-lineage returns, calibrated all-cell Atlas concordance plus QC broad rescue, explicit unknown/OOD detection, immutable state, shallow high-confidence subtypes and audited reports. Seurat/R-first is preferred for full-feature RDS inputs; AnnData/H5AD, SingleCellExperiment, BANKSY, Scanpy/Leiden and external cluster tables are supported. Use for end-to-end or multi-sample annotation, evidence review, unknown-cell discovery, recovery of unresolved observations and release reporting. Never copy parameters or labels from example projects.
---

# Annotate Spatial Transcriptomics

Use one evidence-first state machine for every new project. A first clustering proposes biological questions; it is never the final annotation. The active architecture uses only broad/targeted cohorts, direct returns and terminal residual QC. Retired registries are migration-only provenance.

Read `references/direct-lineage-controller.md`, `references/iterative-controller.md`, `references/quality-standard.md` and `references/state-schema.md` before analysis.

## Default annotation architecture

Execute this order without skipping phases:

1. Discover, inspect and hash inputs; freeze `full_object`, `analysis_set` and `excluded_initial_qc` memberships.
2. Select a whole-tissue broad clustering adaptively. Before showing paper labels or a favored marker interpretation, freeze label-blind positive DEG, anti-DEG, technical flags and simultaneous scores for every eligible broad lineage, including winner, runner-up, margin and contradictions. Validate with `validate_prelabel_broad_evidence.py`.
3. Review every cluster from that frozen open-world matrix. Assign a moderate-or-higher initial broad label only from at least two coherent marker families, a positive winner margin and no unresolved material contradiction; otherwise use `qc_holdout`. Paper markers may explain a frozen result but may not redefine the candidate set after the fact.
   Start the continuous lineage-signal ledger at this boundary. Scan the selected resolution and next two higher available candidates; `below label threshold` means `watch and carry forward`, not absent.
4. Freeze one immutable `broad_class_recluster` cohort for every supported initial broad class and run the complete active-workflow candidate grid. A one-cluster or unsplit result is not penalized: when no stable mixture or real subtype exists, close successfully as `homogeneous_parent_confirmed` and return all observations to the parent. A genuinely underpowered class may remain broad-only only after a hash-bound `not_applicable_reviewed` skip.
   At the selected resolution and next two higher available candidates, rescan every subcluster against the complete catalog plus unexplained programs. The parent broad label is provenance, never a search-space restriction.
5. Write each subcluster directly as a high-confidence shallow fine label, its supported parent broad class, a supported cross-lineage broad/fine return, one decision-relevant targeted cohort, or QC/technical retention.
   Every `watch`, `candidate` or `supported` signal must be resolved with multichannel evidence; it cannot disappear because it failed an earlier naming gate.
6. Create `targeted_recluster` only for a local interpretable mixture, contamination boundary or context-gated identity. It answers one question and cannot become a long-lived catch-all.
   Its subclusters still receive the same full-catalog, selected-plus-two-higher scan.
7. Use RCTD/reference assistance only when applicable. Canonical high plus independent marker/anti-marker, resolution and spatial evidence may support fine; moderate supports broad-only; low enters QC.
8. After every broad/targeted cohort is terminal, freeze the complete residual `qc_holdout` membership once. Do not recluster it. Run one calibrated broad-only Atlas mapping over the complete `analysis_set`, preferably by projecting into a reusable fixed reference representation/index rather than repeating joint integration.
9. Compare every mapped observation with the frozen primary state using `build_global_atlas_concordance.py`. A QC observation with moderate-or-higher Atlas tier, current-query marker/anti-marker support, an independent internal-anchor/spatial channel and no OOD/ontology conflict returns broad-only with `fine_anchor_eligible=false`; lower confidence remains QC. Defined broad labels are never overwritten by mapping alone.
10. Reopen a complete cluster or frozen cohort once when a calibrated alternative differs in a material fraction, a material ontology conflict remains after crosswalk inspection, or a coherent group is out-of-distribution. Use query full-feature positive/anti evidence, pseudobulk, stability, sample consistency, spatial morphology and technical alternatives to retain, supersede, downgrade or mark `Unknown candidate`. Close every trigger with `validate_global_atlas_concordance.py`; do not chase the reference iteratively.
11. Validate evidence content and exact membership closure. Prelabel freezes, cohort outcomes, direct returns, Atlas concordance/reviews and per-label supports must pass their schemas and hashes; empty/status-only evidence fails.
12. Materialize one final annotation: moderate-or-higher broad labels, optional high-confidence fine labels, or explicit retained QC/technical/unknown state.
13. Run `validate_lineage_signal_coverage.py`; pass direct-workflow, state, taxonomy, completion and main-Agent biological-quality gates; generate the lightweight confirmation HTML including label-independent all-cell canonical-marker spatial panels; wait for explicit user confirmation; then produce final DEG, figures and release HTML.

A cohort is a frozen computational query boundary, not a cell type. Direct cross-lineage return preserves source cohort/subcluster/membership/evidence and does not create an intermediate cohort or automatically trigger another target-lineage reclustering.

## Context, inputs and framework

Read `references/context-and-biology.md`. Require and validate a context JSON containing species, tissue, stage/condition, platform, observation unit, biological questions and priority lineages. Run `scripts/resolve_workflow_profile.py`; profiles supply candidate boundaries, marker/anti-marker gates, references and technical constraints, never a label map.

Start every project with:

1. `scripts/discover_inputs.py INPUT_ROOT --out PROJECT_ROOT/input_discovery`
2. `scripts/inspect_r_object.R` for R objects plus `scripts/check_runtime.py` and `scripts/check_r_runtime.R`
3. `scripts/init_annotation_project.py --sample SAMPLE --input-root INPUT_ROOT --project-root PROJECT_ROOT --modality spatial|single-cell`
4. `scripts/register_input_snapshot.py` for every active object and reused clustering table
5. `scripts/audit_feature_scope.R`; if Seurat `Spatial@data` is absent or equals counts, build a separate validation-only LogNormalize object with `scripts/prepare_seurat_full_feature_validation.R`
6. freeze exact analysis membership and SHA256 in `provenance/analysis_scope_policy.json`

When a readable full-feature Seurat RDS exists or R is requested, read `references/r-first-workflow.md` and use Seurat/SCT graphs as the computational backbone. BANKSY is an optional spatial clustering adapter/corroborating view, not the controller. Reuse computation only after object/cell hashes, feature scope, parameters and outputs pass; never reuse historical labels.

When a Seurat RDS was converted from StereoPy `cellbin_PPed`, read `references/seurat-cellbin-preprocessing.md`. Treat it as a raw-count carrier, not an SCT object; never reuse imported StereoPy PCA/UMAP as the R graph.

## Clustering and evidence

Read `references/clustering-selection.md` and `references/evidence-routing.md`. Quantitative ranking only shortlists candidates. Select the complete-grid candidate with the best integrated evidence for the current biological question: preserve stable lineages/true subtypes first, avoid state/technical fragmentation second, and use lower complexity only as a tie-breaker. Generate every cluster/cohort DEG, UMAP, spatial overview and highlight; final positive/negative evidence uses the full-feature object.

Build anchors from coherent multi-gene backbone/support groups plus anti-programs, not marker-hit counts. Use `profile_anchor_programs.R`, `build_program_anchor_membership.R` and `summarize_cluster_programs.R` where applicable. Missing genes in an HVG adapter are unassayed, not negative.

Use a shallow tree. A graph subcluster is not a subtype. ECM, contractility, stress, cell cycle, hypoxia, low RNA, ambient signal, anatomical adjacency or marker intensity remain tags unless an independent functional/lineage program passes the fine-label gate. Zero fine labels is valid.

## Sheep-ovary profile

Sheep/Ovis/ovine/羊 plus ovary/ovarian/卵巢 activates the sheep evidence profile. Read:

- `references/profiles/sheep_ovary_standard_workflow.md`
- `references/profiles/sheep_ovary_candidate_lineage_catalog.json`
- `references/profiles/sheep_ovary_literature_2025_2026.md`
- `references/profiles/sheep_ovary_rfirst_case_reference.md` as a sanitized strategy regression trace only
- `references/profiles/sheep_ovary_oocyte_rfirst_route.md` before any Oocyte decision

The same-batch R-first preset freezes technical preprocessing, not biological outcomes. Use `--strategy-preset sheep_ovary_same_batch_rfirst`, resolve `config/active_strategy_preset.json`, and independently select whole-tissue/cohort resolutions and labels.

For a StereoPy cellbin RDS, activate the same-batch contract only after the workflow profile, explicit preset, object/layer/marker/coordinate audits and conversion provenance are hash-verified. Path names and feature counts are nonbinding hints. Then require the declared `Spatial -> SCT` contract and new PCA/neighbour graph. The formal whole-tissue and cohort resolution grid is `0.1,0.2,0.3,0.4,0.6`. Do not run `0.01`, `0.02` or `0.05`; if 0.1 produces implausible inflation, repair the graph. Final QC is not reclustered.

Run `init_open_world_lineage_audit.py` and `validate_open_world_lineage_audit.py`. Review all candidate boundaries even when absent. The catalog is non-exhaustive: coherent unexplained programs require added review. `Stromal/mesenchymal` is an honest parent; standalone Theca, Smooth muscle, Pericyte/mural, Mesenchymal progenitor-like, luteal, lymphatic, Neural/Schwann/glia or other context-dependent classes require their own complete programs and morphology. Do not use `Theca/follicular wall` or `Stromal/perivascular` as catch-all release labels.

Oocyte uses one contamination-safe targeted cohort containing the complete multi-module recall set. Strict seeds/spatial foci support identity but do not define cohort membership or final census. Require coherent non-ZP identity plus maternal/ooplasm enrichment, somatic anti-program clearance and compatible object morphology. Cortical/section-edge location is never negative evidence and location alone is never positive evidence. Return pregranulosa/granulosa or stromal clusters directly to those lineages with zona/adjacency tags. Report cellbin count separately from putative object count.

Without a usable matched count-level reference, GSE233801 is the primary public adult-sheep somatic Atlas for the terminal all-cell broad concordance pass. It may directly rescue only the frozen residual QC membership; for already-defined cells it is a diagnostic challenger. It cannot rescue Oocyte, Theca or Epithelial/mesothelial automatically. A matched dotplot is marker evidence, not cell-level transfer.

## Atlas and reference calibration

Read `references/matched-single-cell-reference.md`. A count-level matched, stage-compatible reference may be the preferred external channel only after all broad/targeted cohorts and terminal QC membership are frozen. Map the analysis set once for broad concordance, but permit direct writeback only for frozen QC. Current-query full-feature marker/anti-marker and morphology remain authoritative.

Calibration uses disjoint query-like held-out current-query anchors. External reference self-splitting is diagnostic only. Default target-precision tiers are moderate-or-higher 0.90 and high 0.95; these are calibration targets, not universal per-cell score cutoffs. Report mutually exclusive `high`, `moderate_only`, `low_reject` and require `moderate_or_higher_n = high_n + moderate_only_n`.

Use `adjudicate_multichannel_broad_rescue.py` only to decide writeback for frozen terminal QC. Require marker/anti-marker support plus at least one internal-anchor or observed-density spatial channel. Every accepted return participates in final broad DEG/dotplots but cannot seed fine discovery. Precompute the Atlas feature transform, low-dimensional representation and approximate-neighbor index once and validate `atlas_index_manifest.schema.json` with `validate_atlas_index_manifest.py`; never form a full query-by-reference distance matrix or repeat joint Atlas integration per sample. Treat low margin, mixed neighbors and coherent OOD programs as unknown/review signals, not nearest-class assignments.

## State, autonomy and multiple samples

Read `references/autonomous-operation.md` for end-to-end work. Run `autopilot_status.py PROJECT_ROOT` at startup and after every writeback. Run `plan_next_iteration.py` after every biological decision; it consumes structured gap codes rather than English error strings. In Agent mode, `CONTINUE`, `ITERATION_REQUIRED` and `EXPECTED_GATE_BLOCKED` are expected business states, while schema/parse/corruption exceptions are execution failures. Use strict nonzero gap exits only in CI/completion. Validate failed jobs under new run IDs and record every incident; do not overwrite failure evidence.

New projects use `recluster_cohort_registry.tsv`, `direct_return_registry.tsv`, `route_attempt_registry.tsv`, `lineage_signal_boundary_registry.tsv` and `lineage_signal_registry.tsv`. The signal registries keep weak lineage evidence alive across whole-tissue, broad and targeted boundaries until explicitly supported or refuted. Retired registry formats remain readable only through migration tools; never create them in a new project.

Use `scripts/migrate_project_v1_7_to_v1_8.py` to opt an existing v1.7 project into these fail-closed gates. Migration creates empty registries and intentionally blocks completion until the historical boundaries are backfilled and audited; it never fabricates evidence.

For multi-sample work, read `references/multi-sample-agent-orchestration.md`. The main conversation Agent owns progress, user decisions, cross-sample audit and final quality approval. Assign exactly one complete workflow child Agent per sample; parallel workers must meet identical gates and must not be reduced to cluster renaming or audit-only tasks.

## Completion and report

Before confirmation:

1. run `build_final_annotation.py`;
2. run `validate_prelabel_broad_evidence.py`, `validate_lineage_signal_coverage.py`, `validate_global_atlas_concordance.py`, `audit_release_taxonomy.py`, `validate_direct_lineage_workflow.py`, `validate_annotation_support_registry.py` and `audit_annotation_membership_partition.py`;
3. run `validate_state.py`, `plan_next_iteration.py` and `check_completion_gate.py`;
4. populate `annotation_support_registry.tsv` and build frozen lightweight broad spatial/canonical-marker assets plus all-cell canonical-marker spatial panels (no label filtering; fixed point size);
5. request and record one main-Agent biological quality approval only after all annotation phases finish;
6. build the lightweight confirmation HTML and request explicit user confirmation.

Do not spend compute on final assets before confirmation. After confirmation, follow `references/report-contract.md` and produce:

- broad DEG, UMAP/spatial maps and per-node highlights using every final broad member, including direct and Atlas returns;
- high-confidence fine DEG/maps only when real fine labels exist;
- separate canonical and current-data-specific broad tree dotplots, plus fine tree dotplots when applicable, in PNG/PDF/source TSV;
- dotplot absolute values plus display size/color normalized within gene;
- marker spatial maps grouped by supported cell type;
- Chinese/requested-language interactive HTML with annotated overview, expandable tree, node jumps, evidence, workflow timeline and raw state;
- cell/cluster/cohort/direct-return/route/run/incident registries, manifests, session information and checksums.

Run `build_release_manifest.py` and `audit_release.py`. Completion requires `autopilot_status.py` = `COMPLETE` and release audit PASS.

## Execution and testing

Local execution is for discovery and small audits. Submit heavy data operations to compute nodes using `references/job-orchestration.md`. Match requested CPUs to real parallelism; parallelize independent resolutions/jobs rather than reserving many cores for single-thread Leiden or unwrapped `FindAllMarkers`.

Before distribution, run repository validation, Python/R syntax checks, tests, install verification and a fresh-agent forward test. The fresh Agent receives raw inputs and the Skill, not the intended clustering, labels or historical answer.
