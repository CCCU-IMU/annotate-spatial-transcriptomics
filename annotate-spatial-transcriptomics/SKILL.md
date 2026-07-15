---
name: annotate-spatial-transcriptomics
description: Independently annotate spatial transcriptomics or single-cell RNA-seq with adaptive broad clustering, broad-class and targeted reclustering cohorts, direct cross-lineage returns, calibrated residual-QC Atlas/RCTD assistance, immutable state, shallow high-confidence subtypes and audited interactive reports. Seurat/R-first is preferred for full-feature RDS inputs; AnnData/H5AD, SingleCellExperiment, BANKSY, Scanpy/Leiden and external cluster tables are supported. Use for end-to-end or multi-sample annotation, evidence review, recovery of unresolved observations and release reporting. Never copy parameters or labels from example projects.
---

# Annotate Spatial Transcriptomics

Use one evidence-first state machine for every new project. A first clustering proposes biological questions; it is never the final annotation. Persistent biological review pools are not part of the active architecture. Legacy pool registries are migration-only provenance.

Read `references/direct-lineage-controller.md`, `references/iterative-controller.md`, `references/quality-standard.md` and `references/state-schema.md` before analysis.

## Default annotation architecture

Execute this order without skipping phases:

1. Discover, inspect and hash inputs; freeze `full_object`, `analysis_set` and `excluded_initial_qc` memberships.
2. Select a whole-tissue broad clustering adaptively from current DEG, marker/anti-marker, stability, UMAP and spatial evidence.
3. Review every cluster with an open-world lineage catalog. Assign supported moderate-or-higher initial broad labels directly; send low-information, featureless or irreducibly mixed observations to `qc_holdout`. Do not create biological pools.
4. Freeze one immutable `broad_class_recluster` cohort for every supported initial broad class. Recluster and adjudicate every subcluster. A genuinely underpowered class may remain broad-only after a recorded `not_applicable_reviewed` skip.
5. Write each subcluster directly as a high-confidence shallow fine label, its supported parent broad class, a supported cross-lineage broad/fine return, one decision-relevant targeted cohort, or QC/technical retention.
6. Create `targeted_recluster` only for a local interpretable mixture, contamination boundary or context-gated identity. It answers one question and cannot become a long-lived catch-all.
7. Use RCTD/reference assistance only when applicable. Extreme plus independent marker/anti-marker, resolution and spatial evidence may support fine; high supports broad-only; medium/low enters QC.
8. After every broad/targeted cohort is terminal, freeze the complete residual `qc_holdout` membership once. Do not recluster it. Apply calibrated Atlas/internal-anchor/marker/spatial consensus only to that exact membership. Moderate-or-higher returns broad-only with `fine_anchor_eligible=false`; lower confidence remains QC reject.
9. Materialize one final annotation: moderate-or-higher broad labels, optional high-confidence fine labels, or explicit retained QC/technical state.
10. Pass direct-workflow, state, taxonomy, completion and main-Agent biological-quality gates; generate the lightweight confirmation HTML; wait for explicit user confirmation; then produce final DEG, figures and release HTML.

A cohort is a frozen computational query boundary, not a cell type. Direct cross-lineage return preserves source cohort/subcluster/membership/evidence and does not create an intermediate target pool or automatically trigger another target-lineage reclustering.

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

Read `references/clustering-selection.md` and `references/evidence-routing.md`. Quantitative ranking only shortlists candidates. Freeze the lowest-complexity resolution preserving defensible broad lineages. Generate every cluster/cohort DEG, UMAP, spatial overview and highlight; final positive/negative evidence uses the full-feature object.

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

For same-batch StereoPy cellbin RDS, require the declared `Spatial -> SCT` contract and new PCA/neighbour graph. The formal whole-tissue and cohort resolution grid is `0.1,0.2,0.3,0.4,0.6`. Do not run `0.01`, `0.02` or `0.05`; if 0.1 produces implausible inflation, repair the graph. Final QC is not reclustered.

Run `init_open_world_lineage_audit.py` and `validate_open_world_lineage_audit.py`. Review all candidate boundaries even when absent. The catalog is non-exhaustive: coherent unexplained programs require added review. `Stromal/mesenchymal` is an honest parent; standalone Theca, Smooth muscle, Pericyte/mural, Mesenchymal progenitor-like, luteal, lymphatic, Neural/Schwann/glia or other context-dependent classes require their own complete programs and morphology. Do not use `Theca/follicular wall` or `Stromal/perivascular` as catch-all release labels.

Oocyte uses one contamination-safe targeted cohort containing the complete multi-module recall set. Strict seeds/spatial foci support identity but do not define cohort membership or final census. Require coherent non-ZP identity plus maternal/ooplasm enrichment, somatic anti-program clearance and compatible object morphology. Cortical/section-edge location is never negative evidence and location alone is never positive evidence. Return pregranulosa/granulosa or stromal clusters directly to those lineages with zona/adjacency tags. Report cellbin count separately from putative object count.

Without a usable matched count-level reference, GSE233801 is the primary public adult-sheep somatic Atlas only for the terminal residual QC membership. It does not classify the full object and cannot rescue Oocyte, Theca or Epithelial/mesothelial automatically. A matched dotplot is marker evidence, not cell-level transfer.

## Atlas and reference calibration

Read `references/matched-single-cell-reference.md`. A count-level matched, stage-compatible reference may be the preferred external channel only in terminal residual QC. Current-query full-feature marker/anti-marker and morphology remain authoritative.

Calibration uses disjoint query-like held-out current-query anchors. External reference self-splitting is diagnostic only. Default target-precision tiers are moderate-or-higher 0.90 and high 0.95; these are calibration targets, not universal per-cell score cutoffs. Report mutually exclusive `high`, `moderate-only`, `low-reject` and require `moderate_or_higher_n = high_n + moderate_only_n`.

Use `adjudicate_multichannel_broad_rescue.py` only on the frozen terminal QC membership. Require marker/anti-marker support plus at least one internal-anchor or observed-density spatial channel. Every accepted return participates in final broad DEG/dotplots but cannot seed fine discovery.

## State, autonomy and multiple samples

Read `references/autonomous-operation.md` for end-to-end work. Run `autopilot_status.py PROJECT_ROOT` at startup and after every writeback. Run `plan_next_iteration.py` after every biological decision. Validate failed jobs under new run IDs and record every incident; do not overwrite failure evidence.

New projects use `recluster_cohort_registry.tsv`, `direct_return_registry.tsv` and `route_attempt_registry.tsv`. Legacy `pool_registry.tsv`, `pool_snapshot_registry.tsv` and `branch_control_board.tsv` remain readable only for migration. Do not start new biological pools.

For multi-sample work, read `references/multi-sample-agent-orchestration.md`. The main conversation Agent owns progress, user decisions, cross-sample audit and final quality approval. Assign exactly one complete workflow child Agent per sample; parallel workers must meet identical gates and must not be reduced to cluster renaming or audit-only tasks.

## Completion and report

Before confirmation:

1. run `build_final_annotation.py`;
2. run `audit_release_taxonomy.py` and `validate_direct_lineage_workflow.py`;
3. run `validate_state.py`, `plan_next_iteration.py` and `check_completion_gate.py`;
4. populate `annotation_support_registry.tsv` and build frozen lightweight broad spatial/canonical-marker assets;
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
