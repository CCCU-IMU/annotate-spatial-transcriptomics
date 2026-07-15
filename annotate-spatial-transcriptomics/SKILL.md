---
name: annotate-spatial-transcriptomics
description: Independently annotate spatial transcriptomics or single-cell RNA-seq projects with a Seurat/R-first default when a full-feature RDS is available, while supporting AnnData/H5AD, SingleCellExperiment, BANKSY, Scanpy/Leiden and external cluster tables. Use for input discovery, adaptive clustering, broad and optional subtype annotation, full-gene marker/spatial review, iterative unresolved-pool routing, calibrated atlas/RCTD assistance, immutable state tracking, baseline quality comparison and audited HTML reports. Do not copy parameters or labels from example projects.
---

# Annotate Spatial Transcriptomics

Apply an evidence-first, stateful annotation workflow. Treat example projects as strategy demonstrations only. Derive parameters, clusters, thresholds and labels from the current data.

When a readable full-feature Seurat RDS exists, or the user requests an R-centered analysis, read `references/r-first-workflow.md` and use it as the primary execution path. Use SCTransform/Seurat graphs and R pool reclustering as the default computational backbone; borrow the completed Scanpy workflow's iterative decision strategy, not its parameters or labels. Treat BANKSY as an optional spatial clustering adapter or corroborating view, not the default controller. Reuse existing Seurat clustering or pool-reclustering artifacts only after validating object/cell hashes, feature scope, parameters and output completeness; reusing computation never authorizes reusing biological labels.

Read `references/context-and-biology.md` first. Require a biological-context JSON containing species, tissue, stage/condition, platform, observation unit, primary questions and priority lineages. Validate it with `scripts/validate_biological_context.py`. Load a matching tissue profile when available, but never treat it as a label map.

After context and input inspection, run `scripts/resolve_workflow_profile.py`. Sheep/Ovis/ovine/羊 plus ovary/ovarian/卵巢 automatically activates the sheep-ovary evidence profile. If a readable full-feature Seurat RDS exists, use the R-first backbone. If inspection also confirms a `Spatial` raw-count assay and StereoPy `cellbin_PPed` conversion provenance/path, the frozen whole-tissue preprocessing contract is mandatory; numerical deviations require an explicit recorded batch exception. Read `references/profiles/sheep_ovary_literature_2025_2026.md` for the 2025 Science, 2026 Advanced Science and 2025 AJOG boundary evidence. These papers define candidate boundaries and anti-overannotation rules, not labels or subtype quotas.

When the user explicitly asks for the same-batch sheep-ovary R-first standard preset, initialize with `--strategy-preset sheep_ovary_same_batch_rfirst` and resolve with the same option to `config/active_strategy_preset.json`. Read `references/profiles/sheep_ovary_same_batch_rfirst_preset.json` and `sheep_ovary_standard_workflow.md`. Reuse the validated phase order, pool axes, evidence gates, rescue order, Oocyte safeguards, state/report contract and post-completion master quality review. Never reuse the reference case's selected resolution, cluster-to-label mapping, membership, label counts or subtype inventory. The completion gate rejects a requested preset whose bindings or independently verified same-batch preprocessing provenance are missing.

For a sheep-ovary execution, read `references/profiles/sheep_ovary_standard_workflow.md` first and follow its low-freedom phase order. Also read `references/profiles/sheep_ovary_rfirst_case_reference.md` as a sanitized strategy regression trace, never as a label map. Before any Oocyte route, read `references/profiles/sheep_ovary_oocyte_rfirst_route.md`; its two-tier candidate/seed distinction is mandatory.

During sheep-ovary broad annotation, load `references/profiles/sheep_ovary_candidate_lineage_catalog.json` and initialize its ledger-bound scaffold with `scripts/init_open_world_lineage_audit.py`. Review every catalog boundary even when the outcome is negative, then record any coherent unexplained multi-gene program as an additional candidate and run `scripts/validate_open_world_lineage_audit.py`. The catalog is non-exhaustive and never requires a class to be present. Do not collapse this audit to a few example lineages, and do not promote state patterns such as ECM-rich, contractile-only, stress, hypoxia, low-RNA or anatomical adjacency into cell types.

Keep the two sheep profiles distinct. `sheep_ovary_rfirst_profile.json` has role `workflow_preprocessing`; it controls preprocessing, normal resolution grids and reference policy. `sheep_ovary.json` has role `biological_evidence`; it controls taxonomy, marker/anti-marker and multi-route completion. Run `scripts/validate_profile_role.py` and never pass the workflow profile to a biological completion/audit argument.

When a biologically matched single-cell reference is supplied, read `references/matched-single-cell-reference.md`. Harmonize source labels to a shallow candidate broad vocabulary with an explicit crosswalk and transfer ceiling; do not force the spatial query to reproduce the reference taxonomy. A reference dotplot can refine marker hypotheses but cannot perform cell-level transfer without a count-level object. Current-query full-feature marker/anti-marker evidence and spatial morphology remain authoritative when they contradict the reference.

Read `references/iterative-controller.md` before analysis. The workflow is a fail-closed multi-round state machine, not cluster-name assignment. A first pass, a high annotation rate or a rendered report is never completion.

Read `references/quality-standard.md` before freezing labels. For spatial cellbins/spots, accurate broad classes are the primary endpoint; fine subtypes are optional and must never be forced merely to make the tree deeper. Comparable quality across clustering methods means independently supported biological decisions and complete route exhaustion, not matching a previous pipeline's labels.

Read `references/taxonomy-and-pool-design.md` before creating biological pools or release labels. Keep the published candidate-lineage checklist, overinclusive analysis pools and final biological broad classes as three separate layers. A literature class is not a required output, an analysis-pool name is never a release label, and interface/QC states are reported outside the biological broad-class census.

Default to a shallow annotation tree. A graph subcluster is not a biological subtype, and several stable subclusters may merge into one broad or literature-supported functional label. Treat ECM-rich, contractile, stress, cell-cycle, low-RNA, ambient, cortical and follicle-adjacent patterns as state/spatial tags unless an independent lineage or functional program passes the fine-label gate. For sheep ovary, use the limited subtype vocabulary supported by recent sheep studies and the current section; never create one named subtype per Seurat cluster.

Subtype restraint is a release criterion, not a cosmetic preference. The default maximum depth is broad class plus one functional subtype level, and zero subtypes is valid. Before keeping any fine label, explicitly state what biological conclusion would be lost by merging it into the parent; if the answer is only cluster identity, marker intensity, ECM amount, state or location, merge it and retain those properties as tags.

When the user requests end-to-end or minimally supervised execution, read `references/autonomous-operation.md` and keep working through its continuous control loop. Run `scripts/autopilot_status.py PROJECT_ROOT` at startup and after every writeback or release step. Do not pause for routine resolution, pool, label or job-repair approvals. After **all** pool reclustering, open-world/context-specific work, rescue writeback, final-annotation materialization and the completion gate have passed, the main conversation Agent must perform one concise biological quality approval against the validated sheep-ovary R-first reference. This is not a second mechanical completion audit and may pass with documented concerns; it asks whether broad labels, marker/anti-marker and spatial evidence, and context-specific/confounded calls reach comparable quality. It must never run after broad annotation alone. Only after this approval, generate the lightweight confirmation HTML and present the frozen review/census for explicit user confirmation before spending compute on final DEG, full tree dotplots, per-node/per-gene spatial assets and the release HTML.

Read `references/efficient-operation.md` for long-running or multi-sample projects. Its cache, handoff, real-worker and single-report rules reduce token and compute cost without weakening evidence gates.

When the user requests multiple samples or parallel annotation, read `references/multi-sample-agent-orchestration.md`. The main conversation Agent is the only user-facing master and owns progress, biological decision gates, cross-sample audit and release. Assign exactly one full-workflow child Agent per sample, with an isolated project root. Parallel workers must satisfy the same evidence, state and report gates as the master; never create multiple reduced audit Agents for one sample or combine samples because Agent slots are limited. Initialize and audit the cohort control board with `scripts/init_annotation_cohort.py` and `scripts/validate_cohort_state.py`.

Read `references/multi-route-controller.md` before creating unresolved pools. Its ordered controller is the core of this Skill: large usable populations undergo real balanced-anchor reclustering; local interfaces undergo targeted anchors then calibrated deconvolution/reference review when eligible; post-clustering low-information populations undergo a full QC-pool anchor-recluster before atlas rescue. Do not create a generic rare-cell route. First run open-world, literature-informed lineage discovery; Neural/Schwann, luteal, Smooth muscle, Oocyte and any additional coherent program may create sample-specific ordinary biological pools. Oocyte alone retains a dedicated contamination-safe route when a candidate exists because adjacent zona signal is a known confounder. A route name without the required query/anchor boundary and validation artifacts does not count as an attempt.

## Start every project

1. Run `scripts/discover_inputs.py INPUT_ROOT --out PROJECT_ROOT/input_discovery`.
2. Read `references/input-adapters.md` and identify the expression object, every clustering candidate, coordinates, reductions and existing progress/state files.
3. Inspect candidate R objects with `scripts/inspect_r_object.R`, then run `scripts/resolve_workflow_profile.py`; this workflow decision precedes clustering and is not a biological label assignment.
4. Run `scripts/check_runtime.py` and `scripts/check_r_runtime.R` in candidate environments before heavy work; select an existing compatible environment and record it rather than silently installing packages.
5. Run `scripts/init_annotation_project.py --sample SAMPLE --input-root INPUT_ROOT --project-root PROJECT_ROOT --modality spatial|single-cell`.
6. Freeze every active expression object and reused clustering table with `scripts/register_input_snapshot.py` before analysis. Never silently replace an object or clustering table.
7. If existing state exists, validate and resume it. Do not restart completed pools.
8. Distinguish clustering features from validation features. Run `scripts/audit_feature_scope.R` on the expression object used for final marker evidence. HVG/BANKSY features may drive clustering, but final positive/negative marker, open-world discovery and context-specific validation must use a full-feature object when the profile requires it. For Seurat `Spatial`, never assume `@data` is normalized: if `data` is absent or exactly equals `counts`, keep the SCT clustering object immutable and build a separate full-feature validation-only object with `scripts/prepare_seurat_full_feature_validation.R`. Any Wilcoxon DEG/marker job on `Spatial` must pass the resulting manifest to the evidence runner; fail closed otherwise.
9. Separate `full_object`, `analysis_set`, `excluded_initial_qc` and `postcluster_holdout`. Never mix initial QC exclusion with an unresolved biological pool.
10. Freeze the analysis-set membership and SHA256 in `provenance/analysis_scope_policy.json`. Before release, require the cell ledger's scope assignments and the single final annotation to match that exact membership.
11. If a matched single-cell reference exists, register its object/annotation/marker artifacts, copy `assets/matched_reference_crosswalk_template.tsv` to `config/matched_reference_crosswalk.tsv`, preserve every source label verbatim and run `scripts/validate_matched_reference_crosswalk.py`. Use a packaged tissue alias table only as a starting hypothesis.

When a Seurat RDS was batch-converted from StereoPy `cellbin_PPed`, read `references/seurat-cellbin-preprocessing.md` before clustering. Treat the converted RDS as a raw-count input container, not an SCT object. For samples from the same production batch, use `scripts/run_seurat_sct_preprocess.R` and its frozen whole-tissue SCT/PCA/neighbour profile so preprocessing is comparable; keep resolution selection and biological annotation adaptive. Never reuse imported StereoPy PCA/UMAP as the R graph, and never infer a valid SCT run without a matching preprocessing manifest.

## Select clustering adaptively

Read `references/clustering-selection.md`. For an R-first project, inspect or generate a Seurat resolution grid, rank it with the same marker/stability/spatial criteria, and freeze the lowest-complexity resolution that preserves defensible broad lineages. For a BANKSY input, run `scripts/rank_banksy_grid.py GRID_SUMMARY --out PROJECT_ROOT/selection` only as the BANKSY-specific adapter.

Use quantitative ranking only to shortlist candidates. Inspect marker interpretability, cluster-size distribution, adjacent-parameter stability, UMAP structure and spatial morphology before freezing the selected run. Never select by filename, a fixed resolution, or an example's choice.

For sheep ovary, run the complete formal grid `0.1,0.2,0.3,0.4,0.6` for whole tissue and every biological/QC pool. Validate it with `scripts/validate_resolution_grid.py` before submission. Do not run or use `0.01`, `0.02` or `0.05`; if `0.1` produces implausible cluster inflation, repair the graph rather than lowering resolution. For a shortlisted grid, run `scripts/compare_clusterings.py` to quantify ARI/AMI. Retain the all-observation ARI/AMI, including every microcluster. A second `n>=100` macro-restricted score may help rank broad resolutions, but it never removes observations or exempts small clusters from DEG, spatial maps and rare-lineage/technical review.

Record every candidate reviewed and the final rationale in `state/clustering_decision_ledger.tsv`.

## Annotate broad classes first

Read `references/evidence-routing.md`.

1. Generate per-cluster one-vs-rest DEG, UMAP and spatial maps.
2. Examine canonical marker programs, data-specific DEG, absolute detection, spatial location, QC complexity and contradictory programs.
3. Assign only defensible broad classes. Keep uncertainty explicit.
4. Create immutable memberships for defined anchors and unresolved pools.
5. Do not infer a cell identity from one marker, one reference label or spatial proximity alone.

Build anchors from coherent marker groups, not only total marker hits. Use `scripts/profile_anchor_programs.R` to inspect feature coverage and threshold sensitivity, then `scripts/build_program_anchor_membership.R` with required lineage-backbone/support groups, anti-programs, depth limits and minimum anchor counts. After reclustering, use `scripts/summarize_cluster_programs.R` to compare cell-level program and anti-program coherence across every candidate resolution. A missing or tiny anchor class blocks or narrows the next route; it is not permission to lower specificity blindly.

For ovary or another mesenchymal-rich tissue, run a mandatory broad-lineage decomposition audit before closing the stromal background. Test generic stromal/fibroblast, mesenchymal-progenitor-like, smooth-muscle, pericyte/mural, endothelial and steroidogenic-theca programs separately. `ACTA2`/`TAGLN` alone cannot define Smooth muscle; `DCN`/`PDGFRA` alone cannot prove a distinct Mesenchymal class. Use mature contractile backbones, vascular adjacency, steroidogenic anti-programs, stable R pool reclustering and morphology. If an independent class is unsupported, record the negative audit and retain it under the biologically honest parent rather than inventing a missing literature category.

Use axis-covering biological review pools rather than graph-cluster names or a one-pool-per-literature-label design. Preserve `source_key`, original parent/cluster, generation, competing candidate lineages, `state_tags`, `spatial_tags` and `qc_tags`. Pool names must signal routing (`*_review`, `*_candidate`, `*_unresolved`, `*_holdout`) and cannot be copied directly into final labels. ECM-rich, contractile, cortical, ambient or low-RNA are state tags, not pool identities.

For sheep ovary, use `Stromal/mesenchymal` as the honest generic stromal parent. Release standalone `Theca`, `Smooth muscle`, `Pericyte/mural` or `Mesenchymal progenitor-like` only when their independent gates pass. Never publish `Theca/follicular wall` as a broad label: steroidogenic theca must be separated from structural stroma, mature smooth muscle, mural cells and follicular interfaces. A documented negative audit for a literature category is a valid outcome.

For sheep immune review, account for incomplete immunoglobulin gene naming. `LOC101108817`, `LOC101108781`, `LOC121817142` and `LOC114108841` are sheep immunoglobulin constant-region loci. A stable cluster with a coherent multi-locus immunoglobulin program plus an independent B/plasma regulator such as `JCHAIN`, `POU2AF1`, `TENT5C` or `MZB1` may support broad `Immune` even when `PTPRC` is sparse. A single immunoglobulin locus, `JCHAIN`, `CD74` or MHC signal remains insufficient. Keep antibody-secreting/plasma-like identity as a state tag or optional fine label unless its own fine-label gate passes.

Missing genes in an HVG/reduced adapter are unassayed, not negative. Never downgrade a lineage for an absent backbone until the full-feature audit passes.

## Route unresolved populations

Choose the route by failure mode:

- Large continuous populations with usable signal: return to a canonical broad review pool and run balanced reference-only-anchor reclustering with query-only graph/UMAP/DEG.
- Local mixed/interface populations: use targeted anchor reclustering, then run a machine-readable RCTD/reference applicability audit. If eligible, calibrated RCTD assistance is used at low priority; medium/low and other unresolved interface observations enter the frozen post-clustering QC holdout rather than calling Atlas directly. If RCTD is ineligible, record the failing criteria and route unresolved observations to the same appropriate biological/QC successor.
- Low-complexity populations: retain as post-clustering QC holdout, combine all eligible residuals from earlier routes, and run one complete QC-pool anchor reclustering first. Only cells still retained in QC after that reclustering may combine Atlas, internal-anchor, marker and observed-density spatial evidence. Rescue calibrated **moderate-or-high confidence** Atlas labels to broad classes; observations below the calibrated moderate tier remain rejects/review. Never run Atlas routinely over already defined broad classes, high-confidence fine labels or ordinary biological pools.
- Technical or irreducible populations: retain as technical/QC/review; do not force annotation.

Broad-only rescues must have `fine_anchor_eligible=false`. Never use atlas top label or RCTD output as an unvalidated forced classifier.

For sheep ovary without a usable matched count-level reference, GSE233801 is the primary public adult-sheep somatic Atlas **only for the residual post-anchor QC holdout**. Its permitted broad rescue labels are Granulosa, Stromal/mesenchymal, Vascular/endothelial, Pericyte/mural candidates and Immune; it is not an automatic classifier for the full object or biological pools, and it cannot rescue Oocyte, Theca or Epithelial/mesothelial. A matched dotplot is candidate marker evidence only and does not enable transfer. A matched, stage-compatible count object may refine the residual-QC external channel after current-query anchors, while GSE233801 remains an independent sheep reference.

For atlas mapping, distinguish calibration performance from an individual prediction score. By default, calibrate class/route-specific held-out target precision tiers at `moderate-or-higher=0.90` and `high=0.95` (or predeclare justified project-specific alternatives). Use nested cumulative thresholds: every high observation also meets the moderate-or-higher gate. Assign mutually exclusive output labels as `high`, `moderate-only` and `low-reject`, and always report `moderate_or_higher_n = high_n + moderate_only_n`. Both accepted tiers are eligible for broad-only return after independent marker/anti-marker and route/spatial or internal-anchor checks. `0.90/0.95` are calibration targets, not universal per-observation raw-score cutoffs. Final calibration must use disjoint query-like held-out current-query anchors and a machine-readable origin manifest. Splitting one external atlas into training and held-out cells is reference self-classification and diagnostic only; it cannot calibrate query rescue. External Atlas, internal-anchor, marker and observed-density spatial evidence should be combined as a calibrated consensus route when query-like held-out anchors are available; Atlas self-consistency alone cannot define final confidence. The legacy combined `medium_high` calibrator is diagnostic-only and never writeback eligible.

For the frozen residual QC-holdout membership only, materialize the independent channels in one observation-level table and run `scripts/adjudicate_multichannel_broad_rescue.py` with a project-specific copy of `assets/multichannel_broad_rescue_config_template.json`. Bind the input to the validated QC-anchor-recluster outcome and its hash. Calibrate on query-depth-matched held-out anchors, not the query itself. Require current-query marker/anti-marker support and at least one internal-anchor or observed-density spatial channel. Treat its output as a broad-only proposal: inspect the route/space summary, then use an atomic ledger commit. Never apply this step to all broad/fine cells and never copy the template's channel minima blindly.

Treat RCTD as lower-priority assistance. Only an extreme-confidence RCTD call with independent marker, anti-marker, resolution and spatial support may contribute to a fine return; a high-confidence call may return broad-only. Medium/low-confidence observations first enter the post-clustering QC holdout, join its complete QC-anchor reclustering, and only residual rejects may continue to calibrated Atlas/internal-anchor review. A large or spatially diffuse RCTD reject pool is not a local interface and cannot close merely because deconvolution was attempted.

Re-audit large first-pass direct definitions. A graph cluster with a convincing pseudobulk marker can still contain a small true lineage embedded in a large resident population. Require cell-level program purity, anti-programs, spatial continuity and, for large/heterogeneous clusters, anchor-assisted query-only reclustering before freezing it as a lineage anchor.

Every route attempt is recorded in `state/route_attempt_registry.tsv`; every branch and no-repeat status is recorded in `state/branch_control_board.tsv`. `explicitly_retained_closed` is not sufficient evidence for a large or priority pool.

An original graph cluster may be split at observation level. Context-gated identities such as Oocyte must pass the profile's positive, negative, spatial-object and reclustering gates; otherwise roll them back into biologically appropriate somatic/interface pools. For ovarian Oocyte review, cortical, subcortical, peripheral or section-edge location is never negative evidence; small/primordial oocytes may be cortical, while location alone is never sufficient positive evidence.

For sheep-ovary Oocyte review, never confuse recall membership with identity evidence. Every observation passing the predeclared multi-module starting gate enters the full query-only candidate recluster pool. High-specificity marker/anti-program observations and compact spatial foci are **strict seeds/support**, not the candidate census and not the only permitted recluster membership. Retain isolated high-evidence candidates in the full pool even when DBSCAN gives no neighbour; spatial compactness supports object interpretation but cannot delete a small cortical/primordial candidate. After adaptive candidate-pool reclustering, call only cluster(s) with coherent non-ZP identity plus maternal/ooplasm enrichment, somatic anti-program separation and compatible object morphology. Reroute somatic-dominant clusters to Granulosa, Stromal/mesenchymal or follicular-interface review with an ambient/adjacent-oocyte state tag. Never expand Oocyte from zona signal, cortical location or a strict seed alone, and never report cellbin count as biological oocyte count.

## Annotate subtypes

Recluster only biologically coherent broad pools. In an R-first project use `scripts/run_seurat_pool_recluster.R`; use Scanpy or `scripts/run_sce_pool_recluster.R` only when the frozen expression object requires it. For same-batch Seurat cellbin inputs, keep the pool SCT model consistent with `references/seurat-cellbin-preprocessing.md`. For sheep ovary, adapt PCs, k and the selected value but always run the fixed formal grid `0.1,0.2,0.3,0.4,0.6` with `--resolution-contract sheep_ovary`; do not inherit the parent's selected value and do not introduce ultra-low candidates. Use `scripts/rank_pool_resolutions.py` only as a shortlist, then inspect marker/anti-marker completeness, adjacent-resolution migration and space before freezing. Preserve the complete lineage:

`input snapshot -> clustering run -> parent pool -> subcluster -> route -> decision -> final label`.

Close each completed pool. A closed membership cannot re-enter the same annotation round.

After every writeback run `scripts/plan_next_iteration.py --biological-profile ...`. Continue open-world lineage discovery, Route A broad-pool, Route B interface/RCTD applicability, Route C QC-pool then Atlas, and the Oocyte-specific route when a candidate exists until the correct biological profile returns `READY_FOR_COMPLETION_AUDIT`. A first-pass mapping or a completion gate run with the workflow profile is never terminal.

Do not trust an old PASS after changing a ledger or registry. Re-run state validation, the iteration planner and completion gate after the final job status/writeback; `autopilot_status.py` treats stale gates and release artifacts as unfinished work.

Before final assets, build one final annotation with `scripts/build_final_annotation.py`. Release a broad label only at moderate-or-higher confidence and a fine label only at high confidence. Preserve confidence, route and `fine_anchor_eligible` as metadata rather than duplicate strict/inclusive/display labels. Every accepted broad-only rescue enters its final broad-class DEG and dotplots; it cannot enter fine-marker discovery or receive a fabricated subtype.

Run `scripts/audit_release_taxonomy.py` on the current `state/cell_ledger.tsv.gz` with `--broad-column final_broad_label --status-column final_state --pool-column target_pool --out provenance/release_taxonomy_audit.json` and the active biological profile before the completion gate and final confirmation. Its biological broad census and retained-state census must be separate and the audit must pass.

If a prior annotation is supplied as a quality baseline, keep it blinded during fitting and writeback. After freezing the new biological decisions, run the full-feature broad-lineage audit for both versions and compare marker enrichment, anti-program leakage, rare-lineage safety, spatial morphology, unresolved accounting and rescue provenance. Do not call the new run better merely because it labels more observations, has more subtypes or agrees with the baseline. A requested improvement is achieved only when the predeclared class-specific acceptance table passes and no priority lineage materially regresses.

After the correct biological-profile completion gate and zero-open-incident gate pass, populate `state/annotation_support_registry.tsv` and run `scripts/build_confirmation_review_assets.R`. Then run `scripts/request_master_quality_review.py PROJECT_ROOT`; the main conversation Agent reviews the frozen biological result and records `PASS` (optionally with concerns) or `RETURN_FOR_ITERATION` using `scripts/record_master_quality_approval.py`. A return reopens targeted annotation work. A pass authorizes `scripts/build_confirmation_review.py`, followed by `scripts/request_final_annotation_confirmation.py PROJECT_ROOT`. Present the lightweight HTML plus the single final broad/fine census, retained QC/interface counts, rare calls and Atlas/RCTD returns. Do not generate final DEG, full broad/subtype tree dotplots, per-node/per-gene maps or the release HTML before user approval.

## Produce required outputs

Read `references/report-contract.md`. The final release must include all of the following:

- Broad-class DEG and UMAP/spatial maps for every release.
- High-confidence subtype DEG and subtype UMAP/spatial maps only when real fine labels pass; zero subtypes is valid.
- Per-node spatial highlights.
- **A broad-class tree dotplot, plus a separate subtype tree dotplot when high-confidence fine labels exist.**
- For every released dotplot level: canonical and data-specific panels when available, PNG, PDF and source TSV.
- Rebuild the broad data-specific marker panel from the **current final broad DEG** and the subtype panel from the **current high-confidence final fine-label DEG**, never from a prior release. Use `scripts/build_marker_panel_from_deg.py`; fail if any current broad/fine label lacks a marker group.
- Build release metadata with `scripts/prepare_report_metadata.py` and pass only `primary_broad_label` / `primary_subtype_label` to DEG, dotplot and map scripts. Keep `retained_state_display` in the separate uncertainty census; never turn a broad-only return into a synthetic subtype.
- Dotplot source columns retaining absolute detection and average expression.
- Display point size normalized within each gene to 0-100 and color scaled within each gene with a documented clip; never discard absolute values.
- Marker-gene spatial maps grouped by the cell type they support.
- Chinese or requested-language HTML report with expandable annotation tree, workflow, decisions, downloads and audit evidence.
- One final census and overview; route input/threshold/outcome tables; a Chinese phase-by-phase timeline reconstructed from workflow events; and an expandable raw state record at the bottom. Follow the successful Scanpy layout: overall annotated spatial map above the expandable tree, per-node highlight jumps, grouped marker spatial maps and the detailed workflow at the bottom.
- Separate the biological broad-class census/tree/DEG/dotplots from retained anatomical-interface, QC, technical and pending-state censuses; never make technical states look like cell types.
- Cell-level ledger, cluster decision ledger, pool/run registries, session information, manifests and checksums.

After user confirmation, the report and session information are final, run `scripts/build_release_manifest.py RELEASE_ROOT`, then use `scripts/audit_release.py RELEASE_ROOT` as a fail-closed gate. The manifest deliberately excludes its own files and the audit output to avoid circular hashes. Do not call a project complete unless the audit passes.

Run `scripts/check_completion_gate.py --biological-profile ACTIVE_BIOLOGICAL_PROFILE` before requesting confirmation. It must pass before release assets are authorized; a blocked gate returns to the queued iterations and no final HTML is generated.

## State rules

Read `references/state-schema.md` before writing labels.

- Never overwrite provenance columns.
- Version all decisions and preserve `supersedes` links.
- Mark completed pools `closed_and_frozen`.
- Record confidence separately from biological label.
- Distinguish cell/bin counts from inferred biological cell counts.
- Keep `defined_fine`, `defined_broad_only`, `interface_review`, `qc_holdout`, `technical_state`, `pending_review` and `excluded_initial_qc` distinct.
- Preserve legacy strict/inclusive/display fields only when resuming an old project. They are compatibility metadata, not release views. Migrate to one `final_*` annotation and do not display three competing results.
- Update state after every completed branch, not only at report time.

## Execution and testing

Use local execution for discovery and small audits. For submitted analysis read `references/job-orchestration.md` and follow submit → monitor → log inspection → repair/resubmit → artifact validation → state writeback. Before every scheduler submission run `scripts/preflight_generated_job.py` on every generated source and generate the visible name with `scripts/scheduler_job_name.py`. Record every failure in `provenance/incidents/incident_registry.tsv` before accepting its repair; any open incident blocks completion.

Every scheduler request must match measured or explicit algorithmic parallelism. A single Leiden optimization and an unwrapped `FindAllMarkers` loop are single-threaded; never reserve many idle CPUs for them. Parallelize Seurat candidate resolutions through the runner's `--resolution-workers`, or submit independent per-resolution evidence jobs with disjoint outputs and a dependency-gated aggregation step. Record worker count/backend/unit and audit CPU-time versus wall-time before copying any resource request to another sample.

Before distributing an updated Skill, run the standard validator and follow `references/testing.md`. Forward-test with a fresh agent that receives raw inputs and the Skill but not the intended clustering choice or annotation answer.
