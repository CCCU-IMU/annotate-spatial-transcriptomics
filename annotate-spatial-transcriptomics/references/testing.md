# Validation and forward testing

## Static validation

Run the Skill validator, compile Python scripts, parse R scripts, validate templates and run synthetic unit tests.

Run `scripts/audit_release_taxonomy.py` on one passing and one deliberately failing fixture. The failing fixture must include a copied cohort/provenance identifier, a forbidden catch-all broad label or a retained state presented as a biological class.

Run `python -m unittest discover -s tests -v`. The release-contract suite must exercise the sheep-ovary resolver, exact StereoPy cellbin profile, GSE233801 priority, dotplot no-transfer rule, query-like held-out origin, legacy combined-tier rejection and one-worker-per-sample cohort isolation.

For v2, additionally test contract hash freshness, exact upstream BANKSY grid binding, complete all-candidate marker-family coverage, authoritative Atlas classwise calibration, release broad/fine hierarchy, typed residual-QC closure and read-only result-directory consistency. Run the result audit against representative historical projects; old `completion_gate=PASS` is expected to be blocked when the released ledger violates v2 semantics.

## End-to-end test

Use a real input directory without supplying the intended clustering choice or final labels. Require the agent to discover inputs, inventory prior progress, rank candidates, inspect shortlisted outputs, freeze a choice, write state and generate required artifacts.

## Leakage control

Do not pass prior final metadata, mapping tables or expected answer. Case-study documents may be unavailable during the forward test. Compare the agent's process and evidence, not exact agreement with a historical annotation.

## Acceptance gates

- No hard-coded sample paths in reusable scripts.
- No example-specific labels or thresholds in execution defaults.
- Every initial cluster decision has a hash-bound label-blind prelabel evidence artifact. It must compare every declared candidate lineage, bind positive and anti-DEG, record winner/runner-up margin and reject a moderate/high broad call based on fewer than two marker families or unresolved contradictions. Loading paper labels before this freeze fails.
- The same-batch Seurat cellbin test writes the frozen SCT+BANKSY profile to its manifest: entry QC `nCount>=100 AND nFeature>=75`, SCT v2/`glmGamPoi`, 4,000 variable features, at most 50,000 fitting observations, BANKSY `M=0`, `k_geom=30`, `lambda=0.2`, 30 PCs, Leiden `k_neighbors=50`, grid `0.2,0.4,0.6,0.8`, and UMAP `n_neighbors=30,min_dist=0.3,spread=1,n_epochs=300`.
- The Seurat cellbin runner never reuses imported StereoPy PCA/UMAP, never silently falls back from `glmGamPoi`, never auto-selects a final resolution and never auto-reassigns a small cluster.
- Sheep/Ovis/ovine/羊 plus ovary/ovarian/卵巢 and a full-feature Seurat RDS automatically resolves to R-first. A confirmed StereoPy cellbin_PPed conversion activates the exact fixed profile; unrecorded numerical drift and missing SHA256 support fail.
- Broad/targeted cohort SCT tests require query-only graph/Leiden/UMAP/DEG, explicit SCT v2 provenance, adaptive k/PC controls, the complete active grid and a fail-closed `glmGamPoi` dependency.
- Broad dotplots always exist with PNG, PDF and valid source TSV; subtype dotplots are required only when at least one real high-confidence fine label is released.
- Every dotplot asset has normalized and absolute PNG/PDF views from one source table; marker labels containing spaces remain one-column tree-order values.
- Every cell/observation is accounted for exactly once.
- State validator and release audit pass.
- The full release audit must run on every runtime declared compatible by `check_runtime.py`, including Python 3.8; do not rely on newer-only `pathlib` helpers for checksum path containment.
- A fresh agent can explain why it selected a clustering and which alternatives it rejected.
- A deliberately overbroad Oocyte or context-specific label is rejected until positive, negative, spatial-object and reclustering evidence pass; rarity alone never creates a label or route.
- A sheep-ovary Oocyte candidate must not be rejected or downweighted merely because it is cortical, subcortical, peripheral or at the section edge; small/primordial oocytes may be cortical. Location alone must still fail as positive evidence without the molecular, anti-program, spatial-object and reclustering gates.
- A sheep-ovary Oocyte route must send the complete predeclared multi-module starting gate to one query-only targeted cohort. Strict marker/anti-program seeds and compact spatial foci are supporting evidence, not the final census and not the only recluster membership; an isolated starting-gate candidate must remain in the cohort. Zona-only expansion and release without direct somatic return/QC handling must fail.
- In an R-first test, a readable full-feature Seurat RDS is selected as the primary backbone, existing Seurat clustering may be reused only after hash/membership validation, and historical labels remain hidden until the new ledger is frozen.
- An ovarian test cannot close a generic `Stromal/perivascular` super-class until generic stroma, mesenchymal-progenitor-like, mature smooth muscle, pericyte/mural, endothelial and steroidogenic-theca alternatives have machine-readable positive/anti-program and morphology audits. The test must not manufacture a standalone Mesenchymal or Smooth muscle class when their gates fail.
- A requested improvement over a baseline release must use a predeclared, baseline-blinded acceptance table and full-feature evidence comparison. Higher annotation rate, subtype count or baseline agreement alone cannot pass; no priority lineage may materially regress.
- A sheep-ovary R-first forward test must fail one-cluster-one-name output. Every fine label needs a reproducible functional/lineage program beyond its parent; unsupported historical subtypes are merged, and ECM/contractile/stress/low-RNA/anatomical differences remain state tags in the ledger and report.
- Every test must keep three machine-readable layers separate: literature candidate-lineage checklist, computational cohorts/QC state and release labels. Copying a cohort/QC identifier into a biological label fails. Missing a literature class after a documented negative audit passes; lowering its gate to complete a published taxonomy fails.
- When a matched single-cell reference is supplied, the test must preserve source labels, validate an explicit source-to-candidate crosswalk and enforce its transfer ceiling. Exact taxonomy copying fails; a dotplot-only reference claiming cell-level transfer fails; a fine or context-gated rare-lineage prediction that bypasses current-query marker/anti-marker/spatial evidence fails. A calibrated in-scope broad-only return to unlabeled frozen QC follows the separate state-aware rule below.
- A regression fixture containing steroidogenic, mature-contractile, generic-ECM, granulosa and endothelial programs must reject `Theca/follicular wall` as a broad catch-all, recover only coherent steroidogenic observations as `Theca`, and route mature-contractile observations through the Smooth-muscle gate.
- A standalone `Mesenchymal progenitor-like` call must fail when support is generic `DCN/LUM/COL1A1/PDGFRA` without stable S100A4/progenitor-like separation; a `Pericyte/mural` fine call under `Vascular-associated` must fail without its backbone and vascular adjacency. Both negative audits are acceptable terminal evidence.
- An ambiguous endothelial/pericyte or blood/lymphatic split must roll back to `Vascular-associated`; an ambiguous immune split must roll back to `Immune`. The test must prefer the least specific honest release name.
- Biological broad classes and anatomical-interface/QC/technical/pending states must have separate censuses. Retained states in broad biological DEG, canonical/data-specific dotplots or the biological annotation tree fail the release audit.
- A first-pass mapping with unfinished broad/targeted cohorts, interface, QC or open priority-lineage evidence creates a nonempty next-action queue and a blocked completion gate.
- A route named `anchor_assisted` without explicit query/anchor roles and query-only graph evidence is rejected.
- A final residual QC membership that is reclustered before Atlas fails the current direct-workflow contract.
- A v1.7 Atlas query that differs from the complete analysis set fails. Its frozen-QC submembership must equal terminal residual QC cell for cell; only that subset can be directly written back. Mapping a defined cell may create a challenge but cannot overwrite it without a closed orthogonal review.
- All-cell Atlas tests cover concordant, weak-challenge, material-disagreement, ontology-conflict, OOD, QC-writeback and QC-reject outcomes. Every material/OOD queue item must have exactly one evidence-bound decision; Atlas-only superseding decisions fail.
- Efficiency tests forbid dense query-by-reference distances and per-sample joint reference retraining in the default route. The fixed transform/reference/index manifest must be hash-bound and reusable.
- Atlas tests must verify nested cumulative calibration: every high row meets moderate-or-higher; mutually exclusive output counts distinguish `high`, `moderate_only` and `low_reject`; `moderate_or_higher_n = high_n + moderate_only_n`; and both accepted tiers remain broad-only. The default held-out target-precision targets are 0.90 and 0.95, not per-observation confidence cutoffs.
- Reference self-classification must fail as a final query-rescue calibration origin. Only disjoint query-like held-out current-query anchors with a PASS origin manifest can calibrate writeback; the legacy combined `medium_high` route remains diagnostic-only.
- Calibrated Atlas state-routing tests must directly return a moderate/high, non-OOD, ontology-compatible prediction to an unlabeled frozen-QC observation, retain profile-excluded classes, and send defined-label disagreement to review without overwriting the query label. Marker/spatial columns are audit/challenge evidence rather than duplicate per-cell gates. Separate uncalibrated multi-channel rescue tests still require genuinely independent channels.
- Every default broad candidate must expose at least two explicit non-overlapping positive families. A profile that stores all stromal or epithelial markers in one family must fail before annotation. Broad-presence tests must retain an abundant parent supported by absolute full-feature detection/pseudobulk even when centered module scores or one-vs-rest DEG are weak.
- Oocyte regression fixtures must include all canonical members of a passing targeted-recluster cluster even when only one is a strict seed or belongs to a seed-derived spatial object; context-only cells and explicitly hard-excluded cells remain outside.
- Residual QC at least 10% or 50,000 observations must block completion until a hash-bound upstream broad-recall audit closes.
- An interface without a machine-readable RCTD/reference applicability audit remains blocked; `not_applicable` requires a resolvable artifact and failed criteria.
- Low-confidence RCTD rejects are queued into the frozen QC holdout, never directly to Atlas; only canonical high confidence with independent evidence may support fine labels, while moderate confidence is broad-only.
- A QC Atlas attempt fails unless its query membership is exactly the complete terminal residual-QC membership after all broad/targeted cohorts. Reusing any defined broad/fine membership fails.
- The RCTD tier counts partition the complete query; high/moderate observations can return only under their evidence gates and otherwise reroute, every low observation reroutes, and RCTD evidence alone never creates a fine-anchor-eligible cell.
- A large or spatially diffuse terminal interface and a large first-pass direct label without cell-level purity/anchor review both block completion, even when cluster-level DEG or RCTD artifacts exist.
- Continuous lineage-signal coverage fixtures must fail when any whole-tissue/broad/targeted boundary omits a catalog-by-cluster product, when selected-plus-two-higher resolution coverage is incomplete, when an unexplained-program audit is missing, or when a positive-family signal is written as `absent`.
- A lineage below the initial naming threshold must remain `watch`; a later cross-lineage subcluster can support and reconstruct that broad class. Parent cohort labels may never reduce the candidate catalog.
- Project-boundary fixtures must reject project A derived expression as project B query evidence while allowing an explicitly registered Atlas/reference artifact.
- BANKSY broad-resolution fixtures must fail when selection is based on cluster count or omits full-catalog recall, zero-census and large-cluster purity review.
- Large-label dilution fixtures must preserve a spatially coherent 3–5% epithelial program as `candidate` even when its parent-cluster average is weak.
- Zero-census fixtures must block completion when a default sheep-ovary broad lineage lacks a query-derived multichannel negative audit; Atlas non-mapping is insufficient.
- Oocyte context fixtures must prevent context-only observations from entering the Oocyte census and must emit a somatic/pregranulosa candidate when appropriate evidence exists.
- Graph-sensitivity fixtures must reject a claimed rescue when only cluster count increases without improved core/support separation, anti-program clearance and spatial coherence.
- Cross-lineage fixtures must allow a supported epithelial subcluster inside a Stromal parent to create a direct Epithelial/mesothelial broad return.
- The pre-confirmation report must fail under framework v1.8 when any released broad/fine label lacks an all-analysis-set canonical-marker spatial panel, when its denominator differs from the analysis set, or when the panel was label-filtered.
- One final annotation covers the analysis set exactly once; the full-object ledger separately accounts for initial QC exclusions. Broad labels meet moderate-or-higher confidence and fine labels meet high confidence.
- If an analysis-scope policy exists, changing even one cell from `excluded_initial_qc` to `analysis_set`, leaking an excluded label into any view, or registering a stale membership hash must fail state validation. A report-only refresh must preserve pre-adjudicated view fields.
- The final HTML contains a broad/targeted cohort, direct-return, RCTD/Atlas dashboard and a Chinese raw-input-to-release event timeline.
- A one-cluster broad-class cohort after the complete grid is valid evidence and may close as `homogeneous_parent_confirmed`; it is not an execution error or a reason to force higher resolution.
- Synthetic cohort scenarios must cover: homogeneous parent return; a true full-feature/spatially supported subtype; ECM/stress/low-RNA/ribosomal state-only splits merged to parent; and a mixed lineage returned directly to its target.
- Empty JSON, status-only JSON, header-only TSV, stale evidence hashes, incomplete candidate grids and missing adjacent-resolution migration/ARI must fail the evidence schemas.
- Membership-partition fixtures must fail on duplicate direct returns, omitted cells, overlapping outcomes, an inexact targeted successor, a v1.7 Atlas query different from the analysis set, a frozen-QC submembership different from terminal residual QC, or an accepted/rejected writeback mismatch. Legacy v1.6 fixtures retain the residual-QC query contract.
- Repeated ranking with identical seed/config/input must be byte-deterministic for the ranking artifacts.
- Biological benchmarks report macro-F1, per-class precision/recall/F1, rare-type false positives, unresolved fraction and cross-sample stability. Ablations must show lower macro-F1 after removing anti-marker, spatial or broad-cohort review evidence.
- A held-out benchmark dataset is forbidden in reference, route and Atlas registries until annotation is hash-frozen. When GSE233801 is the held-out target, any GSE233801 reference occurrence fails before unblinding.
- Failed scheduler runs and logs remain registered; repaired executions use new run IDs.
- Scheduler regression fixtures reject an invalid stage name before submission, require AIP standard-input submission, and prevent completion-audit self-registration deadlock.
- Query-only cohort fixtures with zero anchors pass without anchor-only artifacts; positive-anchor fixtures still require and validate them.
- Membership readers accept both TSV and TSV.GZ; spatial morphology ignores duplicate coordinate copies in clustering exports; broad-family evidence uses the base sparse-matrix `sum()` generic.
- Atlas routing rejects an unbound combined map and accepts only a hash-bound exact disjoint target-plus-heldout union.
- Release manifests include `review/`; completion and release both reject empty workflow-event provenance, and release audit requires `release_sessionInfo.txt`.
- The HTML is visibly preliminary whenever the completion gate is absent or blocked.
- `autopilot_status.py` remains `CONTINUE` after a successful first pass, a stale completion gate, missing explicit final-annotation user confirmation, missing final assets, a stale report or a missing/full-release audit, and returns `COMPLETE` only after the confirmed current release is audited.
- Before confirmation, `autopilot_status.py` must require the completion gate, frozen broad spatial/canonical-marker assets, then a post-completion main-Agent biological quality approval, and only then a self-contained lightweight HTML. Requesting master approval after broad annotation alone or requesting user confirmation without the master approval fails; pre-confirmation output may not contain final DEG or full release assets.
- `build_report.py` and the full release audit must fail when final user confirmation is missing or when any confirmed ledger/gate/review hash is stale. Final DEG/full-dotplot/per-node/per-gene spatial/report generation is queued only after confirmation.
- The frozen confirmation request includes the single final census and retained QC/interface counts. Recording confirmation binds the exact cell, decision and support ledgers, completion/taxonomy/master-quality gates and lightweight-review hashes.
- A post-confirmation `final_*`/`release_*`/`report_*` run must still trigger `run_control` while nonterminal, but must not reopen state validation or the biological completion gate after it becomes terminal. A post-confirmation non-release run must trigger `confirmation_invalidation`.
- A fresh Agent continues across submitted jobs, failures, broad/targeted cohorts, direct returns, state writebacks and final release without asking for routine intermediate approval.
- A multi-sample fixture assigns exactly one full-workflow logical worker and one isolated root per sample, respects wave-limited parallelism, blocks active double claims without an audited takeover, and requires every sample to pass its own release audit. Worker mode may not reduce report or evidence gates.
