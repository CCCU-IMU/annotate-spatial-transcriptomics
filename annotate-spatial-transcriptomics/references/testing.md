# Validation and forward testing

## Static validation

Run the Skill validator, compile Python scripts, parse R scripts, validate templates and run synthetic unit tests.

Run `scripts/audit_release_taxonomy.py` on one passing and one deliberately failing fixture. The failing fixture must include a copied pool name, a forbidden catch-all broad label or a retained state presented as a biological class.

Run `python -m unittest discover -s tests -v`. The release-contract suite must exercise the sheep-ovary resolver, exact StereoPy cellbin profile, GSE233801 priority, dotplot no-transfer rule, query-like held-out origin, legacy combined-tier rejection and one-worker-per-sample cohort isolation.

## End-to-end test

Use a real input directory without supplying the intended clustering choice or final labels. Require the agent to discover inputs, inventory prior progress, rank candidates, inspect shortlisted outputs, freeze a choice, write state and generate required artifacts.

## Leakage control

Do not pass prior final metadata, mapping tables or expected answer. Case-study documents may be unavailable during the forward test. Compare the agent's process and evidence, not exact agreement with a historical annotation.

## Acceptance gates

- No hard-coded sample paths in reusable scripts.
- No example-specific labels or thresholds in execution defaults.
- The same-batch Seurat cellbin test writes the frozen SCT profile to its manifest: entry QC `nCount>=100 AND nFeature>=75`, SCT v2/`glmGamPoi`, 3,000 variable features, at most 50,000 fitting observations, 50 computed/30 used PCs, k=30 cosine Annoy with 50 trees, and the candidate grid `0.1,0.2,0.3,0.4,0.6`.
- The Seurat cellbin runner never reuses imported StereoPy PCA/UMAP, never silently falls back from `glmGamPoi`, never auto-selects a final resolution and never auto-reassigns a small cluster.
- Sheep/Ovis/ovine/羊 plus ovary/ovarian/卵巢 and a full-feature Seurat RDS automatically resolves to R-first. A confirmed StereoPy cellbin_PPed conversion activates the exact fixed profile; unrecorded numerical drift and missing SHA256 support fail.
- Pool SCT tests require joint query/anchor SCT/PCA with query-only graph/Leiden/UMAP/DEG, explicit SCT v2 provenance, adaptive k/PC/resolution controls and a fail-closed `glmGamPoi` dependency.
- Broad dotplots always exist with PNG, PDF and valid source TSV; subtype dotplots are required only when at least one real high-confidence fine label is released.
- Every cell/observation is accounted for exactly once.
- State validator and release audit pass.
- The full release audit must run on every runtime declared compatible by `check_runtime.py`, including Python 3.8; do not rely on newer-only `pathlib` helpers for checksum path containment.
- A fresh agent can explain why it selected a clustering and which alternatives it rejected.
- A deliberately overbroad Oocyte or context-specific label is rejected until positive, negative, spatial-object and reclustering evidence pass; rarity alone never creates a label or route.
- A sheep-ovary Oocyte candidate must not be rejected or downweighted merely because it is cortical, subcortical, peripheral or at the section edge; small/primordial oocytes may be cortical. Location alone must still fail as positive evidence without the molecular, anti-program, spatial-object and reclustering gates.
- A sheep-ovary Oocyte route must send the complete predeclared multi-module starting gate to query-only candidate-pool reclustering. Strict marker/anti-program seeds and compact spatial foci are supporting evidence, not the final census and not the only recluster membership; an isolated starting-gate candidate must remain in the pool. Zona-only expansion and release without cluster-level somatic rerouting must fail.
- In an R-first test, a readable full-feature Seurat RDS is selected as the primary backbone, existing Seurat clustering may be reused only after hash/membership validation, and historical labels remain hidden until the new ledger is frozen.
- An ovarian test cannot close a generic `Stromal/perivascular` super-class until generic stroma, mesenchymal-progenitor-like, mature smooth muscle, pericyte/mural, endothelial and steroidogenic-theca alternatives have machine-readable positive/anti-program and morphology audits. The test must not manufacture a standalone Mesenchymal or Smooth muscle class when their gates fail.
- A requested improvement over a baseline release must use a predeclared, baseline-blinded acceptance table and full-feature evidence comparison. Higher annotation rate, subtype count or baseline agreement alone cannot pass; no priority lineage may materially regress.
- A sheep-ovary R-first forward test must fail one-cluster-one-name output. Every fine label needs a reproducible functional/lineage program beyond its parent; unsupported historical subtypes are merged, and ECM/contractile/stress/low-RNA/anatomical differences remain state tags in the ledger and report.
- A sheep-ovary test must keep three machine-readable layers separate: literature candidate-lineage checklist, analysis parent pools and release labels. Copying a `_review`, `_candidate`, `_unresolved` or `_holdout` pool name into a biological label fails. Missing a literature class after a documented negative audit passes; lowering its gate to complete a published taxonomy fails.
- When a matched single-cell reference is supplied, the test must preserve source labels, validate an explicit source-to-candidate crosswalk and enforce its transfer ceiling. Exact taxonomy copying fails; a dotplot-only reference claiming cell-level transfer fails; a matched-reference prediction that bypasses current-query marker/anti-marker, spatial or rare-lineage gates fails.
- A regression fixture containing steroidogenic, mature-contractile, generic-ECM, granulosa and endothelial programs must reject `Theca/follicular wall` as a broad catch-all, recover only coherent steroidogenic observations as `Theca`, and route mature-contractile observations through the Smooth-muscle gate.
- A standalone `Mesenchymal progenitor-like` call must fail when support is generic `DCN/LUM/COL1A1/PDGFRA` without stable S100A4/progenitor-like separation; a standalone `Pericyte/mural` call must fail without its backbone and vascular adjacency. Both negative audits are acceptable terminal evidence.
- An ambiguous blood/lymphatic split must roll back to `Vascular/endothelial`; an ambiguous immune split must roll back to `Immune`. The test must prefer the least specific honest release name.
- Biological broad classes and anatomical-interface/QC/technical/pending states must have separate censuses. Retained states in broad biological DEG, canonical/data-specific dotplots or the biological annotation tree fail the release audit.
- A first-pass mapping with broad-only, interface, QC or open priority-lineage pools creates a nonempty next-action queue and a blocked completion gate.
- A route named `anchor_assisted` without explicit query/anchor roles and query-only graph evidence is rejected.
- A large post-clustering QC pool that has atlas mapping but no full QC-pool anchor-recluster remains blocked.
- Atlas tests must verify nested cumulative calibration: every high row meets moderate-or-higher; mutually exclusive output counts distinguish `high`, `moderate-only` and `low-reject`; `moderate_or_higher_n = high_n + moderate_only_n`; and both accepted tiers remain broad-only. The default held-out target-precision targets are 0.90 and 0.95, not per-observation confidence cutoffs.
- Reference self-classification must fail as a final query-rescue calibration origin. Only disjoint query-like held-out current-query anchors with a PASS origin manifest can calibrate writeback; the legacy combined `medium_high` route remains diagnostic-only.
- Multi-channel rescue tests must require current-query marker/anti-marker support plus an independent route/internal-anchor/spatial channel, calibrate only at support counts represented by query-like held-out anchors, reject Atlas-only calls, preserve high as a subset of moderate-or-higher, emit proposal-only broad labels and set `fine_anchor_eligible=false`.
- An interface without a machine-readable RCTD/reference applicability audit remains blocked; `not_applicable` requires a resolvable artifact and failed criteria.
- Medium/low-confidence RCTD rejects are queued into the frozen QC holdout, never directly to Atlas; only extreme confidence with independent evidence may support fine labels, while high confidence is broad-only.
- A QC Atlas attempt fails unless its query membership is exactly the residual-QC membership emitted by a validated earlier complete QC-anchor attempt. Reusing the full QC input snapshot, any biological-pool membership or any defined broad/fine membership fails.
- The RCTD tier counts partition the complete query; extreme/high observations can return only under their evidence gates and otherwise reroute, every medium/low observation reroutes, and the RCTD route itself never creates fine-anchor-eligible cells.
- A large or spatially diffuse terminal interface and a large first-pass direct label without cell-level purity/anchor review both block completion, even when cluster-level DEG or RCTD artifacts exist.
- One final annotation covers the analysis set exactly once; the full-object ledger separately accounts for initial QC exclusions. Broad labels meet moderate-or-higher confidence and fine labels meet high confidence.
- If an analysis-scope policy exists, changing even one cell from `excluded_initial_qc` to `analysis_set`, leaking an excluded label into any view, or registering a stale membership hash must fail state validation. A report-only refresh must preserve pre-adjudicated view fields.
- The final HTML contains a Route A–E dashboard and a Chinese raw-input-to-release event timeline.
- A single-cluster pool at low resolution is valid evidence, not an execution error or a reason to force higher resolution.
- Failed scheduler runs and logs remain registered; repaired executions use new run IDs.
- The HTML is visibly preliminary whenever the completion gate is absent or blocked.
- `autopilot_status.py` remains `CONTINUE` after a successful first pass, a stale completion gate, missing explicit final-annotation user confirmation, missing final assets, a stale report or a missing/full-release audit, and returns `COMPLETE` only after the confirmed current release is audited.
- Before confirmation, `autopilot_status.py` must require the completion gate, frozen broad spatial/canonical-marker assets, then a post-completion main-Agent biological quality approval, and only then a self-contained lightweight HTML. Requesting master approval after broad annotation alone or requesting user confirmation without the master approval fails; pre-confirmation output may not contain final DEG or full release assets.
- `build_report.py` and the full release audit must fail when final user confirmation is missing or when any confirmed ledger/gate/review hash is stale. Final DEG/full-dotplot/per-node/per-gene spatial/report generation is queued only after confirmation.
- The frozen confirmation request includes the single final census and retained QC/interface counts. Recording confirmation binds the exact cell, decision and support ledgers, completion/taxonomy/master-quality gates and lightweight-review hashes.
- A post-confirmation `final_*`/`release_*`/`report_*` run must still trigger `run_control` while nonterminal, but must not reopen state validation or the biological completion gate after it becomes terminal. A post-confirmation non-release run must trigger `confirmation_invalidation`.
- A fresh Agent continues across submitted jobs, failures, child pools, state writebacks and final release without asking for routine intermediate approval.
- A multi-sample fixture assigns exactly one full-workflow logical worker and one isolated root per sample, respects wave-limited parallelism, blocks active double claims without an audited takeover, and requires every sample to pass its own release audit. Worker mode may not reduce report or evidence gates.
