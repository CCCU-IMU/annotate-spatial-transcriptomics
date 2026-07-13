# Validation and forward testing

## Static validation

Run the Skill validator, compile Python scripts, parse R scripts, validate templates and run synthetic unit tests.

## End-to-end test

Use a real input directory without supplying the intended clustering choice or final labels. Require the agent to discover inputs, inventory prior progress, rank candidates, inspect shortlisted outputs, freeze a choice, write state and generate required artifacts.

## Leakage control

Do not pass prior final metadata, mapping tables or expected answer. Case-study documents may be unavailable during the forward test. Compare the agent's process and evidence, not exact agreement with a historical annotation.

## Acceptance gates

- No hard-coded sample paths in reusable scripts.
- No example-specific labels or thresholds in execution defaults.
- Broad and subtype dotplots both exist with PNG, PDF and valid source TSV.
- Every cell/observation is accounted for exactly once.
- State validator and release audit pass.
- The full release audit must run on every runtime declared compatible by `check_runtime.py`, including Python 3.8; do not rely on newer-only `pathlib` helpers for checksum path containment.
- A fresh agent can explain why it selected a clustering and which alternatives it rejected.
- A deliberately overbroad rare-cell label is rejected until positive, negative, spatial-object and reclustering evidence pass.
- A first-pass mapping with broad-only, interface, QC or open priority-lineage pools creates a nonempty next-action queue and a blocked completion gate.
- A route named `anchor_assisted` without explicit query/anchor roles and query-only graph evidence is rejected.
- A large post-clustering QC pool that has atlas mapping but no full QC-pool anchor-recluster remains blocked.
- Atlas tests must verify nested cumulative calibration: every high row meets moderate-or-higher; mutually exclusive output counts distinguish `high`, `moderate-only` and `low-reject`; `moderate_or_higher_n = high_n + moderate_only_n`; and both accepted tiers remain broad-only. The default held-out target-precision targets are 0.90 and 0.95, not per-observation confidence cutoffs.
- Multi-channel rescue tests must require current-query marker/anti-marker support plus an independent route/internal-anchor/spatial channel, calibrate only at support counts represented by query-like held-out anchors, reject Atlas-only calls, preserve high as a subset of moderate-or-higher, emit proposal-only broad labels and set `fine_anchor_eligible=false`.
- An interface without a machine-readable RCTD/reference applicability audit remains blocked; `not_applicable` requires a resolvable artifact and failed criteria.
- Medium/low-confidence RCTD rejects remain queued for calibrated atlas/internal-anchor fallback; only extreme confidence with independent evidence may support fine labels, while high confidence is broad-only.
- The RCTD tier counts partition the complete query; extreme/high observations can return only under their evidence gates and otherwise reroute, every medium/low observation reroutes, and the RCTD route itself never creates fine-anchor-eligible cells.
- A large or spatially diffuse terminal interface and a large first-pass direct label without cell-level purity/anchor review both block completion, even when cluster-level DEG or RCTD artifacts exist.
- Strict/inclusive/display views each cover the analysis set exactly once; the full-object ledger separately accounts for initial QC exclusions.
- If an analysis-scope policy exists, changing even one cell from `excluded_initial_qc` to `analysis_set`, leaking an excluded label into any view, or registering a stale membership hash must fail state validation. A report-only refresh must preserve pre-adjudicated view fields.
- The final HTML contains a Route A–E dashboard and a Chinese raw-input-to-release event timeline.
- A single-cluster pool at low resolution is valid evidence, not an execution error or a reason to force higher resolution.
- Failed scheduler runs and logs remain registered; repaired executions use new run IDs.
- The HTML is visibly preliminary whenever the completion gate is absent or blocked.
- `autopilot_status.py` remains `CONTINUE` after a successful first pass, a stale completion gate, missing explicit final-annotation user confirmation, missing final assets, a stale report or a missing/full-release audit, and returns `COMPLETE` only after the confirmed current release is audited.
- `build_report.py` and the full release audit must fail when final user confirmation is missing or when any confirmed ledger/gate hash is stale. Final DEG/dotplot/spatial/report generation is queued only after confirmation.
- The frozen confirmation request must include strict unresolved/QC/interface counts and strict/inclusive/display broad censuses. Recording confirmation must bind the exact cell-ledger, decision-ledger and completion-gate hashes; changing any one invalidates release generation.
- A post-confirmation `final_*`/`release_*`/`report_*` run must still trigger `run_control` while nonterminal, but must not reopen state validation or the biological completion gate after it becomes terminal. A post-confirmation non-release run must trigger `confirmation_invalidation`.
- A fresh Agent continues across submitted jobs, failures, child pools, state writebacks and final release without asking for routine intermediate approval.
