# Autonomous operation contract

Use this contract when the user asks the Agent to complete an annotation project with little or no intermediate supervision.

## Default authority

Within the supplied input/project scope, autonomously inspect files, choose candidate clusterings and pool-specific resolution grids, submit and monitor compute jobs, inspect logs, repair failed scripts with a new run ID, review biological evidence, write state, and generate the final release. Do not ask the user to approve routine intermediate cluster labels, resolution choices, pool creation, scheduler resubmission or report rendering.

Pause only when a missing fact would materially change the biological question or require new authority: no usable expression object, unknown species/tissue, an external reference or software installation requiring permission, or a destructive/external action outside the supplied scope. Record assumptions when stage/condition is unknown and prefer conservative labels.

## Continuous control loop

1. Run `autopilot_status.py PROJECT_ROOT` at startup and after every atomic writeback, repaired job, report build or release audit.
2. Execute the first emitted phase to completion. A biological phase may require Agent judgment and multiple submitted jobs; a file-existing check alone is insufficient.
3. Validate artifacts before changing the run from `submitted` to `validated_done`.
4. Run `plan_next_iteration.py` after every biological writeback. Execute the entire queue, including child pools created by the current route.
5. Re-run state validation and the completion gate after the last **biological-analysis** run registry update. Never trust a PASS artifact older than a biological ledger or registry before confirmation.
6. Build strict, inclusive and display views. After the biological completion gate passes, freeze the census with `request_final_annotation_confirmation.py` and pause for explicit user approval. Do not generate final DEG, dotplot, spatial or HTML assets before that approval.
7. Bind the approval to the current ledger/gate hashes with `record_final_annotation_confirmation.py`, then generate final assets and the report. Post-confirmation runs must be release-only stages (`final_*`, `release_*` or `report_*`): they remain subject to run monitoring and release checks, but their run-registry entries do not reopen the frozen biological gate. Any later cell/cluster ledger or completion-gate change, or any post-confirmation non-release analysis run, invalidates the approval and returns to step 6.
8. Rebuild the release manifest after every report/state/asset change, then run the full release audit.
9. Stop only when `autopilot_status.py` returns `COMPLETE` and `release_audit.json` is PASS.

## Forbidden early exits

Do not stop because:

- the original graph clusters have tentative labels;
- the broad annotation fraction is high;
- one parent-pool reclustering completed;
- RCTD or atlas mapping produced labels;
- QC holdout shrank but remains unreviewed;
- a report can be opened;
- all scheduler jobs are DONE while their artifacts are unvalidated;
- the user did not explicitly say “continue”.

These are intermediate states. Continue autonomously while an in-scope safe next action exists.

## Biological judgment contract

Use profiles as priors and safety gates, never as an answer key. Derive global and pool resolution from the current data. Use the lowest stable resolution that preserves coherent compartments; roll back state-only fragmentation. Allow broad-only and retained interface/QC outcomes when fine evidence fails. For rare or contamination-prone identities, require the profile’s multi-module, anti-program, spatial-object and strict-recluster gates even if a reference prediction is confident.

## Failure recovery

Preserve failed/cancelled jobs, logs, route attempts and invalid decisions. Repair under a new run/route/decision ID, link `supersedes`, revalidate the query/reference boundary and update the state only after the repaired artifacts pass. Never overwrite failure evidence to make the control board appear complete.

## Final handoff

Return the audited report, cell/cluster ledgers, route/branch/run registries, strict/inclusive/display metadata, DEG and figure source tables, environment/session records, manifest, checksums, release audit and the reusable Skill package. Distinguish biological cell counts from spot/cellbin observations in every summary.
