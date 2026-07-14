# Scheduler and job orchestration

Use local execution only for discovery, schema checks and small summaries. Submit memory- or CPU-heavy expression work through the site's scheduler.

## Required lifecycle

1. Detect scheduler syntax, queues, environment paths and resource limits. Copy and edit a scheduler template; never assume example paths or queues.
2. Create a unique `run_id`, frozen parameter file and output directory before submission.
3. Generate the scheduler-visible name with `scripts/scheduler_job_name.py`, pass that exact name to the scheduler and record it as `scheduler_job_name`. Append the run to `state/run_registry.tsv` with status `submitted`, job ID, environment and logs. Supply a deterministic `work_key`, execution fingerprint, owner assignment and positive attempt. Only one `prepared/submitted/running` run may own a work key.
4. Monitor to a terminal state. Inspect both stdout and stderr and verify expected artifacts, not just scheduler status.
5. On failure, diagnose the concrete cause, preserve the failed run, create a superseding run with the fix and resubmit. Never overwrite failed logs. The superseded run must already be terminally preserved. If a worker disappears while its scheduler job is active, transfer ownership and monitor that job; do not submit a duplicate.
6. On success, validate schemas, dimensions, cell-ID conservation, nonempty figures/tables and session information; then set status `validated_done`.
7. Only validated outputs may be used for biological writeback.

After explicit final-annotation confirmation, continue registering DEG, map, dotplot, report, manifest, audit and package jobs, but name their stages with `final_*`, `release_*` or `report_*`. These release-only records are monitored normally and are included in the release manifest; they do not make the already hash-bound biological completion gate stale. A post-confirmation run outside those release-only stages is fail-closed evidence that biological work resumed and requires renewed validation and user confirmation.

Never start a downstream job while its upstream run is `submitted` or `running`. Large ledgers/metadata are written to a same-filesystem temporary path, checked for decompression, row count, unique IDs and schema, then atomically renamed. A filename appearing is not a completion signal; use `validated_done` plus the artifact checks. Scheduler arrays/parallel jobs may share read-only inputs but never the same output file.

Reference-assisted jobs must validate zero query/reference observation overlap immediately before submission or at job start. Do not assume an anchor table remains independent after a child pool is formed. Use `filter_reference_query_overlap.py` to freeze the repaired reference and its before/excluded/after label counts. Preserve the failed job, route and branch; submit the repaired execution with new IDs.

Resource choices follow object size and algorithm. More cores do not necessarily accelerate sparse R/Seurat steps. Record requested and observed resources. Do not leave memory unspecified unless the scheduler/site contract explicitly grants the intended allocation.

The Agent may repair software/path/read-format problems autonomously inside the authorized project. Biological ambiguity is not a job failure; route it through the iterative controller.

## Mandatory scheduler-visible naming

Every submitted job uses `SAMPLE__Pnn_STAGE[_SCOPE]__Ann`, generated and
validated by `scripts/scheduler_job_name.py` (default maximum 48 characters).
Put sample and stage first so they remain visible when a scheduler UI truncates
the right edge. `SCOPE` is a short pool/release target, not an arbitrary script
name. `Ann` is the execution attempt and increments after a preserved failure.

| Code | Scheduler stage | Meaning |
|---:|---|---|
| P00 | INPUT | input inspection/snapshot |
| P10/P11 | SCT/SCTQC | SCT preprocessing and its validation |
| P20/P21 | RESGRID/RESEVID | resolution grid and candidate evidence |
| P30 | BROAD | broad-class evidence/adjudication compute |
| P40/P41 | POOL/POOLQC | pool reclustering and validation |
| P50/P51 | RCTD/ATLAS | reference-assisted routes |
| P60 | RARE | strict rare-lineage route |
| P70/P75 | WRITEBACK/CONFIRM | atomic writeback and frozen user gate |
| P80/P81/P82 | FINALDEG/DOTPLOT/SPATIAL | confirmed final release assets |
| P90/P99 | REPORT/AUDIT | report assembly and release audit |

Examples: `SAMPLE1__P10_SCT__A01`,
`SAMPLE1__P40_POOL_stromal__A02`, and
`SAMPLE1__P81_DOTPLOT_broad__A01`. Names such as `sct_preprocess_v0`,
`validate_v001`, or a script filename are forbidden because they hide the
workflow phase. The detailed `run_id` may be longer, but it does not replace
the scheduler-visible name.
