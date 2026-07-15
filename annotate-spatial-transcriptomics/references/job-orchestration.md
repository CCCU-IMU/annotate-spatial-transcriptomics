# Scheduler and job orchestration

Use local execution only for discovery, schema checks and small summaries. Submit memory- or CPU-heavy expression work through the site's scheduler.

## Required lifecycle

1. Detect scheduler syntax, queues, environment paths and resource limits. Copy and edit a scheduler template; never assume example paths or queues.
2. Create a unique `run_id`, frozen parameter file and output directory before submission.
3. Generate the scheduler-visible name with `scripts/scheduler_job_name.py`, pass that exact name to the scheduler and record it as `scheduler_job_name`. Before submission run `scripts/preflight_generated_job.py` on the generated script and parameter/profile files; syntax or profile/grid failures stop before `csub`. Append the run to `state/run_registry.tsv` with status `submitted`, job ID, environment and logs. Supply a deterministic `work_key`, execution fingerprint, owner assignment and positive attempt. Only one `prepared/submitted/running` run may own a work key.
4. Monitor to a terminal state. Inspect both stdout and stderr and verify expected artifacts, not just scheduler status.
5. On failure, diagnose the concrete cause, preserve the failed run, and append it to `provenance/incidents/incident_registry.tsv` before creating a superseding run. Record the failure boundary, root cause, reusable artifacts, repair, verification and whether the lesson is generalizable. Never overwrite failed logs. The superseded run must already be terminally preserved. If a worker disappears while its scheduler job is active, transfer ownership and monitor that job; do not submit a duplicate.
6. On success, validate schemas, dimensions, cell-ID conservation, nonempty figures/tables and session information; then set status `validated_done`. Close the linked incident only after the repaired artifact passes its validation; `scripts/validate_incident_registry.py` must pass before completion.
7. Only validated outputs may be used for biological writeback.

After explicit final-annotation confirmation, continue registering DEG, map, dotplot, report, manifest, audit and package jobs, but name their stages with `final_*`, `release_*` or `report_*`. These release-only records are monitored normally and are included in the release manifest; they do not make the already hash-bound biological completion gate stale. A post-confirmation run outside those release-only stages is fail-closed evidence that biological work resumed and requires renewed validation and user confirmation.

Never start a downstream job while its upstream run is `submitted` or `running`. Large ledgers/metadata are written to a same-filesystem temporary path, checked for decompression, row count, unique IDs and schema, then atomically renamed. A filename appearing is not a completion signal; use `validated_done` plus the artifact checks. Scheduler arrays/parallel jobs may share read-only inputs but never the same output file.

Reference-assisted jobs must validate zero query/reference observation overlap immediately before submission or at job start. Do not assume an anchor table remains independent after a child pool is formed. Use `filter_reference_query_overlap.py` to freeze the repaired reference and its before/excluded/after label counts. Preserve the failed job, route and branch; submit the repaired execution with new IDs.

Resource choices follow object size and algorithm. More cores do not necessarily accelerate sparse R/Seurat steps. Record requested and observed resources. Do not leave memory unspecified unless the scheduler/site contract explicitly grants the intended allocation.

## Mandatory CPU-to-parallelism contract

Never request many scheduler CPUs for a stage that is actually single-threaded. Before submission, identify the stage's real parallel unit and record `requested_cpus`, `worker_count`, `parallel_backend` and `parallel_unit` in the parameter/manifest file.

- `leidenbase` does not make one Leiden optimization multithreaded. Seurat 5.5 can parallelize a vector of candidate resolutions through `future`; use `--resolution-workers min(number_of_resolutions, allocated_cpus)` with the Linux `multicore` backend when safe. The scheduler CPU request must equal the worker count plus only demonstrably active helper threads, never an arbitrary 64-core allocation.
- If a frozen graph must be evaluated by separate scripts, submit one read-only job per resolution and run those jobs concurrently. Give each truly single-threaded resolution job one CPU. Use a dependency-gated aggregation job for ARI/AMI and evidence summaries.
- `FindAllMarkers` in current Seurat releases iterates clusters serially. A job that calls it directly is single-threaded unless the wrapper explicitly parallelizes independent clusters/resolutions. Request one CPU, or split independent comparisons into scheduler jobs/fork workers with disjoint outputs and a validated aggregation step.
- `RunUMAP(..., umap.method="uwot")` uses `future::nbrOfWorkers()` as its thread count in Seurat 5.5. Record that count and prevent nested BLAS/OpenMP oversubscription.
- When memory duplication makes process parallelism unsafe, reduce both workers and requested CPUs together and record the reason. Reserving idle CPUs is not an acceptable memory strategy.

After completion, compare observed CPU time with wall time. A multi-core request whose CPU/wall ratio remains close to one is a resource-audit failure: preserve the run, correct future submissions and do not copy that allocation into templates.

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
| P40/P41 | COHORT/TARGET | broad-class or targeted cohort reclustering/validation; legacy `POOL/POOLQC` names remain readable |
| P50/P51 | RCTD/ATLAS | reference-assisted routes |
| P60 | RARE | legacy job-name compatibility only; do not create a generic rare-cell route |
| P61 | CONTEXT | open-world lineage audit or triggered Oocyte/context-specific validation |
| P70/P75 | WRITEBACK/CONFIRM | atomic writeback, lightweight confirmation spatial/dotplot assets, frozen support report and user gate |
| P80/P81/P82 | FINALDEG/DOTPLOT/SPATIAL | confirmed final release assets |
| P90/P99 | REPORT/AUDIT | report assembly and release audit |

Examples: `SAMPLE1__P10_SCT__A01`,
`SAMPLE1__P40_COHORT_stromal__A02`, and
`SAMPLE1__P81_DOTPLOT_broad__A01`. Names such as `sct_preprocess_v0`,
`validate_v001`, or a script filename are forbidden because they hide the
workflow phase. The detailed `run_id` may be longer, but it does not replace
the scheduler-visible name.
