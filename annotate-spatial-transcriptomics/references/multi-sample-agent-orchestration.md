# Master-agent and per-sample worker orchestration

Use this contract when a user requests more than one sample, parallel sample annotation, or a minimally supervised cohort run.

## Roles

- The main conversation Agent is the sole user-facing master. It owns the cohort control board, displays progress, reviews worker evidence, asks only material biological/release questions, performs cross-sample consistency audits and authorizes publication.
- Assign exactly one logical worker Agent to each sample. A worker owns that sample's complete workflow from input discovery through state, iterative routes and report assets. Parallelism never reduces evidence gates, route exhaustion, report content or annotation quality.
- An audit helper is not a sample worker. Do not create several Agents to annotate or audit the same sample unless the user explicitly requests an independent replication.

## Isolation and scheduling

- Every sample has a unique project root under `samples/<sample_id>/` or another explicitly registered absolute root. No worker writes another sample's state.
- Shared references and Skill files are read-only and registered by hash. Biological labels, anchors and thresholds never leak between samples.
- Requested parallelism controls how many sample workers may be active. When scheduler or Agent slots are fewer than samples, assign deterministic waves; never combine samples to fit a slot limit.
- A stopped worker with a live scheduler job is resumed by takeover and monitoring. It is not resubmitted until the existing job is proven terminal. Failed jobs and logs remain registered; repairs get a new run ID with `supersedes`.
- Every sample `run_registry.tsv` uses a unique active `work_key`, execution fingerprint, worker owner, positive attempt number and explicit superseded run. Existing v1.3 projects must run `scripts/migrate_project_v1_3_to_v1_4.py` before worker takeover.

## Lifecycle

`PLANNED -> WORKER_ASSIGNED -> ANALYSIS_RUNNING -> READY_FOR_MASTER_SAMPLE_GATE -> SAMPLE_FROZEN -> CROSS_SAMPLE_AUDIT -> COHORT_CONFIRMATION_PENDING -> RELEASE_RUNNING -> RELEASED`

The sample worker completes all computational/biological gates before `READY_FOR_MASTER_SAMPLE_GATE`, then sends a decision packet containing census, unresolved/interface/QC counts, rare calls, Atlas/RCTD returns, negative audits, route history and missing assets. The master independently verifies the packet and freezes each sample. After every sample is frozen, the master checks taxonomy semantics and evidence quality across samples without forcing identical labels. One aggregate user confirmation can authorize all unchanged frozen samples. Any later ledger/gate change invalidates that sample's confirmation.

After confirmation, the same sample worker may generate final inclusive DEG, both broad and subtype dotplots, maps and audited HTML. Each sample must independently reach `autopilot_status=COMPLETE` and `audit_release=PASS`; a cohort is never released because most samples passed.

## Cohort control files

Initialize with `scripts/init_annotation_cohort.py`, update assignments with `scripts/update_cohort_worker.py`, and validate with `scripts/validate_cohort_state.py`. After every sample is master-frozen and `provenance/cross_sample_audit.json` passes, run `scripts/request_cohort_confirmation.py`; bind the user's actual approval with `scripts/record_cohort_confirmation.py`. Any bound hash change invalidates release.

```text
COHORT_ROOT/
  config/cohort.json
  control/sample_manifest.tsv
  control/worker_registry.tsv
  control/sample_gate_registry.tsv
  control/cohort_run_index.tsv
  control/cohort_event_registry.tsv
  provenance/worker_packets/
  provenance/cross_sample_audit.json
  provenance/cohort_confirmation_request.json
  state/cohort_confirmation.json
  shared/
  samples/<sample_id>/
```

The master is the only writer under `control/`, `provenance/worker_packets/` and cohort-level `state/`. Sample workers write only their registered sample roots and their packet staging file.
