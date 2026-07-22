# Efficient operation without lowering annotation quality

Use deterministic scripts and compact manifests to save Agent context, scheduler time and repeated large-object reads.

- The master reads only cohort/sample checkpoints, censuses, gate JSON and incident summaries. A sample worker alone reads its expression object and evidence tables.
- Maintain one compact `provenance/agent_handoff.json` per sample containing current phase, active job, open incidents, frozen hashes, next gate and user decisions. Do not reconstruct progress from all logs on every turn.
- Run the complete profile grid once. Reuse validated SCT/BANKSY PCA/UMAP and memberships by hash for downstream evidence; never refit them for a plotting, aggregation or syntax-only repair.
- Split compute by real parallel unit. Use one worker per concurrent resolution; give serial Leiden/DEG units one CPU. Aggregate only after dependencies pass.
- Inspect log tails and scheduler metrics first. Read full logs only when the tail cannot identify the failure.
- Before confirmation generate only one broad spatial PNG, one canonical broad marker dotplot PNG and a self-contained lightweight review HTML. Generate final DEG, full dotplots, per-node/per-gene spatial assets and release HTML once, after confirmation. Publish one final annotation rather than three duplicate views.
- Resolve scripts and profiles from the frozen Skill binding. Use the packaged preflights and validators instead of generating project-specific equivalents.
- Record reusable failure lessons in the incident registry; do not repeatedly rediscover an already closed root cause in another sample.
- Precompute and hash the Atlas feature transform, low-dimensional reference matrix/prototypes, broad crosswalk and approximate-neighbor index once. Project every query cell once; do not refit joint integration per sample or materialize dense query-by-reference distances.
- Perform all-cell concordance and cluster aggregation deterministically in tables. Spend full-feature DEG/pseudobulk/spatial compute and Agent interpretation only on frozen-QC writebacks, material broad disagreements and coherent OOD clusters.
