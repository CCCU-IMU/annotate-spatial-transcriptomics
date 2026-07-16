# Annotate Spatial Transcriptomics v1.7.0

This release makes broad annotation evidence-first and extends calibrated Atlas assistance into a safe, efficient all-cell concordance audit.

## Highlights

- Freeze label-blind positive/anti evidence for every candidate broad lineage before paper context or initial labels.
- Map the complete analysis set once after terminal QC freeze, while restricting direct Atlas writeback to the frozen QC subset.
- Detect material broad disagreements, ontology conflicts and coherent out-of-distribution groups without allowing Atlas-only label overwrite.
- Require one query-evidence review for each triggered cluster/cohort and preserve `Unknown candidate` as an honest open-set outcome.
- Reuse a fixed Atlas transform/low-dimensional index; avoid dense pairwise distances, repeated joint integration and whole-object RCTD.
- Validate exact coverage, QC partitioning, review closure, artifact hashes and broad-only/fine-anchor ceilings with new deterministic scripts.

Projects created under v1.6.1 can add the new registry fields with `migrate_project_v1_6_1_to_v1_7_0.py`. The migration intentionally does not invent missing prelabel evidence or Atlas review decisions.
