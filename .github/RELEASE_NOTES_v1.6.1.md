# Annotate Spatial Transcriptomics v1.6.1

This stability release preserves the v1.6 direct-lineage design while making its biological evidence and membership closure machine-verifiable.

Highlights:

- every supported initial broad class still runs one complete broad-class reclustering cohort;
- verified same-batch sheep-ovary StereoPy cellbin R-first projects still run the full `0.1,0.2,0.3,0.4,0.6` grid;
- a homogeneous parent is a successful endpoint and is never forced into two clusters;
- terminal residual QC is frozen once and never reclustered before calibrated consensus rescue;
- cohort, direct-return and per-label evidence must pass content schemas;
- the final annotation, direct returns, cohort successors and Atlas partitions must close exactly cell for cell;
- new-project confidence uses only `low`, `moderate`, `high`; Atlas retains `high`, `moderate_only`, `low_reject`;
- leakage-safe held-out benchmarking reports macro/per-class metrics, unresolved fraction, rare false positives, cross-sample stability and evidence ablations.

Legacy persistent-pool and Route A–E artifacts remain readable only through migration tools under `references/legacy/` and `scripts/legacy/`.
