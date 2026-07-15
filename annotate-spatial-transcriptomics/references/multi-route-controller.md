# Legacy pool-controller migration note

The historical Route A/B/C persistent-pool controller is not used for new projects. New annotation follows `direct-lineage-controller.md`.

When resuming an old project, preserve its pool and branch registries as immutable provenance, migrate current memberships into versioned `broad_class_recluster` or `targeted_recluster` cohorts, and convert terminal pool outcomes into direct-return records. Do not delete historical rows or silently reinterpret old Atlas evidence. If the old Atlas query depended on a QC-anchor child membership, preserve that provenance; new residual QC decisions use the current direct-query contract only after migration is complete.
