# Annotate Spatial Transcriptomics v1.5.1

This release adds a post-completion biological quality gate for every sample.

- Pool reclustering, all applicable Route A–D work, rescue writeback, the single final annotation and the completion gate must finish before review.
- The main conversation Agent then performs a concise quality comparison with the bundled deidentified validated R-first sheep-ovary reference. Exact labels, counts and subtype depth are not answer keys.
- Approval evaluates broad-label reasonableness, marker/anti-marker plus spatial support and literature-candidate/confounded-lineage safety. It may pass with documented concerns; material problems return the sample worker to targeted iteration.
- The approval is hash-bound to the final ledgers, pools, routes, completion evidence, support registry and lightweight evidence assets. User confirmation and final release assets remain blocked until it passes.
- Multi-sample cohorts now expose explicit `READY_FOR_MASTER_QUALITY_GATE` and `MASTER_QUALITY_APPROVED` stages for each independent sample.
- An explicit `sheep_ovary_same_batch_rfirst` preset reuses the deidentified successful decision process while keeping resolutions, memberships and labels adaptive.
- The preset is open-world: a machine-readable, literature-derived catalog audits 14 default, evidence-dependent, stage-dependent and anatomy-dependent boundaries across all major ovarian lineage families. Negative results are valid, and coherent current-query programs outside the catalog create additional candidates.
- New projects do not use a generic rare-cell route. Oocyte alone retains an on-demand contamination-safe route; other supported context-specific candidates use ordinary biological pool validation.
