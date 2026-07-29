# Annotate Spatial Transcriptomics Skill v2.5.0

v2.5.0 is the stable second-round-first, evidence-bound annotation release.

- Uses the first whole-tissue partition only to build one cohort per initial cluster; formal identity is assigned from project-local raw-count second-round SCT/PCA/SNN/Leiden analysis.
- Releases Endothelial, Pericyte/mural and mature nonvascular Smooth muscle as independent competing identities; `Vascular-associated` is legacy-only.
- Runs one calibrated all-cell Atlas route with independent 90% moderate and 95% high precision validation, direct broad-only rescue for eligible unlabeled observations and no silent overwrite of defined labels.
- Reviews every context-evaluable broad as one whole-query precision/recall task with raw-count marker families, full-transcriptome pseudobulk, pairwise competitors and spatial consistency.
- Limits observation splitting to genuine mixed second-round subclusters, audits combined split workload before dispatch and keeps sub-percent stable traces as watch signals.
- Preserves canonical Oocyte and follicle-ROI biological review, exact cell-ID writeback, source-bound repair provenance and a single public `final_cell_type`.
- Stops after the bounded Atlas/unresolved/ROI/per-broad recovery budget: high residual uncertainty becomes a non-release review candidate instead of reopening global residual or QC-anchor reclustering.

The project framework schema remains 2.0.0. The internal controller/artifact protocol remains 2.2.0 for compatibility with frozen projects and reusable derived partitions.
