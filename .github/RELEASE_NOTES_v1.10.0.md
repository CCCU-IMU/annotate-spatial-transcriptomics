# v1.10.0 — Broad-recall and state-aware Atlas correction

This release fixes a broad-annotation regression that could place abundant,
well-supported tissue lineages into QC when their markers were shared across a
continuum and therefore weak in centered module scores or one-vs-rest DEG.

Highlights:

- explicit two-family contracts for every default sheep-ovary broad lineage;
- absolute full-feature detection/pseudobulk as primary broad-presence evidence;
- direct broad-only return for calibrated moderate/high unlabeled QC Atlas calls;
- mandatory upstream recall audit for residual QC >=10% or >=50,000 observations;
- cluster-complete Oocyte materialization after strict cluster adjudication;
- broader epithelial cluster/component recall with two independent families;
- `Vascular-associated` broad parent for endothelial and pericyte/mural cells,
  with mature Smooth muscle retained as a separate broad lineage.
