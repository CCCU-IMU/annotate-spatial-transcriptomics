# v1.4.0 - Sheep-ovary R-first release contract

This release turns the validated sheep-ovary workflow into a deterministic, auditable profile while preserving adaptive biological decisions.

## Highlights

- Automatically resolves sheep/Ovis/ovine/羊 ovary projects with a full-feature Seurat RDS to R-first.
- Applies the frozen StereoPy `cellbin_PPed` SCT/PCA/neighbour profile only when matching conversion provenance is detected; unrecorded parameter drift fails.
- Handles the single-assay scalar JSON emitted by R `jsonlite`, verified against the real Seurat inspection schema rather than only a synthetic assay list.
- Makes GSE233801 the primary public adult-sheep somatic atlas when no usable matched count-level reference exists.
- Requires disjoint query-like held-out current-query anchors for Atlas rescue calibration; reference self-classification and the legacy combined `medium_high` tier cannot write back.
- Incorporates boundary and anti-overannotation evidence from Science 2025, Advanced Science 2026 and the AJOG 2025 adult-ovary review.
- Adds a master Agent plus one complete worker Agent per sample cohort architecture with isolated roots, ownership and release validation.
- Adds release-contract regression tests for workflow routing, fixed preprocessing, reference priority, calibration safety, broad-rescue inclusion and cohort isolation.

The literature and example workflow remain priors and regression evidence, never a sample-specific label map. Broad resolution, pool resolution, labels and optional subtypes remain adaptive to each query.
