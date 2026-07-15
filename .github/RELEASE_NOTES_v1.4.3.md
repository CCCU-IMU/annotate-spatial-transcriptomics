# Annotate Spatial Transcriptomics v1.4.3

This patch release makes the sheep-ovary Seurat/R-first Oocyte route sensitive without widening the final identity gate.

## Two-tier Oocyte route

- Every observation passing the predeclared multi-module starting gate enters the complete query-only candidate recluster pool.
- High-specificity non-ZP/maternal markers, somatic anti-program clearance and compact spatial foci are strict seeds/support, not the final census or the only recluster membership.
- Isolated high-evidence candidates remain in the full pool; cortical, peripheral or section-edge location is never negative evidence.
- Final broad Oocyte requires an enriched recluster cluster with coherent non-ZP/maternal programs, somatic separation and compatible object morphology.
- Granulosa- or stromal-dominant neighbours return to their somatic parent/review pools with ambient/adjacent-oocyte state tags.
- Zona-only, location-only and strict-seed-only calls remain fail closed, and cellbin counts are not reported as histological oocyte counts.

## Regression protection

The rare-cell scripts now emit separate full-candidate and strict-support memberships. A release-contract test verifies that strict seeds cannot shrink the canonical recluster pool.
