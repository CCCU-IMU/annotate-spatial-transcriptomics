# Sheep-ovary R-first Oocyte route

Use this reference whenever sheep/ovine ovarian spatial observations may contain Oocyte. It is a sanitized strategy distilled from a completed Seurat/R-first cellbin analysis; it contains no reusable sample label map or fixed resolution.

## Core distinction

Keep three sets separate:

1. **Starting candidate pool (recall set):** every observation passing the predeclared multi-module raw-count gate. This complete set enters query-only candidate-pool reclustering.
2. **Strict seeds and spatial foci (identity support):** high non-ZP/maternal evidence, low somatic contradiction and/or coherent spatial objects. These help identify an Oocyte-enriched recluster cluster but are not the final census and must not be the only recluster membership.
3. **Final broad Oocyte membership:** observations in a recluster cluster that passes program, anti-program and morphology adjudication. This set may be larger than the strict seeds and smaller than the starting candidate pool.

This separation prevents both failure modes: losing small/isolated true candidates before reclustering and calling zona-contaminated neighbouring somatic cellbins Oocyte.

## Required route

1. Screen raw/full-feature counts for non-ZP identity, maternal/ooplasm and zona modules plus Granulosa, stromal, vascular, epithelial and immune contradictory programs.
2. Build the starting pool with the profile's multi-module gate. Do not enlarge it from `ZP2/ZP3/ZP4`, location or morphology alone.
3. Calculate query-calibrated strict seeds. Spatially group candidates for object-level evidence, but retain isolated high-evidence candidates in the full starting pool.
4. Recluster the **complete starting pool** using a query-only graph. Select resolution adaptively: prefer the lowest complexity that separates an Oocyte-enriched program from coherent somatic alternatives and remains stable in adjacent settings.
5. Compare every candidate cluster for absolute detection and coherence of non-ZP identity/maternal modules, zona support, somatic anti-programs, QC complexity and spatial objects. A seed-rich cluster is still only a candidate until this comparison passes.
6. Freeze broad Oocyte only for passing cluster(s). Reroute Granulosa-dominant, ECM/stromal-dominant or other somatic clusters to their canonical parent/review pools and preserve `oocyte_adjacent`, `zona_ambient` or equivalent state tags.
7. Group spatially contiguous positive cellbins into putative objects for reporting. Report observation count and putative-object count separately; neither proves a histological oocyte count without image review.

## Spatial interpretation

- Cortical, subcortical, peripheral or section-edge location is never negative evidence. Small and primordial oocytes may occur in the ovarian cortex.
- Compact follicular morphology is supporting evidence, not a universal hard filter on candidate-pool membership.
- Location alone is never positive evidence and cannot relax the multi-module or somatic anti-program requirements.
- A singleton with strong multi-module evidence remains in the candidate recluster pool; a singleton does not become Oocyte automatically.

## Regression pattern from the completed R-first route

The completed route deliberately reclustered the full multi-marker candidate set. It separated an Oocyte-enriched cluster from a neighbouring ECM/stromal cluster despite substantial overlap in raw marker-hit counts. The decisive evidence was cluster-level non-ZP/maternal enrichment versus `MFAP4/LTBP4/IGFBP5/DCN`-like somatic coherence, not a higher marker-count cutoff. Therefore marker-hit number, strict-seed membership or spatial-focus membership must never replace full candidate-pool reclustering.

## Fail-closed outcomes

- If the starting pool is too small for defensible reclustering, retain `strict_oocyte_candidate` for manual/image review; do not widen thresholds from zona genes.
- If no cluster clears the somatic anti-program, return all candidates to somatic/interface pools with contamination tags and record a negative Oocyte audit.
- If morphology is unavailable, broad Oocyte may be proposed only when molecular separation is exceptionally coherent; mark morphology pending and require master/user review before release.
- Final Oocyte observations are never eligible to define stage-specific Oocyte subtypes without additional independent stage programs and morphology.
