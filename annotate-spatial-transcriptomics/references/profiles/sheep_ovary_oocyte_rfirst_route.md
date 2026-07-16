# Sheep-ovary R-first Oocyte route

Use this contamination-safe targeted-cohort route whenever sheep ovarian spatial observations may contain Oocyte. It contains no reusable sample label map or fixed selected resolution.

## Three distinct memberships

1. **Starting targeted cohort (recall set):** every observation passing the predeclared multi-module full-count gate. The complete set enters query-only reclustering.
2. **Strict seeds/spatial foci (support):** high non-ZP/maternal evidence, low somatic contradiction and/or coherent spatial objects. These identify an enriched recluster cluster but are not the final census or sole cohort membership.
3. **Final broad Oocyte membership:** observations in recluster clusters passing program, anti-program and morphology adjudication. It may be larger than strict seeds and smaller than the starting cohort.

This prevents both loss of small/isolated true candidates and overcalling zona-contaminated neighbouring somatic cellbins.

## Required route

1. Screen full-feature counts for non-ZP identity, maternal/ooplasm and zona modules plus Granulosa, stromal, vascular, epithelial and immune anti-programs.
2. Build one immutable `oocyte_targeted_recluster` cohort from the multi-module gate. Do not enlarge it from `ZP2/ZP3/ZP4`, location or morphology alone.
3. Calculate query-calibrated strict seeds. Group candidates for object-level evidence, but retain isolated high-evidence candidates in the full cohort.
4. Recluster the complete cohort with a query-only graph. Select the candidate with the best integrated evidence for the Oocyte-versus-somatic question: preserve a stable Oocyte program or coherent somatic alternative first, reject state/technical fragmentation second, and use lower complexity only when candidates are otherwise equivalent.
5. Compare every cluster for absolute detection and coherence of non-ZP identity/maternal modules, zona support, somatic anti-programs, QC complexity and spatial objects. A seed-rich cluster remains only a candidate until this passes.
6. Freeze broad Oocyte only for passing clusters. Return Granulosa-dominant, ECM/stromal-dominant or other coherent somatic clusters directly to their supported broad/fine labels and preserve `oocyte_adjacent`/`zona_ambient` tags. Send non-interpretable remainder to QC. Do not create an intermediate cohort or automatically recluster a direct return again.
7. Group contiguous positive cellbins into putative objects for reporting. Observation count and putative-object count are separate; neither proves a histological oocyte count without image review.

## Spatial interpretation

- Cortical, subcortical, peripheral or section-edge location is never negative evidence; small/primordial oocytes may occur in the cortex.
- Compact follicular morphology supports identity but is not a hard filter on cohort membership.
- Location alone is never positive evidence and cannot relax molecular/anti-program gates.
- A singleton with strong multi-module evidence remains in the cohort but is not automatically Oocyte.

## Regression lesson

The deidentified R-first reference separated an Oocyte-enriched cluster from an ECM/stromal cluster despite overlap in raw marker-hit counts. Cluster-level non-ZP/maternal enrichment versus `MFAP4/LTBP4/IGFBP5/DCN`-like somatic coherence was decisive. Marker count, strict-seed membership or spatial focus must never replace full targeted-cohort reclustering.

## Fail-closed outcomes

- If the cohort is too small for defensible reclustering, record `not_applicable_reviewed` and retain evidence-positive observations for manual/image review; never widen from zona genes.
- If no cluster clears the somatic anti-program, return coherent somatic candidates directly to those lineages, retain the remainder in QC with contamination tags and record a negative Oocyte audit.
- If morphology is unavailable, broad Oocyte needs exceptionally coherent molecular separation and pending morphology must be visible to master/user review.
- Oocyte observations are not eligible for stage-specific subtypes without independent stage programs and morphology.
