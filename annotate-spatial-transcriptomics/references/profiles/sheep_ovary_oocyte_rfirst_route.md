# Sheep-ovary Oocyte canonical-cluster route

Use this rule whenever a second-round sheep-ovary cohort contains a plausible Oocyte program. It contains no sample label map or fixed selected resolution.

## Evidence layers

1. **Non-ZP identity:** ooplasm/maternal and oocyte-regulatory programs, not zona genes alone.
2. **Zona support:** ZP2/ZP3/ZP4 supports a coherent identity but is ambient-prone and never seeds the census alone.
3. **Somatic competitors:** Granulosa, Stromal, Vascular, Epithelial and Immune programs.
4. **Object/ensemble morphology:** compact follicular topology supports a subcluster-level call but is not a per-cell admission gate.
5. **Spatial context window:** nearby observations are evidence for pregranulosa/cumulus or ambient contamination; they cannot become Oocyte by adjacency.

## Canonical second-round rule

At every second-round cohort resolution, scan all subclusters for non-ZP identity, maternal/ooplasm support, zona support, somatic anti-programs and morphology. Select the cohort resolution using the same full-grid and neighboring-resolution rule as other lineages.

Strict seeds/spatial foci are discovery aids and never the final census.

Within the complete cohort, when a second-round Leiden subcluster has a coherent non-ZP/maternal identity, somatic clearance, neighboring-resolution stability and compatible object-level morphology, return the complete subcluster as broad Oocyte except observations with a direct hard contradiction or objective input QC failure. Do not require each member to be a strict seed, a seed-derived spatial object or individually above a two-family threshold. This preserves sparse canonical members.

If the subcluster also contains a separable somatic identity core, trigger the standard local mixed-subcluster route. Return coherent Granulosa/ECM/vascular/epithelial components to those identities and keep `oocyte_adjacent` or `zona_ambient` as state tags. Nonselection never becomes QC automatically.

If no canonical subcluster exists, retain a watch or negative Oocyte audit. Do not create an all-object per-cell Oocyte classifier and do not widen from zona, location or morphology alone.

## Zero-census canonical challenger

When ordinary second-round scanning returns no Oocyte, screen the complete analysis set once with the same label-blind multi-module starting rule. If that screen is nonzero, form exactly one query-only targeted cohort containing **all** starting-gate observations across initial clusters and recluster it from project-local raw counts. Existing broad/fine labels, Atlas labels, prior repair membership and spatial position remain invisible until the candidate cluster and exclusions are frozen.

This route is a cohort-level challenger, not a whole-object per-cell classifier. A passing canonical cluster returns all of its members except objective input-QC failures or direct multigene somatic contradictions spanning independent families. Strict seeds, zona signal, spatial foci and putative object IDs support the cluster decision but cannot shrink the passing cluster. If the complete starting gate is empty or its targeted cohort has no stable canonical cluster, record a justified negative census.

## Pregranulosa/cumulus context

A context window may identify adjacent pregranulosa/cumulus candidates but does not expand Oocyte. Require FOXL2 plus stage-compatible KITL/WNT4/RSPO1/LGR5 evidence, follicular support such as FST/GJA1/CDH2/INHBB/SERPINE2, somatic competitor clearance and compatible topology. Proximity, FOXL2 alone or additional graph clusters are insufficient.

## Reporting

Group contiguous Oocyte-positive cellbins into putative objects only for visualization. Cellbin count and putative-object count are separate and neither equals histological oocyte count without image review. No stage-specific Oocyte subtype is released without an independent stage program and morphology.
