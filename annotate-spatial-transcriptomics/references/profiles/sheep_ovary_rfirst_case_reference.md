# Sanitized sheep-ovary R-first strategy reference

This strategy trace was distilled from a completed adult sheep-ovary spatial cellbin annotation and updated to the current direct-lineage architecture. It contains no sample identifier, private path, observation ID, final mapping table or expected cluster answer. Reproduce its safeguards, never its labels or parameters.

## 1. Input and preprocessing

A full-feature Seurat RDS converted from StereoPy `cellbin_PPed` was treated as a raw-count container. Imported reductions were provenance only. `Spatial` counts were entry-filtered with the same-batch contract, independently SCT-normalized, and used to compute fresh PCA, cosine Annoy neighbours, Leiden candidates and UMAP. The broad resolution was chosen from current adjacent-resolution stability, full-gene marker/anti-marker interpretability, cluster size and morphology.

## 2. Initial broad decisions

Every selected whole-tissue cluster underwent open-world lineage review. Supported clusters received direct moderate-or-higher initial broad labels. Candidate Theca, Smooth muscle, Pericyte/mural, epithelial, immune and Oocyte boundaries were audited independently. Low-information or irreducibly mixed observations entered QC directly; no computational membership name became a label.

Sheep immunoglobulin LOC aliases were resolved. A multi-locus immunoglobulin program plus an independent B/plasma regulator may support broad Immune despite sparse `PTPRC`; a single immunoglobulin locus, `JCHAIN`, CD74 or MHC signal does not.

## 3. Broad and targeted cohorts

Each supported initial broad class received one immutable query-only reclustering cohort. Subclusters were merged unless an independent functional/lineage program would be lost. Outcomes were direct parent-broad return, high-confidence shallow fine label, direct cross-lineage return, one decision-specific targeted cohort, or QC.

The Oocyte targeted cohort used the complete multi-module recall set while high-specificity marker/anti-program observations and spatial foci served only as seeds/support. Final Oocyte required multiple non-zona identity/maternal genes, somatic anti-program clearance and compatible object morphology. Zona-positive coherent Granulosa or stromal neighbours returned directly to those somatic lineages with ambient/adjacent tags.

Local interfaces used targeted reclustering first. RCTD was lower priority: canonical high plus independent evidence could support fine, moderate supported broad-only, and low entered QC. A cross-lineage return was not placed in a new intermediate cohort and was not automatically reclustered again.

## 4. Residual QC rescue

After every broad and targeted cohort closed, residual QC was frozen once and was not reclustered. The external Atlas was mapped over the complete analysis set for a single broad concordance audit. Current-query marker/anti-marker and internal-anchor or observed-density spatial channels gated direct QC writeback; `high` and `moderate_only` returns entered final broad membership with `fine_anchor_eligible=false`, while `low_reject`/OOD calls remained QC. Defined broad labels were challenged but never overwritten by Atlas alone. Atlas self-splitting was never accepted as query calibration.

## 5. Boundary safeguards

- Steroidogenic Theca required coherent enzyme/regulator programs and follicular outer-ring morphology; ECM or `LHCGR` alone remained structural stroma.
- Smooth muscle required a mature contractile core and coherent tracks; `ACTA2/TAGLN` alone was a state.
- Pericyte/mural required a multi-gene mural backbone and vascular adjacency; otherwise a negative audit was correct.
- An isolated neuronal-looking gene did not create Neural/Schwann or Neuroendocrine without full programs and compatible morphology.
- Graph subclusters driven by ECM, cortex, hypoxia, stress, cell cycle, low RNA or marker intensity were merged and retained as tags.

## 6. Release and reporting

The release published one annotation. Broad DEG and broad canonical/data-specific tree dotplots included every accepted observation, including direct and Atlas-rescued broad-only returns. Fine DEG/dotplots included only high-confidence real fine labels. QC/technical states were accounted for separately.

The Chinese HTML contained the annotated section map, expandable broad/fine tree, per-node spatial highlights, broad/fine tree dotplots, marker spatial maps, cohort/direct-return/Atlas evidence, workflow timeline, raw state and checksums. Main-Agent biological quality approval and user confirmation preceded expensive final assets.

## 7. Forward-test acceptance

A fresh Agent passes only when it independently reproduces these safeguards. It fails if it widens Oocyte from zona signal, creates a lineage from one gene, hides steroidogenic/contractile/mural alternatives, skips broad-class cohorts, reintroduces persistent biological pools, reclusters QC before Atlas, calibrates on Atlas self-classification, excludes rescued cells from broad DEG/dotplots, or substitutes many weak subtypes for a reliable shallow tree.
