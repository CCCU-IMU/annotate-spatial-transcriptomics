# Sheep-ovary annotation evidence added in 2025-2026

Use this reference only after the biological context resolves to sheep ovary. It supplements, but never replaces, the current-query full-feature evidence and the sheep-specific sources in `sheep_ovary_evidence.md`.

## Evidence hierarchy

1. Current-query full-feature marker, anti-marker, stability and spatial evidence remains authoritative.
2. A matched count-level sheep single-cell object is the strongest external channel after auditing its labels and provenance.
3. Without such an object, GSE233801 is the primary public adult-sheep somatic atlas. Its scope is mainly Granulosa, Stromal/mesenchymal, Vascular/endothelial, Pericyte/mural candidates and Immune; it is not a sufficient source for Oocyte, steroidogenic Theca or surface Epithelium.
4. Sheep developmental atlases and the 2026 sheep-human atlas provide complementary lineage and state evidence.
5. Human/mouse ovary studies provide conserved boundary, anti-marker and morphology hypotheses only. They cannot transfer cells or require the sheep section to contain the same classes or subtypes.

## Gaylord et al., Science 2025

Gaylord et al., *Comparative analysis of human and mouse ovaries across age*, Science (2025), DOI `10.1126/science.adx0659`.

- The whole-ovary broad taxonomy contains Oocyte, Granulosa, Theca, Fibroblast, Smooth muscle, Pericyte, Epithelia, Endothelia, Immune and Glia; mouse additionally contains Luteal cells in the sampled state.
- Cross-species broad identities are more conserved than fine identities. Granulosa, fibroblast and endothelial subtypes were relatively conserved, whereas theca, pericyte and epithelial subtypes showed stronger species specificity.
- Ovarian glia are biologically plausible but rare and were supported by a coherent `S100B/SOX10`-like glial program and nerve-associated morphology. A neuronal-looking single gene is not sufficient.
- Theca requires a steroidogenic program such as `CYP17A1`, not follicular-wall location or ECM alone. Pericytes require a mural program and neurovascular/vascular adjacency.
- The paper deliberately resolved many subtypes for a comparative aging question. A single sheep spatial section should borrow the broad boundaries, not the paper's full subtype depth.

## Zhao et al., Advanced Science 2026

Zhao et al., *Comparative Single-Cell Transcriptomic Atlas Reveals the Genetic Regulation of Reproductive Traits*, Advanced Science (2026), DOI `10.1002/advs.202517633`.

- The atlas profiles 316,616 sheep cells across 15 reproductive/CNS tissues and integrates 769,901 human cells. It reports strong conservation of many major lineage programs, but also species-specific expression and uneven tissue representation.
- The ovary-focused analysis supports sheep Granulosa and its preantral, antral and atretic trajectories. Those names are candidate functional states, not mandatory release labels for a whole-ovary cellbin section.
- The cross-tissue atlas distinguishes fibroblast, smooth muscle, pericyte, endothelial, epithelial, immune, neuronal and glial categories. Because these are atlas-wide classes, a sheep ovary query must still establish ovary-local marker and spatial evidence before releasing a rare or contractile class.
- Cross-species classifiers in the publication do not license uncalibrated label transfer. For this Skill, any external mapping remains broad-only and must be calibrated on query-like held-out current-query anchors.

## Rooda et al., AJOG 2025 expert review

Rooda et al., *The adult ovary at single cell resolution: an expert review*, American Journal of Obstetrics & Gynecology 232 (2025), DOI `10.1016/j.ajog.2024.05.046`.

- Recurrently expected adult ovarian classes are Oocyte, Granulosa, Theca, Stroma, Endothelial, Smooth muscle, Perivascular and Immune; surface epithelium may be under-captured because it is a thin fragile layer.
- Cortex and medulla contain related broad somatic classes but differ in subtype composition and morphology; cortex is dense and follicle-bearing, while medulla is looser and vascular/nerve-rich.
- Sampling and filtration can remove large growing oocytes, and absence of a lineage in a single-cell reference is not proof of absence in a spatial section.
- The review explicitly warns that one broad type can be divided into effectively endless computational subtypes and that single-marker annotation is weaker than holistic transcriptome evidence. Therefore graph clusters, ECM amount, stress, location and marker intensity remain state tags unless an independent biological program passes.
- Oogonial-stem-cell-like claims require exceptional evidence; reported isolation procedures can instead capture perivascular cells. This Skill does not offer an oogonial stem cell release label.

## Machine-actionable consequences

- Candidate broad checklist: Oocyte, Granulosa, Theca, Stromal/mesenchymal, Smooth muscle, Pericyte/mural, Vascular/endothelial, Immune, Epithelial/mesothelial; Luteal and Neural/Schwann are context-gated.
- This checklist is not a quota. A documented negative audit is a successful result.
- `ACTA2/TAGLN` without `MYH11/CNN1/ACTG2` and coherent tracks remains contractile stroma/myofibroblast state.
- `RGS5/PDGFRB/CSPG4/NOTCH3/MCAM` plus vascular adjacency is required for Pericyte/mural.
- `DLG2`, `RBFOX1`, `TENM3` or another isolated neural gene cannot establish Neural/Schwann or Neuroendocrine. Require a multi-gene resident program, anti-program clearance and track-like morphology; neuroendocrine additionally requires a secretory-neuroendocrine backbone.
- Zona transcripts in Granulosa do not establish Oocyte. Require multiple non-ZP oocyte identity/maternal-ooplasm modules and object-level morphology.
- The primary broad DEG and dotplots use the inclusive biological membership, including every calibrated broad rescue. Strict cells are an additional sensitivity cohort, not a replacement.
