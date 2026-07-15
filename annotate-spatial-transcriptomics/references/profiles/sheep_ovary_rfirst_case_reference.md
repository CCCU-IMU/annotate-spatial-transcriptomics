# Sanitized sheep-ovary R-first forward-test reference

This is a reusable strategy trace distilled from a completed adult sheep-ovary spatial cellbin annotation. It intentionally contains no sample identifier, private path, observation ID, final mapping table or expected cluster answer. Use it to reproduce the decision process, never to copy labels.

## 1. Input and preprocessing

The successful path began from a full-feature Seurat RDS converted from StereoPy `cellbin_PPed`. The converted object was treated as a raw-count container: imported StereoPy reductions were provenance only. The `Spatial` counts were entry-filtered with the frozen same-batch contract, then independently SCT-normalized and used to compute fresh PCA, cosine Annoy neighbours, Leiden candidates and UMAP. The final broad resolution was selected after reviewing adjacent-resolution stability, full-gene marker/anti-marker interpretability, cluster size and spatial morphology; it was not copied from the companion Scanpy or BANKSY run.

## 2. Broad-first annotation and pool construction

High-confidence current-query anchors were built for follicular somatic, stromal/mesenchymal, vascular/endothelial, contractile, immune, epithelial/mesothelial and strict oocyte programs. Candidate Theca, Smooth muscle and Pericyte/mural were audited independently rather than hidden inside a combined follicular-wall or perivascular name.

Sheep immunoglobulin LOC aliases were resolved before closing an apparently immune-negative pool. A stable multi-locus immunoglobulin cluster with an independent B/plasma regulator can support broad Immune despite sparse `PTPRC`; single `JCHAIN`, single immunoglobulin, CD74 or MHC signals still fail closed. Plasma/antibody-secreting character is retained as a state tag unless the optional fine-label gate independently passes.

Uncertain observations were grouped by competing biological axes, not by graph-cluster number:

- follicular somatic review for Granulosa, steroidogenic Theca and follicular stroma;
- stromal/mesenchymal/mural review for generic ECM stroma, mature Smooth muscle and Pericyte/mural;
- vascular/endothelial/mural review for endothelial tracks and surrounding mural/stromal cells;
- immune and epithelial review only where the current query required them;
- a strict oocyte candidate separated from zona-positive adjacent Granulosa;
- a post-clustering QC holdout kept distinct from entry-QC exclusions.

Every pool had immutable membership, source cluster, generation, competing candidate lineages, state/spatial/QC tags, route history and supersession links.

## 3. Iterative routes

Large signal-bearing unresolved pools underwent balanced-anchor query-only reclustering. Local interfaces first underwent targeted anchor review; RCTD was lower-priority assistance and did not force medium/low calls. The complete low-information pool underwent a full QC-pool anchor recluster before remaining rejects were considered for Atlas rescue.

Atlas rescue was broad-only. The external Atlas, current-query marker/anti-marker channel and an internal-anchor or observed-density spatial channel had to agree under thresholds calibrated on disjoint query-like held-out current-query anchors. High and moderate-only returns entered the final broad membership with `fine_anchor_eligible=false`; low calls remained review/QC. A reference self-split was never treated as query calibration.

## 4. Rare and easily confused boundaries

- Oocyte used a two-tier route: the complete multi-module starting candidate pool entered query-only reclustering, while high-specificity marker/anti-program observations and compact spatial foci served only as strict seeds/support. Final Oocyte required a recluster cluster with multiple non-zona identity/maternal-ooplasm genes, somatic anti-program clearance and compatible follicular-object morphology. Zona transcripts in otherwise coherent Granulosa or stromal neighbours were treated as ambient/adjacent signal and rerouted to the appropriate somatic pool.
- Steroidogenic Theca required a coherent enzyme/regulator program and outer-follicular morphology. ECM or `LHCGR` alone remained follicular stroma/interface.
- Smooth muscle required a mature contractile core plus coherent tracks; `ACTA2/TAGLN` alone remained a contractile stromal state.
- Pericyte/mural required a multi-gene mural backbone and vascular adjacency. When no clean independent backbone survived full-feature/reclustering review, the correct output was a negative audit and a stromal/vascular state tag, not a forced class.
- A candidate stromal population with a very strong isolated neuronal gene remained Stromal/mesenchymal because the full neuronal, glial and neuroendocrine programs and nerve-like morphology failed. A single `DLG2`-like signal became a state tag, not a Neural/Schwann or Neuroendocrine label.

## 5. Subtype restraint

Graph subclusters were merged unless an independently supported functional program would be lost. ECM amount, cortical location, contractility without a mature lineage backbone, hypoxia, stress, cell cycle, low RNA and marker intensity were retained as tags. Granulosa stage/state names were used only when a coherent program, negative evidence, stability and morphology passed; zero fine subtypes remained valid.

## 6. Release and reporting

The frozen release published one annotation. The broad DEG and both broad canonical/data-specific tree dotplots used every accepted biological observation, including Atlas/anchor-rescued broad-only observations. Fine DEG/dotplots used only high-confidence observations with real fine labels; no synthetic subtype was created for broad-only rescue.

Entry-QC exclusions, interfaces, remaining QC holdout and technical states were fully accounted for but excluded from the biological broad census and DEG. The final Chinese HTML contained the annotated whole-section map, expandable broad/fine tree, per-node spatial highlights, broad and subtype tree dotplots, marker spatial maps grouped by supported cell type, route/threshold/outcome panels, a Chinese event timeline, raw state records and checksums. After every annotation route and final writeback completed, the main Agent judged biological quality against this process reference; final assets were generated only after that approval and user confirmation bound the frozen ledger and completion gate.

## 7. Forward-test acceptance

A future sheep-ovary Agent passes this reference test only if it can independently reproduce the safeguards above from raw inputs. Agreement with a historical label count is irrelevant. It fails if it widens Oocyte from zona signal, creates a neural class from one gene, hides steroidogenic/contractile/mural alternatives in a catch-all name, skips the full QC-pool recluster, calibrates on Atlas self-classification, excludes rescued cells from final broad DEG/dotplots, or substitutes many weak subtypes for a reliable shallow tree.
