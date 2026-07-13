# Final report contract

## Required sections

1. Data/input audit and selected clustering rationale.
2. Final broad spatial/UMAP overview.
3. Expandable broad-to-subtype annotation tree with node highlights.
4. Broad marker evidence.
5. Subtype marker evidence.
6. DEG downloads at both levels.
7. Spatial marker maps grouped by supported cell type.
8. Unresolved, interface and QC policies with counts.
9. Complete chronological workflow and state provenance.
10. Software/session information, audit results and checksums.
11. Biological context/profile, iteration queue, every pool/run status, rare-cell object audits and completion-gate result.
12. Full-feature audit plus label-level full-feature marker validation; clearly distinguish clustering/HVG evidence from final validation evidence.
13. Strict, inclusive and display census/maps, with an explicit statement of which cohort supplies DEG and marker discovery.
14. Multi-route dashboard: each Route A–E input, applicability, anchors/reference, selected parameters, calibrated thresholds, outcome counts, retained rejects and no-repeat status.
15. A detailed Chinese workflow at the bottom reconstructed from `workflow_event_registry.tsv`, starting from raw-object loading and including every pool generation, scheduler failure/repair, biological review and atomic writeback. Link to it from the top navigation.
16. Source ancestry/control board: `source_key → parent decision → pool snapshot → run/resolution/subcluster → route/action → strict/inclusive/display`.

The annotation tree must provide expand/collapse/search controls and a direct link from every broad/subtype node to its spatial highlight. Put the reviewed annotated spatial overview immediately above the tree. Keep the detailed Chinese raw-input-to-release workflow at the bottom, with a top navigation entry to it.

If the completion gate is absent or blocked, no final report may be built. After it passes, freeze and present `provenance/final_annotation_confirmation_request.json`; do not spend compute on final release assets until the user explicitly approves it. `state/final_annotation_confirmation.json` must hash the current cell ledger, cluster ledger and completion gate. A biologically complete report is built only after every pool is closed with rationale and this confirmation remains current; the release audit separately verifies the figure/table contract.

## Mandatory dotplots

Produce separate broad-class and subtype tree dotplots. For each level produce canonical and data-specific versions when data-specific markers are available.

After all report assets and session files are stable, run `scripts/build_release_manifest.py`, followed by `scripts/audit_release.py --profile full`. The checksum list must cover the report, ledgers, DEG tables, dotplot sources and rendered figures. The manifest, checksum file and release-audit output are intentionally excluded from their own hash set to avoid a circular release dependency.

Each source TSV must include `gene`, `label`, `avg_expression`, `pct_expressed_absolute`, `n_observations`, `marker_group`, `avg_expression_scaled_within_gene`, `pct_expressed_scaled_within_gene`, `analysis_view` and `evidence_cohort`.

Point size uses within-gene normalized detection from 0 to 100. Color uses within-gene scaled average expression with a documented clip. Absolute detection and average expression remain in the source table and report legend.

Render the label dendrogram on the left. Put marker genes on the x axis and facet/group them by the current cell type or program they support. Every strict DEG label, including validated rare labels above the project display minimum, must occur as both a dotplot label and a marker group in canonical and data-specific panels.

Canonical and data-specific marker discovery defaults to the strict, non-QC, non-interface, fine-anchor-eligible cohort. Display/inclusive plots may be supplemental but cannot replace strict evidence plots. DEG must distinguish strict biological evidence from descriptive display-level aggregates.

## Asset formats

All scientific figures require PNG and PDF. Every report link must resolve. Per-node and per-gene assets require an index TSV. The report must state whether observations are cells, nuclei, spots or cellbins.

Use `scripts/build_annotation_maps.R` to generate broad/subtype UMAPs and, when coordinates exist, broad/subtype spatial maps plus per-node highlights. Single-cell projects omit spatial sections explicitly rather than fabricating coordinates.

Use `scripts/build_spatial_gene_maps.R` for spatial projects. It generates one PNG/PDF pair per marker and an index grouped by the biological cell type or program supported by that marker.

Use `scripts/run_final_label_deg.R` after label writeback to generate separate broad and subtype one-vs-rest DEG tables. For SCE it records aggregate normalized mean/detection and fold change; for Seurat it runs sampled Wilcoxon markers. State the method and its limitations in the report.
