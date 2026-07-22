# Final report contract

## Required sections

1. Data/input audit and selected clustering rationale.
2. Final broad spatial/UMAP overview.
3. Expandable broad-to-subtype annotation tree with node highlights.
4. Broad marker evidence.
5. High-confidence subtype marker evidence when at least one fine label passes; otherwise an explicit zero-subtype audit.
6. Broad DEG downloads and, when released, subtype DEG downloads.
7. Spatial marker maps grouped by supported cell type.
8. Unresolved, interface and QC policies with counts.
9. Complete chronological workflow and state provenance.
10. Software/session information, audit results and checksums.
11. Biological context/profile/preset, open-world lineage discovery, iteration queue, every cohort/direct-return/run status, triggered Oocyte/context-specific object audits and completion-gate result.
12. Full-feature audit plus label-level full-feature marker validation; clearly distinguish clustering/HVG evidence from final validation evidence and expose the manifest of the LogNormalize validation-only object when the Seurat raw container had `Spatial@data == counts`.
13. One final census and map. Moderate-or-higher broad assignments supply broad DEG/markers; high-confidence real fine labels supply subtype DEG/markers.
14. Direct-workflow dashboard: each initial broad class, broad/targeted cohort, direct/cross-lineage return, optional RCTD, prelabel evidence freeze, all-cell Atlas concordance, frozen-QC writeback/rejects, OOD census, calibrated thresholds and every closed discrepancy review.
15. A detailed Chinese workflow at the bottom reconstructed from `workflow_event_registry.tsv`, starting from raw-object loading and including every cohort generation/direct return, scheduler failure/repair, biological review and atomic writeback. Link to it from the top navigation.
16. Source ancestry/control board: `source_key → initial broad decision → cohort/run/resolution/subcluster → direct or assisted action → final label/confidence`.
17. Two separately labeled top-level censuses: biological broad classes and retained anatomical-interface/QC/technical/pending states. Only the former enters the broad biological tree, DEG and marker dotplots.
18. Continuous lineage-signal coverage dashboard: required boundaries, selected-plus-two-higher resolutions, candidate coverage, open/closed weak signals, negative closures and any cross-lineage class reconstructed after cohort reclustering.
19. In the pre-confirmation report, canonical marker-expression spatial panels for every broad class and released fine label across all analysis-set observations. These maps do not filter, fade or recolor by assigned cell type; they use fixed point size, a common background and documented per-gene expression scaling.

The annotation tree must provide expand/collapse/search controls and a direct link from every broad/subtype node to its spatial highlight. Put the reviewed annotated spatial overview immediately above the tree. Keep the detailed Chinese raw-input-to-release workflow at the bottom, with a top navigation entry to it.

If the completion gate is absent or blocked, no master quality approval, confirmation or final report may be built. After it passes, generate one high-contrast broad spatial PNG, one canonical broad marker dotplot PNG and compact all-cell canonical-marker spatial panels, freeze the fully annotated snapshot, and obtain the main conversation Agent's concise biological quality approval. Only then build `review/confirmation/index.html` from the approved result and `state/annotation_support_registry.tsv`. This lightweight, self-contained HTML is the only pre-user-approval report and must not run final DEG or full release assets. Then freeze `provenance/final_annotation_confirmation_request.json`; its hashes bind the current cell/cluster/support ledgers, completion/taxonomy/master-quality records and lightweight review. A biologically complete release report is built only after explicit user approval and every cohort/assisted route is terminal with rationale.

## Mandatory dotplots

Always produce broad-class tree dotplots. Produce subtype tree dotplots only when at least one high-confidence fine label is released; zero subtypes is a valid outcome and must not trigger synthetic labels. For every released level produce canonical and data-specific versions when data-specific markers are available.

After all report assets are stable, run `scripts/build_release_session_info.py`, then `scripts/build_release_manifest.py`, followed by `scripts/audit_release.py --profile full`. The checksum list must cover the report, review package, ledgers, DEG tables, dotplot sources and rendered figures. The manifest, checksum file and release-audit output are intentionally excluded from their own hash set to avoid a circular release dependency.

Completion, confirmation/review assets, final HTML and release assets write a sibling `.deps.json` content-hash manifest. The controller may inspect mtime as a cheap scheduling hint for lightweight state, but it must decide whether an expensive asset is stale from the recorded target/dependency hashes. A touched-but-identical dependency does not trigger recomputation; a content change does, even when mtime is preserved.

Each source TSV must include `gene`, `label`, `avg_expression`, `pct_expressed_absolute`, `n_observations`, `marker_group`, `avg_expression_scaled_within_gene`, `pct_expressed_scaled_within_gene`, `analysis_view` and `evidence_cohort`.

Point size uses within-gene normalized detection from 0 to 100. Color uses within-gene scaled average expression with a documented clip. The same report card must switch between this view and an absolute detection/average-expression view built from the same source table; the release audit requires both PNG/PDF pairs.

Render the label dendrogram on the left. Put marker genes on the x axis and facet/group them by the current cell type or program they support. Every final broad label and every high-confidence final fine label must occur as both a dotplot label and a marker group in canonical and data-specific panels.

The broad-class DEG and both dotplots use the final non-QC/non-interface biological cohort: every cellbin formally returned to a broad class participates. Subtype DEG/dotplots use only high-confidence cells with a real fine label; broad-only rescue is excluded and never receives a synthetic subtype. Every source table declares `analysis_view=final` and its evidence cohort.

A cohort identifier is provenance, not a cell type. Show it in ancestry/route tables, never as a biological tree node unless an independently gated biological label happens to use different, approved wording.

## Asset formats

All scientific figures require PNG and PDF. Every report link must resolve. Per-node and per-gene assets require an index TSV. The report must state whether observations are cells, nuclei, spots or cellbins.

Use `scripts/build_annotation_maps.R` to generate broad/subtype UMAPs and, when coordinates exist, broad/subtype spatial maps plus per-node highlights. Single-cell projects omit spatial sections explicitly rather than fabricating coordinates.

First run `scripts/prepare_report_metadata.py`. Pass its `primary_broad_label` to every broad DEG/dotplot/map command and `primary_subtype_label` to every subtype command. Use `retained_state_display` only for the separate QC/interface/technical census. Never use `broad_display` or a retained-state fallback as subtype input, and never manufacture a `Broad only: ...` subtype for broad-only rescue cells.

Use `scripts/build_spatial_gene_maps.R` for spatial projects. It generates one PNG/PDF pair per marker and compact per-marker-group panels, all using every analysis-set observation rather than label-filtered cells. Its indexes record scope, denominator, requested and available markers. Missing markers remain explicit rather than silently disappearing.

Use `scripts/run_final_label_deg.R` after label writeback to generate separate broad and subtype one-vs-rest DEG tables. For SCE it records aggregate normalized mean/detection and fold change; for Seurat it runs sampled Wilcoxon markers. State the method and its limitations in the report.
