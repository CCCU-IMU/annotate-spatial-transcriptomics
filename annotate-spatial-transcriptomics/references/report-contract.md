# Final report contract

## Single public annotation

Every released analysis observation has exactly one `final_cell_type`:

- a high-confidence, parent-locked fine identity replaces its broad parent;
- otherwise the frozen broad identity is used;
- terminal unresolved biological/QC membership is shown as `QC/Unknown`.

`final_broad_label`, `final_fine_label`, candidate IDs, confidence and assignment source remain internal audit provenance. They must not be presented as parallel public annotation systems. Public census, overview, spatial highlights, DEG, canonical/data-specific dotplots and the main downloadable annotation table all use `final_cell_type`.

## Required report sections

1. Input/context audit and selected clustering rationale.
2. One high-contrast whole-tissue `final_cell_type` spatial/UMAP overview.
3. One final-cell-type census, including `QC/Unknown` separately from biological identities.
4. One switchable normalized/absolute canonical marker dotplot grouped by `final_cell_type`.
5. One data-specific marker/DEG result grouped by `final_cell_type`.
6. A fixed-point-size high-contrast spatial highlight and support record for every final cell type.
7. All-cell canonical-marker spatial panels; marker maps never filter observations by assigned identity.
8. Independent state summaries that do not replace cell identity.
9. Atlas, missing-lineage, Oocyte and follicle-ROI biological review summaries.
10. Internal broad/fine/state provenance, cohort routes, checksums and complete workflow in collapsed audit sections.

If formal completion is absent or blocked, do not build a user-approval or final report. A controller result with `PENDING_USER_REVIEW_HIGH_UNRESOLVED` may build an explicitly diagnostic `pending_user_review` report so the unresolved regions and typed reasons can be inspected; that report has no approval or release authority and must not trigger another global residual/QC-anchor loop. For a PASS membership, use `pending_user_review` before user confirmation and `approved_final` after explicit confirmation without changing membership.

## Assets and rendering

Use `scripts/prepare_report_metadata.py` first. Its `primary_final_cell_type` is the only public label column. `primary_broad_label` and `primary_subtype_label` are internal audit columns only.

Use `scripts/build_annotation_maps.R --final-cell-type-col final_cell_type` to generate one final-cell-type UMAP/spatial overview and one per-type highlight. Spatial maps use fixed point diameter, black background and saturated colours; rare groups are never enlarged automatically. Endothelial and Pericyte/mural must use visibly distinct colours. Single-cell projects omit spatial sections rather than fabricating coordinates.

Use `scripts/run_final_label_deg.R --final-cell-type-col final_cell_type` for the public one-vs-rest DEG. Optional broad/fine DEG may be produced only under an explicitly marked internal-audit flag and must not feed the public report.

Use `scripts/build_marker_dotplots.R --final-cell-type-col final_cell_type`. Every non-QC final cell type must appear as both a dotplot label and marker group. Source TSVs contain `gene`, `label`, `avg_expression`, `pct_expressed_absolute`, `n_observations`, `marker_group`, `avg_expression_scaled_within_gene`, `pct_expressed_scaled_within_gene`, `analysis_view` and `evidence_cohort`.

Point size uses within-gene normalized detection from 0 to 100. Colour uses within-gene scaled average expression with a documented clip. The same card switches to absolute detection/mean expression derived from the identical source table. Render PNG and PDF pairs.

Use `scripts/build_spatial_gene_maps.R` for all-cell marker projections. Pass the expected observation count so partial intersections fail. Missing sheep gene symbols remain explicit as unavailable and are not interpreted as biological absence.

Use `scripts/write_frozen_annotations_to_seurat.R`; analysis membership plus optional excluded-initial-QC membership must exactly cover the object. The writer stores `spanno_v2_2_cell_type` and retains broad/fine provenance separately.

Use `scripts/build_frozen_review_report.py` for the sample-agnostic HTML. Every final-cell-type highlight shows route/definition, count, canonical markers, top DEG, competing/anti-program review and spatial support. Internal fine-candidate audits may be available in a collapsed provenance section but are not a second public taxonomy.

## Release audit

Run `scripts/build_release_session_info.py`, `scripts/build_release_manifest.py` and `scripts/audit_release.py --profile full`. The release audit must prove:

- exactly one `final_cell_type` per analysis observation;
- no releasable `Vascular-associated` label;
- no second strict/inclusive/display taxonomy;
- complete final-cell-type census, map, DEG and marker dotplot;
- valid report links and PNG/PDF/source assets;
- content hashes for membership, report, RDS, tables and figures.

A cohort identifier is provenance, never a public cell type.
