# Standard sheep-ovary R-first annotation workflow

Use this contract whenever the resolved context is sheep/Ovis/ovine plus ovary/ovarian and a full-feature Seurat RDS is available. It converts the successful R/Scanpy strategy into a low-freedom execution path. Parameters outside this contract require a recorded user-approved exception.

## Fixed preprocessing and resolution contract

Treat a StereoPy-converted `cellbin_PPed` RDS as a raw-count container. Rebuild `Spatial counts -> SCT v2/glmGamPoi -> PCA50 -> 30-PC SNN(k=30, Annoy/cosine) -> UMAP`, and never reuse imported StereoPy PCA/UMAP.

Use `0.1,0.2,0.3,0.4,0.6` as the complete formal candidate grid for whole-tissue clustering and every broad-class or targeted reclustering cohort. Select adaptively within that grid. Do not recluster the final QC holdout. Do not run, rank, interpret or write back `0.01`, `0.02`, `0.05` or another below-grid resolution. If cluster inflation occurs, audit the graph, input reduction, neighbor construction, connected components and imported metadata; do not lower the resolution to conceal a graph defect.

Run `scripts/validate_resolution_grid.py` before every whole-tissue, broad-class or targeted-cohort submission. Every completed cohort records the exact candidate grid and selected resolution. A missing grid or an out-of-contract value invalidates the result.

## Ordered biological phases

1. **Whole-tissue broad pass.** Generate complete DEG, marker/anti-marker, UMAP, whole-section spatial and per-cluster highlights for the normal grid. Freeze the lowest-complexity resolution that preserves supported broad lineages.
2. **Clusterwise open-world initial broad annotation.** At the selected resolution, run `init_open_world_lineage_audit.py` and review every catalog boundary plus every unexplained multi-gene program. Assign a moderate-or-higher initial broad label directly to each supported cluster. Low-information, featureless or irreducibly mixed clusters enter `qc_holdout`. Do not create biological review pools at this stage.
3. **One reclustering cohort per initial broad class.** Submit an adaptive query-only reclustering job for every initial broad class. A class too small for defensible reclustering closes as `not_applicable_reviewed`, remains broad-only and records the limitation. These are immutable computational cohorts, not biological pools.
4. **Subcluster adjudication.** A coherent high-confidence subtype receives broad plus fine labels. A subcluster that supports only its parent returns directly as broad-only. State/intensity differences remain tags. Do not manufacture subtypes to match a reference.
5. **Direct cross-lineage return.** If a subcluster from one broad-class run coherently supports another lineage—for example a pregranulosa-like subcluster found during Oocyte-candidate reclustering—write it directly to the supported target broad/fine label. Preserve source run, source cluster, cohort, membership hash and evidence. Do not create an intermediate target pool and do not automatically recluster it again with the target class.
6. **Decision-relevant mixtures only.** When sufficient observations contain interpretable competing biological programs, create a temporary `targeted_recluster_cohort`; Oocyte contamination-safe review is one such cohort. If targeted reclustering remains insufficient, calibrated RCTD is lower-priority assistance: extreme plus independent evidence may support a fine call, high may support broad-only, and medium/low returns to QC holdout. A featureless or uninterpretable mixture goes directly to QC holdout.
7. **Residual QC Atlas, without QC reclustering.** After every broad-class and targeted decision is terminal, freeze the complete residual `qc_holdout` membership and map it directly using finalized internal broad anchors plus the declared external reference. Without a matched count-level reference, GSE233801 is primary. Calibrated moderate-or-higher consensus returns broad-only with `fine_anchor_eligible=false`; lower confidence remains `qc_holdout/qc_reject`. Record a zero-query terminal audit when no residual QC exists.
8. **Final single annotation.** Every analysis-set observation receives either one moderate-or-higher `final_broad_label` with an optional high-confidence `final_fine_label`, or an explicit retained QC/technical state. Preserve assignment mode, source lineage, confidence, state/spatial/QC tags and evidence. Do not create strict/inclusive/display copies.
9. **Completion, main-Agent quality approval and lightweight review.** Require exact analysis-set ownership, terminal broad-class/targeted cohorts, validated cross-lineage returns, a terminal residual-QC Atlas phase, current open-world audit, zero open incidents, current profile/preset binding, normal-grid compliance, final taxonomy audit and a fresh direct-workflow completion gate. Generate only the frozen broad spatial and canonical-marker evidence assets. The main Agent then reviews biological quality. Only after a pass build the lightweight user-confirmation HTML; final DEG, full tree dotplots, per-node/per-gene maps and release HTML wait for user confirmation.

## Final evidence and report cohort

All observations formally returned to a biological broad class—including direct cross-lineage returns and calibrated Atlas broad-only rescue—participate in that class's final DEG, canonical marker dotplot, data-specific marker dotplot, UMAP and spatial maps. Only observations with a real high-confidence fine label participate in fine-label DEG/dotplots. Retained QC, technical and initial-QC states remain outside biological DEG.

Publish one annotation tree and one final census. Follow the successful Scanpy report layout: overall annotated spatial map above an expandable tree; per-node jump buttons to spatial highlights; broad and fine marker panels; marker-gene spatial maps grouped by supported cell type; route/threshold/outcome tables; a top link to the workflow; and the detailed Chinese workflow/state reconstruction at the bottom.

## Failure and resource contract

Before submission, parse every generated R, Python, shell and AIP file with `scripts/preflight_generated_job.py`. Record every EXIT, cancellation, fail-closed gate, syntax error, dependency problem and resource mismatch before accepting a repair. A sample cannot complete with an open incident.

Reuse validated immutable SCT/PCA/SNN and per-resolution memberships by hash. Do not reread or recompute large expression objects for a postprocessing-only repair. Match requested CPUs to real workers; one Leiden optimization and unwrapped `FindAllMarkers` are single-threaded. Delay release figures and HTML until the user confirms the frozen first annotation.
