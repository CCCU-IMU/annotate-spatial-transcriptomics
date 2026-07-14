# Evidence and unresolved-population routing

## Minimum evidence

Use coherent marker programs rather than single genes. Retain absolute detection percentages, average expression, DEG effect size, contradictory programs, QC metrics and spatial location.

When matched single-cell data are available, use `matched-single-cell-reference.md`. Query-derived full-feature anchors remain the primary biological evidence; a stage-matched same-project single-cell reference is the preferred external channel, ahead of public same-species and cross-species atlases. Its source labels must pass an explicit crosswalk and transfer ceiling.

## Decision routes

### Direct definition

Use when multiple canonical markers, DEG and morphology agree. Assign broad label first; subtype only if subtype-specific evidence is coherent.

### Broad review pool plus anchors

Use for large unresolved populations with usable transcriptomic structure. Add frozen, R-only or platform-matched trustworthy anchors. Recluster query and anchors, then classify only clusters supported by multiple channels. Broad returns are not automatically subtype anchors.

Anchors participate in joint normalization/PCA only. Query cells alone define graph, clustering, UMAP, DEG and outcome counts. Require `query_or_anchor` and `anchor_label` in the frozen membership and emit source/state/QC composition tables. Otherwise record the run as query-only, not anchor-assisted.

### Targeted RCTD/reference assistance

Use for a local mixed interface, not as a whole-dataset label generator. Build reference types independently, match query depth for validation, compare parameter settings, calibrate thresholds on held-out anchors and preserve rejects. Record weights, margins, precision and coverage.

Recheck the query/reference boundary for every child pool, even when reusing a previously valid reference. Descendant query pools can contain observations that were anchors for an earlier branch. Run `scripts/filter_reference_query_overlap.py` before fitting; require zero final overlap, adequate counts for every retained reference label and a machine-readable exclusion manifest. A boundary failure is a preserved failed run and must be repaired under a new run/route ID.

`scripts/run_rctd_review.R` supplies a versioned local-review adapter. Require explicit query-depth `UMI_min`/`counts_MIN`; its results and weights are evidence artifacts only.

Use a tiered fallback: extreme confidence plus independent current-query evidence may support a subtype; high confidence may support only a broad return; medium/low confidence proceeds to depth-matched atlas/internal-anchor mapping. If atlas also rejects it, retain it explicitly. Never equate an RCTD reject with a terminal anatomical interface without checking size and spatial locality.

### QC holdout plus atlas

First determine whether low information is technical, biological or depth-related. Except for verified zero-count observations, a large post-clustering QC holdout must first undergo complete broad-anchor query-only reclustering. Only the remaining rejects proceed to depth-matched atlas/internal-anchor review. Validate reference performance after downsampling anchors to query depth. Combine atlas, internal anchors, markers and observed-density space. A spatial or atlas top label alone cannot write back. Rescue broad classes at predeclared moderate/high confidence; retain observations below the moderate tier as calibrated rejects with the mapping artifact and rationale.

The default atlas policy is route- and class-specific held-out target precision `moderate-or-higher=0.90`, `high=0.95`. These are calibration precision targets used to derive score/margin thresholds, not requirements that every observation's raw confidence be at least 0.90 or 0.95. Calibrate nested cumulative selections on the full held-out group and require the high selection to be a subset of moderate-or-higher. Output counts are `high`, `moderate-only` and `low-reject`; also report the cumulative `moderate-or-higher` count so a nonzero high tier can never be misdescribed as zero observations meeting the moderate gate. Both accepted tiers may return **broad-only** labels after independent current-query evidence passes, and both must set `fine_anchor_eligible=false`.

Prefer query-like held-out anchors for final transfer calibration. The completed Scanpy-style route used external Atlas, internal anchors, marker programs and observed-density space as four concordance channels; the completed R-style route calibrated those channels on held-out current-query anchors before high/moderate consensus. An external-reference held-out split measures Atlas self-consistency and is useful diagnostically, but it cannot by itself certify transfer confidence in low-depth query observations.

A biologically matched single-cell reference may replace the generic external-atlas channel, but not the internal-anchor, current-query marker/anti-marker or spatial channels. Use its count-level object for depth matching and held-out calibration. If only a marker dotplot is available, use it to refine positive/negative programs and the crosswalk; do not report cell-level mapping confidence from the dotplot.

Use `scripts/adjudicate_multichannel_broad_rescue.py` after assembling the evidence matrix. Configure each channel with its own label and acceptance column, require a current-query marker/anti-marker channel, and require at least one route/internal-anchor/observed-density spatial channel. The script derives a candidate label from concordant channels, calibrates cumulative support-count thresholds separately for each declared route/label group, and refuses to extrapolate to support counts absent from held-out anchors. Its output is proposal-only and never edits the cell ledger.

The default structural minima are two concordant channels for moderate-or-higher and three for high, but held-out precision calibration may select stricter observed thresholds. High must remain a subset of moderate-or-higher. A group without enough query-like held-out support remains rejected; do not borrow a permissive threshold from an unrelated route merely to increase coverage.

For H5AD references, first use `scripts/downsample_anchors_to_query_depth.py` to create disjoint held-out **current-query anchors** matching the QC pool's count-depth distribution. `scripts/run_anchor_knn_mapping.py` creates auditable atlas/internal-anchor predictions for those anchors and the QC query; calibrate nested high/moderate tiers with `scripts/calibrate_tiered_mapping_thresholds.py` and its required origin manifest. `scripts/build_depth_matched_atlas_bundle.py` splits an external reference and is diagnostic reference-self-classification only; it cannot calibrate query rescue. The legacy `calibrate_mapping_thresholds.py` is diagnostic-only and never writeback eligible. Accepted rows remain `defined_broad_only`, `fine_anchor_eligible=false`; rejects remain QC/review. More sophisticated mapping methods may replace kNN but must obey the same disjoint query-like calibration contract.

### Retain uncertainty

Use `interface_review`, `qc_holdout` or `technical_state` when evidence remains insufficient. Annotation rate is not an optimization target.

Before retention, compare the proposed interface/QC fraction with the current analysis-set structure. This is a diagnostic trigger, not a quota: a large or diffuse retained pool must demonstrate that broad anchor reclustering, RCTD tiering and atlas/internal-anchor fallback were all exhausted. Do not preserve a high unresolved fraction simply to maximize specificity.

## Contamination and ambient signal

Require lineage-core programs and spatially plausible focal structure. Neighboring or ambient markers must not override a coherent resident program. Report observation/bin counts separately from biological cell counts.
