# Iterative controller

Iteration occurs inside explicit v2.2 phases. It does not mean repeatedly tuning a first-pass all-object cellbin classifier.

## Round 0: freeze inputs

Bind full object, analysis set, excluded initial QC, project-local raw counts, profiles, catalog, seed and artifact roles. Historical labels and failed diagnostics remain outside runtime evidence.

## Round 1: build cohorts

Select a stable whole-tissue partition and create exactly one cohort for every initial cluster. Record only provisional biology, mixed/unknown status and lineage watches. No formal label or QC membership exists in this round.

## Round 2: annotate second-round subclusters

Recompute SCT/PCA/SNN/Leiden from raw counts independently inside every cohort. Select a resolution from the full grid and nearest neighboring resolutions. Score each subcluster against the full catalog before revealing its provisional parent. Return a subcluster with one clear independent identity wholesale; shared background candidates remain watches. Mark a subcluster mixed only when two material, independently supported direct identity components are separable, or when a specific component is embedded in a supported generic remainder. If high resolution creates many marker-poor unresolved fragments, review the nearest lower stable resolution; do not force a subtype when only broad identity is supported.

## Round 3: resolve local mixtures

Inside a triggered mixed subcluster, first reuse the selected resolution's neighboring Leiden partitions. Only unresolved mixtures proceed to candidate-local spatial observation components. Adjudicate overlaps and rescore the exact local remainder once. Unselected or ambiguous members are not technical QC. A missing broad can be reconstructed here even if it never appeared in Round 1.

## Round 4: merge and freeze

Merge all cohort outcomes into an exact, mutually exclusive analysis-set partition. Freeze broad labels, then parent-lock fine proposals and retain state programs in separate fields.

## Round 5: Atlas and tissue-specific completeness

Map all cells once at broad level. In sheep ovary, the active reference is the immutable contract-bound GSE233801 split-wall v2 bundle; only capability-matrix-supported and current-query-class-calibrated source/release pairs may rescue unlabeled, moderate-or-higher, non-OOD and profile-compatible observations. Challenge-only classes can trigger review but cannot write labels. Review material conflicts with the complete source cohort. Run the tissue-specific Oocyte and follicle-ROI endpoints before broad closure.

## Round 6: strictly serial per-broad whole-query review

Expose exactly one active cell type. Announce `现在开始对 <cell type> 进行专项复核。`, then independently resolve its current-member precision, whole-query recall, molecular identity, whole-section spatial consistency and literature boundary. Apply only an exact cell-ID patch or a bounded source-cohort review when evidence requires it. Close that type before opening the next. Reuse shared raw-count/coordinate caches, but never use one batch calculation to close several biological types.

After all types close, audit present, missing and unmodeled lineages, then convert the remaining unresolved cells to typed QC during final materialization.

An excessive final QC fraction routes to the specific contributing second-round cohort or post-merge unresolved group. It never routes to whole-object per-cellbin threshold tuning. Resume a controlled per-type pause from the frozen transform chain; do not repeat Atlas, rescue, ROI repair or stable clustering.
