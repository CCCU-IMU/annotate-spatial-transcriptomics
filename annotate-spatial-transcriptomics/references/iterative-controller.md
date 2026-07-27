# Iterative controller

Iteration occurs inside explicit v2.2 phases. It does not mean repeatedly tuning a first-pass all-object cellbin classifier.

## Round 0: freeze inputs

Bind full object, analysis set, excluded initial QC, project-local raw counts, profiles, catalog, seed and artifact roles. Historical labels and failed diagnostics remain outside runtime evidence.

## Round 1: build cohorts

Select a stable whole-tissue partition and create exactly one cohort for every initial cluster. Record only provisional biology, mixed/unknown status and lineage watches. No formal label or QC membership exists in this round.

## Round 2: annotate second-round subclusters

Recompute SCT/PCA/SNN/Leiden from raw counts independently inside every cohort. Select a resolution from the full grid and nearest neighboring resolutions. Score each subcluster against the full catalog before revealing its provisional parent. Return high-purity subclusters wholesale; preserve a supported broad when fine evidence is weak; mark genuinely mixed subclusters for local resolution. If high resolution creates many marker-poor unresolved fragments, review the nearest lower stable resolution.

## Round 3: resolve local mixtures

Inside a triggered mixed subcluster, form independent lineage components, adjudicate overlaps and rescore the exact local remainder once. Unselected or ambiguous members are not technical QC. A missing broad can be reconstructed here even if it never appeared in Round 1.

## Round 4: merge and freeze

Merge all cohort outcomes into an exact, mutually exclusive analysis-set partition. Freeze broad labels, then parent-lock fine proposals and retain state programs in separate fields.

## Round 5: Atlas and completeness

Map all cells once at broad level. Directly rescue only unlabeled, moderate-or-higher, non-OOD and profile-compatible observations. Review material conflicts with the complete source cohort. Audit present, missing and unmodeled lineages, then convert the remaining unresolved cells to typed QC during final materialization.

An excessive final QC fraction routes to the specific contributing second-round cohort or post-merge unresolved group. It never routes to whole-object per-cellbin threshold tuning.
