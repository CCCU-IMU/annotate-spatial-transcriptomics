# annotate-spatial-transcriptomics v2.0.4

This release makes observation-level subset writeback the fail-closed path for mixed graph subclusters.

- Whole-subcluster writeback now requires an absolute lineage-support floor and purity margin.
- A supported subset is adjudicated by its own complete-catalog evidence and may differ from the aggregate subcluster winner.
- Multiple disjoint lineage subsets may be returned from one source subcluster; the exact remainder is rescored independently.
- Final broad-class completeness audits every source writeback instead of relying only on a final-label average.
- Every present parent broad class audits its complete machine-actionable fine-candidate catalog; zero fine labels remains valid only after explicit review.
- Smooth muscle still requires a mature MYH11/CNN1/ACTG2-centered program with mural exclusion; ACTA2/TAGLN or location alone cannot drive return.

The project framework schema remains `2.0.0`.
