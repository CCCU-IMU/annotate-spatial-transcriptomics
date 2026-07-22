# annotate-spatial-transcriptomics v2.0.2

This patch closes a biological writeback loophole without replacing the v2 workflow.

- Initial broad clustering now establishes supported resident compartments instead of forcing every final lineage to appear.
- Initial and cohort decisions use one machine-derived lineage table; validators recompute winner, runner-up and margin.
- Whole-subcluster biological returns require two positive families, contradiction clearance and observation-level purity. Parent cohort identity is not evidence.
- Epithelial/mesothelial recall is limited to high-purity reclustered subclusters rather than expanded mixed clusters or spatial components.
- Mature Smooth muscle requires a MYH11/CNN1/ACTG2-centered program after mural exclusion; mural vascular-wall cells remain Vascular-associated.
- One public completion interface retains existing internal validators under four phases.
