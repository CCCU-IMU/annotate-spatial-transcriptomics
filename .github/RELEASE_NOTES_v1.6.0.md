# Annotate Spatial Transcriptomics v1.6.0

This release makes the direct-lineage cohort controller the single architecture for all new projects.

- Initial broad classes are assigned directly after open-world review.
- Every supported broad class receives one immutable reclustering cohort; small classes may close broad-only with a reviewed skip.
- Subclusters return directly to parent/fine/cross-lineage labels, a one-question targeted cohort, or QC.
- Persistent biological pools are retired for new projects.
- Terminal residual QC goes directly to calibrated Atlas review without mandatory QC reclustering.
- RCTD extreme plus independent evidence may support fine, high supports broad-only, and medium/low enters QC.
- State, planner, completion, master-review and report contracts now bind cohort/direct-return provenance.
- Clean GitHub runners now install declared regression-test dependencies before PR validation and Release packaging; PR branches run one validation check instead of duplicate push and pull-request checks.

Historical pool registries remain readable for migration but are not active routing instructions.
