# Staged second-round annotation controller

`run_lineage_controller.py` is the only formal entry and the only chain that may ultimately create release membership. The controller binds its scripts, candidate catalog, profiles, inputs, seed and artifact roles in the annotation contract.

## Authority by phase

1. `whole_tissue_partition` creates only a stable initial-cluster partition and provisional cohort plan.
2. `cluster_cohort_recluster` independently reclusters every initial cluster from project-local non-SCT raw counts and freezes full-catalog subcluster evidence.
3. `local_mixed_subcluster_split` resolves only second-round subclusters with two or more separable competing identity cores.
4. `merge_and_freeze_broad` merges mutually exclusive second-round outcomes and creates formal broad membership.
5. `atlas_and_completeness_review` maps all cells at broad level, rescues only unlabeled cells and audits present, absent and unmodeled lineages.
6. `materialize_final_release` parent-locks fine labels, preserves states and converts terminal unresolved observations to typed QC.

The first phase has no release authority. It cannot emit formal broad, fine, QC or release membership. Artifacts registered as `failed_diagnostic` cannot enter any phase as membership or reference truth.

`cluster_cohort_recluster` may reuse a controller-generated `derived_partition_cache` only when the frozen input, exact cohort membership, resolution grid/contract, seed and clustering-script fingerprints are identical. Cache reuse skips raw-count SCT/PCA/SNN/Leiden fitting but reruns full-catalog scoring, resolution selection and adjudication. The cache carries no label authority and cannot originate under `failed_diagnostic`.

## Whole-tissue partition

Choose a stable resolution that preserves major anatomy without obvious technical fragmentation. The purpose is one initial cluster to one second-round cohort, not complete lineage separation. Store `provisional_broad`, `mixed` or `unknown`, plus candidate/competitor/watch programs. Observation scores may describe those programs but cannot assign formal cells.

All `analysis_set` members must occur exactly once in `whole_tissue_cohort_plan.tsv`. Provisional fields are withheld from the second-round scorer and resolution selector.

## Second-round cohort annotation

Every initial cluster, including an apparently pure one, runs:

`project-local raw counts -> SCT v2/glmGamPoi -> PCA -> query-only SNN -> Leiden grid`

Never run SCTransform on SCT corrected counts. Never substitute whole-tissue BANKSY graph parameters for cohort PCs/k. Evaluate the complete grid, select by identity-program separation, anatomy, neighboring-resolution stability and absence of technical/state fragmentation, and scan every selected subcluster against the complete open catalog plus unexplained programs.

The provisional parent must not narrow the candidates. **Direct cross-lineage return** and missing-broad reconstruction remain available to every second-round subcluster until full-catalog scores are frozen.

Only after scores freeze may the controller read the provisional record and name an outcome `parent_return`, `cross_lineage_return`, `missing_broad_reconstruction`, `fine_candidate`, `state_annotation`, `unmodeled_lineage_candidate` or `unresolved_biological`. An underpowered cohort is explicitly `underpowered_not_evaluable`; it never inherits the first-pass name silently.

A coherent high-purity subcluster may return wholesale. Sparse noncontradictory members inherit its broad identity. Dropout alone is not QC.

Resolve broad identity before fine identity. If no independent child discriminator is evaluable, return the supported parent or cross-lineage broad and leave fine empty. When many small subclusters are marker-poor or unresolved, review the nearest lower resolution before accepting biological uncertainty.

## Local mixed-subcluster resolution

Observation-level splitting is a conditional local operation, not an all-object classifier. Two independently supported identity cores, or one specific core embedded in a generic parent, trigger a bounded separability check. Do not require positive whole-subcluster DEG: a real minority identity can be diluted by the parent. Write a subset only after observation discriminators, candidate-local pseudobulk/DEG, space or neighboring partitions demonstrate separable members; otherwise retain the supported parent and keep the alternative as watch.

Each candidate independently proposes a spatial/expression subset. Resolve overlap with normalized evidence, pairwise discriminators and anti-programs. Specific lineages precede generic Stromal remainder. Ambiguous overlaps remain with a supported common parent or `unresolved_biological`; candidate ordering cannot assign them.

Group thresholds validate a proposed component and never admit individual observations. Exact remainder closure is scoped to the source subcluster, preserves original scores, permits one additional extraction and never converts “not selected” into QC.

## Broad freeze and final review

After every initial-cluster cohort closes, merge whole-subcluster returns, supported local subsets and local parent remainders. Require an exact disjoint cover of `analysis_set`; unresolved conflicts stay biological. Only this phase freezes formal broad membership.

Fine proposals become labels only under the matching frozen broad parent. State annotations remain separate from identity. The single all-cell Atlas pass occurs after broad freeze: it may rescue unlabeled moderate-or-higher, non-OOD, profile-compatible cells, but cannot silently overwrite a defined broad or create fine labels. Material conflicts reopen the complete source cohort once for query-derived evidence review.

Present/zero-census lineages, embedded programs, unmodeled programs and the complete parent-by-fine catalog are audited after merge. A post-merge component that passed its own strict multi-family and spatial validation counts as source-linked local evidence even when the source subcluster aggregate diluted it. An anatomical parent override additionally requires an independently detected target-parent program in that same selected second-round subcluster; neighborhood or morphology alone cannot change identity, and an unsupported override remains `unresolved_biological`. Residual QC limits apply only during final materialization.

For sheep ovary, a post-Atlas follicle-histology failure may trigger one canonical bounded iteration only when all open biological issues are tied to explicit follicle ROI IDs. Rebuild every affected ROI from the selected input's non-SCT raw-count assay through SCT/PCA/SNN/Leiden, score the full catalog, and authorize writeback only for typed failing wall layers. Direct coherent Theca, vascular and mature nonvascular Smooth-muscle identities compete before the generic Stromal remainder. Merge fresh ROI scores into the complete disjoint score ledger and rerun the biological validator on the full membership. A cropped ROI review cannot close Oocyte, other broad classes or whole-section spatial localization.

## Determinism and completion

Use stable ID sorting, a fixed master seed and deterministic per-resolution/per-cohort derived seeds. Identical input, contract and seed must yield identical partitions, membership and semantic hash. Completion additionally requires biological equivalence under neighboring resolution, PCs +/-5 and k +/-20% perturbations.

`validate_lineage_controller_release.py` verifies:

- no first-pass formal membership;
- exactly one second-round cohort per initial cluster;
- project-local raw-count ancestry and full-catalog evidence;
- local splitting only behind a recorded mixed-subcluster trigger;
- exact disjoint broad coverage at merge;
- post-merge-only Atlas, fine and QC authority;
- one raw-count, typed, bounded follicle iteration with full-membership revalidation;
- typed residual reasons, completeness closure and final semantic hash.
