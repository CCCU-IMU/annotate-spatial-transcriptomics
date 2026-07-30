# Staged second-round annotation controller

`run_lineage_controller.py` is the only formal entry and the only chain that may ultimately create release membership. The controller binds its scripts, candidate catalog, profiles, inputs, seed and artifact roles in the annotation contract.

Release-critical scoring, writeback, resolution-selection and residual-QC defaults come from `references/controller_thresholds_v2_2.json`. The project contract freezes and hashes that registry; prose and project-local scripts cannot redefine hidden defaults.

## Authority by phase

1. `whole_tissue_partition` creates only a stable initial-cluster partition and provisional cohort plan.
2. `cluster_cohort_recluster` independently reclusters every initial cluster from project-local non-SCT raw counts and freezes full-catalog subcluster evidence.
3. `local_mixed_subcluster_split` resolves only second-round subclusters with two or more separable competing identity cores.
4. `merge_and_freeze_broad` merges mutually exclusive second-round outcomes and creates formal broad membership.
5. `atlas_and_completeness_review` maps all cells at broad level, lets Atlas rescue only unlabeled cells, performs bounded biological corrections, and then audits every context-evaluable broad on both present-label precision and outside-label recall.
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

Observation-level splitting is a conditional local operation, not an all-object classifier. A candidate being detected or `local_split_worthy` does not make the subcluster mixed. Against a clear specific parent, another specific lineage must have independent subcluster-level multigene, DEG/pseudobulk and neighboring-resolution evidence, and both candidates must contain material pairwise-exclusive direct identity components. A shared or nearly nested program is a watch or catalog-declared contextual lineage-of-origin trace. Inside a supported generic remainder, one specific direct identity component plus a material generic complement may open the check even when the minority whole-subcluster DEG is diluted. Test neighboring Leiden partitions before candidate-local spatial observation components. Two independently strong but inseparable coexpressed programs remain `unresolved_biological`; a single weaker trace does not block a supported parent return. Write a subset only after discriminators, candidate-local pseudobulk/DEG and a validated expression/spatial component demonstrate separable members.

Each candidate independently proposes a spatial/expression subset. Resolve overlap with normalized evidence, pairwise discriminators and anti-programs. Specific lineages precede generic Stromal remainder. Ambiguous overlaps remain with a supported common parent or `unresolved_biological`; candidate ordering cannot assign them.

Group thresholds validate a proposed component and never admit individual observations. Exact remainder closure is scoped to the source subcluster, preserves original scores, permits one additional extraction and never converts “not selected” into QC.

## Broad freeze and final review

After every initial-cluster cohort closes, merge whole-subcluster returns, supported local subsets and local parent remainders. Require an exact disjoint cover of `analysis_set`; unresolved conflicts stay biological. Only this phase freezes formal broad membership.

Fine proposals become labels only under the matching frozen broad parent. State annotations remain separate from identity. The single all-cell Atlas pass occurs after broad freeze: it may rescue unlabeled moderate-or-higher, non-OOD, profile-compatible cells, but cannot silently overwrite a defined broad or create fine labels. Material conflicts reopen the complete source cohort once for query-derived evidence review.

Present/zero-census lineages, embedded programs, unmodeled programs and the complete parent-by-fine catalog are audited after merge. A post-merge component that passed its own strict multi-family and spatial validation counts as source-linked local evidence even when the source subcluster aggregate diluted it. An anatomical parent override additionally requires an independently detected target-parent program in that same selected second-round subcluster; neighborhood or morphology alone cannot change identity, and an unsupported override remains `unresolved_biological`. Residual QC limits apply only during final materialization.

After unresolved-only closure, perform one catalog-wide double-sided review for every context-evaluable broad. One broad cell type is the primary decision unit: its existing membership is checked for precision while all other observations are searched for recall in the same review. Original second-round source subclusters are internal evidence units. Recall uses a direct multi-family spatial component or a source-subcluster group watch when selected-resolution multichannel evidence predicts a lineage deficit that sparse observation scores could not localize. The latter requires targeted raw-count review and cannot assign an entire source group. The audit has no label authority. Exact cell-type decisions may be applied for at most two rounds, and each round must be followed by a complete repeat audit. A retained cell type is reusable only when its reviewed cell-ID scope signature is unchanged. Context-ineligible candidates are `not_evaluable`, never positive completeness evidence or Atlas/fine/final release targets.

The fixed user-facing stage name is **逐大类全样本复核**. “Catalog-wide” is an internal artifact identifier, not a different biological phase; terms such as “目录复核”, “目录审阅决策”, “来源组复核” and “强制复核重放” must not appear as substitute stages. Exactly one broad is active at a time and its evidence packet binds current-member precision, whole-query recall, molecular identity, spatial consistency and literature-boundary consistency to the same membership signature. The packet must include raw-count marker evidence, target-versus-outside and challenger pseudobulk/differential evidence, exact challenger memberships and a spatial plot. Keyword-bearing rationale cannot close a review. A challenger can be refuted only by that bound multichannel evidence; a confirmed precision or recall defect requires an exact patch or bounded targeted review. A patch reopens the changed target and any other broad whose own members/recall/watch signature changed, not unrelated closed types.

For sheep ovary, the external Atlas is the immutable contract-bound GSE233801 split-wall v2 bundle. Its independent Endothelial, Pericyte/mural and Smooth-muscle prototypes were reconstructed from the original reviewed reference clusters rather than obtained by renaming a mixed vascular prototype. The capability matrix and current-query classwise calibration are enforced before routing. Controlled pauses carry an ordered membership-transform ledger and resume state, so subsequent single-type reviews reuse frozen Atlas, unresolved, ROI and clustering products instead of rerunning them; the old merged v1 bundle is accepted only for such an already-frozen resume.

Context is evaluation permission, not identity evidence. Missing, conflicting, `refuted` or `not_evaluable` context removes every release path for that candidate and its context-dependent fine children. Such candidates remain visible in the audit but cannot make a zero census positive. Conversely, an already released context-ineligible broad is a hard inconsistency that must return to its exact source membership; it cannot be excused by fine-program counts.

For sheep ovary, a post-Atlas follicle-histology failure may trigger one canonical bounded iteration only when all open biological issues are tied to explicit follicle ROI IDs. Rebuild every affected ROI from the selected input's non-SCT raw-count assay through SCT/PCA/SNN/Leiden, score the full catalog, and authorize writeback only for typed failing wall layers. Direct coherent Theca, vascular and mature nonvascular Smooth-muscle identities compete before the generic Stromal remainder. Merge fresh ROI scores into the complete disjoint score ledger and rerun the biological validator on the full membership. A cropped ROI review cannot close Oocyte, other broad classes or whole-section spatial localization.

A present broad review cannot close until that one cell type has explicit current-member precision, whole-query recall, molecular identity and whole-section spatial-consistency conclusions. Internal source groups and spatial components provide evidence and patch bounds only.

## Determinism and completion

Use stable ID sorting, a fixed master seed and deterministic per-resolution/per-cohort derived seeds. Identical input, contract and seed must yield identical partitions, membership and semantic hash. Completion additionally requires biological equivalence under neighboring resolution, PCs +/-5 and k +/-20% perturbations.

`validate_lineage_controller_release.py` verifies:

- no first-pass formal membership;
- exactly one second-round cohort per initial cluster;
- project-local raw-count ancestry and full-catalog evidence;
- local splitting only behind a recorded mixed-subcluster trigger;
- exact disjoint broad coverage at merge;
- post-merge-only Atlas, fine and QC authority;
- context-ineligible candidates cannot enter scoring decisions, Atlas rescue, bounded review writeback, fine materialization, completeness positives or final release;
- every context-evaluable broad has a closed catalog-wide precision/recall review after the last membership change;
- one raw-count, typed, bounded follicle iteration with full-membership revalidation;
- typed residual reasons, completeness closure and final semantic hash.
