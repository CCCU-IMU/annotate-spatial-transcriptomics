---
name: annotate-spatial-transcriptomics
description: Independently annotate spatial transcriptomics or single-cell RNA-seq with a provisional whole-tissue partition followed by label-blind second-round cohort reclustering, local mixed-subcluster splitting, post-merge broad freezing, fixed calibrated Atlas review, strictly sequential per-cell-type whole-query review, open-lineage discovery, reliable fine labels and auditable reports. Seurat/R-first is preferred for full-feature RDS inputs; AnnData/H5AD, SingleCellExperiment, BANKSY, Scanpy/Leiden and external cluster tables are supported. Use for end-to-end annotation, evidence review, unknown-lineage discovery, repeatability testing and release reporting. Never copy labels or selected parameters from example projects.
---

# Annotate Spatial Transcriptomics

Use the contract-bound staged controller. The first whole-tissue clustering only partitions the query into second-round cohorts; it never writes a formal biological label. Formal annotation is driven by independent second-round reclustering of every initial cluster.

Read `references/direct-lineage-controller.md`, `references/iterative-controller.md`, `references/quality-standard.md`, `references/state-schema.md` and the applicable biological profile before analysis.

For formal membership, always execute `scripts/run_lineage_controller.py`. Project-local scorers, subset writers and remainder closers are experimental and cannot receive release authority. Artifacts registered as `failed_diagnostic` may be summarized for failure analysis but may never contribute membership, labels, Atlas training, fixture truth or blind-regression input.

## Standard staged workflow

Run these phases in order:

1. `whole_tissue_partition`
2. `cluster_cohort_recluster`
3. `local_mixed_subcluster_split`
4. `merge_and_freeze_broad`
5. `atlas_and_completeness_review`
6. `materialize_final_release`

The annotation contract binds the controller and dependency hashes, candidate catalog, workflow/biological profiles, input snapshot, `full_object`, `analysis_set`, `excluded_initial_qc`, project-local raw-count assay, seed and artifact roles. Historical labels remain invisible until all second-round candidate scores and broad membership are frozen.

### 1. Freeze inputs

- For a preprocessed sheep-ovary SCT+BANKSY Seurat RDS, start with `scripts/bootstrap_sct_banksy_project.py`. It creates the project, infers an unambiguous common BANKSY grid (or requires an explicit mapping), audits the project-local non-SCT raw-count assay and coordinates, freezes the analysis boundary and builds one canonical contract. Do not replace it with project-local audit or contract-adapter scripts.
- Inspect and hash the object, assays/layers, coordinates, observation IDs and upstream clustering grid.
- Freeze exact, disjoint `analysis_set` and `excluded_initial_qc` memberships.
- For SCT+BANKSY input, every second-round cohort must restart from the same project's non-SCT raw counts, preferring `RNA` then `Spatial`. Never SCTransform SCT corrected counts and never reuse a derived expression object from another annotation project.
- Only objectively zero-count or damaged observations may enter `excluded_initial_qc`. Biological uncertainty remains `unresolved_biological` until final closure.

### 2. Whole-tissue partition

Select a stable whole-tissue resolution that preserves the major spatial compartments without obvious technical fragmentation. For sheep ovary SCT+BANKSY, inspect the complete `0.2,0.4,0.6,0.8` grid. Do not optimize this phase to separate every lineage.

For every initial cluster record only:

- exact membership and one second-round cohort ID;
- `provisional_broad`, `mixed` or `unknown`;
- candidate programs, competitors and missing-lineage `watch` signals.

The parent broad label is provenance, never a search-space restriction; always watch and carry forward coherent lineage programs that do not separate during the first partition.

Observation direct/local scores may describe programs and mixtures, but cannot write broad, fine or QC membership here. The only membership output is `whole_tissue_cohort_plan.tsv`. First-pass candidate-local splitting, exact remainder closure, fine labels, residual-QC gates and `release_membership` are forbidden.

### 3. Second-round annotation

Every initial cluster enters its own query-only cohort, including apparently pure clusters. Run:

`project-local raw counts -> SCT v2/glmGamPoi -> PCA -> SNN -> Leiden grid`

For sheep ovary use the complete `0.1,0.2,0.3,0.4,0.6` grid, with PCs and k adapted to current cohort size and complexity. Select the resolution that best separates stable identity programs and anatomy while avoiding depth, state or pure-ECM fragmentation. Compare the nearest lower and higher resolutions; at a grid edge use the two nearest candidates.

Do not rerun SCT/PCA/SNN/Leiden for a scorer, taxonomy, writeback, report or other semantic-only repair when the frozen input SHA, exact cohort-membership SHA, grid, resolution contract, seed and clustering-script SHAs are unchanged. Reuse only the controller-generated `derived_partition_cache` through `--reuse-recluster-manifest`; the controller recomputes full-catalog scoring, resolution selection, adjudication and every downstream membership. Any input, membership, graph/preprocessing parameter, seed or clustering-code change invalidates the cache. A cache has no label authority, and a `failed_diagnostic` cache is forbidden.

Before reading the first-pass provisional record, score every selected-resolution subcluster against the complete open candidate catalog and stable unexplained programs. Freeze those scores, then classify the outcome as:

- `parent_return`;
- `cross_lineage_return`;
- `missing_broad_reconstruction`;
- `fine_candidate`;
- `state_annotation`;
- `unmodeled_lineage_candidate`;
- `unresolved_biological`.

An underpowered cohort is `underpowered_not_evaluable`; it cannot silently inherit the provisional label.

The main annotation unit is the second-round Leiden subcluster. A high-purity subcluster with coherent multi-family identity, DEG/pseudobulk support, reasonable spatial morphology, no independent competitor and neighboring-resolution stability returns wholesale. Sparse, low-RNA, noncontradictory members inherit that broad identity rather than becoming QC through dropout.

Bind exogenous sample context such as reproductive stage in the annotation contract and propagate it automatically through cohort adjudication and every triggered local split. Context only permits evaluation of a context-gated candidate; it never supplies an identity score or membership. For luteal candidates, granulosa- and theca-lineage-of-origin programs are expected competitors/support and must not be treated as hard anti evidence. Stage compatibility still requires an independent luteal program and corpus-luteum-like morphology before release.

Optimize the second round for stable broad identity, not subtype count. When a subcluster lacks an independent child discriminator, keep its supported parent broad and leave the fine label empty. If another broad lineage is coherent, return it cross-lineage at broad level even when no reliable fine subtype exists. A high burden of small marker-poor or unresolved subclusters triggers review of the nearest lower resolution; do not manufacture fine labels to close the cohort.

A shared generic parent/remainder program remains visible in the audit but is not an independent competitor to a specific identity. Conversely, a weak specific signal cannot preempt a well-supported generic parent unless it has its own orthogonal identity evidence. A bounded identity-dominant whole-subcluster route may tolerate modest dropout-associated contradiction only when no independent competing lineage exists; it does not relax local subset validation.

### 4. Evidence and local mixed-subcluster splitting

Use full-feature expression. Scale each gene by its nonzero query 95th percentile; aggregate a family with the strongest two available genes; combine direct/local evidence at `0.35/0.65` by default. One anti-marker or local anti-signal is a penalty. Only a direct, multigene contradiction spanning independent anti-families is hard. Aggregate winners, shared parent markers, generic ECM and `ACTA2/TAGLN` alone cannot suppress or prove a specific lineage.

Load release-critical defaults from `references/controller_thresholds_v2_2.json`; do not duplicate them in project scripts or infer them from prose. The annotation contract binds a frozen copy. Project-specific overrides are permitted only when explicitly recorded in that contract.

Candidate visibility is not biological mixedness. Start one bounded local separability check only after one of two conditions is proven: (a) two specific broad candidates each have an independent subcluster-level identity program and both contain a material pairwise-exclusive direct identity component; or (b) one orthogonally supported specific direct component and a material complement are embedded in a supported generic remainder such as Stromal. Shared steroidogenic, ECM, contractile or locally smoothed signal is a `watch`, not a competitor. When a clear specific parent already exists, a weaker challenger cannot reopen it without its own subcluster-level multigene, DEG/pseudobulk and neighboring-resolution support. A specific minority inside a generic parent may still enter the check when whole-subcluster DEG is diluted, but direct identity—not one background gene—must define its seed. First test the already computed neighboring Leiden partitions; only if they do not resolve the identities may candidate-local spatial observation components be used. If two independently strong programs remain coexpressed without mutually exclusive components, do not split or choose by candidate order: keep the affected observations `unresolved_biological`. If only the parent is independently supported, return the parent and record the alternative as `watch`. Actual writeback still requires observation discriminators, candidate-local DEG/pseudobulk and a validated expression/spatial component. A high contradiction fraction in an aggregate true mixture must not erase split candidates; contradiction is evaluated again on each proposed subset. Generate candidate subsets independently, resolve overlaps with normalized evidence and pairwise discriminators, and resolve specific lineages before a generic Stromal remainder. A stable 0.5%-level trace is not authority to route the subcluster through P41. Before submitting local-split jobs, run `audit_local_split_workload.py` across all completed second-round cohort manifests and pass its PASS artifact through every local-split call. If more than 50% of the analysis set is pending, revisit second-round resolution or mixed-trigger semantics instead of submitting an all-object-equivalent P41 workload.

The `0.70 support / 0.30 competitor margin / <=0.05 contradiction` criteria validate an already formed group; they are never per-observation admission thresholds. Exact remainder closure is local, permits at most one additional extraction, and never equates failure to enter a subset with QC.

After both bounded specific-lineage extractions, return an exact tail to its generic broad parent when the tail itself has at least 70% coherent multigene parent support, no residual separable specific component and a stable parent program. Do not reuse the generic parent's mixed-source anti fraction as a veto after candidate-specific blockers have been removed; only truly ambiguous residual component overlap remains `unresolved_biological`.

When exact-remainder review has already certified an orthogonal cross-lineage identity as a blocker, keep that blocker visible even if the still-mixed remainder has more than 5% aggregate contradiction. The aggregate contradiction ceiling cannot be applied a second time to erase the competing identity and force a wholesale parent return. State programs are recorded separately and do not, by themselves, veto a supported broad parent.

The one additional extraction can only re-evaluate candidate-local components already generated inside that mixed subcluster after accepted members are removed. It cannot scan the whole remainder as a new per-cellbin classifier. Candidate-specific exceptions, such as an Oocyte canonical identity component, must be declared in the catalog/controller and revalidated from raw component evidence; they cannot become general anti-marker bypasses.

### 5. Merge, freeze and review

After all cohorts finish, merge whole-subcluster returns and supported local subsets, require an exact mutually exclusive cover of `analysis_set`, resolve conflicts, then freeze formal broad membership. First-pass provisional labels are not a release source. Fine proposals are materialized only inside the correct frozen broad parent. State annotations occupy a separate column.

A fine candidate may propose broad reconstruction only when its exact proposed members independently retain a parent-broad identity core in at least 25% of members, have no coherent hard contradiction and the fine candidate itself passes its two-family identity/discriminator test. A parent broad candidate that declares required positive families must satisfy every declared family on the exact subset; a parent without such a declaration must still contribute at least one supported broad family. Spatial enrichment is audited but is not required for a widespread true parent background. A synthetic `parent_identity` field, overlapping parent/fine genes or a strong fine discriminator cannot substitute for this exact parent-broad validation. If the parent test fails, reject the promotion and return the members to local remainder review; never expand the surrounding cluster.

Run one calibrated all-cell broad Atlas mapping only after broad freeze. For the sheep-ovary profile, the immutable active reference is the contract-bound GSE233801 split-wall v2 bundle. It reconstructs independent Granulosa, Immune, Stromal/mesenchymal, Endothelial, Pericyte/mural and Smooth-muscle prototypes from the original reviewed reference clusters; Epithelial/mesothelial and Theca are challenge-only, while Oocyte and Luteal are unsupported. `supported` never bypasses current-query classwise calibration. Atlas may directly rescue an unlabeled, non-OOD, profile-compatible observation at moderate-or-higher calibrated confidence; challenge-only classes can reopen biological review but cannot write labels. Atlas cannot silently overwrite an existing broad label or create a fine label. A material conflict reopens the complete source subcluster/cohort once for query-derived biological review. Runtime Atlas substitution, rebuilding or query-reference joint training is forbidden; the deprecated merged v1 bundle is resume-only.

After Atlas routing, first run one bounded post-merge review on observations still marked `unresolved_biological`. It is not a whole-object per-cellbin classifier. Each releasable catalog candidate forms strict direct multi-family seeds inside its original second-round subcluster, expands only through coherent same-candidate support on that subcluster's spatial graph, and requires at least three direct seeds in a five-member component. Resolve specific lineages before a catalog-declared generic remainder; when two specific component proposals remain inseparable, preserve unresolved rather than using candidate order. Catalog-declared anatomical parent overrides may return a component to a stage-defined broad parent—for example, a Theca lineage-of-origin program embedded in a coherent corpus-luteum compartment returns `Luteal`—only when that same selected second-round subcluster independently supports both a steroidogenic core and an independent corpus-luteum identity family. Neighborhood, morphology, PGR/APOD/SFRP4 or generic steroidogenesis alone cannot authorize the override; if the anatomical context passes but target-parent identity support is absent, keep the component `unresolved_biological`. The accepted parent candidate replaces the lineage-of-origin challenger in release provenance so completeness can verify the exact source support. This stage exists to close bounded low-RNA tails and residual mixed components, not to tune the final QC percentage.

After any contract-required tissue-specific biological review and bounded repair, run `audit_catalog_wide_lineage_challengers.py` across every context-evaluable broad lineage, including labels already present and labels currently absent. The **decision, evidence-packet and reporting unit is one broad cell type**. The controller exposes exactly one active target and the Agent must first say `现在开始对 <cell type> 进行专项复核。`; it may not generate or decide packets for the remaining queue in the same invocation. For each annotated broad, one review simultaneously checks **precision** inside the current label and **recall** across all other query observations, then records molecular, spatial and literature-boundary conclusions before producing one retain-or-patch decision. Original second-round source subclusters, direct multi-family spatial components and source-subcluster `group watch` records are child evidence and exact patch bounds, not separate user-facing review tasks. A group watch is raised when selected-resolution DEG/pseudobulk/multichannel evidence predicts more lineage members than the sparse observation scorer recovered; it authorizes targeted raw-count review only and cannot label the whole source group. For a zero-census lineage, direct raw-count support from at least two independent families, at least three genes and a multigene identity-core family creates a bounded source-subcluster challenger even when its cells are too fragmented to form the generic five-member spatial component; this only opens review and never writes labels. Low-RNA members may inherit only inside a supported current-label group, and generic Stromal is not allowed to generate whole-section recall components or group watches. The audit writes no labels. Resolve the active type with query raw-count marker families, target-versus-outside DEG/pseudobulk, pairwise competitors, whole-section spatial distribution, literature-derived alternative identities, over-recall and under-recall; geometry or literature alone has no writeback authority. A raw challenger may be explicitly refuted by the bound multichannel packet; confirmed over-recall or under-recall requires an exact patch or bounded targeted review. The two-round limit applies independently to each cell type. An unrelated type's patch does not reopen a closed type. Moving members into a closed type or creating a new recall/watch challenger reopens it; removing at most 10% of its members can close by deterministic monotonic-subtraction re-audit only when the retained set is an exact subset and no new recall, watch or zero-census challenger appears. Stable partitions are reused; rerun SCT/PCA/SNN/Leiden only when that cell type's existing source units are genuinely not evaluable.

In user-facing progress, this stage has exactly one name: **逐大类全样本复核** (`per-broad whole-query review`). Do not call it “目录复核”, “目录审阅”, “来源组复核”, “强制重放” or any near-synonym. Every queued broad has one controller-built active packet binding current-member precision, whole-query recall, raw-count/full-transcriptome molecular identity, whole-section spatial consistency and the literature boundary/alternative explanation. Free text, keyword presence or a project-local decision script cannot substitute for that packet. If the packet itself provides sufficient multichannel evidence to reject a marker-only, ambient or shared-program question, record that refutation and close without a gratuitous membership patch. After a patch, regenerate the affected type's packet and reopen every other type whose own member/recall/watch signature changed; preserve unrelated closed reviews. For sheep ovary, bind the Oocyte canonical endpoint and follicle-ROI endpoint into the relevant packets before closure.

A present broad review is incomplete unless it explicitly records five conclusions for that one cell type: current-member precision, whole-query recall, molecular identity, whole-section spatial consistency and literature-boundary consistency. Source subclusters and components may justify or bound a patch, but they cannot close these broad-level conclusions on their own.

Every membership-changing step is recorded in `membership_transform_chain.json` with exact source/result semantic hashes, cell universe, delta cell IDs and evidence manifest. A retain/absence decision with zero changed observations is a biological closure, not a membership transform: it reuses the exact source membership and is not appended to the transform chain. Completion reads this generic ordered ledger rather than assuming a fixed list of repair filenames. Use `--resume-review-manifest <previous lineage_controller_manifest.json>` after a controlled per-cell-type pause; this resumes from the frozen membership/transform chain and must not rerun Atlas, unresolved rescue, ROI repair or either clustering round. The controller reuses one contract-bound raw-count/coordinate cache across these pauses and recomputes full-transcriptome pseudobulk only for the active type and its evidenced competitors.

Map the complete `analysis_set`, retain explicit out-of-distribution handling, and never form a full query-by-reference distance matrix; use the contract-bound fixed representation and ANN index.

Then audit present broad labels, zero-census candidate lineages, embedded programs in large labels, unmodeled stable programs and the complete broad-by-fine candidate matrix. A whole-subcluster return must trace to selected-resolution multichannel evidence; a local writeback to canonical `subset_validation.tsv`; an exact parent tail to canonical `exact_remainder_audit.tsv`; a post-merge unresolved return to its component/context audit; and a catalog-wide correction to its exact review queue, evidence decision and bounded membership change. Context-ineligible candidates remain visible as `not_evaluable` but cannot count as positive recall, block a justified zero census, receive Atlas rescue, materialize fine labels or enter final release. Never invalidate a supported local subset merely because its minority program was diluted in the source subcluster aggregate. Stable unmodeled programs require an exact biological review: route them to an existing catalog identity, a state/technical program, an insufficient-identity program, or a genuine novel-lineage candidate. Only the last outcome blocks release for catalog extension and cohort rerun; none may be auto-named from a marker list. Convert remaining unresolved observations to typed QC only during final materialization. `QC <10%` and `QC <50,000` are formal-release requirements, never first-pass tuning targets. After one Atlas route, one bounded post-merge unresolved review, the optional single ROI repair and at most two per-broad decision rounds, stop automatic membership recovery. If either QC bound is still reached, materialize a non-release `PENDING_USER_REVIEW_HIGH_UNRESOLVED` membership/report with typed reasons; do not launch a whole-residual, QC-anchor or repeated unresolved reclustering loop merely to cross the release threshold.

Canonical-cluster identities use one complete detector across second-round competition, zero-census review and release. A weak ordinary/rare program may remain a `watch`, but it cannot make a canonical zero census positive unless the full canonical-cluster rule passes. For Oocyte this means coherent multigene germline/ooplasm evidence, direct core, DEG contrast and object-level spatial review; low-level zona or ambient signal alone cannot block a justified absence call. If ordinary second-round scanning returns zero Oocyte while the label-blind full-analysis-set multi-module starting gate is nonzero, run exactly one query-only canonical targeted cohort containing every starting-gate observation across initial clusters. Freeze its cluster decision before exposing existing labels; return the complete passing cluster except objective input-QC failures or direct multigene somatic hard contradictions. Spatial foci, zona signal and strict seeds never define or shrink this membership.

When a canonical-cluster challenger passes but the census remains zero, generic completeness may mark it `not_evaluable` only if the bound profile requires a downstream object-level biological-quality review. That review must then either justify absence, reconstruct a canonical object, or return the project for iteration. Deferral is not an absence call and is forbidden when no downstream canonical review is contractually required.

For sheep ovary, run `validate_sheep_ovary_biological_quality.py` before the final catalog-wide review closes. Require three compact biological endpoints: complete-label spatial localization, canonical Oocyte credibility and stage-appropriate follicle ROI histology. Read `references/profiles/sheep_ovary_biological_quality_review.md`. A coherent lineage program recurring across a follicle boundary but retained in generic Stromal/unresolved returns only that exact ROI or source subcluster to one bounded raw-count review. Geometry is a challenger and cannot write labels.
When rerunning a triggered follicle ROI, bind its pre-iteration review with `--expected-roi-review`. Freeze layers that already passed; only typed failing layers may receive writeback authority. Reject any candidate that erases a previously confirmed antral cavity or converts it into `not_evaluable`.
When fresh scores cover only the reopened ROI, supply the complete-section x/y ledger through `--coordinate-membership`; it is geometry-only and cannot provide labels or writeback evidence.
Within a reopened follicle wall, distinguish layer detection from writeback: local/coherent signal may prove a layer exists, whereas directly discriminated identity cores or a validated expression subcluster authorize membership. Formal Theca interna is steroidogenic/androgenic; structural/perifollicular Theca remains exploratory. Resolve Theca, Endothelial, Pericyte/mural and mature nonvascular Smooth muscle independently before the generic Stromal remainder, and never demand mutually exclusive labels for the same mixed cellbin. A missing follicle-wall layer is a negative/NOT_EVALUABLE audit; only coherent direct under-recall, unsupported published identity or implausible published morphology reopens an ROI.
The canonical controller may execute exactly one automatic ROI iteration only when every remaining biological-quality problem is a typed follicle ROI problem. It rebuilds each bounded ROI from project-local non-SCT raw counts with SCT/PCA/SNN/Leiden, rescans the full catalog, and writes only directly coherent identity cores within the authorized wall layers. It then replaces ROI scores inside the complete-section score ledger and reruns all three biological endpoints on the full membership. A subset-only PASS has no completion authority.

### 6. Final annotation and reporting

Every analysis-set observation receives exactly one frozen broad label or typed retained state. Fine labels require high confidence, an independent discriminator and a matching frozen parent. Hypoxia, proliferation, atresia, luteinization, stress, cell cycle and similar programs remain state columns unless they form an independently stable identity.

Write the frozen membership back with `write_frozen_annotations_to_seurat.R`; analysis membership plus the optional excluded-initial-QC ledger must exactly cover the Seurat object. Every semantic repair or proposal overlay must first pass through `apply_cell_id_membership_patch.py`: proposal IDs are unique and contained in the frozen base, updated values are joined by `cell_id`, base order is restored and non-proposal members remain byte-equivalent in the updated columns. Positional row assignment is forbidden. Materialize the single public `final_cell_type`: high-confidence parent-locked fine identity overrides broad, otherwise broad is shown, and terminal unresolved is `QC/Unknown`. Build the sample-agnostic HTML with `build_frozen_review_report.py`. Before user confirmation, use `pending_user_review`; after confirmation, use `approved_final` without changing membership. Public census, DEG, dotplot, overview and fixed-point-size highlights use only `final_cell_type`; broad/fine remain internal provenance. For Seurat output, normalized dotplot expression may use the active expression assay, but absolute detection must use the non-SCT raw-count assay bound by `query_input_audit`; never infer the count source at report time.

Treat `current_stage.json` and `next_action_manifest.json` as the only scheduler-facing phase state. `REVIEW_REQUIRED` and `ITERATION_REQUIRED` are controlled, resumable pauses and must exit the wrapper successfully; `FAILED_RUNTIME` is the only failure state. The final broad/fine/state/QC materialization is appended to the ordered membership-transform chain, so terminal QC and public `final_cell_type` cannot bypass provenance.

## Open taxonomy and sheep-ovary rules

Profiles provide candidate lineages, marker families, anti-programs, context requirements and literature provenance; they are not label maps. Every second-round subcluster scans the full catalog, and stable catalog-external programs become `Unmodeled lineage candidate` rather than being forcibly named.

For sheep ovary, read:

- `references/profiles/sheep_ovary_standard_workflow.md`
- `references/profiles/sheep_ovary_candidate_lineage_catalog.json`
- `references/profiles/sheep_ovary_literature_2025_2026.md`
- `references/profiles/sheep_ovary_rfirst_case_reference.md` only as a sanitized failure/regression guide
- `references/profiles/sheep_ovary_oocyte_rfirst_route.md` before an Oocyte decision
- `references/profiles/sheep_ovary_biological_quality_review.md` before biological-quality approval

Release Endothelial, Pericyte/mural and Smooth muscle as independent competing broad identities. Lymphatic endothelial is a high-confidence fine identity under Endothelial and may replace Endothelial only in public `final_cell_type`; the legacy `Vascular-associated` label cannot be released. Endothelial requires a junction backbone plus independent angiovascular support; Pericyte/mural requires a mural backbone plus independent support; Smooth muscle requires a mature nonvascular MYH11/CNN1/ACTG2/SMTN/LMOD1-centered program. ACTA2/TAGLN/MYL9 or vascular proximity alone proves neither mural nor smooth-muscle identity. Epithelial/mesothelial cannot expand from one keratin, one spatial component or a mixed-cluster aggregate winner. Oocyte uses a complete canonical-cluster rule; zona/ambient signal does not expand the census. Granulosa state programs do not substitute for identity subtypes. Discover Theca from its multigene steroidogenic/androgenic program across the complete candidate space; follicle ROI and distance are post hoc challengers, never admission or exclusion gates. Release `Luteal` only when stage, steroidogenic function, an independent corpus-luteum identity family and mass-like morphology agree; generic `steroidogenic` is a program, not a broad label. Mesenchymal progenitor-like and Neuroendocrine-like remain exploratory. No designated regression sample label or membership is a runtime prior or Atlas.

## Execution and validation

Use local execution for discovery and small audits. Submit heavy data work through the user's compute workflow. The contract-bound resource classes are: 1–4 CPUs for state, manifest, hash and preflight work; 8 CPUs by default and no more than 16 for per-broad marker/DEG/pseudobulk evidence; and 64 CPUs for heavy SCT/PCA/SNN/Leiden or large-RDS materialization. Match requests to real parallelism; independent cohorts/resolutions are the primary parallel units. A 64-CPU job must still use `--scoring-workers 1` for any observation-scoring boundary with at least 100,000 observations, because scoring forks replicate large direct/local arrays.

Tests must cover:

- first-pass inability to write formal broad/fine/QC membership;
- exactly one second-round cohort per initial cluster;
- non-SCT raw-count ancestry and provisional-label blindness;
- full-catalog scanning and missing-broad reconstruction;
- local-only observation splitting and non-QC remainder behavior;
- post-freeze fine/Atlas/final authority;
- context-gated candidates cannot leak into resolution selection, unresolved review, completeness, Atlas rescue, fine materialization or final release;
- every context-evaluable broad receives post-Atlas present-label precision and outside-label recall review, with at most two exact decision rounds;
- follicle ROI repair only after a typed post-Atlas review, from non-SCT raw counts, with one bounded iteration and a full-membership re-review;
- failed-diagnostic isolation;
- deterministic repeat runs and parameter robustness;
- Granulosa-rich, Stromal/Smooth/Vascular mixed, rare Oocyte, low-fraction Epithelial and unmodeled-program fixtures.
- sheep-ovary three-endpoint review: diffuse restricted labels fail, Oocyte remains canonical and object-aware, and a multisector Theca/Vascular/Smooth program hidden in a follicle-wall Stromal remainder triggers local iteration without changing membership.

For technical determinism, identical input, contract and seed must yield identical partition, membership and semantic hash. For parameter robustness, preserve every >=1% broad class, broad-census Spearman >=0.95, major spatial Dice >=0.80 and the parent/core program of major fine labels across neighboring resolutions, PCs +/-5 and k +/-20%.

## Conversation output

Report biological decisions rather than gates, hashes, filenames or state-machine chatter. Each update has at most five short bullets:

- whole tissue: selected resolution, major structures, persistent lineage watches and cohort count;
- cohort: selected resolution, major programs and parent/cross-lineage/missing-broad/unresolved outcomes;
- Atlas/completeness: rescue count, material broad conflicts, per-lineage precision/recall questions and whether absent lineages truly lack support;
- final: `final_cell_type` census, spatial structure, retained QC/Unknown and regions needing review.

Only report a validator or execution detail when one concrete failure blocks progress, and include the next action in the same line.
