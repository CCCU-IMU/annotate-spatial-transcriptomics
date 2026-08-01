# Biological quality standard

## Primary endpoint

Optimize correct broad identity, plausible spatial anatomy, faithful uncertainty and reproducible shallow subtypes. Cluster count, subtype count and agreement with an older annotation are not endpoints. Broad-only is correct when the query lacks an independent fine discriminator.

The first whole-tissue pass is judged only as a stable cohort partition. Annotation quality is judged from the second-round subclusters, local mixed-subcluster resolution, merged broad membership and post-merge review.

## Required evidence

Every released broad class must have:

1. coherent multigene marker-family support on the full-feature query;
2. compatible DEG and pseudobulk evidence, including competitor/anti-program review;
3. stability at the selected second-round resolution and its nearest neighbors;
4. plausible spatial morphology and localization;
5. a documented source: whole second-round subcluster, supported local subset, local remainder parent or post-merge unlabeled Atlas rescue.

A fine label additionally requires a matching frozen broad parent and an independent child discriminator. State, depth, ECM, stress, hypoxia, cell cycle or generic contractility alone cannot create a fine identity. Low RNA or weak child markers may make a subtype not evaluable without invalidating a coherent broad identity: return the subcluster to its supported parent broad, or to another supported cross-lineage broad, and leave the fine label empty.

## Method-independent acceptance

- Freeze `full_object`, `analysis_set`, `excluded_initial_qc` and project-local raw-count ancestry.
- Keep historical annotation, repair membership, failed diagnostic artifacts and same-batch labels invisible until new membership freezes.
- Create exactly one second-round cohort per initial cluster, not per provisional broad label.
- Restart each cohort from non-SCT raw counts and scan the complete open candidate catalog.
- Permit the second round to reconstruct a broad lineage absent from the first pass.
- Return a pure second-round subcluster wholesale, including sparse noncontradictory observations.
- Trigger observation-level splitting only for a documented mixed second-round subcluster.
- Recompute exact local remainder without converting nonselection into QC.
- Freeze broad membership only after all cohorts merge into an exact disjoint analysis-set cover.
- Run Atlas and missing-lineage review only after broad freeze. For sheep ovary, new mappings use only the immutable contract-bound GSE233801 split-wall v2 bundle reconstructed from original reference clusters. Supported Granulosa, Immune, Stromal/mesenchymal, Endothelial, Pericyte/mural and Smooth-muscle pairs may rescue an unlabeled broad identity only after current-query classwise calibration; challenge-only Epithelial/mesothelial and Theca can raise review questions but cannot write labels. Atlas cannot overwrite a defined label, create a fine label or substitute another reference. The deprecated merged v1 bundle is resume-only.
- Audit every present broad, every zero-census catalog lineage, every stable unmodeled program and every parent-by-fine candidate.
- After tissue-specific review, audit broad identities strictly one at a time. One invocation exposes one active cell type, one multichannel packet and one decision. The decision must separately close current-member precision, whole-query recall, molecular identity, whole-section spatial consistency and the literature boundary. A patch reopens only types whose own membership/recall/watch signature changed. A small exact monotonic subtraction may close without repeating the full packet only when no new challenger is created; any gain, non-subset change or new challenger reopens the type.
- A zero-census broad cannot close merely because its direct multigene observations are spatially fragmented. A direct two-family/three-gene identity-core challenger triggers bounded source-subcluster review, then requires DEG/pseudobulk and competitor adjudication before any label is written.
- Convert unresolved biological cells to typed QC only at final materialization.

## Failure patterns that block release

- Granulosa, Epithelial, Smooth muscle or another restricted lineage expands diffusely through an anatomically incompatible compartment.
- A mixed subcluster is assigned wholesale from its aggregate winner despite an independent competitor.
- Multiple visible/shared marker programs are treated as mixed without material pairwise-exclusive direct identity components.
- A clear specific parent is reopened by a weak background trace that lacks independent subcluster-level DEG/pseudobulk and stability evidence.
- Observation-level splitting is attempted before the already computed neighboring Leiden partitions are checked.
- Epithelial is inferred from one keratin/surface marker or a single spatial component.
- Smooth muscle is inferred from `ACTA2/TAGLN` without mature nonvascular contractile identity or mural exclusion.
- Endothelial, Pericyte/mural and mature nonvascular Smooth muscle are independently adjudicated; ACTA2/TAGLN-only or irreducibly mixed vascular-wall observations are not forced into any of them.
- Oocyte is expanded by zona/ambient signal instead of a coherent canonical cluster.
- A common ovarian lineage is declared absent despite a repeatable multigene spatial program.
- A stable catalog-external program is silently forced into the nearest known label.
- First-pass provisional labels, failed diagnostics or a historical annotation contribute runtime membership.
- Final QC is at least 10% or 50,000 observations. Before the bounded recovery budget is exhausted this returns only to the contributing second-round/post-merge biological problem; afterward it stops as a non-release high-unresolved review candidate rather than opening another global residual or QC-anchor cycle.

## Stability

Technical determinism requires identical partitions, memberships and semantic hash for identical input, contract and seed.

Parameter robustness is evaluated with neighboring reasonable resolutions, PCs +/-5 and k +/-20%. All broad classes occupying at least 1% must remain present; no new or missing >=1% broad class is allowed; broad census Spearman must be at least 0.95; major broad spatial Dice must be at least 0.80; major fine labels must retain the same parent and core biological program. Rare lineages are reviewed by program and anatomy rather than prevalence alone.

## External comparison

A previous annotation is loaded only after the new result freezes. Compare broad census, follicular/vascular/surface anatomy, major subtype programs and failure patterns. Exact cell counts and cellbin equality are not required. The comparison is an external acceptance test, never a runtime classifier or Atlas.

## Sheep-ovary biological endpoints

For sheep ovary, reduce final biological approval to three required questions:

1. Is the complete membership of every released cell type spatially and molecularly plausible?
2. Is Oocyte supported by a complete canonical group and reconstructable objects without zona/ambient expansion?
3. Where follicles are present, are stage-appropriate follicle ROIs interpretable, and do large/antral candidates show a cavity-bounding Granulosa layer followed by Theca-interna steroidogenesis with interleaved vasculature and an outer fibromuscular/stromal transition?

Run the deterministic sheep-ovary biological-quality review after broad freeze and Atlas. A sample without follicles or without antral follicles may close the corresponding structure as `NOT_EVALUABLE`; absence is not a failure. A coherent multisector lineage program hidden inside generic Stromal/unresolved is an iteration trigger. Reopen only the exact source subcluster or bounded follicle ROI. Spatial shape alone never writes a label.

Only after those three endpoints are reviewed may the controller enter **逐大类全样本复核**. The Agent announces `现在开始对 <cell type> 进行专项复核。`, completes that type and freezes its closure before opening the next type. Shared raw-count and coordinate caches are allowed; batch generation or batch closure of formal cell-type decisions is not.

Runtime efficiency must preserve this biology: reuse one immutable project-local raw-count/coordinate cache, compute full-transcriptome pseudobulk only for the active type and its evidenced competitors, do not rewrite membership for zero-change decisions, and do not append such decisions to the semantic transform chain. Use 1–4 CPUs for state/preflight work, 8–16 for active-type evidence and 64 only for heavy reclustering or large-RDS materialization.

## Main-Agent approval

After computational completion, the main Agent reviews broad spatial reasonableness, marker/anti-marker evidence, missing-lineage audit, fine-parent consistency and retained QC. The user sees a lightweight review report before final release assets are generated.
