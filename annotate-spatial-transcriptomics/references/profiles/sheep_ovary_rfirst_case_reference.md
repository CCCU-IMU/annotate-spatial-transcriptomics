# Sanitized successful sheep-ovary R-first case reference

Use this file only as a regression reference for transferable annotation strategy. It is not a runtime controller, Atlas, label map, parameter preset, or source of expected results. The successful case did **not** perform whole-object observation-level lineage assignment during the first partition.

This sanitized strategy trace records only the algorithmic route to final broad membership; it contains no runtime membership or expected label map.

This reference contains no sample path, identifier, observation membership, historical mapping, selected resolution, expected count, or repair result. If it conflicts with the active annotation contract or standard controller, the contract and controller take precedence.

## Core lesson from the successful case

1. The whole-tissue partition proposed computational cohorts; it did not freeze cell identities.
2. Every initial cluster, including an apparently pure one, underwent query-only second-round reclustering.
3. The second-round Leiden subcluster was the primary annotation unit.
4. Cross-lineage returns, missing broad reconstruction, fine candidates, states, and unmodeled programs were discovered during second-round cohort review.
5. Observation-level splitting and exact remainder closure were used only inside a genuinely mixed second-round subcluster.
6. Broad membership was frozen only after all cohort outcomes were merged; fine labels and Atlas review followed that freeze.

## Phase 1: whole-tissue partition builds cohorts only

- Select a stable whole-tissue resolution from the complete grid using spatial structure, adjacent-resolution stability, biological program preservation, and resistance to technical fragmentation.
- For each initial cluster, record only `provisional_broad`, `mixed`, or `unknown`, plus candidate programs, competitors, and missing-lineage watches.
- Treat `provisional_broad` as a cluster-level planning note. Do not write it to the observation ledger, use it as evidence, or freeze it as membership.
- Precomputed direct or local observation signals may describe within-cluster heterogeneity and watches, but they must not drive formal writeback.
- Assign every analysis observation to exactly one initial-cluster cohort plan.

The first phase must not perform:

- whole-object observation-level lineage classification;
- candidate-local spatial splitting;
- exact remainder closure;
- broad, fine, or QC release writeback;
- residual-QC completion checks;
- release membership or release authority generation.

An apparently pure initial cluster is still provisional. A lineage embedded within it may be recovered in the second round.

## Phase 2: recluster every initial-cluster cohort

- Build one independent cohort for every initial cluster rather than merging clusters by provisional label.
- Read project-local raw counts from the verified non-SCT assay; do not reuse corrected SCT values, a whole-tissue graph, or an expression object from another project.
- Run query-only normalization and graph construction: raw counts, SCT v2 with `glmGamPoi`, PCA, SNN, and a Leiden resolution grid.
- Adapt PCs and neighbourhood size to cohort size and complexity; do not copy whole-tissue BANKSY graph parameters.
- Select a cohort resolution by identity-program separation, DEG and pseudobulk coherence, spatial morphology, adjacent-resolution stability, and avoidance of state- or depth-driven fragmentation.
- Score the complete open candidate catalog for every second-round subcluster while the provisional first-round label remains hidden.
- Read the provisional label only after evidence is frozen, and use it solely to name the outcome.

Allowed second-round outcomes are:

- `parent_return`;
- `cross_lineage_return`;
- `missing_broad_reconstruction`;
- `fine_candidate`;
- `state_annotation`;
- `unmodeled_lineage_candidate`;
- `unresolved_biological`.

A fine candidate is only a parent-bound proposal at this stage. It is not formal fine membership.

## Phase 3: decide subclusters, then split only local mixtures

- Use the second-round Leiden subcluster as the default annotation unit.
- Return a high-purity subcluster as a whole when it has a coherent multigene identity program, compatible DEG and pseudobulk evidence, reasonable spatial morphology, adjacent-resolution stability, and no independent competing lineage.
- Allow sparse, low-detection observations without hard contradiction to inherit the supported broad identity of a high-purity subcluster.
- Do not let an aggregate winner suppress a smaller independent candidate program.
- Treat a single anti-marker or local anti signal as a penalty. Require consistent direct multigene contradictions spanning independent anti-families for hard exclusion.

Trigger observation-level splitting only when one second-round subcluster contains at least two independently supported identities that remain separable by expression or space.

Within that local mixed subcluster:

- generate an identity-core seed and candidate-local spatial components independently for each candidate;
- allow one source subcluster to propose multiple lineage subsets;
- resolve overlaps with normalized evidence, pairwise discriminators, and anti-programs rather than catalog order;
- process specific lineages before a generic Stromal remainder;
- validate support, competitor margin, and contradiction at the proposed-subset level, not as a per-observation admission gate;
- execute exact remainder closure only on the original membership of that local subcluster;
- permit at most one additional subset-extraction pass;
- return a coherent remainder to its supported parent when appropriate;
- retain irreducible biological ambiguity as `unresolved_biological`.

Failure to enter a candidate subset is not a QC definition. Never run this remainder algorithm over the complete first-round object.

## Phase 4: merge and freeze broad membership

- Merge whole-subcluster returns, supported local subsets, and local-remainder parent returns from every cohort.
- Require mutually exclusive membership and exact coverage of the analysis set.
- Recompare evidence for conflicting claims; keep unresolved members biological rather than converting them silently to QC.
- Freeze formal broad membership only after this merge.
- Materialize fine candidates only inside their frozen compatible broad parent.
- Store state programs in a separate field; they do not replace broad or fine identity.

Every released broad member must trace to a second-round whole-subcluster return, supported local subset, local-remainder parent return, or post-merge Atlas rescue. A first-round provisional label is never a release source.

## Phase 5: Atlas challenge and completeness review

- Run one all-cell broad Atlas mapping after broad membership is frozen.
- Directly rescue only unlabeled, moderate- or high-confidence, non-OOD, ontology-compatible observations.
- Use Atlas agreement to close defined labels and Atlas conflict to reopen the complete source subcluster or cohort once; never silently overwrite an existing broad label.
- Do not use Atlas to create fine labels or calibrate the query from its own predictions.
- Audit residual QC reasons, missing and common lineages, embedded programs in large broad classes, spatial plausibility, unmodeled programs, and the complete broad-by-fine candidate matrix.
- Apply residual-QC release thresholds only at the final completion stage. If they fail, return to the implicated second-round cohort, not to first-round global per-observation threshold tuning.

## Transferable ovarian safeguards

- Keep the candidate catalog open. A sample may contain stage-, breed-, or region-specific luteal, neural, glial, epithelial, vascular, or other lineages absent from another sheep ovary.
- Keep blood endothelial, lymphatic endothelial, and pericyte/mural identities under `Vascular-associated`.
- Reserve `Smooth muscle` for a mature nonvascular contractile program. `ACTA2` or `TAGLN` alone does not separate it from vascular-wall or contractile stromal cells.
- Use the canonical Oocyte-cluster route: establish identity with independent maternal and non-zona programs plus somatic contradiction review, then retain compatible canonical-cluster members despite sparse detection. If ordinary second-round scanning is negative but a label-blind full-analysis-set multi-module starting gate is nonzero, run one query-only targeted cohort containing every starting-gate observation before calling absence.
- Never expand Epithelial/mesothelial from one keratin, one surface marker, or one spatial component to an entire mixed subcluster.
- Separate Granulosa identity from fine subtype and state. Hypoxia, proliferation, atresia, and luteinization are states unless an independent stable identity program exists.
- Release steroidogenic or androgenic Theca as broad when supported by its molecular program; use follicle ROI and distance only as post hoc anatomy checks, never as candidate membership gates. Keep structural or perifollicular programs provisional unless they establish an independent identity.
- Release `Luteal` only when stage, a steroidogenic core, an independent corpus-luteum identity family and a mass-like compartment agree. Generic steroidogenesis or PGR/APOD/SFRP4 alone remains a program, not a broad identity.
- Record a stable unmatched multigene, multiresolution, spatially coherent program as `Unmodeled lineage candidate`; do not force it into the nearest known class.

## Historical-result and release boundary

- Keep all historical labels, repair memberships, expected counts, selected resolutions, and same-batch results invisible until the new result is frozen.
- Use the successful historical annotation only for post-freeze biological equivalence review, never as a runtime Atlas or expected membership.
- Do not convert failed diagnostic runs into labels, thresholds, fixtures, or external references.
- Only the contract-bound standard controller may write release membership. Project-local scorers, subset writers, and remainder closers produce experimental evidence only.
- Apply every accepted semantic repair with the canonical cell-ID patch writer; positional row assignment is a failed diagnostic, regardless of apparent table-size agreement.
- Recompute completion from raw evidence and membership products; do not accept self-reported PASS fields.

The transferable success criterion is architectural and biological equivalence: cohort-first discovery, complete second-round lineage review, local handling of genuine mixtures, correct ovarian spatial boundaries, low unresolved burden, and traceable final membership. It is not reproduction of a historical cell count or per-observation label map.
