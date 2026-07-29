# Local mixed-subcluster splitting and exact remainder

## Scope and authority

This route is available only after second-round cohort scoring identifies one exact Leiden subcluster with either two independently supported specific identity cores, or one orthogonally supported specific identity embedded in a supported generic remainder. Whole-subcluster DEG is supporting evidence, not a trigger requirement: a minority Granulosa, Smooth muscle, Endothelial, Pericyte/mural or other lineage can have negative aggregate DEG after dilution by the majority parent. Such a program triggers one bounded candidate-local separability check. It is written back only if observation discriminators, candidate-local pseudobulk/DEG, spatial components or neighboring-resolution partitions show that the signal belongs to a separable member set. The route must never run as a whole-tissue or whole-cohort per-cellbin classifier.

Formal execution is contract-bound through `run_lineage_controller.py`, which may invoke `run_observation_lineage_scoring.R`, `derive_candidate_local_subsets.R` and `close_exact_remainders.py`. Project-local variants are experimental.

## Observation evidence

- Score the complete catalog within the exact source subcluster and nearest-resolution evidence.
- Scale each gene by the current query's nonzero 95th percentile.
- Aggregate a family with the two strongest available genes.
- Combine fixed-kNN direct/local evidence at `0.35/0.65` by default.
- One anti-marker or local anti signal lowers confidence; only direct multigene contradictions spanning independent anti-families are hard.
- Freeze scores before proposal construction and do not read provisional/historical labels.

## Independent candidate components

Generate every candidate's identity-core seed and local spatial/expression components independently. Generic ECM, a shared parent program or `ACTA2/TAGLN` cannot bridge distant identity cores. The aggregate winner cannot suppress another candidate.

One source subcluster may therefore propose multiple subsets. Resolve overlap using normalized evidence, pairwise discriminators and anti-programs. A clear winner receives the cell; unresolved overlap stays with a supported common parent or `unresolved_biological`. Catalog order is never a tie-breaker. Specific lineages are adjudicated before generic Stromal remainder. If simultaneous lineage signals remain coexpressed in the same observations and no expression/spatial partition is validated, do not split or invent a lineage; retain the supported broad parent and keep the alternative as `watch`.

Generic Stromal/ECM coherence is a possible final remainder parent, not an independent competing identity by itself. A specific candidate must demonstrate its own multi-family core plus orthogonal group evidence before it can displace that parent or trigger a split.

Do not use the aggregate subcluster contradiction fraction as a candidate-visibility filter. A true mixed subcluster is expected to be contradictory before its memberships are separated. Candidate visibility requires multi-family direct identity-core, stable and orthogonal group evidence; the contradiction ceiling is applied only after a candidate-local subset has formed, or when deciding whether the entire subcluster is pure enough to inherit one identity.

The default `0.70 support / 0.30 competitor margin / <=0.05 contradiction` thresholds validate an already formed group. They are not per-observation admission gates.

## Whole-subcluster fast path

Do not call the local splitter when a second-round subcluster has a coherent multigene identity, DEG/pseudobulk support, reasonable morphology, neighboring-resolution stability and no independent competitor. Return it wholesale; sparse noncontradictory members inherit broad identity.

## Exact local remainder

After accepting disjoint candidate subsets:

1. remove their exact IDs from the original source subcluster;
2. keep observation scores immutable;
3. recompute candidate prevalence, competition, spatial components and group evidence on the exact remainder;
4. return a reliable common parent when no embedded competitor remains;
5. allow one final candidate extraction if a new coherent component appears;
6. retain genuinely ambiguous members as `unresolved_biological`.

An independently certified cross-lineage identity remains a remainder blocker even when the unsplit aggregate remainder exceeds the ordinary 5% contradiction ceiling. Do not reapply that ceiling to hide the blocker and return the whole remainder to one parent. State-only programs are stored in the state column and are not broad-parent blockers by themselves.

The final extraction may only revalidate residual membership from candidate-local spatial components already generated inside the source mixed subcluster. It cannot rerun an unrestricted per-observation candidate scan. The Oocyte canonical-component route is candidate-specific: it requires a canonical group challenger, an identity-core-only spatial component, two supported families, positive program/direct enrichment and local spatial support. Somatic anti signal remains reported but cannot veto that route by itself; no other candidate inherits the exception.

Failure to enter a subset is never, by itself, a QC reason. Local remainder artifacts cannot read or modify full-object membership outside their source subcluster.

## Required audit

The local artifact records its trigger, source cohort/subcluster and membership hash, each candidate proposal, overlap decision, accepted disjoint subsets, remainder rounds and final coverage. The union of accepted subsets plus parent/unresolved remainder must exactly equal the source subcluster.

Before overlaying any accepted repair on a frozen membership, run `apply_cell_id_membership_patch.py` with the exact columns authorized for change. The base and proposal must contain unique, nonempty `cell_id`; proposal IDs must be a subset of the base; output order must equal base order; and every non-proposal value in an updated column must remain unchanged. Never assign a marker matrix, score table or proposal vector to membership rows by position, even when both tables have the same row count.
