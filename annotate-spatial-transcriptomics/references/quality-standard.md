# Biological quality standard

## Primary endpoint

For spatial spots/cellbins, optimize broad-class correctness, spatial coherence and complete uncertainty handling. Define fine subtypes only when full-feature marker/anti-marker evidence, stable cohort-specific reclustering and morphology agree. A biologically honest broad-only label is a successful endpoint; a weak fine label is not.

Use a shallow-tree default. The number of graph clusters, selected resolutions or DEG tables does not determine the number of biological subtypes. Penalize one-cluster-one-name annotation, literature-name completion and state-only labels presented as lineages. Preserve useful ECM, contractile, stress, low-RNA, ambient and anatomical differences as tags even when several clusters merge into one biological class.

Keep the literature checklist, computational cohorts and release taxonomy distinct. Penalize one-cohort-one-label annotation, copying cohort/QC identifiers into the final tree, and adding a literature class after its query-specific gate fails. The report must separate biological broad classes from retained anatomical-interface, QC, technical and pending states.

## Method-independent acceptance

Do not optimize a BANKSY, Seurat or Scanpy result to reproduce another method's cell-level labels. If a previous annotation exists, keep it blinded during fitting, resolution selection and writeback. Use it only after freezing the new result as a diagnostic comparison.

Comparable annotation quality requires:

1. Exact `full_object`, `analysis_set`, initial-QC and post-clustering-holdout accounting.
2. Freeze a label-blind, all-candidate broad evidence matrix before initial labels: positive DEG, anti-DEG, winner/runner-up, margin, marker-family coverage, contradictions and technical flags. Paper markers cannot narrow the candidate set after seeing the result.
3. Per-class cell/bin-level lineage-core, support and anti-program audits on the full-feature object.
4. Marker, DEG, source/QC composition and spatial morphology review for every large label and every broad/targeted cohort.
5. Continuous full-catalog signal coverage at the whole-tissue, broad-cohort and targeted-cohort levels. Audit the selected resolution and next two higher available resolutions; retain weak coherent signals across boundaries instead of treating failure of the naming threshold as absence.
6. Open-world discovery plus strict context gates for ambient-prone identities; rarity alone is not a class or route. Report spatial objects separately from observation counts.
7. One reclustering cohort for every supported initial broad class, with an underpowered skip allowed only when recorded.
8. Tiered low-priority RCTD evidence: canonical high plus independent evidence may support fine, moderate returns broad-only, and low enters the frozen QC holdout rather than Atlas.
9. After terminal QC is frozen, run one calibrated broad-only Atlas mapping over the complete analysis set. `high`/`moderate_only` plus current-query and independent spatial/internal support may directly rescue only QC; `low_reject`, OOD and ontology-conflicted rows do not write back. Defined cells use the mapping only for concordance challenge. Accepted QC returns cannot seed fine discovery. Default held-out precision targets remain 0.90/0.95.
10. Every material defined-label disagreement and coherent OOD group reopens the complete cluster/cohort exactly once. Atlas alternatives require independent query full-feature evidence to supersede; mixed evidence downgrades or becomes unknown/QC.
11. A small/local retained interface, or a documented irreducible QC/technical remainder after every applicable route. Large/diffuse retention automatically reopens.
12. Broad and subtype evidence assets are separate. The subtype tree may be shallow when the data support only broad identities.
13. A navigable report whose annotation tree, node highlights, DEG, marker dotplots, spatial gene maps, source ancestry and detailed workflow all resolve to audited artifacts.
14. A subtype parsimony audit showing that every fine label adds a reproducible functional or lineage program beyond its parent; unsupported splits are merged and retained as state tags.
15. A taxonomy/cohort audit showing that every release label passed its own gate, every unsupported plausible literature class has a negative audit, and no routing/technical state enters biological DEG or marker discovery.
16. Content-schema validation for every prelabel freeze, cohort outcome, direct return and broad/fine support record. Artifact existence alone, empty evidence and status-only placeholders fail.
17. Exact observation-level partition closure through global Atlas query coverage, frozen-QC accepted/rejected writeback, defined-label review closure and the final annotation.
18. A post-Atlas query-derived broad-class completeness audit: present classes pass full-membership expression/spatial/purity review, while every zero-census default tissue lineage has a multichannel negative audit independent of Atlas.
19. Every expression object used as query evidence is project-local and hash-bound to its raw counts, analysis set and parent artifact; cross-project derived expression is reference-only and explicitly registered.

## Large-label purity trigger

A convincing cluster-level DEG does not validate every observation in that cluster. Reopen a large direct label when a cell-level lineage backbone is sparse, anti-programs dominate, spatial distribution is too broad for the proposed identity, or a small coherent lineage appears embedded in a resident class. Freeze only the coherent cells; return the remainder directly to a supported broad lineage, a targeted cohort, or QC.

The purity trigger and rare-lineage recall use the same persistent signal ledger. Every cohort remains open to cross-lineage evidence. An embedded lineage that was subthreshold at the whole-tissue pass is a tracked hypothesis, not a rejected class; a coherent recurrence in a subcluster must reopen broad-class construction.

Use observation-level core/support/anti hits and local spatial connected components inside large labels. A coherent program occupying only 1–5% of a parent label is not refuted by a weak parent-cluster average. It enters `watch/candidate` and receives one bounded targeted recall. More graph clusters without improved lineage separation and anti-program clearance do not constitute rescue.

## Unresolved-fraction trigger

Do not impose a universal annotation-rate quota. Treat any large or spatially diffuse interface/QC/pending fraction as a mandatory evidence audit. After broad/targeted cohorts are terminal, freeze residual QC once and run the calibrated all-cell broad concordance pass without QC reclustering. Report the complete mapping denominator, frozen-QC accepted/rejected outcomes, defined-label disagreements and OOD groups separately.

## Forward-test target

A fresh Agent receives raw inputs, biological context and this Skill, but not the intended clustering or final labels. It must autonomously discover inputs, select/adapt resolutions, submit and repair jobs, reopen overbroad labels, route uncertainty, write immutable state and build the complete report with only genuinely blocking user questions. Evaluate evidence completeness and biological safety first; historical label agreement is a blind secondary diagnostic.

## Post-completion main-Agent approval

For every sample, the main conversation Agent performs one biological quality approval only after all broad/targeted cohorts, direct returns, all-cell Atlas concordance/OOD review, residual-QC writeback, final label materialization and the completion gate have finished. The completion gate proves workflow closure; the approval does not repeat it. Review broad-label reasonableness, marker/anti-marker plus spatial support, context-gated/confounded lineage safety, and whether the result reaches the evidence quality of the validated deidentified reference strategy. Exact label/count agreement and equally deep subtypes are not required. A result may pass with documented concerns; a material biological error returns the same sample worker to targeted iteration. User confirmation and final assets remain blocked until this approval passes.
