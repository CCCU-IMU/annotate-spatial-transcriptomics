# Biological quality standard

## Primary endpoint

For spatial spots/cellbins, optimize broad-class correctness, spatial coherence and complete uncertainty handling. Define fine subtypes only when full-feature marker/anti-marker evidence, stable pool-specific reclustering and morphology agree. A biologically honest broad-only label is a successful endpoint; a weak fine label is not.

Use a shallow-tree default. The number of graph clusters, selected resolutions or DEG tables does not determine the number of biological subtypes. Penalize one-cluster-one-name annotation, literature-name completion and state-only labels presented as lineages. Preserve useful ECM, contractile, stress, low-RNA, ambient and anatomical differences as tags even when several clusters merge into one biological class.

Keep the literature checklist, analysis pools and release taxonomy distinct. Penalize one-pool-one-label annotation, copying `_review`/`_candidate`/`_unresolved`/`_holdout` names into the final tree, and adding a literature class after its query-specific gate fails. The report must separate biological broad classes from retained anatomical-interface, QC, technical and pending states.

## Method-independent acceptance

Do not optimize a BANKSY, Seurat or Scanpy result to reproduce another method's cell-level labels. If a previous annotation exists, keep it blinded during fitting, resolution selection and writeback. Use it only after freezing the new result as a diagnostic comparison.

Comparable annotation quality requires:

1. Exact `full_object`, `analysis_set`, initial-QC and post-clustering-holdout accounting.
2. Per-class cell/bin-level lineage-core, support and anti-program audits on the full-feature object.
3. Marker, DEG, anchor distance, source/QC composition and spatial morphology review for every large label and every large unresolved pool.
4. Open-world discovery plus strict context gates for ambient-prone identities; rarity alone is not a class or route. Report spatial objects separately from observation counts.
5. Route-A anchor reclustering for large usable uncertainty before supervised assistance.
6. Tiered low-priority RCTD evidence: extreme plus independent evidence may support fine, high returns broad-only, and medium/low enters the frozen QC holdout rather than Atlas.
7. Atlas/internal-anchor rescue only for the residual cells left in QC after the complete QC-holdout anchor reclustering and depth-matched calibration; calibrated moderate-or-high Atlas calls may return broad-only, while lower calls remain rejects. Already defined broad/fine cells and ordinary biological pools are forbidden Atlas queries. Accepted cells cannot seed fine discovery. Default held-out target-precision tiers are moderate 0.90 and high 0.95; they are not raw per-cell score cutoffs.
8. A small/local retained interface, or a documented irreducible QC/technical remainder after every applicable route. Large/diffuse retention automatically reopens.
9. Broad and subtype evidence assets are separate. The subtype tree may be shallow when the data support only broad identities.
10. A navigable report whose annotation tree, node highlights, DEG, marker dotplots, spatial gene maps, source ancestry and detailed workflow all resolve to audited artifacts.
11. A subtype parsimony audit showing that every fine label adds a reproducible functional or lineage program beyond its parent; unsupported splits are merged and retained as state tags.
12. A taxonomy/pool audit showing that every release label passed its own gate, every unsupported plausible literature class has a negative audit, and no routing/technical state enters biological DEG or marker discovery.

## Large-label purity trigger

A convincing cluster-level DEG does not validate every observation in that cluster. Reopen a large direct label when a cell-level lineage backbone is rare, anti-programs dominate, spatial distribution is too broad for the proposed identity, or a small coherent lineage appears embedded in a resident class. Freeze only the coherent cells; return the remainder to the appropriate broad pool.

## Unresolved-fraction trigger

Do not impose a universal annotation-rate quota. Instead, treat any large or spatially diffuse interface/QC/pending fraction as a mandatory route audit. It may close only after the complete frozen membership has undergone the applicable unsupervised anchor, calibrated reference/atlas and technical-retention routes. Report the denominator and every outcome.

## Forward-test target

A fresh Agent receives raw inputs, biological context and this Skill, but not the intended clustering or final labels. It must autonomously discover inputs, select/adapt resolutions, submit and repair jobs, reopen overbroad labels, route uncertainty, write immutable state and build the complete report with only genuinely blocking user questions. Evaluate evidence completeness and biological safety first; historical label agreement is a blind secondary diagnostic.

## Post-completion main-Agent approval

For every sample, the main conversation Agent performs one biological quality approval only after pool reclustering, all applicable multi-route rescue, final label materialization and the completion gate have finished. The completion gate proves workflow closure; the approval does not repeat it. Review broad-label reasonableness, marker/anti-marker plus spatial support, rare or confounded lineage safety, and whether the result reaches the evidence quality of the validated deidentified sheep-ovary R-first reference. Exact label/count agreement and equally deep subtypes are not required. A result may pass with documented concerns; a material biological error returns the same sample worker to targeted iteration. User confirmation and final assets remain blocked until this approval passes.
