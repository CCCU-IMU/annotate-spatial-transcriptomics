# v2.2 testing contract

## Syntax and schema

Run Python compilation, R parsing, JSON parsing, schema validation, repository validation and install verification. Contract/controller/script/profile hashes must match, and formal output must come from the canonical staged controller.

## Architecture tests

- The first phase cannot create formal broad, fine, QC or release membership.
- `whole_tissue_cohort_plan.tsv` exactly partitions `analysis_set`, with one second-round cohort per initial cluster.
- Provisional fields and historical labels cannot enter second-round scoring or resolution selection.
- A cohort must use project-local non-SCT raw counts; SCT corrected counts and cross-project expression fail.
- Every selected second-round subcluster has complete open-catalog evidence.
- A lineage absent in the first pass can become `missing_broad_reconstruction` in the second round.
- Observation-level splitting fails unless an exact mixed second-round subcluster trigger exists.
- Local subset nonmembers remain parent/unresolved rather than becoming QC automatically.
- Local remainder cannot read membership outside its exact source subcluster.
- Formal broad membership cannot exist before all cohorts merge into an exact disjoint cover.
- Fine proposals cannot materialize before broad freeze or outside the matching parent.
- Atlas, completeness and residual-QC completion checks fail before broad merge.
- Atlas may rescue only unlabeled broad membership and cannot overwrite existing broad or create fine.
- `failed_diagnostic` artifacts cannot enter runtime inputs, reference truth or expected fixture labels.
- Identical input, contract and seed reproduce partition, membership and semantic hash.

## Biological integration fixtures

Build fixtures directly from allowed raw input and never from A08/A09 membership or the designated blind-regression sample's historical labels.

1. **Granulosa-rich:** a pure second-round subcluster returns wholesale; sparse tails inherit broad identity; no diffuse stromal expansion.
2. **Stromal/Smooth muscle/Endothelial/Pericyte mixed:** local splitting independently produces mature nonvascular Smooth muscle, Endothelial and Pericyte/mural subsets plus a Stromal remainder; an inseparable complete endothelial–mural mixed cellbin remains unresolved_biological.
3. **Rare Oocyte:** a coherent canonical cluster returns all noncontradictory members; zona/ambient cells and neighboring granulosa stay outside.
4. **Low-fraction Epithelial:** a 3%-5% coherent program is reconstructed without expanding the enclosing stromal population.
5. **Unmodeled program:** a nonstress/noncycle/non-ECM program repeated across neighboring resolutions and space is recorded without automatic naming.

Also test that one anti-marker is soft, independent multigene direct anti-programs can be hard, aggregate winners cannot suppress independent candidates, zero-signal candidates cannot beat real supported candidates, and an empty fine census requires a complete parent-by-candidate audit.

## Robustness tests

Technical repeatability requires exact equality of partition, membership and semantic hash.

Parameter perturbation compares selected and neighboring reasonable resolutions, PCs +/-5 and k +/-20% while holding other contract fields constant. Require:

- every >=1% broad class remains present;
- no >=1% broad class appears or disappears from perturbation alone;
- broad census Spearman >=0.95;
- each major broad spatial Dice >=0.80;
- context-gated candidates receive the same bound exogenous context in second-round and local-split stages, while the context artifact has no label or membership authority;
- granulosa/theca lineage-of-origin programs do not hard-veto a context-compatible Luteal candidate, while generic steroidogenesis without an independent corpus-luteum identity family cannot release Luteal;
- major fine labels retain parent and core identity program.

Rare lineages are reviewed by multigene identity and anatomy rather than prevalence thresholds.

## Independent designated-sample blind regression

After all fixtures pass, create a new project from the fixed SCT+BANKSY input. Do not read v2.1.0 labels, repair membership, A08/A09 outputs or same-batch Atlas before freeze. Run the complete pipeline twice for determinism and run perturbation variants for biological robustness. Only after freeze load v2.1.0 for external comparison.

Acceptance requires biologically equivalent major broad types and follicular/vascular/surface structures, no diffuse Granulosa, correct Endothelial/Pericyte/mural/Smooth-muscle separation, no whole-cluster Epithelial expansion or complete epithelial loss, canonical Oocyte recall, missing-broad reconstruction, multichannel support for every broad and final QC below both 10% and 50,000. Exact cellbin equality is not required.

## Forward test

A fresh Agent receives raw inputs, context and the Skill but not selected resolutions, labels or historical answers. It must autonomously construct cohorts, annotate second-round subclusters, resolve local mixtures, freeze broad membership, run Atlas/completeness, build the report and ask only genuinely blocking questions.
