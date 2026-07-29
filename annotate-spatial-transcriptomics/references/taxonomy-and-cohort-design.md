# Literature taxonomy, computational cohorts and release labels

Keep three layers separate throughout annotation. They answer different questions and must never be copied into one another.

## 1. Literature reference taxonomy

A published atlas supplies a **candidate-lineage checklist**, not a list that the query must reproduce. Record the reference species, age/stage, tissue sampling, observation unit and dissociation/spatial technology before using its labels.

For sheep ovary, the 2025 developmental atlas (PMID 40641558; DOI 10.1016/j.isci.2025.112422) reported nine major types across prenatal, pre-pubertal, post-pubertal and adult samples: stromal, granulosa, oocyte, immune, epithelial, endothelial, theca, mesenchymal and smooth muscle. This is a useful audit checklist, but it is not a nine-label release requirement. Its Mesenchymal population was small and S100A4-associated, and developmental populations may be absent or inseparable in one adult spatial section.

GSE233801 (PMID 37964337) sampled adult Hu sheep within 12 hours of estrus, identified five somatic types and four granulosa subtypes, and retained reference limitations. It is best used for adult somatic and granulosa anchors, not as an exhaustive ovary taxonomy or a source of automatic Oocyte/Theca labels.

The 2024 Tibetan-sheep *Molecular Biology and Evolution* atlas (PMID 38552245; DOI 10.1093/molbev/msae058) provides an important independent counterbalance: 111,548 ovarian cells were summarized as stromal, granulosa, oocyte, endothelial, epithelial, immune and perivascular classes. Its omission of standalone Theca, Mesenchymal and Smooth muscle does not prove those lineages are absent; together with the developmental atlas, it shows that their release level depends on sampling and resolvable programs. A 2026 *Advanced Science* sheep–human reproductive atlas (DOI 10.1002/advs.202517633), 2025 *Science* human–mouse whole-ovary atlas (DOI 10.1126/science.adx0659), 2025 AJOG adult-ovary expert review (DOI 10.1016/j.ajog.2024.05.046) and a 2026 nine-species ovarian integration (PMID 41975518) strengthen conserved lineage/program evidence, but their cross-tissue/cross-species designs still do not impose a fixed adult-section taxonomy. The AJOG review additionally supports explicit sampling limitations and rejects single-marker or unlimited-subtype reasoning.

Targeted high-level sheep studies answer narrower questions. A 2026 *FASEB Journal* whole-ovary study strengthens macrophage–granulosa boundary and crosstalk evidence; purified cumulus, single-oocyte, atresia or maturation multiomics studies refine state modules and strict Oocyte validation. They cannot supply the whole-tissue broad-class tree. Read `profiles/sheep_ovary_evidence.md` for the weighted source table.

Cross-species studies may clarify boundaries but do not override the query. Human theca/stroma work (PMID 36599970) supports a progenitor-to-structural/perifollicular/androgenic continuum. Human morphologically guided spatial work (PMID 38578993) emphasizes oocyte, theca and granulosa programs and cortex/medulla variation while reporting only four major scRNA-seq types. These differences demonstrate why missing literature labels must be audited, not manufactured.

No study defines a fixed sheep-ovary broad list. The catalog therefore starts with common review boundaries but remains open to stage-, breed- and section-specific luteal, neural/glial or other programs. Endothelial, Pericyte/mural and Smooth muscle are independent competing broad identities; Lymphatic endothelial is an optional high-confidence child of Endothelial. The legacy `Vascular-associated` umbrella is an anatomical relationship, not a releasable identity. The catalog is an audit surface and permission model, not a required output list.

## 2. Computational cohorts and unresolved state

New projects use these boundaries:

| Boundary | Membership | Purpose |
|---|---|---|
| `whole_tissue_partition` | every analysis observation in exactly one initial cluster | Build second-round cohorts; record provisional biology only. |
| `initial_cluster_cohort` | exact members of one initial cluster | Perform the main SCT/PCA/SNN/Leiden annotation against the full catalog. |
| `local_mixed_subcluster` | one triggered second-round subcluster | Resolve two or more separable competing identities and its exact remainder. |
| `post_merge_unresolved` | biological members still unlabeled after all cohorts merge | Receive Atlas challenge/rescue and final typed closure. |

A cohort is a computational membership, not a biological category. It is never named from the first-pass provisional record, and cross-lineage or missing-broad returns do not automatically enter another cohort. `unresolved_biological` remains distinct from technical QC until final materialization.

## 3. Release broad classes

Release labels describe biology and require query-specific evidence. Use the following sheep-ovary vocabulary as a naming policy, not a quota.

### Default biological broad classes

- `Granulosa`
- `Stromal/mesenchymal`
- `Endothelial`
- `Pericyte/mural`
- `Smooth muscle`
- `Immune`
- `Epithelial/mesothelial`
- `Oocyte`, only after the strict context gate

### Evidence-dependent standalone broad classes

- `Theca`: reserve for a coherent steroidogenic/androgenic theca program with follicular outer-ring morphology. Do not publish `Theca/follicular wall` as a broad label; structural follicular wall may be stroma, smooth muscle, pericyte or interface.
- `Endothelial`: publish only when an endothelial junction backbone and an independent angiovascular support family are coherent. Branching morphology supports but cannot assign identity.
- `Pericyte/mural`: publish as an independent broad only when the RGS5/PDGFRB/CSPG4/NOTCH3/MCAM/RERGL-like backbone plus independent contractile support separates it from endothelial, mature smooth muscle and generic stroma.
- `Smooth muscle`: publish when a mature MYH11/CNN1/ACTG2/SMTN/LMOD1-centered backbone, independent contractile support, stable separation and coherent nonvascular follicular/hilar/structural tracks pass. This is not synonymous with ACTA2/TAGLN-positive stroma or vascular-wall mural cells.
- `Luteal`: requires bound stage/context, a coherent steroidogenic core, an independent corpus-luteum identity family and a corpus-luteum-like spatial structure; STAR/CYP11A1, PGR/APOD/SFRP4 or a solid component alone is insufficient. Granulosa- and Theca-lineage-of-origin programs are not hard contradictions because both contribute luteal cells. They remain competing origins that must be resolved with luteal-specific programs, DEG/pseudobulk and mass-like rather than thin perifollicular morphology.
- `Glial/Schwann-like`: requires a coherent glial program and nerve-track morphology.

`Mesenchymal progenitor-like`, `Neuroendocrine-like` and structural/perifollicular Theca remain exploratory records. They do not become broad labels without a future catalog revision based on independent sheep-ovary identity evidence.

The final label must be the least specific honest name. Use broad `Endothelial` when blood versus lymphatic separation is ambiguous but endothelial identity is secure; when Endothelial versus Pericyte/mural versus Smooth muscle cannot be resolved, retain `unresolved_biological` rather than inventing a generic vascular label. Use `Immune` when myeloid versus lymphoid support is too shallow and `Stromal/mesenchymal` when a standalone Mesenchymal class is unsupported.

### Non-biological retained states

The following are not biological broad classes and must not be included in the biological broad-class census or broad DEG/dotplot tree:

- `Anatomical interface` or a named resident-resident interface;
- `Low-information/QC holdout`;
- `Technical state`;
- `Pending review`.

Report them in a separate retained-state census and spatial layer. They remain mutually exclusive ledger outcomes and must not be hidden.

## 4. Release decision rules

Before a broad label is frozen from a second-round outcome:

1. Demonstrate at least two explicit independent positive marker families on the full-feature object. Validate the profile contains those families before applying the gate. Use absolute detection/prevalence and pseudobulk for broad presence; centered module scores and one-vs-rest DEG cannot reject a parent program merely because it is shared by several abundant clusters.
2. Quantify major anti-program leakage at observation level, not only cluster-average DEG.
3. Verify spatial morphology when coordinates exist.
4. Review stability across the selected cohort resolution and its nearest neighbors.
5. For a mixed subcluster, use the local candidate-component and exact-remainder route before closure.
6. Record a negative audit for every biologically plausible but unsupported standalone class; never lower its gate merely to match a paper.

If a fine label fails, roll it back to the supported broad parent and retain its ECM, contractile, hypoxic, stress, cycle, ambient or anatomical characteristics as tags. If a broad biological label fails but the population is a genuine interface or irreducible low-information state, retain that state explicitly rather than forcing the closest atlas label.

## 5. Regression lessons from the R-first ovary forward test

The forward test established reusable failure checks:

- A broad `Theca/follicular wall` bucket can absorb mature smooth muscle, generic ECM stroma, granulosa and endothelial cells. Reopen it with separate steroidogenic, contractile, stromal, granulosa and endothelial programs.
- A strong mature-contractile population with ring/track morphology can be hidden inside Theca or Stroma. The smooth-muscle audit is mandatory even when no initial cluster carries that name.
- A stromal cohort containing CDH5/PECAM1/CLDN5/PTPRB/ROBO4/MMRN2-positive branching tracks must return molecularly supported observations directly to broad `Endothelial`; ordinary non-lymphatic endothelium does not need a redundant Blood endothelial fine label.
- Do not create Mesenchymal or Pericyte merely because a reference lists them. A machine-readable negative audit is an acceptable result.
- Zona or other oocyte-adjacent RNA in granulosa does not establish Oocyte. Report cellbin/spot counts separately from inferred biological objects.
- A query cluster with vascular-adjacent markers but dominant granulosa lineage support may remain Granulosa with a spatial/state tag; top DEG alone must not switch its lineage.
- Resolution is selected to separate supported lineages, not to maximize cluster count. Several computational clusters may merge into one release class.
- An ambiguous blood/lymphatic split rolls back to broad Endothelial; an unresolved endothelial–mural–smooth-muscle identity remains unresolved_biological.

These are regression tests for reasoning, not sample-specific label maps.
