# Literature taxonomy, computational cohorts and release labels

Keep three layers separate throughout annotation. They answer different questions and must never be copied into one another.

## 1. Literature reference taxonomy

A published atlas supplies a **candidate-lineage checklist**, not a list that the query must reproduce. Record the reference species, age/stage, tissue sampling, observation unit and dissociation/spatial technology before using its labels.

For sheep ovary, the 2025 developmental atlas (PMID 40641558; DOI 10.1016/j.isci.2025.112422) reported nine major types across prenatal, pre-pubertal, post-pubertal and adult samples: stromal, granulosa, oocyte, immune, epithelial, endothelial, theca, mesenchymal and smooth muscle. This is a useful audit checklist, but it is not a nine-label release requirement. Its Mesenchymal population was small and S100A4-associated, and developmental populations may be absent or inseparable in one adult spatial section.

GSE233801 (PMID 37964337) sampled adult Hu sheep within 12 hours of estrus, identified five somatic types and four granulosa subtypes, and retained reference limitations. It is best used for adult somatic and granulosa anchors, not as an exhaustive ovary taxonomy or a source of automatic Oocyte/Theca labels.

The 2024 Tibetan-sheep *Molecular Biology and Evolution* atlas (PMID 38552245; DOI 10.1093/molbev/msae058) provides an important independent counterbalance: 111,548 ovarian cells were summarized as stromal, granulosa, oocyte, endothelial, epithelial, immune and perivascular classes. Its omission of standalone Theca, Mesenchymal and Smooth muscle does not prove those lineages are absent; together with the developmental atlas, it shows that their release level depends on sampling and resolvable programs. A 2026 *Advanced Science* sheep–human reproductive atlas (DOI 10.1002/advs.202517633), 2025 *Science* human–mouse whole-ovary atlas (DOI 10.1126/science.adx0659), 2025 AJOG adult-ovary expert review (DOI 10.1016/j.ajog.2024.05.046) and a 2026 nine-species ovarian integration (PMID 41975518) strengthen conserved lineage/program evidence, but their cross-tissue/cross-species designs still do not impose a fixed adult-section taxonomy. The AJOG review additionally supports explicit sampling limitations and rejects single-marker or unlimited-subtype reasoning.

Targeted high-level sheep studies answer narrower questions. A 2026 *FASEB Journal* whole-ovary study strengthens macrophage–granulosa boundary and crosstalk evidence; purified cumulus, single-oocyte, atresia or maturation multiomics studies refine state modules and strict Oocyte validation. They cannot supply the whole-tissue broad-class tree. Read `profiles/sheep_ovary_evidence.md` for the weighted source table.

Cross-species studies may clarify boundaries but do not override the query. Human theca/stroma work (PMID 36599970) supports a progenitor-to-structural/perifollicular/androgenic continuum. Human morphologically guided spatial work (PMID 38578993) emphasizes oocyte, theca and granulosa programs and cortex/medulla variation while reporting only four major scRNA-seq types. These differences demonstrate why missing literature labels must be audited, not manufactured.

The cross-study adult sheep candidate backbone is therefore `Granulosa`, `Stromal/mesenchymal`, `Vascular/endothelial`, `Immune`, `Epithelial/mesothelial` and strict `Oocyte`. `profiles/sheep_ovary_candidate_lineage_catalog.json` expands the mandatory boundary audit to evidence-, stage- and anatomy-dependent alternatives including steroidogenic Theca/luteal, mature Smooth muscle, Pericyte/mural, Mesenchymal progenitor-like, lymphatic endothelial and neural/glial/neuroendocrine programs. It is non-exhaustive and is still an audit surface, not a requirement that every label appear.

## 2. Computational cohorts and QC state

New projects do not create persistent biological pools. They use three auditable boundaries:

| Boundary | Membership | Purpose |
|---|---|---|
| `broad_class_recluster` | all observations assigned to one supported initial broad class | Test shallow subtypes and purity within the parent class. |
| `targeted_recluster` | only observations needed to answer one local mixture, contamination or context-gated question | Resolve that question once; it cannot become a long-lived catch-all. |
| `qc_holdout` | all final low-information or irreducibly mixed observations after broad/targeted decisions | Exact query boundary for the terminal calibrated Atlas review; it is not reclustered. |

A cohort is a computational membership, not a biological category. Every subcluster returns directly to a supported broad/fine label, crosses directly to another lineage, enters one targeted question, or remains QC/technical. A cross-lineage return never creates an intermediate target pool and does not automatically trigger another target-lineage reclustering.

## 3. Release broad classes

Release labels describe biology and require query-specific evidence. Use the following sheep-ovary vocabulary as a naming policy, not a quota.

### Default biological broad classes

- `Granulosa`
- `Stromal/mesenchymal`
- `Vascular/endothelial`
- `Immune`
- `Epithelial/mesothelial`
- `Oocyte`, only after the strict context gate

### Evidence-dependent standalone broad classes

- `Theca`: reserve for a coherent steroidogenic/androgenic theca program with follicular outer-ring morphology. Do not publish `Theca/follicular wall` as a broad label; structural follicular wall may be stroma, smooth muscle, pericyte or interface.
- `Smooth muscle`: publish when a mature contractile backbone, stable separation and coherent vessel-wall/hilar/structural tracks pass. This is not synonymous with ACTA2/TAGLN-positive stroma.
- `Pericyte/mural`: publish when the RGS5/PDGFRB/CSPG4/NOTCH3/MCAM backbone separates it from endothelial, smooth muscle and generic stroma. Otherwise retain a mural state tag under `Stromal/mesenchymal`.
- `Mesenchymal progenitor-like`: publish separately from `Stromal/mesenchymal` only when a stable S100A4/progenitor-like program and morphology pass the profile gate. Absence after a documented negative audit is valid.
- `Luteal steroidogenic`: requires stage/context plus a coherent corpus-luteum-like spatial structure; STAR/CYP11A1 alone is insufficient.
- `Neural/Schwann`: requires a coherent glial program and nerve-track morphology.
- `Neuroendocrine`: requires a secretory-neuroendocrine core (`CHGA/CHGB/SYP/SCG/INSM1`), neuronal support and coherent rare-focus morphology. `DLG2/RBFOX1/TENM3` alone remains a state tag under the resident class.

The final label must be the least specific honest name. For example, use `Vascular/endothelial` when blood versus lymphatic separation is ambiguous; use `Immune` when myeloid versus lymphoid support is too shallow; use `Stromal/mesenchymal` when a standalone Mesenchymal class is unsupported.

### Non-biological retained states

The following are not biological broad classes and must not be included in the biological broad-class census or broad DEG/dotplot tree:

- `Anatomical interface` or a named resident-resident interface;
- `Low-information/QC holdout`;
- `Technical state`;
- `Pending review`.

Report them in a separate retained-state census and spatial layer. They remain mutually exclusive ledger outcomes and must not be hidden.

## 4. Release decision rules

Before a broad label is frozen:

1. Demonstrate at least two independent positive marker families on the full-feature object.
2. Quantify major anti-program leakage at observation level, not only cluster-average DEG.
3. Verify spatial morphology when coordinates exist.
4. Review stability across adjacent whole-tissue or cohort resolutions.
5. For a large heterogeneous label, use its broad-class or targeted query-only cohort before closure.
6. Record a negative audit for every biologically plausible but unsupported standalone class; never lower its gate merely to match a paper.

If a fine label fails, roll it back to the supported broad parent and retain its ECM, contractile, hypoxic, stress, cycle, ambient or anatomical characteristics as tags. If a broad biological label fails but the population is a genuine interface or irreducible low-information state, retain that state explicitly rather than forcing the closest atlas label.

## 5. Regression lessons from the R-first ovary forward test

The forward test established reusable failure checks:

- A broad `Theca/follicular wall` bucket can absorb mature smooth muscle, generic ECM stroma, granulosa and endothelial cells. Reopen it with separate steroidogenic, contractile, stromal, granulosa and endothelial programs.
- A strong mature-contractile population with ring/track morphology can be hidden inside Theca or Stroma. The smooth-muscle audit is mandatory even when no initial cluster carries that name.
- A stromal cohort containing CDH5/PECAM1/CLDN5/PTPRB/ROBO4/MMRN2-positive branching tracks must return those observations directly to `Vascular/endothelial` after evidence review.
- Do not create Mesenchymal or Pericyte merely because a reference lists them. A machine-readable negative audit is an acceptable result.
- Zona or other oocyte-adjacent RNA in granulosa does not establish Oocyte. Report cellbin/spot counts separately from inferred biological objects.
- A query cluster with vascular-adjacent markers but dominant granulosa lineage support may remain Granulosa with a spatial/state tag; top DEG alone must not switch its lineage.
- Resolution is selected to separate supported lineages, not to maximize cluster count. Several computational clusters may merge into one release class.
- An ambiguous blood/lymphatic or other fine split must roll back to the broader honest label.

These are regression tests for reasoning, not sample-specific label maps.
