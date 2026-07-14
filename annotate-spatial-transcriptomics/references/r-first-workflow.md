# Seurat/R-first annotation workflow

Use this path whenever a readable full-feature Seurat RDS is available or the user selects R as the primary framework. The workflow remains method-independent biologically, but Seurat is the computational backbone.

## 1. Freeze and inspect the R object

Run `scripts/inspect_r_object.R`, `scripts/check_r_runtime.R` and `scripts/audit_feature_scope.R`. Record object hash, Seurat version, assays/layers, raw-count assay, normalized assay, reductions, metadata columns, coordinates, observation IDs and analysis-set membership. Use RNA/Spatial raw counts for absolute detection and anti-marker gates. SCT residuals may drive PCA/clustering but cannot establish marker absence.

Prefer the full-feature Seurat RDS over a reduced/HVG object. If both exist, freeze both roles explicitly: `clustering_object` and `validation_object`. A full feature list does not prove that `Spatial@data` is normalized. Before any Wilcoxon evidence, verify `data != counts`; when it fails, create a separate manifest-bound LogNormalize validation object with `scripts/prepare_seurat_full_feature_validation.R` and leave the SCT clustering object unchanged.

If a matched single-cell Seurat/AnnData object is available, freeze it as a separate `matched_reference_object`; never merge its cells into the query ledger. Read `matched-single-cell-reference.md`, validate its label crosswalk, and use it as the preferred external reference channel after current-query anchors.

## 2. Reuse computation without leaking labels

Existing whole-tissue or pool reclustering may be reused when all conditions pass:

1. Its source-object hash or exact observation membership is recoverable.
2. It used the current analysis set or has an explicit exclusion ledger.
3. Candidate resolutions, cluster memberships, DEG, UMAP and spatial outputs are complete.
4. Validation can be repeated on a verified full-feature normalized layer; a raw converted `Spatial@data == counts` object must first use the separate validation-object route.
5. Historical annotation columns are hidden until new decisions are frozen.

Register reused artifacts as `validated_reuse`, preserve their original path/hash and create a new decision ledger. Never copy their labels or confidence fields.

## 3. Whole-tissue broad pass

If reusable clustering is unavailable, run SCTransform with `return.only.var.genes=FALSE`, PCA and a coarse Seurat resolution grid appropriate to the current observation count and signal. Do not hard-code an example grid. Shortlist the lowest resolutions that preserve granulosa/follicular somatic, steroidogenic theca, stromal/mesenchymal, mural/smooth muscle, endothelial, immune, epithelial and strict rare candidates.

Generate one-vs-rest DEG, canonical and data-specific marker summaries, UMAP, whole-section spatial maps and per-cluster highlights for every shortlisted resolution. Select the lowest-complexity candidate that preserves supported broad lineages without state-only fragmentation.

## 4. Broad anchors and canonical parent pools

Freeze only full-gene, cell-level coherent anchors. Follow `taxonomy-and-pool-design.md`: a literature checklist, an analysis pool and a release label are different objects. Create the smallest set of axis-covering review pools needed by the actual uncertainty, rather than one graph-cluster pool or one pool per desired cell type:

- `follicular_somatic_review`: granulosa, steroidogenic theca and structural follicular-stromal alternatives;
- `stromal_mesenchymal_mural_review`: generic stroma, progenitor-like mesenchyme, mature smooth muscle, pericyte/mural and structural follicular wall;
- `vascular_endothelial_mural_review`: blood/lymphatic endothelium, mural, smooth-muscle and vascular-adjacent stromal alternatives;
- `immune_review` and `epithelial_mesothelial_review` when those lineages remain unresolved;
- `strict_oocyte_candidate`, local `anatomical_interface_review` and `postcluster_qc_holdout` only when their route definitions apply.

Every pool has immutable membership and explicit competing hypotheses. Its name cannot be copied into a final label. Every large direct label must pass cell-level purity and anti-program review. Split heterogeneous source clusters at observation level when necessary.

## 5. R pool controller

Use `scripts/run_seurat_pool_recluster.R` with frozen query/anchor memberships. Fit normalization/PCA jointly when anchors are used, but construct graph, clusters, UMAP, DEG and outcome counts from query observations only. Test a pool-specific resolution grid and select it independently from the whole-tissue resolution.

Apply the mandatory order:

1. Large usable uncertainty: balanced-anchor Seurat reclustering.
2. Local mixed interface: targeted reclustering, then calibrated RCTD only if applicable.
3. Medium/low RCTD or low-information rejects: full QC-pool anchor reclustering, then calibrated atlas/internal-anchor/marker/observed-density spatial consensus.
4. Rare identity: full-feature positive/negative modules, spatial objects and strict-candidate reclustering.
5. Irreducible technical or biological uncertainty: explicit retained state.

RCTD is low-priority assistance. Atlas moderate-or-higher may rescue broad-only after independent evidence; neither route creates fine anchors.

When a matched single-cell reference exists, use its harmonized broad labels before a generic atlas for the external-reference channel, provided stage, donor composition, raw-count availability and source marker programs pass audit. The matched prediction still has a broad-only default ceiling and cannot bypass strict Oocyte, Theca, mural/smooth-muscle or rare-lineage gates.

## 6. Ovary stromal-lineage decomposition

Do not close a single `Stromal/perivascular` super-class before testing:

- generic ovarian stromal/fibroblast;
- mesenchymal-progenitor-like;
- mature smooth muscle;
- pericyte/mural;
- blood and lymphatic endothelial;
- structural/follicular-wall theca versus androgenic steroidogenic theca.

Use the tissue profile's program families and anti-programs. A mature Smooth muscle call requires a coherent mature-contractile family such as `MYH11/CNN1/ACTG2` plus `TAGLN/ACTA2/MYL9/DES`, stable pool separation and vessel/follicular-wall morphology. Pericyte requires `RGS5/PDGFRB/CSPG4/NOTCH3/MCAM` with vascular adjacency. ECM-rich `ACTA2/TAGLN` cells without mature contractile backbone remain myofibroblast/structural stroma or follicular wall.

`Mesenchymal` is not automatically separate from ovarian stroma. Define a standalone `Mesenchymal progenitor-like` broad class only when a stable population has coherent `S100A4/PDGFRA/TCF21/NR2F2/CD34/PI16/COL14A1/DPT/CFD`-like progenitor/fibroblast support, lacks mature contractile, endothelial, immune, epithelial, granulosa and steroidogenic programs, and has reproducible morphology. Otherwise retain `Stromal/mesenchymal` as the parent label and record that no independent mesenchymal class was supported.

Release `Theca` only for coherent steroidogenic/androgenic theca. Never use `Theca/follicular wall` as a broad release label: return structural collagen/contractile wall cells to stromal, smooth-muscle, mural or interface review as supported. A negative Mesenchymal or Pericyte audit is valid; do not complete the literature taxonomy by lowering gates.

## 7. Shallow subtype policy

Broad classes are the default endpoint for cellbin spatial data. Do not translate every Seurat subcluster into a named biological subtype. Merge subclusters when their separation is driven mainly by depth, LOC/ribosomal/mitochondrial expression, ECM amount, stress, cell cycle, generic contractility or spatial position without an independent functional program.

For adult sheep ovary, permit only a small evidence-supported functional vocabulary by default:

- granulosa: early/undifferentiated, mural/estrogenic, cumulus/oocyte-adjacent, or luteinizing/steroidogenic only when the corresponding multi-gene program and follicular morphology are present;
- follicular somatic: androgenic steroidogenic `Theca` may be retained; structural/follicular-wall stroma remains under its supported stromal/mural parent or as an interface state;
- vascular/mural: blood endothelial, lymphatic endothelial, pericyte and mature smooth muscle only when each backbone passes;
- immune: broad myeloid, lymphoid or plasma categories only when depth supports them;
- stromal/mesenchymal: keep most variation as state/spatial tags unless a stable lineage-level separation passes the profile gate.

The sheep GSE233801 study supports a limited number of granulosa functional subtypes, not arbitrary cluster-by-cluster naming. Developmental sheep atlases contain stage-specific populations that must not be imported into an adult sample without query evidence. A shallower tree with stronger parent labels is preferable to a detailed but weak tree.

For sheep ovary, first run `scripts/resolve_workflow_profile.py`. A readable full-feature Seurat RDS makes this R-first workflow the default. A confirmed StereoPy `cellbin_PPed` conversion additionally activates the frozen whole-tissue preprocessing profile in `seurat-cellbin-preprocessing.md`; this specialized profile overrides any more general SCT examples in this document. Resolution and labels remain adaptive.

The 2025 Science human/mouse ovary atlas, 2026 Advanced Science sheep-human reproductive atlas and 2025 AJOG expert review support a stable shallow checklist and stronger negative gates. They do not require every class to appear. In particular, a neural-looking single gene cannot create Neural/Schwann or Neuroendocrine, and graph-derived granulosa states do not automatically become follicle-stage labels. Read `profiles/sheep_ovary_literature_2025_2026.md`.

## 8. Baseline-independent quality gate

When improvement over an earlier R release is requested, predeclare a comparison table before unblinding it. At minimum require:

- granulosa: no loss of multi-gene lineage enrichment or follicular-wall coherence;
- oocyte: equal or better non-ZP identity/maternal-core purity and lower resident-program leakage;
- theca: steroidogenic core is not diluted by ECM-only follicular wall;
- endothelial/immune/epithelial: no severe collapse in supported rare-lineage coverage;
- stromal: generic stroma, mesenchymal, smooth muscle, pericyte and endothelial alternatives were explicitly audited;
- uncertainty: every observation is accounted for and large rescued populations expose their direct versus transferred evidence;
- report: broad and subtype assets, spatial tree navigation, source ancestry and Chinese detailed workflow pass the release audit.
- subtype restraint: the new tree contains no cluster-renaming-only labels; unsupported historical subtypes are merged while their ECM/contractile/stress/spatial characteristics remain available as state tags.

Unblind the baseline only after the new ledger is frozen. Compare evidence, not label agreement. If the new run fails any priority-lineage gate, reopen the responsible pool and iterate.
