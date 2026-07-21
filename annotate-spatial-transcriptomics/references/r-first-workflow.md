# Seurat/R-first annotation workflow

Use this path whenever a readable full-feature Seurat RDS is available or R is selected as the primary framework. Seurat is the computational backbone; the biological controller is `direct-lineage-controller.md`.

## 1. Freeze and inspect the object

Run `inspect_r_object.R`, `check_r_runtime.R` and `audit_feature_scope.R`. Record object hash, Seurat version, assays/layers, raw-count and normalized layers, reductions, coordinates, observation IDs and analysis-set membership. SCT residuals may drive PCA/clustering but cannot establish marker absence.

Prefer the full-feature RDS. A full feature list does not prove `Spatial@data` is normalized. If `data` is absent or identical to `counts`, leave the SCT clustering object immutable and build a manifest-bound full-feature LogNormalize validation object with `prepare_seurat_full_feature_validation.R` before Wilcoxon DEG or marker/anti-marker decisions.

A matched single-cell object is a separate reference and is never merged into the query ledger.

## 2. Reuse computation without leaking labels

Existing whole-tissue or cohort clustering is reusable only when source-object and membership hashes, analysis scope, candidate resolutions, cluster memberships, DEG, UMAP, spatial outputs and full-feature validation are complete. Register it as `validated_reuse`, hide historical annotation columns and create a new decision ledger. Reusing computation never authorizes copying labels.

## 3. Whole-tissue broad pass

When no valid clustering is reusable, run SCTransform, PCA and an adaptive Seurat resolution grid. Sheep ovary uses the formal grid `0.1,0.2,0.3,0.4,0.6`; never lower below 0.1 to compensate for a defective graph. Other tissues use a declared context-appropriate grid.

Generate per-resolution one-vs-rest DEG, canonical and data-specific marker summaries, UMAP, whole-section spatial maps and per-cluster highlights. Select the integrated-evidence optimum that preserves supported lineages before avoiding state-only fragmentation; use lower complexity only when evidence is otherwise equivalent.

At that resolution, perform open-world lineage review. Supported clusters receive a moderate-or-higher `initial_broad_label` directly. Low-information, featureless or irreducibly mixed observations enter `qc_holdout`; do not create intermediate biological memberships.

## 4. Broad-class and targeted cohorts

Create one immutable `broad_class_recluster` cohort for every supported initial broad class. A genuinely underpowered class may remain broad-only only after a recorded `not_applicable_reviewed` decision. Use `run_seurat_cohort_recluster.R` and register its output in `recluster_cohort_registry.tsv`.

Fit normalization/PCA jointly with frozen anchors when they are needed, while constructing graph, clusters, UMAP, DEG and outcome counts from query observations only. Select cohort PCs, k and resolution from the current membership; for sheep ovary run the full formal grid for every broad or targeted cohort.

Adjudicate every subcluster as exactly one of:

- high-confidence shallow fine label;
- direct return to the parent broad class;
- direct cross-lineage broad/fine return with source provenance;
- one decision-relevant `targeted_recluster` cohort;
- residual QC/technical retention.

Do not create an intermediate cohort. A direct cross-lineage return does not automatically enter the target class's reclustering cohort again.

## 5. Assistance and residual QC

Use a targeted cohort only for a local interpretable mixture, contamination boundary or context-gated identity. RCTD is lower-priority assistance: canonical high confidence plus independent marker/anti-marker, resolution and spatial evidence may support fine; moderate supports broad-only; low enters the final QC holdout.

After all broad and targeted cohorts are terminal, freeze residual QC. Do not recluster it. Run one calibrated broad-only Atlas mapping over the complete analysis set. Apply marker/internal-anchor/observed-density consensus only when writing back the frozen-QC subset; compare defined labels against the same mapping and reopen only material broad disagreement or coherent OOD once. Moderate-or-higher QC returns are broad-only with `fine_anchor_eligible=false`; lower confidence remains QC reject.

For sheep ovary without a usable matched count-level reference, GSE233801 is the primary public adult-sheep somatic Atlas for this terminal residual-QC step. It cannot bypass strict Oocyte, Theca or epithelial evidence gates.

## 6. Ovary stromal-lineage decomposition

Before closing a generic stromal parent, independently test generic stromal/fibroblast, mesenchymal-progenitor-like, mature smooth muscle, pericyte/mural, blood/lymphatic endothelial and steroidogenic-theca programs.

- Mature Smooth muscle requires a coherent `MYH11/CNN1/ACTG2` plus `TAGLN/ACTA2/MYL9/DES` backbone and compatible tracks.
- Pericyte/mural requires `RGS5/PDGFRB/CSPG4/NOTCH3/MCAM`-like support and vascular adjacency and is released as a fine child of `Vascular-associated`.
- Steroidogenic Theca requires a coherent steroidogenic/androgenic program; ECM or `LHCGR` alone is insufficient.
- ECM-rich `ACTA2/TAGLN` without a mature backbone remains structural/contractile stroma.
- `Mesenchymal progenitor-like` is separate only when a stable progenitor/fibroblast program and morphology pass; otherwise use `Stromal/mesenchymal`.

Negative audits are valid outcomes and must not be repaired by lowering gates.

## 7. Shallow subtype policy

Broad classes are the default endpoint for cellbin spatial data. Merge subclusters driven mainly by depth, LOC/ribosomal/mitochondrial expression, ECM amount, stress, cell cycle, generic contractility or spatial position. Keep those characteristics as tags.

Adult sheep ovary permits only evidence-supported shallow functional labels: selected Granulosa states, steroidogenic Theca, blood/lymphatic endothelial and Pericyte/mural under `Vascular-associated`, mature Smooth muscle, and depth-supported immune divisions. Zero fine labels is valid.

## 8. Quality and release

Build one final annotation: broad labels require moderate-or-higher confidence; fine labels require high confidence. Every direct and Atlas-rescued broad observation participates in final broad DEG and dotplots. Atlas broad-only cells and RCTD-only cells cannot seed fine marker discovery.

Run `validate_direct_lineage_workflow.py`, state validation, taxonomy audit and completion gate before main-Agent biological quality approval. User confirmation follows the lightweight report; only then generate final DEG, tree dotplots, spatial assets and release HTML.
