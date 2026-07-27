# R-first annotation workflow

Use this path when a readable full-feature Seurat RDS is available. For verified same-batch StereoPy cellbin inputs, Seurat carries expression and BANKSY may supply the whole-tissue partition.

## 1. Freeze and inspect

Run the object/runtime/feature-scope audits. Record object and observation hashes, assays/layers, raw counts, normalized data, reductions and coordinates. Freeze exact analysis and excluded-initial-QC memberships. A full feature list does not prove that `data` is normalized; create a project-local full-feature LogNormalize validation object when necessary.

Existing clustering can be reused only after input, membership, grid and artifact validation. Hide all historical annotation columns. Reusing computation never authorizes labels.

## 2. Whole-tissue partition

For a fresh verified same-batch sheep-ovary cellbin input, use SCT v2/glmGamPoi with 4,000 HVGs followed by BANKSY `M=0`, `k_geom=30`, `lambda=0.2`, 30 PCs, Leiden k=50 and `0.2,0.4,0.6,0.8`. A reusable input keeps its complete bound grid.

Select a stable spatial/anatomical partition without obvious technical fragmentation. Write only exact initial-cluster membership, one cohort ID per cluster, provisional program summary and lineage watches. Do not write initial biological labels or QC.

## 3. Second-round cohorts

Every initial cluster receives its own immutable cohort. For each cohort extract the same project's non-SCT raw counts, preferring `RNA` then `Spatial`, then run:

`CreateSeuratObject(raw counts) -> SCTransform(v2, glmGamPoi) -> PCA -> FindNeighbors(SNN) -> Leiden grid`

Never use the SCT assay as input to a new SCTransform and never reuse the whole-tissue BANKSY graph. For sheep ovary evaluate `0.1,0.2,0.3,0.4,0.6`, adapting PCs and k. Compute DEG, pseudobulk, spatial morphology and full-catalog evidence at every resolution. Select by stable identity separation and anatomy; use nearest neighboring resolutions for stability.

Score every selected subcluster against the full open catalog while its provisional first-pass record is hidden. After score freeze, return a pure subcluster wholesale, reconstruct a missing/cross-lineage broad, propose fine/state/unmodeled programs, trigger local mixed-subcluster splitting, or retain unresolved biology. An underpowered cohort cannot silently inherit the provisional name.

## 4. Local mixtures and broad freeze

Only a second-round subcluster with two or more separable competing identity cores enters observation-level splitting. Independently form candidate components, adjudicate overlaps and rescore the exact local remainder once. Unselected cells do not become QC.

After all cohorts finish, merge an exact disjoint analysis-set cover and freeze broad membership. Fine labels are materialized only inside the matching frozen broad parent; state annotations remain separate.

## 5. Atlas, completeness and release

Run one all-cell broad Atlas mapping after broad freeze. Rescue only unlabeled moderate-or-higher, non-OOD and profile-compatible observations. Existing broad labels are comparison-only unless a material conflict triggers one source-cohort biological review. Atlas never writes fine labels.

Audit every present and zero-census broad, embedded/unmodeled programs and every broad-by-fine candidate. Convert terminal unresolved members to typed QC only during final materialization. The final QC census must be below both 10% and 50,000.

## 6. Ovary-specific boundaries

- `Vascular-associated` contains endothelial and pericyte/mural children.
- Mature nonvascular Smooth muscle requires MYH11/CNN1/ACTG2-centered identity and mural exclusion; ACTA2/TAGLN alone is insufficient.
- Epithelial/mesothelial cannot expand from one keratin, one component or an aggregate mixed-cluster score.
- Oocyte uses a complete canonical-cluster rule and excludes zona/ambient expansion.
- Theca and luteal identities require independent programs and compatible anatomy/context.
- Granulosa hypoxia, proliferation, atresia and luteinization remain states unless an independent stable identity exists.
- Catalog-external stable programs are recorded as unmodeled candidates rather than forced into known labels.
