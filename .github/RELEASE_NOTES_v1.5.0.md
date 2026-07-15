# Annotate Spatial Transcriptomics v1.5.0

This release makes one evidence-gated annotation the only publication surface and tightens the sheep-ovary unresolved-cell controller.

## One final annotation

- Broad labels require calibrated moderate-or-higher confidence.
- Fine labels require high confidence and a valid fine-anchor route.
- Broad-only rescues are included in their final broad DEG, dotplots, UMAP and spatial statistics but excluded from fine-marker discovery.
- Release metadata binds broad/fine assets to `primary_broad_label` and `primary_subtype_label`; broad-only or retained-QC observations can never be fabricated into subtype labels.
- The HTML follows the validated Scanpy-style layout and does not publish competing strict/inclusive/display views.
- A mandatory lightweight pre-confirmation HTML now shows support/anti-marker reasons, a distinguishable broad spatial projection and a canonical broad marker dotplot. Full DEG, tree dotplots and release assets remain post-confirmation only.

## Residual-QC-only Atlas rescue

- RCTD medium/low observations enter the frozen QC holdout; they cannot call Atlas directly.
- The complete QC holdout first undergoes balanced-anchor query-only reclustering on the formal profile grid.
- Only the exact residual-QC child membership can enter Atlas/internal-anchor review.
- The route registry verifies the preceding QC-anchor route ID, outcome hash, parent/child pool snapshots and identical residual membership hash.
- GSE233801 remains the primary public adult-sheep somatic Atlas when appropriate, but it never runs routinely over the full object, defined broad/fine cells or ordinary biological pools.

## Reproducibility and orchestration

- Sheep-ovary whole-tissue and pool candidates use the complete `0.1,0.2,0.3,0.4,0.6` grid and reject ultra-low shortcuts.
- Generated jobs receive syntax/parameter preflight, failures enter an incident registry and open incidents block release.
- One user-facing master Agent coordinates the cohort; one full-workflow worker Agent owns each sample with the same completion and release gates.
