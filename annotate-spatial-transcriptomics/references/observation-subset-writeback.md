# Observation-level subset writeback

## Principle

A Leiden/Louvain/BANKSY cluster is a computational query boundary, not an indivisible biological unit. Whole-subcluster writeback is only a high-purity fast path. When a subcluster contains two coherent resident programs, freeze disjoint observation-level subsets and adjudicate the exact remainder independently.

This permits, for example, a Stromal-dominant follicular-wall subcluster to return a mature Smooth muscle subset while the remainder returns to Stromal/mesenchymal. The aggregate subcluster winner does not constrain a supported subset; the subset's own complete-catalog numerical evidence is authoritative.

## Required execution

1. Score every observation against the complete broad candidate catalog on the project-local full-feature object. Use lineage-specific identity/core and an independent support family, plus anti-programs. ACTA2/TAGLN, one keratin, one ECM marker, spatial location or one paper marker cannot define a subset.
2. Propose a target subset only when its per-observation winner is stable and its exact membership has at least two positive families, low hard anti-program burden, a positive score margin and coherent spatial structure.
3. Validate the frozen subset with `observation_subset_evidence.schema.json`. The default gate is at least 0.70 lineage-supported observations, at least 0.30 support advantage over the strongest competitor and at most 0.05 contradiction fraction. Thresholds are project-bound in `observation_writeback_policy`.
4. Allow multiple disjoint biological returns from the same source subcluster. Their memberships plus targeted/QC successors must exactly partition the source membership.
5. Remove accepted subsets and rescore the exact remainder against the complete catalog. The remainder may return to any supported broad class or remain QC; parent identity, aggregate winner and spatial adjacency are not priors.
6. Preserve shared programs such as contractility, ECM, hypoxia, stress and follicle adjacency as state tags. A double-positive interface remains broad-only or QC when no subset passes.

## Whole-subcluster gate

Whole-subcluster writeback requires all of the following: complete-catalog positive-margin winner, at least two independent positive families, no unresolved contradiction, raw two-family support of at least 0.40 and a raw support advantage of at least 0.20. Every other eligible broad lineage with raw support at least 0.35 triggers an embedded-component audit; the whole return passes only when that competitor is shown to be spatially incoherent or contradicted. A coherent competitor requires `supported_subset` extraction or a targeted successor. This avoids treating ubiquitous ECM in Granulosa cellbins as true Stroma while still reopening the large follicular-ring Smooth muscle signal. When a legacy table lacks raw support, the conservative fallback is effective support at least 0.25 with margin at least 0.10. Relative winner status alone never passes, and project-specific calibrated values are frozen in the annotation contract.

## Final audit

The post-Atlas broad-class completeness artifact must enumerate every source writeback that contributes to the final class. Query-derived whole/subset returns retain their own purity metrics and membership hashes; calibrated Atlas and canonical Oocyte returns retain their route-specific validation. Counts must exactly reproduce the final class. This prevents a strong majority return from hiding one weak whole-cluster expansion.

## Fine candidates

After broad labels are locked, audit the complete parent-specific fine-candidate catalog, not only vascular children. Every optional child is recorded as `supported`, `refuted` or `not_evaluable`; zero fine labels is valid only after this audit. Fine labels never repair an incorrect broad label, and state-only splits remain tags.
