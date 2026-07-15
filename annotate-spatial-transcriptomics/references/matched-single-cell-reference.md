# Matched single-cell reference policy

A same-project or biologically matched single-cell dataset is the preferred external reference for a spatial query because it is closer in species, tissue, stage and laboratory context than a public atlas. It is still evidence, not ground truth. Harmonize the **meaning of broad classes**, not the number of labels or the final cell assignments.

## Required inputs and registry

When a matched reference is available, register all available forms separately:

- full count-level reference object and stable cell IDs;
- source annotation column and its original vocabulary;
- per-source-label cell counts and donor/sample composition;
- one-vs-rest DEG or pseudobulk marker table;
- canonical marker dotplot, if supplied;
- species, tissue, stage, condition, platform and dissociation protocol;
- whether the reference is from the same animal, same cohort, same stage or only the same tissue.

A dotplot alone can improve candidate marker programs and label crosswalks, but it cannot support cell-level transfer, RCTD or atlas-style calibration. Those routes require an expression matrix and source membership.

Create `config/matched_reference_crosswalk.tsv` from `assets/matched_reference_crosswalk_template.tsv` and validate it with `scripts/validate_matched_reference_crosswalk.py`. Preserve the source label verbatim. Never overwrite it with the harmonized name.

## Reference priority

Use the following evidence order for broad-class adjudication:

1. coherent full-feature evidence and high-confidence anchors from the current query;
2. a same-project/same-species/same-tissue, stage-matched single-cell reference;
3. an independent same-species/same-tissue atlas;
4. a developmental, cross-tissue or cross-species reference.

The matched reference is normally the strongest **external** channel, but it does not outrank contradictory current-query marker/anti-marker evidence or impossible spatial morphology. A same-project reference may share ambient RNA, dissociation bias or annotation assumptions, so it is not an independent second vote unless its marker program and provenance are audited.

## Harmonization rules

- Map detailed reference labels to a shallow candidate broad vocabulary before transfer. Many single-cell immune labels may map to `Immune`; blood and lymphatic endothelial labels may map to `Vascular/endothelial` while retaining subtype candidates.
- Keep evidence-dependent lineages such as `Theca`, `Smooth muscle`, `Pericyte/mural`, `Neural/Schwann` and strict `Oocyte` behind their tissue-profile gates even when the reference contains those names.
- A combined source label such as `APCs&B` has a broad-only ceiling until APC and B-cell programs separate in the current query.
- A source label whose displayed markers do not support its name is a review hypothesis, not a transferable identity. For example, neuronal genes without `CHGA/CHGB/SYP/SCG/INSM1` do not establish a neuroendocrine lineage.
- Preserve absent reference classes as negative audits. Do not lower a query gate to reproduce the reference taxonomy.
- Do not merge reference cells into the spatial query for final DEG, detection percentages or marker dotplots. Final broad DEG/dotplots use every spatial observation formally assigned to that broad class; subtype assets use validated high-confidence fine-label cells only.

## Transfer and calibration

Build depth-matched held-out reference cells only to diagnose reference separability. Final rescue thresholds require disjoint held-out current-query anchors with frozen truth, target membership and an origin manifest. Derive score and margin thresholds per route and candidate broad class from those query-like anchors. Reference self-classification is never an acceptable final calibration target.

For sheep ovary, a count-level, stage-compatible matched reference is the preferred external channel after current-query anchors **when the residual-QC Atlas/reference-rescue route is actually reached**. It may provide marker/anti-marker evidence elsewhere, but cell-level transfer is not run routinely on defined broad/fine cells or ordinary biological pools. If the matched artifact is only a dotplot, it remains marker/anti-marker evidence and GSE233801 is the primary public adult-sheep somatic Atlas only for the residual post-QC-anchor holdout. GSE233801 does not automatically rescue Oocyte, Theca or Epithelial/mesothelial.

The default crosswalk ceiling is `broad_only_after_calibration`. A reference prediction may return a broad label only when:

- its calibrated tier is moderate-or-higher;
- current-query marker/anti-marker evidence supports the same broad class;
- at least one independent query/internal-anchor or observed-density spatial channel agrees; and
- the result does not violate a rare/context-specific tissue gate.

Transferred observations set `fine_anchor_eligible=false`. A fine label requires current-query full-feature evidence, pool-specific stability and morphology; reference agreement is supporting evidence only.

## Reporting

The report must show:

- the source-to-harmonized crosswalk and transfer ceiling;
- source reference counts and marker evidence;
- calibrated high, moderate-only and low-reject counts by source/candidate label;
- discordant current-query marker or spatial evidence;
- broad labels rescued with matched-reference support;
- reference candidates rejected or retained as negative audits.

Display the matched single-cell dotplot as a reference panel, not as the spatial result. The final spatial broad DEG and canonical/data-specific marker dotplots must be recomputed from the complete accepted final spatial broad memberships.
