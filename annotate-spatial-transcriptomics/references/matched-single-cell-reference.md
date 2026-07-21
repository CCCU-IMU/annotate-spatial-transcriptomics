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

- Map detailed reference labels to a shallow candidate broad vocabulary before transfer. Many single-cell immune labels may map to `Immune`; blood endothelial, lymphatic endothelial and pericyte/mural labels map to the release broad parent `Vascular-associated` while retaining their fine candidates.
- Keep evidence-dependent lineages such as `Theca`, `Smooth muscle`, `Neural/Schwann` and strict `Oocyte` behind their tissue-profile gates even when the reference contains those names. Pericyte/mural remains evidence-gated as a fine identity under `Vascular-associated`, not as a separate release broad class.
- A combined source label such as `APCs&B` has a broad-only ceiling until APC and B-cell programs separate in the current query.
- A source label whose displayed markers do not support its name is a review hypothesis, not a transferable identity. For example, neuronal genes without `CHGA/CHGB/SYP/SCG/INSM1` do not establish a neuroendocrine lineage.
- Preserve absent reference classes as negative audits. Do not lower a query gate to reproduce the reference taxonomy.
- Do not merge reference cells into the spatial query for final DEG, detection percentages or marker dotplots. Final broad DEG/dotplots use every spatial observation formally assigned to that broad class; subtype assets use validated high-confidence fine-label cells only.

## Transfer and calibration

Build depth-matched held-out reference cells only to diagnose reference separability. Final rescue thresholds require disjoint held-out current-query anchors with frozen truth, target membership and an origin manifest. Derive score and margin thresholds per route and candidate broad class from those query-like anchors. Reference self-classification is never an acceptable final calibration target.

For sheep ovary, a count-level, stage-compatible matched reference is the preferred external channel **only after every broad/targeted cohort is terminal and residual QC is frozen**. Run one broad-only mapping over the complete analysis set. Unlabeled frozen QC is eligible for calibrated broad-only rescue; defined broad/fine cells use the same mapping only for concordance challenge. A dotplot-only artifact remains marker evidence and cannot support all-cell mapping. Otherwise GSE233801 is the primary public adult-sheep somatic Atlas. It does not automatically rescue Oocyte, Theca or Epithelial/mesothelial.

The default crosswalk ceiling is `broad_only_after_calibration`. A calibrated Atlas prediction may return a broad label to an unlabeled frozen-QC observation only when:

- its calibrated tier is moderate-or-higher;
- it is not OOD or ontology-conflicted;
- the candidate is within the declared Atlas scope and does not violate a rare/context-specific tissue gate.

Current-query marker/anti-marker, internal-anchor and spatial evidence remain mandatory audit columns and can trigger a group-level veto/review when they show coherent contradiction. They are not repeated as per-cell prerequisites after the Atlas tier has already been calibrated on disjoint query-like anchors.

Transferred observations set `fine_anchor_eligible=false`. A fine label requires current-query full-feature evidence, cohort-specific stability and morphology; reference agreement is supporting evidence only.

Treat mapping as open-set. Store maximum support, margin, neighbor composition and OOD status. A coherent low-support or mixed-neighbor cluster is an unknown-lineage hypothesis; never force it to the nearest reference class. Precompute and hash a fixed feature transform, low-dimensional reference representation and approximate-neighbor index for reuse. Query projection and all-cell label comparison are cheap; per-sample reference retraining, dense pairwise distances and whole-object RCTD are forbidden defaults.

## Reporting

The report must show:

- the source-to-harmonized crosswalk and transfer ceiling;
- source reference counts and marker evidence;
- calibrated `high`, `moderate_only` and `low_reject` counts by source/candidate label;
- discordant current-query marker or spatial evidence;
- broad labels rescued with matched-reference support;
- reference candidates rejected or retained as negative audits.

Display the matched single-cell dotplot as a reference panel, not as the spatial result. The final spatial broad DEG and canonical/data-specific marker dotplots must be recomputed from the complete accepted final spatial broad memberships.
