# Biological context and tissue profiles

Collect species, tissue, stage/condition, platform, observation unit, section anatomy, biological questions and reference datasets before choosing pools or subtype resolutions. Missing age/cycle information must reduce confidence for stage-specific ovarian labels.

After input inspection, run `scripts/resolve_workflow_profile.py`. Sheep/Ovis/ovine/羊 plus ovary/ovarian/卵巢 selects the sheep-ovary profile; a readable full-feature Seurat RDS selects R-first, and confirmed StereoPy `cellbin_PPed` provenance plus a full-feature `Spatial` assay activates the fixed batch preprocessing contract. Resolver output is workflow provenance only and cannot assign cell identities.

If the user explicitly requests the same-batch standard process, pass `--strategy-preset sheep_ovary_same_batch_rfirst` and write the resolver result to `config/active_strategy_preset.json`. This binds a deidentified successful decision process, not a sample answer: phase order and safeguards are fixed, while resolutions, memberships and biological labels remain current-query decisions.

A context profile supplies marker programs, contradictory programs, candidate audit axes and safety alerts. It is not a mapping table. Derive final labels from the current query. Follow `taxonomy-and-pool-design.md`: published taxonomies are checklists, analysis pools are routing containers and release classes are independently gated biological conclusions.

If the project has corresponding single-cell data, read `matched-single-cell-reference.md`. Record whether it is the same animal, cohort, stage or merely the same tissue. A matched reference is the preferred external reference channel, but it never replaces current-query evidence. Preserve the original reference vocabulary and create an auditable crosswalk rather than renaming the reference in place.

For sheep ovary, load the packaged `references/profiles/sheep_ovary.json`, `profiles/sheep_ovary_rfirst_profile.json` and `profiles/sheep_ovary_literature_2025_2026.md`. Treat granulosa/theca/follicular stroma as a review axis before release decisions; keep generic stroma, mature smooth muscle, pericyte/mural and vascular lineages separable; use immune and epithelial review pools only when needed; and handle Oocyte through a strict candidate route. Do not require all classes reported by a multi-stage or cross-species atlas to occur in one adult spatial section.

For a sheep-ovary matched reference, `profiles/sheep_ovary_reference_aliases.tsv` provides candidate aliases and evidence ceilings. It deliberately maps detailed immune and endothelial source labels to shallow spatial parents and flags weak or combined labels for review. Validate a project-specific copy and update it from the actual reference object, DEG and stage metadata.

## Oocyte safety gate

Never write an entire cluster as Oocyte from ZP2/ZP3/ZP4 alone. Require cell-level co-detection of a multi-gene program including at least one identity/maternal-core marker, inspect contradictory granulosa/stromal/vascular/epithelial/immune programs, verify compact follicular spatial foci and recluster the strict candidate pool. Calibrate numerical thresholds from the query and known anchors. A prevalence above the profile alert is a mandatory review trigger, not automatic rejection.

Do not use cortical, subcortical, peripheral or section-edge location as negative evidence against Oocyte. Small and primordial oocytes may occur in the ovarian cortex, so an edge-localized compact object must remain eligible for the same molecular and object-level review as an internal focus. Location is also not sufficient positive evidence and does not relax the multi-module, somatic anti-program or strict-reclustering gates.

Split an overbroad candidate into strict Oocyte candidates, oocyte-adjacent follicular somatic/stromal observations, and unresolved interface observations. Observation counts are not inferred oocyte counts.

For Oocyte candidate validation, use the legacy-named `scripts/screen_rare_cell_programs.R` only as an Oocyte multi-module/contradictory-program screen on the full-feature object, then `scripts/screen_spatial_foci.py` for compact-focus evidence and multi-bin object grouping. Its filename is retained for project compatibility and does not authorize a generic rare-cell route. When depth varies, `scripts/calibrate_rare_cell_candidates.py` derives positive/contradictory thresholds from the query distribution while retaining every observation that passes the predeclared multi-module starting gate. Render the whole-section audit with `scripts/plot_rare_cell_foci.py`. Feed the complete calibrated starting pool—not only strict seeds or spatial foci—into a new frozen query-only pool and recluster it. Neither the screen, a spatial focus nor a single-cluster separation is authorized to assign the final label by itself; compare cluster-level non-ZP/maternal programs, somatic anti-programs, adjacent/background observations and morphology.

Use biological questions to set prioritization and stopping criteria. A project prioritizing folliculogenesis must not close a mixed granulosa/theca pool at broad level, while a broad tissue census may retain lower-priority rare states with explicit uncertainty.
