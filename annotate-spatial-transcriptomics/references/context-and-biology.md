# Biological context and tissue profiles

Collect species, tissue, stage/condition, platform, observation unit, section anatomy, biological questions and reference datasets before choosing pools or subtype resolutions. Missing age/cycle information must reduce confidence for stage-specific ovarian labels.

A context profile supplies marker programs, contradictory programs, candidate audit axes and safety alerts. It is not a mapping table. Derive final labels from the current query. Follow `taxonomy-and-pool-design.md`: published taxonomies are checklists, analysis pools are routing containers and release classes are independently gated biological conclusions.

For sheep ovary, load the packaged `references/profiles/sheep_ovary.json`. Treat granulosa/theca/follicular stroma as a review axis before release decisions; keep generic stroma, mature smooth muscle, pericyte/mural and vascular lineages separable; use immune and epithelial review pools only when needed; and handle Oocyte through a strict candidate route. Do not require all nine classes reported by a multi-stage sheep atlas to occur in one adult spatial section.

## Oocyte safety gate

Never write an entire cluster as Oocyte from ZP2/ZP3/ZP4 alone. Require cell-level co-detection of a multi-gene program including at least one identity/maternal-core marker, inspect contradictory granulosa/stromal/vascular/epithelial/immune programs, verify compact follicular spatial foci and recluster the strict candidate pool. Calibrate numerical thresholds from the query and known anchors. A prevalence above the profile alert is a mandatory review trigger, not automatic rejection.

Split an overbroad candidate into strict Oocyte candidates, oocyte-adjacent follicular somatic/stromal observations, and unresolved interface observations. Observation counts are not inferred oocyte counts.

Use `scripts/screen_rare_cell_programs.R` on the full-feature object for cell-level positive and contradictory program screening, then `scripts/screen_spatial_foci.py` for compact-focus evidence and multi-bin object grouping. When depth varies, `scripts/calibrate_rare_cell_candidates.py` derives positive/contradictory thresholds from the query distribution and expands only focus groups containing strict seeds. Render the whole-section audit with `scripts/plot_rare_cell_foci.py`. Feed only the calibrated strict focus set into a new frozen pool and recluster it. Neither screen, a spatial focus nor a single-cluster reclustering result is authorized to assign the final label by itself; compare it against adjacent/background observations and review morphology.

Use biological questions to set prioritization and stopping criteria. A project prioritizing folliculogenesis must not close a mixed granulosa/theca pool at broad level, while a broad tissue census may retain lower-priority rare states with explicit uncertainty.
