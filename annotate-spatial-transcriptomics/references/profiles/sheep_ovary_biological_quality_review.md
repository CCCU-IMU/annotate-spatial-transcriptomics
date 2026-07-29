# Sheep-ovary biological quality review

Run this review after second-round broad membership and the one all-cell Atlas pass are frozen, but before final materialization and user confirmation. It is a query-derived challenger and never a label source.

## Three biological endpoints

1. **Spatial localization of every released broad class.** Audit complete membership against its full-feature multigene identity, hard contradictions, connected morphology and profile-compatible anatomy. Treat branching vasculature, broad ovarian stroma and repeated follicles according to their own morphology; do not require one compact component.
2. **Oocyte annotation quality.** Require at least one complete canonical Oocyte subcluster with non-zona/maternal identity, somatic contradiction clearance and reconstructable spatial objects. Count cellbins and spatial objects separately. Cortical or section-edge location is neither negative nor sufficient positive evidence. A coherent zero-census Oocyte group reopens the canonical cohort; zona/adjacency never expands it.
   If an exact label-blind canonical targeted cohort has already been frozen, pass its adjudication manifest with `--canonical-oocyte-review`. The validator accepts it only when the reviewed membership contains the identical Oocyte cell-ID set, the cohort is cross-resolution stable, location was not used for admission, zona-only admission was forbidden and independent non-zona DEG evidence passed. Stale ordinary second-round Oocyte scores cannot veto that exact canonical result.
3. **Follicle ROI histology.** Detect Granulosa spatial components from the frozen result, assign nearby observations to follicle ROIs and review the complete stage-appropriate structure rather than one lineage in isolation. Large/antral candidates must show a coherent Granulosa boundary, a clear low-cell-density cavity, Theca interna with interleaved vascular cells, an outer mature nonvascular contractile layer and then ovarian stroma. Audit presence, angular continuity, radial order, label/program concordance, boundary blurring, missing layers and abnormal expansion. Small/preantral follicles are not required to have a resolved cavity or every mature wall layer.

## Follicle ROI evidence

- Build ROIs from project-local coordinates and frozen Granulosa membership. Granulosa labels locate the question; full-catalog query scores adjudicate it.
- Measure annularity, cavity density, angular coverage and signed distance from the local Granulosa boundary. Review Granulosa, Theca, endothelial/pericyte-mural, mature Smooth muscle and Stromal programs jointly across radial shells and angular sectors.
- A coherent specific program recurring in multiple sectors but retained mainly as generic Stromal or unresolved is `ITERATION_REQUIRED`, even when direct cellbin seeds fragment into small components.
- Separate **program detection** from **label recall**. Direct plus local signal can establish that a layer exists, but only directly discriminated identity cores (or a validated expression subcluster inheriting its sparse tail) form the denominator for required writeback. One mixed cellbin is never required to carry several mutually exclusive broad labels merely because its neighbours support several layers.
- Theca interna and capillary/mural cells may be interleaved. Do not require a complete pure vascular ring.
- Formal Theca interna uses the steroidogenic/androgenic Theca identity. A structural/perifollicular ring is an exploratory topology program and cannot be counted or written as Theca interna without that independent identity.
- Whole-section compactness is not an admission or exclusion rule for molecularly supported Theca: multiple follicles can produce many small, separated Theca foci. Review its recurrence and radial role inside follicle ROIs instead.
- Split an outer fibromuscular region only when `MYH11/CNN1/ACTG2` supports mature nonvascular Smooth muscle after mural exclusion. Keep the ECM/`TCF21/PDGFRA` remainder Stromal/mesenchymal with an optional contractile or theca-externa-like state.
- Resolve specific wall identities before the generic Stromal remainder. Compare directly supported Theca, Endothelial, Pericyte/mural and mature nonvascular Smooth-muscle identities independently; preserve unresolved specific-specific ties; only observations with no specific identity return to Stromal. Ovarian stroma is a broad background rather than a compulsory thin ring, so require that it extends beyond an observed contractile wall in the ROI, not that every follicle must display all layers or that its global median radius exceed every Smooth-muscle observation.
- The apparent basement membrane is assessed only as a Granulosa–Theca identity-boundary proxy because it is not itself a cell lineage. Geometry may trigger one raw-count anatomy-conditioned targeted cohort. Geometry, a ring, adjacency, ECM or `ACTA2/TAGLN` cannot assign a label.

## Conditional interpretation

- If no coherent Granulosa component exists because the sample contains no follicles, record follicle histology as `NOT_EVALUABLE`; do not fail a luteal or nonfollicular sample.
- If Granulosa is released but remains spatially diffuse without a coherent component, return to the contributing second-round cohort.
- If no large/antral ROI is detected, report cavity clarity as `NOT_EVALUABLE`; do not require an antral cavity from a sample containing only small follicles.
- A review problem returns the exact contributing second-round subcluster or one bounded follicle ROI to local analysis. It never reopens first-pass whole-object cellbin classification.
- A targeted repair binds the pre-iteration ROI review. A previously confirmed large/antral cavity and every layer that already passed are immutable anchors unless that same layer had a typed review failure. A repair that makes a mature ROI disappear is rejected, even if re-detection would otherwise call the sample `NOT_EVALUABLE`.
- A bounded ROI re-review may score only the reopened ROI, but spatial geometry must still come from the complete-section coordinate membership via `--coordinate-membership`. The coordinate artifact supplies no labels and has no writeback authority; it prevents artificial component fragmentation caused by cropping the coordinate universe.

Run `scripts/validate_sheep_ovary_biological_quality.py` on the post-Atlas membership plus every disjoint second-round observation-score table. Preserve its ROI membership and numerical tables for the lightweight review report.
