# Iterative annotation controller

The first clustering is a proposal generator, not the final annotation. Execute this state machine until the completion gate passes.

## Round 0: context and input

Validate the biological context and matching profile. Discover every expression object, clustering, coordinate set and prior state. Freeze hashes. If state exists, resume open pools; never re-run a closed membership in the same decision version.

## Round 1: broad anchors and review pools

Select a broad clustering adaptively. Generate DEG, marker programs, anti-programs, UMAP and space. Write high-confidence broad anchors. Create immutable pools for coherent broad lineages, interfaces, QC holdout, technical observations and context-gated rare candidates.

An original cluster can be split at cell/bin level. Never force all observations in a cluster to inherit one biological label when positive and contradictory programs coexist.

## Round 2+: pool-specific analysis

For every open biological pool:

1. Combine queries with frozen relevant anchors when needed.
2. Test several pool-appropriate resolutions; shortlist quantitatively and choose the lowest resolution preserving stable biological compartments.
3. Generate per-resolution DEG, UMAP, spatial maps and cluster migration/stability evidence.
4. Decide each subcluster as defined, returned to a named parent pool, interface, QC/technical retention, RCTD/reference-assisted review or triggered context-specific validation. Do not route by rarity alone.
5. Write membership, route, run ID, iteration and evidence before closing the run.
6. Roll back splits driven only by ECM, stress, cell cycle, ribosomal/mitochondrial load or graph fragmentation.

Use `multi-route-controller.md` for the mandatory route sequence. For large post-clustering QC, Route C begins with a complete QC-pool anchor reclustering; atlas mapping is the second-stage rescue of remaining rejects. For interfaces, route applicability is machine-readable and cannot be waived by a narrative judgment.

Do not demand a fine subtype when biology or depth supports only a broad type. A broad-only closure needs a rationale and is never a fine anchor.

## Failure-mode routing

- Large usable unresolved population: anchor-assisted pool reclustering.
- Mixed local interface: targeted reclustering, then optional calibrated RCTD.
- Low-depth/QC holdout: complete broad-anchor QC-pool reclustering first; only its frozen residual QC proceeds to depth-matched Atlas/internal-anchor rescue; retain rejects.
- Rare identity: positive modules, negative programs, spatial objects, morphology and strict-pool reclustering.
- Irreducible/technical: explicit retained state, not a forced label.

## Mandatory loop

Run `plan_next_iteration.py` after every writeback. Submit its queued routes, inspect outputs and update state. A non-empty queue means the report is preliminary. Run `check_completion_gate.py` only after the queue is empty; a blocked gate requires another round.

The plan and completion gate must inspect route attempts, branch no-repeat state, the single final annotation and the distinction between full-object, analysis-set, initial-QC and post-clustering holdout ledgers. “All decisions closed” alone is not completion.

The Agent must narrate biological conclusions separately from automated scores, cite evidence artifacts, and record rejected alternatives. Numerical ranks and reference predictions are decision support, not label writers.

After RCTD, route extreme/high/medium-low confidence separately. Medium/low confidence enters the post-clustering QC holdout and cannot call Atlas directly. After the complete QC-holdout anchor reclustering, only residual QC may continue to calibrated Atlas/internal-anchor review. Reopen a terminal interface when it is too large or spatially diffuse to be anatomically local. Also reopen large direct first-pass labels when cell-level lineage purity and spatial continuity were never tested, even if their cluster-level DEG looked canonical.
