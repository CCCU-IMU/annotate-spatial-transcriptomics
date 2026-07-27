#!/usr/bin/env python3
"""Select whole-tissue or cohort resolution by one executable biology score."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from lineage_controller_lib import number, read_tsv, sha256, truth, write_manifest, write_tsv


METRIC_WEIGHTS = {
    "whole_tissue_cohort_partition": {
        "catalog_recall": 0.10,
        "embedded_program_separation": 0.05,
        "deg_antideg_coherence": 0.10,
        "pseudobulk_coherence": 0.05,
        "spatial_morphology_coherence": 0.20,
        "adjacent_resolution_stability": 0.30,
        "observation_support_coherence": 0.10,
    },
    "cohort_identity_resolution": {
        "catalog_recall": 0.05,
        "embedded_program_separation": 0.05,
        "deg_antideg_coherence": 0.05,
        "pseudobulk_coherence": 0.05,
        "spatial_morphology_coherence": 0.05,
        "adjacent_resolution_stability": 0.10,
        "observation_support_coherence": 0.05,
        "directly_resolved_observation_fraction": 0.40,
        "mean_identity_margin": 0.20,
    },
}
COMPOSITE_EQUIVALENCE_TOLERANCE = 0.01
PER_METRIC_EQUIVALENCE_TOLERANCE = 0.02


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grid-evidence", required=True, type=Path)
    parser.add_argument(
        "--selection-purpose", required=True,
        choices=sorted(METRIC_WEIGHTS),
    )
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    rows = read_tsv(args.grid_evidence)
    required = {
        "resolution", "complete_catalog_scanned", "zero_census_audited",
        "technical_fragmentation", "state_overfragmentation", "complexity",
        "selection_purpose", *METRIC_WEIGHTS[args.selection_purpose],
    }
    if args.selection_purpose == "cohort_identity_resolution":
        required.update({
            "mixed_observation_fraction", "unresolved_observation_fraction",
        })
    if not rows or not required.issubset(rows[0]):
        raise SystemExit(
            "resolution evidence is empty or lacks: "
            + ", ".join(sorted(required - set(rows[0] if rows else {})))
        )
    evaluated = []
    metrics = METRIC_WEIGHTS[args.selection_purpose]
    for row in rows:
        if row.get("selection_purpose") != args.selection_purpose:
            raise SystemExit("resolution evidence purpose differs from selector purpose")
        eligible = truth(row["complete_catalog_scanned"]) and truth(row["zero_census_audited"])
        biological = sum(number(row[key]) * weight for key, weight in metrics.items())
        penalty = (
            (0.20 if args.selection_purpose == "whole_tissue_cohort_partition" else 0.10)
            * number(row["technical_fragmentation"])
            + (0.10 if args.selection_purpose == "whole_tissue_cohort_partition" else 0.15)
            * number(row["state_overfragmentation"])
        )
        if args.selection_purpose == "cohort_identity_resolution":
            penalty += (
                0.20 * number(row["mixed_observation_fraction"])
                + 0.30 * number(row["unresolved_observation_fraction"])
            )
        evaluated.append({
            **row,
            "resolution": number(row["resolution"]),
            "eligible": eligible,
            "biological_score": biological,
            "penalty": penalty,
            "selection_score": biological - penalty,
            "complexity": number(row["complexity"]),
        })
    resolution_values = sorted({number(row["resolution"]) for row in evaluated})
    eligible_rows = [row for row in evaluated if row["eligible"]]
    if not eligible_rows:
        raise SystemExit(
            "no resolution completed the catalog and zero-census audit"
        )
    eligible_rows.sort(
        key=lambda row: (-row["selection_score"], row["complexity"], row["resolution"])
    )
    best = eligible_rows[0]
    best_score = best["selection_score"]
    equivalent = [
        row for row in eligible_rows
        if (
            best_score - row["selection_score"]
            <= COMPOSITE_EQUIVALENCE_TOLERANCE
            and all(
                number(best[metric]) - number(row[metric])
                <= PER_METRIC_EQUIVALENCE_TOLERANCE
                for metric in metrics
            )
        )
    ]
    selected = min(equivalent, key=lambda row: (row["complexity"], row["resolution"]))
    selected_resolution = number(selected["resolution"])
    lower = [value for value in resolution_values if value < selected_resolution]
    higher = [value for value in resolution_values if value > selected_resolution]
    neighbors: list[float] = []
    if lower:
        neighbors.append(lower[-1])
    if higher:
        neighbors.append(higher[0])
    for value in sorted(
        (value for value in resolution_values if value != selected_resolution),
        key=lambda value: (abs(value - selected_resolution), value),
    ):
        if value not in neighbors:
            neighbors.append(value)
        if len(neighbors) == 2:
            break
    if len(neighbors) != 2:
        raise SystemExit("resolution stability review requires at least three grid candidates")
    result = {
        "status": "PASS",
        "schema_version": "2.2",
        "selector_version": "2.2.0",
        "selection_purpose": args.selection_purpose,
        "selected_resolution": selected_resolution,
        "selected_and_neighbors": [selected_resolution, *neighbors],
        "neighbor_strategy": "nearest_lower_and_higher_else_two_nearest",
        "selection_score": selected["selection_score"],
        "lower_complexity_used_only_as_tiebreaker": len(equivalent) > 1,
        "equivalence_policy": {
            "composite_tolerance": COMPOSITE_EQUIVALENCE_TOLERANCE,
            "per_metric_tolerance": PER_METRIC_EQUIVALENCE_TOLERANCE,
        },
        "grid_evidence": {
            "path": str(args.grid_evidence.resolve()),
            "sha256": sha256(args.grid_evidence),
        },
    }
    args.out.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out / "resolution_ranking.tsv", evaluated)
    write_manifest(args.out / "resolution_selection.json", result)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
