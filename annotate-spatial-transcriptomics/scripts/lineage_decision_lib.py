#!/usr/bin/env python3
"""Deterministic helpers for cluster-by-lineage decision tables.

The table is the numerical authority for both initial broad decisions and
reclustered-subcluster returns.  Human/Agent prose may explain a decision but
cannot replace or override the ranking encoded here.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from typing import Iterable


REQUIRED_COLUMNS = {
    "cluster",
    "candidate_id",
    "candidate_broad_lineage",
    "program_score",
    "positive_family_count",
    "positive_families",
    "anti_program_burden",
    "contradiction_count",
    "eligible",
}

PURITY_COLUMNS = {
    "purity_status",
    "lineage_supported_fraction",
    "strongest_competing_fraction",
}


def truth(value: object) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "pass", "passed"}


def split_families(value: object) -> list[str]:
    return sorted({item.strip() for item in re.split(r"[;,|]+", str(value or "")) if item.strip()})


def _number(row: dict[str, str], key: str) -> float:
    value = float(row.get(key, "nan"))
    if not math.isfinite(value):
        raise ValueError(key)
    return value


def _integer(row: dict[str, str], key: str) -> int:
    value = _number(row, key)
    if value != int(value):
        raise ValueError(key)
    return int(value)


def validate_rows(
    rows: list[dict[str, str]],
    *,
    require_purity: bool = False,
    expected_clusters: Iterable[str] | None = None,
    expected_candidates: Iterable[str] | None = None,
) -> list[str]:
    errors: list[str] = []
    if not rows:
        return ["lineage decision table is empty"]
    fields = set(rows[0])
    required = REQUIRED_COLUMNS | (PURITY_COLUMNS if require_purity else set())
    if missing := sorted(required - fields):
        return ["lineage decision table lacks columns: " + ", ".join(missing)]
    keys = [(str(row.get("cluster", "")).strip(), str(row.get("candidate_id", "")).strip()) for row in rows]
    duplicates = sorted(key for key, count in Counter(keys).items() if count != 1)
    if duplicates:
        errors.append(f"lineage decision table has duplicate cluster/candidate rows: {duplicates[:5]}")
    if any(not cluster or not candidate for cluster, candidate in keys):
        errors.append("lineage decision table has an empty cluster or candidate_id")
    clusters = {cluster for cluster, _ in keys if cluster}
    candidates_by_cluster = {
        cluster: {candidate for observed_cluster, candidate in keys if observed_cluster == cluster}
        for cluster in clusters
    }
    if expected_clusters is not None and clusters != {str(value) for value in expected_clusters}:
        errors.append("lineage decision table cluster universe differs from the selected partition")
    if expected_candidates is not None:
        expected = {str(value) for value in expected_candidates}
        for cluster, observed in candidates_by_cluster.items():
            if observed != expected:
                errors.append(f"cluster {cluster} does not cover the complete candidate catalog")
    elif candidates_by_cluster:
        first = next(iter(candidates_by_cluster.values()))
        for cluster, observed in candidates_by_cluster.items():
            if observed != first:
                errors.append(f"cluster {cluster} has a different candidate universe")
    for line, row in enumerate(rows, 2):
        try:
            _number(row, "program_score")
            anti = _number(row, "anti_program_burden")
            positive_n = _integer(row, "positive_family_count")
            contradiction_n = _integer(row, "contradiction_count")
            families = split_families(row.get("positive_families", ""))
            if anti < 0 or positive_n < 0 or contradiction_n < 0:
                errors.append(f"lineage decision row {line} contains a negative evidence count/burden")
            if positive_n != len(families):
                errors.append(f"lineage decision row {line} positive-family count is inconsistent")
            if require_purity:
                supported = _number(row, "lineage_supported_fraction")
                competing = _number(row, "strongest_competing_fraction")
                if not 0 <= supported <= 1 or not 0 <= competing <= 1:
                    errors.append(f"lineage decision row {line} purity fractions are outside [0,1]")
                if str(row.get("purity_status", "")).strip().lower() not in {"pass", "mixed", "fail", "not_applicable"}:
                    errors.append(f"lineage decision row {line} has invalid purity_status")
        except (TypeError, ValueError):
            errors.append(f"lineage decision row {line} contains invalid numeric evidence")
    return errors


def cluster_rows(rows: list[dict[str, str]], cluster: object) -> list[dict[str, str]]:
    key = str(cluster)
    return [row for row in rows if str(row.get("cluster", "")) == key]


def rank(rows: list[dict[str, str]]) -> dict[str, object]:
    if len(rows) < 2:
        raise ValueError("at least two candidate rows are required for deterministic ranking")
    ordered = sorted(rows, key=lambda row: (-_number(row, "program_score"), str(row.get("candidate_id", ""))))
    winner, runner = ordered[:2]
    return {
        "winner": str(winner["candidate_id"]),
        "runner_up": str(runner["candidate_id"]),
        "winning_margin": _number(winner, "program_score") - _number(runner, "program_score"),
        "winner_row": winner,
        "runner_up_row": runner,
    }


def eligible(row: dict[str, str], minimum_families: int = 2) -> bool:
    return (
        truth(row.get("eligible"))
        and _integer(row, "positive_family_count") >= minimum_families
        and _integer(row, "contradiction_count") == 0
    )


def purity_passes(row: dict[str, str]) -> bool:
    return (
        str(row.get("purity_status", "")).strip().lower() == "pass"
        and _number(row, "lineage_supported_fraction") > _number(row, "strongest_competing_fraction")
    )


def absolute_supported_families(
    rows: list[dict[str, str]],
    *,
    cluster: object,
    candidate_id: str,
    thresholds: dict,
) -> list[str]:
    minimum_detected = int(thresholds["minimum_detected_marker_count"])
    minimum_fraction = float(thresholds["minimum_any_detection_fraction"])
    minimum_expression = float(thresholds["minimum_mean_normalized_expression"])
    return sorted({
        row["family_name"]
        for row in rows
        if str(row.get("cluster", "")) == str(cluster)
        and row.get("candidate_id") == candidate_id
        and int(float(row.get("detected_marker_count", 0))) >= minimum_detected
        and float(row.get("any_detection_fraction", 0)) >= minimum_fraction
        and float(row.get("mean_normalized_expression", 0)) >= minimum_expression
    })
