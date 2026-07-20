#!/usr/bin/env python3
"""Find small coherent lineage programs embedded inside large broad labels."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


def quantile(values: list[float], q: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.inf
    position = (len(ordered) - 1) * min(max(q, 0.0), 1.0)
    low = int(math.floor(position)); high = int(math.ceil(position))
    if low == high:
        return ordered[low]
    return ordered[low] * (high - position) + ordered[high] * (position - low)


def spatial_components(rows: list[dict[str, object]], radius: float) -> list[list[int]]:
    if not rows:
        return []
    parent = list(range(len(rows)))
    buckets: dict[tuple[int, int], list[int]] = defaultdict(list)

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        a, b = find(left), find(right)
        if a != b:
            parent[b] = a

    for index, row in enumerate(rows):
        x, y = float(row["x"]), float(row["y"])
        key = (math.floor(x / radius), math.floor(y / radius))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for other in buckets.get((key[0] + dx, key[1] + dy), []):
                    ox, oy = float(rows[other]["x"]), float(rows[other]["y"])
                    if (x - ox) ** 2 + (y - oy) ** 2 <= radius ** 2:
                        union(index, other)
        buckets[key].append(index)
    groups: dict[int, list[int]] = defaultdict(list)
    for index in range(len(rows)):
        groups[find(index)].append(index)
    return list(groups.values())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("scores", type=Path, help="TSV: cell_id,parent_label,candidate_lineage,program_score,core_hits,support_hits,anti_hits,x,y")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--positive-quantile", type=float, default=0.95)
    parser.add_argument("--radius", type=float, required=True)
    parser.add_argument("--min-component", type=int, default=25)
    parser.add_argument("--max-parent-fraction", type=float, default=0.10)
    args = parser.parse_args()
    with args.scores.open(newline="", encoding="utf-8") as handle:
        raw = list(csv.DictReader(handle, delimiter="\t"))
    required = {"cell_id", "parent_label", "candidate_lineage", "program_score", "core_hits", "support_hits", "anti_hits", "x", "y"}
    if not raw or required - set(raw[0]):
        raise SystemExit(f"missing columns: {sorted(required - (set(raw[0]) if raw else set()))}")
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in raw:
        parsed: dict[str, object] = dict(row)
        for key in ("program_score", "core_hits", "support_hits", "anti_hits", "x", "y"):
            parsed[key] = float(row[key])
        grouped[(row["parent_label"], row["candidate_lineage"])].append(parsed)
    result_rows: list[dict[str, object]] = []
    for (parent_label, candidate), group in sorted(grouped.items()):
        threshold = quantile([float(row["program_score"]) for row in group], args.positive_quantile)
        positive = [row for row in group if float(row["program_score"]) >= threshold
                    and float(row["core_hits"]) >= 1 and float(row["support_hits"]) >= 1 and float(row["anti_hits"]) == 0]
        for component_id, indices in enumerate(spatial_components(positive, args.radius), start=1):
            members = [positive[index] for index in indices]
            fraction = len(members) / len(group)
            coherent = len(members) >= args.min_component
            embedded = coherent and fraction <= args.max_parent_fraction
            result_rows.append({
                "parent_label": parent_label, "candidate_lineage": candidate, "component_id": component_id,
                "n_component": len(members), "n_parent": len(group), "parent_fraction": fraction,
                "score_threshold": threshold,
                "mean_program_score": sum(float(row["program_score"]) for row in members) / len(members),
                "coherent": str(coherent).lower(), "embedded_candidate": str(embedded).lower(),
                "required_action": "targeted_recall" if embedded else "none",
            })
    args.out.mkdir(parents=True, exist_ok=True)
    dest = args.out / "embedded_lineage_components.tsv"
    fields = ["parent_label", "candidate_lineage", "component_id", "n_component", "n_parent", "parent_fraction",
              "score_threshold", "mean_program_score", "coherent", "embedded_candidate", "required_action"]
    with dest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t"); writer.writeheader(); writer.writerows(result_rows)
    n_candidates = sum(row["embedded_candidate"] == "true" for row in result_rows)
    manifest = {"status": "PASS", "n_rows": len(raw), "n_components": len(result_rows),
                "n_embedded_candidates": n_candidates, "policy": "observation_level_core_support_anti_plus_spatial_connectivity"}
    (args.out / "embedded_lineage_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
