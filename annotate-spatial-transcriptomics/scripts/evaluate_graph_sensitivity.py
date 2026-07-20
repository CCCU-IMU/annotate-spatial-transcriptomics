#!/usr/bin/env python3
"""Decide whether a bounded graph-sensitivity audit rescued biological separation."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def number(row: dict[str, str], key: str) -> float:
    try:
        return float(row.get(key, "0") or 0)
    except ValueError:
        return 0.0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("summary", type=Path, help="TSV with graph_id, role, k, n_clusters and biological metrics")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--min-separation-gain", type=float, default=0.10)
    args = parser.parse_args()
    with args.summary.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    default = next((x for x in rows if x.get("role") == "default"), None)
    local = next((x for x in rows if x.get("role") == "local"), None)
    if default is None or local is None:
        result = {"status": "BLOCKED", "rescued": False, "reason": "default and local graph rows are required"}
    else:
        cluster_gain = number(local, "n_clusters") - number(default, "n_clusters")
        separation_gain = number(local, "core_support_separation") - number(default, "core_support_separation")
        anti_gain = number(local, "anti_clearance") - number(default, "anti_clearance")
        spatial_gain = number(local, "spatial_coherence") - number(default, "spatial_coherence")
        rescued = separation_gain >= args.min_separation_gain and anti_gain > 0 and spatial_gain >= 0
        result = {
            "status": "PASS", "rescued": rescued, "cluster_count_increase": cluster_gain,
            "core_support_separation_gain": separation_gain, "anti_clearance_gain": anti_gain,
            "spatial_coherence_gain": spatial_gain,
            "decision": "biological_rescue" if rescued else "not_rescued_do_not_name_small_clusters",
            "reason": "more clusters alone are insufficient" if cluster_gain > 0 and not rescued else "multichannel comparison completed",
        }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
