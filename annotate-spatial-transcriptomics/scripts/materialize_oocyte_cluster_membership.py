#!/usr/bin/env python3
"""Materialize Oocyte membership from passing targeted-recluster clusters.

Strict seeds and spatial-object IDs are intentionally ignored here.  They are
cluster-adjudication evidence, not observation-level inclusion filters.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--canonical-membership", required=True, type=Path)
    ap.add_argument("--passing-clusters", required=True, type=Path)
    ap.add_argument("--explicit-exclusions", type=Path)
    ap.add_argument("--cell-id-column", default="cell_id")
    ap.add_argument("--cluster-column", default="recluster_cluster")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    canonical = read_rows(args.canonical_membership)
    passing_rows = read_rows(args.passing_clusters)
    if not canonical:
        raise SystemExit("canonical Oocyte membership is empty")
    required = {args.cell_id_column, args.cluster_column}
    missing = required.difference(canonical[0])
    if missing:
        raise SystemExit(f"canonical membership lacks columns: {sorted(missing)}")
    passing = {
        str(row.get(args.cluster_column, "")).strip()
        for row in passing_rows
        if str(row.get(args.cluster_column, "")).strip()
        and str(row.get("adjudication_status", "pass")).strip().lower()
        in {"pass", "passed", "supported", "oocyte"}
    }
    if not passing:
        raise SystemExit("no passing Oocyte recluster cluster was supplied")

    exclusions: dict[str, str] = {}
    if args.explicit_exclusions:
        for row in read_rows(args.explicit_exclusions):
            cell = str(row.get(args.cell_id_column, "")).strip()
            reason = str(row.get("exclusion_reason", "")).strip()
            if not cell or not reason:
                raise SystemExit("explicit exclusions require cell ID and exclusion_reason")
            exclusions[cell] = reason

    seen: set[str] = set()
    materialized: list[dict[str, str]] = []
    eligible_n = 0
    for row in canonical:
        cell = str(row[args.cell_id_column]).strip()
        cluster = str(row[args.cluster_column]).strip()
        if not cell or cell in seen:
            raise SystemExit("canonical membership cell IDs must be nonempty and unique")
        seen.add(cell)
        if cluster not in passing:
            continue
        eligible_n += 1
        if cell in exclusions:
            continue
        materialized.append(
            {
                args.cell_id_column: cell,
                args.cluster_column: cluster,
                "final_broad_label": "Oocyte",
                "decision_basis": "canonical_member_of_passing_oocyte_recluster_cluster",
                "fine_anchor_eligible": "false",
            }
        )

    args.out.mkdir(parents=True, exist_ok=True)
    table = args.out / "materialized_oocyte_membership.tsv"
    fields = [
        args.cell_id_column,
        args.cluster_column,
        "final_broad_label",
        "decision_basis",
        "fine_anchor_eligible",
    ]
    with table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(materialized)
    manifest = {
        "status": "PASS",
        "canonical_n": len(canonical),
        "passing_clusters": sorted(passing),
        "eligible_canonical_members_n": eligible_n,
        "explicit_exclusion_n": sum(cell in exclusions for cell in seen),
        "final_oocyte_n": len(materialized),
        "strict_seed_used_as_membership_filter": False,
        "spatial_object_used_as_membership_filter": False,
    }
    (args.out / "materialized_oocyte_membership_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
