#!/usr/bin/env python3
"""Keep canonical Oocyte recall membership distinct from spatial context discovery."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def ids(path: Path, column: str) -> set[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {str(row[column]) for row in csv.DictReader(handle, delimiter="\t")}


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--canonical-membership", required=True, type=Path)
    parser.add_argument("--context-membership", required=True, type=Path)
    parser.add_argument("--outcomes", required=True, type=Path, help="TSV: cell_id, final_broad_label")
    parser.add_argument("--cell-id-column", default="cell_id")
    parser.add_argument("--cluster-column", default="recluster_cluster")
    parser.add_argument("--passing-clusters", type=Path,
                        help="Optional TSV binding passing Oocyte clusters to the complete canonical membership")
    parser.add_argument("--explicit-exclusions", type=Path,
                        help="Optional TSV of cell_id,exclusion_reason hard exclusions")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    canonical = ids(args.canonical_membership, args.cell_id_column)
    context = ids(args.context_membership, args.cell_id_column)
    errors: list[str] = []
    if not canonical.issubset(context):
        errors.append("spatial context must contain the complete canonical Oocyte recall set")
    outcome_rows = rows(args.outcomes)
    oocyte_outcomes = {
        str(row.get(args.cell_id_column, ""))
        for row in outcome_rows
        if row.get("final_broad_label", "").lower() in {"oocyte", "germ cell", "germ_cell"}
    }
    context_only_oocyte = sorted((context - canonical).intersection(oocyte_outcomes))
    if context_only_oocyte:
        errors.append(f"{len(context_only_oocyte)} context-only observations were written into the Oocyte census")
    expected_cluster_members: set[str] = set()
    missing_cluster_members: list[str] = []
    if args.passing_clusters:
        canonical_rows = rows(args.canonical_membership)
        if canonical_rows and args.cluster_column not in canonical_rows[0]:
            errors.append(f"canonical membership lacks cluster column: {args.cluster_column}")
        else:
            passing = {
                str(row.get(args.cluster_column, "")).strip()
                for row in rows(args.passing_clusters)
                if str(row.get(args.cluster_column, "")).strip()
                and str(row.get("adjudication_status", "pass")).strip().lower()
                in {"pass", "passed", "supported", "oocyte"}
            }
            exclusions = ids(args.explicit_exclusions, args.cell_id_column) if args.explicit_exclusions else set()
            expected_cluster_members = {
                str(row.get(args.cell_id_column, ""))
                for row in canonical_rows
                if str(row.get(args.cluster_column, "")).strip() in passing
                and str(row.get(args.cell_id_column, "")) not in exclusions
            }
            missing_cluster_members = sorted(expected_cluster_members.difference(oocyte_outcomes))
            if missing_cluster_members:
                errors.append(
                    f"{len(missing_cluster_members)} canonical members of passing Oocyte clusters were omitted; "
                    "strict-seed/object membership cannot filter the final census"
                )
    result = {"status": "PASS" if not errors else "BLOCKED", "canonical_n": len(canonical),
              "context_n": len(context), "context_only_n": len(context - canonical),
              "context_only_oocyte_n": len(context_only_oocyte),
              "expected_passing_cluster_member_n": len(expected_cluster_members),
              "missing_passing_cluster_member_n": len(missing_cluster_members), "errors": errors}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
