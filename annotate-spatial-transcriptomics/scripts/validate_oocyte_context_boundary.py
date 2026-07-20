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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--canonical-membership", required=True, type=Path)
    parser.add_argument("--context-membership", required=True, type=Path)
    parser.add_argument("--outcomes", required=True, type=Path, help="TSV: cell_id, final_broad_label")
    parser.add_argument("--cell-id-column", default="cell_id")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    canonical = ids(args.canonical_membership, args.cell_id_column)
    context = ids(args.context_membership, args.cell_id_column)
    errors: list[str] = []
    if not canonical.issubset(context):
        errors.append("spatial context must contain the complete canonical Oocyte recall set")
    context_only_oocyte: list[str] = []
    with args.outcomes.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            cell = str(row.get(args.cell_id_column, ""))
            label = row.get("final_broad_label", "").lower()
            if cell in context - canonical and label in {"oocyte", "germ cell", "germ_cell"}:
                context_only_oocyte.append(cell)
    if context_only_oocyte:
        errors.append(f"{len(context_only_oocyte)} context-only observations were written into the Oocyte census")
    result = {"status": "PASS" if not errors else "BLOCKED", "canonical_n": len(canonical),
              "context_n": len(context), "context_only_n": len(context - canonical),
              "context_only_oocyte_n": len(context_only_oocyte), "errors": errors}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
