#!/usr/bin/env python3
"""Fail closed on inconsistent or taxonomy-invalid final fine labels."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path


TRUE = {"1", "true", "yes", "high"}


def opener(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(encoding="utf-8", newline="")


def catalog_map(path: Path) -> dict[str, str]:
    doc = json.loads(path.read_text(encoding="utf-8"))
    result: dict[str, str] = {}
    for items in doc.get("machine_actionable_fine_candidate_catalog", {}).values():
        for item in items if isinstance(items, list) else []:
            if isinstance(item, dict) and item.get("release_label"):
                result[str(item["release_label"])] = str(item.get("parent_release_label", ""))
    return result


def validate(ledger: Path, catalog: Path) -> dict[str, object]:
    labels = catalog_map(catalog)
    errors: list[str] = []
    fine_census: Counter[str] = Counter()
    with opener(ledger) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"cell_id", "final_state", "final_broad_label", "final_fine_label", "final_fine_confidence", "final_fine_eligible", "fine_anchor_eligible"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            errors.append("ledger lacks required final fine-state columns")
        else:
            for row in reader:
                fine = str(row.get("final_fine_label", "")).strip()
                if not fine:
                    if row.get("final_state") == "defined_fine" or str(row.get("final_fine_eligible", "")).lower() in TRUE:
                        errors.append(f"{row['cell_id']}: empty fine label has a fine terminal state/eligibility")
                    continue
                fine_census[fine] += 1
                parent = labels.get(fine)
                if not parent:
                    errors.append(f"{row['cell_id']}: fine label is absent from the bound catalog: {fine}")
                elif parent != str(row.get("final_broad_label", "")).strip():
                    errors.append(f"{row['cell_id']}: fine label crosses parent: {row.get('final_broad_label')} -> {fine}")
                if row.get("final_state") != "defined_fine":
                    errors.append(f"{row['cell_id']}: nonempty fine label is not defined_fine")
                if str(row.get("final_fine_confidence", "")).lower() != "high":
                    errors.append(f"{row['cell_id']}: released fine label is not high-confidence")
                if str(row.get("final_fine_eligible", "")).lower() not in TRUE or str(row.get("fine_anchor_eligible", "")).lower() not in TRUE:
                    errors.append(f"{row['cell_id']}: released fine label lacks fine eligibility")
                if len(errors) >= 100:
                    errors.append("error list truncated at 100 rows")
                    break
    return {"status": "PASS" if not errors else "BLOCKED", "fine_census": dict(sorted(fine_census.items())), "errors": errors}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("ledger", type=Path)
    parser.add_argument("--catalog", required=True, type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = validate(args.ledger, args.catalog)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
