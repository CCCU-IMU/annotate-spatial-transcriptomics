#!/usr/bin/env python3
"""Validate a source-label to candidate-taxonomy crosswalk without assigning cells."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


REQUIRED = (
    "source_label",
    "canonical_review_axis",
    "candidate_release_broad",
    "candidate_fine_label",
    "transfer_ceiling",
    "positive_markers_in_reference",
    "required_query_evidence",
    "rationale",
)

ALLOWED_CEILINGS = {
    "candidate_only",
    "broad_only_after_calibration",
    "strict_context_gate",
    "fine_candidate_after_query_validation",
    "retain_source_only",
}

FORBIDDEN_BROAD_SUFFIXES = ("_review", "_candidate", "_unresolved", "_holdout")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("crosswalk", type=Path)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    errors: list[str] = []
    warnings: list[str] = []
    rows: list[dict[str, str]] = []
    with args.crosswalk.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [name for name in REQUIRED if name not in (reader.fieldnames or [])]
        if missing:
            errors.append("missing columns: " + ", ".join(missing))
        else:
            rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]

    seen: set[str] = set()
    for line_no, row in enumerate(rows, start=2):
        source = row["source_label"]
        if not source:
            errors.append(f"line {line_no}: source_label is empty")
        elif source in seen:
            errors.append(f"line {line_no}: duplicate source_label {source!r}")
        seen.add(source)

        ceiling = row["transfer_ceiling"]
        if ceiling not in ALLOWED_CEILINGS:
            errors.append(f"line {line_no}: unsupported transfer_ceiling {ceiling!r}")
        broad = row["candidate_release_broad"]
        broad_norm = broad.lower().replace(" ", "_")
        if any(broad_norm.endswith(suffix) for suffix in FORBIDDEN_BROAD_SUFFIXES):
            errors.append(f"line {line_no}: analysis-pool name copied into candidate_release_broad")
        if not broad and ceiling not in {"candidate_only", "retain_source_only"}:
            errors.append(f"line {line_no}: blank candidate_release_broad requires candidate_only or retain_source_only")
        if row["candidate_fine_label"] and ceiling == "broad_only_after_calibration":
            warnings.append(f"line {line_no}: fine label is provenance-only under a broad-only ceiling")
        if not row["positive_markers_in_reference"]:
            warnings.append(f"line {line_no}: no positive reference markers recorded")
        if not row["required_query_evidence"]:
            errors.append(f"line {line_no}: required_query_evidence is empty")
        if source == broad and broad:
            warnings.append(f"line {line_no}: source and candidate broad names match; still preserve the crosswalk decision")

    result = {
        "status": "PASS" if not errors else "FAIL",
        "crosswalk": str(args.crosswalk.resolve()),
        "n_source_labels": len(rows),
        "transfer_ceiling_counts": {
            ceiling: sum(row.get("transfer_ceiling") == ceiling for row in rows)
            for ceiling in sorted(ALLOWED_CEILINGS)
        },
        "errors": errors,
        "warnings": warnings,
    }
    out = args.out or args.crosswalk.with_suffix(".validation.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
