#!/usr/bin/env python3
"""Require complete, hash-bound support evidence for every final broad/fine label."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from evidence_schema_lib import active_registry_rows, path_at, read_tsv, sha256, validate_evidence_artifact, validate_json_against_schema, write_result


def validate(root: Path, registry_path: Path, cell_ledger_path: Path) -> dict:
    rows = active_registry_rows(read_tsv(registry_path), "support_id")
    cells = read_tsv(cell_ledger_path)
    errors: list[str] = []
    if not rows:
        errors.append("annotation support registry is missing or empty")
    if not cells:
        errors.append("final cell ledger is missing or empty")
    schema = Path(__file__).resolve().parents[1] / "schemas/annotation_support.schema.json"
    validated_broad: set[str] = set()
    validated_fine: set[tuple[str, str]] = set()
    for line, row in enumerate(rows, 2):
        if row.get("status") != "validated":
            continue
        artifact_value = row.get("support_artifact", "")
        artifact_hash = row.get("support_artifact_sha256", "")
        if not artifact_value:
            errors.append(f"support row {line} lacks support_artifact")
            continue
        artifact_path = path_at(root, artifact_value)
        support, support_errors = validate_json_against_schema(artifact_path, schema)
        errors.extend(f"support row {line}: {message}" for message in support_errors)
        if support_errors:
            continue
        if not artifact_hash or sha256(artifact_path) != artifact_hash:
            errors.append(f"support row {line} has a stale support artifact hash")
        for index, artifact in enumerate(support["evidence_artifacts"]):
            _, artifact_errors = validate_evidence_artifact(root, artifact, f"support row {line} evidence artifact {index + 1}")
            errors.extend(artifact_errors)
        _, full_feature_errors = validate_evidence_artifact(root, support["full_feature_evidence"]["artifact"], f"support row {line} full-feature evidence")
        errors.extend(full_feature_errors)
        for field in ("positive_marker_evidence", "anti_marker_evidence", "resolution_evidence", "spatial_evidence", "provenance_evidence"):
            if set(support[field]) <= {"status", "created_at", "completed_at"}:
                errors.append(f"support row {line} {field} contains status metadata but no evidence")
        comparisons = {
            "support_id": support["support_id"],
            "label_level": support["label_level"],
            "broad_label": support["broad_label"],
            "fine_label": support["fine_label"],
            "confidence": support["confidence"],
        }
        for key, expected in comparisons.items():
            if row.get(key, "") != str(expected):
                errors.append(f"support row {line} disagrees with support artifact: {key}")
        if support["label_level"] == "broad":
            if support["fine_label"]:
                errors.append(f"support row {line} broad support contains a fine label")
            validated_broad.add(support["broad_label"])
        else:
            if not support["fine_label"] or support["confidence"] != "high":
                errors.append(f"support row {line} fine support lacks a high-confidence fine label")
            validated_fine.add((support["broad_label"], support["fine_label"]))
        if any(not str(value).strip() for value in support["alternative_hypotheses"]):
            errors.append(f"support row {line} contains an empty alternative hypothesis")
    analysis = [row for row in cells if row.get("analysis_scope") == "analysis_set"]
    final_broad = {row.get("final_broad_label", "").strip() for row in analysis if row.get("final_broad_label", "").strip()}
    final_fine = {
        (row.get("final_broad_label", "").strip(), row.get("final_fine_label", "").strip())
        for row in analysis if row.get("final_fine_label", "").strip()
    }
    canonical = {"low", "moderate", "high"}
    for row in analysis:
        confidence = row.get("final_confidence", "").strip()
        if row.get("final_broad_label", "").strip() and confidence not in {"moderate", "high"}:
            errors.append(f"final broad label for {row.get('cell_id','')} lacks canonical moderate/high confidence")
        if row.get("final_fine_label", "").strip():
            if confidence != "high" or str(row.get("fine_anchor_eligible", "")).strip().lower() not in {"1", "true", "yes"}:
                errors.append(f"final fine label for {row.get('cell_id','')} is not high-confidence and fine-anchor eligible")
        elif confidence and confidence not in canonical:
            errors.append(f"noncanonical final confidence for {row.get('cell_id','')}: {confidence}")
    if final_broad != validated_broad:
        errors.append(f"final broad labels differ from validated broad supports: final={sorted(final_broad)} support={sorted(validated_broad)}")
    if final_fine != validated_fine:
        errors.append(f"final fine labels differ from validated fine supports: final={sorted(final_fine)} support={sorted(validated_fine)}")
    return {
        "status": "PASS" if not errors else "FAIL",
        "registry": str(registry_path),
        "registry_sha256": sha256(registry_path) if registry_path.is_file() else None,
        "cell_ledger": str(cell_ledger_path),
        "cell_ledger_sha256": sha256(cell_ledger_path) if cell_ledger_path.is_file() else None,
        "validated_broad_labels": sorted(validated_broad),
        "validated_fine_labels": [list(value) for value in sorted(validated_fine)],
        "errors": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--registry", type=Path)
    parser.add_argument("--cell-ledger", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    registry = args.registry or root / "state/annotation_support_registry.tsv"
    ledger = args.cell_ledger or (root / "state/cell_ledger.tsv.gz")
    if not ledger.exists():
        ledger = root / "state/cell_ledger.tsv"
    result = validate(root, registry, ledger)
    out = args.out or root / "provenance/annotation_support_validation.json"
    write_result(out, result)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
