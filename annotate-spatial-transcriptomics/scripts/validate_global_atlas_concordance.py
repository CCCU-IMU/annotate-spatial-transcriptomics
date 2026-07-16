#!/usr/bin/env python3
"""Validate exact all-cell Atlas coverage and closure of every triggered review."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

from evidence_schema_lib import validate_evidence_artifact, validate_json_against_schema


ALLOWED_OUTCOMES = {"retain_primary", "supersede_broad", "downgrade_to_broad", "return_to_qc", "unknown_candidate"}


def open_text(path: Path, mode: str):
    return gzip.open(path, mode, newline="", encoding="utf-8") if path.suffix == ".gz" else path.open(mode, newline="", encoding="utf-8")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with open_text(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit-root", required=True, type=Path)
    parser.add_argument("--decisions", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    manifest_path = args.audit_root / "global_atlas_concordance_manifest.json"
    errors: list[str] = []
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"cannot load global Atlas manifest: {exc}")
    artifacts = manifest.get("artifacts", {})
    for name, record in artifacts.items():
        path = Path(record.get("path", ""))
        if not path.is_file() or not record.get("sha256") or sha256(path) != record["sha256"]:
            errors.append(f"artifact is missing or stale: {name}")
    all_cells = read_tsv(Path(artifacts.get("all_cell_concordance", {}).get("path", "")))
    ids = [row.get("cell_id", "") for row in all_cells]
    if len(ids) != int(manifest.get("n_analysis_set", -1)) or not ids or "" in ids or len(ids) != len(set(ids)):
        errors.append("all-cell concordance does not exactly contain unique analysis-set cells")
    if any(row.get("comparison_status") == "qc_writeback_candidate" and row.get("atlas_tier") not in {"high", "moderate_only"} for row in all_cells):
        errors.append("QC writeback contains a below-moderate Atlas result")
    if any(row.get("comparison_status") == "qc_writeback_candidate" and row.get("consensus_pass") != "true" for row in all_cells):
        errors.append("QC writeback bypasses multichannel consensus")
    if any(row.get("comparison_status") == "qc_writeback_candidate" and (row.get("out_of_distribution") == "true" or row.get("ontology_conflict") == "true") for row in all_cells):
        errors.append("QC writeback contains OOD or ontology-conflicted observations")
    queue = read_tsv(Path(artifacts.get("discrepancy_review_queue", {}).get("path", "")))
    queue_ids = {row.get("review_id", "") for row in queue}
    if "" in queue_ids or len(queue_ids) != len(queue):
        errors.append("discrepancy review queue IDs are empty or duplicated")
    decision_rows = read_tsv(args.decisions) if args.decisions else []
    decision_ids = [row.get("review_id", "") for row in decision_rows]
    if queue_ids:
        if set(decision_ids) != queue_ids or len(decision_ids) != len(set(decision_ids)):
            errors.append("discrepancy decisions do not close every queued review exactly once")
        for line, row in enumerate(decision_rows, 2):
            if row.get("outcome") not in ALLOWED_OUTCOMES:
                errors.append(f"decision row {line} has invalid outcome")
            if len(row.get("rationale", "").strip()) < 20:
                errors.append(f"decision row {line} rationale is too short")
            evidence = Path(row.get("evidence_artifact", ""))
            if not evidence.is_absolute() and args.decisions:
                evidence = args.decisions.parent / evidence
            if not evidence.is_file() or not row.get("evidence_sha256") or sha256(evidence) != row.get("evidence_sha256"):
                errors.append(f"decision row {line} evidence artifact is missing or stale")
                continue
            schema = Path(__file__).resolve().parents[1] / "schemas/atlas_discrepancy_decision.schema.json"
            decision, decision_errors = validate_json_against_schema(evidence, schema)
            errors.extend(f"decision row {line}: {message}" for message in decision_errors)
            if decision_errors:
                continue
            if decision.get("review_id") != row.get("review_id") or decision.get("outcome") != row.get("outcome"):
                errors.append(f"decision row {line} disagrees with its orthogonal evidence artifact")
            for index, artifact in enumerate(decision.get("evidence_artifacts", []), 1):
                _, artifact_errors = validate_evidence_artifact(evidence.parent, artifact, f"decision row {line} query evidence {index}")
                errors.extend(artifact_errors)
            if decision.get("atlas_only") is not False or row.get("atlas_only", "").strip().lower() in {"1", "true", "yes"}:
                errors.append(f"decision row {line} attempts an Atlas-only biological overwrite")
    elif decision_rows:
        errors.append("discrepancy decisions exist although the review queue is empty")
    result = {
        "status": "PASS" if not errors else "FAIL",
        "audit_manifest": str(manifest_path.resolve()),
        "audit_manifest_sha256": sha256(manifest_path),
        "n_analysis_set": len(all_cells),
        "review_queue_n": len(queue),
        "review_decisions_n": len(decision_rows),
        "decisions": str(args.decisions.resolve()) if args.decisions else "",
        "decisions_sha256": sha256(args.decisions) if args.decisions and args.decisions.is_file() else "",
        "fine_anchor_eligible": False,
        "errors": errors,
    }
    output = args.out or args.audit_root / "global_atlas_concordance_validation.json"
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
