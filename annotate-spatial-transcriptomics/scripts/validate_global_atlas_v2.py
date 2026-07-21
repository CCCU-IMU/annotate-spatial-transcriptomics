#!/usr/bin/env python3
"""Validate v2 all-cell Atlas routing and closure of every material review."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from evidence_schema_lib import read_tsv, sha256, validate_artifact_ref, validate_evidence_artifact, validate_json_against_schema


def truth(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "pass", "passed"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--routing-manifest", required=True, type=Path)
    ap.add_argument("--decisions", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    manifest = json.loads(args.routing_manifest.read_text(encoding="utf-8"))
    errors: list[str] = []
    if manifest.get("schema_version") != "2.0" or manifest.get("authoritative_router") != "route_global_atlas_v2.py":
        errors.append("routing manifest is not authoritative v2 output")
    if manifest.get("ledger_writeback_performed") is not False or manifest.get("fine_anchor_eligible") is not False:
        errors.append("Atlas routing must remain a broad-only atomic-commit proposal")
    resolved: dict[str, Path] = {}
    for key in ("cell_ledger", "atlas_mapping", "calibration_manifest", "workflow_profile"):
        path, artifact_errors = validate_artifact_ref(args.routing_manifest.parent, manifest.get(key, {}), key)
        errors.extend(artifact_errors)
        if path:
            resolved[key] = path
    for key, record in manifest.get("artifacts", {}).items():
        path, artifact_errors = validate_artifact_ref(args.routing_manifest.parent, record, key)
        errors.extend(artifact_errors)
        if path:
            resolved[key] = path
    routes = read_tsv(resolved.get("routing", Path())) if "routing" in resolved else []
    if len(routes) != int(manifest.get("n_analysis_set", -1)):
        errors.append("routing rows do not equal n_analysis_set")
    ids = [row.get("cell_id", "") for row in routes]
    if not ids or "" in ids or len(ids) != len(set(ids)):
        errors.append("routing cell IDs are empty or duplicated")
    for line, row in enumerate(routes, 2):
        route = row.get("atlas_state_route", "")
        if row.get("writeback_status") != "proposal_only_requires_atomic_commit" or row.get("fine_anchor_eligible") != "false":
            errors.append(f"route row {line} violates proposal/fine-anchor semantics")
        if route == "direct_qc_broad_return":
            if row.get("primary_broad") or row.get("primary_state") not in {"qc_holdout", "low_information_qc_holdout", "pending_qc"}:
                errors.append(f"route row {line} writes outside unlabeled frozen QC")
            if row.get("atlas_tier") not in {"high", "moderate_only"} or row.get("atlas_class_calibrated") != "true" or row.get("atlas_scope_pass") != "true":
                errors.append(f"route row {line} bypasses calibrated class/scope gates")
            if truth(row.get("out_of_distribution", "")) or truth(row.get("ontology_conflict", "")):
                errors.append(f"route row {line} writes an OOD/ontology-conflicted result")
        if row.get("primary_broad") and row.get("proposed_broad_label") != row.get("primary_broad"):
            errors.append(f"route row {line} silently overwrites a defined primary label")
    queue = read_tsv(resolved.get("review_queue", Path())) if "review_queue" in resolved else []
    queue_ids = {row.get("review_id", "") for row in queue}
    if "" in queue_ids or len(queue_ids) != len(queue):
        errors.append("review queue IDs are empty or duplicated")
    routed_review_ids = [row.get("review_id", "") for row in routes if truth(row.get("review_required", ""))]
    if any(not value for value in routed_review_ids):
        errors.append("a review-required routing row lacks review_id")
    if set(routed_review_ids) != queue_ids:
        errors.append("routing review IDs do not equal the material review queue")
    if any(row.get("review_id") and not truth(row.get("review_required", "")) for row in routes):
        errors.append("a routing row carries review_id without review_required=true")
    expected_manifest_status = "REVIEW_REQUIRED" if queue_ids else "PASS_NO_REVIEW"
    if manifest.get("status") != expected_manifest_status:
        errors.append("routing manifest status disagrees with the material review queue")
    decisions = read_tsv(args.decisions) if args.decisions else []
    if queue_ids:
        decision_ids = [row.get("review_id", "") for row in decisions]
        if set(decision_ids) != queue_ids or len(decision_ids) != len(set(decision_ids)):
            errors.append("review decisions do not close the queue exactly once")
        schema = Path(__file__).resolve().parents[1] / "schemas/atlas_discrepancy_decision.schema.json"
        for line, row in enumerate(decisions, 2):
            evidence = Path(row.get("evidence_artifact", ""))
            if not evidence.is_absolute() and args.decisions:
                evidence = args.decisions.parent / evidence
            document, document_errors = validate_json_against_schema(evidence, schema)
            errors.extend(f"decision row {line}: {message}" for message in document_errors)
            if document_errors:
                continue
            if document.get("review_id") != row.get("review_id") or document.get("atlas_only") is not False:
                errors.append(f"decision row {line} is stale or Atlas-only")
            for index, artifact in enumerate(document.get("evidence_artifacts", []), 1):
                _, artifact_errors = validate_evidence_artifact(evidence.parent, artifact, f"decision row {line} evidence {index}")
                errors.extend(artifact_errors)
    elif decisions:
        errors.append("decisions exist although the review queue is empty")
    result = {
        "status": "PASS" if not errors else "BLOCKED",
        "schema_version": "2.0",
        "routing_manifest": str(args.routing_manifest.resolve()),
        "routing_manifest_sha256": sha256(args.routing_manifest),
        "n_analysis_set": len(routes),
        "review_queue_n": len(queue),
        "review_decisions_n": len(decisions),
        "errors": errors,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
