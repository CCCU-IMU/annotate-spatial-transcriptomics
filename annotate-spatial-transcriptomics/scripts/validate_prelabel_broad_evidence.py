#!/usr/bin/env python3
"""Fail closed when initial broad decisions were not preceded by label-blind evidence."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from evidence_schema_lib import (
    active_registry_rows,
    path_at,
    read_tsv,
    sha256,
    validate_evidence_artifact,
    validate_json_against_schema,
    write_result,
)


def truth(value: object) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "pass", "passed"}


def validate(root: Path, ledger_path: Path) -> dict:
    rows = active_registry_rows(read_tsv(ledger_path), "decision_id")
    schema = Path(__file__).resolve().parents[1] / "schemas/prelabel_broad_evidence.schema.json"
    errors: list[str] = []
    validated = 0
    project_path = root / "config/project.json"
    try:
        framework = json.loads(project_path.read_text(encoding="utf-8")).get("framework_version", "0")
    except (OSError, json.JSONDecodeError):
        framework = "0"
    try:
        framework_tuple = tuple(int(part) for part in str(framework).split(".")[:3])
    except ValueError:
        framework_tuple = (0,)
    absolute_broad_contract = framework_tuple >= (1, 10, 0)
    for line, row in enumerate(rows, 2):
        cluster = row.get("source_cluster", "").strip()
        artifact_value = row.get("prelabel_evidence_artifact", "").strip()
        expected_hash = row.get("prelabel_evidence_sha256", "").strip()
        if not artifact_value:
            errors.append(f"cluster row {line} ({cluster}) lacks prelabel_evidence_artifact")
            continue
        artifact = path_at(root, artifact_value)
        document, document_errors = validate_json_against_schema(artifact, schema)
        errors.extend(f"cluster row {line} ({cluster}): {message}" for message in document_errors)
        if document_errors:
            continue
        if not expected_hash or sha256(artifact) != expected_hash:
            errors.append(f"cluster row {line} ({cluster}) has a stale prelabel evidence hash")
        if document["source_cluster"] != cluster:
            errors.append(f"cluster row {line} disagrees with prelabel artifact source_cluster")
        try:
            if int(float(row.get("n_observations", 0) or 0)) != document["n_observations"]:
                errors.append(f"cluster row {line} disagrees with prelabel artifact n_observations")
        except ValueError:
            errors.append(f"cluster row {line} has invalid n_observations")
        candidates = document["candidate_lineages"]
        hypotheses = document["lineage_hypotheses"]
        names = [item["candidate_id"] for item in hypotheses]
        if len(names) != len(set(names)) or set(names) != set(candidates):
            errors.append(f"cluster row {line} hypotheses do not exactly cover candidate_lineages")
        catalog_path, catalog_errors = validate_evidence_artifact(root, document["candidate_catalog_artifact"], f"cluster row {line} candidate catalog")
        errors.extend(catalog_errors)
        if catalog_path and not catalog_errors:
            try:
                catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                catalog = {}
            boundaries = catalog.get("candidate_boundaries", []) if isinstance(catalog, dict) else []
            required_entries = [
                item for item in boundaries
                if isinstance(item, dict) and item.get("review_required") is True
            ]
            required_ids = {str(item.get("candidate_id", "")).strip() for item in required_entries}
            required_ids.discard("")
            if not required_ids:
                errors.append(f"cluster row {line} candidate catalog has no required boundary reviews")
            if len(required_ids) != len(required_entries):
                errors.append(f"cluster row {line} candidate catalog has duplicate or empty required candidate IDs")
            if set(candidates) != required_ids:
                errors.append(f"cluster row {line} candidate set differs from the bound catalog")
        winner = document["winner"]
        runner_up = document["runner_up"]
        if winner == runner_up or winner not in candidates or runner_up not in candidates:
            errors.append(f"cluster row {line} winner/runner_up contract failed")
        hypothesis_by_name = {item["candidate_id"]: item for item in hypotheses}
        winning = hypothesis_by_name.get(winner, {})
        initial = row.get("initial_broad_label", "").strip()
        confidence = row.get("confidence", "").strip().lower()
        if initial:
            if initial != winning.get("candidate_broad_lineage", ""):
                errors.append(f"cluster row {line} initial broad label differs from frozen winner")
            if confidence not in {"moderate", "high"}:
                errors.append(f"cluster row {line} assigned broad label below moderate confidence")
            if document["winning_margin"] <= 0:
                errors.append(f"cluster row {line} assigned broad label without a positive winner margin")
            if int(winning.get("positive_marker_family_count", 0)) < 2:
                errors.append(f"cluster row {line} assigned broad label from fewer than two marker families")
            if absolute_broad_contract:
                absolute_families = winning.get("absolute_supported_families", [])
                if int(winning.get("absolute_family_support_count", 0)) < 2 or len(absolute_families) < 2:
                    errors.append(f"cluster row {line} assigned broad label without two absolute-expression-supported families")
                if int(winning.get("absolute_family_support_count", 0)) != len(absolute_families):
                    errors.append(f"cluster row {line} absolute family support count is inconsistent")
                if winning.get("comparative_score_role") != "comparative_not_absence_gate":
                    errors.append(f"cluster row {line} does not declare centered/DEG scores as comparative-only")
                absolute = winning.get("absolute_detection_evidence")
                if not isinstance(absolute, dict):
                    errors.append(f"cluster row {line} lacks absolute detection/pseudobulk evidence")
                else:
                    _, absolute_errors = validate_evidence_artifact(
                        root, absolute, f"cluster row {line} absolute broad-presence evidence"
                    )
                    errors.extend(absolute_errors)
            if int(winning.get("contradiction_count", 0)) > 0:
                errors.append(f"cluster row {line} assigned broad label despite unresolved contradictions")
        if row.get("prelabel_winner", "") != winner or row.get("prelabel_runner_up", "") != runner_up:
            errors.append(f"cluster row {line} winner/runner_up registry fields are stale")
        try:
            if abs(float(row.get("prelabel_winning_margin", "nan")) - float(document["winning_margin"])) > 1e-12:
                errors.append(f"cluster row {line} winning margin registry field is stale")
        except ValueError:
            errors.append(f"cluster row {line} has invalid prelabel_winning_margin")
        for label, evidence in (("positive DEG", document["positive_deg_artifact"]), ("anti DEG", document["anti_deg_artifact"])):
            _, evidence_errors = validate_evidence_artifact(root, evidence, f"cluster row {line} {label}")
            errors.extend(evidence_errors)
        for hypothesis in hypotheses:
            if hypothesis["positive_marker_family_count"] != len(hypothesis["positive_marker_families"]):
                errors.append(f"cluster row {line} marker-family count is inconsistent for {hypothesis['candidate_broad_lineage']}")
            for index, evidence in enumerate(hypothesis["evidence_artifacts"], 1):
                _, evidence_errors = validate_evidence_artifact(
                    root, evidence, f"cluster row {line} {hypothesis['candidate_id']} evidence {index}"
                )
                errors.extend(evidence_errors)
        if not truth(row.get("prelabel_evidence_frozen", "")):
            errors.append(f"cluster row {line} is not marked prelabel_evidence_frozen")
        validated += 1
    if not rows:
        errors.append("cluster decision ledger is missing or empty")
    return {
        "status": "PASS" if not errors else "FAIL",
        "ledger": str(ledger_path),
        "ledger_sha256": sha256(ledger_path) if ledger_path.is_file() else None,
        "active_cluster_decisions": len(rows),
        "validated_prelabel_evidence": validated,
        "errors": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--ledger", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    result = validate(root, args.ledger or root / "state/cluster_decision_ledger.tsv")
    write_result(args.out or root / "provenance/prelabel_broad_evidence_validation.json", result)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
