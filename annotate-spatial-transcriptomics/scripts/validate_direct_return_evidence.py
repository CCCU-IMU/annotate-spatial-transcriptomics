#!/usr/bin/env python3
"""Validate a direct-return evidence bundle and its registry binding."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from evidence_schema_lib import membership_ids, read_tsv, sha256, validate_evidence_artifact, validate_json_against_schema, write_result


def validate(root: Path, evidence_path: Path, registry_path: Path | None = None) -> dict:
    schema = Path(__file__).resolve().parents[1] / "schemas/direct_return_evidence.schema.json"
    evidence, errors = validate_json_against_schema(evidence_path, schema)
    if errors:
        return {"status": "FAIL", "artifact": str(evidence_path), "errors": errors}
    members, membership_errors = membership_ids(root, evidence["membership"], "direct return membership")
    errors.extend(membership_errors)
    for index, artifact in enumerate(evidence["evidence_artifacts"]):
        _, artifact_errors = validate_evidence_artifact(root, artifact, f"evidence artifact {index + 1}")
        errors.extend(artifact_errors)
    _, full_feature_errors = validate_evidence_artifact(root, evidence["full_feature_scope"]["artifact"], "direct return full-feature evidence")
    errors.extend(full_feature_errors)
    for field in ("anti_marker_results", "spatial_evidence"):
        if set(evidence[field]) <= {"status", "created_at", "completed_at"}:
            errors.append(f"{field} contains status metadata but no biological evidence")
    if evidence["target_fine_label"] and evidence["confidence"] != "high":
        errors.append("fine direct return must have canonical high confidence")
    if any(not str(value).strip() for value in evidence["positive_marker_families"]):
        errors.append("positive marker families cannot contain empty entries")
    if any(not str(value).strip() for value in evidence["alternative_hypotheses"]):
        errors.append("alternative hypotheses cannot contain empty entries")
    if registry_path and registry_path.is_file():
        matches = [row for row in read_tsv(registry_path) if row.get("return_id") == evidence["return_id"]]
        if len(matches) != 1:
            errors.append("direct return evidence must match exactly one registry row")
        else:
            registry = matches[0]
            comparisons = {
                "source_cohort_id": evidence["source_cohort_id"],
                "source_cluster": evidence["source_subcluster"],
                "target_broad_label": evidence["target_broad_label"],
                "target_fine_label": evidence["target_fine_label"],
                "membership_sha256": evidence["membership"]["sha256"],
                "confidence": evidence["confidence"],
            }
            for key, expected in comparisons.items():
                if registry.get(key, "") != str(expected):
                    errors.append(f"direct return registry mismatch: {key}")
            if registry.get("evidence_artifact_sha256") and registry["evidence_artifact_sha256"] != sha256(evidence_path):
                errors.append("direct return evidence hash differs from registry")
    return {
        "status": "PASS" if not errors else "FAIL",
        "artifact": str(evidence_path),
        "artifact_sha256": sha256(evidence_path),
        "return_id": evidence["return_id"],
        "n_observations": len(members),
        "errors": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("evidence", type=Path)
    parser.add_argument("--registry", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = validate(args.project_root.resolve(), args.evidence.resolve(), args.registry)
    write_result(args.out, result)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
