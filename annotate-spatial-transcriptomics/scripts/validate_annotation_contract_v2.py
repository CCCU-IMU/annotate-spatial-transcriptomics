#!/usr/bin/env python3
"""Fail closed when a v2 annotation contract or any bound artifact is stale."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from evidence_schema_lib import sha256, validate_artifact_ref, validate_json_against_schema


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("contract", type=Path)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    schema = Path(__file__).resolve().parents[1] / "schemas/annotation_contract.schema.json"
    contract, errors = validate_json_against_schema(args.contract, schema)
    root = args.contract.parent
    resolved = {}
    for key in ("project_config", "workflow_profile", "biological_profile", "candidate_catalog"):
        path, artifact_errors = validate_artifact_ref(root, contract.get(key, {}), key)
        errors.extend(artifact_errors)
        if path:
            resolved[key] = path
    if "project_config" in resolved:
        project = json.loads(resolved["project_config"].read_text(encoding="utf-8"))
        for key in ("project_id", "sample_id", "framework_version"):
            if project.get(key) != contract.get(key):
                errors.append(f"project config disagrees with contract: {key}")
    snapshot_registry = Path(contract.get("selected_input_snapshot", {}).get("registry_path", ""))
    if snapshot_registry.is_file():
        with snapshot_registry.open(newline="", encoding="utf-8") as handle:
            rows = [row for row in csv.DictReader(handle, delimiter="\t") if row.get("snapshot_id") == contract.get("selected_input_snapshot", {}).get("snapshot_id")]
        if len(rows) != 1 or rows[0].get("sha256") != contract.get("selected_input_snapshot", {}).get("sha256"):
            errors.append("selected input snapshot is missing, duplicated or changed")
    else:
        errors.append("selected input snapshot registry is missing")
    whole = contract.get("whole_tissue_partition", {})
    query = contract.get("query_reclustering", {})
    if whole.get("method") == "BANKSY" and whole.get("candidate_grid_source") not in {"bound_upstream_input", "fresh_project_computation"}:
        errors.append("BANKSY whole-tissue selection requires a declared fresh or bound-upstream grid source")
    if whole.get("candidate_grid_source") == "bound_upstream_input" or whole.get("method") == "BANKSY":
        grid_path, grid_errors = validate_artifact_ref(root, whole.get("grid_artifact", {}), "whole-tissue grid artifact")
        errors.extend(grid_errors)
        if grid_path:
            try:
                grid_payload = json.loads(grid_path.read_text(encoding="utf-8"))
                grid = grid_payload.get("candidate_resolutions", grid_payload.get("resolutions", []))
                if sorted({float(value) for value in grid}) != whole.get("candidate_resolutions"):
                    errors.append("whole-tissue resolutions differ from the bound upstream grid artifact")
            except (OSError, UnicodeDecodeError, json.JSONDecodeError, TypeError, ValueError):
                errors.append("whole-tissue grid artifact is not valid resolution JSON")
    if whole.get("candidate_grid_source") == "fresh_project_computation" and "biological_profile" in resolved:
        biological = json.loads(resolved["biological_profile"].read_text(encoding="utf-8"))
        fresh = biological.get("resolution_policy", {}).get("fresh_sct_banksy_whole_tissue_candidate_resolutions", [])
        if fresh and sorted({float(value) for value in fresh}) != whole.get("candidate_resolutions"):
            errors.append("fresh whole-tissue resolutions differ from the bound biological profile")
    if whole.get("method") == "BANKSY" and whole.get("candidate_resolutions") == query.get("candidate_resolutions"):
        # Equality is possible, but it must not arise by silently substituting
        # the query grid. The bound-upstream source above is the decisive gate.
        pass
    result = {
        "status": "PASS" if not errors else "BLOCKED",
        "schema_version": "2.0",
        "contract": str(args.contract.resolve()),
        "contract_sha256": sha256(args.contract),
        "errors": errors,
    }
    out = args.out or args.contract.parent.parent / "provenance/annotation_contract_validation.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
