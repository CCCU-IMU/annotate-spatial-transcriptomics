#!/usr/bin/env python3
"""Validate a reusable low-cost Atlas representation and forbid costly defaults."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from evidence_schema_lib import validate_artifact_ref, validate_json_against_schema, write_result


def validate(root: Path, manifest_path: Path) -> dict:
    schema = Path(__file__).resolve().parents[1] / "schemas/atlas_index_manifest.schema.json"
    document, errors = validate_json_against_schema(manifest_path, schema)
    if not errors:
        for field in ("feature_transform", "reference_representation", "broad_crosswalk"):
            _, field_errors = validate_artifact_ref(root, document[field], f"Atlas index {field}")
            errors.extend(field_errors)
        engine = document["mapping_engine"]
        if engine == "fixed_projection_ann":
            if "neighbor_index" not in document:
                errors.append("fixed_projection_ann requires neighbor_index")
            else:
                _, field_errors = validate_artifact_ref(root, document["neighbor_index"], "Atlas neighbor index")
                errors.extend(field_errors)
        if engine == "diagnostic_exact_knn_small_reference" and document["n_reference"] > 50000:
            errors.append("exact kNN is restricted to references with at most 50000 observations")
        if document["n_dimensions"] > 256:
            errors.append("Atlas mapping dimensions exceed the efficient default ceiling of 256")
    return {
        "status": "PASS" if not errors else "FAIL",
        "index_manifest": str(manifest_path.resolve()),
        "reference_id": document.get("reference_id", ""),
        "mapping_engine": document.get("mapping_engine", ""),
        "reusable_across_queries": document.get("reusable_across_queries", False),
        "dense_pairwise_matrix": document.get("dense_pairwise_matrix", True),
        "query_reference_joint_retraining": document.get("query_reference_joint_retraining", True),
        "whole_object_rctd": document.get("whole_object_rctd", True),
        "errors": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = validate(args.project_root.resolve(), args.manifest)
    write_result(args.out, result)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
