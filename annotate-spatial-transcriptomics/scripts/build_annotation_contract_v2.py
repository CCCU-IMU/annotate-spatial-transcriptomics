#!/usr/bin/env python3
"""Freeze the v2 project/profile/input/resolution contract before annotation."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from evidence_schema_lib import sha256


def artifact(path: Path) -> dict[str, str]:
    if not path.is_file():
        raise SystemExit(f"missing contract artifact: {path}")
    return {"path": str(path.resolve()), "sha256": sha256(path)}


def resolutions(value: str) -> list[float]:
    try:
        result = sorted({float(x) for x in value.replace(";", ",").split(",") if x.strip()})
    except ValueError as exc:
        raise SystemExit(f"invalid resolution grid: {exc}")
    if len(result) < 3:
        raise SystemExit("v2 whole-tissue contract requires at least three candidate resolutions")
    return result


def grid_artifact(path: Path, expected: list[float]) -> dict[str, str]:
    record = artifact(path)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        observed = payload.get("candidate_resolutions", payload.get("resolutions", []))
        observed = sorted({float(value) for value in observed})
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, TypeError, ValueError):
        raise SystemExit("--whole-tissue-grid-artifact must be JSON with candidate_resolutions/resolutions")
    if observed != expected:
        raise SystemExit("whole-tissue grid differs from the bound upstream grid artifact")
    return record


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("project_root", type=Path)
    ap.add_argument("--workflow-profile", required=True, type=Path)
    ap.add_argument("--biological-profile", required=True, type=Path)
    ap.add_argument("--candidate-catalog", required=True, type=Path)
    ap.add_argument("--snapshot-id", required=True)
    ap.add_argument("--whole-tissue-method", choices=["BANKSY", "Seurat", "Scanpy", "external"], required=True)
    ap.add_argument("--whole-tissue-grid", required=True)
    ap.add_argument("--whole-tissue-grid-artifact", type=Path)
    ap.add_argument("--grid-source", choices=["bound_upstream_input", "fresh_project_computation"], required=True)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    root = args.project_root.resolve()
    project_path = root / "config/project.json"
    project = json.loads(project_path.read_text(encoding="utf-8"))
    if project.get("framework_version") != "2.0.0":
        raise SystemExit("project framework_version must be 2.0.0")
    snapshot_registry = root / "state/input_snapshot_registry.tsv"
    with snapshot_registry.open(newline="", encoding="utf-8") as handle:
        matches = [row for row in csv.DictReader(handle, delimiter="\t") if row.get("snapshot_id") == args.snapshot_id and row.get("status") in {"frozen", "validated", "active"}]
    if len(matches) != 1 or len(matches[0].get("sha256", "")) != 64:
        raise SystemExit("selected snapshot is not one unique frozen/validated registry row")
    biological = json.loads(args.biological_profile.read_text(encoding="utf-8"))
    whole_grid = resolutions(args.whole_tissue_grid)
    query_grid = biological.get("resolution_policy", {}).get("query_reclustering_candidate_resolutions", [])
    if len(query_grid) < 3:
        raise SystemExit("biological profile lacks the v2 query-reclustering grid")
    if args.grid_source == "bound_upstream_input" and not args.whole_tissue_grid_artifact:
        raise SystemExit("bound_upstream_input requires --whole-tissue-grid-artifact")
    fresh_grid = biological.get("resolution_policy", {}).get("fresh_r_first_whole_tissue_candidate_resolutions", [])
    if args.grid_source == "fresh_project_computation" and fresh_grid and whole_grid != sorted({float(value) for value in fresh_grid}):
        raise SystemExit("fresh whole-tissue grid differs from the bound biological profile")
    contract = {
        "schema_version": "2.0",
        "framework_version": "2.0.0",
        "project_id": project["project_id"],
        "sample_id": project["sample_id"],
        "project_config": artifact(project_path),
        "workflow_profile": artifact(args.workflow_profile),
        "biological_profile": artifact(args.biological_profile),
        "candidate_catalog": artifact(args.candidate_catalog),
        "selected_input_snapshot": {"registry_path": str(snapshot_registry.resolve()), "snapshot_id": args.snapshot_id, "sha256": matches[0]["sha256"]},
        "expression_ancestry_policy": {
            "project_local_query_evidence": True,
            "cross_project_expression_requires_external_reference": True,
        },
        "whole_tissue_partition": {
            "method": args.whole_tissue_method,
            "candidate_grid_source": args.grid_source,
            "candidate_resolutions": whole_grid,
            "selection_endpoint": "broad_recall_and_purity",
        },
        "query_reclustering": {
            "candidate_resolutions": query_grid,
            "separate_from_upstream_banksy_grid": True,
        },
        "broad_family_evidence": {
            "required": True,
            "feature_scope": "full_feature",
            "complete_cartesian_product": True,
        },
        "atlas_routing": {
            "authoritative_router": "route_global_atlas_v2.py",
            "mapping_scope": "complete_analysis_set",
            "writeback_scope": "unlabeled_frozen_qc_only",
            "fine_anchor_eligible": False,
        },
        "release_taxonomy": {
            "vascular_parent": "Vascular-associated",
            "hierarchy_validation_required": True,
        },
    }
    if args.whole_tissue_grid_artifact:
        contract["whole_tissue_partition"]["grid_artifact"] = grid_artifact(args.whole_tissue_grid_artifact, whole_grid)
    out = args.out or root / "config/annotation_contract.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(contract, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(contract, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
