#!/usr/bin/env python3
"""Fail-closed validation of annotation state registries and cell ledger."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path


ALLOWED = {"defined_fine", "defined_broad_only", "interface_review", "qc_holdout", "technical_state", "pending_review", "excluded_initial_qc", "closed_and_frozen"}


def rows(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def require(ok: bool, message: str, errors: list[str]) -> None:
    if not ok: errors.append(message)

def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("project_root", type=Path)
    ap.add_argument("--cell-ledger", type=Path)
    args = ap.parse_args()
    state = args.project_root / "state"
    config_path = args.project_root / "config/project.json"
    config = json.loads(config_path.read_text()) if config_path.exists() else {}
    multi_required = bool(config.get("multi_route_completion_required"))
    cell_ledger = args.cell_ledger
    if cell_ledger is None:
        for candidate in (state / "cell_ledger.tsv.gz", state / "cell_ledger.tsv"):
            if candidate.exists():
                cell_ledger = candidate
                break
    errors = []
    scope_policy_path = args.project_root / "provenance/analysis_scope_policy.json"
    expected_analysis_ids = None
    scope_policy = None
    if scope_policy_path.exists():
        try:
            scope_policy = json.loads(scope_policy_path.read_text(encoding="utf-8"))
            require(scope_policy.get("status") == "PASS", "analysis-scope policy status is not PASS", errors)
            membership_value = scope_policy.get("membership_path", "")
            membership_path = Path(membership_value)
            membership_path = membership_path if membership_path.is_absolute() else args.project_root / membership_path
            require(membership_path.exists(), "analysis-scope membership artifact is missing", errors)
            if membership_path.exists():
                expected_hash = scope_policy.get("membership_sha256", "")
                require(bool(expected_hash), "analysis-scope membership SHA256 is missing", errors)
                if expected_hash:
                    require(sha256(membership_path) == expected_hash, "analysis-scope membership SHA256 is stale", errors)
                membership_rows = list(rows(membership_path))
                require(bool(membership_rows) and "cell_id" in membership_rows[0], "analysis-scope membership lacks cell_id", errors)
                expected_analysis_ids = {row.get("cell_id", "") for row in membership_rows}
                require("" not in expected_analysis_ids, "analysis-scope membership has empty cell_id", errors)
                require(len(expected_analysis_ids) == len(membership_rows), "analysis-scope membership has duplicate cell_id", errors)
                try:
                    require(
                        len(expected_analysis_ids) == int(scope_policy.get("analysis_set_n", -1)),
                        "analysis-scope policy count disagrees with frozen membership",
                        errors,
                    )
                except (TypeError, ValueError):
                    require(False, "analysis-scope policy count is invalid", errors)
        except (OSError, json.JSONDecodeError) as exc:
            require(False, f"analysis-scope policy is unreadable: {exc}", errors)
    required = ["input_snapshot_registry.tsv", "clustering_decision_ledger.tsv", "cluster_decision_ledger.tsv", "pool_registry.tsv", "run_registry.tsv"]
    if multi_required: required += ["pool_snapshot_registry.tsv", "route_attempt_registry.tsv", "branch_control_board.tsv", "workflow_event_registry.tsv", "annotation_view_registry.tsv"]
    for f in required: require((state / f).exists(), f"missing registry: {f}", errors)
    cluster_path = state / "cluster_decision_ledger.tsv"
    active_ids = set()
    if cluster_path.exists():
        cluster_rows=list(rows(cluster_path));ids=[]
        for row in cluster_rows:ids.append(row.get("decision_id") or f"{row.get('decision_version','')}:{row.get('source_run_id','')}:{row.get('source_cluster','')}")
        require(len(ids)==len(set(ids)),"duplicate cluster decision_id",errors);known=set(ids);superseded=set()
        for row in cluster_rows:superseded.update(x for x in re.split(r"[;,\s]+",row.get("supersedes","").strip()) if x)
        active_ids=known-superseded
        for i, row in enumerate(cluster_rows, 2):
            if not row.get("source_cluster"): continue
            rid=ids[i-2];sup=[x for x in re.split(r"[;,\s]+",row.get("supersedes","").strip()) if x]
            require(rid not in sup,f"line {i}: decision supersedes itself",errors)
            for x in sup:require(x in known,f"line {i}: supersedes unknown decision {x}",errors)
            require(row.get("state") in ALLOWED, f"cluster ledger line {i}: invalid state", errors)
            if row.get("state") == "defined_fine":
                require(bool(row.get("broad_label") and row.get("fine_label")), f"line {i}: fine label lacks hierarchy", errors)
            if row.get("state") == "defined_broad_only":
                require(row.get("fine_anchor_eligible", "").lower() in {"false", "0", "no"}, f"line {i}: broad-only is fine anchor", errors)
            if rid not in superseded and row.get("closed","").lower() in {"true","1","yes"}:require(bool(row.get("closure_rationale")),f"line {i}: active closed decision lacks closure rationale",errors)
            if rid not in superseded and row.get("validation_status")=="passed":
                vp=Path(row.get("validation_artifact",""));vp=vp if vp.is_absolute() else args.project_root/vp
                require(vp.exists(),f"line {i}: passed validation artifact missing",errors)
    if cell_ledger:
        seen = set(); n = 0; it=iter(rows(cell_ledger));first=next(it,None);required_fields={"sample_id","cell_id","decision_id","state","broad_label","fine_label","validation_feature_scope"}
        if multi_required: required_fields.update({"analysis_scope","source_key","parent_decision_id","pool_snapshot_id","generation","route_attempt_id","state_tags","spatial_tags","qc_tags","candidate_lineages","strict_state","strict_broad_label","strict_fine_label","inclusive_state","inclusive_broad_label","inclusive_fine_label","display_state","display_broad_label","display_fine_label","display_policy"})
        fields=set(first or {});missing_fields=required_fields-fields;require(not missing_fields,f"cell ledger missing fields: {sorted(missing_fields)}",errors);oocyte_object_missing=False;analysis_n=0;analysis_ids_seen=set();excluded_n=0
        import itertools
        for i, row in enumerate(itertools.chain([first] if first else [],it), 2):
            key = (row.get("sample_id"), row.get("cell_id")); n += 1
            require(key not in seen, f"cell ledger duplicate at line {i}: {key}", errors); seen.add(key)
            require(row.get("state") in ALLOWED, f"cell ledger line {i}: invalid state", errors)
            if row.get("state") == "defined_broad_only":
                require(row.get("fine_anchor_eligible", "").lower() in {"false", "0", "no"}, f"cell line {i}: broad-only is fine anchor", errors)
            if "decision_id" in fields:require(row.get("decision_id") in active_ids,f"cell line {i}: decision_id is not active",errors)
            if "oocyte" in (row.get("broad_label","")+row.get("fine_label","")).lower() and not row.get("spatial_object_id"):oocyte_object_missing=True
            if multi_required:
                scope=row.get("analysis_scope","");require(scope in {"analysis_set","excluded_initial_qc"},f"cell line {i}: invalid analysis_scope",errors)
                if scope=="analysis_set":
                    analysis_n+=1
                    analysis_ids_seen.add(row.get("cell_id", ""))
                    for view in ["strict","inclusive","display"]:require(bool(row.get(f"{view}_state")),f"cell line {i}: {view} view is empty",errors)
                elif scope=="excluded_initial_qc":
                    excluded_n+=1
                    for view in ["strict", "inclusive", "display"]:
                        require(row.get(f"{view}_state") == "excluded_initial_qc", f"cell line {i}: excluded observation leaked into {view} state", errors)
                        require(not row.get(f"{view}_broad_label"), f"cell line {i}: excluded observation has {view} broad label", errors)
                        require(not row.get(f"{view}_fine_label"), f"cell line {i}: excluded observation has {view} fine label", errors)
                if row.get("route") in {"depth_matched_atlas_anchor_mapping_fullfeature_rescue","qc_anchor_recluster_broad_rescue","interface_rctd_broad_return"}:
                    require(row.get("fine_anchor_eligible","").lower() in {"false","0","no"},f"cell line {i}: rescue route is a fine anchor",errors)
                    require(not row.get("strict_broad_label"),f"cell line {i}: calibrated broad rescue leaked into strict biological view",errors)
        require(not oocyte_object_missing,"Oocyte-associated cellbins lack spatial_object_id",errors)
        if multi_required:
            if expected_analysis_ids is not None:
                require(analysis_ids_seen == expected_analysis_ids, "cell-ledger analysis_set differs from frozen analysis-scope membership", errors)
                try:
                    expected_excluded = int(scope_policy.get("excluded_initial_qc_n", -1)) if scope_policy else -1
                    require(excluded_n == expected_excluded, "cell-ledger excluded_initial_qc count disagrees with policy", errors)
                except (TypeError, ValueError):
                    require(False, "analysis-scope excluded count is invalid", errors)
            view_rows=list(rows(state/"annotation_view_registry.tsv"));validated={x.get("view"):x for x in view_rows if x.get("status")=="validated"}
            for view in ["strict","inclusive","display"]:
                require(view in validated,f"missing validated {view} annotation view",errors)
                if view in validated:
                    try:require(int(float(validated[view].get("n_observations",0)))==analysis_n,f"{view} view count differs from analysis set",errors)
                    except ValueError:require(False,f"{view} view count is invalid",errors)
    validated={str(cluster_path.relative_to(args.project_root)):sha256(cluster_path)} if cluster_path.exists() else {}
    if cell_ledger and cell_ledger.exists():validated[str(cell_ledger.resolve().relative_to(args.project_root.resolve()))]=sha256(cell_ledger)
    result = {"status": "PASS" if not errors else "FAIL", "errors": errors,"validated_files":validated}
    out = args.project_root / "provenance/state_validation.json"; out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
