#!/usr/bin/env python3
"""Build the next-action queue for the direct cohort annotation architecture."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

from validate_direct_lineage_workflow import audit as audit_direct
from validate_profile_role import load_profile


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def active_decisions(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    def decision_id(row: dict[str, str]) -> str:
        return row.get("decision_id") or f"{row.get('decision_version','')}:{row.get('source_run_id','')}:{row.get('source_cluster','')}"

    superseded: set[str] = set()
    for row in rows:
        superseded.update(item for item in re.split(r"[;,\s]+", row.get("supersedes", "").strip()) if item)
    return [row for row in rows if decision_id(row) not in superseded]


def classify_error(message: str) -> tuple[int, str, str, str]:
    if "initial broad class lacks" in message:
        label = message.rsplit(":", 1)[-1].strip()
        return 1, "broad_class_recluster_or_underpowered_skip", label, "terminal broad-class cohort or recorded underpowered skip"
    if "recluster cohort" in message or message.startswith("cohort "):
        cohort = message.split(":", 1)[0].replace("cohort ", "").strip()
        return 1, "repair_or_finish_recluster_cohort", cohort, "validated terminal cohort artifact and writeback"
    if "direct return" in message:
        return 1, "repair_direct_return", message.split(":", 1)[0], "validated membership, evidence and confidence"
    if "Atlas" in message or "residual QC Atlas" in message:
        return 2, "residual_qc_calibrated_atlas_review", "final_qc_holdout", "terminal calibrated broad-only rescue/QC partition"
    if "final annotation" in message or "final" in message and "view" in message:
        return 4, "build_final_annotation", "analysis_set", "single final broad/fine/QC annotation"
    if "cell ledger" in message:
        return 0, "write_or_repair_cell_ledger", "analysis_set", "valid cell-level state"
    return 2, "repair_annotation_workflow_gap", "project", "validated state writeback"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--context", required=True, type=Path)
    parser.add_argument("--biological-profile", type=Path)
    parser.add_argument("--profile", type=Path, help="deprecated biological-profile alias")
    args = parser.parse_args()
    root = args.project_root.resolve()
    context = json.loads(args.context.read_text(encoding="utf-8"))
    profile_path = args.biological_profile or args.profile
    profile = load_profile(profile_path, "biological_evidence") if profile_path else {}
    project_path = root / "config/project.json"
    project = json.loads(project_path.read_text(encoding="utf-8")) if project_path.is_file() else {}

    if project.get("routing_model", "direct_cross_lineage_recluster_cohorts") not in {"direct_cross_lineage_recluster_cohorts", "direct_cross_branch_recluster_cohorts"}:
        raise SystemExit("legacy pool project detected; migrate it before using the current planner")

    cluster_rows = active_decisions(read_tsv(root / "state/cluster_decision_ledger.tsv"))
    result = audit_direct(root)
    rows: list[dict[str, object]] = []

    if profile.get("final_validation"):
        full_feature = root / "provenance/full_feature_validation.json"
        passed = False
        if full_feature.is_file():
            try:
                passed = json.loads(full_feature.read_text(encoding="utf-8")).get("status") == "PASS"
            except json.JSONDecodeError:
                passed = False
        if not passed:
            rows.append({
                "priority": 0, "source_run_id": "project", "source_cluster": "__FULL_FEATURE_VALIDATION__",
                "n_observations": 0, "current_state": "project_gate", "broad_label": "", "fine_label": "",
                "required_route": "full_feature_marker_validation",
                "reason": "final marker/anti-marker and context-specific evidence require a validated full-feature layer",
                "target_scope": "whole_tissue_and_all_recluster_cohorts",
                "blocked_until": "full-feature audit PASS and evidence writeback",
            })

    seen: set[tuple[str, str]] = set()
    for message in result.get("errors", []):
        priority, route, target, blocked = classify_error(message)
        key = (route, target)
        if key in seen:
            continue
        seen.add(key)
        rows.append({
            "priority": priority, "source_run_id": "project", "source_cluster": target,
            "n_observations": 0, "current_state": "annotation_workflow_gap", "broad_label": target if route.startswith("broad_class") else "",
            "fine_label": "", "required_route": route, "reason": message,
            "target_scope": target, "blocked_until": blocked,
        })

    for row in cluster_rows:
        if row.get("state") in {"defined_fine", "defined_broad_only"} and row.get("validation_feature_scope") != "full_feature" and profile.get("final_validation"):
            key = ("full_feature_marker_validation", row.get("source_cluster", ""))
            if key not in seen:
                seen.add(key)
                rows.append({
                    "priority": 1, "source_run_id": row.get("source_run_id", ""), "source_cluster": row.get("source_cluster", ""),
                    "n_observations": int(float(row.get("n_observations", 0) or 0)), "current_state": row.get("state", ""),
                    "broad_label": row.get("broad_label", ""), "fine_label": row.get("fine_label", ""),
                    "required_route": "full_feature_marker_validation", "reason": "biological decision lacks full-feature validation",
                    "target_scope": row.get("recluster_cohort_id", ""), "blocked_until": "full-feature positive/anti-marker evidence",
                })

    rows.sort(key=lambda row: (int(row["priority"]), -int(row["n_observations"])))
    queue_path = root / "state/next_action_queue.tsv"
    queue_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "priority", "source_run_id", "source_cluster", "n_observations", "current_state", "broad_label",
        "fine_label", "required_route", "reason", "target_scope", "blocked_until",
    ]
    with queue_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    plan = {
        "status": "ITERATION_REQUIRED" if rows else "READY_FOR_COMPLETION_AUDIT",
        "workflow_model": result.get("workflow_model"),
        "queued_actions": len(rows),
        "queued_observations": sum(int(row["n_observations"]) for row in rows),
        "active_decisions": len(cluster_rows),
        "workflow_audit_status": result.get("status"),
        "context": context,
    }
    provenance = root / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    (provenance / "direct_lineage_workflow_audit.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (provenance / "iteration_plan.json").write_text(json.dumps(plan, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(plan, ensure_ascii=False, indent=2))
    return 2 if rows else 0


if __name__ == "__main__":
    raise SystemExit(main())
