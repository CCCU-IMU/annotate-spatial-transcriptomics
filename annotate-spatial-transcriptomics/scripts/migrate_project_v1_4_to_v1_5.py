#!/usr/bin/env python3
"""Add v1.5 single-final-annotation, incident and resolution-contract fields."""

from __future__ import annotations

import argparse
import csv
import json
import os
import tempfile
from pathlib import Path


def add_column(path: Path, field: str) -> None:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)
    if field in fields:
        return
    fields.append(field)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader(); writer.writerows(rows); handle.flush(); os.fsync(handle.fileno())
        Path(name).replace(path)
    finally:
        Path(name).unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    for field in [
        "candidate_resolutions", "parent_pool_snapshot_id", "source_state",
        "query_membership_artifact", "residual_qc_membership_artifact",
        "residual_qc_membership_sha256", "prerequisite_route_attempt_id",
        "prerequisite_outcome_sha256", "prerequisite_residual_qc_sha256",
        "successor_state", "reroute_membership_artifact", "reroute_membership_sha256",
    ]:
        add_column(root / "state/route_attempt_registry.tsv", field)
    incident_dir = root / "provenance/incidents"; incident_dir.mkdir(parents=True, exist_ok=True)
    incident = incident_dir / "incident_registry.tsv"
    if not incident.exists():
        with incident.open("w", newline="", encoding="utf-8") as handle:
            csv.writer(handle, delimiter="\t").writerow([
                "incident_id", "detected_at_utc", "run_id", "scheduler_job_id", "scheduler_terminal",
                "failure_class", "failure_stage", "symptom", "root_cause", "failure_boundary",
                "accepted_prior_artifacts", "repair_run_id", "repair_job_id", "repair_action",
                "repair_verification", "state_mutated", "biological_labels_changed",
                "skill_prevention_candidate", "regression_test_candidate", "status", "evidence_paths",
            ])
    support = root / "state/annotation_support_registry.tsv"
    if not support.exists():
        with support.open("w", newline="", encoding="utf-8") as handle:
            csv.writer(handle, delimiter="\t").writerow([
                "support_id", "sample_id", "label_level", "broad_label", "fine_label",
                "n_observations", "confidence", "positive_marker_evidence", "anti_marker_evidence",
                "resolution_evidence", "spatial_evidence", "literature_context", "route_summary",
                "source_decision_ids", "validation_artifacts", "status", "supersedes", "created_at",
            ])
    project_path = root / "config/project.json"
    project = json.loads(project_path.read_text(encoding="utf-8"))
    project["framework_version"] = "1.5.1"
    project["required_annotation_views"] = ["final"]
    project["required_dotplot_levels"] = ["broad"]
    project["subtype_assets_required_if_final_fine_present"] = True
    project["final_broad_minimum_confidence"] = "moderate"
    project["final_fine_minimum_confidence"] = "high"
    project["preconfirmation_lightweight_review_required"] = True
    project["post_completion_master_quality_approval_required"] = True
    project_path.write_text(json.dumps(project, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "project_root": str(root), "framework_version": "1.5.1"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
