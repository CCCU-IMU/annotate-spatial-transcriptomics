#!/usr/bin/env python3
"""Add v1.7 evidence-freeze and global-Atlas fields without inventing evidence."""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


CLUSTER_FIELDS = [
    "prelabel_evidence_artifact", "prelabel_evidence_sha256", "prelabel_evidence_frozen",
    "prelabel_winner", "prelabel_runner_up", "prelabel_winning_margin",
]
ROUTE_FIELDS = [
    "qc_membership_artifact", "qc_membership_sha256", "qc_membership_n_observations",
    "concordance_artifact", "concordance_artifact_sha256",
    "cluster_concordance_artifact", "cluster_concordance_artifact_sha256",
    "discrepancy_review_artifact", "discrepancy_review_artifact_sha256",
    "n_defined_concordant", "n_defined_weak_challenge", "n_defined_discordant",
    "n_defined_ood", "n_ontology_conflict",
]


def add_fields(path: Path, additions: list[str], backup: Path) -> None:
    if not path.is_file():
        return
    shutil.copy2(path, backup / path.name)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    fields.extend(field for field in additions if field not in fields)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    project_path = root / "config/project.json"
    if not project_path.is_file():
        raise SystemExit("missing config/project.json")
    project = json.loads(project_path.read_text(encoding="utf-8"))
    if project.get("framework_version") not in {"1.6.0", "1.6.1", "1.7.0"}:
        raise SystemExit(f"unsupported source framework version: {project.get('framework_version')}")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup = root / "provenance/migration_backups" / f"v1.6.1_to_v1.7.0_{stamp}"
    backup.mkdir(parents=True, exist_ok=False)
    shutil.copy2(project_path, backup / "project.json")
    add_fields(root / "state/cluster_decision_ledger.tsv", CLUSTER_FIELDS, backup)
    add_fields(root / "state/route_attempt_registry.tsv", ROUTE_FIELDS, backup)
    project.update({
        "framework_version": "1.7.0",
        "routing_model": "direct_cross_lineage_recluster_cohorts_global_atlas",
        "prelabel_evidence_freeze_required": True,
        "global_atlas_concordance_required_when_reference_applicable": True,
        "global_atlas_mapping_scope": "complete_analysis_set",
        "global_atlas_mapping_ceiling": "broad_only",
        "global_atlas_ood_rejection_required": True,
        "migrated_at_utc": datetime.now(timezone.utc).isoformat(),
    })
    project_path.write_text(json.dumps(project, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: migrated registry schemas to v1.7.0; evidence remains intentionally unfilled; backup: {backup}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
