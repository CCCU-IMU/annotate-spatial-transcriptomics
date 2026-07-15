#!/usr/bin/env python3
"""Freeze the annotation census that must be shown to the user before release assets."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from master_quality_lib import validate_master_approval


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--review-report", default="review/confirmation/index.html")
    parser.add_argument("--review-manifest", default="provenance/confirmation_review_manifest.json")
    args = parser.parse_args()
    root = args.project_root.resolve()
    project = json.loads((root / "config/project.json").read_text(encoding="utf-8"))
    cell_ledger = root / "state/cell_ledger.tsv.gz"
    cluster_ledger = root / "state/cluster_decision_ledger.tsv"
    completion_gate = root / "provenance/completion_gate.json"
    taxonomy_audit = root / "provenance/release_taxonomy_audit.json"
    review_report = root / args.review_report
    review_manifest = root / args.review_manifest
    master_approval_path = root / "state/master_quality_approval.json"
    if not all(path.is_file() for path in (cell_ledger, cluster_ledger, completion_gate, taxonomy_audit, master_approval_path, review_report, review_manifest)):
        raise SystemExit("cell/cluster ledgers, audits and the lightweight confirmation review are required")
    completion = json.loads(completion_gate.read_text(encoding="utf-8"))
    if completion.get("status") != "PASS":
        raise SystemExit("completion gate must pass before requesting final user confirmation")
    taxonomy = json.loads(taxonomy_audit.read_text(encoding="utf-8"))
    if taxonomy.get("pass") is not True or taxonomy.get("metadata_sha256") != sha256(cell_ledger):
        raise SystemExit("release taxonomy audit must pass on the current cell ledger")
    master_ok, master_errors, _ = validate_master_approval(root, master_approval_path)
    if not master_ok:
        raise SystemExit("main-Agent annotation-quality approval is missing or stale: " + "; ".join(master_errors))
    review = json.loads(review_manifest.read_text(encoding="utf-8"))
    expected_review_hashes = {
        "cell_ledger_sha256": sha256(cell_ledger),
        "cluster_ledger_sha256": sha256(cluster_ledger),
        "completion_gate_sha256": sha256(completion_gate),
        "report_sha256": sha256(review_report),
        "master_quality_approval_sha256": sha256(master_approval_path),
    }
    if review.get("status") != "PASS" or any(review.get(key) != value for key, value in expected_review_hashes.items()):
        raise SystemExit("lightweight confirmation review is failed or stale for the frozen annotation")
    support_registry = Path(str(review.get("support_registry", "")))
    if not support_registry.is_file() or review.get("support_registry_sha256") != sha256(support_registry):
        raise SystemExit("annotation support reasons changed after the lightweight review")

    state_counts: Counter[str] = Counter()
    final_state_counts: Counter[str] = Counter()
    final_broad_counts: Counter[str] = Counter()
    final_fine_counts: Counter[str] = Counter()
    retained_state_counts: Counter[str] = Counter()
    analysis_n = 0
    excluded_n = 0
    final_unresolved_n = 0
    with gzip.open(cell_ledger, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"analysis_scope", "state", "final_state", "final_broad_label", "final_fine_label"}
        if not required.issubset(reader.fieldnames or []):
            raise SystemExit("cell ledger lacks single-final-annotation fields")
        for row in reader:
            if row["analysis_scope"] == "analysis_set":
                analysis_n += 1
                state_counts[row["state"] or "blank"] += 1
                final_state_counts[row["final_state"] or "blank"] += 1
                final_broad_counts[row["final_broad_label"] or "unresolved"] += 1
                if row["final_fine_label"]:
                    final_fine_counts[row["final_fine_label"]] += 1
                if not row["final_broad_label"]:
                    final_unresolved_n += 1
                    retained_state_counts[row["final_state"] or "blank"] += 1
            else:
                excluded_n += 1

    request = {
        "status": "AWAITING_USER_CONFIRMATION",
        "sample_id": project.get("sample_id"),
        "decision_version": project.get("decision_version"),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "cell_ledger": "state/cell_ledger.tsv.gz",
        "cell_ledger_sha256": sha256(cell_ledger),
        "cluster_ledger": "state/cluster_decision_ledger.tsv",
        "cluster_ledger_sha256": sha256(cluster_ledger),
        "completion_gate": "provenance/completion_gate.json",
        "completion_gate_sha256": sha256(completion_gate),
        "release_taxonomy_audit": "provenance/release_taxonomy_audit.json",
        "release_taxonomy_audit_sha256": sha256(taxonomy_audit),
        "master_quality_approval": "state/master_quality_approval.json",
        "master_quality_approval_sha256": sha256(master_approval_path),
        "confirmation_review_report": args.review_report,
        "confirmation_review_report_sha256": sha256(review_report),
        "confirmation_review_manifest": args.review_manifest,
        "confirmation_review_manifest_sha256": sha256(review_manifest),
        "annotation_support_registry": str(support_registry),
        "annotation_support_registry_sha256": sha256(support_registry),
        "biological_broad_census": taxonomy.get("biological_broad_census", {}),
        "retained_state_census": taxonomy.get("retained_state_census", {}),
        "full_object_n": analysis_n + excluded_n,
        "analysis_set_n": analysis_n,
        "excluded_initial_qc_n": excluded_n,
        "final_unresolved_n": final_unresolved_n,
        "final_unresolved_fraction": final_unresolved_n / analysis_n if analysis_n else None,
        "analysis_set_state_counts": dict(state_counts),
        "final_state_counts": dict(final_state_counts),
        "final_broad_counts": dict(final_broad_counts),
        "final_fine_counts": dict(final_fine_counts),
        "retained_state_counts": dict(retained_state_counts),
        "user_must_review": [
            "broad-class census and unresolved/QC/interface counts",
            "rare-cell and context-sensitive calls",
            "atlas/RCTD broad-only returns and retained rejects",
            "one final annotation: moderate-or-higher broad labels and high-confidence fine labels",
            "separate biological broad-class and retained anatomical/QC/technical-state censuses",
            "lightweight review HTML: support/anti-marker reasons, distinguishable broad spatial colors and canonical broad marker dotplot",
            "main-Agent quality approval performed only after all annotation routes, final writeback and completion gate passed",
        ],
        "release_rule": "Only the lightweight confirmation HTML/assets are allowed before approval. Do not compute final DEG, full tree dotplots, per-node/per-gene spatial assets or the release HTML until the user confirms this frozen snapshot.",
    }
    out = root / "provenance/final_annotation_confirmation_request.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(request, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(request, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
