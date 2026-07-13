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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    project = json.loads((root / "config/project.json").read_text(encoding="utf-8"))
    cell_ledger = root / "state/cell_ledger.tsv.gz"
    cluster_ledger = root / "state/cluster_decision_ledger.tsv"
    completion_gate = root / "provenance/completion_gate.json"
    if not all(path.is_file() for path in (cell_ledger, cluster_ledger, completion_gate)):
        raise SystemExit("cell ledger, cluster ledger and completion gate are required")
    completion = json.loads(completion_gate.read_text(encoding="utf-8"))
    if completion.get("status") != "PASS":
        raise SystemExit("completion gate must pass before requesting final user confirmation")

    state_counts: Counter[str] = Counter()
    strict_counts: Counter[str] = Counter()
    strict_broad_counts: Counter[str] = Counter()
    inclusive_broad_counts: Counter[str] = Counter()
    display_broad_counts: Counter[str] = Counter()
    strict_unresolved_state_counts: Counter[str] = Counter()
    analysis_n = 0
    excluded_n = 0
    strict_unresolved_n = 0
    with gzip.open(cell_ledger, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "analysis_scope", "state", "strict_state", "strict_broad_label",
            "inclusive_broad_label", "display_broad_label",
        }
        if not required.issubset(reader.fieldnames or []):
            raise SystemExit("cell ledger lacks analysis-scope/strict annotation fields")
        for row in reader:
            if row["analysis_scope"] == "analysis_set":
                analysis_n += 1
                state_counts[row["state"] or "blank"] += 1
                strict_counts[row["strict_state"] or "blank"] += 1
                strict_broad_counts[row["strict_broad_label"] or "unresolved"] += 1
                inclusive_broad_counts[row["inclusive_broad_label"] or "unresolved"] += 1
                display_broad_counts[row["display_broad_label"] or "unresolved"] += 1
                if not row["strict_broad_label"]:
                    strict_unresolved_n += 1
                    strict_unresolved_state_counts[row["strict_state"] or "blank"] += 1
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
        "full_object_n": analysis_n + excluded_n,
        "analysis_set_n": analysis_n,
        "excluded_initial_qc_n": excluded_n,
        "strict_unresolved_n": strict_unresolved_n,
        "strict_unresolved_fraction": strict_unresolved_n / analysis_n if analysis_n else None,
        "analysis_set_state_counts": dict(state_counts),
        "strict_state_counts": dict(strict_counts),
        "strict_broad_counts": dict(strict_broad_counts),
        "inclusive_broad_counts": dict(inclusive_broad_counts),
        "display_broad_counts": dict(display_broad_counts),
        "strict_unresolved_state_counts": dict(strict_unresolved_state_counts),
        "user_must_review": [
            "broad-class census and unresolved/QC/interface counts",
            "rare-cell and context-sensitive calls",
            "atlas/RCTD broad-only returns and retained rejects",
            "strict versus inclusive/display interpretation",
        ],
        "release_rule": "Do not compute final DEG/dotplot/spatial/report assets and do not write a confirmation artifact until the user explicitly confirms this frozen annotation snapshot.",
    }
    out = root / "provenance/final_annotation_confirmation_request.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(request, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(request, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
