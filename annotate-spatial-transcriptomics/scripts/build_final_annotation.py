#!/usr/bin/env python3
"""Build one final annotation: moderate+ broad labels and high-confidence fine labels."""

from __future__ import annotations

import argparse
import csv
import fcntl
import gzip
import hashlib
import os
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


RANK = {
    "low": 0, "low_confidence": 0, "medium_low": 0,
    "moderate": 1, "moderate_only": 1, "medium": 1, "medium_high": 1, "moderate_high": 1,
    "high": 2, "high_confidence": 2, "very_high": 2,
}
BIOLOGICAL = {"defined_fine", "defined_broad_only"}


def truth(value: str) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y", "pass", "passed"}


def confidence(value: str) -> tuple[str, int]:
    normalized = str(value).strip().lower().replace("-", "_").replace(" ", "_")
    return normalized, RANK.get(normalized, -1)


def open_table(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, newline="", encoding="utf-8")
    return path.open(mode.replace("t", ""), newline="", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--cell-ledger", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--version", default="v001")
    args = parser.parse_args()
    final_fields = [
        "analysis_scope", "final_state", "final_broad_label", "final_fine_label", "final_confidence",
        "final_broad_confidence", "final_fine_confidence",
        "final_assignment_tier", "final_broad_eligible", "final_fine_eligible",
    ]
    counts: Counter[tuple[str, str, str, str]] = Counter()
    analysis_n = broad_n = fine_n = 0
    args.out.parent.mkdir(parents=True, exist_ok=True)
    temp = args.out.with_name(args.out.stem + ".tmp" + args.out.suffix)
    with open_table(args.cell_ledger, "rt") as source:
        reader = csv.DictReader(source, delimiter="\t")
        fields = reader.fieldnames or []
        required = {"cell_id", "state", "broad_label", "fine_label", "confidence", "fine_anchor_eligible"}
        missing = sorted(required - set(fields))
        if missing:
            raise SystemExit("cell ledger lacks final-annotation fields: " + ",".join(missing))
        output_fields = fields + [field for field in final_fields if field not in fields]
        with open_table(temp, "wt") as target:
            writer = csv.DictWriter(target, fieldnames=output_fields, delimiter="\t")
            writer.writeheader()
            for row in reader:
                scope = row.get("analysis_scope") or ("excluded_initial_qc" if row.get("state") == "excluded_initial_qc" else "analysis_set")
                row["analysis_scope"] = scope
                state = row.get("state", "")
                broad_confidence, broad_rank = confidence(row.get("broad_confidence") or row.get("confidence", ""))
                fine_confidence, fine_rank = confidence(row.get("fine_confidence") or row.get("confidence", ""))
                route_text = " ".join([row.get("route", ""), row.get("evidence_status", "")]).lower()
                legacy_mapping_tier = any(token in route_text for token in ["atlas", "calibrated", "mapping"]) and broad_confidence in {"medium_high", "moderate_high"}
                if legacy_mapping_tier:
                    broad_rank = -1
                if scope != "analysis_set":
                    final_state, broad, fine, tier = "excluded_initial_qc", "", "", "excluded"
                elif state in BIOLOGICAL and row.get("broad_label") and broad_rank >= 1:
                    broad = row.get("broad_label", "")
                    fine_ok = state == "defined_fine" and bool(row.get("fine_label")) and fine_rank >= 2 and truth(row.get("fine_anchor_eligible", ""))
                    fine = row.get("fine_label", "") if fine_ok else ""
                    final_state = "defined_fine" if fine else "defined_broad_only"
                    tier = "high_fine" if fine else ("high_broad" if broad_rank >= 2 else "moderate_broad")
                else:
                    retained_state = "pending_review" if state in BIOLOGICAL else (state or "pending_review")
                    final_state, broad, fine, tier = retained_state, "", "", "retained_or_unresolved"
                row.update({
                    "final_state": final_state, "final_broad_label": broad, "final_fine_label": fine,
                    "final_confidence": broad_confidence,
                    "final_broad_confidence": broad_confidence if broad else "",
                    "final_fine_confidence": fine_confidence if fine else "",
                    "final_assignment_tier": tier,
                    "final_broad_eligible": "true" if broad else "false",
                    "final_fine_eligible": "true" if fine else "false",
                })
                if scope == "analysis_set":
                    analysis_n += 1; broad_n += bool(broad); fine_n += bool(fine)
                    counts[(final_state, broad, fine, tier)] += 1
                writer.writerow({field: row.get(field, "") for field in output_fields})
    temp.replace(args.out)
    tables = args.project_root / "tables"; tables.mkdir(parents=True, exist_ok=True)
    census = tables / "final_annotation_census.tsv"
    with census.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["state", "broad_label", "fine_label", "assignment_tier", "n_observations"])
        for key, n in counts.most_common():
            writer.writerow([*key, n])
    registry = args.project_root / "state/annotation_view_registry.tsv"
    with registry.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t"); registry_fields = reader.fieldnames or []; rows = [row for row in reader if row.get("view") != "final"]
    rows.append({
        "view_id": f"{args.sample}_final_{args.version}", "sample_id": args.sample, "view": "final",
        "membership_path": str(args.out.resolve()), "membership_sha256": sha256(args.out),
        "n_observations": analysis_n,
        "policy": "single final annotation: moderate-or-higher broad; high-confidence evidence-gated fine",
        "marker_deg_eligible": "true", "status": "validated", "artifact": str(census.resolve()),
        "created_at": datetime.now(timezone.utc).isoformat(),
    })
    lock = registry.with_suffix(registry.suffix + ".lock")
    with lock.open("w") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        fd, name = tempfile.mkstemp(prefix=registry.name + ".", suffix=".tmp", dir=registry.parent)
        try:
            with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=registry_fields, delimiter="\t")
                writer.writeheader(); writer.writerows(rows); handle.flush(); os.fsync(handle.fileno())
            Path(name).replace(registry)
        finally:
            Path(name).unlink(missing_ok=True)
    print(f"PASS analysis={analysis_n} final_broad={broad_n} final_fine={fine_n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
