#!/usr/bin/env python3
"""Freeze all master-approved sample packets into one cohort confirmation request."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("cohort_root", type=Path)
    args = parser.parse_args()
    root = args.cohort_root.resolve()
    with (root / "control/sample_gate_registry.tsv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or any(row.get("status") != "SAMPLE_FROZEN" for row in rows):
        raise SystemExit("every sample must be SAMPLE_FROZEN before cohort confirmation")
    samples = []
    for row in rows:
        item = {"sample_id": row["sample_id"]}
        if row.get("master_decision") != "APPROVED":
            raise SystemExit(f"sample lacks main-Agent quality approval: {row['sample_id']}")
        for prefix in ["packet", "completion_gate", "master_quality_approval", "confirmation_review"]:
            path = Path(row.get(f"{prefix}_path", ""))
            if not path.is_absolute():
                path = root / path
            if not path.is_file():
                raise SystemExit(f"missing {prefix} for {row['sample_id']}")
            observed = sha(path)
            if observed != row.get(f"{prefix}_sha256"):
                raise SystemExit(f"stale {prefix} hash for {row['sample_id']}")
            item[f"{prefix}_path"] = str(path.resolve())
            item[f"{prefix}_sha256"] = observed
        samples.append(item)
    cross_path = root / "provenance/cross_sample_audit.json"
    cross = json.loads(cross_path.read_text(encoding="utf-8"))
    if cross.get("status") != "PASS":
        raise SystemExit("cross-sample audit is absent or not PASS")
    request = {"status": "COHORT_CONFIRMATION_PENDING", "n_samples": len(samples), "samples": samples, "cross_sample_audit_path": str(cross_path.resolve()), "cross_sample_audit_sha256": sha(cross_path), "created_at": datetime.now(timezone.utc).isoformat(), "warning": "Any bound packet, gate or cross-sample audit change invalidates this request."}
    out = root / "provenance/cohort_confirmation_request.json"
    out.write_text(json.dumps(request, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(request, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
