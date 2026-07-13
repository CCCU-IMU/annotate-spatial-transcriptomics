#!/usr/bin/env python3
"""Audit route applicability, execution and final annotation views."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from multiroute_lib import audit_multiroute


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--context", required=True, type=Path)
    parser.add_argument("--profile", required=True, type=Path)
    args = parser.parse_args()
    result = audit_multiroute(
        args.project_root,
        json.loads(args.context.read_text()),
        json.loads(args.profile.read_text()),
    )
    out = args.project_root / "provenance/multiroute_audit.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    queue = args.project_root / "state/multiroute_gap_queue.tsv"
    fields = ["decision_id", "n_observations", "required_route", "reason"]
    with queue.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(result["gaps"])
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
