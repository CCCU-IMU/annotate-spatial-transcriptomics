#!/usr/bin/env python3
"""Append one immutable chronological workflow event from JSON."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from registry_io import locked_tsv_update


def main() -> int:
    parser=argparse.ArgumentParser();parser.add_argument("project_root",type=Path);parser.add_argument("--record",required=True,type=Path);args=parser.parse_args()
    path=args.project_root/"state/workflow_event_registry.tsv";record=json.loads(args.record.read_text());required=["event_id","sample_id","phase","action","status","decision_summary_zh"]
    missing=[x for x in required if not record.get(x)];
    if missing:raise SystemExit(f"missing workflow event fields: {missing}")
    record.setdefault("timestamp",datetime.now(timezone.utc).isoformat())
    def mutate(rows,fields):
        if any(x.get("event_id")==record["event_id"] for x in rows):raise SystemExit("event_id already exists; workflow events are immutable")
        rows.append({x:record.get(x,"") for x in fields});return rows
    locked_tsv_update(path,mutate)
    print(record["event_id"]);return 0


if __name__=="__main__":raise SystemExit(main())
