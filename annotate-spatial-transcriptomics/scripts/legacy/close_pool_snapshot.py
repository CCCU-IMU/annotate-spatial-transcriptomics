#!/usr/bin/env python3
"""Legacy pool closure for migration-only projects."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
from registry_io import locked_tsv_update


def main() -> int:
    ap=argparse.ArgumentParser();ap.add_argument("project_root",type=Path);ap.add_argument("--snapshot-id",required=True);ap.add_argument("--closure-artifact",required=True);args=ap.parse_args()
    registry=args.project_root/"state/pool_snapshot_registry.tsv";artifact=Path(args.closure_artifact);artifact=artifact if artifact.is_absolute() else args.project_root/artifact
    if not artifact.exists():raise SystemExit("closure artifact does not exist")
    changed={"value":False}
    def mutate(rows,fields):
        hits=[row for row in rows if row.get("snapshot_id")==args.snapshot_id]
        if len(hits)!=1:raise SystemExit("snapshot_id must exist exactly once")
        row=hits[0]
        if row.get("status")=="closed_and_frozen":return rows
        row["status"]="closed_and_frozen";row["closed_at"]=datetime.now(timezone.utc).isoformat();changed["value"]=True;return rows
    locked_tsv_update(registry,mutate);print(args.snapshot_id,"closed_and_frozen" if changed["value"] else "already_closed");return 0


if __name__=="__main__":raise SystemExit(main())
