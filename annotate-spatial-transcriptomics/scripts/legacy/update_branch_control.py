#!/usr/bin/env python3
"""Legacy branch-board updater for migration-only projects."""

from __future__ import annotations

import argparse,csv,json
from datetime import datetime,timezone
from pathlib import Path
from registry_io import locked_tsv_update

def main()->int:
    p=argparse.ArgumentParser();p.add_argument("project_root",type=Path);p.add_argument("--record",required=True,type=Path);a=p.parse_args();path=a.project_root/"state/branch_control_board.tsv";record=json.loads(a.record.read_text());required=["branch_id","sample_id","parent_decision_id","pool_snapshot_id","generation","current_state","recluster_policy","terminal","next_action","authoritative_artifact"];missing=[x for x in required if str(record.get(x,"")).strip()==""]
    if missing:raise SystemExit(f"missing branch fields: {missing}")
    record.setdefault("updated_at",datetime.now(timezone.utc).isoformat())
    def mutate(rows,fields):
        old=[x for x in rows if x.get("branch_id")==record["branch_id"]]
        if old and str(old[0].get("terminal","")).lower() in {"true","1","yes"}:raise SystemExit("terminal branch is immutable; create a superseding generation")
        rows=[x for x in rows if x.get("branch_id")!=record["branch_id"]];rows.append({x:record.get(x,"") for x in fields});return rows
    locked_tsv_update(path,mutate)
    print(record["branch_id"],record["current_state"]);return 0
if __name__=="__main__":raise SystemExit(main())
