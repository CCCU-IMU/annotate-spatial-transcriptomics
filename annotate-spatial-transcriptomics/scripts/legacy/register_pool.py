#!/usr/bin/env python3
"""Legacy pool registration for migration-only projects."""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json
from datetime import datetime,timezone
from pathlib import Path

FIELDS=["pool_id","sample_id","parent_pool_id","membership_path","membership_sha256","n_observations","purpose","status","decision_version","created_at","closed_at"]
def open_text(p):return gzip.open(p,"rt",encoding="utf-8",newline="") if p.suffix==".gz" else p.open(newline="",encoding="utf-8")
def main():
    p=argparse.ArgumentParser();p.add_argument("project_root",type=Path);p.add_argument("--pool-id",required=True);p.add_argument("--membership",required=True,type=Path);p.add_argument("--cell-id-col",default="cell_id");p.add_argument("--parent-pool-id",required=True);p.add_argument("--purpose",required=True);p.add_argument("--status",default="open_pending_review",choices=["open_pending_review","running","closed_and_frozen","explicitly_retained_closed"]);p.add_argument("--decision-version",required=True);a=p.parse_args()
    with open_text(a.membership) as h:ids=[r[a.cell_id_col] for r in csv.DictReader(h,delimiter="\t")]
    if not ids or len(ids)!=len(set(ids)):raise SystemExit("membership must be nonempty with unique IDs")
    sha=hashlib.sha256(a.membership.read_bytes()).hexdigest();cfg=json.loads((a.project_root/"config/project.json").read_text());path=a.project_root/"state/pool_registry.tsv";rows=[]
    if path.exists():
        with path.open(newline="",encoding="utf-8") as h:rows=list(csv.DictReader(h,delimiter="\t"))
    old=[r for r in rows if r.get("pool_id")==a.pool_id]
    if old:
        if old[0].get("membership_sha256")!=sha:raise SystemExit("pool_id already exists with different membership; create a new versioned pool_id")
        if old[0].get("status") in {"closed_and_frozen","explicitly_retained_closed"} and a.status!=old[0].get("status"):raise SystemExit("closed pool is immutable")
        rows=[r for r in rows if r.get("pool_id")!=a.pool_id]
    now=datetime.now(timezone.utc).isoformat();row={"pool_id":a.pool_id,"sample_id":cfg["sample_id"],"parent_pool_id":a.parent_pool_id,"membership_path":str(a.membership.resolve()),"membership_sha256":sha,"n_observations":len(ids),"purpose":a.purpose,"status":a.status,"decision_version":a.decision_version,"created_at":old[0].get("created_at",now) if old else now,"closed_at":now if a.status in {"closed_and_frozen","explicitly_retained_closed"} else ""};rows.append(row)
    with path.open("w",newline="",encoding="utf-8") as h:w=csv.DictWriter(h,fieldnames=FIELDS,delimiter="\t");w.writeheader();w.writerows(rows)
    print(json.dumps(row,ensure_ascii=False))
if __name__=="__main__":main()
