#!/usr/bin/env python3
"""Legacy pool-snapshot registration for migration-only projects."""

from __future__ import annotations

import argparse,csv,gzip,hashlib
from datetime import datetime,timezone
from pathlib import Path
from registry_io import locked_tsv_update

REQUIRED={"cell_id","query_or_anchor","source_key","state_tags","spatial_tags","qc_tags","candidate_lineages"}
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
    return h.hexdigest()
def main()->int:
    p=argparse.ArgumentParser();p.add_argument("project_root",type=Path);p.add_argument("--snapshot-id",required=True);p.add_argument("--pool-id",required=True);p.add_argument("--sample",required=True);p.add_argument("--generation",required=True,type=int);p.add_argument("--membership",required=True,type=Path);p.add_argument("--parent-decision-ids",required=True);p.add_argument("--supersedes-snapshot-id",default="");p.add_argument("--status",choices=["frozen_open","closed_and_frozen"],default="frozen_open");a=p.parse_args()
    op=gzip.open if a.membership.suffix==".gz" else open
    with op(a.membership,"rt",newline="",encoding="utf-8") as h:r=csv.DictReader(h,delimiter="\t");fields=set(r.fieldnames or []);missing=REQUIRED-fields
    if missing:raise SystemExit(f"membership missing: {sorted(missing)}")
    with op(a.membership,"rt",newline="",encoding="utf-8") as h:rows=list(csv.DictReader(h,delimiter="\t"))
    ids=[x["cell_id"] for x in rows]
    if len(ids)!=len(set(ids)):raise SystemExit("duplicate membership cell IDs")
    roles=[x["query_or_anchor"] for x in rows]
    if any(x not in {"query","anchor"} for x in roles):raise SystemExit("invalid query_or_anchor value")
    if any(x["query_or_anchor"]=="anchor" and not x.get("anchor_label","") for x in rows):raise SystemExit("anchor rows require anchor_label")
    registry=a.project_root/"state/pool_snapshot_registry.tsv"
    digest=sha(a.membership);now=datetime.now(timezone.utc).isoformat();row={"snapshot_id":a.snapshot_id,"pool_id":a.pool_id,"sample_id":a.sample,"generation":a.generation,"membership_path":str(a.membership.resolve()),"membership_sha256":digest,"n_query":roles.count("query"),"n_anchors":roles.count("anchor"),"parent_decision_ids":a.parent_decision_ids,"supersedes_snapshot_id":a.supersedes_snapshot_id,"status":a.status,"created_at":now,"closed_at":now if a.status=="closed_and_frozen" else ""}
    already={"value":False}
    def mutate(old,cols):
        prior=[x for x in old if x.get("snapshot_id")==a.snapshot_id]
        if prior:
            if prior[0].get("membership_sha256")!=digest:raise SystemExit("snapshot ID already exists with a different membership")
            already["value"]=True;return old
        old.append(row);return old
    locked_tsv_update(registry,mutate)
    if already["value"]:print(a.snapshot_id,"already_registered");return 0
    print(a.snapshot_id,row["n_query"],row["n_anchors"]);return 0
if __name__=="__main__":raise SystemExit(main())
