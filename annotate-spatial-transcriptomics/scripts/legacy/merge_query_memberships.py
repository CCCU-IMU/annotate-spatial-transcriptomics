#!/usr/bin/env python3
"""Legacy cross-pool merger for migration-only projects."""

from __future__ import annotations

import argparse,csv,gzip,hashlib,json
from pathlib import Path

def read(path):
    op=gzip.open if path.suffix==".gz" else open
    with op(path,"rt",newline="",encoding="utf-8") as h:return list(csv.DictReader(h,delimiter="\t"))
def main():
    p=argparse.ArgumentParser();p.add_argument("--inputs",nargs="+",type=Path,required=True);p.add_argument("--out",type=Path,required=True);p.add_argument("--pool-id",required=True);p.add_argument("--generation",type=int,required=True);a=p.parse_args();rows=[];seen=set();source_counts={}
    for path in a.inputs:
        part=read(path);source_counts[str(path.resolve())]=len(part)
        for row in part:
            cid=row.get("cell_id","")
            if not cid or cid in seen:raise SystemExit(f"missing or overlapping cell_id: {cid}")
            seen.add(cid);row["query_or_anchor"]="query";row["anchor_label"]="";row["parent_pool_id"]=a.pool_id;row["pool_snapshot_id"]=f"{a.pool_id}_g{a.generation}";row["generation"]=str(a.generation);rows.append(row)
    fields=["cell_id","query_or_anchor","anchor_label","source_key","parent_decision_id","parent_pool_id","pool_snapshot_id","generation","state_tags","spatial_tags","qc_tags","candidate_lineages"]
    a.out.parent.mkdir(parents=True,exist_ok=True);op=gzip.open if a.out.suffix==".gz" else open
    with op(a.out,"wt",newline="",encoding="utf-8") as h:w=csv.DictWriter(h,fieldnames=fields,delimiter="\t");w.writeheader();w.writerows([{k:r.get(k,"") for k in fields} for r in rows])
    hh=hashlib.sha256(a.out.read_bytes()).hexdigest();result={"status":"PASS","pool_id":a.pool_id,"generation":a.generation,"n_query":len(rows),"source_counts":source_counts,"sha256":hh};Path(str(a.out)+".manifest.json").write_text(json.dumps(result,indent=2)+"\n");print(json.dumps(result,indent=2));return 0
if __name__=="__main__":raise SystemExit(main())
