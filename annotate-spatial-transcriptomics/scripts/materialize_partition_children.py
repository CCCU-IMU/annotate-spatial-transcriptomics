#!/usr/bin/env python3
"""Materialize open child-pool query memberships from a conserved partition."""

from __future__ import annotations

import argparse,csv,gzip,hashlib,json,re
from pathlib import Path

def read(path):
    op=gzip.open if path.suffix==".gz" else open
    with op(path,"rt",newline="",encoding="utf-8") as h:return list(csv.DictReader(h,delimiter="\t"))
def write(path,rows,fields):
    op=gzip.open if path.suffix==".gz" else open
    with op(path,"wt",newline="",encoding="utf-8") as h:w=csv.DictWriter(h,fieldnames=fields,delimiter="\t");w.writeheader();w.writerows(rows)
def sha(path):
    h=hashlib.sha256();
    with path.open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
    return h.hexdigest()
def safe(x):return re.sub(r"[^A-Za-z0-9_.-]+","_",x).strip("_")
def main():
    p=argparse.ArgumentParser();p.add_argument("partition",type=Path);p.add_argument("--target-pools",required=True);p.add_argument("--out-dir",type=Path,required=True);p.add_argument("--generation",type=int,required=True);a=p.parse_args();rows=read(a.partition);targets=[x.strip() for x in a.target_pools.split(",") if x.strip()];a.out_dir.mkdir(parents=True,exist_ok=True);index=[]
    fields=["cell_id","query_or_anchor","anchor_label","source_key","parent_decision_id","parent_pool_id","pool_snapshot_id","generation","state_tags","spatial_tags","qc_tags","candidate_lineages"]
    for target in targets:
        subset=[x for x in rows if x.get("target_pool")==target]
        if not subset:raise SystemExit(f"target pool has no cells: {target}")
        out=[]
        for x in subset:out.append({"cell_id":x["cell_id"],"query_or_anchor":"query","anchor_label":"","source_key":x.get("source_key",""),"parent_decision_id":x.get("parent_decision_id",""),"parent_pool_id":target,"pool_snapshot_id":f"{safe(target)}_g{a.generation}","generation":a.generation,"state_tags":x.get("state_tags",""),"spatial_tags":x.get("spatial_tags",""),"qc_tags":x.get("qc_tags",""),"candidate_lineages":x.get("provisional_broad_label",x.get("candidate_lineages",""))})
        path=a.out_dir/f"{safe(target)}_g{a.generation}.tsv.gz";write(path,out,fields);index.append({"target_pool":target,"generation":a.generation,"n_query":len(out),"membership":str(path.resolve()),"sha256":sha(path)})
    idx=a.out_dir/f"child_pool_index_g{a.generation}.tsv";write(idx,index,list(index[0]));print(json.dumps({"status":"PASS","n_children":len(index),"n_query":sum(x["n_query"] for x in index),"index":str(idx.resolve())},indent=2));return 0
if __name__=="__main__":raise SystemExit(main())
