#!/usr/bin/env python3
"""Legacy pool-query builder for migration-only projects."""

from __future__ import annotations

import argparse,csv,gzip,hashlib,json
from pathlib import Path

def open_text(path,mode):return gzip.open(path,mode,newline="",encoding="utf-8") if path.suffix==".gz" else path.open(mode.replace("t",""),newline="",encoding="utf-8")
def main()->int:
    p=argparse.ArgumentParser();p.add_argument("--cell-ledger",required=True,type=Path);p.add_argument("--decision-ids",required=True,type=Path);p.add_argument("--out",required=True,type=Path);p.add_argument("--pool-id",required=True);p.add_argument("--candidate-lineages",required=True);p.add_argument("--state-tags",default="");p.add_argument("--spatial-tags",default="");a=p.parse_args();ids={x.strip() for x in a.decision_ids.read_text().splitlines() if x.strip() and not x.startswith("#")};fields=["cell_id","query_or_anchor","anchor_label","source_key","parent_decision_id","parent_pool_id","pool_snapshot_id","generation","state_tags","spatial_tags","qc_tags","candidate_lineages"]
    a.out.parent.mkdir(parents=True,exist_ok=True);tmp=a.out.with_name(a.out.stem+".tmp"+a.out.suffix);n=0;seen=set();source_counts={}
    with open_text(a.cell_ledger,"rt") as inp,open_text(tmp,"wt") as out:
        r=csv.DictReader(inp,delimiter="\t");w=csv.DictWriter(out,fieldnames=fields,delimiter="\t");w.writeheader()
        for row in r:
            if row.get("decision_id") not in ids:continue
            cid=row["cell_id"]
            if cid in seen:raise SystemExit(f"duplicate cell {cid}")
            seen.add(cid);source=row.get("source_key") or f"{row.get('source_run_id','')}|{row.get('source_cluster','')}";source_counts[source]=source_counts.get(source,0)+1
            w.writerow({"cell_id":cid,"query_or_anchor":"query","anchor_label":"","source_key":source,"parent_decision_id":row.get("decision_id",""),"parent_pool_id":row.get("parent_pool_id",""),"pool_snapshot_id":a.pool_id,"generation":str(int(float(row.get("generation",row.get("iteration",1)) or 1))+1),"state_tags":"|".join(x for x in [row.get("state_tags",""),a.state_tags] if x),"spatial_tags":"|".join(x for x in [row.get("spatial_tags",""),a.spatial_tags] if x),"qc_tags":row.get("qc_tags",""),"candidate_lineages":a.candidate_lineages});n+=1
    tmp.replace(a.out);digest=hashlib.sha256()
    with a.out.open("rb") as handle:
        for block in iter(lambda:handle.read(8*1024*1024),b""):digest.update(block)
    manifest={"status":"PASS","pool_id":a.pool_id,"decision_ids":sorted(ids),"n_query":n,"source_counts":source_counts,"membership_sha256":digest.hexdigest(),"membership":str(a.out.resolve())};Path(str(a.out)+".manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n");print(json.dumps(manifest,ensure_ascii=False,indent=2));return 0
if __name__=="__main__":raise SystemExit(main())
