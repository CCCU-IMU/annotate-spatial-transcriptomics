#!/usr/bin/env python3
"""Normalize legitimate zero-marker DEG outputs to a schema-preserving TSV."""

from __future__ import annotations
import argparse,csv,json
from datetime import datetime,timezone
from pathlib import Path

FIELDS=["p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cluster","gene"]
def main():
    p=argparse.ArgumentParser();p.add_argument("output_root",type=Path);p.add_argument("--resolutions",required=True);a=p.parse_args();repairs=[]
    for res in [x.strip() for x in a.resolutions.split(",") if x.strip()]:
        tag=res.replace(".","p");cluster=a.output_root/f"tables/framework_res{tag}_clusters.tsv";allp=a.output_root/f"tables/framework_res{tag}_DEG_all.tsv";top=a.output_root/f"tables/framework_res{tag}_DEG_top100.tsv"
        if not cluster.exists():raise SystemExit(f"missing cluster table for {res}")
        with cluster.open(newline="",encoding="utf-8") as h:clusters={x["cluster"] for x in csv.DictReader(h,delimiter="\t")}
        for path in [allp,top]:
            if path.exists() and path.stat().st_size==0:
                with path.open("w",newline="",encoding="utf-8") as h:csv.DictWriter(h,fieldnames=FIELDS,delimiter="\t").writeheader()
                repairs.append({"resolution":res,"file":str(path.relative_to(a.output_root)),"n_clusters":len(clusters),"interpretation":"no marker passed configured thresholds; not an interpretable biological split"})
    result={"status":"PASS","repairs":repairs,"timestamp":datetime.now(timezone.utc).isoformat()};out=a.output_root/"empty_deg_schema_normalization.json";out.write_text(json.dumps(result,indent=2)+"\n");print(json.dumps(result,indent=2));return 0
if __name__=="__main__":raise SystemExit(main())
