#!/usr/bin/env python3
"""Append or update a scheduler/local analysis run without losing failed history."""
from __future__ import annotations
import argparse,csv
from datetime import datetime,timezone
from pathlib import Path
from registry_io import locked_tsv_update

FIELDS=["run_id","sample_id","stage","script","parameters_path","environment","scheduler_job_id","status","output_root","started_at","finished_at"]
TERMINAL={"validated_done","failed_preserved","cancelled_preserved"}
def main():
    p=argparse.ArgumentParser();p.add_argument("project_root",type=Path);p.add_argument("--run-id",required=True);p.add_argument("--stage",required=True);p.add_argument("--script",required=True);p.add_argument("--parameters-path",default="");p.add_argument("--environment",required=True);p.add_argument("--job-id",default="");p.add_argument("--status",required=True,choices=["prepared","submitted","running","failed_preserved","validated_done","cancelled_preserved"]);p.add_argument("--output-root",required=True);a=p.parse_args()
    path=a.project_root/"state/run_registry.tsv";cfg=__import__("json").loads((a.project_root/"config/project.json").read_text());now=datetime.now(timezone.utc).isoformat()
    def mutate(rows,fields):
        hit=[i for i,r in enumerate(rows) if r.get("run_id")==a.run_id]
        if len(hit)>1:raise SystemExit("duplicate run_id in registry")
        old=rows[hit[0]] if hit else {}
        if old.get("status") in TERMINAL and old.get("status")!=a.status:raise SystemExit("terminal run is immutable; create a superseding run_id")
        row={"run_id":a.run_id,"sample_id":cfg["sample_id"],"stage":a.stage,"script":a.script,"parameters_path":a.parameters_path,"environment":a.environment,"scheduler_job_id":a.job_id or old.get("scheduler_job_id",""),"status":a.status,"output_root":a.output_root,"started_at":old.get("started_at") or now,"finished_at":now if a.status in TERMINAL else ""}
        if hit:rows[hit[0]]=row
        else:rows.append(row)
        return rows
    locked_tsv_update(path,mutate)
    print(f"{a.run_id}\t{a.status}")
if __name__=="__main__":main()
