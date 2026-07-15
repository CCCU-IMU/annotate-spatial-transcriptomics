#!/usr/bin/env python3
"""Register one immutable broad-class or targeted reclustering cohort."""

from __future__ import annotations
import argparse, csv, hashlib, json
from datetime import datetime, timezone
from pathlib import Path

FIELDS = ["cohort_id","sample_id","cohort_type","source_broad_label","source_run_ids","source_cluster_ids","membership_path","membership_sha256","n_observations","purpose","competing_hypotheses","candidate_resolutions","selected_resolution","applicability","applicability_rationale","outcome_artifact","status","created_at","closed_at"]

def digest(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(8*1024*1024),b""): h.update(block)
    return h.hexdigest()

def main() -> int:
    p=argparse.ArgumentParser();p.add_argument("project_root",type=Path);p.add_argument("--cohort-id",required=True);p.add_argument("--cohort-type",required=True,choices=["broad_class_recluster","targeted_recluster","oocyte_targeted_recluster"]);p.add_argument("--source-broad-label",default="");p.add_argument("--source-run-ids",required=True);p.add_argument("--source-cluster-ids",required=True);p.add_argument("--membership",type=Path);p.add_argument("--purpose",required=True);p.add_argument("--competing-hypotheses",default="");p.add_argument("--candidate-resolutions",default="");p.add_argument("--selected-resolution",default="");p.add_argument("--applicability",choices=["applicable","not_applicable"],required=True);p.add_argument("--applicability-rationale",required=True);p.add_argument("--outcome-artifact",type=Path);p.add_argument("--status",choices=["validated_done","not_applicable_reviewed"],required=True);a=p.parse_args()
    root=a.project_root.resolve();cfg=json.loads((root/"config/project.json").read_text());membership_path="";membership_sha="";n=0
    if a.applicability=="applicable":
        if not a.membership or not a.membership.is_file(): raise SystemExit("applicable cohort requires a membership artifact")
        with a.membership.open(newline="",encoding="utf-8") as h: rows=list(csv.DictReader(h,delimiter="\t"))
        ids=[r.get("cell_id","") for r in rows]
        if not rows or "cell_id" not in rows[0] or "" in ids or len(ids)!=len(set(ids)): raise SystemExit("membership must contain unique nonempty cell_id")
        if not a.outcome_artifact or not a.outcome_artifact.is_file(): raise SystemExit("applicable cohort requires a validated outcome artifact")
        membership_path=str(a.membership.resolve());membership_sha=digest(a.membership);n=len(ids)
    elif a.status!="not_applicable_reviewed" or len(a.applicability_rationale.strip())<20:
        raise SystemExit("not-applicable cohort requires reviewed terminal status and substantive rationale")
    now=datetime.now(timezone.utc).isoformat();registry=root/"state/recluster_cohort_registry.tsv";old=[]
    if registry.is_file():
        with registry.open(newline="",encoding="utf-8") as h: old=list(csv.DictReader(h,delimiter="\t"))
    existing=[r for r in old if r.get("cohort_id")==a.cohort_id]
    if existing: raise SystemExit("cohort_id already exists; create a versioned successor")
    row={"cohort_id":a.cohort_id,"sample_id":cfg["sample_id"],"cohort_type":a.cohort_type,"source_broad_label":a.source_broad_label,"source_run_ids":a.source_run_ids,"source_cluster_ids":a.source_cluster_ids,"membership_path":membership_path,"membership_sha256":membership_sha,"n_observations":n,"purpose":a.purpose,"competing_hypotheses":a.competing_hypotheses,"candidate_resolutions":a.candidate_resolutions,"selected_resolution":a.selected_resolution,"applicability":a.applicability,"applicability_rationale":a.applicability_rationale,"outcome_artifact":str(a.outcome_artifact.resolve()) if a.outcome_artifact else "","status":a.status,"created_at":now,"closed_at":now}
    with registry.open("w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=FIELDS,delimiter="\t",extrasaction="ignore");w.writeheader();w.writerows(old+[row])
    print(json.dumps(row,indent=2));return 0
if __name__=="__main__": raise SystemExit(main())
