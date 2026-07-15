#!/usr/bin/env python3
"""Register a direct parent-broad, fine or cross-lineage observation return."""

from __future__ import annotations
import argparse,csv,hashlib,json
from datetime import datetime,timezone
from pathlib import Path

FIELDS=["return_id","sample_id","source_cohort_id","source_run_id","source_cluster","membership_path","membership_sha256","n_observations","target_broad_label","target_fine_label","confidence","assignment_mode","rctd_tier","independent_evidence","evidence_artifact","fine_anchor_eligible","status","created_at"]
def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(8*1024*1024),b""):h.update(block)
    return h.hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument("project_root",type=Path);p.add_argument("--return-id",required=True);p.add_argument("--source-cohort-id",required=True);p.add_argument("--source-run-id",required=True);p.add_argument("--source-cluster",required=True);p.add_argument("--membership",required=True,type=Path);p.add_argument("--target-broad-label",required=True);p.add_argument("--target-fine-label",default="");p.add_argument("--confidence",required=True);p.add_argument("--assignment-mode",required=True,choices=["parent_broad_direct","fine_direct","cross_lineage_direct","rctd_assisted","atlas_broad_rescue"]);p.add_argument("--rctd-tier",default="");p.add_argument("--independent-evidence",default="false");p.add_argument("--evidence-artifact",required=True,type=Path);p.add_argument("--fine-anchor-eligible",default="false");a=p.parse_args();root=a.project_root.resolve()
    if not a.membership.is_file() or not a.evidence_artifact.is_file():raise SystemExit("membership and evidence artifacts must exist")
    with a.membership.open(newline="",encoding="utf-8") as h:members=list(csv.DictReader(h,delimiter="\t"))
    ids=[r.get("cell_id","") for r in members]
    if not members or "cell_id" not in members[0] or "" in ids or len(ids)!=len(set(ids)):raise SystemExit("membership must contain unique nonempty cell_id")
    cfg=json.loads((root/"config/project.json").read_text());registry=root/"state/direct_return_registry.tsv";old=[]
    if registry.is_file():
        with registry.open(newline="",encoding="utf-8") as h:old=list(csv.DictReader(h,delimiter="\t"))
    if any(r.get("return_id")==a.return_id for r in old):raise SystemExit("return_id already exists; create a versioned successor")
    row={"return_id":a.return_id,"sample_id":cfg["sample_id"],"source_cohort_id":a.source_cohort_id,"source_run_id":a.source_run_id,"source_cluster":a.source_cluster,"membership_path":str(a.membership.resolve()),"membership_sha256":digest(a.membership),"n_observations":len(ids),"target_broad_label":a.target_broad_label,"target_fine_label":a.target_fine_label,"confidence":a.confidence,"assignment_mode":a.assignment_mode,"rctd_tier":a.rctd_tier,"independent_evidence":a.independent_evidence,"evidence_artifact":str(a.evidence_artifact.resolve()),"fine_anchor_eligible":a.fine_anchor_eligible,"status":"validated_done","created_at":datetime.now(timezone.utc).isoformat()}
    with registry.open("w",newline="",encoding="utf-8") as h:w=csv.DictWriter(h,fieldnames=FIELDS,delimiter="\t");w.writeheader();w.writerows(old+[row])
    print(json.dumps(row,indent=2));return 0
if __name__=="__main__":raise SystemExit(main())
