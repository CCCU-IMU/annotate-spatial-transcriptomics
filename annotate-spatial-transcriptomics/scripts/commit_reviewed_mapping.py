#!/usr/bin/env python3
"""Commit reviewed cluster mappings and candidate decisions into project state."""

from __future__ import annotations
import argparse,csv
from datetime import datetime,timezone
from pathlib import Path

ALLOWED_STATES = {"defined_fine", "defined_broad_only", "interface_review", "qc_holdout", "technical_state", "pending_review", "unknown_candidate", "excluded_initial_qc", "closed_and_frozen"}

def canonical_state(row):
    state = row.get("state", "").strip()
    aliases = {
        "technical_or_low_information": "technical_state",
        "candidate_defined_pending_context_safe_recluster": "pending_review",
    }
    if state == "defined":
        state = "defined_fine" if row.get("fine_label", "").strip() else "defined_broad_only"
    state = aliases.get(state, state)
    if state not in ALLOWED_STATES:
        raise SystemExit(f"invalid annotation state for cluster {row.get('source_cluster', '')}: {row.get('state', '')}")
    return state

def read(path,delimiter="\t"):
    with path.open(newline="",encoding="utf-8") as h:return list(csv.DictReader(h,delimiter=delimiter))

def append_rows(path,fields,new_rows):
    old=read(path) if path.exists() and path.stat().st_size else []
    with path.open("w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter="\t",extrasaction="ignore");w.writeheader()
        for r in old:w.writerow({k:r.get(k,"") for k in fields})
        for r in new_rows:w.writerow(r)

def main():
    p=argparse.ArgumentParser();p.add_argument("project_root",type=Path);p.add_argument("--mapping",required=True,type=Path);p.add_argument("--counts",required=True,type=Path);p.add_argument("--selected-run",required=True);p.add_argument("--method",required=True);p.add_argument("--ranking",type=Path);p.add_argument("--sample",required=True);p.add_argument("--decision-version",default="v001");p.add_argument("--selection-rationale",required=True);a=p.parse_args()
    state=a.project_root/"state";state.mkdir(parents=True,exist_ok=True);now=datetime.now(timezone.utc).isoformat();maps=read(a.mapping);counts={str(r.get("cluster",r.get("source_cluster"))):r.get("n_observations",r.get("n_cells","")) for r in read(a.counts)}
    cfields=["decision_version","decision_id","sample_id","source_run_id","source_cluster","n_observations","spatial_object_count","count_interpretation","prelabel_evidence_artifact","prelabel_evidence_sha256","prelabel_evidence_frozen","prelabel_winner","prelabel_runner_up","prelabel_winning_margin","initial_broad_label","broad_label","fine_label","state","confidence","evidence_status","validation_status","validation_artifact","validation_feature_scope","route","route_run_id","recluster_cohort_id","assignment_mode","cross_lineage_target","iteration","fine_anchor_eligible","next_action","closure_rationale","supersedes","closed","created_at"]
    existing=read(state/"cluster_decision_ledger.tsv") if (state/"cluster_decision_ledger.tsv").exists() else [];known={r.get("decision_id") for r in existing if r.get("decision_id")};new=[]
    for r in maps:
        did=r.get("decision_id") or f"{a.decision_version}:{a.selected_run}:{r['source_cluster']}"
        if did in known:raise SystemExit(f"decision_id already exists: {did}")
        known.add(did);new.append({"decision_version":a.decision_version,"decision_id":did,"sample_id":a.sample,"source_run_id":a.selected_run,"source_cluster":r["source_cluster"],"n_observations":counts.get(str(r["source_cluster"]),""),"spatial_object_count":r.get("spatial_object_count",""),"count_interpretation":r.get("count_interpretation","observations_not_inferred_cells"),"prelabel_evidence_artifact":r.get("prelabel_evidence_artifact",""),"prelabel_evidence_sha256":r.get("prelabel_evidence_sha256",""),"prelabel_evidence_frozen":r.get("prelabel_evidence_frozen","false"),"prelabel_winner":r.get("prelabel_winner",r.get("broad_label","")),"prelabel_runner_up":r.get("prelabel_runner_up",""),"prelabel_winning_margin":r.get("prelabel_winning_margin",""),"initial_broad_label":r.get("initial_broad_label",r["broad_label"]),"broad_label":r["broad_label"],"fine_label":r["fine_label"],"state":canonical_state(r),"confidence":r["confidence"],"evidence_status":r["evidence_status"],"validation_status":r.get("validation_status","not_required"),"validation_artifact":r.get("validation_artifact",""),"validation_feature_scope":r.get("validation_feature_scope","unknown"),"route":r["route"],"route_run_id":r.get("route_run_id",""),"recluster_cohort_id":r.get("recluster_cohort_id",""),"assignment_mode":r.get("assignment_mode","initial_broad_direct"),"cross_lineage_target":r.get("cross_lineage_target",r.get("cross_branch_target","")),"iteration":r.get("iteration","1"),"fine_anchor_eligible":r["fine_anchor_eligible"],"next_action":r["next_action"],"closure_rationale":r.get("closure_rationale",""),"supersedes":r.get("supersedes",""),"closed":r["closed"],"created_at":now})
    append_rows(state/"cluster_decision_ledger.tsv",cfields,new)
    dfields=["decision_version","sample_id","method","run_id","parameters","n_clusters","quantitative_rank","marker_review","spatial_review","decision","rationale","created_at"]
    ranks=read(a.ranking) if a.ranking else []
    drows=[]
    for r in ranks:
        run=r["run_id"];drows.append({"decision_version":a.decision_version,"sample_id":a.sample,"method":a.method,"run_id":run,"parameters":f"resolution={r.get('resolution','')};k={r.get('k_neighbors','')}","n_clusters":r.get("n_clusters",""),"quantitative_rank":r.get("quantitative_rank",""),"marker_review":"reviewed" if run==a.selected_run else "candidate_screen","spatial_review":"reviewed" if run==a.selected_run else "candidate_screen","decision":"selected" if run==a.selected_run else "rejected_or_deferred","rationale":a.selection_rationale if run==a.selected_run else "not selected after comparative shortlist review","created_at":now})
    append_rows(state/"clustering_decision_ledger.tsv",dfields,drows)
    print(f"committed {len(maps)} clusters")
if __name__=="__main__":main()
