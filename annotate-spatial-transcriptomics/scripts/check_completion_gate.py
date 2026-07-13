#!/usr/bin/env python3
"""Block final release until all iterative biological and provenance gates close."""
from __future__ import annotations
import argparse,csv,hashlib,json,re
from pathlib import Path
from multiroute_lib import audit_multiroute

def read(path):
    if not path.exists():return []
    with path.open(newline="",encoding="utf-8") as h:return list(csv.DictReader(h,delimiter="\t"))
def active_decisions(rows):
    def rid(x):return x.get("decision_id") or f"{x.get('decision_version','')}:{x.get('source_run_id','')}:{x.get('source_cluster','')}"
    superseded=set()
    for x in rows:superseded.update(v for v in re.split(r"[;,\s]+",x.get("supersedes","").strip()) if v)
    return [x for x in rows if rid(x) not in superseded]
def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument("project_root",type=Path);p.add_argument("--context",type=Path);p.add_argument("--profile",type=Path);a=p.parse_args();r=a.project_root;errors=[];all_clusters=read(r/"state/cluster_decision_ledger.tsv");clusters=active_decisions(all_clusters);pools=read(r/"state/pool_registry.tsv");runs=read(r/"state/run_registry.tsv");queue=read(r/"state/next_action_queue.tsv")
    context_path=a.context or r/"config/biological_context.json";context={}
    if not context_path.exists():errors.append("missing biological context")
    else:context=json.loads(context_path.read_text())
    cv=r/"provenance/biological_context_validation.json"
    if not cv.exists() or json.loads(cv.read_text()).get("status")!="PASS":errors.append("biological-context validation has not passed")
    if context.get("profile"):
        ff=r/"provenance/full_feature_validation.json"
        if not ff.exists() or json.loads(ff.read_text()).get("status")!="PASS":errors.append("full-feature marker validation has not passed")
    if queue:errors.append(f"next-action queue still contains {len(queue)} items")
    review_gate_path=r/"provenance/annotation_review_gate.json"
    if review_gate_path.exists():
        review_gate=json.loads(review_gate_path.read_text())
        if review_gate.get("status") not in {"PASS","CLOSED"}:
            errors.append(
                "annotation review gate remains open: "
                + str(review_gate.get("required_next_action", review_gate.get("reason", "unspecified review")))
            )
    if a.profile:
        profile=json.loads(a.profile.read_text());multi=audit_multiroute(r,context,profile);gap_routes={}
        for gap in multi.get("gaps",[]):gap_routes.setdefault(gap.get("decision_id",""),set()).add(gap.get("required_route",""))
    else:multi={"status":"NOT_RUN"};gap_routes={}
    for x in clusters:
        did=x.get("decision_id") or f"{x.get('decision_version','')}:{x.get('source_run_id','')}:{x.get('source_cluster','')}";closed=str(x.get("closed","")).lower() in {"true","1","yes"}
        artifact=x.get("validation_artifact","");ap=Path(artifact) if artifact else None;artifact_ok=bool(ap and (ap if ap.is_absolute() else r/ap).exists());evidence=(x.get("route","")+" "+x.get("evidence_status","")+" "+x.get("validation_status","")).lower()
        if not closed:errors.append(f"cluster {x.get('source_cluster')} is not closed")
        if closed and not x.get("closure_rationale"):errors.append(f"cluster {x.get('source_cluster')} is closed without rationale")
        if closed and x.get("state") in {"defined_fine","defined_broad_only","interface_review"} and not artifact_ok:errors.append(f"cluster {x.get('source_cluster')} lacks a resolvable validation artifact")
        if context.get("profile"):
            if x.get("state") in {"defined_fine","defined_broad_only"} and x.get("validation_feature_scope")!="full_feature":errors.append(f"biological cluster {x.get('source_cluster')} lacks full-feature validation scope")
            if "oocyte" in evidence and x.get("validation_feature_scope")!="full_feature":errors.append(f"Oocyte-related decision {x.get('source_cluster')} lacks full-feature validation scope")
        if closed and x.get("state")=="qc_holdout" and "zero_count" not in evidence and "qc_atlas_review" in gap_routes.get(did,set()):errors.append(f"non-zero QC cluster {x.get('source_cluster')} lacks a calibrated atlas/internal-anchor attempt")
        if x.get("state")=="defined_broad_only" and str(x.get("fine_anchor_eligible","")).lower() in {"true","1","yes"}:errors.append(f"broad-only cluster {x.get('source_cluster')} is a fine anchor")
        if any(k in (x.get("broad_label","")+x.get("fine_label","")).lower() for k in ["oocyte","germ cell"]):
            if x.get("validation_status")!="passed":errors.append(f"rare-cell cluster {x.get('source_cluster')} lacks passed strict validation")
            if context.get("observation_unit") in {"cellbin","spot"}:
                try:object_ok=float(x.get("spatial_object_count",0) or 0)>0 and "not" in x.get("count_interpretation","").lower()
                except ValueError:object_ok=False
                if not object_ok:errors.append(f"rare-cell cluster {x.get('source_cluster')} lacks object-level count interpretation")
    for x in pools:
        if x.get("status") not in {"closed_and_frozen","explicitly_retained_closed"}:errors.append(f"pool {x.get('pool_id')} remains {x.get('status')}")
    for x in runs:
        if x.get("status") not in {"validated_done","failed_preserved","cancelled_preserved"}:errors.append(f"run {x.get('run_id')} remains {x.get('status')}")
    required=["provenance/state_validation.json"]
    for f in required:
        if not (r/f).exists():errors.append(f"missing {f}")
    sv=r/"provenance/state_validation.json"
    if sv.exists():
        val=json.loads(sv.read_text())
        if val.get("status")!="PASS":errors.append("state validation did not pass")
        validated=val.get("validated_files",{})
        if "state/cluster_decision_ledger.tsv" not in validated:errors.append("state validation did not hash cluster decision ledger")
        if not any(k.endswith("state/cell_ledger.tsv.gz") or k.endswith("state/cell_ledger.tsv") for k in validated):errors.append("state validation did not hash current cell ledger")
        for rel,expected in validated.items():
            pth=r/rel
            if not pth.exists() or sha256(pth)!=expected:errors.append(f"state validation is stale for {rel}")
    if a.profile:
        (r/"provenance/multiroute_audit.json").write_text(json.dumps(multi,ensure_ascii=False,indent=2)+"\n")
        if multi.get("status")!="PASS":errors.append(f"multi-route completion is blocked: {len(multi.get('gaps',[]))} route gaps, {len(multi.get('invalid_attempts',[]))} invalid attempts, missing views={multi.get('missing_views',[])}")
    result={"status":"PASS" if not errors else "BLOCKED","errors":errors,"active_decisions":len(clusters),"historical_decisions":len(all_clusters),"pools":len(pools),"runs":len(runs),"multiroute_status":multi.get("status"),"context":context};(r/"provenance/completion_gate.json").parent.mkdir(parents=True,exist_ok=True);(r/"provenance/completion_gate.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n");print(json.dumps(result,ensure_ascii=False,indent=2));return 0 if not errors else 2
if __name__=="__main__":raise SystemExit(main())
