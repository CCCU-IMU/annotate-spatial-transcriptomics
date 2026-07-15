#!/usr/bin/env python3
"""Create a fail-closed next-action queue from cluster/pool state and biological context."""
from __future__ import annotations
import argparse,csv,json,re
from pathlib import Path
from multiroute_lib import audit_multiroute
from validate_profile_role import load_profile

def read(path):
    if not path.exists(): return []
    with path.open(newline="",encoding="utf-8") as h:return list(csv.DictReader(h,delimiter="\t"))

def active_decisions(rows):
    def rid(r): return r.get("decision_id") or f"{r.get('decision_version','')}:{r.get('source_run_id','')}:{r.get('source_cluster','')}"
    superseded=set()
    for r in rows:
        superseded.update(x for x in re.split(r"[;,\s]+",r.get("supersedes","").strip()) if x)
    return [r for r in rows if rid(r) not in superseded]

def main():
    p=argparse.ArgumentParser();p.add_argument("project_root",type=Path);p.add_argument("--context",required=True,type=Path);p.add_argument("--biological-profile",type=Path);p.add_argument("--profile",type=Path,help="deprecated biological-profile alias");a=p.parse_args();root=a.project_root;ctx=json.loads(a.context.read_text());profile_path=a.biological_profile or a.profile;profile=load_profile(profile_path,"biological_evidence") if profile_path else {};all_clusters=read(root/"state/cluster_decision_ledger.tsv");clusters=active_decisions(all_clusters);pools=read(root/"state/pool_registry.tsv");pool_status={r.get("pool_id"):r.get("status") for r in pools};rows=[]
    context_rules=profile.get("context_specific_identity_rules",profile.get("rare_cell_rules",{}));priority_terms=[str(x).lower() for x in ctx.get("priority_lineages",[])];multiroute=audit_multiroute(root,ctx,profile);gap_routes={}
    for gap in multiroute.get("gaps",[]):gap_routes.setdefault(gap.get("decision_id",""),set()).add(gap.get("required_route",""))
    if profile.get("final_validation"):
        ff=root/"provenance/full_feature_validation.json";ok=ff.exists() and json.loads(ff.read_text()).get("status")=="PASS"
        if not ok:rows.append({"priority":0,"source_run_id":"project","source_cluster":"__FULL_FEATURE_VALIDATION__","n_observations":0,"current_state":"project_gate","broad_label":"","fine_label":"","required_route":"full_feature_marker_validation","reason":"HVG/reduced features cannot provide final positive/negative marker evidence","target_pool":"whole_tissue_and_context_gated_pools","blocked_until":"full-feature audit PASS + final lineage/rare-cell evidence writeback"})
    for r in clusters:
        did=r.get("decision_id") or f"{r.get('decision_version','')}:{r.get('source_run_id','')}:{r.get('source_cluster','')}";n=int(float(r.get("n_observations",0) or 0));state=r.get("state","");broad=r.get("broad_label","");fine=r.get("fine_label","");text=(broad+" "+fine).lower();action=r.get("next_action","");evidence=(r.get("route","")+" "+r.get("evidence_status","")+" "+r.get("validation_status","")).lower();artifact=r.get("validation_artifact","");ap=Path(artifact) if artifact else None;artifact_ok=bool(ap and (ap if ap.is_absolute() else root/ap).exists());feature_scope=r.get("validation_feature_scope","")
        closed=str(r.get("closed","")).lower() in {"true","1","yes"};closure=bool(r.get("closure_rationale",""));iteration=int(float(r.get("iteration",1) or 1));validated=r.get("validation_status","")=="passed";object_ok=True
        if ctx.get("observation_unit") in {"cellbin","spot"}:
            try:object_ok=float(r.get("spatial_object_count",0) or 0)>0 and "not" in r.get("count_interpretation","").lower()
            except ValueError:object_ok=False
        if state=="defined_broad_only" and not (closed and closure): route="broad_pool_anchor_recluster"; reason="broad-only labels require pool-level subtype review or explicit low-priority closure"
        elif state=="interface_review" and not (closed and closure and artifact_ok): route="targeted_interface_recluster_then_optional_RCTD";reason="mixed/interface state requires a targeted-review artifact before retained closure"
        elif state=="qc_holdout" and "zero_count" not in evidence and "qc_anchor_recluster" in gap_routes.get(did,set()): route="qc_anchor_recluster";reason="the complete frozen QC holdout must undergo balanced-anchor reclustering before any Atlas route"
        elif state=="qc_holdout" and "zero_count" not in evidence and "qc_atlas_review" in gap_routes.get(did,set()): route="depth_matched_atlas_internal_anchor_review";reason="only the residual QC child left after the complete QC-anchor route may enter calibrated Atlas/internal-anchor review"
        elif state=="qc_holdout" and not (closed and closure): route="qc_anchor_recluster";reason="QC holdout must first undergo complete anchor reclustering; only its residual child may later reach Atlas"
        elif state in {"pending_review","technical_state"} and not (closed and closure): route="targeted_review_or_explicit_retention";reason="unresolved state needs evidence-backed closure"
        else: route="none";reason="currently defined"
        for name,rule in context_rules.items():
            if any(re.search(pat,text,re.I) for pat in rule.get("label_patterns",[])):
                if not (validated and closed and closure and artifact_ok and feature_scope=="full_feature" and object_ok): route="context_specific_identity_validation";reason=f"{name} requires full-feature core/negative/spatial/recluster validation and object-level counting"
                break
        if route=="none" and iteration<2 and any(term in text for term in priority_terms):
            route="priority_lineage_pool_recluster";reason="priority biological lineage requires at least one pool-level adaptive reclustering/validation round"
        if route=="none" and state in {"defined_fine","defined_broad_only"} and profile.get("final_validation") and feature_scope!="full_feature":
            route="full_feature_marker_validation";reason="final biological labels require full-feature positive/negative marker validation"
        if route=="none" and "oocyte" in evidence and profile.get("final_validation") and feature_scope!="full_feature":
            route="full_feature_marker_validation";reason="context-gated rare/interface decision was made from reduced features"
        if route=="none" and (not closed or not closure or (state in {"defined_fine","defined_broad_only"} and not artifact_ok)):
            route="closure_evidence_review";reason="decision lacks an explicit rationale or a resolvable evidence artifact"
        if route!="none":
            rows.append({"priority":1 if route=="context_specific_identity_validation" else 2 if "interface" in route else 3,"source_run_id":r.get("source_run_id",""),"source_cluster":r.get("source_cluster",""),"n_observations":n,"current_state":state,"broad_label":broad,"fine_label":fine,"required_route":route,"reason":reason,"target_pool":r.get("target_pool",action),"blocked_until":"validation artifact + state writeback + pool closure"})
    (root/"provenance/multiroute_audit.json").write_text(json.dumps(multiroute,ensure_ascii=False,indent=2)+"\n")
    existing={(x.get("source_cluster"),x.get("required_route")) for x in rows}
    decision_lookup={r.get("decision_id") or f"{r.get('decision_version','')}:{r.get('source_run_id','')}:{r.get('source_cluster','')}":r for r in clusters}
    for gap in multiroute.get("gaps",[]):
        decision=decision_lookup.get(gap["decision_id"],{});key=(decision.get("source_cluster",gap["decision_id"]),gap["required_route"])
        if key in existing:continue
        rows.append({"priority":1 if gap["required_route"]=="context_specific_identity_review" else 2,"source_run_id":decision.get("source_run_id",""),"source_cluster":decision.get("source_cluster",gap["decision_id"]),"n_observations":gap["n_observations"],"current_state":decision.get("state","multiroute_gap"),"broad_label":decision.get("broad_label",""),"fine_label":decision.get("fine_label",""),"required_route":gap["required_route"],"reason":gap["reason"],"target_pool":decision.get("target_pool",""),"blocked_until":"validated route attempt + atomic writeback"})
    for view in multiroute.get("missing_views",[]):rows.append({"priority":4,"source_run_id":"project","source_cluster":f"__ANNOTATION_VIEW_{view.upper()}__","n_observations":0,"current_state":"project_gate","broad_label":"","fine_label":"","required_route":"build_final_annotation","reason":f"missing validated {view} annotation view","target_pool":"analysis_set","blocked_until":"single final census + cell-level fields + artifact validation"})
    for invalid in multiroute.get("invalid_attempts",[]):rows.append({"priority":1,"source_run_id":"project","source_cluster":invalid.get("route_attempt_id","__INVALID_ROUTE__"),"n_observations":0,"current_state":"invalid_route_attempt","broad_label":"","fine_label":"","required_route":"repair_route_attempt","reason":"; ".join(invalid.get("errors",[])),"target_pool":"route_attempt_registry","blocked_until":"valid terminal route record"})
    rows.sort(key=lambda x:(x["priority"],-x["n_observations"]));out=root/"state/next_action_queue.tsv";out.parent.mkdir(parents=True,exist_ok=True)
    fields=list(rows[0]) if rows else ["priority","source_cluster","required_route"]
    with out.open("w",newline="",encoding="utf-8") as h:w=csv.DictWriter(h,fieldnames=fields,delimiter="\t");w.writeheader();w.writerows(rows)
    action_observations=sum(r["n_observations"] for r in rows);unique_decisions={}
    for item in rows:
        key=(item.get("source_run_id",""),item.get("source_cluster",""));unique_decisions[key]=max(unique_decisions.get(key,0),int(item.get("n_observations",0) or 0))
    result={"status":"ITERATION_REQUIRED" if rows else "READY_FOR_COMPLETION_AUDIT","queued_actions":len(rows),"queued_observations":sum(unique_decisions.values()),"queued_action_observations":action_observations,"queued_observation_note":"queued_observations counts each active decision once; queued_action_observations can exceed sample size when one decision requires multiple routes","active_decisions":len(clusters),"historical_decisions":len(all_clusters),"multiroute_status":multiroute.get("status"),"context":ctx}
    (root/"provenance/iteration_plan.json").parent.mkdir(parents=True,exist_ok=True);(root/"provenance/iteration_plan.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n");print(json.dumps(result,ensure_ascii=False,indent=2));return 2 if rows else 0
if __name__=="__main__":raise SystemExit(main())
