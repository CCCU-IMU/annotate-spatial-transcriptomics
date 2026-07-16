#!/usr/bin/env python3
"""Block final release until all iterative biological and provenance gates close."""
from __future__ import annotations
import argparse,csv,hashlib,json,re
from pathlib import Path
from validate_direct_lineage_workflow import audit as audit_direct_lineage
from validate_incident_registry import validate as validate_incidents
from validate_profile_role import load_profile
from audit_annotation_membership_partition import audit as audit_membership_partition
from validate_annotation_support_registry import validate as validate_support_registry
from dependency_manifest import build as build_dependency_manifest
from evidence_schema_lib import active_registry_rows

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
def referenced_files(root, rows):
    """Collect hash dependencies declared by active registries and evidence JSON."""
    found=set()
    pending=[]
    for row in rows:
        for key,value in row.items():
            if not value or not (key.endswith("_path") or key.endswith("_artifact") or key in {"membership_path"}):continue
            path=Path(value);path=path if path.is_absolute() else root/path
            if path.is_file():pending.append(path.resolve())
    while pending:
        path=pending.pop()
        if path in found:continue
        found.add(path)
        if path.suffix.lower()!=".json":continue
        try:payload=json.loads(path.read_text(encoding="utf-8"))
        except (OSError,json.JSONDecodeError,UnicodeDecodeError):continue
        stack=[payload]
        while stack:
            value=stack.pop()
            if isinstance(value,dict):
                for key,item in value.items():
                    if isinstance(item,(dict,list)):stack.append(item)
                    elif isinstance(item,str) and (key=="path" or key.endswith("_path") or key.endswith("_artifact")):
                        child=Path(item);child=child if child.is_absolute() else root/child
                        if child.is_file() and child.resolve() not in found:pending.append(child.resolve())
            elif isinstance(value,list):stack.extend(value)
    return sorted(found)
def main():
    p=argparse.ArgumentParser();p.add_argument("project_root",type=Path);p.add_argument("--context",type=Path);p.add_argument("--biological-profile",type=Path);p.add_argument("--profile",type=Path,help="deprecated alias; must still be a biological_evidence profile");a=p.parse_args();r=a.project_root;errors=[];all_clusters=read(r/"state/cluster_decision_ledger.tsv");clusters=active_decisions(all_clusters);cohorts=active_registry_rows(read(r/"state/recluster_cohort_registry.tsv"),"cohort_id");returns=active_registry_rows(read(r/"state/direct_return_registry.tsv"),"return_id");runs=read(r/"state/run_registry.tsv");queue=read(r/"state/next_action_queue.tsv")
    project_path=r/"config/project.json";project=json.loads(project_path.read_text()) if project_path.exists() else {};preset_requested=project.get("strategy_preset_requested","");preset_record={}
    if preset_requested:
        preset_path=r/"config/active_strategy_preset.json"
        if not preset_path.exists():errors.append("requested strategy preset is not activated in config/active_strategy_preset.json")
        else:
            preset_record=json.loads(preset_path.read_text())
            if preset_record.get("strategy_preset_status")!="ACTIVE" or preset_record.get("strategy_preset_id")!=preset_requested:errors.append("active strategy preset does not match the requested preset")
            if preset_record.get("strategy_preset_preprocessing_mode")!="fixed_verified_same_batch_contract" or preset_record.get("fixed_cellbin_preprocessing_required") is not True:errors.append("same-batch R-first preset lacks independently verified fixed cellbin preprocessing provenance")
            bindings=preset_record.get("strategy_preset_bindings") or {}
            for key,expected in bindings.items():
                if not key.endswith("_sha256"):continue
                source_key=key[:-7];value=bindings.get(source_key,"");path=Path(value)
                if not value or not path.is_file() or sha256(path)!=expected:errors.append(f"strategy preset binding is missing or stale: {source_key}")
            discovery_audit=r/"provenance/open_world_lineage_audit_validation.json"
            if not discovery_audit.exists():errors.append("same-batch strategy preset lacks the open-world lineage discovery audit")
            else:
                discovery=json.loads(discovery_audit.read_text());ledger=r/"state/cluster_decision_ledger.tsv"
                if discovery.get("status")!="PASS" or not ledger.exists() or discovery.get("cluster_ledger_sha256")!=sha256(ledger):errors.append("open-world lineage discovery audit is failed or stale")
                for source_key,hash_key in (("audit","audit_sha256"),("candidate_catalog","candidate_catalog_sha256"),("biological_profile","biological_profile_sha256")):
                    value=str(discovery.get(source_key,"") or "");path=Path(value) if value else Path()
                    if not value or not path.is_file() or discovery.get(hash_key)!=sha256(path):errors.append(f"open-world lineage discovery binding is missing or stale: {source_key}")
    context_path=a.context or r/"config/biological_context.json";context={}
    if not context_path.exists():errors.append("missing biological context")
    else:context=json.loads(context_path.read_text())
    cv=r/"provenance/biological_context_validation.json"
    if not cv.exists() or json.loads(cv.read_text()).get("status")!="PASS":errors.append("biological-context validation has not passed")
    if queue:errors.append(f"next-action queue still contains {len(queue)} items")
    review_gate_path=r/"provenance/annotation_review_gate.json"
    if review_gate_path.exists():
        review_gate=json.loads(review_gate_path.read_text())
        if review_gate.get("status") not in {"PASS","CLOSED"}:
            errors.append(
                "annotation review gate remains open: "
                + str(review_gate.get("required_next_action", review_gate.get("reason", "unspecified review")))
            )
    profile_path=a.biological_profile or a.profile;profile={}
    if not profile_path:
        open_world_path=r/"provenance/open_world_lineage_audit_validation.json"
        if open_world_path.is_file():
            try:
                bound_profile=str(json.loads(open_world_path.read_text()).get("biological_profile","") or "")
                if bound_profile:profile_path=Path(bound_profile)
            except (OSError,json.JSONDecodeError):pass
    if not profile_path and preset_record:
        bound_profile=str((preset_record.get("strategy_preset_bindings") or {}).get("biological_profile","") or "")
        if bound_profile:profile_path=Path(bound_profile)
    direct_workflow = project.get("routing_model", "direct_cross_lineage_recluster_cohorts") in {"direct_cross_lineage_recluster_cohorts", "direct_cross_branch_recluster_cohorts", "direct_cross_lineage_recluster_cohorts_global_atlas"}
    if not direct_workflow:errors.append("legacy pool-routing project must be migrated to the direct-lineage cohort architecture before completion")
    if profile_path:
        try:profile=load_profile(profile_path,"biological_evidence")
        except Exception as exc:errors.append(f"invalid biological profile binding: {exc}")
    requires_full_feature=bool(profile.get("final_validation"))
    if requires_full_feature:
        ff=r/"provenance/full_feature_validation.json"
        if not ff.exists() or json.loads(ff.read_text()).get("status")!="PASS":errors.append("active biological profile requires full-feature validation, but the gate has not passed")
    multi=audit_direct_lineage(r)
    if preset_requested and not profile_path:errors.append("active strategy preset lacks a resolvable biological profile binding")
    if profile.get("annotation_workflow_policy",{}).get("incident_registry_required",False):
        incident_path=r/"provenance/incidents/incident_registry.tsv"
        incident=validate_incidents(incident_path)
        (r/"provenance/incident_registry_validation.json").write_text(json.dumps(incident,ensure_ascii=False,indent=2)+"\n")
        if incident.get("status")!="PASS":errors.append(f"incident registry blocks completion: {len(incident.get('open_incidents',[]))} open; {'; '.join(incident.get('errors',[]))}")
    for x in clusters:
        did=x.get("decision_id") or f"{x.get('decision_version','')}:{x.get('source_run_id','')}:{x.get('source_cluster','')}";closed=str(x.get("closed","")).lower() in {"true","1","yes"}
        artifact=x.get("validation_artifact","");ap=Path(artifact) if artifact else None;artifact_ok=bool(ap and (ap if ap.is_absolute() else r/ap).exists());evidence=(x.get("route","")+" "+x.get("evidence_status","")+" "+x.get("validation_status","")).lower()
        if not closed:errors.append(f"cluster {x.get('source_cluster')} is not closed")
        if closed and not x.get("closure_rationale"):errors.append(f"cluster {x.get('source_cluster')} is closed without rationale")
        if closed and x.get("state") in {"defined_fine","defined_broad_only","interface_review"} and not artifact_ok:errors.append(f"cluster {x.get('source_cluster')} lacks a resolvable validation artifact")
        if requires_full_feature:
            if x.get("state") in {"defined_fine","defined_broad_only"} and x.get("validation_feature_scope")!="full_feature":errors.append(f"biological cluster {x.get('source_cluster')} lacks full-feature validation scope")
            if "oocyte" in evidence and x.get("validation_feature_scope")!="full_feature":errors.append(f"Oocyte-related decision {x.get('source_cluster')} lacks full-feature validation scope")
        # Residual QC is gated once as the exact terminal membership.  A
        # cluster-level Atlas prerequisite would revive the retired QC-routing
        # model and can disagree with the final cell-ledger partition.
        if x.get("state")=="defined_broad_only" and str(x.get("fine_anchor_eligible","")).lower() in {"true","1","yes"}:errors.append(f"broad-only cluster {x.get('source_cluster')} is a fine anchor")
        if any(k in (x.get("broad_label","")+x.get("fine_label","")).lower() for k in ["oocyte","germ cell"]):
            if x.get("validation_status")!="passed":errors.append(f"rare-cell cluster {x.get('source_cluster')} lacks passed strict validation")
            if context.get("observation_unit") in {"cellbin","spot"}:
                try:object_ok=float(x.get("spatial_object_count",0) or 0)>0 and "not" in x.get("count_interpretation","").lower()
                except ValueError:object_ok=False
                if not object_ok:errors.append(f"rare-cell cluster {x.get('source_cluster')} lacks object-level count interpretation")
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
    if profile_path:
        taxonomy_path=r/"provenance/release_taxonomy_audit.json"
        cell_ledger=r/"state/cell_ledger.tsv.gz"
        if not taxonomy_path.exists():errors.append("missing provenance/release_taxonomy_audit.json")
        else:
            taxonomy=json.loads(taxonomy_path.read_text())
            if taxonomy.get("pass") is not True:errors.append("release taxonomy audit did not pass")
            try:taxonomy_metadata=Path(taxonomy.get("metadata","")).resolve()
            except Exception:taxonomy_metadata=Path()
            if taxonomy_metadata!=cell_ledger.resolve():errors.append("release taxonomy audit was not run on the current cell ledger")
            if not cell_ledger.exists() or taxonomy.get("metadata_sha256")!=sha256(cell_ledger):errors.append("release taxonomy audit is stale for state/cell_ledger.tsv.gz")
            if taxonomy.get("profile_sha256")!=sha256(profile_path):errors.append("release taxonomy audit is stale for the active biological profile")
    partition=audit_membership_partition(r)
    (r/"provenance/annotation_membership_partition_audit.json").write_text(json.dumps(partition,ensure_ascii=False,indent=2)+"\n")
    if partition.get("status")!="PASS":errors.append(f"annotation membership partition audit failed: {len(partition.get('errors',[]))} errors")
    cell_ledger=r/"state/cell_ledger.tsv.gz"
    if not cell_ledger.exists():cell_ledger=r/"state/cell_ledger.tsv"
    support=validate_support_registry(r,r/"state/annotation_support_registry.tsv",cell_ledger)
    (r/"provenance/annotation_support_validation.json").write_text(json.dumps(support,ensure_ascii=False,indent=2)+"\n")
    if support.get("status")!="PASS":errors.append(f"annotation support coverage validation failed: {len(support.get('errors',[]))} errors")
    (r/"provenance/direct_lineage_workflow_audit.json").write_text(json.dumps(multi,ensure_ascii=False,indent=2)+"\n")
    if multi.get("status")!="PASS":errors.append(f"annotation workflow completion is blocked: {len(multi.get('gaps',[]))} gaps, {len(multi.get('invalid_attempts',[]))} invalid attempts, missing views={multi.get('missing_views',[])}")
    result={"status":"PASS" if not errors else "BLOCKED","errors":errors,"active_decisions":len(clusters),"historical_decisions":len(all_clusters),"recluster_cohorts":len(cohorts),"direct_returns":len(returns),"persistent_biological_pools":False,"runs":len(runs),"annotation_workflow_status":multi.get("status"),"strategy_preset_requested":preset_requested,"strategy_preset_id":preset_record.get("strategy_preset_id"),"context":context}
    completion_path=r/"provenance/completion_gate.json";completion_path.parent.mkdir(parents=True,exist_ok=True);completion_path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    dependencies=[path for path in [cell_ledger,r/"state/cluster_decision_ledger.tsv",r/"state/recluster_cohort_registry.tsv",r/"state/direct_return_registry.tsv",r/"state/route_attempt_registry.tsv",r/"state/annotation_support_registry.tsv",r/"provenance/prelabel_broad_evidence_validation.json",r/"provenance/annotation_membership_partition_audit.json",r/"provenance/annotation_support_validation.json",r/"provenance/direct_lineage_workflow_audit.json",r/"provenance/whole_tissue_resolution_grid_validation.json",r/"config/active_workflow_profile.json",r/"config/active_strategy_preset.json",r/"provenance/input_contract_validation.json"] if path.is_file()]
    dependencies.extend(path for path in referenced_files(r,clusters+cohorts+returns+read(r/"state/route_attempt_registry.tsv")+read(r/"state/annotation_support_registry.tsv")) if path not in dependencies)
    if profile_path and profile_path.is_file() and profile_path not in dependencies:dependencies.append(profile_path)
    build_dependency_manifest(completion_path,dependencies,{"gate":"completion"})
    print(json.dumps(result,ensure_ascii=False,indent=2));return 0 if not errors else 2
if __name__=="__main__":raise SystemExit(main())
