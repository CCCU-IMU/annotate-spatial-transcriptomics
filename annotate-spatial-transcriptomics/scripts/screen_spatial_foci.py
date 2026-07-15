#!/usr/bin/env python3
"""Add spatial-focus evidence to a rare-cell marker screen without assigning final labels."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors

def main():
    p=argparse.ArgumentParser();p.add_argument("--screen",required=True,type=Path);p.add_argument("--coordinates",required=True,type=Path);p.add_argument("--out",required=True,type=Path);p.add_argument("--cell-id-col",default="cell_id");p.add_argument("--x-col",default="sdimx");p.add_argument("--y-col",default="sdimy");p.add_argument("--eps",type=float);p.add_argument("--min-samples",type=int,default=2);p.add_argument("--max-contradictory-hits",type=int);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    s=pd.read_csv(a.screen,sep="\t",dtype={a.cell_id_col:str})
    # The program-screen membership may already carry coordinates.  Avoid a
    # merge that silently creates ``sdimx_x/sdimx_y`` and then loses the
    # requested column names.
    if a.x_col in s.columns and a.y_col in s.columns:
        d=s
    else:
        c=pd.read_csv(a.coordinates,sep=None,engine="python",dtype={a.cell_id_col:str},usecols=[a.cell_id_col,a.x_col,a.y_col])
        d=s.merge(c,on=a.cell_id_col,how="left",validate="one_to_one")
    if d[[a.x_col,a.y_col]].isna().any().any():
        raise ValueError("Missing coordinates after cell-ID alignment")
    gate=d.starting_marker_gate.astype(str).str.lower().isin(["true","1"]);xy=d.loc[gate,[a.x_col,a.y_col]].to_numpy(float)
    if len(xy)>=2:
        if a.eps is None:
            nn=NearestNeighbors(n_neighbors=2).fit(xy);dist=nn.kneighbors(xy,return_distance=True)[0][:,1];eps=float(max(np.median(dist)*3,np.quantile(dist,.75)))
        else:eps=a.eps
        lab=DBSCAN(eps=eps,min_samples=a.min_samples).fit_predict(xy)
        d["spatial_focus_id"]=-2;d.loc[gate,"spatial_focus_id"]=lab
    else:eps=a.eps or 0.0;d["spatial_focus_id"]=-2
    d["spatial_focus_supported"]=gate & d.spatial_focus_id.ge(0)
    if a.max_contradictory_hits is None:d["contradiction_gate_review"]=True
    else:d["contradiction_gate_review"]=d.contradictory_somatic_hits.le(a.max_contradictory_hits)
    d["strict_candidate_for_recluster"]=d.spatial_focus_supported & d.contradiction_gate_review
    seed_foci=set(d.loc[d.strict_candidate_for_recluster & d.spatial_focus_id.ge(0),"spatial_focus_id"])
    # Preserve the complete multi-bin spatial object around a strict seed.  The
    # seed is evidence; the expanded object remains a review/recluster pool and
    # is not a final rare-cell call.
    d["strict_focus_for_recluster"]=d.spatial_focus_id.isin(seed_foci) & gate
    # R-first two-tier route: spatial foci are supporting evidence.  The full
    # multi-module starting gate remains the canonical query-only recluster
    # pool, including isolated high-evidence candidates.
    d["full_candidate_recluster_member"]=gate
    out=a.out/"rare_cell_spatial_focus_screen.tsv.gz";d.to_csv(out,sep="\t",index=False,compression="gzip")
    focus=d.loc[gate & d.spatial_focus_id.ge(0)].groupby("spatial_focus_id",as_index=False).agg(
        n_observations=(a.cell_id_col,"size"),x_min=(a.x_col,"min"),x_max=(a.x_col,"max"),
        y_min=(a.y_col,"min"),y_max=(a.y_col,"max"),strict_seed_n=("strict_candidate_for_recluster","sum"),
        median_oocyte_hits=("total_oocyte_program_hits","median"),
        median_contradictory_hits=("contradictory_somatic_hits","median"))
    focus["strict_focus_for_recluster"]=focus.spatial_focus_id.isin(seed_foci)
    focus["morphology_status"]="pending_image_or_histology_review"
    focus.to_csv(a.out/"rare_cell_focus_objects.tsv",sep="\t",index=False)
    membership_cols=[a.cell_id_col,a.x_col,a.y_col,"spatial_focus_id","total_oocyte_program_hits","contradictory_somatic_hits","strict_candidate_for_recluster"]
    d.loc[d.strict_focus_for_recluster,membership_cols].to_csv(a.out/"strict_focus_recluster_membership.tsv.gz",sep="\t",index=False,compression="gzip")
    full_membership_cols=membership_cols+["spatial_focus_supported","strict_focus_for_recluster","full_candidate_recluster_member"]
    d.loc[gate,full_membership_cols].to_csv(a.out/"full_candidate_recluster_membership.tsv.gz",sep="\t",index=False,compression="gzip")
    result={"n_screened":len(d),"starting_marker_gate":int(gate.sum()),"full_candidate_recluster_pool":int(gate.sum()),"eps":eps,"min_samples":a.min_samples,"spatial_focus_supported":int(d.spatial_focus_supported.sum()),"strict_candidate_for_recluster":int(d.strict_candidate_for_recluster.sum()),"strict_focus_support_pool":int(d.strict_focus_for_recluster.sum()),"spatial_focus_groups":int(d.loc[d.spatial_focus_id.ge(0),"spatial_focus_id"].nunique()),"strict_focus_groups":len(seed_foci),"canonical_recluster_membership":"full_candidate_recluster_membership.tsv.gz","warning":"The complete starting marker gate is the canonical query-only recluster pool. Strict seeds/foci are supporting evidence and this script does not assign Oocyte."};(a.out/"rare_cell_spatial_focus_summary.json").write_text(json.dumps(result,indent=2)+"\n");print(json.dumps(result,indent=2))
if __name__=="__main__":main()
