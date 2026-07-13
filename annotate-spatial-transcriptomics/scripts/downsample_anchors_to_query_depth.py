#!/usr/bin/env python3
"""Create held-out anchors downsampled to a low-depth query distribution."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import anndata as ad
import numpy as np,pandas as pd
from scipy import sparse

def main():
    p=argparse.ArgumentParser();p.add_argument("--h5ad",required=True,type=Path);p.add_argument("--anchor-membership",required=True,type=Path);p.add_argument("--target-membership",required=True,type=Path);p.add_argument("--out",required=True,type=Path);p.add_argument("--cell-id-col",default="cell_id");p.add_argument("--truth-col",default="true_label");p.add_argument("--layer");p.add_argument("--seed",type=int,default=20260713);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True);x=ad.read_h5ad(a.h5ad);x.obs_names=x.obs_names.astype(str);am=pd.read_csv(a.anchor_membership,sep="\t",dtype={a.cell_id_col:str});tm=pd.read_csv(a.target_membership,sep="\t",dtype={a.cell_id_col:str});
    for name,m in [("anchor",am),("target",tm)]:
        miss=pd.Index(m[a.cell_id_col]).difference(x.obs_names)
        if len(miss):raise SystemExit(f"{len(miss)} {name} IDs absent")
    mat=x.layers[a.layer] if a.layer else x.X;mat=mat.tocsr() if sparse.issparse(mat) else sparse.csr_matrix(mat);pos={v:i for i,v in enumerate(x.obs_names)};ai=np.array([pos[v] for v in am[a.cell_id_col]]);ti=np.array([pos[v] for v in tm[a.cell_id_col]]);target_tot=np.asarray(mat[ti].sum(1)).ravel().astype(int);target_tot=target_tot[target_tot>0]
    if not len(target_tot):raise SystemExit("target membership has no positive-count observations")
    rng=np.random.default_rng(a.seed);rows=[];cols=[];vals=[];orig=[];used=[]
    for j,i in enumerate(ai):
        row=mat.getrow(i);counts=np.rint(row.data).astype(int);total=int(counts.sum());depth=min(total,int(rng.choice(target_tot)));orig.append(total);used.append(depth)
        if depth>0 and total>0:
            draw=rng.multinomial(depth,counts/counts.sum());nz=draw>0;rows.extend([j]*int(nz.sum()));cols.extend(row.indices[nz]);vals.extend(draw[nz])
    outmat=sparse.csr_matrix((vals,(rows,cols)),shape=(len(ai),x.n_vars));obs=am.set_index(a.cell_id_col).copy();obs["original_counts"]=orig;obs["downsampled_counts"]=used;out=ad.AnnData(outmat,obs=obs,var=x.var.copy());out.layers["counts"]=outmat.copy();out.write_h5ad(a.out/"depth_matched_anchors.h5ad",compression="gzip");summary={"n_anchors":len(ai),"n_target":len(ti),"target_depth_median":float(np.median(target_tot)),"anchor_depth_median_before":float(np.median(orig)),"anchor_depth_median_after":float(np.median(used)),"truth_col":a.truth_col,"warning":"Use only as held-out calibration evidence, never as a reference-label source for the same observations."};(a.out/"depth_matching_summary.json").write_text(json.dumps(summary,indent=2)+"\n");print(json.dumps(summary,indent=2))
if __name__=="__main__":main()
