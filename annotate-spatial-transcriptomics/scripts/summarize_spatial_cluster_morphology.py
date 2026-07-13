#!/usr/bin/env python3
"""Quantify spatial fragmentation/compactness across candidate resolutions."""
from __future__ import annotations
import argparse,glob,json,re
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
def rtag(p):
    m=re.search(r"res(?:olution)?([0-9]+p?[0-9]*)",Path(p).stem,re.I);return m.group(1).replace("p",".") if m else Path(p).stem
def main():
    p=argparse.ArgumentParser();p.add_argument("--cluster-glob",required=True);p.add_argument("--coordinates",required=True,type=Path);p.add_argument("--out",required=True,type=Path);p.add_argument("--cell-id-col",default="cell_id");p.add_argument("--x-col",default="sdimx");p.add_argument("--y-col",default="sdimy");p.add_argument("--eps",type=float);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True);co=pd.read_csv(a.coordinates,sep=None,engine="python",dtype={a.cell_id_col:str},usecols=[a.cell_id_col,a.x_col,a.y_col]);xy=co[[a.x_col,a.y_col]].to_numpy(float)
    if a.eps is None:
        sample=xy if len(xy)<=50000 else xy[np.linspace(0,len(xy)-1,50000,dtype=int)];dist=NearestNeighbors(n_neighbors=2).fit(sample).kneighbors(sample,return_distance=True)[0][:,1];eps=float(max(np.median(dist)*3,np.quantile(dist,.75)))
    else:eps=a.eps
    rows=[]
    for f in glob.glob(a.cluster_glob):
        c=pd.read_csv(f,sep=None,engine="python",dtype={a.cell_id_col:str,"cluster":str});d=c.merge(co,on=a.cell_id_col,validate="one_to_one");res=rtag(f)
        for cl,x in d.groupby("cluster"):
            pts=x[[a.x_col,a.y_col]].to_numpy(float);lab=DBSCAN(eps=eps,min_samples=2).fit_predict(pts) if len(pts)>1 else np.array([-1]);valid=lab[lab>=0];sizes=pd.Series(valid).value_counts();largest=int(sizes.max()) if len(sizes) else 1;span=np.ptp(pts,axis=0) if len(pts)>1 else np.array([0.,0.]);rows.append({"resolution":res,"cluster":cl,"n_observations":len(x),"spatial_components":int(len(sizes)),"noise_fraction":float((lab<0).mean()),"largest_component_fraction":largest/len(x),"x_span":span[0],"y_span":span[1],"bbox_aspect_ratio":float((max(span)+1)/(min(span)+1))})
    detail=pd.DataFrame(rows);detail.to_csv(a.out/"spatial_cluster_morphology.tsv",sep="\t",index=False);summary=detail.groupby("resolution").apply(lambda x:pd.Series({"n_clusters":len(x),"weighted_largest_component_fraction":np.average(x.largest_component_fraction,weights=x.n_observations),"weighted_noise_fraction":np.average(x.noise_fraction,weights=x.n_observations),"median_components":x.spatial_components.median()})).reset_index();summary.to_csv(a.out/"spatial_resolution_morphology_summary.tsv",sep="\t",index=False);manifest={"status":"EVIDENCE_ONLY","eps":eps,"n_resolutions":len(summary),"warning":"Fragmentation is interpreted against expected lineage morphology; vascular networks and repeated follicles need not be compact."};(a.out/"spatial_morphology_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n");print(json.dumps(manifest,indent=2))
if __name__=="__main__":main()
