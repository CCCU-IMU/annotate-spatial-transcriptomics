#!/usr/bin/env python3
"""Render auditable whole-section rare-cell screening and focus labels."""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

def main():
    p=argparse.ArgumentParser();p.add_argument("--screen",required=True,type=Path);p.add_argument("--background",required=True,type=Path);p.add_argument("--candidate-membership",type=Path);p.add_argument("--out",required=True,type=Path);p.add_argument("--cell-id-col",default="cell_id");p.add_argument("--x-col",default="sdimx");p.add_argument("--y-col",default="sdimy");p.add_argument("--title",default="Rare-cell strict validation");a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    s=pd.read_csv(a.screen,sep="\t",dtype={a.cell_id_col:str});b=pd.read_csv(a.background,sep=None,engine="python",dtype={a.cell_id_col:str},usecols=[a.cell_id_col,a.x_col,a.y_col])
    if a.x_col not in s or a.y_col not in s:s=s.merge(b,on=a.cell_id_col,how="left",validate="one_to_one")
    fig,ax=plt.subplots(figsize=(11,9));ax.scatter(b[a.x_col],b[a.y_col],s=.08,c="#d9d9d9",rasterized=True,label="all observations")
    g=s[s.starting_marker_gate.astype(str).str.lower().isin(["true","1"])];ax.scatter(g[a.x_col],g[a.y_col],s=4,c="#f4a261",alpha=.75,label="multi-marker screen")
    if a.candidate_membership:
        q=pd.read_csv(a.candidate_membership,sep="\t",dtype={a.cell_id_col:str});
        if a.x_col not in q or a.y_col not in q:q=q.merge(b,on=a.cell_id_col,how="left",validate="one_to_one")
    else:q=s[s.strict_focus_for_recluster.astype(str).str.lower().isin(["true","1"])]
    ax.scatter(q[a.x_col],q[a.y_col],s=9,c="#d62828",alpha=.9,label="calibrated focus recluster pool")
    if "strict_seed_calibrated" in q:
        seed=q[q.strict_seed_calibrated.astype(str).str.lower().isin(["true","1"])];ax.scatter(seed[a.x_col],seed[a.y_col],s=13,c="#370617",alpha=.95,label="calibrated strict seeds")
    for fid,x in q.groupby("spatial_focus_id"):
        ax.text(x[a.x_col].mean(),x[a.y_col].mean(),str(int(fid)),fontsize=6,color="#6a040f")
    ax.invert_yaxis();ax.set_aspect("equal");ax.set_axis_off();ax.set_title(a.title);ax.legend(loc="best",markerscale=3,frameon=False)
    fig.tight_layout();fig.savefig(a.out/"rare_cell_focus_whole_section.png",dpi=400,bbox_inches="tight");fig.savefig(a.out/"rare_cell_focus_whole_section.pdf",bbox_inches="tight");plt.close(fig)
    print(a.out/"rare_cell_focus_whole_section.png")
if __name__=="__main__":main()
