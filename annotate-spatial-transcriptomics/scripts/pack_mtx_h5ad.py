#!/usr/bin/env python3
"""Pack the framework R sparse-export contract into a standard H5AD."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import anndata as ad
import pandas as pd
from scipy.io import mmread
from scipy import sparse
def main():
    p=argparse.ArgumentParser();p.add_argument("export_dir",type=Path);p.add_argument("--out",required=True,type=Path);p.add_argument("--cell-id-col",default="cell_id");a=p.parse_args();e=a.export_dir;obs=pd.read_csv(e/"observations.tsv",sep="\t",dtype={a.cell_id_col:str});var=pd.read_csv(e/"features.tsv",sep="\t",dtype=str);mat=mmread(e/"counts_genes_by_observations.mtx").tocsr().T
    if mat.shape!=(len(obs),len(var)):raise SystemExit("matrix/metadata dimension mismatch")
    if obs[a.cell_id_col].duplicated().any() or var.gene.duplicated().any():raise SystemExit("duplicate observation or feature IDs")
    obs=obs.set_index(a.cell_id_col);var=var.set_index("gene");x=ad.AnnData(sparse.csr_matrix(mat),obs=obs,var=var);x.layers["counts"]=x.X.copy();a.out.parent.mkdir(parents=True,exist_ok=True);x.write_h5ad(a.out,compression="gzip");manifest={"status":"STANDARD_H5AD_ADAPTER","n_observations":x.n_obs,"n_genes":x.n_vars,"source_export":str(e.resolve()),"output":str(a.out.resolve())};a.out.with_suffix(".manifest.json").write_text(json.dumps(manifest,indent=2)+"\n");print(json.dumps(manifest,indent=2))
if __name__=="__main__":main()
