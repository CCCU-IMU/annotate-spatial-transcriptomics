#!/usr/bin/env python3
"""Calibrate strict rare-cell seeds from the query, then expand spatial objects."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import pandas as pd
def main():
    p=argparse.ArgumentParser();p.add_argument("--screen",required=True,type=Path);p.add_argument("--out",required=True,type=Path);p.add_argument("--cell-id-col",default="cell_id");p.add_argument("--positive-quantile",type=float,default=.75);p.add_argument("--contradictory-quantile",type=float,default=.25);p.add_argument("--modules-required",type=int,default=3);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True);d=pd.read_csv(a.screen,sep="\t",dtype={a.cell_id_col:str});g=d[d.starting_marker_gate.astype(str).str.lower().isin(["true","1"])].copy()
    need={"spatial_focus_supported","spatial_focus_id","total_oocyte_program_hits","identity_core_hits","modules_detected","contradictory_somatic_hits"};missing=need-set(g.columns)
    if missing:raise SystemExit(f"screen lacks columns: {sorted(missing)}")
    if len(g)<10:raise SystemExit("too few starting-gate observations for query-calibrated quantiles")
    tt=float(g.total_oocyte_program_hits.quantile(a.positive_quantile));ti=float(g.identity_core_hits.quantile(a.positive_quantile));tc=float(g.contradictory_somatic_hits.quantile(a.contradictory_quantile));supported=g.spatial_focus_supported.astype(str).str.lower().isin(["true","1"]);seed=supported & g.total_oocyte_program_hits.ge(tt)&g.identity_core_hits.ge(ti)&g.modules_detected.ge(a.modules_required)&g.contradictory_somatic_hits.le(tc);groups=set(g.loc[seed,"spatial_focus_id"]);pool=g[g.spatial_focus_id.isin(groups)].copy();pool["strict_seed_calibrated"]=pool[a.cell_id_col].isin(set(g.loc[seed,a.cell_id_col]));pool.to_csv(a.out/"calibrated_rare_focus_pool.tsv.gz",sep="\t",index=False,compression="gzip");summary={"status":"CALIBRATED_CANDIDATE_POOL","starting_gate":len(g),"spatial_supported":int(supported.sum()),"positive_quantile":a.positive_quantile,"contradictory_quantile":a.contradictory_quantile,"thresholds":{"total_program_hits":tt,"identity_core_hits":ti,"contradictory_hits":tc,"modules_required":a.modules_required},"strict_seeds":int(seed.sum()),"strict_seed_groups":len(groups),"expanded_focus_pool":len(pool),"warning":"Query-calibrated seeds and expanded spatial objects remain candidates, not final rare-cell labels; observations are not biological-cell counts."};(a.out/"calibrated_rare_focus_summary.json").write_text(json.dumps(summary,indent=2)+"\n");print(json.dumps(summary,indent=2))
if __name__=="__main__":main()
