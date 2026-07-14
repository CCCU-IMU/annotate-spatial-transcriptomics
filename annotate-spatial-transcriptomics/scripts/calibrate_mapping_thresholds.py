#!/usr/bin/env python3
"""Calibrate mapping confidence/margin on held-out anchors and retain rejects."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
def main():
    p=argparse.ArgumentParser(description="Legacy one-tier diagnostic calibration; never eligible for writeback");p.add_argument("--predictions",required=True,type=Path);p.add_argument("--truth",required=True,type=Path);p.add_argument("--query-predictions",required=True,type=Path);p.add_argument("--out",required=True,type=Path);p.add_argument("--cell-id-col",default="cell_id");p.add_argument("--truth-col",default="true_label");p.add_argument("--target-precision",type=float,default=.9);p.add_argument("--min-support",type=int,default=20);p.add_argument("--diagnostic-only",action="store_true");a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    if not a.diagnostic_only:
        raise SystemExit("Legacy one-tier calibration is disabled for writeback. Use calibrate_tiered_mapping_thresholds.py with a query-like held-out origin manifest.")
    pr=pd.read_csv(a.predictions,sep="\t",dtype={a.cell_id_col:str});tr=pd.read_csv(a.truth,sep="\t",dtype={a.cell_id_col:str});d=pr.merge(tr,on=a.cell_id_col,validate="one_to_one");rows=[]
    for lab,x in d.groupby("predicted_label"):
        best=None
        cgrid=np.unique(np.quantile(x.confidence,np.linspace(0,1,min(31,len(x)))))
        mgrid=np.unique(np.quantile(x.margin,np.linspace(0,1,min(31,len(x)))))
        for c in cgrid:
            for m in mgrid:
                y=x[(x.confidence>=c)&(x.margin>=m)];n=len(y)
                if n<a.min_support:continue
                precision=float((y.predicted_label==y[a.truth_col]).mean());coverage=n/len(x)
                if precision>=a.target_precision and (best is None or coverage>best["calibration_coverage"]):best={"predicted_label":lab,"confidence_threshold":float(c),"margin_threshold":float(m),"calibration_precision":precision,"calibration_coverage":coverage,"calibration_support":n}
        if best:rows.append(best)
    cols=["predicted_label","confidence_threshold","margin_threshold","calibration_precision","calibration_coverage","calibration_support"];thresholds=pd.DataFrame(rows,columns=cols);thresholds.to_csv(a.out/"mapping_thresholds.tsv",sep="\t",index=False);q=pd.read_csv(a.query_predictions,sep="\t",dtype={a.cell_id_col:str});q=q.merge(thresholds,on="predicted_label",how="left");q["legacy_threshold_pass"]=q.confidence.ge(q.confidence_threshold)&q.margin.ge(q.margin_threshold);q["mapping_status"]="legacy_diagnostic_only_not_writeback_eligible";q["writeback_eligible"]=False;q["fine_anchor_eligible"]=False;q.to_csv(a.out/"calibrated_query_mapping.tsv.gz",sep="\t",index=False,compression="gzip");result={"status":"LEGACY_DIAGNOSTIC_ONLY","target_precision":a.target_precision,"labels_with_thresholds":len(thresholds),"diagnostic_threshold_pass":int(q.legacy_threshold_pass.sum()),"writeback_eligible":0,"warning":"This legacy combined tier is diagnostic only. Use tiered calibration on query-like held-out current-query anchors."};(a.out/"calibration_summary.json").write_text(json.dumps(result,indent=2)+"\n");print(json.dumps(result,indent=2))
if __name__=="__main__":main()
