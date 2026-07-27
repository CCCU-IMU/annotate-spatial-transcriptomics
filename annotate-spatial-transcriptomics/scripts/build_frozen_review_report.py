#!/usr/bin/env python3
"""Build one portable, sample-agnostic HTML review/final annotation report."""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import os
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import quote


def open_text(path: Path):
    return gzip.open(path, "rt", newline="") if path.suffix == ".gz" else path.open(newline="")


def read_tsv(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.is_file():
        return []
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_json(path: Path | None) -> dict:
    if path is None or not path.is_file():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def rel_url(path: str | Path, out: Path) -> str:
    if not path:
        return ""
    target = Path(path).resolve()
    return quote(os.path.relpath(target, out.parent).replace(os.sep, "/"))


def image(path: str | Path, out: Path, alt: str, wide: bool = False) -> str:
    if not path or not Path(path).is_file():
        return '<p class="muted">图像不可用。</p>'
    link = rel_url(path, out)
    klass = "wide-figure" if wide else "figure"
    return f'<a href="{link}"><img class="{klass}" src="{link}" alt="{html.escape(alt)}"></a>'


def table(rows: list[dict[str, object]], columns: list[tuple[str, str]]) -> str:
    if not rows:
        return '<p class="muted">无记录。</p>'
    head = "".join(f"<th>{html.escape(title)}</th>" for _, title in columns)
    body = []
    for row in rows:
        cells = "".join(f"<td>{html.escape(str(row.get(key, '')))}</td>" for key, _ in columns)
        body.append(f"<tr>{cells}</tr>")
    return f'<div class="table-wrap"><table><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'


def top_deg(rows: list[dict[str, str]], n: int = 5) -> dict[str, list[str]]:
    ranked: dict[str, list[tuple[float, float, str]]] = defaultdict(list)
    for row in rows:
        label = row.get("label", "")
        gene = row.get("gene", "") or row.get("rownames", "")
        if not label or not gene or gene.startswith("LOC") or "-r0" in gene:
            continue
        try:
            pct = float(row.get("pct.1", row.get("pct_expressed_absolute", "1")) or 0)
            lfc = float(row.get("avg_log2FC", row.get("avg_logFC", "0")) or 0)
            p_adj = float(row.get("p_val_adj", "1") or 1)
        except ValueError:
            continue
        if pct >= 0.05 and lfc > 0:
            ranked[label].append((lfc, -p_adj, gene))
    result: dict[str, list[str]] = {}
    for label, values in ranked.items():
        genes: list[str] = []
        for _, _, gene in sorted(values, reverse=True):
            if gene not in genes:
                genes.append(gene)
            if len(genes) == n:
                break
        result[label] = genes
    return result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--biological-context", default="")
    ap.add_argument("--observation-unit", default="cellbin")
    ap.add_argument("--release-status", choices=["pending_user_review", "approved_final"], default="pending_user_review")
    ap.add_argument("--membership", required=True, type=Path)
    ap.add_argument("--release-manifest", required=True, type=Path)
    ap.add_argument("--support", type=Path)
    ap.add_argument("--maps-index", required=True, type=Path)
    ap.add_argument("--maps-dir", required=True, type=Path)
    ap.add_argument("--canonical-dotplots", required=True, type=Path)
    ap.add_argument("--state-dotplots", type=Path)
    ap.add_argument("--canonical-gene-panels", action="append", type=Path, default=[])
    ap.add_argument("--state-gene-panels", type=Path)
    ap.add_argument("--broad-deg", type=Path)
    ap.add_argument("--fine-deg", type=Path)
    ap.add_argument("--fine-audit", type=Path)
    ap.add_argument("--atlas-review", type=Path)
    ap.add_argument("--zero-census", type=Path)
    ap.add_argument("--biological-quality-review", type=Path)
    ap.add_argument("--annotated-rds", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    release = read_json(args.release_manifest)
    membership_rows = read_tsv(args.membership)
    total = int(release.get("n_analysis_set", release.get("n_observations", len(membership_rows))))
    qc_n = int(release.get("residual_qc_n", 0))
    qc_fraction = float(release.get("residual_qc_fraction", release.get("residual_unresolved_or_qc_fraction", qc_n / max(total, 1))))
    broad_census = {str(k): int(v) for k, v in release.get("broad_census", {}).items()}
    if qc_n:
        broad_census["QC/Unknown"] = qc_n
    fine_census = {str(k): int(v) for k, v in release.get("fine_census", {}).items()}
    state_census = {str(k): int(v) for k, v in release.get("state_census", {}).items()}
    membership_record = release.get("membership", {})
    semantic_hash = str(membership_record.get("semantic_sha256", release.get("semantic_hash", "")))
    membership_sha = str(membership_record.get("sha256", release.get("membership_sha256", "")))

    support_rows = read_tsv(args.support)
    support = {(row.get("level", ""), row.get("label", "")): row for row in support_rows}
    nodes = read_tsv(args.maps_index)
    node_by_key = {(row.get("level", ""), row.get("label", "")): row for row in nodes}
    gene_panels = {
        row.get("marker_group", ""): row
        for path in args.canonical_gene_panels for row in read_tsv(path)
    }
    state_panels = {row.get("marker_group", ""): row for row in read_tsv(args.state_gene_panels)}
    broad_deg = top_deg(read_tsv(args.broad_deg))
    fine_deg = top_deg(read_tsv(args.fine_deg))
    source_counts: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    for row in membership_rows:
        broad = row.get("final_broad_label", "") or "QC/Unknown"
        source_counts[("broad", broad)][row.get("broad_freeze_source", "") or row.get("assignment_origin", "") or row.get("final_state", "") or "unresolved"] += 1
        fine = row.get("final_fine_label", "")
        if fine:
            source_counts[("subtype", fine)][row.get("final_fine_assignment_source", "") or "second-round fine evidence"] += 1
        for state in (row.get("state_annotations", "") or row.get("final_state_annotation", "")).split(";"):
            if state:
                source_counts[("state", state)]["second-round state program"] += 1

    dot_rows = read_tsv(args.canonical_dotplots) + read_tsv(args.state_dotplots)

    def dotplot_toggle(level: str, panel: str, title: str) -> str:
        row = next((x for x in dot_rows if x.get("level") == level and x.get("panel") == panel), None)
        if not row:
            return '<p class="muted">该层级无可发布 dotplot。</p>'
        ident = f"dot_{level}_{panel}".replace("-", "_")
        normal = rel_url(row.get("png", ""), args.out)
        absolute = rel_url(row.get("absolute_png", ""), args.out)
        return (
            f'<div class="plot-block"><h3>{html.escape(title)}</h3><div class="toggle">'
            f'<button class="active" onclick="swapPlot(\'{ident}\',\'{normal}\',this)">基因内归一化</button>'
            f'<button onclick="swapPlot(\'{ident}\',\'{absolute}\',this)">绝对表达/检出率</button></div>'
            f'<a id="{ident}_link" href="{normal}"><img id="{ident}" class="wide-figure" src="{normal}" alt="{html.escape(title)}"></a></div>'
        )

    def source_text(level: str, label: str) -> str:
        return "; ".join(f"{name}: {count:,}" for name, count in source_counts[(level, label)].most_common(4))

    def cards(level: str, census: dict[str, int], panels: dict[str, dict[str, str]], degs: dict[str, list[str]]) -> str:
        blocks = []
        for label, count in sorted(census.items(), key=lambda item: (-item[1], item[0])):
            info = support.get((level, label), {})
            node = node_by_key.get((level, label), {})
            panel = panels.get(label, {})
            parent = info.get("parent_label", "") or node.get("parent_label", "")
            title = f"{parent} → {label}" if parent and level != "broad" else label
            markers = info.get("canonical_markers", "") or panel.get("available_genes", "")
            blocks.append(
                '<article class="annotation-card">'
                f'<h3>{html.escape(title)} <span class="n">n={count:,} · {100*count/max(total,1):.2f}%</span></h3>'
                '<div class="card-grid"><div>' + image(node.get("png", ""), args.out, f"{title} spatial") + '</div><div class="evidence">'
                f'<p><b>定义来源</b><br>{html.escape(info.get("evidence_source", "第二轮 query-only 证据与最终生物学复核"))}</p>'
                f'<p><b>成员来源</b><br>{html.escape(source_text(level, label) or "未单独登记")}</p>'
                f'<p><b>典型 marker</b><br>{html.escape(markers or "见 marker dotplot")}</p>'
                f'<p><b>Top DEG</b><br>{html.escape(";".join(degs.get(label, [])) or "未单独计算/未达到阈值")}</p>'
                f'<p><b>空间支撑</b><br>{html.escape(info.get("spatial_support", "见完整成员空间高亮与生物学复核"))}</p>'
                f'<p><b>结论</b><br>{html.escape(info.get("release_interpretation", "保留当前最具体且可靠的身份层级"))}</p></div></div>'
                + (f'<details><summary>所有 {args.observation_unit} 的 marker 空间投影</summary>{image(panel.get("png", ""), args.out, title, True)}</details>' if panel.get("png") else "")
                + '</article>'
            )
        return "".join(blocks)

    quality = read_json(args.biological_quality_review)
    endpoints = quality.get("quality_endpoints", {})
    follicle = endpoints.get("follicle_roi_histology", {})
    status_text = "已同意审阅 · 最终注释" if args.release_status == "approved_final" else "待用户审阅 · 冻结候选"
    broad_overview = args.maps_dir / "figures/final_broad_spatial.png"
    fine_overview = args.maps_dir / "figures/final_subtype_spatial.png"
    state_overview = args.maps_dir / "figures/final_state_spatial.png"
    css = """
    :root{--ink:#18212b;--muted:#66717d;--line:#dfe5eb;--accent:#a51c30;--bg:#f4f6f8;--card:#fff}*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font-family:Inter,"Noto Sans SC","Microsoft YaHei",Arial,sans-serif;line-height:1.55}main{max-width:1500px;margin:auto;padding:28px}.hero{background:linear-gradient(135deg,#17202b,#313d4c);color:white;padding:32px;border-radius:18px}.hero h1{margin:0 0 8px;font-size:30px}.hero p{margin:4px 0;color:#dce4ec}.status{display:inline-block;background:#f4b942;color:#17202b;padding:4px 10px;border-radius:999px;font-weight:700}.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px;margin:18px 0}.metric,.panel,.annotation-card{background:var(--card);border:1px solid var(--line);border-radius:14px;box-shadow:0 3px 12px #15202b0d}.metric{padding:16px}.metric b{display:block;font-size:24px;color:var(--accent)}section{margin:30px 0}.panel{padding:20px;margin:14px 0}.overview-grid,.card-grid{display:grid;grid-template-columns:minmax(0,1.2fr) minmax(320px,.8fr);gap:18px;align-items:start}.annotation-card{padding:18px;margin:16px 0}.annotation-card h3{margin-top:0}.n{font-size:14px;color:var(--muted);font-weight:500}.figure{width:100%;max-height:760px;object-fit:contain;background:white}.wide-figure{display:block;width:100%;max-height:1000px;object-fit:contain;background:white}.evidence p{margin:0 0 12px}.muted{color:var(--muted)}.toggle{display:flex;gap:8px;margin:8px 0}.toggle button{border:1px solid var(--line);background:#fff;padding:7px 12px;border-radius:7px;cursor:pointer}.toggle button.active{background:var(--accent);border-color:var(--accent);color:#fff}.plot-block{margin:18px 0}.table-wrap{overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:8px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}th{background:#eef2f5}.hash{font-family:ui-monospace,monospace;word-break:break-all;font-size:12px}@media(max-width:900px){main{padding:12px}.overview-grid,.card-grid{grid-template-columns:1fr}}
    """
    script = "function swapPlot(id,src,button){const img=document.getElementById(id);img.src=src;document.getElementById(id+'_link').href=src;button.parentElement.querySelectorAll('button').forEach(b=>b.classList.remove('active'));button.classList.add('active')}"
    report = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(args.sample_id)} v2.2.0 注释报告</title><style>{css}</style></head><body><main>
    <header class="hero"><span class="status">{status_text}</span><h1>{html.escape(args.sample_id)} · annotate-spatial-transcriptomics v2.2.0</h1><p>{html.escape(args.biological_context)}</p><p>结果来自当前样本 query-only 第二轮重聚类；历史标签与同批样本注释未参与 membership 冻结。</p><p class="hash">semantic hash: {html.escape(semantic_hash)}</p></header>
    <div class="metrics"><div class="metric"><b>{total:,}</b>分析 {html.escape(args.observation_unit)}</div><div class="metric"><b>{len([x for x in broad_census if x!='QC/Unknown'])}</b>生物学 broad</div><div class="metric"><b>{len(fine_census)}</b>可靠 fine</div><div class="metric"><b>{len(state_census)}</b>独立状态</div><div class="metric"><b>{qc_fraction*100:.2f}%</b>QC/Unknown</div><div class="metric"><b>{broad_census.get('Oocyte',0):,}</b>canonical Oocyte</div></div>
    <section><h2>1. 全组织空间结果</h2><div class="overview-grid"><div class="panel"><h3>Broad + QC/Unknown</h3>{image(broad_overview,args.out,'broad spatial',True)}</div><div class="panel"><h3>可靠 fine</h3>{image(fine_overview,args.out,'fine spatial',True)}</div></div>{('<div class="panel"><h3>状态</h3>'+image(state_overview,args.out,'state spatial',True)+'</div>') if state_overview.is_file() else ''}</section>
    <section><h2>2. 典型 marker dotplot</h2><p class="muted">同一位置切换基因内归一化与绝对表达/检出率，两者来自同一数值表。</p>{dotplot_toggle('broad','canonical','Broad 典型 marker（每类 5 个）')}{dotplot_toggle('subtype','canonical','可靠 fine 典型 marker（每类 5 个）')}{dotplot_toggle('broad','state','状态 marker')}</section>
    <section><h2>3. Broad 注释与逐类支撑</h2>{cards('broad',broad_census,gene_panels,broad_deg)}</section>
    <section><h2>4. 可靠 fine 注释</h2><p class="muted">证据不足时保留 broad 背景，不强求亚型。</p>{cards('subtype',fine_census,gene_panels,fine_deg)}</section>
    <section><h2>5. 独立状态</h2>{cards('state',state_census,state_panels,{})}</section>
    <section><h2>6. Fine 候选完整审计</h2>{table(read_tsv(args.fine_audit),[('parent_broad_label','Broad parent'),('release_label','候选'),('status','结论'),('supported_subcluster_n','支持亚簇'),('rationale','理由')])}</section>
    <section><h2>7. Atlas、缺失谱系与卵泡组织学复核</h2><div class="panel"><p><b>生物学复核：</b>{html.escape(str(quality.get('status','未提供')))}</p><p><b>Oocyte：</b>{html.escape(str(endpoints.get('oocyte_annotation_quality',{}).get('status','未提供')))}</p><p><b>卵泡 ROI：</b>{html.escape(str(follicle.get('status','未提供')))}；antral ROI {html.escape(str(follicle.get('antral_roi_n','NA')))}；腔体 {html.escape(str(follicle.get('antral_cavity_status','NA')))}</p><p><b>层次：</b>{html.escape(str(follicle.get('histological_sequence','未评估')))}</p></div>{table(read_tsv(args.zero_census),[('release_broad_label','候选谱系'),('decision','复核结论'),('rationale','解释')])}</section>
    <section><h2>8. 交付物</h2><div class="panel"><p><b>冻结 membership：</b><span class="hash">{html.escape(str(args.membership.resolve()))}</span></p><p><b>写回注释 RDS：</b><a href="{rel_url(args.annotated_rds,args.out)}">{html.escape(args.annotated_rds.name)}</a></p><p><b>membership SHA256：</b><span class="hash">{html.escape(membership_sha)}</span></p></div></section>
    </main><script>{script}</script></body></html>'''
    args.out.write_text(report, encoding="utf-8")
    print(json.dumps({"status": "PASS", "report": str(args.out.resolve()), "sample_id": args.sample_id, "n_observations": total}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
