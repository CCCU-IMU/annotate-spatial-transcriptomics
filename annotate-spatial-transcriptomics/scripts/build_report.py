#!/usr/bin/env python3
"""Build a portable static HTML report from framework assets and state."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
from pathlib import Path
from dependency_manifest import build as build_dependency_manifest
from evidence_schema_lib import active_registry_rows


def read_tsv(path: Path):
    if not path.exists(): return []
    with path.open(newline="", encoding="utf-8") as h: return list(csv.DictReader(h, delimiter="\t"))


def rel(path: str, report_dir: Path) -> str:
    try: return str(Path(path).resolve().relative_to(report_dir.parent.resolve()))
    except ValueError: return Path(path).resolve().as_uri()


def asset_path(path: str, root: Path) -> Path:
    p = Path(path)
    return p if p.is_absolute() else root / p


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_valid_confirmation(root: Path) -> dict:
    path = root / "state/final_annotation_confirmation.json"
    if not path.is_file():
        raise SystemExit(
            "final report blocked: request and obtain explicit user confirmation first"
        )
    record = json.loads(path.read_text(encoding="utf-8"))
    if record.get("status") != "CONFIRMED":
        raise SystemExit("final report blocked: user confirmation status is not CONFIRMED")
    for key, hash_key in (
        ("cell_ledger", "cell_ledger_sha256"),
        ("cluster_ledger", "cluster_ledger_sha256"),
        ("completion_gate", "completion_gate_sha256"),
        ("master_quality_approval", "master_quality_approval_sha256"),
        ("confirmation_review_report", "confirmation_review_report_sha256"),
        ("confirmation_review_manifest", "confirmation_review_manifest_sha256"),
        ("annotation_support_registry", "annotation_support_registry_sha256"),
    ):
        target = root / str(record.get(key, ""))
        if not target.is_file() or sha256(target) != record.get(hash_key):
            raise SystemExit(
                f"final report blocked: confirmed annotation snapshot is stale ({key})"
            )
    return record


def active_decisions(rows):
    def rid(x): return x.get("decision_id") or f"{x.get('decision_version','')}:{x.get('source_run_id','')}:{x.get('source_cluster','')}"
    superseded=set()
    for x in rows: superseded.update(v for v in re.split(r"[;,\s]+",x.get("supersedes","").strip()) if v)
    return [x for x in rows if rid(x) not in superseded]


def node_id(level: str, label: str, parent: str = "") -> str:
    raw = f"{level}-{parent}-{label}".lower()
    clean = re.sub(r"[^a-z0-9\u4e00-\u9fff]+", "-", raw).strip("-")
    return "node-" + (clean or "unlabeled")


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("project_root", type=Path); ap.add_argument("--language", default="zh"); args = ap.parse_args()
    root = args.project_root.resolve(); cfg = json.loads((root / "config/project.json").read_text())
    confirmation = load_valid_confirmation(root)
    report_dir = root / "report"; report_dir.mkdir(parents=True, exist_ok=True)
    dot = read_tsv(root / "figures/marker_dotplots/marker_dotplot_asset_index.tsv")
    state = (root / "state/annotation_state.md").read_text(encoding="utf-8") if (root / "state/annotation_state.md").exists() else ""
    clustering = read_tsv(root / "state/clustering_decision_ledger.tsv")
    clusters = read_tsv(root / "state/cluster_decision_ledger.tsv")
    active_clusters = active_decisions(clusters)
    cohorts = active_registry_rows(read_tsv(root / "state/recluster_cohort_registry.tsv"), "cohort_id")
    direct_returns = active_registry_rows(read_tsv(root / "state/direct_return_registry.tsv"), "return_id")
    runs = read_tsv(root / "state/run_registry.tsv")
    queue = read_tsv(root / "state/next_action_queue.tsv")
    route_attempts = active_registry_rows(read_tsv(root / "state/route_attempt_registry.tsv"), "route_attempt_id")
    workflow_events = read_tsv(root / "state/workflow_event_registry.tsv")
    annotation_views = [row for row in read_tsv(root / "state/annotation_view_registry.tsv") if row.get("view") == "final"]
    workflow_audit_path = root / "provenance/direct_lineage_workflow_audit.json"
    workflow_audit = json.loads(workflow_audit_path.read_text()) if workflow_audit_path.exists() else {"status":"NOT_RUN","gaps":[]}
    context_path = root / "config/biological_context.json"
    context = json.loads(context_path.read_text(encoding="utf-8")) if context_path.exists() else {}
    gate_path = root / "provenance/completion_gate.json"
    gate = json.loads(gate_path.read_text(encoding="utf-8")) if gate_path.exists() else {"status":"NOT_RUN","errors":["completion gate has not run"]}
    master_quality_path = root / "state/master_quality_approval.json"
    master_quality = json.loads(master_quality_path.read_text(encoding="utf-8")) if master_quality_path.exists() else {"status":"MISSING"}
    feature_path = root / "provenance/full_feature_validation.json"
    feature_audit = json.loads(feature_path.read_text(encoding="utf-8")) if feature_path.exists() else {"status":"NOT_RUN"}
    scope_path = root / "provenance/analysis_scope_policy.json"
    scope_policy = json.loads(scope_path.read_text(encoding="utf-8")) if scope_path.exists() else {}
    atlas_policy_path = root / "provenance/atlas_broad_rescue_policy.json"
    if not atlas_policy_path.exists():
        atlas_policy_path = root / "provenance/atlas_broad_rescue_policy_quality_reopen.json"  # legacy project compatibility
    atlas_policy = json.loads(atlas_policy_path.read_text(encoding="utf-8")) if atlas_policy_path.exists() else {}
    open_world_validation_path = root / "provenance/open_world_lineage_audit_validation.json"
    open_world_validation = json.loads(open_world_validation_path.read_text(encoding="utf-8")) if open_world_validation_path.exists() else {"status": "NOT_RUN"}
    open_world_source = asset_path(str(open_world_validation.get("audit", "")), root) if open_world_validation.get("audit") else None
    open_world_catalog = asset_path(str(open_world_validation.get("candidate_catalog", "")), root) if open_world_validation.get("candidate_catalog") else None
    open_world = json.loads(open_world_source.read_text(encoding="utf-8")) if open_world_source and open_world_source.is_file() else {}
    final_census = read_tsv(root / "tables/final_annotation_census.tsv")
    final_status = "最终版（completion gate 已通过）" if gate.get("status") == "PASS" else "初步版（仍需迭代，禁止作为最终注释）"
    nodes = read_tsv(root / "tables/spatial_node_asset_index.tsv")
    genes = read_tsv(root / "tables/spatial_gene_asset_index.tsv")
    def table(rows, cols):
        if not rows: return "<p>暂无记录</p>"
        return "<table><thead><tr>" + "".join(f"<th>{html.escape(c)}</th>" for c in cols) + "</tr></thead><tbody>" + "".join(
            "<tr>" + "".join(f"<td>{html.escape(str(r.get(c,'')))}</td>" for c in cols) + "</tr>" for r in rows
        ) + "</tbody></table>"
    overview_cards = []
    for filename, title in [("final_broad_spatial.png","大类空间图"),("final_subtype_spatial.png","亚群空间图"),("final_broad_UMAP.png","大类 UMAP"),("final_subtype_UMAP.png","亚群 UMAP")]:
        p = root / "figures" / filename
        if p.exists(): overview_cards.append(f"<article class='card'><h3>{title}</h3><img src='../figures/{filename}'></article>")
    open_world_rows = open_world.get("candidate_reviews", []) if isinstance(open_world.get("candidate_reviews", []), list) else []
    if open_world_validation.get("status") != "NOT_RUN":
        overview_cards.append(
            "<article class='card'><h3>开放式谱系发现</h3>"
            f"<p>验证状态：{html.escape(str(open_world_validation.get('status','')))}；已审查家族：{html.escape('、'.join(map(str,open_world_validation.get('reviewed_families',[]))))}。候选目录不是封闭标签表，也不要求每类存在。</p>"
            + table(open_world_rows, ["candidate_id", "outcome", "evidence_summary", "action", "limitation"])
            + "</article>"
        )
    node_cards = []
    for r in nodes:
        p = asset_path(r.get("png", ""), root)
        parent = r.get("parent_label", "")
        path_label = f"{parent} → {r.get('label','')}" if r.get("level") == "subtype" and parent else r.get("label", "")
        nid = node_id(r.get("level", ""), r.get("label", ""), parent)
        if p.exists(): node_cards.append(f"<details id='{html.escape(nid)}' class='node-card' data-search='{html.escape(path_label.lower())}'><summary>{html.escape(r.get('level',''))} → {html.escape(path_label)} (n={html.escape(r.get('n_observations',''))})</summary><p><a class='button' href='#tree'>返回注释树</a></p><img src='../{html.escape(rel(str(p), report_dir))}'></details>")
    fine_parent = {}
    for r in active_clusters:
        if r.get("broad_label") and r.get("fine_label"): fine_parent.setdefault(r["fine_label"], set()).add(r["broad_label"])
    broad_nodes = {r.get("label"): r for r in nodes if r.get("level") == "broad"}
    subtype_nodes = [r for r in nodes if r.get("level") == "subtype"]
    tree_parts = []
    for broad, br in sorted(broad_nodes.items()):
        bp = asset_path(br.get("png", ""), root); children = []; broad_id = node_id("broad", broad or "未定义")
        for sr in subtype_nodes:
            parents = set(filter(None, sr.get("parent_label", "").split(" | "))) or fine_parent.get(sr.get("label"), set())
            if broad not in parents: continue
            sp = asset_path(sr.get("png", ""), root)
            subtype_id = node_id("subtype", sr.get("label", ""), broad or "")
            children.append(f"<details class='tree-node' data-search='{html.escape((broad+' '+sr.get('label','')).lower())}'><summary>{html.escape(sr.get('label',''))} (n={html.escape(sr.get('n_observations',''))})</summary><p><a class='button' href='#{html.escape(subtype_id)}'>跳转到该亚群空间高亮</a></p>" + (f"<img src='../{html.escape(rel(str(sp),report_dir))}'>" if sp.exists() else "") + "</details>")
        tree_parts.append(f"<details class='tree-node' data-search='{html.escape((broad or '未定义').lower())}'><summary><b>{html.escape(broad or '未定义')}</b> (n={html.escape(br.get('n_observations',''))})</summary><p><a class='button' href='#{html.escape(broad_id)}'>跳转到该大类空间高亮</a></p>" + (f"<img src='../{html.escape(rel(str(bp),report_dir))}'>" if bp.exists() else "") + "".join(children) + "</details>")
    gene_groups = {}
    for r in genes:
        p = asset_path(r.get("png", ""), root)
        if p.exists():
            group=r.get('marker_group',r.get('marker_source','未分组'));gene_groups.setdefault(group,[]).append(f"<article class='card'><h3>{html.escape(r.get('gene',''))}</h3><img src='../{html.escape(rel(str(p), report_dir))}'></article>")
    gene_sections="".join(f"<details open><summary><b>{html.escape(group)}</b>（{len(cards)} genes）</summary><div class='grid'>{''.join(cards)}</div></details>" for group,cards in sorted(gene_groups.items()))
    rare_cards = []
    rare_paths=sorted(root.glob("runs/**/rare_cell_focus_whole_section.png"));fullgene_paths=[p for p in rare_paths if "fullgene" in str(p).lower()]
    for p in (fullgene_paths or rare_paths):
        rare_cards.append(f"<article class='card'><h3>{html.escape(p.parent.name)}</h3><img src='../{html.escape(rel(str(p),report_dir))}'></article>")
    state_summary = {}
    unresolved = []
    for r in active_clusters:
        try: n = int(float(r.get("n_observations", 0) or 0))
        except ValueError: n = 0
        state_summary[r.get("state", "unknown")] = state_summary.get(r.get("state", "unknown"), 0) + n
        if r.get("state") not in {"defined_fine", "defined_broad_only", "closed_and_frozen"}: unresolved.append(r)
    summary_rows = [{"state": k, "n_observations": v} for k, v in sorted(state_summary.items(), key=lambda x: -x[1])]
    downloads = []
    download_candidates = [
        ("Final broad DEG", root/"tables/final_broad_DEG_one_vs_rest_all.tsv"), ("Final subtype DEG", root/"tables/final_subtype_DEG_one_vs_rest_all.tsv"),
        ("Broad DEG", root/"tables/broad_DEG_one_vs_rest_all.tsv"), ("Subtype DEG", root/"tables/subtype_DEG_one_vs_rest_all.tsv"),
        ("Broad DEG top", root/"tables/broad_DEG_top100.tsv"), ("Subtype DEG top", root/"tables/subtype_DEG_top100.tsv"),
        ("Cell ledger", root/"state/cell_ledger.tsv.gz"), ("Cell ledger", root/"tables/cell_ledger.tsv.gz"),
        ("Cluster decisions", root/"state/cluster_decision_ledger.tsv"), ("Clustering decisions", root/"state/clustering_decision_ledger.tsv"),
        ("Run registry", root/"state/run_registry.tsv"), ("Recluster cohort registry", root/"state/recluster_cohort_registry.tsv"),
        ("Direct return registry", root/"state/direct_return_registry.tsv"),
        ("State validation", root/"provenance/state_validation.json"), ("Release audit", root/"provenance/release_audit.json"),
        ("Biological context", context_path), ("Next-action queue", root/"state/next_action_queue.tsv"),
        ("Full-feature validation", feature_path),
        ("Analysis-scope policy", scope_path),
        ("Atlas broad-rescue policy", atlas_policy_path),
        ("Open-world lineage audit validation", open_world_validation_path),
        ("Final annotation user confirmation", root/"state/final_annotation_confirmation.json"),
        ("Main-Agent annotation-quality approval", master_quality_path),
        ("Main-Agent quality review request", root/"provenance/master_quality_review_request.json"),
        ("Pre-confirmation lightweight review", root/"review/confirmation/index.html"),
        ("Annotation support registry", root/"state/annotation_support_registry.tsv"),
        ("Confirmation review manifest", root/"provenance/confirmation_review_manifest.json"),
        ("Completion gate", gate_path), ("Iteration plan", root/"provenance/iteration_plan.json"),
        ("Direct-lineage workflow audit", workflow_audit_path), ("Assisted route attempts", root/"state/route_attempt_registry.tsv"),
        ("Workflow events", root/"state/workflow_event_registry.tsv"),
        ("Annotation views", root/"state/annotation_view_registry.tsv"),
        ("Checksums", root/"provenance/checksums.sha256"), ("Session info", root/"provenance/release_sessionInfo.txt"),
    ]
    if open_world_source and open_world_source.is_file():
        download_candidates.append(("Open-world lineage audit source", open_world_source))
    if open_world_catalog and open_world_catalog.is_file():
        download_candidates.append(("Sheep-ovary candidate lineage catalog", open_world_catalog))
    seen_downloads = set()
    for label, p in download_candidates:
        if p.exists() and p not in seen_downloads:
            seen_downloads.add(p); downloads.append(f"<li><a href='../{html.escape(rel(str(p), report_dir))}'>{html.escape(label)}</a></li>")
    sections = []
    for level, title in [("broad", "大类 marker 带树点图"), ("subtype", "亚群 marker 带树点图")]:
        cards = []
        for r in dot:
            if r.get("level") != level: continue
            png = rel(str(asset_path(r["png"], root)), report_dir); pdf = rel(str(asset_path(r["pdf"], root)), report_dir); src = rel(str(asset_path(r["source"], root)), report_dir)
            cards.append(f"<article class='card'><h3>{html.escape(r.get('panel',''))}</h3><p>view={html.escape(r.get('analysis_view','未记录'))}；cohort={html.escape(r.get('evidence_cohort','未记录'))}</p><img src='../{html.escape(png)}'><p><a href='../{html.escape(pdf)}'>PDF</a> · <a href='../{html.escape(src)}'>源表</a></p></article>")
        sections.append(f"<section id='dot-{level}'><h2>{title}</h2><p>marker 基因位于 X 轴并按其支持的细胞类型/程序分栏；左侧树聚类当前标签。点大小和颜色均按基因内部归一化，绝对检出率与绝对平均表达保存在源表。</p><div class='grid'>{''.join(cards) or '<p>缺失：发布审计将失败。</p>'}</div></section>")
    tree_overview_path = root / "figures/final_broad_spatial.png"
    tree_overview = (
        "<article class='card'><h3>注释后大类空间总图</h3>"
        "<img src='../figures/final_broad_spatial.png'></article>"
        if tree_overview_path.exists() else
        "<p>注释后空间总图尚未生成；发布审计将失败。</p>"
    )
    atlas_external = atlas_policy.get("external_atlas_result", atlas_policy.get("current_external_atlas_result_g067", {}))
    atlas_internal = atlas_policy.get("internal_anchor_result", atlas_policy.get("current_internal_anchor_result_g067", {}))
    atlas_rows = []
    if atlas_external:
        atlas_rows.append({
            "channel": "GSE233801 external Atlas", "n_query": atlas_external.get("n_query", ""),
            "high_n": atlas_external.get("high_n", ""), "moderate_only_n": atlas_external.get("moderate_only_n", ""),
            "moderate_or_higher_n": atlas_external.get("moderate_or_higher_n", ""),
            "low_reject_n": atlas_external.get("low_reject_n", ""),
            "stromal_moderate_or_higher_n": atlas_external.get("moderate_or_higher_stromal_n", ""),
            "direct_writeback": atlas_external.get("direct_label_writeback_performed", ""),
        })
    if atlas_internal:
        atlas_rows.append({
            "channel": "internal anchors", "n_query": "route-constrained",
            "high_n": atlas_internal.get("high_n", ""), "moderate_only_n": atlas_internal.get("moderate_only_n", ""),
            "moderate_or_higher_n": atlas_internal.get("moderate_or_higher_n", ""),
            "low_reject_n": "", "stromal_moderate_or_higher_n": "", "direct_writeback": False,
        })
    writeback_counts = atlas_policy.get("writeback_counts", atlas_policy.get("nested_recalibration_writeback_counts", {}))
    atlas_writeback_rows = [writeback_counts] if writeback_counts else []
    body = f"""<!doctype html><html lang='zh'><head><meta charset='utf-8'><title>{html.escape(cfg['sample_id'])} 注释报告</title>
<style>body{{font-family:Arial,'Noto Sans CJK SC',sans-serif;margin:0;background:#f5f7fa;color:#20242a}}header,main{{max-width:1500px;margin:auto;padding:24px}}header{{background:#17324d;color:white;max-width:none}}nav a{{color:#d8efff;margin-right:18px}}section{{background:white;margin:18px 0;padding:20px;border-radius:10px;scroll-margin-top:12px}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(520px,1fr));gap:16px}}.card{{border:1px solid #d9e0e8;padding:12px;border-radius:8px}}img{{width:100%;height:auto}}table{{border-collapse:collapse;width:100%;font-size:12px}}th,td{{border:1px solid #ccd4dd;padding:5px;text-align:left}}pre{{white-space:pre-wrap}}details{{margin:8px 0;padding:6px;border-left:3px solid #d9e0e8}}summary{{cursor:pointer}}.button,button{{display:inline-block;background:#176b87;color:white!important;border:0;border-radius:5px;padding:7px 10px;text-decoration:none;cursor:pointer;margin:3px}}input[type=search]{{padding:8px;min-width:320px;border:1px solid #9ca9b5;border-radius:5px}}</style></head>
<body><header><h1>{html.escape(cfg['sample_id'])} 可追溯注释报告</h1><p>{html.escape(final_status)} · {html.escape(cfg['modality'])} · observation={html.escape(cfg['observation_unit'])} · framework={html.escape(cfg['framework_version'])}</p><nav><a href='#context'>生物学背景</a><a href='#overview'>注释空间图</a><a href='#final'>最终注释</a><a href='#selection'>聚类选择</a><a href='#tree'>注释树</a><a href='#routes'>cohort 与回归裁决</a><a href='#iterations'>迭代状态</a><a href='#uncertainty'>保留状态</a><a href='#dot-broad'>大类点图</a><a href='#dot-subtype'>亚群点图</a><a href='#rare-audit'>Oocyte/特异候选审计</a><a href='#nodes'>空间节点</a><a href='#genes'>空间基因</a><a href='#downloads'>下载与审计</a><a href='#workflow'>全流程详细版</a></nav></header><main>
<section id='context'><h2>生物学背景与完成门</h2><p><b>报告状态：</b>{html.escape(final_status)}</p><p><b>主 Agent 注释质量审批：</b>{html.escape(str(master_quality.get('status')))}；{html.escape(str(master_quality.get('rationale','')))}</p><p><b>用户终审：</b>{html.escape(str(confirmation.get('status')))}；确认时间={html.escape(str(confirmation.get('confirmed_at','')))}；decision version={html.escape(str(confirmation.get('decision_version','')))}。</p><p><b>物种/组织：</b>{html.escape(str(context.get('species','未记录')))} / {html.escape(str(context.get('tissue','未记录')))}；<b>阶段：</b>{html.escape(str(context.get('developmental_or_reproductive_stage','未记录')))}；<b>平台/单位：</b>{html.escape(str(context.get('platform','未记录')))} / {html.escape(str(context.get('observation_unit',cfg.get('observation_unit',''))))}</p><p><b>主要问题：</b>{html.escape('；'.join(map(str,context.get('primary_questions',[]))))}</p><p><b>分析范围：</b>full object={html.escape(str(scope_policy.get('full_object_n','未记录')))}；analysis set={html.escape(str(scope_policy.get('analysis_set_n','未记录')))}；excluded initial QC={html.escape(str(scope_policy.get('excluded_initial_qc_n','未记录')))}。初始 QC 排除仍保留在完整账本和空间图，但不进入 DEG、marker、锚点或生物学比例。</p><p><b>全特征验证：</b>{html.escape(str(feature_audit.get('status')))}；features={html.escape(str(feature_audit.get('n_features','NA')))}；profile marker coverage={html.escape(str(feature_audit.get('profile_marker_coverage','NA')))}</p><p><b>completion gate：</b>{html.escape(str(gate.get('status')))}；{html.escape('；'.join(map(str,gate.get('errors',[]))))}</p></section>
<section id='overview'><h2>{'最终' if gate.get('status') == 'PASS' else '当前'}注释总览</h2><div class='grid'>{''.join(overview_cards) or '<p>当前阶段尚未生成总览图。</p>'}</div></section>
<section id='final'><h2>单一最终注释</h2><p>本报告只发布一套标签：大类至少为中等置信度，亚群必须为高置信度。经审核的 broad-only 救回参与所属大类的最终 DEG、dotplot、UMAP 和空间统计，但不能成为亚群锚点或被强行赋予 fine label。界面、QC holdout、技术状态和初始 QC 排除不进入生物学 DEG。</p>{table(annotation_views,['view','n_observations','policy','marker_deg_eligible','status','artifact'])}<h3>最终注释 census</h3>{table(final_census,['state','broad_label','fine_label','assignment_tier','n_observations'])}</section>
<section id='selection'><h2>聚类选择与候选审计</h2>{table(clustering,['run_id','parameters','n_clusters','quantitative_rank','decision','rationale'])}</section>
<section id='tree'><h2>可展开注释树与决策</h2>{tree_overview}<p><button onclick="toggleDetails('tree',true)">全部展开</button><button onclick="toggleDetails('tree',false)">全部折叠</button><input id='tree-search' type='search' placeholder='搜索大类或亚群' oninput="filterNodes(this.value,'.tree-node')"></p>{''.join(tree_parts) or '<p>注释树将在最终空间节点资产生成后显示。</p>'}<h3>逐簇/cohort 决策账本（含历史与 supersedes）</h3>{table(clusters,['decision_id','source_run_id','source_cluster','n_observations','spatial_object_count','count_interpretation','broad_label','fine_label','state','confidence','validation_feature_scope','iteration','route','supersedes','next_action'])}</section>
<section id='routes'><h2>大类 cohort、跨谱系回归与残余 QC 裁决</h2><p><b>工作流完成状态：</b>{html.escape(str(workflow_audit.get('status','NOT_RUN')))}。确认门位于全组织开放谱系审查、每个可信初始大类的 cohort 重聚类、所有直接/跨谱系回归、必要的定向 cohort/RCTD、完整残余 QC Atlas 复核和单一最终注释之后。新流程不建立持久化生物学池，也不在 Atlas 前重聚类 QC。</p><p><b>Atlas broad rescue：</b>Atlas 只处理所有大类/定向判断完成后冻结的完整残余 QC membership。RCTD low 汇入该 membership。held-out target precision 默认 moderate-or-higher={html.escape(str(atlas_policy.get('moderate_target_precision',0.90)))}、high={html.escape(str(atlas_policy.get('high_target_precision',0.95)))}；中等及以上只能 broad-only 回归，低于中等保留 reject/QC，且 <code>fine_anchor_eligible=false</code>。</p><h3>大类与定向重聚类 cohort</h3>{table(cohorts,['cohort_id','cohort_type','source_broad_label','n_observations','candidate_resolutions','selected_resolution','terminal_outcome','applicability','status','outcome_artifact'])}<h3>直接与跨谱系回归</h3>{table(direct_returns,['return_id','source_cohort_id','source_cluster','n_observations','target_broad_label','target_fine_label','confidence','assignment_mode','rctd_tier','status'])}<h3>残余 QC 的 Atlas / 内部锚点校准结果</h3>{table(atlas_rows,['channel','n_query','high_n','moderate_only_n','moderate_or_higher_n','low_reject_n','stromal_moderate_or_higher_n','direct_writeback'])}<h3>辅助路线记录</h3>{table(route_attempts,['route_attempt_id','decision_id','source_state','route_class','failure_mode','applicability','reference_id','n_query','query_membership_sha256','calibration_origin','status','rctd_high_n','rctd_moderate_n','rctd_low_n','rctd_fine_return_n','rctd_broad_return_n','n_defined_broad_only','n_qc_retained','validation_artifact'])}<h3>尚未闭环的工作流问题</h3>{table(workflow_audit.get('gaps',[]),['code','entity_type','entity_id','required_action','detail','blocking'])}</section>
<section id='iterations'><h2>多轮 cohort、作业与下一动作</h2><h3>分析运行</h3>{table(runs,['run_id','stage','environment','scheduler_job_name','scheduler_job_id','status','output_root'])}<h3>强制下一动作队列</h3>{table(queue,['priority','source_cluster','n_observations','current_state','required_route','reason','target_scope'])}</section>
<section id='uncertainty'><h2>保留状态与未解决群体</h2><p>最终 census 只统计冻结 analysis set。没有达到中等大类置信度的 observation 保留为界面、QC、技术或待审查状态，不伪装成细胞类型。</p>{table(final_census,['state','broad_label','fine_label','assignment_tier','n_observations'])}<h3>Full-object active decision 状态</h3>{table(summary_rows,['state','n_observations'])}<h3>保留审查/排除状态的历史决策</h3>{table(unresolved,['source_cluster','n_observations','state','evidence_status','route','next_action'])}</section>
{''.join(sections)}<section id='rare-audit'><h2>Oocyte 与样本特异候选空间对象审计</h2><p>仅展示实际触发的候选验证；不存在按“稀有”强制生成的细胞类型。对象数与 cellbin 数均不等同于生物学细胞计数。</p><div class='grid'>{''.join(rare_cards) or '<p>当前项目没有需要专门空间对象审计的候选。</p>'}</div></section><section id='nodes'><h2>逐节点空间高亮</h2><p><button onclick="toggleDetails('nodes',true)">全部展开</button><button onclick="toggleDetails('nodes',false)">全部折叠</button><input type='search' placeholder='搜索节点空间图' oninput="filterNodes(this.value,'.node-card')"></p>{''.join(node_cards) or '<p>单细胞项目或当前阶段无节点空间图。</p>'}</section><section id='genes'><h2>按细胞类型分组的空间 marker</h2>{gene_sections or '<p>当前阶段尚未生成空间基因资产。</p>'}</section><section id='downloads'><h2>结果下载、运行记录与审计</h2><ul>{''.join(downloads) or '<li>当前阶段尚无发布文件。</li>'}</ul></section><section id='workflow'><h2>中文全流程详细版：从原始输入到最终发布</h2><p>下表由 workflow event registry 还原，包含输入、参数、调度作业、失败修复、生物学裁决和原子写回。最底部保留状态文件原文。</p>{table(workflow_events,['timestamp','phase','branch_id','action','input_scope','parameters','scheduler_job_id','status','decision_summary_zh','artifact','supersedes_event_id'])}<h3>状态文件原文折叠备份</h3><details><summary>展开 annotation_state.md</summary><pre>{html.escape(state)}</pre></details></section></main><script>function toggleDetails(id,open){{document.querySelectorAll('#'+id+' details').forEach(x=>x.open=open);}}function filterNodes(value,selector){{const q=value.trim().toLowerCase();document.querySelectorAll(selector).forEach(x=>{{const s=(x.dataset.search||x.textContent).toLowerCase();x.style.display=!q||s.includes(q)?'block':'none';}});}}</script></body></html>"""
    out = report_dir / "index.html"; out.write_text(body, encoding="utf-8")
    dependency_keys=["cell_ledger","cluster_ledger","completion_gate","master_quality_approval","confirmation_review_report","confirmation_review_manifest","annotation_support_registry"]
    dependencies=[root/str(confirmation[key]) for key in dependency_keys]
    dependencies.append(root/"state/final_annotation_confirmation.json")
    build_dependency_manifest(out,dependencies,{"asset_class":"final_html_report"})
    print(out); return 0


if __name__ == "__main__": raise SystemExit(main())
