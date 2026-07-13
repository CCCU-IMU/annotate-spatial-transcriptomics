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
    pools = read_tsv(root / "state/pool_registry.tsv")
    runs = read_tsv(root / "state/run_registry.tsv")
    queue = read_tsv(root / "state/next_action_queue.tsv")
    route_attempts = read_tsv(root / "state/route_attempt_registry.tsv")
    branches = read_tsv(root / "state/branch_control_board.tsv")
    workflow_events = read_tsv(root / "state/workflow_event_registry.tsv")
    annotation_views = read_tsv(root / "state/annotation_view_registry.tsv")
    multi_path = root / "provenance/multiroute_audit.json"
    multi_audit = json.loads(multi_path.read_text()) if multi_path.exists() else {"status":"NOT_RUN","gaps":[]}
    context_path = root / "config/biological_context.json"
    context = json.loads(context_path.read_text(encoding="utf-8")) if context_path.exists() else {}
    gate_path = root / "provenance/completion_gate.json"
    gate = json.loads(gate_path.read_text(encoding="utf-8")) if gate_path.exists() else {"status":"NOT_RUN","errors":["completion gate has not run"]}
    feature_path = root / "provenance/full_feature_validation.json"
    feature_audit = json.loads(feature_path.read_text(encoding="utf-8")) if feature_path.exists() else {"status":"NOT_RUN"}
    scope_path = root / "provenance/analysis_scope_policy.json"
    scope_policy = json.loads(scope_path.read_text(encoding="utf-8")) if scope_path.exists() else {}
    atlas_policy_path = root / "provenance/atlas_broad_rescue_policy_quality_reopen.json"
    atlas_policy = json.loads(atlas_policy_path.read_text(encoding="utf-8")) if atlas_policy_path.exists() else {}
    strict_census = read_tsv(root / "tables/strict_annotation_census.tsv")
    inclusive_census = read_tsv(root / "tables/inclusive_annotation_census.tsv")
    display_census = read_tsv(root / "tables/display_annotation_census.tsv")
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
    view_overview_sections = []
    for view in ["strict","inclusive","display"]:
        cards = []
        for kind,title in [
            ("broad_spatial","大类空间图"),("subtype_spatial","亚群空间图"),
            ("broad_UMAP","大类 UMAP"),("subtype_UMAP","亚群 UMAP"),
        ]:
            p=root/"figures"/f"{view}_{kind}.png"
            if p.exists():cards.append(f"<article class='card'><h3>{view} {title}</h3><img src='../figures/{p.name}'></article>")
        view_overview_sections.append(
            f"<details {'open' if view == 'display' else ''}><summary><b>{view}</b> 注释总览</summary>"
            f"<div class='grid'>{''.join(cards)}</div></details>"
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
        ("Strict broad DEG", root/"tables/strict_broad_DEG_one_vs_rest_all.tsv"), ("Strict subtype DEG", root/"tables/strict_subtype_DEG_one_vs_rest_all.tsv"),
        ("Display broad DEG", root/"tables/display_broad_DEG_one_vs_rest_all.tsv"), ("Display subtype DEG", root/"tables/display_subtype_DEG_one_vs_rest_all.tsv"),
        ("Broad DEG", root/"tables/broad_DEG_one_vs_rest_all.tsv"), ("Subtype DEG", root/"tables/subtype_DEG_one_vs_rest_all.tsv"),
        ("Broad DEG top", root/"tables/broad_DEG_top100.tsv"), ("Subtype DEG top", root/"tables/subtype_DEG_top100.tsv"),
        ("Cell ledger", root/"state/cell_ledger.tsv.gz"), ("Cell ledger", root/"tables/cell_ledger.tsv.gz"),
        ("Cluster decisions", root/"state/cluster_decision_ledger.tsv"), ("Clustering decisions", root/"state/clustering_decision_ledger.tsv"),
        ("Run registry", root/"state/run_registry.tsv"), ("Pool registry", root/"state/pool_registry.tsv"),
        ("State validation", root/"provenance/state_validation.json"), ("Release audit", root/"provenance/release_audit.json"),
        ("Biological context", context_path), ("Next-action queue", root/"state/next_action_queue.tsv"),
        ("Full-feature validation", feature_path),
        ("Analysis-scope policy", scope_path),
        ("Atlas broad-rescue policy", atlas_policy_path),
        ("Final annotation user confirmation", root/"state/final_annotation_confirmation.json"),
        ("Completion gate", gate_path), ("Iteration plan", root/"provenance/iteration_plan.json"),
        ("Multi-route audit", multi_path), ("Route attempts", root/"state/route_attempt_registry.tsv"),
        ("Branch control board", root/"state/branch_control_board.tsv"), ("Workflow events", root/"state/workflow_event_registry.tsv"),
        ("Annotation views", root/"state/annotation_view_registry.tsv"),
        ("Checksums", root/"provenance/checksums.sha256"), ("Session info", root/"provenance/release_sessionInfo.txt"),
    ]
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
    atlas_external = atlas_policy.get("current_external_atlas_result_g067", {})
    atlas_internal = atlas_policy.get("current_internal_anchor_result_g067", {})
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
    atlas_writeback_rows = [atlas_policy.get("nested_recalibration_writeback_counts", {})] if atlas_policy.get("nested_recalibration_writeback_counts") else []
    body = f"""<!doctype html><html lang='zh'><head><meta charset='utf-8'><title>{html.escape(cfg['sample_id'])} 注释报告</title>
<style>body{{font-family:Arial,'Noto Sans CJK SC',sans-serif;margin:0;background:#f5f7fa;color:#20242a}}header,main{{max-width:1500px;margin:auto;padding:24px}}header{{background:#17324d;color:white;max-width:none}}nav a{{color:#d8efff;margin-right:18px}}section{{background:white;margin:18px 0;padding:20px;border-radius:10px;scroll-margin-top:12px}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(520px,1fr));gap:16px}}.card{{border:1px solid #d9e0e8;padding:12px;border-radius:8px}}img{{width:100%;height:auto}}table{{border-collapse:collapse;width:100%;font-size:12px}}th,td{{border:1px solid #ccd4dd;padding:5px;text-align:left}}pre{{white-space:pre-wrap}}details{{margin:8px 0;padding:6px;border-left:3px solid #d9e0e8}}summary{{cursor:pointer}}.button,button{{display:inline-block;background:#176b87;color:white!important;border:0;border-radius:5px;padding:7px 10px;text-decoration:none;cursor:pointer;margin:3px}}input[type=search]{{padding:8px;min-width:320px;border:1px solid #9ca9b5;border-radius:5px}}</style></head>
<body><header><h1>{html.escape(cfg['sample_id'])} 可追溯注释报告</h1><p>{html.escape(final_status)} · {html.escape(cfg['modality'])} · observation={html.escape(cfg['observation_unit'])} · framework={html.escape(cfg['framework_version'])}</p><nav><a href='#context'>生物学背景</a><a href='#overview'>全局图</a><a href='#views'>三层注释</a><a href='#selection'>聚类选择</a><a href='#tree'>注释决策</a><a href='#routes'>多路线裁决</a><a href='#iterations'>多轮池状态</a><a href='#uncertainty'>不确定群体</a><a href='#dot-broad'>大类点图</a><a href='#dot-subtype'>亚群点图</a><a href='#rare-audit'>稀有细胞审计</a><a href='#nodes'>空间节点</a><a href='#genes'>空间基因</a><a href='#downloads'>下载与审计</a><a href='#workflow'>全流程详细版</a></nav></header><main>
<section id='context'><h2>生物学背景与完成门</h2><p><b>报告状态：</b>{html.escape(final_status)}</p><p><b>用户终审：</b>{html.escape(str(confirmation.get('status')))}；确认时间={html.escape(str(confirmation.get('confirmed_at','')))}；decision version={html.escape(str(confirmation.get('decision_version','')))}。</p><p><b>物种/组织：</b>{html.escape(str(context.get('species','未记录')))} / {html.escape(str(context.get('tissue','未记录')))}；<b>阶段：</b>{html.escape(str(context.get('developmental_or_reproductive_stage','未记录')))}；<b>平台/单位：</b>{html.escape(str(context.get('platform','未记录')))} / {html.escape(str(context.get('observation_unit',cfg.get('observation_unit',''))))}</p><p><b>主要问题：</b>{html.escape('；'.join(map(str,context.get('primary_questions',[]))))}</p><p><b>分析范围：</b>full object={html.escape(str(scope_policy.get('full_object_n','未记录')))}；analysis set={html.escape(str(scope_policy.get('analysis_set_n','未记录')))}；excluded initial QC={html.escape(str(scope_policy.get('excluded_initial_qc_n','未记录')))}。初始 QC 排除仍保留在完整账本和空间图，但不进入 DEG、marker、锚点或生物学比例。</p><p><b>全特征验证：</b>{html.escape(str(feature_audit.get('status')))}；features={html.escape(str(feature_audit.get('n_features','NA')))}；profile marker coverage={html.escape(str(feature_audit.get('profile_marker_coverage','NA')))}</p><p><b>completion gate：</b>{html.escape(str(gate.get('status')))}；{html.escape('；'.join(map(str,gate.get('errors',[]))))}</p></section>
<section id='overview'><h2>{'最终' if gate.get('status') == 'PASS' else '当前'}注释总览</h2><div class='grid'>{''.join(overview_cards) or '<p>当前阶段尚未生成总览图。</p>'}</div></section>
<section id='views'><h2>Strict / Inclusive / Display 三层注释</h2><p>Strict 用于保守生物学定义、锚点和 marker/DEG 证据；Inclusive 接收经校准的 broad-only rescue；Display 是人工裁决后的主报告层。三层必须分别覆盖 analysis set，初始 QC 排除另存于 full-object ledger。</p>{table(annotation_views,['view','n_observations','policy','marker_deg_eligible','status','artifact'])}{''.join(view_overview_sections)}<h3>Strict 证据层 census</h3>{table(strict_census,['state','broad_label','fine_label','n_observations'])}<h3>Inclusive 主解释层 census</h3>{table(inclusive_census,['state','broad_label','fine_label','n_observations'])}</section>
<section id='selection'><h2>聚类选择与候选审计</h2>{table(clustering,['run_id','parameters','n_clusters','quantitative_rank','decision','rationale'])}</section>
<section id='tree'><h2>可展开注释树与决策</h2>{tree_overview}<p><button onclick="toggleDetails('tree',true)">全部展开</button><button onclick="toggleDetails('tree',false)">全部折叠</button><input id='tree-search' type='search' placeholder='搜索大类或亚群' oninput="filterNodes(this.value,'.tree-node')"></p>{''.join(tree_parts) or '<p>注释树将在最终空间节点资产生成后显示。</p>'}<h3>逐簇/逐池决策账本（含历史与 supersedes）</h3>{table(clusters,['decision_id','source_run_id','source_cluster','n_observations','spatial_object_count','count_interpretation','broad_label','fine_label','state','confidence','validation_feature_scope','iteration','route','supersedes','next_action'])}</section>
<section id='routes'><h2>多路线待定义细胞裁决</h2><p><b>多路线完成状态：</b>{html.escape(str(multi_audit.get('status','NOT_RUN')))}。大型生物学池必须完成 balanced-anchor query-only 重聚类；局部界面必须先完成无监督复核，再将 RCTD 作为低优先级分层证据；大型聚类后 QC 必须先完成完整 QC 大池锚点重聚类，再处理 atlas rejects。</p><p><b>Atlas broad rescue：</b>held-out target precision 默认 moderate-or-higher={html.escape(str(atlas_policy.get('moderate_target_precision',0.90)))}、high={html.escape(str(atlas_policy.get('high_target_precision',0.95)))}，并使用嵌套阈值；每个 high 都同时达到 moderate-or-higher 门槛。互斥结果应分别报告 high、moderate-only 和 low-reject，同时报告 high + moderate-only 的累计数。中等及以上均可成为 broad-only 回归候选，低于中等保留 reject/review。上述数值是校准精度目标，不是单个 observation 的原始 confidence 阈值；外部 Atlas 必须与本样本 marker/anti-marker、路线及空间或内部锚点证据联合，所有回归均 <code>fine_anchor_eligible=false</code>。</p><h3>嵌套 Atlas / 内部锚点校准结果</h3>{table(atlas_rows,['channel','n_query','high_n','moderate_only_n','moderate_or_higher_n','low_reject_n','stromal_moderate_or_higher_n','direct_writeback'])}<h3>v073 联合裁决实际写回</h3>{table(atlas_writeback_rows,['granulosa_broad_only','theca_broad_only','strict_changes','vascular_candidates_rejected','legacy_stromal_labels_retained'])}<h3>Route attempts</h3>{table(route_attempts,['route_attempt_id','decision_id','pool_snapshot_id','route_class','failure_mode','applicability','n_query','n_anchors','query_only_graph','depth_matched_validation','observed_density_spatial_prior','selected_resolution','status','rctd_extreme_n','rctd_high_n','rctd_medium_low_n','rctd_fine_return_n','rctd_broad_return_n','fallback_route_attempt_id','n_rerouted','n_interface_retained','n_qc_retained','validation_artifact'])}<h3>尚未闭环的路线</h3>{table(multi_audit.get('gaps',[]),['decision_id','n_observations','required_route','reason'])}<h3>Branch control / no-repeat</h3>{table(branches,['branch_id','parent_decision_id','pool_snapshot_id','generation','run_id','selected_resolution','n_query','current_state','recluster_policy','terminal','next_action','authoritative_artifact'])}</section>
<section id='iterations'><h2>多轮池、作业与下一动作</h2><h3>池注册表</h3>{table(pools,['pool_id','parent_pool_id','n_observations','purpose','status','decision_version'])}<h3>分析运行</h3>{table(runs,['run_id','stage','environment','scheduler_job_id','status','output_root'])}<h3>强制下一动作队列</h3>{table(queue,['priority','source_cluster','n_observations','current_state','required_route','reason','target_pool'])}</section>
<section id='uncertainty'><h2>状态数量与未解决群体</h2><p>下方 display census 只统计冻结的 analysis set；随后 active-decision 表保留 full-object 历史状态，因此可能包含已划入 excluded-initial-QC 的 observation，二者不可混为生物学未注释率。</p>{table(display_census,['state','broad_label','fine_label','n_observations'])}<h3>Full-object active decision 状态</h3>{table(summary_rows,['state','n_observations'])}<h3>保留审查/排除状态的历史决策</h3>{table(unresolved,['source_cluster','n_observations','state','evidence_status','route','next_action'])}</section>
{''.join(sections)}<section id='rare-audit'><h2>稀有细胞空间对象审计</h2><p>此处对象数与 cellbin 数均不等同于生物学细胞计数；只有通过 completion gate 的标签才进入最终注释。</p><div class='grid'>{''.join(rare_cards) or '<p>当前项目无稀有细胞对象审计图。</p>'}</div></section><section id='nodes'><h2>逐节点空间高亮</h2><p><button onclick="toggleDetails('nodes',true)">全部展开</button><button onclick="toggleDetails('nodes',false)">全部折叠</button><input type='search' placeholder='搜索节点空间图' oninput="filterNodes(this.value,'.node-card')"></p>{''.join(node_cards) or '<p>单细胞项目或当前阶段无节点空间图。</p>'}</section><section id='genes'><h2>按细胞类型分组的空间 marker</h2>{gene_sections or '<p>当前阶段尚未生成空间基因资产。</p>'}</section><section id='downloads'><h2>结果下载、运行记录与审计</h2><ul>{''.join(downloads) or '<li>当前阶段尚无发布文件。</li>'}</ul></section><section id='workflow'><h2>中文全流程详细版：从原始输入到最终发布</h2><p>下表由 workflow event registry 还原，包含输入、参数、调度作业、失败修复、生物学裁决和原子写回。最底部保留状态文件原文。</p>{table(workflow_events,['timestamp','phase','branch_id','action','input_scope','parameters','scheduler_job_id','status','decision_summary_zh','artifact','supersedes_event_id'])}<h3>状态文件原文折叠备份</h3><details><summary>展开 annotation_state.md</summary><pre>{html.escape(state)}</pre></details></section></main><script>function toggleDetails(id,open){{document.querySelectorAll('#'+id+' details').forEach(x=>x.open=open);}}function filterNodes(value,selector){{const q=value.trim().toLowerCase();document.querySelectorAll(selector).forEach(x=>{{const s=(x.dataset.search||x.textContent).toLowerCase();x.style.display=!q||s.includes(q)?'block':'none';}});}}</script></body></html>"""
    out = report_dir / "index.html"; out.write_text(body, encoding="utf-8"); print(out); return 0


if __name__ == "__main__": raise SystemExit(main())
