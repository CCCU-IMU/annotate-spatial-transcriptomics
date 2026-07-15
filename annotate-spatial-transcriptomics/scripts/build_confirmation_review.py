#!/usr/bin/env python3
"""Build the lightweight evidence report shown before final user confirmation."""

from __future__ import annotations

import argparse
import base64
import csv
import gzip
import hashlib
import html
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from master_quality_lib import sha256 as master_sha256, validate_master_approval


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def resolve(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def active_rows(rows: list[dict[str, str]], id_field: str) -> list[dict[str, str]]:
    superseded: set[str] = set()
    for row in rows:
        superseded.update(x for x in row.get("supersedes", "").replace(",", ";").split(";") if x)
    return [row for row in rows if row.get(id_field, "") not in superseded]


def table(rows: list[dict[str, object]], columns: list[str]) -> str:
    if not rows:
        return "<p>无。</p>"
    head = "".join(f"<th>{html.escape(column)}</th>" for column in columns)
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(row.get(column, '')))}</td>" for column in columns) + "</tr>"
        for row in rows
    )
    return f"<div class='scroll'><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--asset-manifest", type=Path)
    parser.add_argument("--support-registry", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    cell_ledger = root / "state/cell_ledger.tsv.gz"
    cluster_ledger = root / "state/cluster_decision_ledger.tsv"
    completion_gate = root / "provenance/completion_gate.json"
    master_approval_path = root / "state/master_quality_approval.json"
    support_registry = args.support_registry or root / "state/annotation_support_registry.tsv"
    asset_manifest_path = args.asset_manifest or root / "review/confirmation/assets/review_asset_manifest.json"
    output = args.out or root / "review/confirmation/index.html"
    required = [cell_ledger, cluster_ledger, completion_gate, support_registry, asset_manifest_path, master_approval_path]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise SystemExit("confirmation review inputs are missing: " + ", ".join(missing))
    completion = json.loads(completion_gate.read_text(encoding="utf-8"))
    if completion.get("status") != "PASS":
        raise SystemExit("completion gate must pass before building the confirmation review")
    master_ok, master_errors, master_approval = validate_master_approval(root, master_approval_path)
    if not master_ok:
        raise SystemExit("main-Agent annotation-quality approval is required: " + "; ".join(master_errors))

    cells = read_tsv(cell_ledger)
    broad_counts: Counter[str] = Counter()
    fine_counts: Counter[tuple[str, str]] = Counter()
    retained_counts: Counter[str] = Counter()
    for row in cells:
        if row.get("analysis_scope") != "analysis_set":
            continue
        broad = row.get("final_broad_label", "")
        fine = row.get("final_fine_label", "")
        if broad:
            broad_counts[broad] += 1
            if fine:
                fine_counts[(broad, fine)] += 1
        else:
            retained_counts[row.get("final_state", "") or "pending_review"] += 1
    if not broad_counts:
        raise SystemExit("confirmation review has no final broad labels")

    support_rows = active_rows(read_tsv(support_registry), "support_id")
    support_rows = [row for row in support_rows if row.get("status") == "validated"]
    broad_support = {row.get("broad_label", ""): row for row in support_rows if row.get("label_level") == "broad"}
    fine_support = {
        (row.get("broad_label", ""), row.get("fine_label", "")): row
        for row in support_rows if row.get("label_level") == "fine"
    }
    missing_broad = sorted(set(broad_counts) - set(broad_support))
    missing_fine = sorted(set(fine_counts) - set(fine_support))
    if missing_broad or missing_fine:
        raise SystemExit(f"validated support reasons are incomplete: broad={missing_broad}; fine={missing_fine}")

    asset_manifest = json.loads(asset_manifest_path.read_text(encoding="utf-8"))
    if asset_manifest.get("status") != "PASS" or asset_manifest.get("label_column") != "primary_broad_label":
        raise SystemExit("review assets did not pass or were not built from primary_broad_label")
    asset_keys = ["spatial_png", "dotplot_png", "dotplot_source", "palette_tsv"]
    assets: dict[str, Path] = {}
    for key in asset_keys:
        path = resolve(root, str(asset_manifest.get(key, "")))
        if not path.is_file() or sha256(path) != asset_manifest.get(f"{key}_sha256"):
            raise SystemExit(f"review asset is missing or stale: {key}")
        assets[key] = path
    palette_rows = read_tsv(assets["palette_tsv"])
    palette = {row.get("label", ""): row.get("color", "") for row in palette_rows}
    if set(palette) != set(broad_counts):
        raise SystemExit("review palette labels do not exactly match final broad labels")
    colors = [value.upper() for value in palette.values()]
    if len(colors) != len(set(colors)) or any(len(value) != 7 or not value.startswith("#") for value in colors):
        raise SystemExit("review palette colors must be unique #RRGGBB values")

    def png_data_uri(path: Path) -> str:
        return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode("ascii")

    broad_rows = []
    for label, count in broad_counts.most_common():
        support = broad_support[label]
        broad_rows.append({
            "大类": label, "数量": count, "置信度": support.get("confidence", ""),
            "正向 marker": support.get("positive_marker_evidence", ""),
            "anti-marker": support.get("anti_marker_evidence", ""),
            "分辨率/稳定性": support.get("resolution_evidence", ""),
            "空间依据": support.get("spatial_evidence", ""),
            "路线": support.get("route_summary", ""),
        })
    fine_rows = []
    for (broad, fine), count in sorted(fine_counts.items(), key=lambda item: -item[1]):
        support = fine_support[(broad, fine)]
        fine_rows.append({
            "大类": broad, "高置信亚群": fine, "数量": count,
            "正向 marker": support.get("positive_marker_evidence", ""),
            "anti-marker": support.get("anti_marker_evidence", ""),
            "分辨率/稳定性": support.get("resolution_evidence", ""),
            "空间依据": support.get("spatial_evidence", ""),
        })
    palette_html = "".join(
        f"<span class='chip'><i style='background:{html.escape(palette[label])}'></i>{html.escape(label)} ({broad_counts[label]:,})</span>"
        for label in broad_counts
    )
    quality_rows = []
    for item, entry in master_approval.get("checklist", {}).items():
        quality_rows.append({"质量维度": item, "结论": entry.get("status", ""), "备注": entry.get("note", entry.get("evidence", ""))})
    review_html = f"""<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>确认前轻量注释审阅</title>
<style>body{{font-family:Arial,'Noto Sans CJK SC',sans-serif;background:#f4f6f8;color:#20242a;margin:0}}header,main{{max-width:1500px;margin:auto;padding:22px}}header{{max-width:none;background:#17324d;color:#fff}}section{{background:#fff;margin:16px 0;padding:18px;border-radius:10px}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(520px,1fr));gap:15px}}img{{max-width:100%;height:auto;border:1px solid #d8dee6}}table{{border-collapse:collapse;width:100%;font-size:12px}}th,td{{border:1px solid #cdd5df;padding:6px;text-align:left;vertical-align:top}}th{{background:#eef3f7}}.scroll{{overflow:auto}}.chip{{display:inline-flex;align-items:center;margin:4px 10px 4px 0;padding:5px 8px;border:1px solid #ccd4dd;border-radius:16px}}.chip i{{width:16px;height:16px;border-radius:3px;margin-right:6px;border:1px solid #555}}.warn{{background:#fff4ce;border-left:5px solid #d99b00;padding:12px}}</style></head><body>
<header><h1>确认前轻量注释审阅</h1><p>这是完成门通过后的冻结候选注释，仅供人工确认；尚未生成最终 DEG、完整双层树点图、逐节点空间图或发布 HTML。</p></header><main>
<section class='warn'><b>确认目标：</b>中等及以上置信大类 + 仅高置信亚群。请重点检查大类数量、支持/反证、空间分布及典型 marker；若需修改，返回迭代路线而不是直接确认。</section>
<section><h2>主 Agent 注释质量审批</h2><p><b>结论：</b>{html.escape(master_approval.get('status',''))}　<b>与验证参考流程的比较：</b>{html.escape(master_approval.get('comparison_to_reference',''))}</p><p><b>审批说明：</b>{html.escape(master_approval.get('rationale',''))}</p>{table(quality_rows, ['质量维度','结论','备注'])}</section>
<section><h2>大类配色与数量</h2>{palette_html}</section>
<section class='grid'><article><h2>大类空间投影</h2><img src='{png_data_uri(assets['spatial_png'])}'></article><article><h2>大类典型 marker dotplot</h2><img src='{png_data_uri(assets['dotplot_png'])}'></article></section>
<section><h2>大类注释支持原因</h2>{table(broad_rows, ['大类','数量','置信度','正向 marker','anti-marker','分辨率/稳定性','空间依据','路线'])}</section>
<section><h2>高置信亚群支持原因</h2>{table(fine_rows, ['大类','高置信亚群','数量','正向 marker','anti-marker','分辨率/稳定性','空间依据'])}</section>
<section><h2>保留的非生物学状态</h2>{table([{'状态': key, '数量': value} for key, value in retained_counts.most_common()], ['状态','数量'])}</section>
</main></body></html>"""
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(review_html, encoding="utf-8")
    manifest = {
        "status": "PASS", "created_at": datetime.now(timezone.utc).isoformat(),
        "report": str(output.resolve()), "report_sha256": sha256(output),
        "cell_ledger": str(cell_ledger.resolve()), "cell_ledger_sha256": sha256(cell_ledger),
        "cluster_ledger": str(cluster_ledger.resolve()), "cluster_ledger_sha256": sha256(cluster_ledger),
        "completion_gate": str(completion_gate.resolve()), "completion_gate_sha256": sha256(completion_gate),
        "master_quality_approval": str(master_approval_path.resolve()), "master_quality_approval_sha256": master_sha256(master_approval_path),
        "support_registry": str(support_registry.resolve()), "support_registry_sha256": sha256(support_registry),
        "asset_manifest": str(asset_manifest_path.resolve()), "asset_manifest_sha256": sha256(asset_manifest_path),
        "broad_labels": dict(broad_counts),
        "fine_labels": {f"{broad} -> {fine}": count for (broad, fine), count in fine_counts.items()},
        "retained_states": dict(retained_counts),
        "scope": "preconfirmation_lightweight_only_no_final_deg_or_release_assets",
    }
    manifest_path = root / "provenance/confirmation_review_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
