#!/usr/bin/env python3
"""Fail-closed audit for a framework release directory."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from html.parser import HTMLParser
from pathlib import Path


def read_tsv(path: Path):
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", newline="", encoding="utf-8") as handle: return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def is_within(path: Path, parent: Path) -> bool:
    """Python 3.8-compatible equivalent of Path.is_relative_to()."""
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


class Links(HTMLParser):
    def __init__(self): super().__init__(); self.paths = []
    def handle_starttag(self, tag, attrs):
        key = "src" if tag == "img" else "href" if tag == "a" else None
        if key:
            value = dict(attrs).get(key)
            if value and not value.startswith(("#", "http://", "https://", "file://")): self.paths.append(value)


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("release_root", type=Path); ap.add_argument("--profile", choices=["core", "full"], default="full"); args = ap.parse_args()
    root = args.release_root.resolve(); errors = []
    def req(ok, msg):
        if not ok: errors.append(msg)
    report = root / "report/index.html"
    req(report.exists() and report.stat().st_size > 2000, "missing or undersized HTML report")
    index = root / "figures/marker_dotplots/marker_dotplot_asset_index.tsv"
    req(index.exists(), "missing dotplot asset index")
    assets = read_tsv(index) if index.exists() else []
    req({"broad", "subtype"}.issubset({r.get("level") for r in assets}), "both broad and subtype dotplots are mandatory")
    combos = {(r.get("level"), r.get("panel")) for r in assets}
    req({("broad", "canonical"), ("subtype", "canonical")}.issubset(combos), "canonical broad and subtype dotplots are mandatory")
    for r in assets:
        for key in ["png", "pdf", "source"]:
            p = Path(r.get(key, "")); p = p if p.is_absolute() else root / p
            req(p.exists() and p.stat().st_size > (5000 if key != "source" else 100), f"invalid dotplot {key}: {p}")
        src = Path(r.get("source", "")); src = src if src.is_absolute() else root / src
        if src.exists():
            rows = read_tsv(src)
            required = {"gene", "label", "avg_expression", "pct_expressed_absolute", "n_observations", "marker_group", "avg_expression_scaled_within_gene", "pct_expressed_scaled_within_gene", "analysis_view", "evidence_cohort"}
            req(bool(rows) and required.issubset(rows[0]), f"dotplot source schema invalid: {src}")
            for x in rows:
                try:
                    vals = [float(x[k]) for k in ["avg_expression", "pct_expressed_absolute", "avg_expression_scaled_within_gene", "pct_expressed_scaled_within_gene"]]
                    req(all(math.isfinite(v) for v in vals), f"non-finite dotplot value: {src}")
                    req(0 <= vals[1] <= 100 and 0 <= vals[3] <= 100, f"dotplot percentage outside 0..100: {src}")
                except Exception: req(False, f"unparseable dotplot values: {src}"); break
    for name in ["state/annotation_state.md", "state/clustering_decision_ledger.tsv", "state/cluster_decision_ledger.tsv", "provenance/state_validation.json"]:
        req((root / name).exists(), f"missing state artifact: {name}")
    if args.profile == "full":
        config_path = root / "config/project.json"
        req(config_path.exists(), "missing project config")
        config = json.loads(config_path.read_text()) if config_path.exists() else {}
        def any_file(pattern, min_bytes=100):
            return any(p.is_file() and p.stat().st_size >= min_bytes for p in root.glob(pattern))
        context_path = root / "config/biological_context.json"
        context_validation = root / "provenance/biological_context_validation.json"
        completion_path = root / "provenance/completion_gate.json"
        confirmation_path = root / "state/final_annotation_confirmation.json"
        req(context_path.exists(), "missing biological context")
        req(context_validation.exists(), "missing biological-context validation")
        if context_validation.exists():
            req(json.loads(context_validation.read_text()).get("status") == "PASS", "biological-context validation did not pass")
        req(completion_path.exists(), "missing completion gate")
        if completion_path.exists():
            req(json.loads(completion_path.read_text()).get("status") == "PASS", "completion gate did not pass")
        req(confirmation_path.exists(), "missing explicit final-annotation user confirmation")
        if confirmation_path.exists():
            confirmation = json.loads(confirmation_path.read_text())
            req(confirmation.get("status") == "CONFIRMED", "final annotation was not explicitly confirmed by the user")
            for key, hash_key in (
                ("cell_ledger", "cell_ledger_sha256"),
                ("cluster_ledger", "cluster_ledger_sha256"),
                ("completion_gate", "completion_gate_sha256"),
            ):
                target = root / str(confirmation.get(key, ""))
                req(target.is_file(), f"confirmed snapshot target is missing: {key}")
                if target.is_file():
                    req(sha256(target) == confirmation.get(hash_key), f"confirmed snapshot is stale: {key}")
        if config.get("multi_route_completion_required"):
            multi_path = root / "provenance/multiroute_audit.json"
            req(multi_path.exists(), "missing multi-route audit")
            if multi_path.exists(): req(json.loads(multi_path.read_text()).get("status") == "PASS", "multi-route audit did not pass")
            for name in ["state/route_attempt_registry.tsv", "state/branch_control_board.tsv", "state/workflow_event_registry.tsv", "state/annotation_view_registry.tsv"]:
                req((root/name).exists(), f"missing multi-route registry: {name}")
            route_rows = read_tsv(root/"state/route_attempt_registry.tsv") if (root/"state/route_attempt_registry.tsv").exists() else []
            event_rows = read_tsv(root/"state/workflow_event_registry.tsv") if (root/"state/workflow_event_registry.tsv").exists() else []
            view_rows = read_tsv(root/"state/annotation_view_registry.tsv") if (root/"state/annotation_view_registry.tsv").exists() else []
            req(bool(route_rows), "multi-route release has no route attempts")
            req(bool(event_rows), "multi-route release has no workflow events")
            req({"strict","inclusive","display"}.issubset({x.get("view") for x in view_rows if x.get("status")=="validated"}), "strict/inclusive/display views are not all validated")
            req(all(x.get("analysis_view")=="strict" for x in assets), "mandatory marker dotplots must use the strict evidence view")
            scope_path = root / "provenance/analysis_scope_policy.json"
            req(scope_path.exists(), "missing frozen analysis-scope policy")
            if scope_path.exists():
                scope = json.loads(scope_path.read_text())
                req(scope.get("status") == "PASS", "analysis-scope policy did not pass")
                membership_value = scope.get("membership_path", "")
                membership = Path(membership_value)
                membership = membership if membership.is_absolute() else root / membership
                req(membership.is_file(), "analysis-scope membership artifact is missing")
                if membership.is_file():
                    req(sha256(membership) == scope.get("membership_sha256"), "analysis-scope membership hash is stale")
            for view in ["strict","inclusive","display"]:
                for level in ["broad","subtype"]:
                    stem=f"{view}_{level}_UMAP";req((root/"figures"/f"{stem}.png").exists() and (root/"figures"/f"{stem}.pdf").exists(),f"missing annotation-view overview: {stem}")
                    if config.get("modality")=="spatial":
                        stem=f"{view}_{level}_spatial";req((root/"figures"/f"{stem}.png").exists() and (root/"figures"/f"{stem}.pdf").exists(),f"missing annotation-view spatial overview: {stem}")
            for level in ["broad","subtype"]:
                named=any_file(f"tables/strict_{level}_DEG_one_vs_rest_all.tsv*")
                generic=list(root.glob(f"tables/{level}_DEG_one_vs_rest_all.tsv*"))
                generic_strict=False
                for candidate in generic:
                    if not candidate.is_file() or candidate.stat().st_size < 100: continue
                    try:
                        rows=read_tsv(candidate)
                        generic_strict=bool(rows) and all(row.get("analysis_view")=="strict" for row in rows[:min(100,len(rows))])
                    except (OSError, UnicodeDecodeError):
                        generic_strict=False
                    if generic_strict: break
                req(named or generic_strict,f"missing strict {level} DEG with explicit analysis_view provenance")
            expected_labels = {}
            for level in ["broad", "subtype"]:
                candidates = [root/f"tables/strict_{level}_DEG_one_vs_rest_all.tsv", root/f"tables/{level}_DEG_one_vs_rest_all.tsv"]
                deg_path = next((path for path in candidates if path.is_file()), None)
                if deg_path:
                    deg_rows = read_tsv(deg_path)
                    expected_labels[level] = {row.get("label", "") for row in deg_rows if row.get("label", "") and row.get("analysis_view") == "strict"}
            for asset in assets:
                level = asset.get("level", "")
                expected = expected_labels.get(level, set())
                source = Path(asset.get("source", "")); source = source if source.is_absolute() else root/source
                if not expected or not source.is_file():
                    continue
                source_rows = read_tsv(source)
                observed_labels = {row.get("label", "") for row in source_rows if row.get("label", "")}
                marker_groups = {row.get("marker_group", "") for row in source_rows if row.get("marker_group", "")}
                req(observed_labels == expected, f"{level}/{asset.get('panel')} dotplot does not cover every strict DEG label")
                req(expected.issubset(marker_groups), f"{level}/{asset.get('panel')} marker groups do not cover every current label")
                try:
                    req(int(asset.get("n_labels", 0)) == len(expected), f"{level}/{asset.get('panel')} n_labels disagrees with strict DEG")
                except (TypeError, ValueError):
                    req(False, f"{level}/{asset.get('panel')} n_labels is invalid")
        req({("broad", "canonical"), ("broad", "data_specific"), ("subtype", "canonical"), ("subtype", "data_specific")}.issubset(combos), "full release requires canonical and data-specific dotplots at broad and subtype levels")
        if not config.get("multi_route_completion_required"):
            req(any_file("tables/broad_DEG_one_vs_rest_all.tsv*"), "missing broad one-vs-rest DEG")
            req(any_file("tables/subtype_DEG_one_vs_rest_all.tsv*"), "missing subtype one-vs-rest DEG")
        req(any_file("state/cell_ledger.tsv*", 1000) or any_file("tables/cell_ledger.tsv*", 1000) or any_file("tables/cell_metadata*.tsv*", 1000), "missing cell-level ledger")
        for name in ["state/pool_registry.tsv","state/run_registry.tsv","state/next_action_queue.tsv","state/final_annotation_confirmation.json","provenance/full_feature_validation.json","provenance/release_manifest.tsv","provenance/checksums.sha256","provenance/release_sessionInfo.txt"]:
            req((root/name).exists(),f"missing full-release provenance: {name}")
        feature_validation = root / "provenance/full_feature_validation.json"
        if feature_validation.exists():
            req(json.loads(feature_validation.read_text()).get("status") == "PASS", "full-feature validation did not pass")
        session = root / "provenance/release_sessionInfo.txt"
        if session.exists():
            req(session.stat().st_size > 100, "release session info is empty or undersized")
            release_evidence=[root/"tables/broad_DEG_one_vs_rest_all.tsv",root/"tables/subtype_DEG_one_vs_rest_all.tsv",root/"figures/marker_dotplots/marker_dotplot_asset_index.tsv"]
            req(session.stat().st_mtime >= max((p.stat().st_mtime for p in release_evidence if p.exists()),default=0), "release session info is older than final DEG/dotplot evidence")
        checksums = root / "provenance/checksums.sha256"
        if checksums.exists():
            checksum_rows = []
            for number, line in enumerate(checksums.read_text(encoding="utf-8").splitlines(), 1):
                if not line.strip(): continue
                parts = line.split(maxsplit=1)
                req(len(parts) == 2, f"malformed checksum line {number}")
                if len(parts) != 2: continue
                expected, relative = parts[0], parts[1].lstrip("*")
                target = (root / relative).resolve()
                req(len(expected) == 64 and all(c in "0123456789abcdefABCDEF" for c in expected), f"invalid SHA256 at line {number}")
                req(is_within(target, root), f"checksum path escapes release root: {relative}")
                req(target.is_file(), f"checksum target missing: {relative}")
                if target.is_file() and is_within(target, root):
                    req(sha256(target).lower() == expected.lower(), f"checksum mismatch: {relative}")
                checksum_rows.append(relative)
            req(len(checksum_rows) >= 20, "checksum manifest is unexpectedly small")
            checksum_set = set(checksum_rows)
            req("report/index.html" in checksum_set, "checksums do not cover the HTML report")
            req("state/cluster_decision_ledger.tsv" in checksum_set, "checksums do not cover the cluster ledger")
            req(any(x.startswith("state/cell_ledger.tsv") for x in checksum_set), "checksums do not cover the cell ledger")
            broad_deg_prefixes=["tables/strict_broad_DEG_one_vs_rest_all.tsv","tables/broad_DEG_one_vs_rest_all.tsv"] if config.get("multi_route_completion_required") else ["tables/broad_DEG_one_vs_rest_all.tsv"]
            subtype_deg_prefixes=["tables/strict_subtype_DEG_one_vs_rest_all.tsv","tables/subtype_DEG_one_vs_rest_all.tsv"] if config.get("multi_route_completion_required") else ["tables/subtype_DEG_one_vs_rest_all.tsv"]
            req(any(any(x.startswith(prefix) for prefix in broad_deg_prefixes) for x in checksum_set), "checksums do not cover broad DEG")
            req(any(any(x.startswith(prefix) for prefix in subtype_deg_prefixes) for x in checksum_set), "checksums do not cover subtype DEG")
            req(any(x.endswith("_source.tsv") and "marker_dotplots/" in x for x in checksum_set), "checksums do not cover dotplot source tables")
        manifest = root / "provenance/release_manifest.tsv"
        if manifest.exists() and checksums.exists():
            manifest_rows = read_tsv(manifest)
            manifest_map = {row.get("relative_path", ""): (row.get("size_bytes", ""), row.get("sha256", "")) for row in manifest_rows}
            req(set(manifest_map) == set(checksum_rows), "release manifest and checksum file list differ")
            for relative, (size, value) in manifest_map.items():
                target = root / relative
                req(target.is_file(), f"release manifest target missing: {relative}")
                if target.is_file():
                    req(str(target.stat().st_size) == str(size), f"release manifest size mismatch: {relative}")
                    req(sha256(target).lower() == value.lower(), f"release manifest SHA256 mismatch: {relative}")
        for stem in ["final_broad_UMAP", "final_subtype_UMAP"]:
            req((root/"figures"/f"{stem}.png").exists() and (root/"figures"/f"{stem}.pdf").exists(), f"missing overview pair: {stem}")
        if config.get("modality") == "spatial":
            for stem in ["final_broad_spatial", "final_subtype_spatial"]:
                req((root/"figures"/f"{stem}.png").exists() and (root/"figures"/f"{stem}.pdf").exists(), f"missing spatial overview pair: {stem}")
            req((root/"tables/spatial_node_asset_index.tsv").exists(), "missing spatial node index")
            req((root/"tables/spatial_gene_asset_index.tsv").exists(), "missing spatial gene index")
    if report.exists():
        report_text=report.read_text(encoding="utf-8")
        if args.profile=="full" and (root/"config/project.json").exists() and json.loads((root/"config/project.json").read_text()).get("multi_route_completion_required"):
            for section in ["id='views'","id='routes'","id='workflow'"]:req(section in report_text,f"report lacks required multi-route section {section}")
        parser = Links(); parser.feed(report.read_text(encoding="utf-8"))
        for value in parser.paths:
            req((report.parent / value).resolve().exists(), f"broken report link: {value}")
    result = {"status": "PASS" if not errors else "FAIL", "profile": args.profile, "errors": sorted(set(errors)), "dotplot_assets": len(assets)}
    out = root / "provenance/release_audit.json"; out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if not errors else 1


if __name__ == "__main__": raise SystemExit(main())
