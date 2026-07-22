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

from audit_annotation_membership_partition import audit as audit_membership_partition
from dependency_manifest import build as build_dependency_manifest
from validate_annotation_support_registry import validate as validate_support_registry


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
    config_path = root / "config/project.json"
    config = json.loads(config_path.read_text()) if config_path.exists() else {}
    workflow_required = config.get(
        "annotation_workflow_completion_required",
        config.get("multi_route_completion_required", False),
    )
    final_census_path = root / "tables/final_annotation_census.tsv"
    final_census = read_tsv(final_census_path) if final_census_path.exists() else []
    has_final_fine = any(row.get("fine_label", "") and int(float(row.get("n_observations", 0) or 0)) > 0 for row in final_census)
    subtype_required = has_final_fine if workflow_required else True
    def req(ok, msg):
        if not ok: errors.append(msg)
    report = root / "report/index.html"
    req(report.exists() and report.stat().st_size > 2000, "missing or undersized HTML report")
    index = root / "figures/marker_dotplots/marker_dotplot_asset_index.tsv"
    req(index.exists(), "missing dotplot asset index")
    assets = read_tsv(index) if index.exists() else []
    required_levels = {"broad"} | ({"subtype"} if subtype_required else set())
    req(required_levels.issubset({r.get("level") for r in assets}), "required broad/high-confidence subtype dotplots are missing")
    combos = {(r.get("level"), r.get("panel")) for r in assets}
    required_canonical = {("broad", "canonical")} | ({("subtype", "canonical")} if subtype_required else set())
    req(required_canonical.issubset(combos), "required canonical broad/high-confidence subtype dotplots are missing")
    for r in assets:
        for key in ["png", "pdf", "absolute_png", "absolute_pdf", "source"]:
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
        req(config_path.exists(), "missing project config")
        def any_file(pattern, min_bytes=100):
            return any(p.is_file() and p.stat().st_size >= min_bytes for p in root.glob(pattern))
        context_path = root / "config/biological_context.json"
        context_validation = root / "provenance/biological_context_validation.json"
        completion_path = root / "provenance/completion_gate.json"
        confirmation_path = root / "state/final_annotation_confirmation.json"
        master_approval_path = root / "state/master_quality_approval.json"
        req(context_path.exists(), "missing biological context")
        req(context_validation.exists(), "missing biological-context validation")
        if context_validation.exists():
            req(json.loads(context_validation.read_text()).get("status") == "PASS", "biological-context validation did not pass")
        req(completion_path.exists(), "missing completion gate")
        if completion_path.exists():
            req(json.loads(completion_path.read_text()).get("status") == "PASS", "completion gate did not pass")
        req(confirmation_path.exists(), "missing explicit final-annotation user confirmation")
        req(master_approval_path.exists(), "missing post-completion main-Agent annotation-quality approval")
        if master_approval_path.exists():
            from master_quality_lib import validate_master_approval
            master_ok, master_errors, _ = validate_master_approval(root, master_approval_path)
            req(master_ok, "invalid main-Agent annotation-quality approval: " + "; ".join(master_errors))
        if confirmation_path.exists():
            confirmation = json.loads(confirmation_path.read_text())
            req(confirmation.get("status") == "CONFIRMED", "final annotation was not explicitly confirmed by the user")
            for key, hash_key in (
                ("cell_ledger", "cell_ledger_sha256"),
                ("cluster_ledger", "cluster_ledger_sha256"),
                ("completion_gate", "completion_gate_sha256"),
                ("confirmation_review_report", "confirmation_review_report_sha256"),
                ("confirmation_review_manifest", "confirmation_review_manifest_sha256"),
                ("annotation_support_registry", "annotation_support_registry_sha256"),
                ("release_taxonomy_audit", "release_taxonomy_audit_sha256"),
                ("master_quality_approval", "master_quality_approval_sha256"),
            ):
                target = root / str(confirmation.get(key, ""))
                req(target.is_file(), f"confirmed snapshot target is missing: {key}")
                if target.is_file():
                    req(sha256(target) == confirmation.get(hash_key), f"confirmed snapshot is stale: {key}")
        if workflow_required:
            req(final_census_path.exists() and bool(final_census), "missing or empty single-final-annotation census")
            workflow_path = root / "provenance/direct_lineage_workflow_audit.json"
            req(workflow_path.exists(), "missing direct-lineage workflow audit")
            if workflow_path.exists(): req(json.loads(workflow_path.read_text()).get("status") == "PASS", "direct-lineage workflow audit did not pass")
            partition_result=audit_membership_partition(root)
            req(partition_result.get("status")=="PASS","annotation membership partition audit failed: "+"; ".join(partition_result.get("errors",[])[:5]))
            support_registry=root/"state/annotation_support_registry.tsv"
            current_ledger=root/"state/cell_ledger.tsv.gz"
            if not current_ledger.exists():current_ledger=root/"state/cell_ledger.tsv"
            support_result=validate_support_registry(root,support_registry,current_ledger)
            req(support_result.get("status")=="PASS","per-label annotation support validation failed: "+"; ".join(support_result.get("errors",[])[:5]))
            for name in ["state/recluster_cohort_registry.tsv", "state/direct_return_registry.tsv", "state/route_attempt_registry.tsv", "state/workflow_event_registry.tsv", "state/annotation_view_registry.tsv", "state/annotation_support_registry.tsv"]:
                req((root/name).exists(), f"missing annotation workflow registry: {name}")
            route_rows = read_tsv(root/"state/route_attempt_registry.tsv") if (root/"state/route_attempt_registry.tsv").exists() else []
            event_rows = read_tsv(root/"state/workflow_event_registry.tsv") if (root/"state/workflow_event_registry.tsv").exists() else []
            view_rows = read_tsv(root/"state/annotation_view_registry.tsv") if (root/"state/annotation_view_registry.tsv").exists() else []
            req(bool(route_rows), "annotation release has no terminal assisted-route record, including zero-query Atlas audit")
            req(bool(event_rows), "annotation release has no workflow events")
            req("final" in {x.get("view") for x in view_rows if x.get("status")=="validated"}, "single final annotation is not validated")
            req(all(x.get("analysis_view")=="final" for x in assets), "mandatory marker dotplots must use the final annotation")
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
            for level in ["broad"] + (["subtype"] if subtype_required else []):
                stem=f"final_{level}_UMAP";req((root/"figures"/f"{stem}.png").exists() and (root/"figures"/f"{stem}.pdf").exists(),f"missing final annotation overview: {stem}")
                if config.get("modality")=="spatial":
                    stem=f"final_{level}_spatial";req((root/"figures"/f"{stem}.png").exists() and (root/"figures"/f"{stem}.pdf").exists(),f"missing final annotation spatial overview: {stem}")
            for level in ["broad"] + (["subtype"] if subtype_required else []):
                named=any_file(f"tables/final_{level}_DEG_one_vs_rest_all.tsv*")
                generic=list(root.glob(f"tables/{level}_DEG_one_vs_rest_all.tsv*"))
                generic_final=False
                for candidate in generic:
                    if not candidate.is_file() or candidate.stat().st_size < 100: continue
                    try:
                        rows=read_tsv(candidate)
                        generic_final=bool(rows) and all(row.get("analysis_view")=="final" for row in rows[:min(100,len(rows))])
                    except (OSError, UnicodeDecodeError):
                        generic_final=False
                    if generic_final: break
                req(named or generic_final,f"missing final {level} DEG with explicit analysis_view provenance")
            expected_labels = {}
            for level in ["broad", "subtype"]:
                candidates = [root/f"tables/final_{level}_DEG_one_vs_rest_all.tsv", root/f"tables/{level}_DEG_one_vs_rest_all.tsv"]
                deg_path = next((path for path in candidates if path.is_file()), None)
                if deg_path:
                    deg_rows = read_tsv(deg_path)
                    expected_labels[level] = {row.get("label", "") for row in deg_rows if row.get("label", "") and row.get("analysis_view") == "final"}
            for asset in assets:
                level = asset.get("level", "")
                expected = expected_labels.get(level, set())
                source = Path(asset.get("source", "")); source = source if source.is_absolute() else root/source
                if not expected or not source.is_file():
                    continue
                source_rows = read_tsv(source)
                observed_labels = {row.get("label", "") for row in source_rows if row.get("label", "")}
                marker_groups = {row.get("marker_group", "") for row in source_rows if row.get("marker_group", "")}
                req(observed_labels == expected, f"{level}/{asset.get('panel')} dotplot does not cover every final DEG label")
                req(expected.issubset(marker_groups), f"{level}/{asset.get('panel')} marker groups do not cover every current label")
                try:
                    req(int(asset.get("n_labels", 0)) == len(expected), f"{level}/{asset.get('panel')} n_labels disagrees with final DEG")
                except (TypeError, ValueError):
                    req(False, f"{level}/{asset.get('panel')} n_labels is invalid")
        required_combos = {("broad", "canonical"), ("broad", "data_specific")}
        if subtype_required:
            required_combos.update({("subtype", "canonical"), ("subtype", "data_specific")})
        req(required_combos.issubset(combos), "full release requires canonical and data-specific dotplots for broad and every released high-confidence subtype level")
        if not workflow_required:
            req(any_file("tables/broad_DEG_one_vs_rest_all.tsv*"), "missing broad one-vs-rest DEG")
            req(any_file("tables/subtype_DEG_one_vs_rest_all.tsv*"), "missing subtype one-vs-rest DEG")
        req(any_file("state/cell_ledger.tsv*", 1000) or any_file("tables/cell_ledger.tsv*", 1000) or any_file("tables/cell_metadata*.tsv*", 1000), "missing cell-level ledger")
        for name in ["state/recluster_cohort_registry.tsv","state/direct_return_registry.tsv","state/run_registry.tsv","state/next_action_queue.tsv","state/master_quality_approval.json","state/final_annotation_confirmation.json","state/annotation_support_registry.tsv","review/confirmation/index.html","provenance/confirmation_review_manifest.json","provenance/full_feature_validation.json","provenance/release_manifest.tsv","provenance/checksums.sha256","provenance/release_sessionInfo.txt"]:
            req((root/name).exists(),f"missing full-release provenance: {name}")
        feature_validation = root / "provenance/full_feature_validation.json"
        if feature_validation.exists():
            req(json.loads(feature_validation.read_text()).get("status") == "PASS", "full-feature validation did not pass")
        session = root / "provenance/release_sessionInfo.txt"
        if session.exists():
            req(session.stat().st_size > 100, "release session info is empty or undersized")
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
            req("state/annotation_support_registry.tsv" in checksum_set, "checksums do not cover annotation support reasons")
            req("review/confirmation/index.html" in checksum_set, "checksums do not cover the pre-confirmation review")
            req(any(x.startswith("state/cell_ledger.tsv") for x in checksum_set), "checksums do not cover the cell ledger")
            broad_deg_prefixes=["tables/final_broad_DEG_one_vs_rest_all.tsv","tables/broad_DEG_one_vs_rest_all.tsv"] if workflow_required else ["tables/broad_DEG_one_vs_rest_all.tsv"]
            subtype_deg_prefixes=["tables/final_subtype_DEG_one_vs_rest_all.tsv","tables/subtype_DEG_one_vs_rest_all.tsv"] if workflow_required else ["tables/subtype_DEG_one_vs_rest_all.tsv"]
            req(any(any(x.startswith(prefix) for prefix in broad_deg_prefixes) for x in checksum_set), "checksums do not cover broad DEG")
            if subtype_required:
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
        for stem in ["final_broad_UMAP"] + (["final_subtype_UMAP"] if subtype_required else []):
            req((root/"figures"/f"{stem}.png").exists() and (root/"figures"/f"{stem}.pdf").exists(), f"missing overview pair: {stem}")
        if config.get("modality") == "spatial":
            for stem in ["final_broad_spatial"] + (["final_subtype_spatial"] if subtype_required else []):
                req((root/"figures"/f"{stem}.png").exists() and (root/"figures"/f"{stem}.pdf").exists(), f"missing spatial overview pair: {stem}")
            req((root/"tables/spatial_node_asset_index.tsv").exists(), "missing spatial node index")
            req((root/"tables/spatial_gene_asset_index.tsv").exists(), "missing spatial gene index")
    if report.exists():
        report_text=report.read_text(encoding="utf-8")
        if args.profile=="full" and workflow_required:
            for section in ["id='final'","id='routes'","id='workflow'"]:req(section in report_text,f"report lacks required annotation-workflow section {section}")
        parser = Links(); parser.feed(report.read_text(encoding="utf-8"))
        for value in parser.paths:
            req((report.parent / value).resolve().exists(), f"broken report link: {value}")
    result = {"status": "PASS" if not errors else "FAIL", "profile": args.profile, "errors": sorted(set(errors)), "dotplot_assets": len(assets)}
    out = root / "provenance/release_audit.json"; out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    dependencies = [path for path in [report, root/"provenance/completion_gate.json", root/"provenance/release_manifest.tsv", root/"provenance/checksums.sha256"] if path.is_file()]
    build_dependency_manifest(out, dependencies, {"asset_class": "release_audit"})
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if not errors else 1


if __name__ == "__main__": raise SystemExit(main())
