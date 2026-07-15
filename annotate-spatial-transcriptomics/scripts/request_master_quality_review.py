#!/usr/bin/env python3
"""Freeze the fully annotated sample for independent main-Agent quality review."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

from master_quality_lib import BOUND_ARTIFACTS, CHECKLIST_IDS, sha256


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def rel_or_abs(root: Path, path: Path) -> str:
    try:
        return str(path.resolve().relative_to(root))
    except ValueError:
        return str(path.resolve())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--quality-reference", type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    default_reference = Path(__file__).resolve().parents[1] / "references/profiles/sheep_ovary_rfirst_case_reference.md"
    paths = {
        "cell_ledger": root / "state/cell_ledger.tsv.gz",
        "cluster_ledger": root / "state/cluster_decision_ledger.tsv",
        "completion_gate": root / "provenance/completion_gate.json",
        "release_taxonomy_audit": root / "provenance/release_taxonomy_audit.json",
        "multiroute_audit": root / "provenance/multiroute_audit.json",
        "state_validation": root / "provenance/state_validation.json",
        "iteration_plan": root / "provenance/iteration_plan.json",
        "next_action_queue": root / "state/next_action_queue.tsv",
        "pool_registry": root / "state/pool_registry.tsv",
        "route_attempt_registry": root / "state/route_attempt_registry.tsv",
        "run_registry": root / "state/run_registry.tsv",
        "final_annotation_census": root / "tables/final_annotation_census.tsv",
        "annotation_view_registry": root / "state/annotation_view_registry.tsv",
        "annotation_support_registry": root / "state/annotation_support_registry.tsv",
        "review_asset_manifest": root / "review/confirmation/assets/review_asset_manifest.json",
        "quality_reference": (args.quality_reference or default_reference).resolve(),
    }
    missing = [key for key in BOUND_ARTIFACTS if not paths[key].is_file()]
    if missing:
        raise SystemExit("master review cannot start; required final artifacts are missing: " + ", ".join(missing))

    completion = json.loads(paths["completion_gate"].read_text(encoding="utf-8"))
    taxonomy = json.loads(paths["release_taxonomy_audit"].read_text(encoding="utf-8"))
    multiroute = json.loads(paths["multiroute_audit"].read_text(encoding="utf-8"))
    state_validation = json.loads(paths["state_validation"].read_text(encoding="utf-8"))
    asset_manifest = json.loads(paths["review_asset_manifest"].read_text(encoding="utf-8"))
    errors: list[str] = []
    if completion.get("status") != "PASS":
        errors.append("completion gate is not PASS")
    if taxonomy.get("pass") is not True or taxonomy.get("metadata_sha256") != sha256(paths["cell_ledger"]):
        errors.append("release taxonomy audit is not PASS/current")
    if multiroute.get("status") != "PASS":
        errors.append("Route A/B/C/D completion audit is not PASS")
    if state_validation.get("status") != "PASS":
        errors.append("state validation is not PASS")
    # The completion gate already owns detailed route/pool/run closure.  The
    # master gate must not make the reviewer repeat that mechanical audit.
    views = read_tsv(paths["annotation_view_registry"])
    if not any(row.get("view") == "final" and row.get("status") == "validated" for row in views):
        errors.append("the single final annotation view is not validated")
    if not read_tsv(paths["final_annotation_census"]):
        errors.append("final annotation census is empty")
    if not any(row.get("status") == "validated" for row in read_tsv(paths["annotation_support_registry"])):
        errors.append("validated annotation support reasons are absent")
    if asset_manifest.get("status") != "PASS" or asset_manifest.get("label_column") != "primary_broad_label":
        errors.append("lightweight broad spatial/dotplot evidence assets are not PASS")
    for key in ("spatial_png", "dotplot_png", "dotplot_source", "palette_tsv"):
        value = str(asset_manifest.get(key, ""))
        artifact = Path(value)
        if value and not artifact.is_absolute():
            artifact = root / artifact
        if not value or not artifact.is_file() or sha256(artifact) != asset_manifest.get(f"{key}_sha256"):
            errors.append(f"lightweight quality-review evidence is missing or stale: {key}")
    if errors:
        raise SystemExit("master review is forbidden before all annotation is complete: " + "; ".join(errors))

    project = json.loads((root / "config/project.json").read_text(encoding="utf-8"))
    request = {
        "status": "AWAITING_MASTER_QUALITY_APPROVAL",
        "sample_id": project.get("sample_id"),
        "decision_version": project.get("decision_version"),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "stage_contract": "post_pool_reclustering_post_multiroute_rescue_post_final_writeback_post_completion_gate",
        "reviewer_role_required": "main_conversation_agent",
        "comparison_semantics": "Compare biological reasonableness and evidence quality with the validated sheep-ovary R-first reference; exact labels, counts and subtype depth need not match.",
        "required_verdict": "PASS (including PASS with documented concerns) when annotation quality is comparable to the validated reference; otherwise RETURN_FOR_ITERATION.",
        "checklist_ids": list(CHECKLIST_IDS),
    }
    for key in BOUND_ARTIFACTS:
        request[key] = rel_or_abs(root, paths[key])
        request[f"{key}_sha256"] = sha256(paths[key])
    if project.get("strategy_preset_requested"):
        active_preset = root / "config/active_strategy_preset.json"
        open_world_audit = root / "provenance/open_world_lineage_audit_validation.json"
        if not active_preset.is_file():
            raise SystemExit("requested strategy preset record is missing at the master quality gate")
        if not open_world_audit.is_file():
            raise SystemExit("requested strategy preset lacks the open-world lineage audit at the master quality gate")
        open_world_validation = json.loads(open_world_audit.read_text(encoding="utf-8"))
        if open_world_validation.get("status") != "PASS":
            raise SystemExit("requested strategy preset has a failed open-world lineage audit")
        open_world_source = Path(str(open_world_validation.get("audit", "")))
        candidate_catalog = Path(str(open_world_validation.get("candidate_catalog", "")))
        open_world_biological_profile = Path(str(open_world_validation.get("biological_profile", "")))
        if not open_world_source.is_file() or sha256(open_world_source) != open_world_validation.get("audit_sha256"):
            raise SystemExit("open-world lineage audit source is missing or stale")
        if not candidate_catalog.is_file() or sha256(candidate_catalog) != open_world_validation.get("candidate_catalog_sha256"):
            raise SystemExit("open-world candidate lineage catalog is missing or stale")
        if not open_world_biological_profile.is_file() or sha256(open_world_biological_profile) != open_world_validation.get("biological_profile_sha256"):
            raise SystemExit("open-world biological profile is missing or stale")
        request["strategy_preset_record"] = rel_or_abs(root, active_preset)
        request["strategy_preset_record_sha256"] = sha256(active_preset)
        request["open_world_lineage_audit"] = rel_or_abs(root, open_world_audit)
        request["open_world_lineage_audit_sha256"] = sha256(open_world_audit)
        request["open_world_lineage_audit_source"] = rel_or_abs(root, open_world_source)
        request["open_world_lineage_audit_source_sha256"] = sha256(open_world_source)
        request["candidate_lineage_catalog"] = rel_or_abs(root, candidate_catalog)
        request["candidate_lineage_catalog_sha256"] = sha256(candidate_catalog)
        request["open_world_biological_profile"] = rel_or_abs(root, open_world_biological_profile)
        request["open_world_biological_profile_sha256"] = sha256(open_world_biological_profile)
    out = root / "provenance/master_quality_review_request.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(request, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    template = {
        "reviewer_role": "main_conversation_agent",
        "reviewer_id": "",
        "reviewed_after_all_annotation_complete": True,
        "verdict": "PASS_or_RETURN_FOR_ITERATION",
        "comparison_to_reference": "",
        "rationale": "",
        "biological_concerns": [],
        "requested_revisions": [],
        "checklist": {item: {"status": "PASS_CONCERN_or_BLOCK", "note": ""} for item in CHECKLIST_IDS},
    }
    template_path = root / "provenance/master_quality_assessment.template.json"
    template_path.write_text(json.dumps(template, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({**request, "assessment_template": str(template_path)}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
