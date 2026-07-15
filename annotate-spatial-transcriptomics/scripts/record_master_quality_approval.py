#!/usr/bin/env python3
"""Record a main-Agent verdict for the fully completed annotation snapshot."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from master_quality_lib import BOUND_ARTIFACTS, CHECKLIST_IDS, sha256, validate_bound_record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--assessment", required=True, type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    request_path = root / "provenance/master_quality_review_request.json"
    if not request_path.is_file():
        raise SystemExit("request_master_quality_review.py must run after the completion gate")
    request = json.loads(request_path.read_text(encoding="utf-8"))
    if request.get("status") != "AWAITING_MASTER_QUALITY_APPROVAL":
        raise SystemExit("master quality review request has an invalid status")
    binding_errors = validate_bound_record(root, request)
    if binding_errors:
        raise SystemExit("master quality review request is stale: " + "; ".join(binding_errors))
    assessment = json.loads(args.assessment.read_text(encoding="utf-8"))
    verdict = str(assessment.get("verdict", "")).upper()
    if assessment.get("reviewer_role") != "main_conversation_agent" or not str(assessment.get("reviewer_id", "")).strip():
        raise SystemExit("only an identified main conversation Agent may record this gate")
    if assessment.get("reviewed_after_all_annotation_complete") is not True:
        raise SystemExit("assessment must explicitly verify that all annotation routes and final writeback were complete")
    if verdict not in {"PASS", "RETURN_FOR_ITERATION"}:
        raise SystemExit("verdict must be PASS or RETURN_FOR_ITERATION")
    for field in ("comparison_to_reference", "rationale"):
        if len(str(assessment.get(field, "")).strip()) < 20:
            raise SystemExit(f"master assessment requires a substantive {field}")
    checklist = assessment.get("checklist", {})
    blocks: list[str] = []
    for item in CHECKLIST_IDS:
        entry = checklist.get(item, {}) if isinstance(checklist, dict) else {}
        if entry.get("status") not in {"PASS", "CONCERN", "BLOCK"}:
            raise SystemExit(f"master assessment checklist is incomplete: {item}")
        if entry.get("status") == "BLOCK":
            blocks.append(item)
    if verdict == "PASS" and blocks:
        raise SystemExit("PASS is forbidden while quality dimensions are blocked: " + ", ".join(blocks))
    if verdict == "RETURN_FOR_ITERATION" and not assessment.get("requested_revisions"):
        raise SystemExit("RETURN_FOR_ITERATION requires requested_revisions")

    status = "PASS" if verdict == "PASS" else "RETURN_FOR_ITERATION"
    approval = {
        "status": status,
        "reviewed_at": datetime.now(timezone.utc).isoformat(),
        "reviewer_role": assessment["reviewer_role"],
        "reviewer_id": assessment["reviewer_id"],
        "reviewed_after_all_annotation_complete": True,
        "comparison_semantics": request["comparison_semantics"],
        "comparison_to_reference": assessment["comparison_to_reference"],
        "rationale": assessment["rationale"],
        "biological_concerns": assessment.get("biological_concerns", []),
        "requested_revisions": assessment.get("requested_revisions", []),
        "checklist": checklist,
        "request_path": "provenance/master_quality_review_request.json",
        "request_sha256": sha256(request_path),
    }
    for key in BOUND_ARTIFACTS:
        approval[key] = request[key]
        approval[f"{key}_sha256"] = request[f"{key}_sha256"]
    if request.get("strategy_preset_record"):
        approval["strategy_preset_record"] = request["strategy_preset_record"]
        approval["strategy_preset_record_sha256"] = request["strategy_preset_record_sha256"]
        approval["open_world_lineage_audit"] = request["open_world_lineage_audit"]
        approval["open_world_lineage_audit_sha256"] = request["open_world_lineage_audit_sha256"]
        approval["open_world_lineage_audit_source"] = request["open_world_lineage_audit_source"]
        approval["open_world_lineage_audit_source_sha256"] = request["open_world_lineage_audit_source_sha256"]
        approval["candidate_lineage_catalog"] = request["candidate_lineage_catalog"]
        approval["candidate_lineage_catalog_sha256"] = request["candidate_lineage_catalog_sha256"]
        approval["open_world_biological_profile"] = request["open_world_biological_profile"]
        approval["open_world_biological_profile_sha256"] = request["open_world_biological_profile_sha256"]
    out = root / "state/master_quality_approval.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(approval, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(approval, ensure_ascii=False, indent=2))
    return 0 if status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
