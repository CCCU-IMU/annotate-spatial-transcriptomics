#!/usr/bin/env python3
"""Validate the immutable GSE233801 sheep-ovary Atlas descriptor and assets."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from evidence_schema_lib import sha256


ACTIVE_BUNDLE_ID = "sheep_ovary_GSE233801_split_wall_v2"
LEGACY_BUNDLE_ID = "sheep_ovary_GSE233801_v1"
BUNDLE_ID = ACTIVE_BUNDLE_ID
REFERENCE_IDS = {
    ACTIVE_BUNDLE_ID: "GSE233801_independent_R_res0p4_split_wall_v002",
    LEGACY_BUNDLE_ID: "GSE233801_independent_R_reviewed_broad_v001",
}
CAPABILITIES = {"supported", "challenge_only", "not_evaluable", "unsupported"}


def validate(document: dict, asset_root: Path | None = None) -> list[str]:
    errors: list[str] = []
    bundle_id = str(document.get("bundle_id", ""))
    expected_schema = "2.0" if bundle_id == ACTIVE_BUNDLE_ID else "1.0"
    if document.get("schema_version") != expected_schema:
        errors.append("fixed Atlas schema version differs")
    if bundle_id not in REFERENCE_IDS:
        errors.append("sheep ovary Atlas bundle ID is not an approved GSE233801 bundle")
    if document.get("reference_id") != REFERENCE_IDS.get(bundle_id):
        errors.append("fixed Atlas reference ID differs")
    if document.get("immutable") is not True or document.get("reusable_across_queries") is not True:
        errors.append("fixed Atlas is not declared immutable and reusable")
    if document.get("query_reference_joint_retraining") is not False:
        errors.append("fixed Atlas permits query-reference joint retraining")
    source = document.get("source_label_capability", {})
    release = document.get("release_broad_capability", {})
    if not source or not release:
        errors.append("fixed Atlas capability matrix is empty")
    if any(value not in CAPABILITIES for value in [*source.values(), *release.values()]):
        errors.append("fixed Atlas capability matrix contains an invalid state")
    if source.get("Vascular/endothelial") == "supported":
        errors.append("legacy mixed vascular label cannot have release capability")
    if bundle_id == ACTIVE_BUNDLE_ID:
        for broad in (
            "Granulosa", "Stromal/mesenchymal", "Endothelial",
            "Pericyte/mural", "Smooth muscle", "Immune",
        ):
            if release.get(broad) != "supported":
                errors.append(f"split-wall GSE233801 lacks supported capability for {broad}")
        for broad in ("Epithelial/mesothelial", "Theca"):
            if release.get(broad) != "challenge_only":
                errors.append(f"split-wall GSE233801 must keep {broad} challenge-only")
        for broad in ("Oocyte", "Luteal"):
            if release.get(broad) != "unsupported":
                errors.append(f"split-wall GSE233801 cannot adjudicate {broad}")
        if (
            document.get("new_routing_authority") is not True
            or document.get("legacy_resume_only") is not False
        ):
            errors.append("active split-wall Atlas routing authority is invalid")
    elif bundle_id == LEGACY_BUNDLE_ID:
        for broad in (
            "Endothelial", "Pericyte/mural", "Smooth muscle",
            "Epithelial/mesothelial", "Theca", "Oocyte", "Luteal",
        ):
            if release.get(broad) == "supported":
                errors.append(f"legacy GSE233801 cannot formally adjudicate {broad}")
        if (
            document.get("new_routing_authority") is not False
            or document.get("legacy_resume_only") is not True
        ):
            errors.append("legacy Atlas must be resume-only")
    runtime = document.get("runtime_policy", {})
    if (
        runtime.get("writeback_scope") != "unlabeled_after_broad_merge_only"
        or runtime.get("defined_label_conflict_action")
        != "biological_review_only_no_overwrite"
        or runtime.get("fine_label_authority") is not False
    ):
        errors.append("fixed Atlas runtime authority is too broad")
    if asset_root is not None:
        for name, expected in document.get("asset_hashes", {}).items():
            path = asset_root / name
            if not path.is_file() or sha256(path) != expected:
                errors.append(f"fixed Atlas asset is missing or stale: {name}")
    return errors


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle-manifest", required=True, type=Path)
    ap.add_argument("--asset-root", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    document = json.loads(args.bundle_manifest.read_text(encoding="utf-8"))
    errors = validate(document, args.asset_root)
    result = {
        "status": "PASS" if not errors else "BLOCKED",
        "bundle_id": document.get("bundle_id", ""),
        "reference_id": document.get("reference_id", ""),
        "bundle_manifest": {
            "path": str(args.bundle_manifest.resolve()),
            "sha256": sha256(args.bundle_manifest),
        },
        "asset_root": str(args.asset_root.resolve()) if args.asset_root else "",
        "errors": errors,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
