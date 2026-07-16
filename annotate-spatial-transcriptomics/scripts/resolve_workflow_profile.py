#!/usr/bin/env python3
"""Resolve a biological workflow profile without assigning biological labels."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path


def load(path: Path | None) -> dict:
    return json.loads(path.read_text(encoding="utf-8")) if path else {}


def normalized(value: object) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip().lower())


def alias_hit(value: object, aliases: list[str]) -> bool:
    text = normalized(value)
    return any(normalized(alias) in text for alias in aliases)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--context", required=True, type=Path)
    parser.add_argument("--discovery", type=Path)
    parser.add_argument("--r-inspection", type=Path)
    parser.add_argument("--input-path")
    parser.add_argument("--conversion-source", default="")
    parser.add_argument("--conversion-manifest", type=Path)
    parser.add_argument("--layer-audit", type=Path)
    parser.add_argument("--coordinate-audit", type=Path)
    parser.add_argument("--marker-coverage-audit", type=Path)
    parser.add_argument("--profile", type=Path)
    parser.add_argument("--strategy-preset", choices=["none", "sheep_ovary_same_batch_rfirst"], default="none")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--workflow-record-out", type=Path)
    parser.add_argument("--input-contract-out", type=Path)
    args = parser.parse_args()

    skill_root = Path(__file__).resolve().parents[1]
    profile_path = args.profile or skill_root / "references/profiles/sheep_ovary_rfirst_profile.json"
    profile = load(profile_path)
    context = load(args.context)
    discovery = load(args.discovery)
    inspection = load(args.r_inspection)

    species = context.get("species", "")
    tissue = context.get("tissue", "")
    sheep_ovary = alias_hit(species, profile["species_aliases"]) and alias_hit(
        tissue, profile["tissue_aliases"]
    )
    inspected_seurat = normalized(inspection.get("type")) == "seurat"
    raw_assays = inspection.get("assays", [])
    if isinstance(raw_assays, str):
        raw_assays = [raw_assays]
    elif not isinstance(raw_assays, list):
        raw_assays = list(raw_assays or [])
    assays = {normalized(item) for item in raw_assays}
    full_feature = int(inspection.get("n_features", 0) or 0) >= 10000
    spatial_assay = "spatial" in assays

    paths: list[str] = []
    if args.input_path:
        paths.append(args.input_path)
    for row in discovery.get("r_objects", []):
        if isinstance(row, dict) and row.get("path"):
            paths.append(str(row["path"]))
    source_text = " ".join(paths + [args.conversion_source])
    stereopy_cellbin_hint = bool(re.search(r"stereopy|cellbin[_-]?pped", source_text, re.I))
    conversion = load(args.conversion_manifest)
    layer_audit = load(args.layer_audit)
    coordinate_audit = load(args.coordinate_audit)
    marker_coverage = load(args.marker_coverage_audit)
    checks = {
        "object_inspection_passed": inspected_seurat and spatial_assay and full_feature,
        "counts_data_layer_audit_passed": layer_audit.get("status") == "PASS" and layer_audit.get("raw_counts_verified") is True,
        "profile_marker_coverage_passed": marker_coverage.get("status") == "PASS" and marker_coverage.get("coverage_sufficient") is True,
        "cell_coordinate_consistency_passed": coordinate_audit.get("status") == "PASS" and coordinate_audit.get("cell_coordinate_ids_match") is True,
        "conversion_provenance_verified": conversion.get("status") == "PASS" and conversion.get("source_platform") == "StereoPy" and conversion.get("cellbin_pped_conversion") is True and conversion.get("same_batch") is True,
    }
    fixed_contract = sheep_ovary and args.strategy_preset == "sheep_ovary_same_batch_rfirst" and all(checks.values())

    preset = None
    preset_bindings = None
    if args.strategy_preset != "none":
        if not sheep_ovary:
            raise SystemExit("the sheep-ovary same-batch R-first preset requires verified sheep ovary context")
        if not (inspected_seurat and full_feature):
            raise SystemExit("the sheep-ovary same-batch R-first preset requires a full-feature Seurat RDS inspection")
        preset_path = skill_root / "references/profiles/sheep_ovary_same_batch_rfirst_preset.json"
        preset = load(preset_path)
        reference_root = preset_path.parent
        preset_bindings = {"preset": str(preset_path.resolve()), "preset_sha256": sha256(preset_path)}
        for key, filename in preset["bound_references"].items():
            path = reference_root / filename
            if not path.is_file():
                raise SystemExit(f"strategy preset binding is missing: {key}")
            preset_bindings[key] = str(path.resolve())
            preset_bindings[f"{key}_sha256"] = sha256(path)

    if sheep_ovary and inspected_seurat and full_feature:
        backbone = "seurat_r_first"
        status = "RESOLVED"
    elif sheep_ovary:
        backbone = "inspect_inputs_then_prefer_r_first_if_full_feature_rds_exists"
        status = "NEEDS_INPUT_INSPECTION"
    else:
        backbone = "generic_context_adaptive"
        status = "GENERIC_PROFILE"

    result = {
        "status": status,
        "profile_id": profile["profile_id"] if sheep_ovary else "generic",
        "sheep_ovary_context": sheep_ovary,
        "preferred_backbone": backbone,
        "full_feature_seurat_detected": inspected_seurat and full_feature,
        "stereopy_cellbin_path_or_feature_hint": stereopy_cellbin_hint,
        "stereopy_cellbin_pped_detected": fixed_contract,
        "fixed_cellbin_preprocessing_required": fixed_contract,
        "fixed_cellbin_profile": profile["stereopy_cellbin_pped_contract"] if fixed_contract else None,
        "strategy_preset_status": "ACTIVE" if preset else "NOT_REQUESTED",
        "strategy_preset_id": preset["preset_id"] if preset else None,
        "strategy_preset": preset,
        "strategy_preset_bindings": preset_bindings,
        "strategy_preset_preprocessing_mode": ("fixed_verified_same_batch_contract" if fixed_contract else "strategy_only_pending_cellbin_provenance") if preset else None,
        "primary_public_atlas": profile["external_reference_policy"]["primary_public_atlas"] if sheep_ovary else None,
        "atlas_calibration_origin_required": profile["external_reference_policy"]["calibration_origin"] if sheep_ovary else None,
        "matched_reference_count_object_changes_external_priority": True,
        "dotplot_only_reference_allows_cell_transfer": False,
        "adaptive_fields": ["selected broad resolution", "cohort membership", "cohort PCs/k/resolution", "biological labels", "optional subtypes"],
        "evidence_only": True,
        "warning": "Profile/preset resolution selects a workflow and safety contract; it never assigns cells or imports example labels, resolutions, memberships or proportions."
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    workflow_record = {
        "status": "ACTIVE" if fixed_contract else "PENDING_PROVENANCE",
        "workflow_profile_id": "seurat_stereopy_cellbin_same_batch_rfirst" if fixed_contract else profile.get("profile_id"),
        "workflow_profile": str(profile_path.resolve()),
        "workflow_profile_sha256": sha256(profile_path),
        "candidate_resolution_grid": profile["stereopy_cellbin_pped_contract"]["clustering"]["candidate_resolutions"] if fixed_contract else [],
    }
    input_contract = {
        "status": "PASS" if fixed_contract else "BLOCKED",
        "verified_stereopy_cellbin_same_batch": fixed_contract,
        "checks": checks,
        "path_or_feature_hint_only": stereopy_cellbin_hint,
        "source_artifacts": {
            "r_inspection": str(args.r_inspection.resolve()) if args.r_inspection else "",
            "conversion_manifest": str(args.conversion_manifest.resolve()) if args.conversion_manifest else "",
            "layer_audit": str(args.layer_audit.resolve()) if args.layer_audit else "",
            "coordinate_audit": str(args.coordinate_audit.resolve()) if args.coordinate_audit else "",
            "marker_coverage_audit": str(args.marker_coverage_audit.resolve()) if args.marker_coverage_audit else ""
        }
    }
    if args.workflow_record_out:
        args.workflow_record_out.parent.mkdir(parents=True, exist_ok=True)
        args.workflow_record_out.write_text(json.dumps(workflow_record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if args.input_contract_out:
        args.input_contract_out.parent.mkdir(parents=True, exist_ok=True)
        args.input_contract_out.write_text(json.dumps(input_contract, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
