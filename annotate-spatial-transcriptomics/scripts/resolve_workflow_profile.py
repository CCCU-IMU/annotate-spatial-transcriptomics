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
    parser.add_argument("--profile", type=Path)
    parser.add_argument("--strategy-preset", choices=["none", "sheep_ovary_same_batch_rfirst"], default="none")
    parser.add_argument("--out", required=True, type=Path)
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
    stereopy_cellbin = bool(re.search(r"stereopy|cellbin[_-]?pped", source_text, re.I))
    fixed_contract = sheep_ovary and inspected_seurat and full_feature and spatial_assay and stereopy_cellbin

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
        "stereopy_cellbin_pped_detected": stereopy_cellbin,
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
        "adaptive_fields": ["selected broad resolution", "pool membership", "pool PCs/k/resolution", "biological labels", "optional subtypes"],
        "evidence_only": True,
        "warning": "Profile/preset resolution selects a workflow and safety contract; it never assigns cells or imports example labels, resolutions, memberships or proportions."
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
