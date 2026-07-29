#!/usr/bin/env python3
"""Ensure every default broad candidate can satisfy the broad-family gate."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from lineage_controller_lib import catalog_candidates


def resolve(payload: dict, dotted: str):
    value = payload
    for key in dotted.split("."):
        value = value[key]
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    profile = json.loads(args.profile.read_text(encoding="utf-8"))
    catalog = json.loads(args.catalog.read_text(encoding="utf-8"))
    errors: list[str] = []
    rows: list[dict] = []
    candidates = list(catalog_candidates(catalog).values())
    for candidate in candidates:
        if candidate.get("review_required") is not True:
            continue
        candidate_id = candidate.get("candidate_id", "")
        required_fields = {
            "candidate_role", "release_broad_label", "release_fine_label",
            "parent_broad_label", "writeback_strategy", "specificity_priority",
            "hard_anti_families", "soft_anti_families", "context_requirements",
        }
        missing_fields = sorted(required_fields - set(candidate))
        if missing_fields:
            errors.append(
                f"{candidate_id}: v2.2 candidate taxonomy lacks {', '.join(missing_fields)}"
            )
        role = candidate.get("candidate_role")
        if role not in {"broad", "fine", "state", "exploratory"}:
            errors.append(f"{candidate_id}: invalid candidate_role {role!r}")
        if role in {"state", "exploratory"} and (
            candidate.get("release_broad_label") or candidate.get("release_fine_label")
        ):
            errors.append(f"{candidate_id}: non-release role carries a release label")
        if candidate_id == "pericyte_mural" and (
            role != "broad"
            or candidate.get("release_broad_label") != "Pericyte/mural"
            or candidate.get("release_fine_label")
            or candidate.get("parent_broad_label")
        ):
            errors.append("pericyte_mural: must be an independent broad lineage")
        if candidate_id == "vascular_endothelial" and (
            role != "broad"
            or candidate.get("release_broad_label") != "Endothelial"
            or candidate.get("release_fine_label")
            or candidate.get("parent_broad_label")
        ):
            errors.append("vascular_endothelial: must release independent broad Endothelial")
        if candidate_id == "lymphatic_endothelial" and (
            role != "fine"
            or candidate.get("release_broad_label") != "Endothelial"
            or candidate.get("parent_broad_label") != "Endothelial"
            or candidate.get("release_fine_label") != "Lymphatic endothelial"
            or candidate.get("parent_broad_reconstruction_allowed") is not True
        ):
            errors.append("lymphatic_endothelial: must be parent-locked to Endothelial")
        forbidden = set(catalog.get("taxonomy_policy", {}).get(
            "forbidden_runtime_release_labels", []
        ))
        if "Vascular-associated" not in forbidden:
            errors.append("catalog does not forbid legacy Vascular-associated release")
        if candidate.get("release_broad_label") == "Vascular-associated":
            errors.append(f"{candidate_id}: releases forbidden Vascular-associated")
        unit_anti = candidate.get("hard_anti_families_by_observation_unit", {})
        if unit_anti:
            if set(unit_anti) != {"cell", "nucleus", "cellbin", "spot"}:
                errors.append(f"{candidate_id}: incomplete observation-unit anti policy")
            if not all(isinstance(value, list) for value in unit_anti.values()):
                errors.append(f"{candidate_id}: invalid observation-unit anti policy")
        path = candidate.get("profile_program", "")
        try:
            program = resolve(profile, path)
        except (KeyError, TypeError):
            errors.append(f"{candidate_id}: profile program cannot be resolved: {path}")
            continue
        families = program.get("positive_families", {}) if isinstance(program, dict) else {}
        if not families and candidate.get("candidate_role") == "state":
            families = candidate.get("state_positive_families", {})
        if not families and candidate.get("candidate_role") == "fine":
            parent = candidate.get("parent_broad_label")
            parent_candidates = [
                row for row in candidates
                if row.get("candidate_role") == "broad"
                and row.get("release_broad_label") == parent
            ]
            if not parent_candidates:
                errors.append(f"{candidate_id}: fine candidate lacks a resolvable broad parent")
                continue
            try:
                parent_program = resolve(profile, parent_candidates[0]["profile_program"])
            except (KeyError, TypeError):
                errors.append(f"{candidate_id}: parent profile program cannot be resolved")
                continue
            subtype_genes = {
                str(gene).strip() for gene in (program if isinstance(program, list) else [])
                if str(gene).strip()
            }
            parent_genes = {
                str(gene).strip()
                for genes in parent_program.get("positive_families", {}).values()
                for gene in genes
                if str(gene).strip()
            } - subtype_genes
            families = {
                "parent_identity": sorted(parent_genes),
                "fine_discriminator": sorted(subtype_genes),
            }
        valid = {
            name: sorted({str(gene).strip() for gene in genes if str(gene).strip()})
            for name, genes in families.items()
            if isinstance(genes, list)
        }
        valid = {name: genes for name, genes in valid.items() if genes}
        overlaps = []
        names = sorted(valid)
        for i, left in enumerate(names):
            for right in names[i + 1 :]:
                shared = sorted(set(valid[left]).intersection(valid[right]))
                if shared:
                    overlaps.append({"left": left, "right": right, "shared": shared})
        if len(valid) < 2:
            errors.append(
                f"{candidate_id}: open-world broad scan requires >=2 explicit positive_families, found {len(valid)}"
            )
        if overlaps:
            errors.append(f"{candidate_id}: positive families overlap and are not independent: {overlaps}")
        rows.append(
            {
                "candidate_id": candidate_id,
                "profile_program": path,
                "positive_family_n": len(valid),
                "positive_families": valid,
                "overlaps": overlaps,
            }
        )
    result = {
        "status": "PASS" if not errors else "BLOCKED",
        "schema_version": "2.0",
        "profile": str(args.profile.resolve()),
        "profile_sha256": sha256(args.profile),
        "catalog": str(args.catalog.resolve()),
        "catalog_sha256": sha256(args.catalog),
        "review_required_candidates_checked": len(rows),
        "default_candidates_checked": sum(
            candidate.get("release_level") in {"default_broad_candidate", "context_specific_broad_candidate"}
            for candidate in candidates
            if candidate.get("review_required") is True
        ),
        "candidates": rows,
        "errors": errors,
        "broad_presence_rule": "absolute full-feature detection/prevalence and pseudobulk; centered scores are comparative only",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
