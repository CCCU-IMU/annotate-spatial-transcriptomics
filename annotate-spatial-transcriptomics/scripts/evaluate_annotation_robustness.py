#!/usr/bin/env python3
"""Evaluate biological equivalence across reasonable v2.2 parameter perturbations."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

from lineage_controller_lib import (
    deterministic_membership_hash, read_tsv, sha256, write_tsv,
)


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: values[index])
    result = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + end - 1) / 2 + 1
        for index in order[start:end]:
            result[index] = rank
        start = end
    return result


def spearman(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or not left:
        return float("nan")
    a, b = ranks(left), ranks(right)
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    numerator = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    denominator = math.sqrt(
        sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)
    )
    return numerator / denominator if denominator else 1.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", required=True, type=Path)
    ap.add_argument("--variant", required=True, action="append", type=Path)
    ap.add_argument("--deterministic-replicate", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--major-fraction", type=float, default=0.01)
    ap.add_argument("--minimum-spearman", type=float, default=0.95)
    ap.add_argument("--minimum-dice", type=float, default=0.80)
    args = ap.parse_args()

    baseline_rows = read_tsv(args.baseline)
    baseline = {str(row.get("cell_id", "")): row for row in baseline_rows}
    if not baseline or "" in baseline or len(baseline) != len(baseline_rows):
        raise SystemExit("baseline membership is empty or duplicated")
    total = len(baseline)
    base_census = Counter(
        str(row.get("final_broad_label", "")) for row in baseline_rows
        if row.get("final_broad_label")
    )
    major = {
        label for label, count in base_census.items()
        if count / total >= args.major_fraction
    }
    base_fine_census = Counter(
        str(row.get("final_fine_label", "")) for row in baseline_rows
        if row.get("final_fine_label")
    )
    major_fine = {
        label for label, count in base_fine_census.items()
        if count / total >= args.major_fraction
    }
    base_fine_signature: dict[str, tuple[str, str]] = {}
    for label in sorted(major_fine):
        signatures = {
            (
                str(row.get("final_broad_label", "")),
                str(row.get("final_fine_candidate_id", "")),
            )
            for row in baseline_rows
            if row.get("final_fine_label") == label
        }
        if len(signatures) != 1 or not next(iter(signatures))[0]:
            raise SystemExit(
                f"baseline major fine label lacks one stable parent/program signature: {label}"
            )
        base_fine_signature[label] = next(iter(signatures))
    rows: list[dict[str, object]] = []
    errors: list[str] = []
    deterministic_status = "NOT_PROVIDED"
    if args.deterministic_replicate:
        replicate_rows = read_tsv(args.deterministic_replicate)
        if (
            replicate_rows == baseline_rows
            and deterministic_membership_hash(replicate_rows)
            == deterministic_membership_hash(baseline_rows)
        ):
            deterministic_status = "PASS"
        else:
            deterministic_status = "BLOCKED"
            errors.append(
                "deterministic replicate differs in membership order/content or semantic hash"
            )
    for path in args.variant:
        variant_rows = read_tsv(path)
        variant = {str(row.get("cell_id", "")): row for row in variant_rows}
        if set(variant) != set(baseline) or len(variant) != len(variant_rows):
            errors.append(f"{path}: variant does not cover baseline observations exactly")
            continue
        census = Counter(
            str(row.get("final_broad_label", "")) for row in variant_rows
            if row.get("final_broad_label")
        )
        variant_major = {
            label for label, count in census.items()
            if count / total >= args.major_fraction
        }
        labels = sorted(set(base_census) | set(census))
        rho = spearman(
            [base_census[label] for label in labels],
            [census[label] for label in labels],
        )
        dice_by_label: dict[str, float] = {}
        for label in sorted(major):
            left = {cell for cell, row in baseline.items() if row.get("final_broad_label") == label}
            right = {cell for cell, row in variant.items() if row.get("final_broad_label") == label}
            dice_by_label[label] = 2 * len(left & right) / (len(left) + len(right)) if left or right else 1.0
        missing = sorted(major - variant_major)
        added = sorted(variant_major - major)
        failed_dice = sorted(label for label, value in dice_by_label.items() if value < args.minimum_dice)
        missing_major_fine: list[str] = []
        changed_fine_signature: list[str] = []
        for label in sorted(major_fine):
            fine_rows = [row for row in variant_rows if row.get("final_fine_label") == label]
            if not fine_rows:
                missing_major_fine.append(label)
                continue
            signatures = {
                (
                    str(row.get("final_broad_label", "")),
                    str(row.get("final_fine_candidate_id", "")),
                )
                for row in fine_rows
            }
            if signatures != {base_fine_signature[label]}:
                changed_fine_signature.append(label)
        status = "PASS" if (
            not missing and not added and rho >= args.minimum_spearman
            and not failed_dice and not missing_major_fine
            and not changed_fine_signature
        ) else "BLOCKED"
        if status != "PASS":
            errors.append(f"{path}: parameter perturbation is not biologically equivalent")
        rows.append({
            "variant": str(path.resolve()),
            "variant_sha256": sha256(path),
            "broad_census_spearman": rho,
            "missing_major_broads": ";".join(missing),
            "added_major_broads": ";".join(added),
            "failed_dice_broads": ";".join(failed_dice),
            "missing_major_fine_labels": ";".join(missing_major_fine),
            "changed_fine_parent_or_program": ";".join(changed_fine_signature),
            "minimum_major_broad_dice": min(dice_by_label.values(), default=1.0),
            "status": status,
        })

    args.out.mkdir(parents=True, exist_ok=True)
    table = args.out / "parameter_robustness.tsv"
    write_tsv(table, rows)
    manifest = {
        "status": "PASS" if not errors else "BLOCKED",
        "schema_version": "2.2",
        "baseline": {"path": str(args.baseline.resolve()), "sha256": sha256(args.baseline)},
        "variant_n": len(args.variant),
        "technical_determinism": deterministic_status,
        "major_fraction": args.major_fraction,
        "minimum_spearman": args.minimum_spearman,
        "minimum_dice": args.minimum_dice,
        "results": {"path": str(table.resolve()), "sha256": sha256(table)},
        "errors": errors,
    }
    manifest_path = args.out / "parameter_robustness_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
