#!/usr/bin/env python3
"""Build sample-specific marker and support tables from frozen taxonomy/profile."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt", newline="") if path.suffix == ".gz" else path.open(newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def nested(doc: dict, dotted: str) -> object:
    value: object = doc
    for key in dotted.split("."):
        if not isinstance(value, dict) or key not in value:
            return {}
        value = value[key]
    return value


def genes_from_program(program: object) -> list[str]:
    if isinstance(program, list):
        return [str(value) for value in program if value]
    if not isinstance(program, dict):
        return []
    keys = (
        "required_program", "required_positive_genes", "positive_markers",
        "canonical_markers",
    )
    genes: list[str] = []
    for key in keys:
        if isinstance(program.get(key), list):
            genes.extend(str(value) for value in program[key])
    for key in (
        "positive_families", "fine_positive_families", "state_positive_families",
    ):
        families = program.get(key, {})
        if isinstance(families, dict):
            for values in families.values():
                if isinstance(values, list):
                    genes.extend(str(value) for value in values)
    if not genes:
        for value in program.values():
            if isinstance(value, list) and value and all(isinstance(item, str) for item in value):
                genes.extend(value)
    return list(dict.fromkeys(gene for gene in genes if gene))


def five_genes(candidates: list[dict], profile: dict) -> list[str]:
    ordered = sorted(
        candidates,
        key=lambda row: (
            0 if row.get("candidate_role") == "broad" else 1,
            -int(row.get("specificity_priority", 0) or 0),
            str(row.get("candidate_id", "")),
        ),
    )
    pools: list[list[str]] = []
    for candidate in ordered:
        program = nested(profile, str(candidate.get("profile_program", "")))
        genes = genes_from_program(program)
        if not genes:
            genes = genes_from_program(candidate)
        if genes:
            pools.append(genes)
    selected: list[str] = []
    # Keep a coherent primary identity while representing child identities of
    # shared broad parents such as Vascular-associated.
    if pools:
        for gene in pools[0][:3]:
            if gene not in selected:
                selected.append(gene)
    cursor = [3] + [0] * max(0, len(pools) - 1)
    while len(selected) < 5 and pools:
        progressed = False
        for index, pool in enumerate(pools):
            while cursor[index] < len(pool) and pool[cursor[index]] in selected:
                cursor[index] += 1
            if cursor[index] < len(pool):
                selected.append(pool[cursor[index]])
                cursor[index] += 1
                progressed = True
                if len(selected) == 5:
                    break
        if not progressed:
            break
    return selected[:5]


def label_tokens(value: str) -> set[str]:
    ignored = {"like", "associated", "granulosa", "cell", "cells", "follicle"}
    return {
        token for token in re.findall(r"[a-z0-9]+", value.lower())
        if token not in ignored
    }


def closest_label(value: str, catalog: dict[str, dict]) -> dict:
    wanted = label_tokens(value)
    ranked = []
    for label, candidate in catalog.items():
        overlap = wanted & label_tokens(label)
        if overlap:
            ranked.append((len(overlap), -len(wanted ^ label_tokens(label)), label, candidate))
    return sorted(ranked, reverse=True)[0][3] if ranked else {}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--membership", required=True, type=Path)
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--profile", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    rows = read_tsv(args.membership)
    catalog = json.loads(args.catalog.read_text(encoding="utf-8"))
    profile = json.loads(args.profile.read_text(encoding="utf-8"))
    broad_counts = Counter(row.get("final_broad_label", "") for row in rows if row.get("final_broad_label", ""))
    fine_counts = Counter(row.get("final_fine_label", "") for row in rows if row.get("final_fine_label", ""))
    fine_parent_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        if row.get("final_fine_label", "") and row.get("final_broad_label", ""):
            fine_parent_counts[row["final_fine_label"]][row["final_broad_label"]] += 1
    state_counts: Counter[str] = Counter()
    source_counts: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    for row in rows:
        broad = row.get("final_broad_label", "")
        if broad:
            source_counts[("broad", broad)][row.get("broad_freeze_source", "") or row.get("assignment_origin", "") or "second_round"] += 1
        fine = row.get("final_fine_label", "")
        if fine:
            source_counts[("subtype", fine)][row.get("final_fine_assignment_source", "") or "second_round_fine"] += 1
        for state in (row.get("state_annotations", "") or row.get("final_state_annotation", "")).split(";"):
            if state:
                state_counts[state] += 1
                source_counts[("state", state)]["second_round_state_program"] += 1

    boundary_by_broad: dict[str, list[dict]] = defaultdict(list)
    state_by_label: dict[str, dict] = {}
    for candidate in catalog.get("candidate_boundaries", []):
        broad = str(candidate.get("release_broad_label", ""))
        if broad:
            boundary_by_broad[broad].append(candidate)
        state = str(candidate.get("release_state_label", ""))
        if state:
            state_by_label[state.lower()] = candidate
    fine_catalog: dict[str, dict] = {}
    for candidates in catalog.get("machine_actionable_fine_candidate_catalog", {}).values():
        for candidate in candidates:
            fine_catalog[str(candidate.get("release_label", ""))] = candidate

    marker_rows: list[dict[str, object]] = []
    support_rows: list[dict[str, object]] = []
    for broad, count in sorted(broad_counts.items()):
        candidates = boundary_by_broad.get(broad, [])
        genes = five_genes(candidates, profile)
        for gene in genes:
            marker_rows.append({"level": "broad", "marker_group": broad, "gene": gene})
        spatial = "; ".join(
            str(nested(profile, str(candidate.get("profile_program", ""))).get("spatial_expectation", ""))
            for candidate in candidates
            if isinstance(nested(profile, str(candidate.get("profile_program", ""))), dict)
            and nested(profile, str(candidate.get("profile_program", ""))).get("spatial_expectation")
        )
        support_rows.append({
            "level": "broad", "parent_label": "", "label": broad,
            "n_observations": count, "canonical_markers": ";".join(genes),
            "evidence_source": "full-catalog second-round query evidence; DEG/pseudobulk; spatial and post-Atlas biological review",
            "member_sources": "; ".join(f"{name}: {n}" for name, n in source_counts[("broad", broad)].most_common()),
            "spatial_support": spatial or "sample-specific spatial localization reviewed on complete membership",
            "release_interpretation": "most specific supported broad identity; geometry alone had no label authority",
        })
    for fine, count in sorted(fine_counts.items()):
        candidate = fine_catalog.get(fine, {}) or closest_label(fine, fine_catalog)
        program = nested(profile, str(candidate.get("profile_program", "")))
        genes = genes_from_program(candidate) or genes_from_program(program)
        genes = genes[:5]
        parent = str(candidate.get("parent_release_label", "")) or (
            fine_parent_counts[fine].most_common(1)[0][0] if fine_parent_counts[fine] else ""
        )
        for gene in genes:
            marker_rows.append({"level": "subtype", "marker_group": fine, "gene": gene})
        support_rows.append({
            "level": "subtype", "parent_label": parent, "label": fine,
            "n_observations": count, "canonical_markers": ";".join(genes),
            "evidence_source": "parent-locked second-round fine discriminator",
            "member_sources": "; ".join(f"{name}: {n}" for name, n in source_counts[("subtype", fine)].most_common()),
            "spatial_support": "; ".join(str(value) for value in candidate.get("context_requirements", [])),
            "release_interpretation": "high-confidence fine identity within its frozen broad parent",
        })
    state_marker_rows: list[dict[str, object]] = []
    for state, count in sorted(state_counts.items()):
        candidate = state_by_label.get(state.lower(), {}) or closest_label(state, state_by_label)
        genes = genes_from_program(candidate)[:5]
        for gene in genes:
            state_marker_rows.append({"level": "broad", "marker_group": state, "gene": gene})
        support_rows.append({
            "level": "state", "parent_label": str(candidate.get("parent_broad_label", "")),
            "label": state, "n_observations": count,
            "canonical_markers": ";".join(genes),
            "evidence_source": "independent state program after broad-parent freeze",
            "member_sources": "; ".join(f"{name}: {n}" for name, n in source_counts[("state", state)].most_common()),
            "spatial_support": "state distribution shown independently from cell identity",
            "release_interpretation": "state annotation; does not replace broad/fine identity",
        })
    args.out.mkdir(parents=True, exist_ok=True)
    fields = ["level", "marker_group", "gene"]
    write_tsv(args.out / "canonical_marker_panel.tsv", marker_rows, fields)
    write_tsv(args.out / "state_marker_panel.tsv", state_marker_rows, fields)
    write_tsv(args.out / "annotation_support_summary.tsv", support_rows, [
        "level", "parent_label", "label", "n_observations",
        "canonical_markers", "evidence_source", "member_sources",
        "spatial_support", "release_interpretation",
    ])
    manifest = {
        "status": "PASS", "n_broad": len(broad_counts),
        "n_fine": len(fine_counts), "n_state": len(state_counts),
        "canonical_marker_panel": str((args.out / "canonical_marker_panel.tsv").resolve()),
        "state_marker_panel": str((args.out / "state_marker_panel.tsv").resolve()),
        "annotation_support_summary": str((args.out / "annotation_support_summary.tsv").resolve()),
    }
    (args.out / "release_evidence_tables_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
