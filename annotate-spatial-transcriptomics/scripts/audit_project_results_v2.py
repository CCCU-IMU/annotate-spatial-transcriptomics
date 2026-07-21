#!/usr/bin/env python3
"""Read-only audit of an annotation project's state and results surfaces."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path

from evidence_schema_lib import sha256


BIOLOGICAL = {"defined_broad_only", "defined_fine"}
QC = {"qc_holdout", "low_information_qc_holdout", "pending_qc", "unknown_candidate"}
DEFAULT_ALIASES = {"Vascular/endothelial", "Vascular/perivascular", "Pericyte/mural"}
DEFAULT_FORBIDDEN_FINE = [
    "low-information", "low information", "qc", "holdout", "pending", "unknown",
    "unresolved", "candidate", "review", "technical", "ambient-only",
]


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return {}


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(encoding="utf-8", newline="")


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def choose_ledger(root: Path) -> tuple[Path | None, str]:
    for rel, role in (
        ("state/cell_ledger.tsv.gz", "committed"),
        ("state/cell_ledger.tsv", "committed"),
    ):
        path = root / rel
        if path.is_file():
            return path, role
    candidates = sorted(
        {
            path
            for pattern in ("*final*ledger*.tsv*", "*candidate*ledger*.tsv*")
            for path in (root / "results").rglob(pattern)
            if ".summary." not in path.name and "census" not in path.name
        },
        key=lambda path: ("candidate" not in path.name.lower(), path.stat().st_mtime),
        reverse=True,
    ) if (root / "results").is_dir() else []
    return (candidates[0], "uncommitted_candidate") if candidates else (None, "missing")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("project_root", type=Path)
    ap.add_argument("--biological-profile", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    root = args.project_root.resolve()
    blockers: list[str] = []
    warnings: list[str] = []
    findings: list[dict] = []

    project_path = root / "config/project.json"
    project = read_json(project_path)
    if not project:
        blockers.append("missing_or_unreadable_project_config")
    framework = str(project.get("framework_version", ""))
    if framework and not framework.startswith("2."):
        warnings.append(f"legacy_framework:{framework}")

    profile = read_json(args.biological_profile.resolve()) if args.biological_profile else {}
    taxonomy = profile.get("release_taxonomy", {}) if isinstance(profile.get("release_taxonomy", {}), dict) else {}
    aliases = set(taxonomy.get("broad_aliases", {}).keys()) or DEFAULT_ALIASES
    forbidden_fine = taxonomy.get("forbidden_fine_semantics", []) or DEFAULT_FORBIDDEN_FINE
    forbidden_fine_pattern = re.compile("|".join(re.escape(str(value)) for value in forbidden_fine), re.I)

    ledger, ledger_role = choose_ledger(root)
    census: Counter[tuple[str, str, str]] = Counter()
    analysis_n = qc_n = 0
    fields: list[str] = []
    if ledger is None:
        blockers.append("no_committed_or_candidate_cell_ledger")
    else:
        fields, rows = read_tsv(ledger)
        state_col = "final_state" if "final_state" in fields else "state"
        broad_col = "final_broad_label" if "final_broad_label" in fields else "broad_label"
        fine_col = "final_fine_label" if "final_fine_label" in fields else "fine_label"
        for row in rows:
            if row.get("analysis_scope", "analysis_set") == "excluded_initial_qc":
                continue
            analysis_n += 1
            state, broad, fine = row.get(state_col, ""), row.get(broad_col, ""), row.get(fine_col, "")
            census[(state, broad, fine)] += 1
            if state in QC:
                qc_n += 1
        if ledger_role != "committed":
            blockers.append("latest_annotation_ledger_is_uncommitted_candidate")
        alias_counts = Counter()
        for (_, broad, _), count in census.items():
            if broad in aliases:
                alias_counts[broad] += count
        if alias_counts:
            blockers.append("legacy_or_child_label_used_as_release_broad")
            findings.append({"finding": "legacy_broad_aliases", "counts": dict(alias_counts)})
        bad_fine: Counter[str] = Counter()
        for (_, _, fine), count in census.items():
            if fine and forbidden_fine_pattern.search(fine):
                bad_fine[fine] += count
        if bad_fine:
            blockers.append("non_biological_semantics_published_as_fine_labels")
            findings.append({"finding": "non_biological_fine_labels", "counts": dict(bad_fine)})
        noncanonical = Counter()
        for (state, _, _), count in census.items():
            if state not in BIOLOGICAL | QC:
                noncanonical[state] += count
        if noncanonical:
            warnings.append("noncanonical_final_states")
            findings.append({"finding": "noncanonical_final_states", "counts": dict(noncanonical)})

    completion_path = root / "provenance/completion_gate.json"
    completion = read_json(completion_path)
    completion_status = completion.get("status", "MISSING")
    if completion_status == "PASS" and blockers:
        blockers.append("completion_gate_pass_conflicts_with_result_audit")
    if framework.startswith("2."):
        for rel in (
            "config/annotation_contract.json",
            "provenance/broad_family_evidence_validation.json",
            "provenance/release_taxonomy_audit.json",
        ):
            if not (root / rel).is_file():
                blockers.append("missing_v2_gate:" + rel)

    state_validation = read_json(root / "provenance/state_validation.json")
    if ledger_role == "committed":
        validated = state_validation.get("validated_files", {})
        relative = str(ledger.relative_to(root))
        if state_validation.get("status") != "PASS" or validated.get(relative) != sha256(ledger):
            blockers.append("state_validation_missing_or_stale_for_committed_cell_ledger")

    derived_path = root / "state/derived_expression_registry.tsv"
    if derived_path.is_file():
        _, derived = read_tsv(derived_path)
        project_id = str(project.get("project_id", ""))
        leaking = [row.get("artifact_id", "") for row in derived if row.get("external_reference", "").lower() not in {"true", "1", "yes"} and row.get("project_id", "") != project_id]
        if leaking:
            blockers.append("cross_project_query_expression_artifacts")
            findings.append({"finding": "cross_project_query_expression_artifacts", "artifact_ids": leaking})

    run_registry = root / "state/run_registry.tsv"
    if run_registry.is_file():
        _, runs = read_tsv(run_registry)
        nonterminal = [
            {"run_id": row.get("run_id", ""), "stage": row.get("stage", ""), "status": row.get("status", "")}
            for row in runs
            if row.get("status", "") not in {"validated_done", "failed_preserved", "cancelled_preserved"}
        ]
        if nonterminal:
            blockers.append("run_registry_contains_nonterminal_runs")
            findings.append({"finding": "nonterminal_runs", "runs": nonterminal})
    incident_registry = root / "provenance/incidents/incident_registry.tsv"
    if incident_registry.is_file():
        _, incidents = read_tsv(incident_registry)
        closed = {"repaired", "repaired_validated", "resolved", "closed", "waived_with_rationale"}
        open_incidents = [row.get("incident_id", "") for row in incidents if row.get("status", "").strip().lower() not in closed]
        if open_incidents:
            blockers.append("incident_registry_contains_open_incidents")
            findings.append({"finding": "open_incidents", "incident_ids": open_incidents})

    empty_results = []
    result_files = []
    if (root / "results").is_dir():
        for path in (root / "results").rglob("*"):
            if path.is_file():
                result_files.append(path)
                if path.stat().st_size == 0:
                    empty_results.append(str(path.relative_to(root)))
    if empty_results:
        warnings.append(f"empty_result_files:{len(empty_results)}")

    result = {
        "schema_version": "2.0",
        "status": "PASS" if not blockers else "BLOCKED",
        "project_root": str(root),
        "framework_version": framework,
        "biological_profile": str(args.biological_profile.resolve()) if args.biological_profile else "",
        "biological_profile_sha256": sha256(args.biological_profile.resolve()) if args.biological_profile else "",
        "completion_gate_status": completion_status,
        "ledger": str(ledger) if ledger else "",
        "ledger_role": ledger_role,
        "ledger_sha256": sha256(ledger) if ledger else "",
        "analysis_n": analysis_n,
        "residual_qc_n": qc_n,
        "residual_qc_fraction": qc_n / analysis_n if analysis_n else None,
        "annotation_census": [
            {"state": state, "broad_label": broad, "fine_label": fine, "n": count}
            for (state, broad, fine), count in census.most_common()
        ],
        "result_file_n": len(result_files),
        "blockers": sorted(set(blockers)),
        "warnings": sorted(set(warnings)),
        "findings": findings,
        "read_only_audit": True,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if not blockers else 2


if __name__ == "__main__":
    raise SystemExit(main())
