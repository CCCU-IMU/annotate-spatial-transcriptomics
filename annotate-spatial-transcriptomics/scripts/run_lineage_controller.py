#!/usr/bin/env python3
"""Run the phase-authorized v2.2 annotation controller."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from evidence_schema_lib import sha256


PHASES = (
    "whole_tissue_partition",
    "cluster_cohort_recluster",
    "local_mixed_subcluster_split",
    "merge_and_freeze_broad",
    "atlas_and_completeness_review",
    "materialize_final_release",
)
WHOLE_TISSUE_FORK_WORKER_CAP = 64
CANONICAL_SCRIPTS = (
    "run_observation_lineage_scoring.R",
    "derive_candidate_local_subsets.R",
    "close_exact_remainders.py",
    "build_whole_tissue_cohort_plan.py",
    "adjudicate_second_round_subclusters.py",
    "build_candidate_context_evidence.py",
    "merge_and_freeze_broad_membership.py",
    "route_global_atlas_v2.py",
    "validate_global_atlas_v2.py",
    "apply_post_merge_atlas_routing.py",
    "review_post_merge_unresolved_components.py",
    "audit_post_merge_completeness.py",
    "validate_sheep_ovary_biological_quality.py",
    "apply_sheep_ovary_follicle_roi_repair.py",
    "screen_rare_cell_programs.R",
    "screen_spatial_foci.py",
    "materialize_oocyte_cluster_membership.py",
    "apply_cell_id_membership_patch.py",
    "materialize_parent_locked_fine_proposals.py",
    "materialize_final_release_v2_2.py",
    "evaluate_annotation_robustness.py",
    "run_lineage_controller.py",
)
FORBIDDEN_SCORING_COLUMNS = {
    "provisional_broad", "provisional_broad_after_score_freeze",
    "initial_broad_label", "broad_label", "fine_label", "final_broad_label",
    "final_fine_label", "celltype", "cell_type", "annotation",
    "historical_label", "repair_membership", "atlas_label",
}


def artifact(path: Path, role: str = "runtime_input") -> dict[str, str]:
    if not path.is_file():
        raise RuntimeError(f"missing artifact: {path}")
    return {
        "path": str(path.resolve()), "sha256": sha256(path),
        "artifact_role": role,
    }


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def bound_artifact_paths(records: object, label: str) -> list[Path]:
    if not isinstance(records, list):
        raise RuntimeError(f"{label} registry is missing")
    paths: list[Path] = []
    for index, record in enumerate(records, 1):
        if not isinstance(record, dict):
            raise RuntimeError(f"{label} artifact {index} is malformed")
        path = Path(str(record.get("path", "")))
        if not path.is_file() or record.get("sha256") != sha256(path):
            raise RuntimeError(f"{label} artifact {index} is missing or stale")
        paths.append(path)
    return paths


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def run(command: list[str], log: Path, allowed_codes: tuple[int, ...] = (0,)) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        result = subprocess.run(
            command, stdout=handle, stderr=subprocess.STDOUT,
            text=True, check=False,
        )
    if result.returncode not in allowed_codes:
        raise RuntimeError(
            f"controller command failed ({result.returncode}); inspect {log}"
        )


def resolve_bound(contract_path: Path, record: dict, label: str) -> Path:
    path = Path(str(record.get("path", "")))
    if not path.is_absolute():
        path = (contract_path.parent / path).resolve()
    if (
        record.get("artifact_role") not in {"runtime_input", "external_reference"}
        or not path.is_file()
        or record.get("sha256") != sha256(path)
    ):
        raise RuntimeError(f"bound {label} is missing, stale or role-forbidden")
    return path


def validate_contract(contract_path: Path) -> tuple[dict, dict[str, Path]]:
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    controller = contract.get("canonical_lineage_controller", {})
    if contract.get("schema_version") != "2.0":
        raise RuntimeError("annotation contract is not schema 2.0")
    if controller.get("controller_version") != "2.2.0":
        raise RuntimeError("annotation contract does not bind controller v2.2.0")
    if controller.get("phase_order") != list(PHASES):
        raise RuntimeError("annotation contract does not bind the staged v2.2 architecture")
    script_dir = Path(__file__).resolve().parent
    paths = {name: script_dir / name for name in CANONICAL_SCRIPTS}
    for name, path in paths.items():
        record = controller.get("scripts", {}).get(name, {})
        if (
            Path(str(record.get("path", ""))).resolve() != path.resolve()
            or not path.is_file()
            or record.get("sha256") != sha256(path)
        ):
            raise RuntimeError(f"contract does not bind installed canonical {name}")
    separately_bound = {
        "select_lineage_resolution.py": controller.get("resolution_selector", {}),
        "build_resolution_grid_evidence.py": controller.get("resolution_evidence_builder", {}),
        "discover_unmodeled_lineages.py": controller.get("unmodeled_discovery", {}),
        "run_seurat_cohort_recluster.R": controller.get("cohort_reclustering", {}).get("entrypoint", {}),
        "run_seurat_cohort_recluster_impl.R": controller.get("cohort_reclustering", {}).get("implementation", {}),
    }
    for name, record in separately_bound.items():
        path = script_dir / name
        if (
            Path(str(record.get("path", ""))).resolve() != path.resolve()
            or not path.is_file()
            or record.get("sha256") != sha256(path)
        ):
            raise RuntimeError(f"contract does not bind installed canonical {name}")
        paths[name] = path
    dependency = script_dir / "lineage_controller_lib.py"
    record = controller.get("dependencies", {}).get("lineage_controller_lib.py", {})
    if (
        Path(str(record.get("path", ""))).resolve() != dependency.resolve()
        or not dependency.is_file()
        or record.get("sha256") != sha256(dependency)
    ):
        raise RuntimeError("contract does not bind lineage_controller_lib.py")
    policy = contract.get("artifact_role_policy", {})
    if (
        policy.get("runtime_allowed_roles") != ["runtime_input", "external_reference"]
        or policy.get("failed_diagnostic_runtime_forbidden") is not True
    ):
        raise RuntimeError("contract does not forbid failed diagnostic runtime inputs")
    selected = contract.get("selected_input_snapshot", {})
    if selected.get("artifact_role") != "runtime_input":
        raise RuntimeError("selected input snapshot is not runtime_input")
    selected_path = Path(str(selected.get("path", ""))).resolve()
    if not selected_path.is_file() or len(str(selected.get("sha256", ""))) != 64:
        raise RuntimeError("selected runtime input snapshot is missing or unbound")
    diagnostic = policy.get("diagnostic_registry")
    if diagnostic:
        diagnostic_path = Path(str(diagnostic.get("path", "")))
        if not diagnostic_path.is_file() or diagnostic.get("sha256") != sha256(diagnostic_path):
            raise RuntimeError("failed diagnostic registry is missing or stale")
        failed = [
            Path(row.get("artifact_root", "")).resolve()
            for row in read_tsv(diagnostic_path)
            if row.get("artifact_role") == "failed_diagnostic"
            and row.get("artifact_root")
        ]
        if any(selected_path == root or root in selected_path.parents for root in failed):
            raise RuntimeError("selected input is inside a failed diagnostic artifact")
    paths["lineage_controller_lib.py"] = dependency
    paths["profile"] = resolve_bound(
        contract_path, contract.get("biological_profile", {}), "biological profile"
    )
    paths["workflow_profile"] = resolve_bound(
        contract_path, contract.get("workflow_profile", {}), "workflow profile"
    )
    paths["catalog"] = resolve_bound(
        contract_path, contract.get("candidate_catalog", {}), "candidate catalog"
    )
    paths["selected_input"] = selected_path
    input_scope = contract.get("input_scope", {})
    paths["analysis_set"] = resolve_bound(
        contract_path, input_scope.get("analysis_set", {}), "analysis_set"
    )
    paths["excluded_initial_qc"] = resolve_bound(
        contract_path, input_scope.get("excluded_initial_qc", {}),
        "excluded_initial_qc",
    )
    paths["whole_tissue_partitions"] = resolve_bound(
        contract_path,
        contract.get("whole_tissue_partition", {}).get("partition_grid", {}),
        "whole-tissue partition grid",
    )
    full_record = input_scope.get("full_object", {})
    full_object = Path(str(full_record.get("path", ""))).resolve()
    if (
        full_record.get("artifact_role") != "runtime_input"
        or full_object != selected_path
        or full_record.get("sha256") != selected.get("sha256")
    ):
        raise RuntimeError("input scope full_object differs from selected snapshot")
    return contract, paths


def validate_selected_input(path: Path, paths: dict[str, Path]) -> None:
    if path.resolve() != paths["selected_input"]:
        raise RuntimeError(
            "controller input differs from the contract-bound runtime snapshot"
        )


def runner_parameters(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    if not rows or not {"parameter", "value"}.issubset(rows[0]):
        raise RuntimeError("cohort runner manifest is empty or malformed")
    return {str(row["parameter"]): str(row["value"]) for row in rows}


def recluster_cache_fingerprint(
    contract: dict, paths: dict[str, Path], args, grid: list[float]
) -> tuple[str, dict[str, object]]:
    """Bind only deterministic clustering inputs, never annotation semantics."""
    payload: dict[str, object] = {
        "schema_version": "2.2",
        "selected_input_sha256": str(
            contract["selected_input_snapshot"]["sha256"]
        ),
        "query_membership_sha256": sha256(args.membership),
        "cohort_id": str(args.cohort_id),
        "source_initial_cluster": str(args.source_initial_cluster),
        "candidate_resolutions": [float(value) for value in grid],
        "resolution_contract": str(args.resolution_contract),
        "seed": int(args.seed),
        "normalization": "SCT",
        "sct_method": "glmGamPoi",
        "clustering_path": "raw_counts_SCTv2_PCA_SNN_Leiden",
        "clustering_scripts": {
            name: sha256(paths[name])
            for name in (
                "run_seurat_cohort_recluster.R",
                "run_seurat_cohort_recluster_impl.R",
            )
        },
    }
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest(), payload


def failed_diagnostic_roots(contract: dict) -> list[Path]:
    record = contract.get("artifact_role_policy", {}).get("diagnostic_registry")
    if not isinstance(record, dict):
        return []
    path = Path(str(record.get("path", "")))
    if not path.is_file() or record.get("sha256") != sha256(path):
        raise RuntimeError("failed diagnostic registry is missing or stale")
    return [
        Path(str(row.get("artifact_root", ""))).resolve()
        for row in read_tsv(path)
        if row.get("artifact_role") == "failed_diagnostic"
        and row.get("artifact_root")
    ]


def validate_recluster_cache(
    manifest_path: Path,
    expected_fingerprint: str,
    contract: dict,
) -> tuple[Path, Path, Path, Path]:
    manifest_path = manifest_path.resolve()
    if any(
        manifest_path == root or root in manifest_path.parents
        for root in failed_diagnostic_roots(contract)
    ):
        raise RuntimeError("failed_diagnostic partition cache is forbidden")
    doc = json.loads(manifest_path.read_text(encoding="utf-8"))
    cached_payload = doc.get("execution_fingerprint", {})
    cached_payload_sha = hashlib.sha256(json.dumps(
        cached_payload, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")).hexdigest()
    if (
        doc.get("status") != "PASS"
        or doc.get("artifact_role") != "derived_partition_cache"
        or doc.get("label_authority") is not False
        or doc.get("execution_fingerprint_sha256") != expected_fingerprint
        or cached_payload_sha != expected_fingerprint
    ):
        raise RuntimeError("recluster partition cache fingerprint is incompatible")
    artifacts = doc.get("artifacts", {})
    required = ("query_rds", "run_manifest", "zero_count", "partition_grid")
    resolved: dict[str, Path] = {}
    for name in required:
        record = artifacts.get(name, {})
        path = Path(str(record.get("path", "")))
        if (
            record.get("artifact_role") != "derived_partition_cache"
            or not path.is_file()
            or record.get("sha256") != sha256(path)
        ):
            raise RuntimeError(f"recluster cache artifact is missing or stale: {name}")
        resolved[name] = path.resolve()
    resolution_records = artifacts.get("resolution_memberships", [])
    if len(resolution_records) != len(
        cached_payload.get("candidate_resolutions", [])
    ):
        raise RuntimeError("recluster cache resolution membership grid is incomplete")
    for record in resolution_records:
        path = Path(str(record.get("path", "")))
        if (
            record.get("artifact_role") != "derived_partition_cache"
            or not path.is_file()
            or record.get("sha256") != sha256(path)
        ):
            raise RuntimeError("recluster cache resolution membership is stale")
    return (
        resolved["query_rds"], resolved["run_manifest"],
        resolved["zero_count"], resolved["partition_grid"],
    )


def write_recluster_cache_manifest(
    path: Path,
    fingerprint: str,
    fingerprint_payload: dict[str, object],
    query_rds: Path,
    run_manifest: Path,
    zero_count: Path,
    grid_partitions: Path,
    resolution_memberships: list[Path],
) -> None:
    cache_record = lambda value: artifact(value, "derived_partition_cache")
    doc = {
        "status": "PASS",
        "schema_version": "2.2",
        "artifact_role": "derived_partition_cache",
        "label_authority": False,
        "annotation_fields_used_for_partition_or_scoring": False,
        "execution_fingerprint_sha256": fingerprint,
        "execution_fingerprint": fingerprint_payload,
        "artifacts": {
            "query_rds": cache_record(query_rds),
            "run_manifest": cache_record(run_manifest),
            "zero_count": cache_record(zero_count),
            "partition_grid": cache_record(grid_partitions),
            "resolution_memberships": [
                cache_record(value) for value in resolution_memberships
            ],
        },
    }
    path.write_text(
        json.dumps(doc, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def validate_blinded_table(path: Path) -> None:
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = {str(value).lower() for value in (reader.fieldnames or [])}
    leaked = sorted(fields.intersection(FORBIDDEN_SCORING_COLUMNS))
    if leaked:
        raise RuntimeError(
            "scoring input exposes provisional/historical annotation columns: "
            + ", ".join(leaked)
        )


def validate_cohort_plan_binding(args) -> tuple[dict, dict[str, str]]:
    """Bind one second-round run to one exact first-round cluster cohort."""
    whole = json.loads(args.whole_manifest.read_text(encoding="utf-8"))
    if not controller_manifest_bound(
        whole, args.contract, "whole_tissue_partition"
    ):
        raise RuntimeError(
            "second-round cohort requires the bound whole-tissue controller manifest"
        )
    plan_record = whole.get("cohort_plan", {})
    plan_path = Path(str(plan_record.get("path", "")))
    if (
        not plan_path.is_file()
        or plan_record.get("sha256") != sha256(plan_path)
    ):
        raise RuntimeError("whole-tissue cohort plan is missing or stale")
    matches = [
        row for row in read_tsv(plan_path)
        if str(row.get("cohort_id", "")) == args.cohort_id
    ]
    if len(matches) != 1:
        raise RuntimeError("cohort_id is not unique in the whole-tissue plan")
    plan_row = matches[0]
    if (
        str(plan_row.get("source_initial_cluster", ""))
        != args.source_initial_cluster
        or str(plan_row.get("provisional_status", ""))
        != args.provisional_status
        or str(plan_row.get("provisional_broad_after_score_freeze", ""))
        != args.provisional_broad
    ):
        raise RuntimeError(
            "second-round cohort metadata differs from the frozen provisional plan"
        )
    planned_membership = Path(str(plan_row.get("membership_path", "")))
    if (
        not planned_membership.is_file()
        or plan_row.get("membership_sha256") != sha256(planned_membership)
    ):
        raise RuntimeError("planned initial-cluster membership is missing or stale")
    planned_ids = {
        str(row.get("cell_id", "")) for row in read_tsv(planned_membership)
    }
    supplied_ids = {
        str(row.get("cell_id", "")) for row in read_tsv(args.membership)
    }
    if (
        not planned_ids
        or "" in planned_ids
        or planned_ids != supplied_ids
        or sha256(args.membership) != plan_row.get("membership_sha256")
    ):
        raise RuntimeError(
            "second-round membership is not the exact bound initial-cluster cohort"
        )
    return whole, plan_row


def materialize_selected_neighbors(
    source: Path, selection: dict, destination: Path
) -> None:
    selected = float(selection["selected_resolution"])
    values = [float(value) for value in selection["selected_and_neighbors"]]
    if len(values) != 3 or values[0] != selected or len(set(values)) != 3:
        raise RuntimeError("selector did not provide selected plus two neighbors")
    roles = {values[0]: "selected", values[1]: "neighbor_1", values[2]: "neighbor_2"}
    rows = read_tsv(source)
    output: list[dict[str, object]] = []
    coverage: dict[float, set[str]] = {value: set() for value in values}
    for row in rows:
        try:
            resolution = float(row.get("resolution", ""))
        except ValueError:
            continue
        if resolution not in roles:
            continue
        cell = str(row.get("cell_id", ""))
        if not cell or cell in coverage[resolution]:
            raise RuntimeError("partition has an empty or duplicate observation")
        coverage[resolution].add(cell)
        output.append({
            "cell_id": cell,
            "boundary_id": row.get("boundary_id", "whole_tissue"),
            "resolution": resolution,
            "cluster": row.get("cluster", ""),
            "resolution_role": roles[resolution],
        })
    baseline = coverage[selected]
    if not baseline or any(value != baseline for value in coverage.values()):
        raise RuntimeError("selected/neighbor partitions do not cover identical observations")
    write_tsv(
        destination, output,
        ["cell_id", "boundary_id", "resolution", "cluster", "resolution_role"],
    )


def materialize_selected_grid_evidence(
    grid_scoring: Path, selection: dict, destination: Path
) -> tuple[Path, Path]:
    """Reuse label-blind grid evidence without a second whole-object scorer run."""
    selected = float(selection["selected_resolution"])
    neighbor_values = {
        float(value) for value in selection["selected_and_neighbors"]
    }
    source_evidence = (
        grid_scoring / "tables/cluster_candidate_multichannel_evidence.tsv.gz"
    )
    evidence_rows = [
        row for row in read_tsv(source_evidence)
        if float(row.get("resolution", "nan")) == selected
    ]
    if not evidence_rows:
        raise RuntimeError("selected resolution has no frozen grid evidence")
    for row in evidence_rows:
        row["resolution_role"] = "selected"
    destination.mkdir(parents=True, exist_ok=True)
    selected_evidence = destination / "tables/cluster_candidate_multichannel_evidence.tsv.gz"
    write_tsv(selected_evidence, evidence_rows, list(evidence_rows[0]))

    source_programs = (
        grid_scoring / "tables/resolution_deg_coexpression_programs.tsv.gz"
    )
    selected_programs = destination / "tables/resolution_deg_coexpression_programs.tsv.gz"
    program_rows = [
        row for row in read_tsv(source_programs)
        if float(row.get("resolution", "nan")) in neighbor_values
    ]
    program_fields = list(program_rows[0]) if program_rows else [
        "program_id", "resolution", "source_boundary", "source_cluster",
        "n_observations", "genes", "coexpressed_gene_count",
        "mean_top_log2fc", "mean_detection_difference",
        "catalog_marker_overlap_fraction", "spatially_coherent",
        "excluded_program_classes", "candidate_status",
    ]
    write_tsv(selected_programs, program_rows, program_fields)
    manifest = {
        "status": "PASS", "schema_version": "2.2",
        "stage": "label_blind_selected_grid_evidence_view",
        "formal_membership_written": False,
        "observation_scores_written": False,
        "historical_labels_read": False,
        "selected_resolution": selected,
        "selected_and_neighbors": sorted(neighbor_values),
        "source_full_grid_scoring": artifact(
            grid_scoring / "observation_scoring_manifest.json"
        ),
        "cluster_multichannel_evidence": artifact(selected_evidence),
        "resolution_deg_coexpression_programs": artifact(selected_programs),
    }
    manifest_path = destination / "selected_grid_evidence_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return selected_evidence, manifest_path


def build_cohort_grid(output: Path, cohort_id: str, grid: list[float]) -> Path:
    rows: list[dict[str, object]] = []
    baseline: set[str] | None = None
    for index, resolution in enumerate(grid):
        tag = str(resolution).replace(".", "p")
        source = output / "tables" / f"framework_res{tag}_clusters.tsv"
        records = read_tsv(source)
        ids = {str(row.get("cell_id", "")) for row in records}
        if not ids or "" in ids or (baseline is not None and ids != baseline):
            raise RuntimeError("cohort resolution grids do not cover identical observations")
        baseline = ids
        role = ("selected", "neighbor_1", "neighbor_2")[index] if index < 3 else "grid"
        for row in records:
            rows.append({
                "cell_id": row["cell_id"], "boundary_id": cohort_id,
                "resolution": resolution, "cluster": row["cluster"],
                "resolution_role": role,
            })
    path = output / "tables" / "cohort_partition_grid.tsv.gz"
    write_tsv(
        path, rows,
        ["cell_id", "boundary_id", "resolution", "cluster", "resolution_role"],
    )
    return path


def stage_authority(
    phase: str, contract_path: Path, paths: dict[str, Path], output: Path,
    **artifacts: Path | list[Path] | None,
) -> Path:
    record = {
        "schema_version": "2.2",
        "mode": "stage_authority",
        "phase": phase,
        "annotation_contract": str(contract_path.resolve()),
        "annotation_contract_sha256": sha256(contract_path),
        "scripts": {
            name: artifact(paths[name]) for name in CANONICAL_SCRIPTS
        },
        "dependencies": {
            "lineage_controller_lib.py": artifact(paths["lineage_controller_lib.py"])
        },
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    for key, path in artifacts.items():
        if isinstance(path, list):
            record[key] = [artifact(item) for item in path]
        else:
            record[key] = artifact(path) if path else None
    authority_path = output / "stage_authority.json"
    authority_path.parent.mkdir(parents=True, exist_ok=True)
    authority_path.write_text(
        json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return authority_path


def run_scorer(
    args, paths: dict[str, Path], rds: Path, partitions: Path, output: Path,
    grid_only: bool = False, analysis_membership: Path | None = None,
    assay: str = "",
) -> None:
    scoring_workers = (
        min(args.scoring_workers, WHOLE_TISSUE_FORK_WORKER_CAP)
        if analysis_membership else args.scoring_workers
    )
    command = [
        args.rscript, str(paths["run_observation_lineage_scoring.R"]),
        "--rds", str(rds.resolve()), "--profile", str(paths["profile"]),
        "--catalog", str(paths["catalog"]), "--partitions", str(partitions.resolve()),
        "--out", str(output), "--seed", str(args.seed),
        "--workers", str(scoring_workers),
    ]
    if grid_only:
        command.extend(["--grid-evidence-only", "true"])
    if analysis_membership:
        command.extend(["--analysis-membership", str(analysis_membership.resolve())])
    if assay:
        command.extend(["--assay", assay])
    run(command, output.parent / "logs" / f"{output.name}.log")


def phase_whole(args, contract: dict, paths: dict[str, Path]) -> dict:
    validate_selected_input(args.rds, paths)
    validate_blinded_table(args.partitions)
    if (
        args.partitions.resolve() != paths["whole_tissue_partitions"]
        or sha256(args.partitions) != sha256(paths["whole_tissue_partitions"])
    ):
        raise RuntimeError("whole-tissue partitions differ from the contract-bound grid")
    analysis_ids = {row.get("cell_id", "") for row in read_tsv(paths["analysis_set"])}
    partition_rows = read_tsv(args.partitions)
    coverage: dict[str, set[str]] = {}
    for row in partition_rows:
        resolution = str(row.get("resolution", ""))
        cell = str(row.get("cell_id", ""))
        cluster = str(row.get("cluster", ""))
        if not resolution or not cell or not cluster:
            raise RuntimeError("whole-tissue grid contains an incomplete partition row")
        if cell in coverage.setdefault(resolution, set()):
            raise RuntimeError("whole-tissue grid duplicates an observation within a resolution")
        coverage[resolution].add(cell)
    if not coverage or any(ids != analysis_ids for ids in coverage.values()):
        raise RuntimeError("every whole-tissue grid partition must exactly equal analysis_set")
    output = args.out.resolve()
    grid_scoring = output / "00_full_grid_scoring"
    run_scorer(
        args, paths, args.rds, args.partitions, grid_scoring,
        grid_only=True, analysis_membership=paths["analysis_set"],
    )
    evidence_out = output / "01_resolution_evidence"
    run([
        sys.executable, str(paths["build_resolution_grid_evidence.py"]),
        "--scoring-output", str(grid_scoring),
        "--catalog", str(paths["catalog"]),
        "--selection-purpose", "whole_tissue_cohort_partition",
        "--out", str(evidence_out),
    ], output / "logs/01_resolution_evidence.log")
    selector_out = output / "02_resolution_selection"
    run([
        sys.executable, str(paths["select_lineage_resolution.py"]),
        "--grid-evidence", str(evidence_out / "resolution_grid_evidence.tsv"),
        "--selection-purpose", "whole_tissue_cohort_partition",
        "--out", str(selector_out),
    ], output / "logs/02_resolution_selection.log")
    selection = json.loads((selector_out / "resolution_selection.json").read_text(encoding="utf-8"))
    selected_partitions = selector_out / "selected_neighbor_partitions.tsv.gz"
    materialize_selected_neighbors(args.partitions, selection, selected_partitions)
    selected_ids = {
        row.get("cell_id", "") for row in read_tsv(selected_partitions)
        if row.get("resolution_role") == "selected"
    }
    if selected_ids != analysis_ids:
        raise RuntimeError("selected whole-tissue partition differs from analysis_set")
    scoring = output / "03_provisional_evidence"
    selected_evidence, selected_evidence_manifest = materialize_selected_grid_evidence(
        grid_scoring, selection, scoring
    )
    plan_out = output / "04_cohort_plan"
    run([
        sys.executable, str(paths["build_whole_tissue_cohort_plan.py"]),
        "--partitions", str(selected_partitions),
        "--cluster-evidence", str(selected_evidence),
        "--catalog", str(paths["catalog"]), "--out", str(plan_out),
    ], output / "logs/04_cohort_plan.log")
    plan_manifest = json.loads(
        (plan_out / "whole_tissue_cohort_plan_manifest.json").read_text(encoding="utf-8")
    )
    return {
        "status": "PASS", "phase": "whole_tissue_partition",
        "formal_membership_written": False, "release_authority_written": False,
        "full_grid_scoring": artifact(grid_scoring / "observation_scoring_manifest.json"),
        "resolution_evidence": artifact(evidence_out / "resolution_grid_evidence_manifest.json"),
        "selection": artifact(selector_out / "resolution_selection.json"),
        "provisional_evidence": artifact(selected_evidence_manifest),
        "cohort_plan": plan_manifest["cohort_plan"],
        "cluster_membership": plan_manifest["cluster_membership"],
    }


def materialize_underpowered_cohort(
    args, contract: dict, paths: dict[str, Path], membership_ids: set[str]
) -> dict:
    """Record a mathematically non-evaluable tiny cohort without label inheritance."""
    output = args.out.resolve()
    underpowered = output / "00_underpowered_not_evaluable"
    underpowered.mkdir(parents=True, exist_ok=True)
    catalog_doc = json.loads(paths["catalog"].read_text(encoding="utf-8"))
    candidates = sorted(
        catalog_doc.get("candidate_boundaries", []),
        key=lambda row: str(row.get("candidate_id", "")),
    )
    evidence_path = underpowered / "cluster_candidate_multichannel_evidence.tsv.gz"
    write_tsv(
        evidence_path,
        [
            {
                "resolution": "", "resolution_role": "underpowered_not_evaluable",
                "source_boundary": args.cohort_id,
                "source_cluster": "underpowered",
                "candidate_id": str(candidate.get("candidate_id", "")),
                "evaluation_status": "underpowered_not_evaluable",
                "n_observations": len(membership_ids),
            }
            for candidate in candidates
        ],
        [
            "resolution", "resolution_role", "source_boundary",
            "source_cluster", "candidate_id", "evaluation_status",
            "n_observations",
        ],
    )
    ancestry_path = underpowered / "raw_count_ancestry.json"
    ancestry = {
        "schema_version": "2.2", "status": "UNDERPOWERED_NOT_EVALUABLE",
        "source_runtime_snapshot": {
            "path": str(paths["selected_input"]),
            "sha256": contract["selected_input_snapshot"]["sha256"],
            "artifact_role": "runtime_input",
        },
        "query_membership": artifact(args.membership),
        "raw_count_assay": "not_evaluated_underpowered",
        "clustering_path": "not_run_fewer_than_three_observations",
    }
    ancestry_path.write_text(
        json.dumps(ancestry, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    scoring_path = underpowered / "underpowered_scoring_manifest.json"
    scoring = {
        "status": "UNDERPOWERED_NOT_EVALUABLE",
        "schema_version": "2.2", "controller_version": "2.2.0",
        "stage": "label_blind_observation_lineage_scoring",
        "scorer": artifact(paths["run_observation_lineage_scoring.R"]),
        "historical_labels_read": False,
        "candidate_universe": [
            str(candidate.get("candidate_id", "")) for candidate in candidates
        ],
        "n_observations": len(membership_ids),
        "cluster_multichannel_evidence": artifact(evidence_path),
    }
    scoring_path.write_text(
        json.dumps(scoring, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    candidate_path = underpowered / "base_candidate_membership.tsv.gz"
    write_tsv(
        candidate_path,
        [
            {
                "cell_id": cell, "source_boundary": args.cohort_id,
                "source_cluster": "underpowered", "candidate_id": "",
                "proposed_state": "unresolved_biological",
                "proposed_broad_label": "", "confidence": "low",
                "assignment_origin": "second_round_underpowered",
                "unresolved_reason": "underpowered_not_evaluable",
            }
            for cell in sorted(membership_ids)
        ],
        [
            "cell_id", "source_boundary", "source_cluster", "candidate_id",
            "proposed_state", "proposed_broad_label", "confidence",
            "assignment_origin", "unresolved_reason",
        ],
    )
    pending_path = underpowered / "pending_local_split_membership.tsv.gz"
    write_tsv(
        pending_path, [],
        ["cell_id", "source_boundary", "source_cluster", "pending_reason"],
    )
    fine_path = underpowered / "fine_candidate_proposals.tsv"
    write_tsv(
        fine_path,
        [
            {
                "cohort_id": args.cohort_id, "subcluster_id": "underpowered",
                "candidate_id": str(candidate.get("candidate_id", "")),
                "parent_broad_label": str(candidate.get("release_broad_label", "")),
                "proposed_fine_label": str(candidate.get("release_fine_label", "")),
                "status": "not_evaluable", "release_candidate": "false",
                "lineage_supported_fraction": "",
                "strongest_competing_fraction": "",
                "contradiction_fraction": "",
                "reason": "underpowered_not_evaluable",
            }
            for candidate in candidates
            if str(candidate.get("candidate_role", "")) == "fine"
            and str(candidate.get("release_broad_label", ""))
        ],
        [
            "cohort_id", "subcluster_id", "candidate_id",
            "parent_broad_label", "proposed_fine_label", "status",
            "release_candidate", "lineage_supported_fraction",
            "strongest_competing_fraction", "contradiction_fraction", "reason",
        ],
    )
    state_path = underpowered / "state_annotation_proposals.tsv"
    write_tsv(
        state_path, [],
        [
            "cohort_id", "subcluster_id", "candidate_id", "state_annotation",
            "lineage_supported_fraction", "contradiction_fraction",
            "assignment_scope",
        ],
    )
    unmodeled_rows = underpowered / "unmodeled_lineage_candidates.tsv"
    write_tsv(unmodeled_rows, [], ["program_id", "status"])
    unmodeled_path = underpowered / "unmodeled_discovery_manifest.json"
    unmodeled_path.write_text(json.dumps({
        "status": "UNDERPOWERED_NOT_EVALUABLE", "schema_version": "2.2",
        "accepted_program_n": 0,
        "accepted_programs": artifact(unmodeled_rows),
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    adjudication_path = underpowered / "second_round_adjudication_manifest.json"
    adjudication = {
        "status": "UNDERPOWERED_NOT_EVALUABLE", "schema_version": "2.2",
        "stage": "cluster_cohort_recluster", "cohort_id": args.cohort_id,
        "source_initial_cluster": args.source_initial_cluster,
        "formal_membership_written": False, "full_catalog_scan": True,
        "all_candidates_not_evaluable": True,
        "provisional_broad_visible_during_scoring": False,
        "n_observations": len(membership_ids), "n_subclusters": 0,
        "n_pending_local_split": 0,
        "base_candidate_membership": artifact(candidate_path),
        "pending_local_split_membership": artifact(pending_path),
    }
    adjudication_path.write_text(
        json.dumps(adjudication, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    outcome = {
        "schema_version": "2.2", "cohort_id": args.cohort_id,
        "cohort_type": "initial_cluster_recluster",
        "question_mode": "open_world_identity",
        "source_initial_cluster": args.source_initial_cluster,
        "provisional_broad_after_score_freeze": args.provisional_broad,
        "raw_count_assay": "not_evaluated_underpowered",
        "raw_count_ancestry": artifact(ancestry_path),
        "full_catalog_scan": True, "formal_membership_written": False,
        "whole_tissue_manifest": artifact(args.whole_manifest),
        "context_evidence": (
            artifact(args.context_evidence) if args.context_evidence else None
        ),
        "query_membership": {
            "path": str(args.membership.resolve()),
            "sha256": sha256(args.membership),
            "n_observations": len(membership_ids),
        },
        "candidate_grid": [
            float(value) for value in
            contract["query_reclustering"]["candidate_resolutions"]
        ],
        "second_round_adjudication": artifact(adjudication_path),
        "selected_scoring": artifact(scoring_path),
        "selected_cluster_evidence": artifact(evidence_path),
        "fine_candidate_proposals": artifact(fine_path),
        "state_annotation_proposals": artifact(state_path),
        "unmodeled_discovery": artifact(unmodeled_path),
        "local_split_required": False, "n_pending_local_split": 0,
        "terminal_outcome": "underpowered_not_evaluable",
    }
    outcome_path = output / "cohort_outcome.json"
    outcome_path.write_text(
        json.dumps(outcome, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return {
        "status": "PASS", "cohort_status": "UNDERPOWERED_NOT_EVALUABLE",
        "phase": "cluster_cohort_recluster", "cohort_id": args.cohort_id,
        "source_initial_cluster": args.source_initial_cluster,
        "formal_membership_written": False,
        "whole_tissue_manifest": artifact(args.whole_manifest),
        "raw_count_reclustering": "not_run_fewer_than_three_observations",
        "raw_count_assay": "not_evaluated_underpowered",
        "raw_count_ancestry": artifact(ancestry_path),
        "selected_scoring": artifact(scoring_path),
        "selected_cluster_evidence": artifact(evidence_path),
        "adjudication": artifact(adjudication_path),
        "base_candidate_membership": artifact(candidate_path),
        "pending_local_split_membership": artifact(pending_path),
        "n_pending_local_split": 0,
        "unmodeled": artifact(unmodeled_path),
        "fine_candidate_proposals": artifact(fine_path),
        "state_annotation_proposals": artifact(state_path),
        "cohort_outcome": artifact(outcome_path),
        "query_membership": outcome["query_membership"],
        "context_evidence": outcome["context_evidence"],
    }


def phase_cohort(args, contract: dict, paths: dict[str, Path]) -> dict:
    validate_selected_input(args.rds, paths)
    validate_blinded_table(args.membership)
    whole, plan_row = validate_cohort_plan_binding(args)
    if not args.context_evidence:
        context_record = contract.get("candidate_context_evidence") or {}
        context_path = Path(str(context_record.get("path", "")))
        if context_path.is_file():
            if context_record.get("sha256") != sha256(context_path):
                raise RuntimeError("contract-bound context evidence is stale")
            args.context_evidence = context_path
    if args.context_evidence:
        validate_blinded_table(args.context_evidence)
    membership_ids = {row.get("cell_id", "") for row in read_tsv(args.membership)}
    analysis_ids = {row.get("cell_id", "") for row in read_tsv(paths["analysis_set"])}
    if not membership_ids or not membership_ids <= analysis_ids:
        raise RuntimeError("cohort membership is empty or outside analysis_set")
    if len(membership_ids) < 3:
        return materialize_underpowered_cohort(
            args, contract, paths, membership_ids
        )
    output = args.out.resolve()
    grid = [float(value) for value in contract["query_reclustering"]["candidate_resolutions"]]
    workers = min(len(grid), args.resolution_workers)
    cache_fingerprint, cache_payload = recluster_cache_fingerprint(
        contract, paths, args, grid
    )
    cache_reused = bool(args.reuse_recluster_manifest)
    if cache_reused:
        query_rds, runner_manifest_path, zero_count_path, grid_partitions = (
            validate_recluster_cache(
                args.reuse_recluster_manifest, cache_fingerprint, contract
            )
        )
        runner_out = runner_manifest_path.parent
    else:
        runner_out = output / "00_recluster"
        run([
            args.rscript, str(paths["run_seurat_cohort_recluster.R"]),
            "--rds", str(args.rds.resolve()), "--membership", str(args.membership.resolve()),
            "--out", str(runner_out), "--cell-id-col", "cell_id",
            "--resolutions", ",".join(str(value) for value in grid),
            "--resolution-contract", args.resolution_contract,
            "--normalization", "SCT", "--sct-method", "glmGamPoi",
            "--seed", str(args.seed), "--resolution-workers", str(workers),
        ], output / "logs/00_recluster.log")
        query_rds = runner_out / "cohort_reclustered_query_seurat.rds"
        runner_manifest_path = runner_out / "run_manifest.tsv"
        zero_count_path = runner_out / "tables/zero_count_observations.tsv"
        grid_partitions = build_cohort_grid(runner_out, args.cohort_id, grid)
        resolution_memberships = [
            runner_out / "tables" / (
                "framework_res" + str(resolution).replace(".", "p")
                + "_clusters.tsv"
            )
            for resolution in grid
        ]
    cache_manifest_path = (
        args.reuse_recluster_manifest.resolve() if cache_reused
        else (runner_out / "recluster_partition_cache_manifest.json").resolve()
    )
    zero_count = read_tsv(zero_count_path)
    if zero_count:
        raise RuntimeError(
            "cohort contains raw-count-zero observations; move them to "
            "excluded_initial_qc and refreeze the analysis_set"
        )
    parameters = runner_parameters(runner_manifest_path)
    raw_assay = parameters.get("full_feature_deg_assay", "")
    if not raw_assay or raw_assay == "SCT" or parameters.get("normalization") != "SCT":
        raise RuntimeError("cohort did not restart SCT from a non-SCT raw-count assay")
    resolution_seeds = [
        value for value in parameters.get("resolution_seeds", "").split(",")
        if value
    ]
    if (
        parameters.get("seed") != str(args.seed)
        or len(resolution_seeds) != len(grid)
        or len(set(resolution_seeds)) != len(grid)
    ):
        raise RuntimeError("cohort runner did not record deterministic per-resolution seeds")
    if not cache_reused:
        write_recluster_cache_manifest(
            cache_manifest_path, cache_fingerprint, cache_payload, query_rds,
            runner_manifest_path, zero_count_path, grid_partitions,
            resolution_memberships,
        )
    ancestry_root = (
        output / "00_recluster_cache_binding" if cache_reused else runner_out
    )
    ancestry_root.mkdir(parents=True, exist_ok=True)
    ancestry_path = ancestry_root / "raw_count_ancestry.json"
    ancestry = {
        "schema_version": "2.2",
        "status": "PASS",
        "source_runtime_snapshot": {
            "path": str(paths["selected_input"]),
            "sha256": contract["selected_input_snapshot"]["sha256"],
            "artifact_role": "runtime_input",
        },
        "query_membership": artifact(args.membership),
        "raw_count_assay": raw_assay,
        "sct_method": parameters.get("sct_method", ""),
        "sct_vst_flavor": parameters.get("sct_vst_flavor", ""),
        "clustering_path": "raw_counts_SCTv2_PCA_SNN_Leiden",
        "partition_cache_reused": cache_reused,
        "partition_cache_manifest": (
            str(args.reuse_recluster_manifest.resolve())
            if cache_reused else str(
                (runner_out / "recluster_partition_cache_manifest.json").resolve()
            )
        ),
    }
    ancestry_path.write_text(
        json.dumps(ancestry, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    grid_scoring = output / "01_full_grid_scoring"
    run_scorer(
        args, paths, query_rds, grid_partitions, grid_scoring,
        grid_only=True, assay=raw_assay,
    )
    evidence_out = output / "02_resolution_evidence"
    run([
        sys.executable, str(paths["build_resolution_grid_evidence.py"]),
        "--scoring-output", str(grid_scoring), "--catalog", str(paths["catalog"]),
        "--selection-purpose", "cohort_identity_resolution",
        "--out", str(evidence_out),
    ], output / "logs/02_resolution_evidence.log")
    selection_out = output / "03_resolution_selection"
    run([
        sys.executable, str(paths["select_lineage_resolution.py"]),
        "--grid-evidence", str(evidence_out / "resolution_grid_evidence.tsv"),
        "--selection-purpose", "cohort_identity_resolution",
        "--out", str(selection_out),
    ], output / "logs/03_resolution_selection.log")
    selection = json.loads(
        (selection_out / "resolution_selection.json").read_text(encoding="utf-8")
    )
    selected_partitions = selection_out / "selected_neighbor_partitions.tsv.gz"
    materialize_selected_neighbors(grid_partitions, selection, selected_partitions)
    scoring = output / "04_selected_scoring"
    run_scorer(args, paths, query_rds, selected_partitions, scoring, assay=raw_assay)
    unmodeled = output / "05_unmodeled"
    run([
        sys.executable, str(paths["discover_unmodeled_lineages.py"]),
        "--programs", str(scoring / "tables/resolution_deg_coexpression_programs.tsv.gz"),
        "--out", str(unmodeled),
    ], output / "logs/05_unmodeled.log")
    adjudication = output / "06_adjudication"
    adjudication_command = [
        sys.executable, str(paths["adjudicate_second_round_subclusters.py"]),
        "--partitions", str(selected_partitions),
        "--cluster-evidence", str(scoring / "tables/cluster_candidate_multichannel_evidence.tsv.gz"),
        "--scores", str(scoring / "tables/observation_lineage_scores.tsv.gz"),
        "--catalog", str(paths["catalog"]), "--contract", str(args.contract),
        "--cohort-id", args.cohort_id,
        "--source-initial-cluster", args.source_initial_cluster,
        "--provisional-status", args.provisional_status,
        "--provisional-broad", args.provisional_broad,
        "--out", str(adjudication),
    ]
    if args.context_evidence:
        adjudication_command.extend(
            ["--context-evidence", str(args.context_evidence.resolve())]
        )
    run(adjudication_command, output / "logs/06_adjudication.log")
    adjudication_manifest = json.loads(
        (adjudication / "second_round_adjudication_manifest.json").read_text(encoding="utf-8")
    )
    resolution_records = []
    for resolution in grid:
        tag = str(resolution).replace(".", "p")
        membership_path = runner_out / "tables" / f"framework_res{tag}_clusters.tsv"
        membership_rows = read_tsv(membership_path)
        resolution_records.append({
            "resolution": resolution,
            "membership": {
                "path": str(membership_path.resolve()),
                "sha256": sha256(membership_path),
                "n_observations": len(membership_rows),
            },
            "cluster_count": len({row.get("cluster", "") for row in membership_rows}),
            "evidence_index": artifact(
                grid_scoring / "observation_scoring_manifest.json"
            ),
        })
    selected_values = [float(value) for value in selection["selected_and_neighbors"]]
    cohort_outcome = {
        "schema_version": "2.2",
        "cohort_id": args.cohort_id,
        "cohort_type": "initial_cluster_recluster",
        "question_mode": "open_world_identity",
        "source_initial_cluster": args.source_initial_cluster,
        "provisional_broad_after_score_freeze": args.provisional_broad,
        "raw_count_assay": raw_assay,
        "raw_count_ancestry": artifact(ancestry_path),
        "partition_cache": artifact(
            cache_manifest_path, "derived_partition_cache"
        ),
        "partition_cache_reused": cache_reused,
        "full_catalog_scan": True,
        "formal_membership_written": False,
        "whole_tissue_manifest": artifact(args.whole_manifest),
        "query_membership": {
            "path": str(args.membership.resolve()),
            "sha256": sha256(args.membership),
            "n_observations": len(membership_ids),
        },
        "candidate_grid": grid,
        "resolutions": resolution_records,
        "selected_resolution": selected_values[0],
        "resolution_neighbors": selected_values[1:],
        "second_round_adjudication": artifact(
            adjudication / "second_round_adjudication_manifest.json"
        ),
        "selected_scoring": artifact(
            scoring / "observation_scoring_manifest.json"
        ),
        "selected_cluster_evidence": artifact(
            scoring / "tables/cluster_candidate_multichannel_evidence.tsv.gz"
        ),
        "fine_candidate_proposals": artifact(
            adjudication / "fine_candidate_proposals.tsv"
        ),
        "state_annotation_proposals": artifact(
            adjudication / "state_annotation_proposals.tsv"
        ),
        "unmodeled_discovery": artifact(
            unmodeled / "unmodeled_discovery_manifest.json"
        ),
        "local_split_required": bool(adjudication_manifest["n_pending_local_split"]),
        "n_pending_local_split": adjudication_manifest["n_pending_local_split"],
        "terminal_outcome": (
            "local_split_required"
            if adjudication_manifest["n_pending_local_split"]
            else "candidate_partition_complete"
        ),
        "context_evidence": (
            artifact(args.context_evidence)
            if args.context_evidence else None
        ),
    }
    cohort_outcome_path = output / "cohort_outcome.json"
    cohort_outcome_path.write_text(
        json.dumps(cohort_outcome, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return {
        "status": "PASS",
        "cohort_status": adjudication_manifest["status"],
        "phase": "cluster_cohort_recluster", "cohort_id": args.cohort_id,
        "source_initial_cluster": args.source_initial_cluster,
        "formal_membership_written": False,
        "whole_tissue_manifest": artifact(args.whole_manifest),
        "raw_count_reclustering": "SCTv2_PCA_SNN_Leiden",
        "raw_count_assay": raw_assay,
        "raw_count_ancestry": artifact(ancestry_path),
        "partition_cache": artifact(
            cache_manifest_path, "derived_partition_cache"
        ),
        "partition_cache_reused": cache_reused,
        "selection": artifact(selection_out / "resolution_selection.json"),
        "selected_scoring": artifact(scoring / "observation_scoring_manifest.json"),
        "selected_cluster_evidence": artifact(
            scoring / "tables/cluster_candidate_multichannel_evidence.tsv.gz"
        ),
        "adjudication": artifact(adjudication / "second_round_adjudication_manifest.json"),
        "base_candidate_membership": adjudication_manifest["base_candidate_membership"],
        "pending_local_split_membership": adjudication_manifest["pending_local_split_membership"],
        "n_pending_local_split": adjudication_manifest["n_pending_local_split"],
        "unmodeled": artifact(unmodeled / "unmodeled_discovery_manifest.json"),
        "fine_candidate_proposals": artifact(
            adjudication / "fine_candidate_proposals.tsv"
        ),
        "state_annotation_proposals": artifact(
            adjudication / "state_annotation_proposals.tsv"
        ),
        "cohort_outcome": artifact(cohort_outcome_path),
        "query_membership": {
            "path": str(args.membership.resolve()),
            "sha256": sha256(args.membership),
            "n_observations": len(membership_ids),
        },
        "context_evidence": (
            artifact(args.context_evidence)
            if args.context_evidence else None
        ),
    }


def phase_local(args, contract: dict, paths: dict[str, Path]) -> dict:
    output = args.out.resolve()
    scores = args.scoring_output / "tables/observation_lineage_scores.tsv.gz"
    evidence = args.scoring_output / "tables/cluster_candidate_multichannel_evidence.tsv.gz"
    trigger = json.loads(args.trigger_manifest.read_text(encoding="utf-8"))
    if not controller_manifest_bound(
        trigger, args.contract, "cluster_cohort_recluster"
    ):
        raise RuntimeError("local split trigger is not a bound second-round controller manifest")
    if trigger.get("cohort_id") != args.source_boundary:
        raise RuntimeError("local split source boundary differs from trigger cohort")
    selected_scoring = trigger.get("selected_scoring", {})
    expected_scoring_manifest = args.scoring_output / "observation_scoring_manifest.json"
    if (
        Path(str(selected_scoring.get("path", ""))).resolve()
        != expected_scoring_manifest.resolve()
        or selected_scoring.get("sha256") != sha256(expected_scoring_manifest)
    ):
        raise RuntimeError("local split scoring output differs from its trigger cohort")
    adjudication_path = Path(str(trigger.get("adjudication", {}).get("path", "")))
    if (
        not adjudication_path.is_file()
        or trigger.get("adjudication", {}).get("sha256") != sha256(adjudication_path)
    ):
        raise RuntimeError("local split trigger adjudication is missing or stale")
    adjudication = json.loads(adjudication_path.read_text(encoding="utf-8"))
    pending_record = adjudication.get("pending_local_split_membership", {})
    pending_path = Path(str(pending_record.get("path", "")))
    if (
        not pending_path.is_file()
        or pending_record.get("sha256") != sha256(pending_path)
    ):
        raise RuntimeError("local split trigger membership is missing or stale")
    pending_ids = {
        str(row.get("cell_id", ""))
        for row in read_tsv(pending_path)
        if row.get("source_boundary") == args.source_boundary
        and row.get("source_cluster") == args.source_cluster
    }
    score_ids = {
        str(row.get("cell_id", ""))
        for row in read_tsv(scores)
        if row.get("source_boundary") == args.source_boundary
        and row.get("source_cluster") == args.source_cluster
    }
    if not pending_ids or pending_ids != score_ids:
        raise RuntimeError("local split source is not one exact triggered mixed subcluster")
    context_path = None
    context_record = trigger.get("context_evidence") or {}
    if context_record:
        candidate = Path(str(context_record.get("path", "")))
        if (
            not candidate.is_file()
            or context_record.get("sha256") != sha256(candidate)
        ):
            raise RuntimeError("local split trigger context evidence is missing or stale")
        validate_blinded_table(candidate)
        context_path = candidate
    authority = stage_authority(
        "local_mixed_subcluster_split", args.contract, paths, output,
        cluster_evidence=evidence, trigger_manifest=args.trigger_manifest,
        trigger_membership=pending_path,
        context_evidence=context_path,
    )
    subset_out = output / "00_candidate_subsets"
    subset_command = [
        args.rscript, str(paths["derive_candidate_local_subsets.R"]),
        "--scores", str(scores), "--cluster-evidence", str(evidence),
        "--catalog", str(paths["catalog"]), "--release-level", "broad",
        "--source-boundary", args.source_boundary,
        "--source-cluster", args.source_cluster,
        "--workers", str(args.subset_workers), "--out", str(subset_out),
    ]
    if context_path:
        subset_command.extend(["--context-evidence", str(context_path)])
    run(subset_command, output / "logs/00_candidate_subsets.log")
    closure = output / "01_local_remainder"
    closure_command = [
        sys.executable, str(paths["close_exact_remainders.py"]),
        "--scores", str(scores), "--cluster-evidence", str(evidence),
        "--subset-membership", str(subset_out / "candidate_subset_membership.tsv.gz"),
        "--subset-evidence", str(subset_out / "candidate_subset_evidence.tsv"),
        "--catalog", str(paths["catalog"]), "--contract", str(args.contract.resolve()),
        "--stage-authority", str(authority), "--scope", "local_mixed_subcluster",
        "--source-boundary", args.source_boundary,
        "--source-cluster", args.source_cluster,
        "--release-level", "broad", "--out", str(closure),
    ]
    if context_path:
        closure_command.extend(["--context-evidence", str(context_path)])
    run(closure_command, output / "logs/01_local_remainder.log")
    manifest = json.loads(
        (closure / "exact_remainder_closure_manifest.json").read_text(encoding="utf-8")
    )
    return {
        "status": "PASS", "phase": "local_mixed_subcluster_split",
        "formal_membership_written": False,
        "source_boundary": args.source_boundary,
        "source_cluster": args.source_cluster,
        "trigger_manifest": artifact(args.trigger_manifest),
        "trigger_membership": artifact(pending_path),
        "candidate_membership": manifest["candidate_membership"],
        "local_subset_validation": artifact(
            closure / "subset_validation.tsv"
        ),
        "local_remainder_audit": artifact(
            closure / "exact_remainder_audit.tsv"
        ),
        "unresolved_biological_n": manifest["unresolved_biological_n"],
    }


def phase_merge(args, contract: dict, paths: dict[str, Path]) -> dict:
    output = args.out.resolve()
    if (
        args.analysis_membership.resolve() != paths["analysis_set"]
        or sha256(args.analysis_membership) != sha256(paths["analysis_set"])
    ):
        raise RuntimeError("merge analysis membership differs from frozen analysis_set")
    whole = json.loads(args.whole_manifest.read_text(encoding="utf-8"))
    if not controller_manifest_bound(whole, args.contract, "whole_tissue_partition"):
        raise RuntimeError("broad merge requires the bound whole-tissue cohort plan")
    plan_record = whole.get("cohort_plan", {})
    plan_path = Path(str(plan_record.get("path", "")))
    if not plan_path.is_file() or plan_record.get("sha256") != sha256(plan_path):
        raise RuntimeError("whole-tissue cohort plan is missing or stale")
    plan_rows = read_tsv(plan_path)
    plan_by_cohort = {
        str(row.get("cohort_id", "")): row for row in plan_rows
    }
    planned = set(plan_by_cohort)
    if not planned or "" in planned or len(planned) != len(plan_rows):
        raise RuntimeError("whole-tissue cohort plan has invalid cohort IDs")

    candidate_paths: list[Path] = []
    source_manifests: list[Path] = []
    pending_groups: set[tuple[str, str]] = set()
    seen_cohorts: set[str] = set()
    cluster_evidence_paths: list[Path] = []
    fine_proposal_paths: list[Path] = []
    state_proposal_paths: list[Path] = []
    unmodeled_manifest_paths: list[Path] = []
    observation_score_paths: list[Path] = []
    local_subset_validation_paths: list[Path] = []
    local_remainder_audit_paths: list[Path] = []
    for source_path in args.cohort_manifest:
        manifest = json.loads(source_path.read_text(encoding="utf-8"))
        if not controller_manifest_bound(
            manifest, args.contract, "cluster_cohort_recluster"
        ):
            raise RuntimeError(f"invalid second-round cohort source: {source_path}")
        cohort_id = str(manifest.get("cohort_id", ""))
        if not cohort_id or cohort_id in seen_cohorts:
            raise RuntimeError("duplicate or empty second-round cohort source")
        if cohort_id not in plan_by_cohort:
            raise RuntimeError("second-round cohort is absent from the whole-tissue plan")
        seen_cohorts.add(cohort_id)
        plan_row = plan_by_cohort[cohort_id]
        whole_record = manifest.get("whole_tissue_manifest", {})
        if (
            Path(str(whole_record.get("path", ""))).resolve()
            != args.whole_manifest.resolve()
            or whole_record.get("sha256") != sha256(args.whole_manifest)
        ):
            raise RuntimeError(
                f"{cohort_id}: second-round cohort is not bound to this whole-tissue plan"
            )
        query_record = manifest.get("query_membership", {})
        if (
            str(manifest.get("source_initial_cluster", ""))
            != str(plan_row.get("source_initial_cluster", ""))
            or query_record.get("sha256") != plan_row.get("membership_sha256")
            or int(query_record.get("n_observations", -1))
            != int(plan_row.get("n_observations", -2))
        ):
            raise RuntimeError(
                f"{cohort_id}: second-round query differs from the exact initial-cluster cohort"
            )
        candidate_record = manifest.get("base_candidate_membership", {})
        candidate_path = Path(str(candidate_record.get("path", "")))
        if (
            not candidate_path.is_file()
            or candidate_record.get("sha256") != sha256(candidate_path)
        ):
            raise RuntimeError(f"{cohort_id}: base candidate membership is stale")
        candidate_paths.append(candidate_path)
        source_manifests.append(source_path)
        for key, destination in (
            ("selected_cluster_evidence", cluster_evidence_paths),
            ("fine_candidate_proposals", fine_proposal_paths),
            ("state_annotation_proposals", state_proposal_paths),
            ("unmodeled", unmodeled_manifest_paths),
        ):
            record = manifest.get(key, {})
            artifact_path = Path(str(record.get("path", "")))
            if (
                not artifact_path.is_file()
                or record.get("sha256") != sha256(artifact_path)
            ):
                raise RuntimeError(f"{cohort_id}: {key} is missing or stale")
            destination.append(artifact_path)
        scoring_record = manifest.get("selected_scoring", {})
        scoring_path = Path(str(scoring_record.get("path", "")))
        if (
            not scoring_path.is_file()
            or scoring_record.get("sha256") != sha256(scoring_path)
        ):
            raise RuntimeError(f"{cohort_id}: selected scoring manifest is missing or stale")
        scoring_manifest = json.loads(scoring_path.read_text(encoding="utf-8"))
        score_record = scoring_manifest.get("outputs", {}).get("observation_scores")
        if manifest.get("cohort_status") == "UNDERPOWERED_NOT_EVALUABLE":
            if score_record:
                raise RuntimeError(
                    f"{cohort_id}: underpowered cohort unexpectedly exposes observation scores"
                )
        else:
            score_path = Path(
                str(score_record.get("path", ""))
                if isinstance(score_record, dict) else str(score_record or "")
            )
            if not score_path.is_file():
                raise RuntimeError(f"{cohort_id}: observation scores are missing")
            if (
                isinstance(score_record, dict)
                and score_record.get("sha256")
                and score_record.get("sha256") != sha256(score_path)
            ):
                raise RuntimeError(f"{cohort_id}: observation scores are stale")
            observation_score_paths.append(score_path)
        pending_record = manifest.get("pending_local_split_membership", {})
        pending_path = Path(str(pending_record.get("path", "")))
        if (
            not pending_path.is_file()
            or pending_record.get("sha256") != sha256(pending_path)
        ):
            raise RuntimeError(f"{cohort_id}: pending local-split membership is stale")
        pending_groups.update(
            (
                str(row.get("source_boundary", "")),
                str(row.get("source_cluster", "")),
            )
            for row in read_tsv(pending_path)
        )
    if seen_cohorts != planned:
        raise RuntimeError("every initial cluster must have exactly one second-round cohort before merge")
    if any(not boundary or not cluster for boundary, cluster in pending_groups):
        raise RuntimeError("pending local-split registry contains an invalid source group")

    local_groups: set[tuple[str, str]] = set()
    for source_path in args.local_manifest:
        manifest = json.loads(source_path.read_text(encoding="utf-8"))
        if not controller_manifest_bound(
            manifest, args.contract, "local_mixed_subcluster_split"
        ):
            raise RuntimeError(f"invalid local mixed-subcluster source: {source_path}")
        key = (
            str(manifest.get("source_boundary", "")),
            str(manifest.get("source_cluster", "")),
        )
        if not all(key) or key in local_groups:
            raise RuntimeError("duplicate or empty local mixed-subcluster source")
        local_groups.add(key)
        candidate_record = manifest.get("candidate_membership", {})
        candidate_path = Path(str(candidate_record.get("path", "")))
        if (
            not candidate_path.is_file()
            or candidate_record.get("sha256") != sha256(candidate_path)
        ):
            raise RuntimeError(f"{key}: local candidate membership is stale")
        candidate_paths.append(candidate_path)
        source_manifests.append(source_path)
        for artifact_key, destination in (
            ("local_subset_validation", local_subset_validation_paths),
            ("local_remainder_audit", local_remainder_audit_paths),
        ):
            record = manifest.get(artifact_key, {})
            artifact_path = Path(str(record.get("path", "")))
            if (
                not artifact_path.is_file()
                or record.get("sha256") != sha256(artifact_path)
            ):
                raise RuntimeError(
                    f"{artifact_key} is missing or stale: {source_path}"
                )
            destination.append(artifact_path)
    if local_groups != pending_groups:
        raise RuntimeError("local split sources do not exactly close all mixed-subcluster triggers")

    authority = stage_authority(
        "merge_and_freeze_broad", args.contract, paths, output,
        whole_manifest=args.whole_manifest,
        candidate_source_manifests=source_manifests,
    )
    command = [
        sys.executable, str(paths["merge_and_freeze_broad_membership.py"]),
        "--contract", str(args.contract.resolve()),
        "--stage-authority", str(authority),
        "--analysis-membership", str(args.analysis_membership.resolve()),
        "--out", str(output / "00_broad_freeze"),
    ]
    for membership_path, source_path in zip(candidate_paths, source_manifests):
        command.extend(["--candidate-membership", str(membership_path.resolve())])
        command.extend(["--candidate-source-manifest", str(source_path.resolve())])
    run(command, output / "logs/00_broad_freeze.log")
    manifest_path = output / "00_broad_freeze/broad_freeze_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    return {
        "status": "PASS", "phase": "merge_and_freeze_broad",
        "formal_broad_membership_written": True,
        "broad_freeze": artifact(manifest_path),
        "membership": manifest["membership"],
        "cluster_evidence": [artifact(path) for path in cluster_evidence_paths],
        "fine_candidate_proposals": [artifact(path) for path in fine_proposal_paths],
        "state_annotation_proposals": [artifact(path) for path in state_proposal_paths],
        "unmodeled_manifests": [artifact(path) for path in unmodeled_manifest_paths],
        "observation_scores": [artifact(path) for path in observation_score_paths],
        "local_subset_validations": [
            artifact(path) for path in local_subset_validation_paths
        ],
        "local_remainder_audits": [
            artifact(path) for path in local_remainder_audit_paths
        ],
    }


def validate_prerequisite(args, required_stage: str) -> dict:
    prerequisite = json.loads(args.prerequisite_manifest.read_text(encoding="utf-8"))
    if not controller_manifest_bound(prerequisite, args.contract, required_stage):
        raise RuntimeError(
            f"{args.phase} requires a PASS {required_stage} manifest"
        )
    return prerequisite


def controller_manifest_bound(manifest: dict, contract_path: Path, phase: str) -> bool:
    return (
        manifest.get("status") == "PASS"
        and manifest.get("controller_version") == "2.2.0"
        and manifest.get("phase") == phase
        and manifest.get("annotation_contract", {}).get("sha256")
        == sha256(contract_path)
    )


def run_targeted_follicle_roi_iteration(
    args, contract: dict, paths: dict[str, Path], output: Path,
    post_membership: Path, quality_path: Path,
    observation_score_paths: list[Path],
) -> dict | None:
    """Run one bounded raw-count follicle-wall iteration when all failures are ROI-local."""
    quality = json.loads(quality_path.read_text(encoding="utf-8"))
    if quality.get("status") != "ITERATION_REQUIRED":
        return None
    action_record = quality.get("required_next_actions", {})
    action_path = Path(str(action_record.get("path", "")))
    if (
        not action_path.is_file()
        or action_record.get("sha256") != sha256(action_path)
    ):
        raise RuntimeError("biological-quality next actions are missing or stale")
    actions = read_tsv(action_path)
    target_rows = [
        row for row in actions
        if row.get("endpoint") == "follicle_roi_histology"
        and re.fullmatch(r"F\d+", str(row.get("scope_id", "")))
    ]
    if not target_rows or len(target_rows) != len(actions):
        return None
    target_rois = sorted({str(row["scope_id"]) for row in target_rows})
    follicle = quality.get("quality_endpoints", {}).get(
        "follicle_roi_histology", {}
    )
    roi_record = follicle.get("roi_membership", {})
    roi_path = Path(str(roi_record.get("path", "")))
    expected_record = follicle.get("roi_review", {})
    expected_roi_path = Path(str(expected_record.get("path", "")))
    if (
        not roi_path.is_file() or roi_record.get("sha256") != sha256(roi_path)
        or not expected_roi_path.is_file()
        or expected_record.get("sha256") != sha256(expected_roi_path)
    ):
        raise RuntimeError("follicle ROI membership/review is missing or stale")
    roi_rows = read_tsv(roi_path)
    output_root = output / "05_follicle_roi_repair"
    grid = [
        float(value)
        for value in contract["query_reclustering"]["candidate_resolutions"]
    ]
    workers = min(len(grid), args.resolution_workers)
    repair_scores: dict[str, Path] = {}
    ancestry_paths: dict[str, Path] = {}
    for roi_id in target_rois:
        roi_output = output_root / "roi_cohorts" / roi_id
        membership_path = roi_output / "roi_membership.tsv.gz"
        members = sorted({
            str(row.get("cell_id", "")) for row in roi_rows
            if row.get("follicle_roi_id") == roi_id and row.get("cell_id")
        })
        if len(members) < 3:
            raise RuntimeError(f"{roi_id} is too small for a raw-count ROI recluster")
        write_tsv(
            membership_path,
            [{"cell_id": cell} for cell in members],
            ["cell_id"],
        )
        runner_out = roi_output / "00_recluster"
        run([
            args.rscript, str(paths["run_seurat_cohort_recluster.R"]),
            "--rds", str(paths["selected_input"]),
            "--membership", str(membership_path),
            "--out", str(runner_out), "--cell-id-col", "cell_id",
            "--resolutions", ",".join(str(value) for value in grid),
            "--resolution-contract", "sheep_ovary",
            "--normalization", "SCT", "--sct-method", "glmGamPoi",
            "--seed", str(args.seed), "--resolution-workers", str(workers),
        ], roi_output / "logs/00_recluster.log")
        if read_tsv(runner_out / "tables/zero_count_observations.tsv"):
            raise RuntimeError(
                f"{roi_id} contains raw-count-zero observations; ROI repair cannot relabel them"
            )
        parameters = runner_parameters(runner_out / "run_manifest.tsv")
        raw_assay = parameters.get("full_feature_deg_assay", "")
        if (
            not raw_assay or raw_assay == "SCT"
            or parameters.get("normalization") != "SCT"
        ):
            raise RuntimeError(
                f"{roi_id} did not restart from a project-local non-SCT raw-count assay"
            )
        ancestry_path = runner_out / "raw_count_ancestry.json"
        ancestry_path.write_text(json.dumps({
            "schema_version": "2.2", "status": "PASS",
            "source_runtime_snapshot": {
                "path": str(paths["selected_input"]),
                "sha256": contract["selected_input_snapshot"]["sha256"],
                "artifact_role": "runtime_input",
            },
            "query_membership": artifact(membership_path),
            "raw_count_assay": raw_assay,
            "sct_method": parameters.get("sct_method", ""),
            "sct_vst_flavor": parameters.get("sct_vst_flavor", ""),
            "clustering_path": "raw_counts_SCTv2_PCA_SNN_Leiden_follicle_ROI",
        }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        grid_partitions = build_cohort_grid(
            runner_out, f"follicle_roi_{roi_id}", grid
        )
        query_rds = runner_out / "cohort_reclustered_query_seurat.rds"
        grid_scoring = roi_output / "01_full_grid_scoring"
        run_scorer(
            args, paths, query_rds, grid_partitions, grid_scoring,
            grid_only=True, assay=raw_assay,
        )
        evidence_out = roi_output / "02_resolution_evidence"
        run([
            sys.executable, str(paths["build_resolution_grid_evidence.py"]),
            "--scoring-output", str(grid_scoring),
            "--catalog", str(paths["catalog"]),
            "--selection-purpose", "cohort_identity_resolution",
            "--out", str(evidence_out),
        ], roi_output / "logs/02_resolution_evidence.log")
        selection_out = roi_output / "03_resolution_selection"
        run([
            sys.executable, str(paths["select_lineage_resolution.py"]),
            "--grid-evidence", str(evidence_out / "resolution_grid_evidence.tsv"),
            "--selection-purpose", "cohort_identity_resolution",
            "--out", str(selection_out),
        ], roi_output / "logs/03_resolution_selection.log")
        selection = json.loads(
            (selection_out / "resolution_selection.json").read_text(encoding="utf-8")
        )
        selected_partitions = selection_out / "selected_neighbor_partitions.tsv.gz"
        materialize_selected_neighbors(
            grid_partitions, selection, selected_partitions
        )
        scoring = roi_output / "04_selected_scoring"
        run_scorer(
            args, paths, query_rds, selected_partitions, scoring,
            assay=raw_assay,
        )
        repair_scores[roi_id] = (
            scoring / "tables/observation_lineage_scores.tsv.gz"
        )
        ancestry_paths[roi_id] = ancestry_path

    repair_authority = stage_authority(
        "atlas_and_completeness_review", args.contract, paths, output_root,
        pre_repair_membership=post_membership,
        pre_repair_biological_quality=quality_path,
        base_scores=observation_score_paths,
        repair_scores=list(repair_scores.values()),
        repair_ancestry=list(ancestry_paths.values()),
    )
    materialized = output_root / "materialized_repair"
    command = [
        sys.executable,
        str(paths["apply_sheep_ovary_follicle_roi_repair.py"]),
        "--contract", str(args.contract),
        "--stage-authority", str(repair_authority),
        "--membership", str(post_membership),
        "--quality-review", str(quality_path),
        "--out", str(materialized),
    ]
    for path in observation_score_paths:
        command.extend(["--base-scores", str(path)])
    for roi_id in target_rois:
        command.extend(["--repair-score", f"{roi_id}={repair_scores[roi_id]}"])
        command.extend(["--repair-ancestry", f"{roi_id}={ancestry_paths[roi_id]}"])
    run(command, output_root / "logs/05_materialize_repair.log")
    repair_manifest_path = materialized / "follicle_roi_repair_manifest.json"
    repair_manifest = json.loads(
        repair_manifest_path.read_text(encoding="utf-8")
    )
    repaired_membership = Path(
        repair_manifest["repaired_membership"]["path"]
    )
    combined_scores = Path(
        repair_manifest["combined_observation_scores"]["path"]
    )
    post_review = output_root / "post_repair_biological_quality"
    run([
        sys.executable,
        str(paths["validate_sheep_ovary_biological_quality.py"]),
        "--membership", str(repaired_membership),
        "--coordinate-membership", str(repaired_membership),
        "--catalog", str(paths["catalog"]),
        "--scores", str(combined_scores),
        "--expected-roi-review", str(expected_roi_path),
        "--out", str(post_review),
    ], output_root / "logs/06_post_repair_biological_quality.log", allowed_codes=(0, 2))
    post_quality_path = post_review / "sheep_ovary_biological_quality_review.json"
    post_quality = json.loads(post_quality_path.read_text(encoding="utf-8"))
    return {
        "status": post_quality.get("status", ""),
        "membership": repaired_membership,
        "scores": combined_scores,
        "quality_review": post_quality_path,
        "repair_manifest": repair_manifest_path,
        "target_rois": target_rois,
    }


def phase_atlas(args, contract: dict, paths: dict[str, Path]) -> dict:
    prerequisite = validate_prerequisite(args, "merge_and_freeze_broad")
    cluster_evidence_paths = bound_artifact_paths(
        prerequisite.get("cluster_evidence"), "second-round cluster evidence"
    )
    fine_proposal_paths = bound_artifact_paths(
        prerequisite.get("fine_candidate_proposals"), "fine proposal"
    )
    state_proposal_paths = bound_artifact_paths(
        prerequisite.get("state_annotation_proposals"), "state proposal"
    )
    unmodeled_manifest_paths = bound_artifact_paths(
        prerequisite.get("unmodeled_manifests"), "unmodeled discovery"
    )
    observation_score_paths = bound_artifact_paths(
        prerequisite.get("observation_scores"), "second-round observation score"
    )
    local_subset_validation_paths = bound_artifact_paths(
        prerequisite.get("local_subset_validations"),
        "local subset validation",
    )
    local_remainder_audit_paths = bound_artifact_paths(
        prerequisite.get("local_remainder_audits"),
        "local remainder audit",
    )
    membership_record = prerequisite.get("membership", {})
    if (
        Path(str(membership_record.get("path", ""))).resolve()
        != args.frozen_broad_membership.resolve()
        or membership_record.get("sha256") != sha256(args.frozen_broad_membership)
    ):
        raise RuntimeError("Atlas phase membership differs from broad freeze")
    output = args.out.resolve()
    authority = stage_authority(
        "atlas_and_completeness_review", args.contract, paths, output,
        frozen_broad=args.frozen_broad_membership,
        atlas_mapping=args.atlas_mapping,
        calibration_manifest=args.calibration_manifest,
        atlas_decisions=args.atlas_decisions,
        unmodeled_decisions=args.unmodeled_decisions,
        cluster_evidence=cluster_evidence_paths,
        fine_candidate_proposals=fine_proposal_paths,
        state_annotation_proposals=state_proposal_paths,
        unmodeled_manifests=unmodeled_manifest_paths,
        observation_scores=observation_score_paths,
        local_subset_validations=local_subset_validation_paths,
        local_remainder_audits=local_remainder_audit_paths,
    )
    routing_out = output / "00_atlas_routing"
    run([
        sys.executable, str(paths["route_global_atlas_v2.py"]),
        "--cell-ledger", str(args.frozen_broad_membership.resolve()),
        "--atlas-mapping", str(args.atlas_mapping.resolve()),
        "--calibration-manifest", str(args.calibration_manifest.resolve()),
        "--workflow-profile", str(paths["workflow_profile"]),
        "--out", str(routing_out),
    ], output / "logs/00_atlas_routing.log", allowed_codes=(0, 2))
    validation_path = output / "01_atlas_validation.json"
    validation_command = [
        sys.executable, str(paths["validate_global_atlas_v2.py"]),
        "--routing-manifest", str(routing_out / "atlas_state_routing_manifest.json"),
        "--out", str(validation_path),
    ]
    if args.atlas_decisions:
        validation_command.extend(["--decisions", str(args.atlas_decisions.resolve())])
    run(validation_command, output / "logs/01_atlas_validation.log")
    apply_out = output / "02_post_atlas_membership"
    run([
        sys.executable, str(paths["apply_post_merge_atlas_routing.py"]),
        "--frozen-broad", str(args.frozen_broad_membership.resolve()),
        "--routing", str(routing_out / "atlas_state_routing.tsv.gz"),
        "--atlas-validation", str(validation_path),
        "--out", str(apply_out),
    ], output / "logs/02_post_atlas_membership.log")
    applied = json.loads(
        (apply_out / "post_atlas_membership_manifest.json").read_text(encoding="utf-8")
    )
    post_membership = Path(applied["membership"]["path"])
    profile = json.loads(paths["profile"].read_text(encoding="utf-8"))
    quality_required = bool(
        profile.get("biological_quality_endpoints", {}).get("required", False)
    )
    unresolved_review_authority = stage_authority(
        "atlas_and_completeness_review", args.contract, paths, output,
        post_atlas_membership=post_membership,
        observation_scores=observation_score_paths,
        cluster_evidence=cluster_evidence_paths,
    )
    unresolved_review_out = output / "03_post_merge_unresolved_review"
    unresolved_review_command = [
        sys.executable,
        str(paths["review_post_merge_unresolved_components.py"]),
        "--contract", str(args.contract),
        "--stage-authority", str(unresolved_review_authority),
        "--membership", str(post_membership),
        "--catalog", str(paths["catalog"]),
        "--workers", str(max(1, int(args.scoring_workers))),
        "--out", str(unresolved_review_out),
    ]
    for path in observation_score_paths:
        unresolved_review_command.extend(["--scores", str(path.resolve())])
    for path in cluster_evidence_paths:
        unresolved_review_command.extend(
            ["--cluster-evidence", str(path.resolve())]
        )
    run(
        unresolved_review_command,
        output / "logs/03_post_merge_unresolved_review.log",
    )
    unresolved_review_manifest_path = (
        unresolved_review_out / "post_merge_unresolved_review_manifest.json"
    )
    unresolved_review_manifest = json.loads(
        unresolved_review_manifest_path.read_text(encoding="utf-8")
    )
    post_membership = Path(unresolved_review_manifest["membership"]["path"])
    completeness_out = output / "04_completeness"
    command = [
        sys.executable, str(paths["audit_post_merge_completeness.py"]),
        "--membership", str(post_membership),
        "--catalog", str(paths["catalog"]),
        "--post-merge-review-manifest", str(unresolved_review_manifest_path),
        "--out", str(completeness_out),
    ]
    for path in cluster_evidence_paths:
        command.extend(["--cluster-evidence", str(path.resolve())])
    for path in fine_proposal_paths:
        command.extend(["--fine-audit", str(path.resolve())])
    for path in unmodeled_manifest_paths:
        command.extend(["--unmodeled", str(path.resolve())])
    if args.unmodeled_decisions:
        command.extend([
            "--unmodeled-review", str(args.unmodeled_decisions.resolve())
        ])
    for path in local_subset_validation_paths:
        command.extend(["--local-subset-validation", str(path.resolve())])
    for path in local_remainder_audit_paths:
        command.extend(["--local-remainder-audit", str(path.resolve())])
    if quality_required:
        command.append("--defer-canonical-zero-to-biological-review")
    run(command, output / "logs/04_completeness.log")
    completeness_manifest = completeness_out / "post_merge_completeness_manifest.json"
    quality_manifest: dict[str, str] | None = None
    quality_status = "NOT_REQUIRED"
    follicle_repair_manifest: dict[str, str] | None = None
    membership_output = unresolved_review_manifest["membership"]
    n_unresolved_biological = int(
        unresolved_review_manifest["n_remaining_unresolved"]
    )
    if quality_required:
        if not observation_score_paths:
            raise RuntimeError(
                "sheep-ovary biological review requires second-round observation scores"
            )
        quality_out = output / "05_sheep_ovary_biological_quality"
        quality_command = [
            sys.executable,
            str(paths["validate_sheep_ovary_biological_quality.py"]),
            "--membership", str(post_membership.resolve()),
            "--catalog", str(paths["catalog"]),
            "--out", str(quality_out),
        ]
        for path in observation_score_paths:
            quality_command.extend(["--scores", str(path.resolve())])
        run(
            quality_command,
            output / "logs/05_sheep_ovary_biological_quality.log",
            allowed_codes=(0, 2),
        )
        quality_path = quality_out / "sheep_ovary_biological_quality_review.json"
        quality_doc = json.loads(quality_path.read_text(encoding="utf-8"))
        quality_status = str(quality_doc.get("status", ""))
        if quality_status not in {"PASS", "ITERATION_REQUIRED"}:
            raise RuntimeError("invalid sheep-ovary biological quality status")
        quality_manifest = artifact(quality_path)
        if quality_status == "ITERATION_REQUIRED":
            repair = run_targeted_follicle_roi_iteration(
                args, contract, paths, output, post_membership, quality_path,
                observation_score_paths,
            )
            if repair is not None:
                post_membership = repair["membership"]
                quality_path = repair["quality_review"]
                quality_status = repair["status"]
                if quality_status not in {"PASS", "ITERATION_REQUIRED"}:
                    raise RuntimeError("invalid post-repair sheep-ovary quality status")
                quality_manifest = artifact(quality_path)
                follicle_repair_manifest = artifact(repair["repair_manifest"])
                membership_output = artifact(post_membership)
                post_rows = read_tsv(post_membership)
                n_unresolved_biological = sum(
                    not str(row.get("final_broad_label", row.get("broad_label", "")))
                    for row in post_rows
                )
                completeness_out = output / "06_post_repair_completeness"
                command = [
                    sys.executable,
                    str(paths["audit_post_merge_completeness.py"]),
                    "--membership", str(post_membership),
                    "--catalog", str(paths["catalog"]),
                    "--out", str(completeness_out),
                ]
                for path in cluster_evidence_paths:
                    command.extend(["--cluster-evidence", str(path.resolve())])
                for path in fine_proposal_paths:
                    command.extend(["--fine-audit", str(path.resolve())])
                for path in unmodeled_manifest_paths:
                    command.extend(["--unmodeled", str(path.resolve())])
                if args.unmodeled_decisions:
                    command.extend([
                        "--unmodeled-review",
                        str(args.unmodeled_decisions.resolve()),
                    ])
                for path in local_subset_validation_paths:
                    command.extend([
                        "--local-subset-validation", str(path.resolve())
                    ])
                for path in local_remainder_audit_paths:
                    command.extend([
                        "--local-remainder-audit", str(path.resolve())
                    ])
                if quality_required:
                    command.append(
                        "--defer-canonical-zero-to-biological-review"
                    )
                run(command, output / "logs/06_post_repair_completeness.log")
                completeness_manifest = (
                    completeness_out / "post_merge_completeness_manifest.json"
                )
    return {
        "status": (
            "ITERATION_REQUIRED"
            if quality_status == "ITERATION_REQUIRED" else "PASS"
        ),
        "phase": "atlas_and_completeness_review",
        "stage_authority": artifact(unresolved_review_authority),
        "prerequisite": artifact(args.prerequisite_manifest),
        "formal_membership_written": False,
        "membership": membership_output,
        "atlas_validation": artifact(validation_path),
        "post_merge_unresolved_review": artifact(
            unresolved_review_manifest_path
        ),
        "completeness": artifact(completeness_manifest),
        "biological_quality_status": quality_status,
        "biological_quality_review": quality_manifest,
        "targeted_follicle_roi_repair": follicle_repair_manifest,
        "unlabeled_broad_rescue_n": applied["unlabeled_broad_rescue_n"],
        "n_unresolved_biological": n_unresolved_biological,
        "fine_candidate_proposals": [artifact(path) for path in fine_proposal_paths],
        "state_annotation_proposals": [artifact(path) for path in state_proposal_paths],
    }


def phase_final(args, contract: dict, paths: dict[str, Path]) -> dict:
    prerequisite = validate_prerequisite(args, "atlas_and_completeness_review")
    membership_record = prerequisite.get("membership", {})
    if (
        Path(str(membership_record.get("path", ""))).resolve()
        != args.post_atlas_membership.resolve()
        or membership_record.get("sha256") != sha256(args.post_atlas_membership)
    ):
        raise RuntimeError(
            "final phase membership differs from the reviewed post-Atlas membership"
        )
    fine_proposal_paths = bound_artifact_paths(
        prerequisite.get("fine_candidate_proposals"), "fine proposal"
    )
    state_proposal_paths = bound_artifact_paths(
        prerequisite.get("state_annotation_proposals"), "state proposal"
    )
    output = args.out.resolve()
    fine_out = output / "00_fine_materialization"
    fine_command = [
        sys.executable, str(paths["materialize_parent_locked_fine_proposals.py"]),
        "--membership", str(args.post_atlas_membership.resolve()),
        "--catalog", str(paths["catalog"]), "--out", str(fine_out),
    ]
    for path in fine_proposal_paths:
        fine_command.extend(["--fine-audit", str(path.resolve())])
    run(fine_command, output / "logs/00_fine_materialization.log")
    fine_assignments = fine_out / "parent_locked_fine_assignments.tsv.gz"
    authority = stage_authority(
        "materialize_final_release", args.contract, paths, output,
        post_atlas_membership=args.post_atlas_membership,
        prerequisite_manifest=args.prerequisite_manifest,
        fine_assignments=fine_assignments,
        state_annotation_proposals=state_proposal_paths,
    )
    final_out = output / "01_final_release"
    final_command = [
        sys.executable, str(paths["materialize_final_release_v2_2.py"]),
        "--contract", str(args.contract), "--stage-authority", str(authority),
        "--post-atlas-membership", str(args.post_atlas_membership.resolve()),
        "--atlas-completeness-manifest", str(args.prerequisite_manifest.resolve()),
        "--fine-assignments", str(fine_assignments),
        "--out", str(final_out),
    ]
    for path in state_proposal_paths:
        final_command.extend(["--state-proposals", str(path.resolve())])
    run(final_command, output / "logs/01_final_release.log")
    final_manifest = json.loads(
        (final_out / "final_release_manifest.json").read_text(encoding="utf-8")
    )
    return {
        "status": "PASS", "phase": "materialize_final_release",
        "stage_authority": artifact(authority),
        "prerequisite": artifact(args.prerequisite_manifest),
        "formal_membership_written": True,
        "membership": final_manifest["membership"],
        "broad_census": final_manifest["broad_census"],
        "fine_census": final_manifest["fine_census"],
        "state_census": final_manifest["state_census"],
        "residual_qc_n": final_manifest["residual_qc_n"],
        "residual_qc_fraction": final_manifest["residual_qc_fraction"],
    }


def add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=2200)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument(
        "--scoring-workers", type=int,
        default=1,
        help=(
            "marker-scoring fork count; default 1 because each worker holds "
            "full observation-scale arrays. Scheduler CPUs remain available "
            "to SCT, graph construction and resolution evaluation."
        ),
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="phase", required=True)
    whole = sub.add_parser("whole_tissue_partition")
    add_common(whole)
    whole.add_argument("--rds", required=True, type=Path)
    whole.add_argument("--partitions", required=True, type=Path)

    cohort = sub.add_parser("cluster_cohort_recluster")
    add_common(cohort)
    cohort.add_argument("--rds", required=True, type=Path)
    cohort.add_argument("--membership", required=True, type=Path)
    cohort.add_argument("--whole-manifest", required=True, type=Path)
    cohort.add_argument("--cohort-id", required=True)
    cohort.add_argument("--source-initial-cluster", required=True)
    cohort.add_argument(
        "--provisional-status", required=True,
        choices=["provisional_broad", "mixed", "unknown"],
    )
    cohort.add_argument("--provisional-broad", default="")
    cohort.add_argument("--context-evidence", type=Path)
    cohort.add_argument(
        "--reuse-recluster-manifest", type=Path,
        help=(
            "reuse only a validated raw-count SCT/PCA/SNN/Leiden partition "
            "cache with an identical input/membership/parameter fingerprint; "
            "all scoring, selection and writeback are recomputed"
        ),
    )
    cohort.add_argument("--resolution-contract", default="sheep_ovary", choices=["generic", "sheep_ovary"])
    cohort.add_argument("--resolution-workers", type=int, default=5)

    local = sub.add_parser("local_mixed_subcluster_split")
    add_common(local)
    local.add_argument("--trigger-manifest", required=True, type=Path)
    local.add_argument("--scoring-output", required=True, type=Path)
    local.add_argument("--source-boundary", required=True)
    local.add_argument("--source-cluster", required=True)
    local.add_argument("--subset-workers", type=int, default=max(1, int(os.environ.get("LSB_DJOB_NUMPROC", "1"))))

    merge = sub.add_parser("merge_and_freeze_broad")
    add_common(merge)
    merge.add_argument("--analysis-membership", required=True, type=Path)
    merge.add_argument("--whole-manifest", required=True, type=Path)
    merge.add_argument("--cohort-manifest", required=True, type=Path, action="append")
    merge.add_argument("--local-manifest", type=Path, action="append", default=[])

    atlas = sub.add_parser("atlas_and_completeness_review")
    add_common(atlas)
    atlas.add_argument("--prerequisite-manifest", required=True, type=Path)
    atlas.add_argument("--frozen-broad-membership", required=True, type=Path)
    atlas.add_argument("--atlas-mapping", required=True, type=Path)
    atlas.add_argument("--calibration-manifest", required=True, type=Path)
    atlas.add_argument("--atlas-decisions", type=Path)
    atlas.add_argument("--unmodeled-decisions", type=Path)
    atlas.add_argument("--resolution-workers", type=int, default=5)

    final = sub.add_parser("materialize_final_release")
    add_common(final)
    final.add_argument("--prerequisite-manifest", required=True, type=Path)
    final.add_argument("--post-atlas-membership", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.contract = args.contract.resolve()
    contract, paths = validate_contract(args.contract)
    bound_seed = int(
        contract.get("canonical_lineage_controller", {}).get("random_seed", -1)
    )
    if args.seed != bound_seed:
        raise RuntimeError(
            f"controller seed {args.seed} differs from contract-bound seed {bound_seed}"
        )
    if args.scoring_workers < 1:
        raise RuntimeError("scoring workers must be >= 1")
    scoring_n = None
    if args.phase == "whole_tissue_partition":
        scoring_n = int(contract.get("input_scope", {}).get("analysis_set_n", 0))
    elif args.phase == "cluster_cohort_recluster":
        scoring_n = len(read_tsv(args.membership))
    if scoring_n is not None and scoring_n >= 100000 and args.scoring_workers > 1:
        raise RuntimeError(
            "observation-scale lineage scoring with >=100,000 observations "
            "must use --scoring-workers 1; forked workers replicate large "
            "direct/local arrays and can exceed node memory"
        )
    scheduler_limits = []
    for name in ("LSB_DJOB_NUMPROC", "SLURM_CPUS_PER_TASK", "NSLOTS", "AIP_CPUS"):
        try:
            value = int(os.environ.get(name, ""))
        except ValueError:
            continue
        if value > 0:
            scheduler_limits.append(value)
    if scheduler_limits and args.scoring_workers > max(scheduler_limits):
        raise RuntimeError("scoring workers exceed the scheduler CPU allocation")
    if args.phase == "whole_tissue_partition":
        result = phase_whole(args, contract, paths)
    elif args.phase == "cluster_cohort_recluster":
        result = phase_cohort(args, contract, paths)
    elif args.phase == "local_mixed_subcluster_split":
        result = phase_local(args, contract, paths)
    elif args.phase == "merge_and_freeze_broad":
        result = phase_merge(args, contract, paths)
    elif args.phase == "atlas_and_completeness_review":
        result = phase_atlas(args, contract, paths)
    else:
        result = phase_final(args, contract, paths)
    result.update({
        "schema_version": "2.2", "controller_version": "2.2.0",
        "annotation_contract": artifact(args.contract),
        "seed": args.seed, "finished_at_utc": datetime.now(timezone.utc).isoformat(),
    })
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = args.out / "lineage_controller_manifest.json"
    manifest.write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
