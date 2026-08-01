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
from controller_runtime_state import materialize_runtime_state
from lineage_controller_lib import deterministic_cell_id_set_hash
from membership_transform_lib import load_and_validate_chain
from runtime_dependency_registry import (
    CANONICAL_DEPENDENCIES,
    CANONICAL_SCRIPTS,
    PHASE_ORDER,
    REGISTRY_VERSION,
)
from validate_fixed_atlas_bundle import (
    ACTIVE_BUNDLE_ID,
    validate as validate_fixed_atlas,
)


PHASES = PHASE_ORDER
WHOLE_TISSUE_FORK_WORKER_CAP = 64
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


def initialize_membership_transform_chain(
    paths: dict[str, Path], membership: Path, output: Path, log_root: Path,
) -> Path:
    out = output / "T0000_initial"
    run([
        sys.executable, str(paths["manage_membership_transform_chain.py"]),
        "init", "--membership", str(membership.resolve()), "--out", str(out),
    ], log_root / "membership_transform_T0000.log")
    return out / "membership_transform_chain.json"


def append_membership_transform(
    paths: dict[str, Path], chain: Path, operation: str, source: Path,
    result: Path, evidence_manifest: Path, output: Path, index: int,
    log_root: Path, target_cell_type: str = "",
) -> Path:
    out = output / f"T{index:04d}_{operation}"
    command = [
        sys.executable, str(paths["manage_membership_transform_chain.py"]),
        "append", "--chain", str(chain.resolve()), "--operation", operation,
        "--source", str(source.resolve()), "--result", str(result.resolve()),
        "--evidence-manifest", str(evidence_manifest.resolve()),
        "--out", str(out),
    ]
    if target_cell_type:
        command.extend(["--target-cell-type", target_cell_type])
    run(command, log_root / f"membership_transform_T{index:04d}.log")
    return out / "membership_transform_chain.json"


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
    if controller.get("runtime_dependency_registry_version") != REGISTRY_VERSION:
        raise RuntimeError("annotation contract binds a stale runtime dependency registry")
    observation_unit = str(contract.get("observation_unit", "")).strip().lower()
    if observation_unit not in {"cell", "nucleus", "cellbin", "spot"}:
        raise RuntimeError("annotation contract has no supported observation_unit")
    release_taxonomy = contract.get("release_taxonomy", {})
    if (
        release_taxonomy.get("independent_vascular_lineages")
        != ["Endothelial", "Pericyte/mural", "Smooth muscle"]
        or release_taxonomy.get("lymphatic_parent") != "Endothelial"
        or release_taxonomy.get("legacy_vascular_associated_release_forbidden") is not True
        or release_taxonomy.get("single_public_annotation_column") != "final_cell_type"
    ):
        raise RuntimeError("annotation contract does not bind the v2.2 vascular taxonomy")
    script_dir = Path(__file__).resolve().parent
    if set(controller.get("scripts", {})) != set(CANONICAL_SCRIPTS):
        raise RuntimeError("contract script set differs from the canonical registry")
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
    if set(controller.get("dependencies", {})) != set(CANONICAL_DEPENDENCIES):
        raise RuntimeError("contract dependency set differs from the canonical registry")
    for dependency_name in CANONICAL_DEPENDENCIES:
        dependency = script_dir / dependency_name
        record = controller.get("dependencies", {}).get(dependency_name, {})
        if (
            Path(str(record.get("path", ""))).resolve() != dependency.resolve()
            or not dependency.is_file()
            or record.get("sha256") != sha256(dependency)
        ):
            raise RuntimeError(
                f"contract does not bind {dependency_name}"
            )
        paths[dependency_name] = dependency
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
    paths["threshold_registry"] = resolve_bound(
        contract_path, contract.get("threshold_registry", {}),
        "controller threshold registry",
    )
    paths["profile"] = resolve_bound(
        contract_path, contract.get("biological_profile", {}), "biological profile"
    )
    paths["workflow_profile"] = resolve_bound(
        contract_path, contract.get("workflow_profile", {}), "workflow profile"
    )
    paths["catalog"] = resolve_bound(
        contract_path, contract.get("candidate_catalog", {}), "candidate catalog"
    )
    paths["atlas_bundle"] = resolve_bound(
        contract_path, contract.get("atlas_bundle", {}), "fixed Atlas bundle"
    )
    atlas_bundle = json.loads(paths["atlas_bundle"].read_text(encoding="utf-8"))
    atlas_errors = validate_fixed_atlas(atlas_bundle)
    if atlas_errors:
        raise RuntimeError(
            "fixed GSE233801 Atlas descriptor is invalid: "
            + "; ".join(atlas_errors)
        )
    atlas_policy = contract.get("atlas_routing", {})
    if not (
        atlas_policy.get("bundle_id") == atlas_bundle.get("bundle_id")
        and atlas_policy.get("capability_matrix_enforced") is True
        and atlas_policy.get("runtime_atlas_substitution_forbidden") is True
    ):
        raise RuntimeError("contract does not enforce the fixed GSE233801 Atlas")
    context_record = contract.get("candidate_context_evidence")
    if context_record:
        paths["context_evidence"] = resolve_bound(
            contract_path, context_record, "candidate context evidence"
        )
    input_audit_record = contract.get("query_input_audit")
    if input_audit_record:
        paths["query_input_audit"] = resolve_bound(
            contract_path, input_audit_record, "query input audit"
        )
        paths["raw_count_assay"] = str(
            input_audit_record.get("raw_count_assay", "")
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
        "query_membership_semantic_sha256": deterministic_cell_id_set_hash(
            read_tsv(args.membership)
        ),
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
            name: artifact(paths[name])
            for name in (
                "controller_thresholds.py", "lineage_controller_lib.py"
            )
        },
        "threshold_registry": artifact(paths["threshold_registry"]),
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
    assay: str = "", observation_unit: str = "",
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
        "--threshold-registry", str(paths["threshold_registry"]),
        "--observation-unit", observation_unit,
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
        observation_unit=contract["observation_unit"],
    )
    evidence_out = output / "01_resolution_evidence"
    run([
        sys.executable, str(paths["build_resolution_grid_evidence.py"]),
        "--scoring-output", str(grid_scoring),
        "--catalog", str(paths["catalog"]),
        "--selection-purpose", "whole_tissue_cohort_partition",
        "--out", str(evidence_out),
        *(
            ["--context-evidence", str(paths["context_evidence"])]
            if "context_evidence" in paths else []
        ),
    ], output / "logs/01_resolution_evidence.log")
    selector_out = output / "02_resolution_selection"
    run([
        sys.executable, str(paths["select_lineage_resolution.py"]),
        "--grid-evidence", str(evidence_out / "resolution_grid_evidence.tsv"),
        "--selection-purpose", "whole_tissue_cohort_partition",
        "--threshold-registry", str(paths["threshold_registry"]),
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
    resource_out = output / "00_resource_plan"
    resource_command = [
        sys.executable, str(paths["plan_cohort_resources.py"]),
        "--rds", str(args.rds.resolve()),
        "--rds-sha256", str(
            contract["selected_input_snapshot"]["sha256"]
        ),
        "--membership", str(args.membership.resolve()),
        "--resolution-workers", str(args.resolution_workers),
        "--grid-size", str(len(grid)), "--out", str(resource_out),
    ]
    if args.allocated_memory_gb is not None:
        resource_command.extend([
            "--allocated-memory-gb", str(args.allocated_memory_gb)
        ])
    run(resource_command, output / "logs/00_resource_plan.log")
    resource_plan_path = resource_out / "cohort_resource_plan.json"
    resource_plan = json.loads(resource_plan_path.read_text(encoding="utf-8"))
    workers = int(resource_plan["effective_resolution_workers"])
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
        observation_unit=contract["observation_unit"],
    )
    evidence_out = output / "02_resolution_evidence"
    run([
        sys.executable, str(paths["build_resolution_grid_evidence.py"]),
        "--scoring-output", str(grid_scoring), "--catalog", str(paths["catalog"]),
        "--selection-purpose", "cohort_identity_resolution",
        "--out", str(evidence_out),
        *(
            ["--context-evidence", str(paths["context_evidence"])]
            if "context_evidence" in paths else []
        ),
    ], output / "logs/02_resolution_evidence.log")
    selection_out = output / "03_resolution_selection"
    run([
        sys.executable, str(paths["select_lineage_resolution.py"]),
        "--grid-evidence", str(evidence_out / "resolution_grid_evidence.tsv"),
        "--selection-purpose", "cohort_identity_resolution",
        "--threshold-registry", str(paths["threshold_registry"]),
        "--out", str(selection_out),
    ], output / "logs/03_resolution_selection.log")
    selection = json.loads(
        (selection_out / "resolution_selection.json").read_text(encoding="utf-8")
    )
    selected_partitions = selection_out / "selected_neighbor_partitions.tsv.gz"
    materialize_selected_neighbors(grid_partitions, selection, selected_partitions)
    scoring = output / "04_selected_scoring"
    run_scorer(
        args, paths, query_rds, selected_partitions, scoring, assay=raw_assay,
        observation_unit=contract["observation_unit"],
    )
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
        "resource_plan": artifact(resource_plan_path),
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
        "resource_plan": artifact(resource_plan_path),
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
        "n_observations": len(membership_ids),
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
    workload = json.loads(args.workload_audit.read_text(encoding="utf-8"))
    if (
        workload.get("stage") != "pre_local_split_workload_audit"
        or workload.get("status") != "PASS"
    ):
        raise RuntimeError(
            "local split requires a PASS combined-cohort workload audit"
        )
    threshold_record = workload.get("threshold_registry", {})
    if (
        Path(str(threshold_record.get("path", ""))).resolve()
        != paths["threshold_registry"].resolve()
        or threshold_record.get("sha256")
        != sha256(paths["threshold_registry"])
    ):
        raise RuntimeError("local split workload audit uses stale thresholds")
    bound_triggers = {
        (
            Path(str(record.get("path", ""))).resolve(),
            str(record.get("sha256", "")),
        )
        for record in workload.get("cohort_manifests", [])
    }
    if (
        args.trigger_manifest.resolve(), sha256(args.trigger_manifest)
    ) not in bound_triggers:
        raise RuntimeError(
            "local split trigger is absent from the combined-cohort workload audit"
        )
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
        local_split_workload_audit=args.workload_audit,
        context_evidence=context_path,
    )
    subset_out = output / "00_candidate_subsets"
    subset_command = [
        args.rscript, str(paths["derive_candidate_local_subsets.R"]),
        "--scores", str(scores), "--cluster-evidence", str(evidence),
        "--catalog", str(paths["catalog"]), "--release-level", "broad",
        "--source-boundary", args.source_boundary,
        "--source-cluster", args.source_cluster,
        "--threshold-registry", str(paths["threshold_registry"]),
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
        "local_split_workload_audit": artifact(args.workload_audit),
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
        candidate_catalog=paths["catalog"],
        context_evidence=paths.get("context_evidence"),
    )
    command = [
        sys.executable, str(paths["merge_and_freeze_broad_membership.py"]),
        "--contract", str(args.contract.resolve()),
        "--stage-authority", str(authority),
        "--analysis-membership", str(args.analysis_membership.resolve()),
        "--catalog", str(paths["catalog"]),
        "--out", str(output / "00_broad_freeze"),
    ]
    if "context_evidence" in paths:
        command.extend([
            "--context-evidence", str(paths["context_evidence"])
        ])
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
    contract_sha = sha256(contract_path)
    manifest_contract = manifest.get("annotation_contract", {})
    compatible = False
    if manifest_contract.get("sha256") != contract_sha:
        contract = json.loads(contract_path.read_text(encoding="utf-8"))
        for record in contract.get("resume_compatible_contracts", []):
            path = Path(str((record or {}).get("path", "")))
            if (
                path.is_file()
                and (record or {}).get("sha256") == sha256(path)
                and manifest_contract.get("sha256") == sha256(path)
                and Path(str(manifest_contract.get("path", ""))).resolve()
                == path.resolve()
            ):
                compatible = True
                break
    return (
        manifest.get("status") == "PASS"
        and manifest.get("controller_version") == "2.2.0"
        and manifest.get("phase") == phase
        and (manifest_contract.get("sha256") == contract_sha or compatible)
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
    all_roi_rows = [
        row for row in actions
        if row.get("endpoint") == "follicle_roi_histology"
        and re.fullmatch(r"F\d+", str(row.get("scope_id", "")))
    ]
    if not all_roi_rows or len(all_roi_rows) != len(actions):
        return None
    # The bounded ROI repair writer is deliberately restricted to follicle-wall
    # identities.  A Granulosa boundary whose released label lacks a matching
    # program must proceed to the formal Granulosa whole-query review; it cannot
    # be silently treated as a wall-layer repair target.
    repairable_tokens = (
        "theca", "endothelial", "vascular", "pericyte", "mural",
        "lymphatic", "smooth", "contractile", "stromal",
        "order_inverted", "not_interleaved", "background_not_resolved",
        "identity_boundary_blurred", "complete_follicle_wall",
    )
    target_rows = [
        row for row in all_roi_rows
        if any(
            token in str(row.get("issue_code", ""))
            for token in repairable_tokens
        )
    ]
    if not target_rows:
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
            observation_unit=contract["observation_unit"],
        )
        evidence_out = roi_output / "02_resolution_evidence"
        run([
            sys.executable, str(paths["build_resolution_grid_evidence.py"]),
            "--scoring-output", str(grid_scoring),
            "--catalog", str(paths["catalog"]),
            "--selection-purpose", "cohort_identity_resolution",
            "--out", str(evidence_out),
            *(
                ["--context-evidence", str(paths["context_evidence"])]
                if "context_evidence" in paths else []
            ),
        ], roi_output / "logs/02_resolution_evidence.log")
        selection_out = roi_output / "03_resolution_selection"
        run([
            sys.executable, str(paths["select_lineage_resolution.py"]),
            "--grid-evidence", str(evidence_out / "resolution_grid_evidence.tsv"),
            "--selection-purpose", "cohort_identity_resolution",
            "--threshold-registry", str(paths["threshold_registry"]),
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
            assay=raw_assay, observation_unit=contract["observation_unit"],
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
        candidate_catalog=paths["catalog"],
        context_evidence=paths.get("context_evidence"),
    )
    materialized = output_root / "materialized_repair"
    command = [
        sys.executable,
        str(paths["apply_sheep_ovary_follicle_roi_repair.py"]),
        "--contract", str(args.contract),
        "--stage-authority", str(repair_authority),
        "--membership", str(post_membership),
        "--quality-review", str(quality_path),
        "--catalog", str(paths["catalog"]),
        "--out", str(materialized),
    ]
    if "context_evidence" in paths:
        command.extend(["--context-evidence", str(paths["context_evidence"])])
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
    coordinate_membership = Path(
        repair_manifest["coordinate_membership"]["path"]
    )
    post_review = output_root / "post_repair_biological_quality"
    run([
        sys.executable,
        str(paths["validate_sheep_ovary_biological_quality.py"]),
        "--membership", str(repaired_membership),
        "--coordinate-membership", str(coordinate_membership),
        "--catalog", str(paths["catalog"]),
        "--scores", str(combined_scores),
        "--expected-roi-review", str(expected_roi_path),
        *(
            ["--canonical-oocyte-review", str(args.canonical_oocyte_review)]
            if getattr(args, "canonical_oocyte_review", None) else []
        ),
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
        "coordinate_membership": coordinate_membership,
        "target_rois": target_rois,
    }


def prepare_catalog_review_static_evidence(
    args, paths: dict[str, Path], output: Path, membership: Path,
) -> tuple[Path, Path, Path]:
    """Build or reuse immutable marker/count evidence for serial review."""
    static_out = output / "07_catalog_wide_review_evidence_static"
    manifest_path = static_out / "catalog_review_static_evidence_manifest.json"
    supplied = getattr(args, "lineage_review_static_evidence_manifest", None)
    candidate = supplied.resolve() if supplied else manifest_path
    cell_set_sha = deterministic_cell_id_set_hash(read_tsv(membership))

    if candidate.is_file():
        document = json.loads(candidate.read_text(encoding="utf-8"))
        expected = {
            "selected_input": paths["selected_input"],
            "profile": paths["profile"],
            "candidate_catalog": paths["catalog"],
            "threshold_registry": paths["threshold_registry"],
        }
        compatible = (
            document.get("status") == "PASS"
            and document.get("artifact_role")
            == "catalog_review_static_raw_count_evidence"
            and document.get("analysis_cell_id_set_sha256") == cell_set_sha
        )
        for key, expected_path in expected.items():
            record = document.get(key, {})
            compatible = compatible and (
                Path(str(record.get("path", ""))).resolve()
                == expected_path.resolve()
                and record.get("sha256") == sha256(expected_path)
            )
        artifacts = document.get("artifacts", {})
        marker_path = Path(str(artifacts.get("marker_manifest", {}).get("path", "")))
        count_manifest = Path(str(artifacts.get("count_export_manifest", {}).get("path", "")))
        for key, path in (
            ("marker_manifest", marker_path),
            ("count_export_manifest", count_manifest),
        ):
            record = artifacts.get(key, {})
            compatible = compatible and path.is_file() and record.get("sha256") == sha256(path)
        count_document = (
            json.loads(count_manifest.read_text(encoding="utf-8"))
            if compatible else {}
        )
        for key in (
            "marker_matrix", "full_count_matrix", "gene_map", "cells",
            "library_size", "coordinates",
        ):
            path = Path(str(count_document.get(key, "")))
            record = artifacts.get(f"count_export_{key}", {})
            compatible = compatible and path.is_file() and record.get("sha256") == sha256(path)
        if compatible:
            return candidate, marker_path, count_manifest.parent
        if supplied:
            raise RuntimeError("supplied catalog review static evidence is stale or incompatible")

    marker_out = static_out / "00_marker_manifest"
    run([
        sys.executable,
        str(paths["build_cell_type_review_marker_manifest.py"]),
        "--profile", str(paths["profile"]),
        "--catalog", str(paths["catalog"]),
        "--out", str(marker_out),
    ], output / "logs/07_catalog_wide_marker_manifest.log")
    marker_path = marker_out / "cell_type_review_marker_manifest.tsv"
    count_root = static_out / "01_raw_count_export"
    run([
        args.rscript, str(paths["export_cell_type_review_counts.R"]),
        "--rds", str(paths["selected_input"]),
        "--analysis-membership", str(membership.resolve()),
        "--marker-manifest", str(marker_path),
        "--out", str(count_root),
    ], output / "logs/07_catalog_wide_raw_count_export.log")
    count_manifest = count_root / "cell_type_review_count_export_manifest.json"
    count_document = json.loads(count_manifest.read_text(encoding="utf-8"))
    document = {
        "schema_version": "2.5",
        "status": "PASS",
        "artifact_role": "catalog_review_static_raw_count_evidence",
        "analysis_cell_id_set_sha256": cell_set_sha,
        "selected_input": artifact(paths["selected_input"]),
        "profile": artifact(paths["profile"]),
        "candidate_catalog": artifact(paths["catalog"]),
        "threshold_registry": artifact(paths["threshold_registry"]),
        "artifacts": {
            "marker_manifest": artifact(marker_path),
            "marker_manifest_summary": artifact(
                marker_out / "cell_type_review_marker_manifest.json"
            ),
            "count_export_manifest": artifact(count_manifest),
            **{
                f"count_export_{key}": artifact(Path(str(count_document[key])))
                for key in (
                    "marker_matrix", "full_count_matrix", "gene_map", "cells",
                    "library_size", "coordinates",
                )
            },
        },
    }
    static_out.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest_path, marker_path, count_root


def run_catalog_wide_review_iterations(
    args, paths: dict[str, Path], output: Path, membership: Path,
    observation_score_paths: list[Path], cluster_evidence_paths: list[Path],
    biological_quality_required: bool,
    membership_already_includes_prior_applies: bool = False,
) -> dict:
    """Run exactly one formal per-cell-type review step without reclustering."""
    policy = json.loads(
        paths["threshold_registry"].read_text(encoding="utf-8")
    )["catalog_wide_lineage_review_policy"]
    maximum_decisions = int(policy["maximum_decision_rounds"])
    supplied_decisions = list(args.lineage_review_decisions or [])
    if len(supplied_decisions) > 1:
        raise RuntimeError(
            "one controller invocation may decide only one active cell type"
        )
    current_membership = membership
    prior_validations: list[Path] = list(
        args.lineage_review_prior_validation or []
    )
    apply_manifests: list[Path] = list(
        args.lineage_review_prior_apply or []
    )
    evidence_packet_manifests: list[Path] = list(
        args.lineage_review_prior_packet or []
    )
    manual_adjudications: list[Path] = list(
        args.lineage_review_manual_adjudication or []
    )
    previous_review: Path | None = args.lineage_review_previous_review
    last_review: Path | None = None

    prior_replay_membership = current_membership
    for index, manifest_path in enumerate(apply_manifests, 1):
        document = json.loads(manifest_path.read_text(encoding="utf-8"))
        source = document.get("source_membership", {})
        result = document.get("membership", {})
        source_path = Path(str(source.get("path", "")))
        result_path = Path(str(result.get("path", "")))
        if (
            document.get("stage") != "catalog_wide_lineage_review_apply"
            or document.get("formal_batch_closure_performed") is not False
            or not source_path.is_file()
            or source.get("sha256") != sha256(source_path)
            or not result_path.is_file()
            or result.get("sha256") != sha256(result_path)
        ):
            raise RuntimeError(
                f"prior cell-type review apply {index} is stale or non-contiguous"
            )
        if not membership_already_includes_prior_applies:
            if source_path.resolve() != current_membership.resolve():
                raise RuntimeError(
                    f"prior cell-type review apply {index} is non-contiguous"
                )
            current_membership = result_path
        elif index > 1:
            previous = json.loads(
                apply_manifests[index - 2].read_text(encoding="utf-8")
            )
            previous_result = Path(str(previous["membership"]["path"]))
            if source_path.resolve() != previous_result.resolve():
                raise RuntimeError(
                    f"prior cell-type review apply {index} is non-contiguous"
                )
    if membership_already_includes_prior_applies and apply_manifests:
        final_prior = json.loads(
            apply_manifests[-1].read_text(encoding="utf-8")
        )
        final_prior_path = Path(str(final_prior["membership"]["path"]))
        if final_prior_path.resolve() != prior_replay_membership.resolve():
            raise RuntimeError(
                "resume membership differs from the last prior cell-type review apply"
            )

    static_evidence_manifest, marker_manifest, count_out = (
        prepare_catalog_review_static_evidence(args, paths, output, membership)
    )
    for review_round in range(len(prior_validations) + 1, len(prior_validations) + 2):
        review_out = output / f"07_catalog_wide_review_round_{review_round:02d}"
        zero_out = review_out / "00_zero_census_direct_challengers"
        zero_command = [
            sys.executable,
            str(paths["build_zero_census_direct_challengers.py"]),
            "--count-export", str(count_out),
            "--marker-manifest", str(marker_manifest),
            "--membership", str(current_membership.resolve()),
            "--catalog", str(paths["catalog"]),
            "--threshold-registry", str(paths["threshold_registry"]),
            "--out", str(zero_out),
        ]
        if "context_evidence" in paths:
            zero_command.extend([
                "--context-evidence", str(paths["context_evidence"])
            ])
        run(
            zero_command,
            output
            / f"logs/07_zero_census_challengers_round_{review_round:02d}.log",
        )
        zero_challenger_manifest = (
            zero_out / "zero_census_direct_challenger_manifest.json"
        )
        review_authority = stage_authority(
            "atlas_and_completeness_review", args.contract, paths, review_out,
            post_atlas_membership=current_membership,
            candidate_catalog=paths["catalog"],
            threshold_registry=paths["threshold_registry"],
            context_evidence=paths.get("context_evidence"),
            observation_scores=observation_score_paths,
            cluster_evidence=cluster_evidence_paths,
            zero_census_direct_challenger=zero_challenger_manifest,
        )
        command = [
            sys.executable,
            str(paths["audit_catalog_wide_lineage_challengers.py"]),
            "--contract", str(args.contract),
            "--stage-authority", str(review_authority),
            "--membership", str(current_membership.resolve()),
            "--catalog", str(paths["catalog"]),
            "--threshold-registry", str(paths["threshold_registry"]),
            "--round-index", str(review_round),
            "--zero-census-direct-challenger-manifest",
            str(zero_challenger_manifest),
            "--workers", str(max(1, int(args.resolution_workers))),
            "--out", str(review_out),
        ]
        if "context_evidence" in paths:
            command.extend([
                "--context-evidence", str(paths["context_evidence"])
            ])
        if previous_review:
            command.extend([
                "--previous-review-manifest", str(previous_review.resolve())
            ])
        for path in prior_validations:
            command.extend([
                "--prior-decision-validation", str(path.resolve())
            ])
        for path in manual_adjudications:
            command.extend([
                "--manual-biological-adjudication", str(path.resolve())
            ])
        for path in observation_score_paths:
            command.extend(["--scores", str(path.resolve())])
        for path in cluster_evidence_paths:
            command.extend(["--cluster-evidence", str(path.resolve())])
        run(
            command,
            output / f"logs/07_catalog_wide_review_round_{review_round:02d}.log",
            allowed_codes=(0, 2),
        )
        last_review = review_out / "catalog_wide_lineage_review_manifest.json"
        review = json.loads(last_review.read_text(encoding="utf-8"))
        state_out = review_out / "00_sequential_cell_type_state"
        state_command = [
            sys.executable,
            str(paths["manage_cell_type_review_queue.py"]),
            "--review-manifest", str(last_review.resolve()),
            "--maximum-decisions-per-cell-type", str(maximum_decisions),
            "--out", str(state_out),
        ]
        if args.lineage_review_previous_state:
            state_command.extend([
                "--previous-state",
                str(args.lineage_review_previous_state.resolve()),
            ])
        for path in prior_validations:
            state_command.extend([
                "--prior-decision-validation", str(path.resolve())
            ])
        for path in manual_adjudications:
            state_command.extend([
                "--manual-biological-adjudication", str(path.resolve())
            ])
        run(
            state_command,
            output / f"logs/07_cell_type_state_round_{review_round:02d}.log",
            allowed_codes=(0, 2),
        )
        review_state_path = state_out / "cell_type_review_state.json"
        review_state = json.loads(
            review_state_path.read_text(encoding="utf-8")
        )
        if review.get("status") == "PASS" and review_state.get("status") != "COMPLETE":
            raise RuntimeError(
                "PASS catalog review did not materialize a COMPLETE serial state"
            )
        active = review_state.get("active_cell_type_review") or {}
        active_review_id = str(active.get("review_id", ""))
        if (
            review_state.get("status") == "COMPLETE"
            and review_state.get("next_action") == "all_cell_types_closed"
            and int(review_state.get("active_review_n", -1)) == 0
            and int(review_state.get("queued_review_n", -1)) == 0
            and int(review_state.get("blocked_review_n", -1)) == 0
        ):
            if supplied_decisions:
                raise RuntimeError("unused catalog-wide review decision file remains")
            return {
                "status": "PASS", "membership": current_membership,
                "review_manifest": last_review,
                "decision_validations": prior_validations,
                "apply_manifests": apply_manifests,
                "evidence_packet_manifests": evidence_packet_manifests,
                "manual_adjudications": manual_adjudications,
                "static_evidence_manifest": static_evidence_manifest,
                "review_state": review_state_path,
            }
        if review_state.get("status") == "BLOCKED":
            return {
                "status": "ITERATION_REQUIRED",
                "membership": current_membership,
                "review_manifest": last_review,
                "decision_validations": prior_validations,
                "apply_manifests": apply_manifests,
                "evidence_packet_manifests": evidence_packet_manifests,
                "manual_adjudications": manual_adjudications,
                "static_evidence_manifest": static_evidence_manifest,
                "review_state": review_state_path,
                "active_cell_type": "",
                "required_progress_message": (
                    "该细胞类型的一次专项复核没有形成可验证结论，需要人工生物学裁决。"
                ),
            }
        if review_state.get("status") != "REVIEW_REQUIRED" or not active_review_id:
            raise RuntimeError("open review did not yield exactly one active cell type")
        evidence_root = review_out / "broad_cell_type_evidence"
        biological_quality_path: Path | None = None
        if biological_quality_required:
            biological_quality_out = evidence_root / "00_sheep_ovary_quality"
            quality_command = [
                sys.executable,
                str(paths["validate_sheep_ovary_biological_quality.py"]),
                "--membership", str(current_membership.resolve()),
                "--catalog", str(paths["catalog"]),
                "--out", str(biological_quality_out),
            ]
            for score_path in observation_score_paths:
                quality_command.extend(["--scores", str(score_path.resolve())])
            if args.canonical_oocyte_review:
                quality_command.extend([
                    "--canonical-oocyte-review",
                    str(args.canonical_oocyte_review.resolve()),
                ])
            for path in args.lineage_review_manual_adjudication:
                quality_command.extend([
                    "--manual-biological-adjudication", str(path.resolve())
                ])
            run(
                quality_command,
                output / f"logs/07_catalog_wide_quality_round_{review_round:02d}.log",
                allowed_codes=(0, 2),
            )
            biological_quality_path = (
                biological_quality_out
                / "sheep_ovary_biological_quality_review.json"
            )
        broad_evidence_out = evidence_root / "00_raw_count_marker_spatial"
        broad_command = [
            sys.executable,
            str(paths["build_broad_cell_type_review_evidence.py"]),
            "--count-export", str(count_out),
            "--marker-manifest", str(marker_manifest),
            "--membership", str(current_membership.resolve()),
            "--review-manifest", str(last_review.resolve()),
            "--active-review-id", active_review_id,
            "--catalog", str(paths["catalog"]),
            "--threshold-registry", str(paths["threshold_registry"]),
            "--zero-census-direct-challenger-manifest",
            str(zero_challenger_manifest),
            "--workers", str(max(1, int(args.resolution_workers))),
            "--out", str(broad_evidence_out),
        ]
        if "context_evidence" in paths:
            broad_command.extend([
                "--context-evidence", str(paths["context_evidence"])
            ])
        run(
            broad_command,
            output / f"logs/07_catalog_wide_evidence_round_{review_round:02d}.log",
        )
        broad_evidence_manifest = (
            broad_evidence_out / "broad_cell_type_review_evidence_manifest.json"
        )
        pseudobulk_out = evidence_root / "01_full_transcriptome_pseudobulk"
        count_export_document = json.loads(
            (count_out / "cell_type_review_count_export_manifest.json").read_text(
                encoding="utf-8"
            )
        )
        active_broad = str(active.get("target_broad_label", ""))
        comparison_broads = {active_broad}
        for row in read_tsv(
            broad_evidence_out / "broad_cell_type_current_member_questions.tsv.gz"
        ):
            challenger = str(row.get("challenger_broad_label", ""))
            if challenger:
                comparison_broads.add(challenger)
        for row in read_tsv(
            broad_evidence_out
            / "broad_cell_type_outside_recall_membership.tsv.gz"
        ):
            origin = str(row.get("current_broad_label", ""))
            if origin:
                comparison_broads.add(origin)
        pseudobulk_command = [
            args.rscript,
            str(paths["export_broad_cell_type_review_pseudobulk.R"]),
            "--rds", str(paths["selected_input"]),
            "--count-cache", str(count_export_document["full_count_matrix"]),
            "--raw-count-assay", str(count_export_document["raw_count_assay"]),
            "--membership", str(current_membership.resolve()),
            "--recall-membership", str(
                broad_evidence_out
                / "broad_cell_type_outside_recall_membership.tsv.gz"
            ),
            "--active-broad-label", active_broad,
            "--out", str(pseudobulk_out),
        ]
        for broad in sorted(comparison_broads):
            if broad:
                pseudobulk_command.extend(["--comparison-broad-label", broad])
        run(
            pseudobulk_command,
            output
            / f"logs/07_catalog_wide_pseudobulk_round_{review_round:02d}.log",
        )
        pseudobulk_summary_out = evidence_root / "02_pseudobulk_summary"
        run([
            sys.executable,
            str(paths["summarize_broad_cell_type_review_pseudobulk.py"]),
            "--pseudobulk", str(
                pseudobulk_out / "broad_cell_type_review_pseudobulk.tsv.gz"
            ),
            "--pseudobulk-manifest", str(
                pseudobulk_out
                / "broad_cell_type_review_pseudobulk_manifest.json"
            ),
            "--broad-evidence-manifest", str(broad_evidence_manifest),
            "--marker-manifest", str(marker_manifest),
            "--out", str(pseudobulk_summary_out),
        ], output / f"logs/07_catalog_wide_pseudobulk_summary_round_{review_round:02d}.log")
        pseudobulk_summary_manifest = (
            pseudobulk_summary_out
            / "broad_cell_type_review_pseudobulk_summary_manifest.json"
        )
        packet_out = evidence_root / "03_evidence_packets"
        packet_command = [
            sys.executable,
            str(paths["build_broad_cell_type_review_packet_index.py"]),
            "--review-manifest", str(last_review.resolve()),
            "--broad-evidence-manifest", str(broad_evidence_manifest),
            "--pseudobulk-summary-manifest", str(pseudobulk_summary_manifest),
            "--threshold-registry", str(paths["threshold_registry"]),
            "--active-review-id", active_review_id,
            "--out", str(packet_out),
        ]
        if biological_quality_path:
            packet_command.extend([
                "--biological-quality-review",
                str(biological_quality_path.resolve()),
            ])
        run(
            packet_command,
            output / f"logs/07_catalog_wide_packet_round_{review_round:02d}.log",
        )
        evidence_packet_manifest = (
            packet_out / "broad_cell_type_review_packet_manifest.json"
        )
        evidence_packet_manifests.append(evidence_packet_manifest)
        if not supplied_decisions:
            return {
                "status": "REVIEW_REQUIRED", "membership": current_membership,
                "review_manifest": last_review,
                "decision_validations": prior_validations,
                "apply_manifests": apply_manifests,
                "evidence_packet_manifests": evidence_packet_manifests,
                "manual_adjudications": manual_adjudications,
                "static_evidence_manifest": static_evidence_manifest,
                "review_state": review_state_path,
                "active_cell_type": active.get("target_broad_label", ""),
                "required_progress_message": active.get(
                    "required_progress_message", ""
                ),
            }
        decision_out = output / f"07_catalog_wide_decision_round_{review_round:02d}"
        run([
            sys.executable,
            str(paths["validate_catalog_wide_lineage_review_decisions.py"]),
            "--review-manifest", str(last_review),
            "--review-state", str(review_state_path),
            "--evidence-packet-manifest", str(evidence_packet_manifest),
            "--decisions", str(supplied_decisions[0].resolve()),
            "--out", str(decision_out),
        ], output / f"logs/07_catalog_wide_decision_round_{review_round:02d}.log", allowed_codes=(0, 2))
        validation_path = (
            decision_out / "catalog_wide_lineage_decision_validation.json"
        )
        validation = json.loads(validation_path.read_text(encoding="utf-8"))
        if validation.get("status") != "PASS":
            return {
                "status": "ITERATION_REQUIRED", "membership": current_membership,
                "review_manifest": last_review,
                "decision_validations": prior_validations + [validation_path],
                "apply_manifests": apply_manifests,
                "evidence_packet_manifests": evidence_packet_manifests,
                "manual_adjudications": manual_adjudications,
                "static_evidence_manifest": static_evidence_manifest,
                "review_state": review_state_path,
                "active_cell_type": active.get("target_broad_label", ""),
            }
        apply_out = output / f"07_catalog_wide_apply_round_{review_round:02d}"
        apply_authority = stage_authority(
            "atlas_and_completeness_review", args.contract, paths, apply_out,
            post_atlas_membership=current_membership,
            catalog_review_manifest=last_review,
            cell_type_review_state=review_state_path,
            catalog_decision_validation=validation_path,
            candidate_catalog=paths["catalog"],
            context_evidence=paths.get("context_evidence"),
        )
        apply_command = [
            sys.executable, str(paths["apply_catalog_wide_lineage_review.py"]),
            "--contract", str(args.contract),
            "--stage-authority", str(apply_authority),
            "--membership", str(current_membership.resolve()),
            "--review-manifest", str(last_review),
            "--review-state", str(review_state_path),
            "--decision-validation", str(validation_path),
            "--catalog", str(paths["catalog"]),
            "--out", str(apply_out),
        ]
        if "context_evidence" in paths:
            apply_command.extend([
                "--context-evidence", str(paths["context_evidence"])
            ])
        run(
            apply_command,
            output / f"logs/07_catalog_wide_apply_round_{review_round:02d}.log",
        )
        apply_manifest = apply_out / "catalog_wide_lineage_review_apply_manifest.json"
        apply_doc = json.loads(apply_manifest.read_text(encoding="utf-8"))
        current_membership = Path(str(apply_doc["membership"]["path"]))
        prior_validations.append(validation_path)
        apply_manifests.append(apply_manifest)
        return {
            "status": "REVIEW_STEP_APPLIED",
            "membership": current_membership,
            "review_manifest": last_review,
            "decision_validations": prior_validations,
            "apply_manifests": apply_manifests,
            "evidence_packet_manifests": evidence_packet_manifests,
            "manual_adjudications": manual_adjudications,
            "static_evidence_manifest": static_evidence_manifest,
            "review_state": review_state_path,
            "active_cell_type": active.get("target_broad_label", ""),
            "required_progress_message": (
                "当前类型已关闭或修订；重新运行本阶段以激活下一个细胞类型。"
            ),
        }
    raise RuntimeError("catalog-wide review loop exceeded its bounded policy")


def resume_sequential_cell_type_review(
    args, contract: dict, paths: dict[str, Path], prerequisite: dict,
    cluster_evidence_paths: list[Path], fine_proposal_paths: list[Path],
    state_proposal_paths: list[Path], unmodeled_manifest_paths: list[Path],
    observation_score_paths: list[Path],
    local_subset_validation_paths: list[Path],
    local_remainder_audit_paths: list[Path],
) -> dict:
    """Resume only the final serial biological review; never rerun Atlas or clustering."""
    resume_path = args.resume_review_manifest.resolve()
    resume = json.loads(resume_path.read_text(encoding="utf-8"))
    resume_contract_record = resume.get("annotation_contract", {})
    same_contract = resume_contract_record.get("sha256") == sha256(args.contract)
    compatible_contract = False
    compatible_scope = ""
    if not same_contract:
        for record in contract.get("resume_compatible_contracts", []):
            path = Path(str((record or {}).get("path", "")))
            if (
                path.is_file()
                and (record or {}).get("sha256") == sha256(path)
                and resume_contract_record.get("sha256") == sha256(path)
                and Path(str(resume_contract_record.get("path", ""))).resolve()
                == path.resolve()
            ):
                compatible_contract = True
                compatible_scope = str(
                    (record or {}).get("compatibility_scope", "")
                )
                break
    serial_pause = (
        resume.get("status") == "REVIEW_REQUIRED"
        and resume.get("pause_reason") in {
            "single_active_cell_type_review",
            "single_cell_type_review_step_applied",
            "manual_biological_adjudication_required",
        }
    )
    post_completeness_iteration = (
        resume.get("status") == "ITERATION_REQUIRED"
        and resume.get("catalog_wide_lineage_review_status") == "PASS"
        and isinstance(resume.get("completeness"), dict)
        and isinstance(resume.get("membership"), dict)
        and isinstance(resume.get("membership_transform_chain"), dict)
    )
    if (
        resume.get("phase") != "atlas_and_completeness_review"
        or not (serial_pause or post_completeness_iteration)
        or not (same_contract or compatible_contract)
        or (
            compatible_contract
            and not (
                post_completeness_iteration
                or resume.get("pause_reason")
                == "manual_biological_adjudication_required"
                or (
                    compatible_scope
                    == "manual_biological_adjudication_provenance_resume"
                    and resume.get("pause_reason") in {
                        "single_active_cell_type_review",
                        "single_cell_type_review_step_applied",
                    }
                )
            )
        )
        or resume.get("prerequisite", {}).get("sha256")
        != sha256(args.prerequisite_manifest)
    ):
        raise RuntimeError(
            "resume manifest is not a canonical pause from this contract/prerequisite"
        )

    def one(record: dict, label: str) -> Path:
        path = Path(str((record or {}).get("path", "")))
        if not path.is_file() or (record or {}).get("sha256") != sha256(path):
            raise RuntimeError(f"resume {label} is missing or stale")
        return path

    current_membership = one(resume.get("membership", {}), "membership")
    transform_chain_path = one(
        resume.get("membership_transform_chain", {}),
        "membership transform chain",
    )
    transform_document = load_and_validate_chain(
        transform_chain_path, current_membership
    )
    transform_index = int(transform_document.get("transform_n", 0)) + 1
    prior_validations = bound_artifact_paths(
        resume.get("catalog_wide_decision_validations", []),
        "resume cell-type decision validation",
    )
    prior_applies = bound_artifact_paths(
        resume.get("catalog_wide_apply_manifests", []),
        "resume cell-type apply manifest",
    )
    prior_packets = bound_artifact_paths(
        resume.get("catalog_wide_evidence_packet_manifests", []),
        "resume cell-type evidence packet",
    )
    recovered_manual_adjudications = bound_artifact_paths(
        resume.get("catalog_wide_manual_biological_adjudications", []),
        "resume manual biological adjudication",
    )
    static_record = resume.get("catalog_wide_static_evidence")
    recovered_static_evidence = (
        one(static_record, "catalog review static evidence")
        if static_record else None
    )
    if resume.get("pause_reason") == "manual_biological_adjudication_required":
        # A rejected or still-pending decision file is diagnostic, not a
        # completed biological decision. Permit a corrected, packet-bound
        # submission without consuming the type's single formal conclusion.
        retained_validations: list[Path] = []
        for path in prior_validations:
            document = json.loads(path.read_text(encoding="utf-8"))
            if document.get("status") == "PASS":
                retained_validations.append(path)
            elif document.get("status") not in {"BLOCKED", "ITERATION_REQUIRED"}:
                raise RuntimeError("manual retry contains an invalid prior validation")
        prior_validations = retained_validations
        # Failed retries can append the same frozen packet more than once.
        # Retain one occurrence per canonical path, preferring the newest.
        deduplicated_packets: list[Path] = []
        seen_packet_paths: set[Path] = set()
        for path in reversed(prior_packets):
            resolved = path.resolve()
            if resolved not in seen_packet_paths:
                deduplicated_packets.append(path)
                seen_packet_paths.add(resolved)
        prior_packets = list(reversed(deduplicated_packets))
    previous_review = one(
        resume.get("catalog_wide_lineage_review", {}),
        "previous cell-type review",
    )
    previous_state_record = resume.get("cell_type_review_state")
    previous_state = (
        one(previous_state_record, "previous cell-type review state")
        if previous_state_record else None
    )
    for supplied, recovered, label in (
        (args.lineage_review_prior_validation, prior_validations, "validation"),
        (args.lineage_review_prior_apply, prior_applies, "apply"),
        (args.lineage_review_prior_packet, prior_packets, "packet"),
    ):
        if supplied and [path.resolve() for path in supplied] != [
            path.resolve() for path in recovered
        ]:
            raise RuntimeError(f"resume {label} list differs from the pause manifest")
    args.lineage_review_prior_validation = prior_validations
    args.lineage_review_prior_apply = prior_applies
    args.lineage_review_prior_packet = prior_packets
    if recovered_static_evidence is not None:
        supplied_static = getattr(
            args, "lineage_review_static_evidence_manifest", None
        )
        if supplied_static and supplied_static.resolve() != recovered_static_evidence.resolve():
            raise RuntimeError(
                "resume static catalog evidence differs from the pause manifest"
            )
        args.lineage_review_static_evidence_manifest = recovered_static_evidence
    supplied_manual_adjudications = list(
        args.lineage_review_manual_adjudication or []
    )
    if recovered_manual_adjudications:
        if supplied_manual_adjudications and [
            path.resolve() for path in supplied_manual_adjudications
        ] != [path.resolve() for path in recovered_manual_adjudications]:
            raise RuntimeError(
                "resume manual adjudication list differs from the pause manifest"
            )
        args.lineage_review_manual_adjudication = recovered_manual_adjudications
    else:
        if supplied_manual_adjudications and resume.get("pause_reason") \
                != "manual_biological_adjudication_required":
            state_manual = []
            if previous_state is not None:
                state_document = json.loads(
                    previous_state.read_text(encoding="utf-8")
                )
                state_manual = bound_artifact_paths(
                    state_document.get("manual_biological_adjudications", []),
                    "resume state manual biological adjudication",
                )
            if [path.resolve() for path in supplied_manual_adjudications] != [
                path.resolve() for path in state_manual
            ]:
                raise RuntimeError(
                    "a new manual adjudication requires the canonical manual pause"
                )
        args.lineage_review_manual_adjudication = supplied_manual_adjudications
    args.lineage_review_previous_review = previous_review
    args.lineage_review_previous_state = previous_state

    profile = json.loads(paths["profile"].read_text(encoding="utf-8"))
    quality_required = bool(
        profile.get("biological_quality_endpoints", {}).get("required", False)
    )
    prior_apply_n = len(prior_applies)
    supplied_decisions = list(args.lineage_review_decisions or [])
    if supplied_decisions:
        if len(supplied_decisions) != 1:
            raise RuntimeError(
                "one controller invocation may decide only one active cell type"
            )
        if previous_state is None or not prior_packets:
            raise RuntimeError(
                "a resumed decision requires the frozen active state and evidence packet"
            )
        state_document = json.loads(previous_state.read_text(encoding="utf-8"))
        active = state_document.get("active_cell_type_review") or {}
        active_review_id = str(active.get("review_id", ""))
        active_cell_type = str(active.get("target_broad_label", ""))
        if (
            state_document.get("status") != "REVIEW_REQUIRED"
            or not active_review_id
            or not active_cell_type
        ):
            raise RuntimeError("resumed decision lacks one active cell type")
        review_round = len(prior_validations) + 1
        decision_out = (
            args.out.resolve()
            / f"07_catalog_wide_decision_round_{review_round:02d}"
        )
        run([
            sys.executable,
            str(paths["validate_catalog_wide_lineage_review_decisions.py"]),
            "--review-manifest", str(previous_review),
            "--review-state", str(previous_state),
            "--evidence-packet-manifest", str(prior_packets[-1]),
            "--decisions", str(supplied_decisions[0].resolve()),
            "--out", str(decision_out),
        ], args.out.resolve() / f"logs/07_catalog_wide_decision_round_{review_round:02d}.log",
            allowed_codes=(0, 2))
        validation_path = (
            decision_out / "catalog_wide_lineage_decision_validation.json"
        )
        validation = json.loads(validation_path.read_text(encoding="utf-8"))
        if validation.get("status") != "PASS":
            review = {
                "status": "ITERATION_REQUIRED",
                "membership": current_membership,
                "review_manifest": previous_review,
                "decision_validations": prior_validations + [validation_path],
                "apply_manifests": prior_applies,
                "evidence_packet_manifests": prior_packets,
                "manual_adjudications": list(
                    args.lineage_review_manual_adjudication or []
                ),
                "review_state": previous_state,
                "active_cell_type": active_cell_type,
                "required_progress_message": "",
            }
        else:
            apply_out = (
                args.out.resolve()
                / f"07_catalog_wide_apply_round_{review_round:02d}"
            )
            apply_authority = stage_authority(
                "atlas_and_completeness_review", args.contract, paths, apply_out,
                post_atlas_membership=current_membership,
                catalog_review_manifest=previous_review,
                cell_type_review_state=previous_state,
                catalog_decision_validation=validation_path,
                candidate_catalog=paths["catalog"],
                context_evidence=paths.get("context_evidence"),
            )
            apply_command = [
                sys.executable,
                str(paths["apply_catalog_wide_lineage_review.py"]),
                "--contract", str(args.contract),
                "--stage-authority", str(apply_authority),
                "--membership", str(current_membership.resolve()),
                "--review-manifest", str(previous_review),
                "--review-state", str(previous_state),
                "--decision-validation", str(validation_path),
                "--catalog", str(paths["catalog"]),
                "--out", str(apply_out),
            ]
            if "context_evidence" in paths:
                apply_command.extend([
                    "--context-evidence", str(paths["context_evidence"])
                ])
            run(
                apply_command,
                args.out.resolve()
                / f"logs/07_catalog_wide_apply_round_{review_round:02d}.log",
            )
            apply_manifest = (
                apply_out / "catalog_wide_lineage_review_apply_manifest.json"
            )
            apply_document = json.loads(
                apply_manifest.read_text(encoding="utf-8")
            )
            review = {
                "status": "REVIEW_STEP_APPLIED",
                "membership": Path(str(apply_document["membership"]["path"])),
                "review_manifest": previous_review,
                "decision_validations": prior_validations + [validation_path],
                "apply_manifests": prior_applies + [apply_manifest],
                "evidence_packet_manifests": prior_packets,
                "manual_adjudications": list(
                    args.lineage_review_manual_adjudication or []
                ),
                "review_state": previous_state,
                "active_cell_type": active_cell_type,
                "required_progress_message": (
                    "当前类型已关闭或修订；重新运行本阶段以激活下一个细胞类型。"
                ),
            }
    else:
        review = run_catalog_wide_review_iterations(
            args, paths, args.out.resolve(), current_membership,
            observation_score_paths, cluster_evidence_paths, quality_required,
            membership_already_includes_prior_applies=True,
        )
    review_status = str(review["status"])
    current_membership = Path(review["membership"])
    apply_manifests = list(review["apply_manifests"])
    transform_root = args.out.resolve() / "membership_transform_chain"
    for apply_manifest_path in apply_manifests[prior_apply_n:]:
        apply_doc = json.loads(
            apply_manifest_path.read_text(encoding="utf-8")
        )
        if int(apply_doc.get("n_changed_observations", 0)) == 0:
            continue
        transform_chain_path = append_membership_transform(
            paths, transform_chain_path, "cell_type_review_patch",
            Path(str(apply_doc["source_membership"]["path"])),
            Path(str(apply_doc["membership"]["path"])),
            apply_manifest_path, transform_root, transform_index,
            args.out.resolve() / "logs",
            target_cell_type=str(apply_doc.get("active_cell_type", "")),
        )
        transform_index += 1

    common = {
        "phase": "atlas_and_completeness_review",
        "formal_membership_written": False,
        "membership": artifact(current_membership),
        "membership_transform_chain": artifact(transform_chain_path),
        "prerequisite": artifact(args.prerequisite_manifest),
        "resume_source": artifact(resume_path),
        "atlas_validation": resume.get("atlas_validation"),
        "post_merge_unresolved_review": resume.get(
            "post_merge_unresolved_review"
        ),
        "biological_quality_status": resume.get(
            "biological_quality_status", "NOT_REQUIRED"
        ),
        "biological_quality_review": resume.get(
            "biological_quality_review"
        ),
        "targeted_follicle_roi_repair": resume.get(
            "targeted_follicle_roi_repair"
        ),
        "catalog_wide_lineage_review_status": review_status,
        "catalog_wide_lineage_review": artifact(review["review_manifest"]),
        "catalog_wide_static_evidence": (
            artifact(Path(review["static_evidence_manifest"]))
            if review.get("static_evidence_manifest")
            else artifact(recovered_static_evidence)
            if recovered_static_evidence else None
        ),
        "cell_type_review_state": (
            artifact(Path(review["review_state"]))
            if review.get("review_state") else None
        ),
        "active_cell_type_review": review.get("active_cell_type", ""),
        "required_progress_message": review.get(
            "required_progress_message", ""
        ),
        "catalog_wide_decision_validations": [
            artifact(path) for path in review["decision_validations"]
        ],
        "catalog_wide_evidence_packet_manifests": [
            artifact(path) for path in review["evidence_packet_manifests"]
        ],
        "catalog_wide_manual_biological_adjudications": [
            artifact(path) for path in review.get("manual_adjudications", [])
        ],
        "catalog_wide_apply_manifests": [
            artifact(path) for path in apply_manifests
        ],
        "observation_scores": [artifact(path) for path in observation_score_paths],
        "cluster_evidence": [artifact(path) for path in cluster_evidence_paths],
        "fine_candidate_proposals": [artifact(path) for path in fine_proposal_paths],
        "state_annotation_proposals": [artifact(path) for path in state_proposal_paths],
    }
    if review_status != "PASS":
        pause_reason = {
            "REVIEW_REQUIRED": "single_active_cell_type_review",
            "REVIEW_STEP_APPLIED": "single_cell_type_review_step_applied",
            "ITERATION_REQUIRED": "manual_biological_adjudication_required",
        }.get(review_status, "single_active_cell_type_review")
        return {
            **common,
            "status": "REVIEW_REQUIRED",
            "pause_reason": pause_reason,
        }

    completeness_out = args.out.resolve() / "08_post_catalog_completeness"
    command = [
        sys.executable, str(paths["audit_post_merge_completeness.py"]),
        "--membership", str(current_membership),
        "--catalog", str(paths["catalog"]),
        "--membership-transform-chain", str(transform_chain_path),
        "--out", str(completeness_out),
    ]
    for path in cluster_evidence_paths:
        command.extend(["--cluster-evidence", str(path)])
    for path in fine_proposal_paths:
        command.extend(["--fine-audit", str(path)])
    for path in unmodeled_manifest_paths:
        command.extend(["--unmodeled", str(path)])
    if args.unmodeled_decisions:
        command.extend(["--unmodeled-review", str(args.unmodeled_decisions)])
    for path in local_subset_validation_paths:
        command.extend(["--local-subset-validation", str(path)])
    for path in local_remainder_audit_paths:
        command.extend(["--local-remainder-audit", str(path)])
    if "context_evidence" in paths:
        command.extend(["--context-evidence", str(paths["context_evidence"])])
    command.extend([
        "--catalog-wide-review-summary", str(review["review_manifest"])
    ])
    if quality_required:
        command.append("--defer-canonical-zero-to-biological-review")
    run(
        command, args.out.resolve() / "logs/08_post_catalog_completeness.log",
        allowed_codes=(0, 2),
    )
    completeness_manifest = (
        completeness_out / "post_merge_completeness_manifest.json"
    )
    completeness_status = str(json.loads(
        completeness_manifest.read_text(encoding="utf-8")
    ).get("status", ""))
    quality_status = "NOT_REQUIRED"
    quality_manifest = None
    if quality_required:
        quality_out = args.out.resolve() / "09_post_catalog_biological_quality"
        quality_command = [
            sys.executable,
            str(paths["validate_sheep_ovary_biological_quality.py"]),
            "--membership", str(current_membership),
            "--catalog", str(paths["catalog"]), "--out", str(quality_out),
        ]
        for path in observation_score_paths:
            quality_command.extend(["--scores", str(path)])
        if args.canonical_oocyte_review:
            quality_command.extend([
                "--canonical-oocyte-review", str(args.canonical_oocyte_review)
            ])
        for path in args.lineage_review_manual_adjudication:
            quality_command.extend([
                "--manual-biological-adjudication", str(path.resolve())
            ])
        run(
            quality_command,
            args.out.resolve() / "logs/09_post_catalog_biological_quality.log",
            allowed_codes=(0, 2),
        )
        quality_path = (
            quality_out / "sheep_ovary_biological_quality_review.json"
        )
        quality_status = str(json.loads(
            quality_path.read_text(encoding="utf-8")
        ).get("status", ""))
        quality_manifest = artifact(quality_path)
    unresolved_n = sum(
        not str(row.get("final_broad_label", ""))
        for row in read_tsv(current_membership)
    )
    return {
        **common,
        "status": (
            "PASS"
            if completeness_status == "PASS"
            and quality_status in {"PASS", "NOT_REQUIRED"}
            else "ITERATION_REQUIRED"
        ),
        "completeness": artifact(completeness_manifest),
        "biological_quality_status": quality_status,
        "biological_quality_review": quality_manifest,
        "n_unresolved_biological": unresolved_n,
        "unlabeled_broad_rescue_n": int(
            resume.get("unlabeled_broad_rescue_n", 0) or 0
        ),
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
    if args.resume_review_manifest:
        return resume_sequential_cell_type_review(
            args, contract, paths, prerequisite, cluster_evidence_paths,
            fine_proposal_paths, state_proposal_paths,
            unmodeled_manifest_paths, observation_score_paths,
            local_subset_validation_paths, local_remainder_audit_paths,
        )
    atlas_bundle = json.loads(paths["atlas_bundle"].read_text(encoding="utf-8"))
    if atlas_bundle.get("bundle_id") != ACTIVE_BUNDLE_ID:
        raise RuntimeError(
            "new Atlas routing requires the active split-wall GSE233801 bundle; "
            "the merged v1 bundle is resume-only"
        )
    output = args.out.resolve()
    authority = stage_authority(
        "atlas_and_completeness_review", args.contract, paths, output,
        frozen_broad=args.frozen_broad_membership,
        atlas_bundle=paths["atlas_bundle"],
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
        "--atlas-bundle-manifest", str(paths["atlas_bundle"]),
        "--catalog", str(paths["catalog"]),
        "--out", str(routing_out),
        *(
            ["--context-evidence", str(paths["context_evidence"])]
            if "context_evidence" in paths else []
        ),
    ], output / "logs/00_atlas_routing.log", allowed_codes=(0, 2))
    validation_path = output / "01_atlas_validation.json"
    validation_command = [
        sys.executable, str(paths["validate_global_atlas_v2.py"]),
        "--routing-manifest", str(routing_out / "atlas_state_routing_manifest.json"),
        "--out", str(validation_path),
    ]
    if args.atlas_decisions:
        validation_command.extend(["--decisions", str(args.atlas_decisions.resolve())])
    run(
        validation_command, output / "logs/01_atlas_validation.log",
        allowed_codes=(0, 2),
    )
    atlas_validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if atlas_validation.get("status") == "REVIEW_REQUIRED":
        routing_manifest_path = (
            routing_out / "atlas_state_routing_manifest.json"
        )
        routing_manifest = json.loads(
            routing_manifest_path.read_text(encoding="utf-8")
        )
        return {
            "status": "REVIEW_REQUIRED",
            "phase": "atlas_and_completeness_review",
            "pause_reason": "atlas_discrepancy_decisions_required",
            "formal_membership_written": False,
            "frozen_broad_membership": artifact(
                args.frozen_broad_membership
            ),
            "atlas_routing": artifact(routing_manifest_path),
            "atlas_validation": artifact(validation_path),
            "review_queue": routing_manifest.get("artifacts", {}).get(
                "review_queue"
            ),
            "resume_requires": (
                "rerun the same phase with decisions bound to this review queue"
            ),
        }
    if atlas_validation.get("status") != "PASS":
        raise RuntimeError("Atlas validation is blocked by invalid evidence or decisions")
    apply_out = output / "02_post_atlas_membership"
    run([
        sys.executable, str(paths["apply_post_merge_atlas_routing.py"]),
        "--frozen-broad", str(args.frozen_broad_membership.resolve()),
        "--routing", str(routing_out / "atlas_state_routing.tsv.gz"),
        "--atlas-validation", str(validation_path),
        "--catalog", str(paths["catalog"]),
        "--out", str(apply_out),
        *(
            ["--context-evidence", str(paths["context_evidence"])]
            if "context_evidence" in paths else []
        ),
    ], output / "logs/02_post_atlas_membership.log")
    applied = json.loads(
        (apply_out / "post_atlas_membership_manifest.json").read_text(encoding="utf-8")
    )
    post_membership = Path(applied["membership"]["path"])
    transform_root = output / "membership_transform_chain"
    transform_chain_path = initialize_membership_transform_chain(
        paths, args.frozen_broad_membership, transform_root, output / "logs"
    )
    transform_index = 1
    transform_chain_path = append_membership_transform(
        paths, transform_chain_path, "atlas_unlabeled_broad_rescue",
        args.frozen_broad_membership, post_membership,
        apply_out / "post_atlas_membership_manifest.json",
        transform_root, transform_index, output / "logs",
    )
    transform_index += 1
    profile = json.loads(paths["profile"].read_text(encoding="utf-8"))
    quality_required = bool(
        profile.get("biological_quality_endpoints", {}).get("required", False)
    )
    unresolved_review_authority = stage_authority(
        "atlas_and_completeness_review", args.contract, paths, output,
        post_atlas_membership=post_membership,
        observation_scores=observation_score_paths,
        cluster_evidence=cluster_evidence_paths,
        context_evidence=paths.get("context_evidence"),
    )
    unresolved_review_out = output / "03_post_merge_unresolved_review"
    pre_unresolved_membership = post_membership
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
    if "context_evidence" in paths:
        unresolved_review_command.extend([
            "--context-evidence", str(paths["context_evidence"])
        ])
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
    transform_chain_path = append_membership_transform(
        paths, transform_chain_path, "post_merge_unresolved_return",
        pre_unresolved_membership, post_membership,
        unresolved_review_manifest_path, transform_root, transform_index,
        output / "logs",
    )
    transform_index += 1
    completeness_out = output / "04_completeness"
    command = [
        sys.executable, str(paths["audit_post_merge_completeness.py"]),
        "--membership", str(post_membership),
        "--catalog", str(paths["catalog"]),
        "--membership-transform-chain", str(transform_chain_path),
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
    if "context_evidence" in paths:
        command.extend(["--context-evidence", str(paths["context_evidence"])])
    if quality_required:
        command.append("--defer-canonical-zero-to-biological-review")
    run(
        command, output / "logs/04_completeness.log",
        allowed_codes=(0, 2),
    )
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
        if args.canonical_oocyte_review:
            quality_command.extend([
                "--canonical-oocyte-review",
                str(args.canonical_oocyte_review.resolve()),
            ])
        for path in args.lineage_review_manual_adjudication:
            quality_command.extend([
                "--manual-biological-adjudication", str(path.resolve())
            ])
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
            pre_follicle_repair_membership = post_membership
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
                transform_chain_path = append_membership_transform(
                    paths, transform_chain_path, "follicle_roi_reconciliation",
                    pre_follicle_repair_membership, post_membership,
                    repair["repair_manifest"], transform_root, transform_index,
                    output / "logs",
                )
                transform_index += 1
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
                    "--membership-transform-chain", str(transform_chain_path),
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
                if "context_evidence" in paths:
                    command.extend([
                        "--context-evidence", str(paths["context_evidence"])
                    ])
                if quality_required:
                    command.append(
                        "--defer-canonical-zero-to-biological-review"
                    )
                run(
                    command, output / "logs/06_post_repair_completeness.log",
                    allowed_codes=(0, 2),
                )
                completeness_manifest = (
                    completeness_out / "post_merge_completeness_manifest.json"
                )
    catalog_review = run_catalog_wide_review_iterations(
        args, paths, output, post_membership,
        observation_score_paths, cluster_evidence_paths,
        quality_required,
    )
    catalog_review_status = str(catalog_review["status"])
    post_membership = Path(catalog_review["membership"])
    catalog_apply_manifests = list(catalog_review["apply_manifests"])
    for apply_manifest_path in catalog_apply_manifests:
        apply_doc = json.loads(
            apply_manifest_path.read_text(encoding="utf-8")
        )
        source_path = Path(str(apply_doc["source_membership"]["path"]))
        result_path = Path(str(apply_doc["membership"]["path"]))
        if int(apply_doc.get("n_changed_observations", 0)) == 0:
            continue
        transform_chain_path = append_membership_transform(
            paths, transform_chain_path, "cell_type_review_patch",
            source_path, result_path, apply_manifest_path, transform_root,
            transform_index, output / "logs",
            target_cell_type=str(apply_doc.get("active_cell_type", "")),
        )
        transform_index += 1
    if catalog_review_status != "PASS":
        review_state_path = catalog_review.get("review_state")
        pause_reason = {
            "REVIEW_REQUIRED": "single_active_cell_type_review",
            "REVIEW_STEP_APPLIED": "single_cell_type_review_step_applied",
            "ITERATION_REQUIRED": "manual_biological_adjudication_required",
        }.get(catalog_review_status, "single_active_cell_type_review")
        return {
            "status": "REVIEW_REQUIRED",
            "phase": "atlas_and_completeness_review",
            "pause_reason": pause_reason,
            "formal_membership_written": False,
            "membership": artifact(post_membership),
            "membership_transform_chain": artifact(transform_chain_path),
            "stage_authority": artifact(authority),
            "prerequisite": artifact(args.prerequisite_manifest),
            "atlas_validation": artifact(validation_path),
            "post_merge_unresolved_review": artifact(
                unresolved_review_manifest_path
            ),
            "biological_quality_status": quality_status,
            "biological_quality_review": quality_manifest,
            "targeted_follicle_roi_repair": follicle_repair_manifest,
            "catalog_wide_lineage_review_status": catalog_review_status,
            "catalog_wide_lineage_review": artifact(
                catalog_review["review_manifest"]
            ),
            "catalog_wide_static_evidence": artifact(
                Path(catalog_review["static_evidence_manifest"])
            ),
            "cell_type_review_state": (
                artifact(Path(review_state_path)) if review_state_path else None
            ),
            "active_cell_type_review": catalog_review.get(
                "active_cell_type", ""
            ),
            "required_progress_message": catalog_review.get(
                "required_progress_message", ""
            ),
            "catalog_wide_decision_validations": [
                artifact(path) for path in catalog_review["decision_validations"]
            ],
            "catalog_wide_evidence_packet_manifests": [
                artifact(path)
                for path in catalog_review["evidence_packet_manifests"]
            ],
            "catalog_wide_apply_manifests": [
                artifact(path) for path in catalog_apply_manifests
            ],
            "catalog_wide_manual_biological_adjudications": [
                artifact(path)
                for path in catalog_review.get("manual_adjudications", [])
            ],
            "observation_scores": [
                artifact(path) for path in observation_score_paths
            ],
            "cluster_evidence": [
                artifact(path) for path in cluster_evidence_paths
            ],
            "fine_candidate_proposals": [
                artifact(path) for path in fine_proposal_paths
            ],
            "state_annotation_proposals": [
                artifact(path) for path in state_proposal_paths
            ],
            "unlabeled_broad_rescue_n": applied["unlabeled_broad_rescue_n"],
        }
    membership_output = artifact(post_membership)
    post_rows = read_tsv(post_membership)
    n_unresolved_biological = sum(
        not str(row.get("final_broad_label", row.get("broad_label", "")))
        for row in post_rows
    )
    if catalog_apply_manifests:
        completeness_out = output / "08_post_catalog_completeness"
        command = [
            sys.executable, str(paths["audit_post_merge_completeness.py"]),
            "--membership", str(post_membership),
            "--catalog", str(paths["catalog"]),
            "--membership-transform-chain", str(transform_chain_path),
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
        if "context_evidence" in paths:
            command.extend(["--context-evidence", str(paths["context_evidence"])])
        command.extend([
            "--catalog-wide-review-summary",
            str(catalog_review["review_manifest"]),
        ])
        if quality_required:
            command.append("--defer-canonical-zero-to-biological-review")
        run(
            command, output / "logs/08_post_catalog_completeness.log",
            allowed_codes=(0, 2),
        )
        completeness_manifest = (
            completeness_out / "post_merge_completeness_manifest.json"
        )
        if quality_required:
            quality_out = output / "09_post_catalog_biological_quality"
            quality_command = [
                sys.executable,
                str(paths["validate_sheep_ovary_biological_quality.py"]),
                "--membership", str(post_membership.resolve()),
                "--catalog", str(paths["catalog"]),
                "--out", str(quality_out),
            ]
            for path in observation_score_paths:
                quality_command.extend(["--scores", str(path.resolve())])
            if args.canonical_oocyte_review:
                quality_command.extend([
                    "--canonical-oocyte-review",
                    str(args.canonical_oocyte_review.resolve()),
                ])
            for path in args.lineage_review_manual_adjudication:
                quality_command.extend([
                    "--manual-biological-adjudication", str(path.resolve())
                ])
            run(
                quality_command,
                output / "logs/09_post_catalog_biological_quality.log",
                allowed_codes=(0, 2),
            )
            quality_path = (
                quality_out / "sheep_ovary_biological_quality_review.json"
            )
            quality_doc = json.loads(quality_path.read_text(encoding="utf-8"))
            quality_status = str(quality_doc.get("status", ""))
            quality_manifest = artifact(quality_path)
    completeness_status = str(json.loads(
        completeness_manifest.read_text(encoding="utf-8")
    ).get("status", ""))
    return {
        "status": (
            "ITERATION_REQUIRED"
            if quality_status == "ITERATION_REQUIRED"
            or catalog_review_status == "ITERATION_REQUIRED"
            or completeness_status != "PASS"
            else "PASS"
        ),
        "phase": "atlas_and_completeness_review",
        "stage_authority": artifact(authority),
        "prerequisite": artifact(args.prerequisite_manifest),
        "formal_membership_written": False,
        "membership": membership_output,
        "membership_transform_chain": artifact(transform_chain_path),
        "atlas_validation": artifact(validation_path),
        "post_merge_unresolved_review": artifact(
            unresolved_review_manifest_path
        ),
        "completeness": artifact(completeness_manifest),
        "biological_quality_status": quality_status,
        "biological_quality_review": quality_manifest,
        "canonical_oocyte_review": (
            artifact(args.canonical_oocyte_review)
            if args.canonical_oocyte_review else None
        ),
        "targeted_follicle_roi_repair": follicle_repair_manifest,
        "catalog_wide_lineage_review_status": catalog_review_status,
        "catalog_wide_lineage_review": artifact(
            catalog_review["review_manifest"]
        ),
        "catalog_wide_static_evidence": artifact(
            Path(catalog_review["static_evidence_manifest"])
        ),
        "catalog_wide_decision_validations": [
            artifact(path) for path in catalog_review["decision_validations"]
        ],
        "catalog_wide_evidence_packet_manifests": [
            artifact(path)
            for path in catalog_review["evidence_packet_manifests"]
        ],
        "catalog_wide_apply_manifests": [
            artifact(path) for path in catalog_apply_manifests
        ],
        "catalog_wide_manual_biological_adjudications": [
            artifact(path)
            for path in catalog_review.get("manual_adjudications", [])
        ],
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
    transform_record = prerequisite.get("membership_transform_chain", {})
    transform_chain = Path(str(transform_record.get("path", "")))
    if (
        not transform_chain.is_file()
        or transform_record.get("sha256") != sha256(transform_chain)
    ):
        raise RuntimeError("final phase lacks the reviewed membership transform chain")
    chain_document = load_and_validate_chain(
        transform_chain, args.post_atlas_membership
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
        "--catalog", str(paths["catalog"]),
        "--threshold-registry", str(paths["threshold_registry"]),
        "--out", str(fine_out),
    ]
    for path in fine_proposal_paths:
        fine_command.extend(["--fine-audit", str(path.resolve())])
    if "context_evidence" in paths:
        fine_command.extend([
            "--context-evidence", str(paths["context_evidence"])
        ])
    run(fine_command, output / "logs/00_fine_materialization.log")
    fine_assignments = fine_out / "parent_locked_fine_assignments.tsv.gz"
    authority = stage_authority(
        "materialize_final_release", args.contract, paths, output,
        post_atlas_membership=args.post_atlas_membership,
        prerequisite_manifest=args.prerequisite_manifest,
        fine_assignments=fine_assignments,
        state_annotation_proposals=state_proposal_paths,
        context_evidence=paths.get("context_evidence"),
    )
    final_out = output / "01_final_release"
    final_command = [
        sys.executable, str(paths["materialize_final_release_v2_2.py"]),
        "--contract", str(args.contract), "--stage-authority", str(authority),
        "--post-atlas-membership", str(args.post_atlas_membership.resolve()),
        "--atlas-completeness-manifest", str(args.prerequisite_manifest.resolve()),
        "--catalog", str(paths["catalog"]),
        "--fine-assignments", str(fine_assignments),
        "--out", str(final_out),
    ]
    for path in state_proposal_paths:
        final_command.extend(["--state-proposals", str(path.resolve())])
    if "context_evidence" in paths:
        final_command.extend([
            "--context-evidence", str(paths["context_evidence"])
        ])
    run(final_command, output / "logs/01_final_release.log")
    final_manifest = json.loads(
        (final_out / "final_release_manifest.json").read_text(encoding="utf-8")
    )
    final_manifest_path = final_out / "final_release_manifest.json"
    final_status = str(final_manifest.get("status", ""))
    if final_status not in {"PASS", "PENDING_USER_REVIEW_HIGH_UNRESOLVED"}:
        raise RuntimeError("final materializer returned an invalid status")
    final_membership = Path(str(final_manifest["membership"]["path"]))
    final_transform_chain = append_membership_transform(
        paths, transform_chain, "final_release_materialization",
        args.post_atlas_membership, final_membership, final_manifest_path,
        output / "02_final_membership_transform",
        int(chain_document.get("transform_n", 0)) + 1,
        output / "logs",
    )
    deliverables_out = output / "03_final_deliverables"
    deliverables_command = [
        sys.executable, str(paths["materialize_final_deliverables.py"]),
        "--contract", str(args.contract),
        "--membership", str(final_membership),
        "--release-manifest", str(final_manifest_path),
        "--rscript", str(args.rscript),
        "--threads", str(max(1, int(args.writer_threads))),
        "--release-status", str(args.release_status),
        "--biological-context", str(args.biological_context),
        "--project-root", str(args.contract.resolve().parent.parent),
        "--out", str(deliverables_out),
    ]
    run(deliverables_command, output / "logs/03_final_deliverables.log")
    deliverables_manifest = (
        deliverables_out / "final_deliverables_checkpoint.json"
    )
    return {
        "status": final_status, "phase": "materialize_final_release",
        "stage_authority": artifact(authority),
        "prerequisite": artifact(args.prerequisite_manifest),
        "formal_membership_written": bool(final_manifest["release_ready"]),
        "review_candidate_membership_written": bool(
            final_manifest["review_candidate_membership_written"]
        ),
        "required_next_action": final_manifest["required_next_action"],
        "membership": final_manifest["membership"],
        "membership_transform_chain": artifact(final_transform_chain),
        "final_deliverables": artifact(deliverables_manifest),
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


def run_phase_runtime_preflight(
    args: argparse.Namespace, paths: dict[str, Path]
) -> Path:
    """Fail before phase work when the exact runtime or semantic reader is invalid."""
    out = args.out.resolve() / "00_phase_runtime_preflight"
    command = [
        sys.executable, str(paths["validate_phase_runtime.py"]),
        "--python", sys.executable,
        "--rscript", str(args.rscript),
        "--semantic-input", f"annotation_contract:json:{args.contract}",
        "--semantic-input", f"biological_profile:json:{paths['profile']}",
        "--semantic-input", f"candidate_catalog:json:{paths['catalog']}",
        "--semantic-input", f"analysis_set:tsv:{paths['analysis_set']}",
        "--out", str(out),
    ]
    for name in CANONICAL_SCRIPTS:
        if name.endswith(".py"):
            command.extend(["--python-script", str(paths[name])])
        elif name.endswith(".R"):
            command.extend(["--r-script", str(paths[name])])
    python_imports = {
        "whole_tissue_partition": ("numpy", "pandas", "scipy"),
        "cluster_cohort_recluster": ("numpy", "pandas", "scipy"),
        "local_mixed_subcluster_split": ("numpy", "pandas", "scipy"),
        "merge_and_freeze_broad": ("numpy", "pandas"),
        "atlas_and_completeness_review": (
            "numpy", "pandas", "scipy", "sklearn", "matplotlib",
        ),
        "materialize_final_release": ("numpy", "pandas"),
    }
    for module in python_imports[args.phase]:
        command.extend(["--python-import", module])
    r_packages = {
        "whole_tissue_partition": (
            "Seurat", "SeuratObject", "Matrix", "data.table", "jsonlite",
        ),
        "cluster_cohort_recluster": (
            "Seurat", "SeuratObject", "Matrix", "data.table", "glmGamPoi",
            "future", "future.apply",
        ),
        "local_mixed_subcluster_split": ("Matrix", "data.table"),
        "merge_and_freeze_broad": (),
        "atlas_and_completeness_review": (
            "Seurat", "Matrix", "data.table", "jsonlite",
        ),
        "materialize_final_release": (
            "Seurat", "SeuratObject", "Matrix", "data.table", "jsonlite",
            "ggplot2", "scattermore", "patchwork",
        ),
    }
    for package in r_packages[args.phase]:
        command.extend(["--r-package", package])
    if "context_evidence" in paths:
        context_kind = (
            "json" if paths["context_evidence"].suffix.lower() == ".json"
            else "tsv"
        )
        command.extend([
            "--semantic-input",
            f"candidate_context_evidence:{context_kind}:{paths['context_evidence']}",
        ])
    run(command, args.out.resolve() / "logs/00_phase_runtime_preflight.log")
    return out / "phase_runtime_preflight.json"


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
    cohort.add_argument("--allocated-memory-gb", type=float)
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
    cohort.add_argument(
        "--resolution-workers", type=int, default=1,
        help=(
            "Leiden-resolution fork count; default 1 because each worker can "
            "replicate the loaded Seurat/SCT carrier. Parallelize independent "
            "cohorts instead of resolutions for large inputs."
        ),
    )

    local = sub.add_parser("local_mixed_subcluster_split")
    add_common(local)
    local.add_argument("--trigger-manifest", required=True, type=Path)
    local.add_argument(
        "--workload-audit", required=True, type=Path,
        help=(
            "PASS audit_local_split_workload.py result built from every "
            "second-round cohort manifest before any P41 job is submitted"
        ),
    )
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
    atlas.add_argument(
        "--canonical-oocyte-review", type=Path,
        help=(
            "optional frozen label-blind canonical Oocyte adjudication whose "
            "exact released member set is revalidated after catalog review"
        ),
    )
    atlas.add_argument(
        "--lineage-review-decisions", type=Path, action="append", default=[],
        help=(
            "exactly one decision file for the currently active cell type; "
            "batch decisions are forbidden"
        ),
    )
    atlas.add_argument(
        "--lineage-review-prior-validation", type=Path,
        action="append", default=[],
        help="PASS single-cell-type validations from earlier invocations",
    )
    atlas.add_argument(
        "--lineage-review-prior-apply", type=Path,
        action="append", default=[],
        help="contiguous single-cell-type apply manifests from earlier invocations",
    )
    atlas.add_argument(
        "--lineage-review-prior-packet", type=Path,
        action="append", default=[],
        help="earlier single-active-type evidence packet manifests",
    )
    atlas.add_argument(
        "--lineage-review-manual-adjudication", type=Path,
        action="append", default=[],
        help=(
            "legacy/user-authorized non-mutating closure for one exact "
            "blocked review scope"
        ),
    )
    atlas.add_argument("--lineage-review-previous-review", type=Path)
    atlas.add_argument("--lineage-review-previous-state", type=Path)
    atlas.add_argument(
        "--lineage-review-static-evidence-manifest", type=Path,
        help=(
            "optional canonical raw-count/coordinate evidence cache; resume "
            "manifests recover it automatically so the large RDS is not reopened"
        ),
    )
    atlas.add_argument(
        "--resume-review-manifest", type=Path,
        help=(
            "canonical REVIEW_REQUIRED controller manifest; skips Atlas, "
            "unresolved rescue, ROI repair and clustering, then activates only "
            "the next single cell-type review"
        ),
    )
    atlas.add_argument("--resolution-workers", type=int, default=1)

    final = sub.add_parser("materialize_final_release")
    add_common(final)
    final.add_argument("--prerequisite-manifest", required=True, type=Path)
    final.add_argument("--post-atlas-membership", required=True, type=Path)
    final.add_argument(
        "--writer-threads", type=int,
        default=max(1, int(os.environ.get("LSB_DJOB_NUMPROC", "1"))),
    )
    final.add_argument("--biological-context", default="")
    final.add_argument(
        "--release-status",
        choices=["pending_user_review", "approved_final"],
        default="pending_user_review",
    )
    return parser.parse_args()


def execute(args: argparse.Namespace) -> tuple[dict, Path]:
    args.contract = args.contract.resolve()
    contract, paths = validate_contract(args.contract)
    runtime_preflight = run_phase_runtime_preflight(args, paths)
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
        "runtime_preflight": artifact(runtime_preflight),
        "schema_version": "2.2", "controller_version": "2.2.0",
        "annotation_contract": artifact(args.contract),
        "seed": args.seed, "finished_at_utc": datetime.now(timezone.utc).isoformat(),
    })
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = args.out / "lineage_controller_manifest.json"
    manifest.write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return result, manifest


def main() -> int:
    args = parse_args()
    try:
        result, manifest = execute(args)
    except Exception as exc:
        try:
            materialize_runtime_state(
                args.out.resolve(), args.phase, args.contract.resolve(), error=exc
            )
        except Exception:
            pass
        raise
    current_stage, next_action = materialize_runtime_state(
        args.out.resolve(), args.phase, args.contract.resolve(),
        result=result, controller_manifest=manifest,
    )
    result["runtime_state"] = {
        "current_stage": artifact(current_stage),
        "next_action": artifact(next_action),
    }
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
