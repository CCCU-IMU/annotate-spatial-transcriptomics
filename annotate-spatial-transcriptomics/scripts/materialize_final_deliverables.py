#!/usr/bin/env python3
"""Build the canonical, resumable RDS/plot/HTML delivery from final membership."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from evidence_schema_lib import sha256


SCRIPT_DIR = Path(__file__).resolve().parent


def semantic_hash(release: dict) -> str:
    value = str(release.get("membership", {}).get("semantic_sha256", ""))
    if len(value) != 64:
        raise RuntimeError("final release lacks a semantic membership hash")
    return value


def run(command: list[str], log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        result = subprocess.run(
            command, stdout=handle, stderr=subprocess.STDOUT,
            text=True, check=False,
        )
    if result.returncode != 0:
        raise RuntimeError(
            f"final deliverable step failed ({result.returncode}); inspect {log}"
        )


def record(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise RuntimeError(f"final deliverable is missing: {path}")
    return {
        "path": str(path.resolve()), "sha256": sha256(path),
        "n_bytes": path.stat().st_size,
    }


def bound_path(value: object, label: str) -> Path:
    if not isinstance(value, dict):
        raise RuntimeError(f"contract-bound {label} is absent")
    path = Path(str(value.get("path", ""))).resolve()
    if not path.is_file():
        raise RuntimeError(f"contract-bound {label} is missing")
    expected = str(value.get("sha256", ""))
    if len(expected) != 64 or expected != sha256(path):
        raise RuntimeError(f"contract-bound {label} is stale")
    return path


def valid(records: list[dict[str, object]]) -> bool:
    return bool(records) and all(
        (path := Path(str(item.get("path", "")))).is_file()
        and item.get("sha256") == sha256(path)
        for item in records
    )


def write_checkpoint(path: Path, document: dict[str, object]) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(
        json.dumps(document, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    temp.replace(path)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contract", required=True, type=Path)
    ap.add_argument("--membership", required=True, type=Path)
    ap.add_argument("--release-manifest", required=True, type=Path)
    ap.add_argument("--rscript", default="Rscript")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--biological-context", default="")
    ap.add_argument(
        "--project-root", type=Path,
        help="write small canonical latest/final pointers without copying the RDS",
    )
    ap.add_argument(
        "--release-status",
        choices=["pending_user_review", "approved_final"],
        default="pending_user_review",
    )
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    release = json.loads(args.release_manifest.read_text(encoding="utf-8"))
    membership_record = release.get("membership", {})
    if (
        Path(str(membership_record.get("path", ""))).resolve()
        != args.membership.resolve()
        or membership_record.get("sha256") != sha256(args.membership)
    ):
        raise RuntimeError("final deliverables do not bind final release membership")
    selected = bound_path(
        contract.get("selected_input_snapshot"), "selected RDS"
    )
    scope = contract.get("input_scope", {})
    excluded = bound_path(scope.get("excluded_initial_qc"), "excluded membership")
    coordinates = bound_path(scope.get("analysis_set"), "analysis membership")
    profile = bound_path(contract.get("biological_profile"), "biological profile")
    catalog = bound_path(contract.get("candidate_catalog"), "candidate catalog")
    input_audit = contract.get("query_input_audit")
    raw_count_assay = ""
    if input_audit:
        bound_path(input_audit, "query input audit")
        raw_count_assay = str(input_audit.get("raw_count_assay", "")).strip()
        if not raw_count_assay or raw_count_assay.upper() == "SCT":
            raise RuntimeError("query input audit lacks a non-SCT raw-count assay")
    if selected.suffix.lower() != ".rds":
        raise RuntimeError("canonical Seurat final delivery requires an RDS input")

    args.out.mkdir(parents=True, exist_ok=True)
    checkpoint_path = args.out / "final_deliverables_checkpoint.json"
    sem_hash = semantic_hash(release)
    checkpoint: dict[str, object]
    if checkpoint_path.is_file():
        checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        if (
            checkpoint.get("annotation_contract_sha256") != sha256(args.contract)
            or checkpoint.get("membership_semantic_sha256") != sem_hash
        ):
            raise RuntimeError(
                "existing final-deliverables checkpoint belongs to another contract or membership"
            )
    else:
        checkpoint = {
            "schema_version": "1.0", "status": "IN_PROGRESS",
            "stage": "canonical_final_deliverables",
            "annotation_contract_sha256": sha256(args.contract),
            "membership_semantic_sha256": sem_hash,
            "steps": {},
        }
        write_checkpoint(checkpoint_path, checkpoint)

    steps = checkpoint["steps"]

    def step(name: str, command: list[str], outputs: list[Path]) -> None:
        prior = steps.get(name, {}) if isinstance(steps, dict) else {}
        prior_outputs = prior.get("outputs", []) if isinstance(prior, dict) else []
        if prior.get("status") == "PASS" and valid(prior_outputs):
            return
        run(command, args.out / f"logs/{name}.log")
        records = [record(path) for path in outputs]
        steps[name] = {
            "status": "PASS", "outputs": records,
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        }
        write_checkpoint(checkpoint_path, checkpoint)

    evidence = args.out / "01_release_evidence"
    step("01_release_evidence", [
        sys.executable, str(SCRIPT_DIR / "build_release_evidence_tables.py"),
        "--membership", str(args.membership), "--catalog", str(catalog),
        "--profile", str(profile), "--out", str(evidence),
    ], [
        evidence / "release_evidence_tables_manifest.json",
        evidence / "canonical_marker_panel.tsv",
        evidence / "annotation_support_summary.tsv",
    ])

    annotated = args.out / "02_annotated_rds" / f"{contract['sample_id']}_annotated.rds"
    writer_manifest = args.out / "02_annotated_rds/writer_manifest.json"
    step("02_annotated_rds", [
        args.rscript, str(SCRIPT_DIR / "write_frozen_annotations_to_seurat.R"),
        "--rds", str(selected), "--membership", str(args.membership),
        "--excluded-membership", str(excluded), "--out", str(annotated),
        "--semantic-hash", sem_hash, "--sample-id", str(contract["sample_id"]),
        "--threads", str(max(1, args.threads)),
        "--writer-manifest", str(writer_manifest),
    ], [annotated, writer_manifest])

    maps = args.out / "03_annotation_maps"
    step("03_annotation_maps", [
        args.rscript, str(SCRIPT_DIR / "build_annotation_maps.R"),
        "--rds", str(selected), "--metadata", str(args.membership),
        "--coordinates", str(coordinates), "--skip-umap",
        "--cell-id-col", "cell_id", "--final-cell-type-col", "final_cell_type",
        "--out", str(maps),
    ], [
        maps / "figures/final_cell_type_spatial.png",
        maps / "tables/spatial_node_asset_index.tsv",
    ])

    dotplots = args.out / "04_marker_dotplots"
    dotplot_command = [
        args.rscript, str(SCRIPT_DIR / "build_marker_dotplots.R"),
        "--rds", str(selected), "--metadata", str(args.membership),
        "--markers", str(evidence / "canonical_marker_panel.tsv"),
        "--cell-id-col", "cell_id", "--final-cell-type-col", "final_cell_type",
        "--out", str(dotplots),
    ]
    if raw_count_assay:
        dotplot_command.extend(["--count-assay", raw_count_assay])
    step("04_marker_dotplots", dotplot_command, [
        dotplots / "marker_dotplot_asset_index.tsv"
    ])

    genes = args.out / "05_spatial_marker_panels"
    step("05_spatial_marker_panels", [
        args.rscript, str(SCRIPT_DIR / "build_spatial_gene_maps.R"),
        "--rds", str(selected), "--coordinates", str(coordinates),
        "--markers", str(evidence / "canonical_marker_panel.tsv"),
        "--cell-id-col", "cell_id", "--level", "cell_type",
        "--expected-observations", str(release["n_analysis_set"]),
        "--out", str(genes),
    ], [genes / "tables/spatial_gene_group_asset_index.tsv"])

    deg = args.out / "06_final_deg"
    step("06_final_deg", [
        args.rscript, str(SCRIPT_DIR / "run_final_label_deg.R"),
        "--rds", str(selected), "--metadata", str(args.membership),
        "--cell-id-col", "cell_id", "--final-cell-type-col", "final_cell_type",
        "--seed", str(contract["canonical_lineage_controller"]["random_seed"]),
        "--out", str(deg),
    ], [deg / "tables/cell_type_DEG_one_vs_rest_all.tsv"])

    prerequisite_path = bound_path(
        release.get("atlas_completeness_review"),
        "final Atlas/completeness review",
    )
    prerequisite = json.loads(prerequisite_path.read_text(encoding="utf-8"))
    quality = prerequisite.get("biological_quality_review", {})
    atlas = prerequisite.get("atlas_validation", {})
    report = args.out / "07_report" / f"{contract['sample_id']}_annotation_report.html"
    command = [
        sys.executable, str(SCRIPT_DIR / "build_frozen_review_report.py"),
        "--sample-id", str(contract["sample_id"]),
        "--biological-context", args.biological_context,
        "--observation-unit", str(contract["observation_unit"]),
        "--release-status", args.release_status,
        "--membership", str(args.membership),
        "--release-manifest", str(args.release_manifest),
        "--support", str(evidence / "annotation_support_summary.tsv"),
        "--maps-index", str(maps / "tables/spatial_node_asset_index.tsv"),
        "--maps-dir", str(maps),
        "--canonical-dotplots", str(dotplots / "marker_dotplot_asset_index.tsv"),
        "--canonical-gene-panels", str(genes / "tables/spatial_gene_group_asset_index.tsv"),
        "--cell-type-deg", str(deg / "tables/cell_type_DEG_one_vs_rest_all.tsv"),
        "--annotated-rds", str(annotated), "--out", str(report),
    ]
    if isinstance(quality, dict) and quality.get("path"):
        command.extend(["--biological-quality-review", str(quality["path"])])
    if isinstance(atlas, dict) and atlas.get("path"):
        command.extend(["--atlas-review", str(atlas["path"])])
    step("07_report", command, [report])

    checkpoint["status"] = "PASS"
    checkpoint["public_annotation_column"] = "final_cell_type"
    checkpoint["annotated_rds"] = record(annotated)
    checkpoint["report"] = record(report)
    checkpoint["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
    write_checkpoint(checkpoint_path, checkpoint)
    if args.project_root:
        project_root = args.project_root.resolve()
        pointer = {
            "schema_version": "1.0",
            "status": "PASS",
            "artifact_role": "canonical_final_release_pointer",
            "sample_id": str(contract["sample_id"]),
            "release_status": args.release_status,
            "public_annotation_column": "final_cell_type",
            "membership_semantic_sha256": sem_hash,
            "annotation_contract": record(args.contract),
            "final_release_manifest": record(args.release_manifest),
            "final_membership": record(args.membership),
            "final_deliverables_checkpoint": record(checkpoint_path),
            "annotated_rds": checkpoint["annotated_rds"],
            "report": checkpoint["report"],
            "completed_at_utc": checkpoint["completed_at_utc"],
        }
        for destination in (
            project_root / "deliverables/latest.json",
            project_root / "state/final_release_pointer.json",
        ):
            destination.parent.mkdir(parents=True, exist_ok=True)
            write_checkpoint(destination, pointer)
    print(json.dumps(checkpoint, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
