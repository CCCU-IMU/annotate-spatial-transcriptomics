#!/usr/bin/env python3
"""Create one canonical project, frozen SCT+BANKSY boundary and v2 contract."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[1]
SCRIPTS = PACKAGE / "scripts"
PROFILES = PACKAGE / "references/profiles"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(command: list[str], log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    log.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(
            f"canonical bootstrap command failed ({result.returncode}); inspect {log}"
        )


def append_snapshot(project_root: Path, sample: str, rds: Path) -> str:
    registry = project_root / "state/input_snapshot_registry.tsv"
    snapshot_id = f"{sample}_sct_banksy_input_v001"
    with registry.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    if any(row.get("snapshot_id") == snapshot_id for row in rows):
        raise RuntimeError("canonical bootstrap snapshot ID already exists")
    rows.append({
        "snapshot_id": snapshot_id,
        "sample_id": sample,
        "path": str(rds.resolve()),
        "kind": "seurat_sct_banksy_preprocessed_rds",
        "size_bytes": str(rds.stat().st_size),
        "sha256": sha256(rds),
        "status": "frozen",
        "created_at": datetime.now(timezone.utc).isoformat(),
    })
    with registry.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    return snapshot_id


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--rds", required=True, type=Path)
    ap.add_argument("--project-root", required=True, type=Path)
    ap.add_argument(
        "--observation-unit",
        choices=["cell", "nucleus", "spot", "cellbin"],
        default="cellbin",
    )
    ap.add_argument("--cluster-mapping", type=Path)
    ap.add_argument("--preprocess-manifest", type=Path)
    ap.add_argument("--context-evidence", type=Path)
    ap.add_argument("--rscript", default="Rscript")
    ap.add_argument("--seed", type=int, default=2200)
    args = ap.parse_args()

    rds = args.rds.resolve()
    root = args.project_root.resolve()
    if not rds.is_file():
        raise SystemExit("SCT+BANKSY RDS is missing")
    if root.exists() and any(root.iterdir()):
        raise SystemExit("project root already exists and is not empty")
    root.mkdir(parents=True, exist_ok=True)
    logs = root / "logs/bootstrap"

    run([
        sys.executable, str(SCRIPTS / "init_annotation_project.py"),
        "--sample", args.sample,
        "--input-root", str(rds.parent),
        "--project-root", str(root),
        "--modality", "spatial",
        "--observation-unit", args.observation_unit,
        "--strategy-preset", "sheep_ovary_same_batch_rfirst",
    ], logs / "01_init_project.log")

    freeze_root = root / "input/frozen_sct_banksy"
    freeze_root.mkdir(parents=True, exist_ok=True)
    digest = sha256(rds)
    sha_record = freeze_root / "input_rds.sha256"
    sha_record.write_text(f"{digest}  {rds}\n", encoding="utf-8")
    run([
        args.rscript, str(SCRIPTS / "freeze_sct_banksy_input.R"),
        str(rds), str(sha_record),
        str(args.cluster_mapping.resolve()) if args.cluster_mapping else "AUTO",
        (
            str(args.preprocess_manifest.resolve())
            if args.preprocess_manifest else "AUTO"
        ),
        args.sample, str(freeze_root),
    ], logs / "02_freeze_sct_banksy.log")

    snapshot_id = append_snapshot(root, args.sample, rds)
    grid_manifest = freeze_root / "whole_tissue_grid.json"
    grid = json.loads(grid_manifest.read_text(encoding="utf-8"))
    resolutions = ",".join(str(value) for value in grid["candidate_resolutions"])
    contract = root / "config/annotation_contract.json"
    command = [
        sys.executable, str(SCRIPTS / "build_annotation_contract_v2.py"),
        str(root),
        "--workflow-profile", str(PROFILES / "sheep_ovary_rfirst_profile.json"),
        "--biological-profile", str(PROFILES / "sheep_ovary.json"),
        "--candidate-catalog",
        str(PROFILES / "sheep_ovary_candidate_lineage_catalog.json"),
        "--snapshot-id", snapshot_id,
        "--analysis-membership", str(freeze_root / "analysis_membership.tsv.gz"),
        "--excluded-initial-qc", str(freeze_root / "excluded_initial_qc.tsv.gz"),
        "--input-audit-manifest", str(freeze_root / "input_audit_manifest.json"),
        "--whole-tissue-method", "BANKSY",
        "--whole-tissue-grid", resolutions,
        "--whole-tissue-grid-artifact", str(grid_manifest),
        "--whole-tissue-partitions", str(freeze_root / "partition_grid.tsv.gz"),
        "--grid-source", "bound_upstream_input",
        "--contract-id", "bootstrap_v2_5",
        "--seed", str(args.seed),
        "--out", str(contract),
    ]
    if args.context_evidence:
        command.extend(["--context-evidence", str(args.context_evidence.resolve())])
    run(command, logs / "03_build_contract.log")
    validation_root = root / "provenance/bootstrap_contract_validation"
    run([
        sys.executable, str(SCRIPTS / "validate_annotation_contract_v2.py"),
        str(contract), "--out", str(validation_root),
    ], logs / "04_validate_contract.log")

    manifest = {
        "schema_version": "1.0",
        "status": "PASS",
        "stage": "canonical_sct_banksy_bootstrap",
        "sample_id": args.sample,
        "input_rds": {"path": str(rds), "sha256": digest},
        "input_audit": {
            "path": str((freeze_root / "input_audit_manifest.json").resolve()),
            "sha256": sha256(freeze_root / "input_audit_manifest.json"),
        },
        "analysis_membership": {
            "path": str((freeze_root / "analysis_membership.tsv.gz").resolve()),
            "sha256": sha256(freeze_root / "analysis_membership.tsv.gz"),
        },
        "whole_tissue_partition_grid": {
            "path": str((freeze_root / "partition_grid.tsv.gz").resolve()),
            "sha256": sha256(freeze_root / "partition_grid.tsv.gz"),
        },
        "annotation_contract": {
            "path": str(contract.resolve()), "sha256": sha256(contract),
        },
        "next_action": "run whole_tissue_partition with the frozen partition grid",
    }
    manifest_path = root / "provenance/bootstrap_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
