#!/usr/bin/env python3
"""Register one low-priority assisted route under the direct-lineage controller."""

from __future__ import annotations

import argparse
import csv
import fcntl
import gzip
import hashlib
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path


ROUTE_CLASSES = {
    "targeted_rctd_review",
    "residual_qc_atlas_review",
    "contextual_reference_review",
    "explicit_technical_retention",
}
TERMINAL = {"validated_done", "not_applicable_reviewed"}


def truth(value: object) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "pass", "passed"}


def path_at(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def membership_set(path: Path) -> set[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    ids = [row.get("cell_id", "").strip() for row in rows]
    if not rows or "cell_id" not in rows[0] or any(not item for item in ids) or len(ids) != len(set(ids)):
        raise ValueError("membership is not a unique nonempty cell_id table")
    return set(ids)


def membership_count(path: Path) -> int:
    return len(membership_set(path))


def validate_terminal(root: Path, record: dict) -> list[str]:
    errors: list[str] = []
    if record.get("route_class") not in ROUTE_CLASSES:
        errors.append(f"invalid route_class: {record.get('route_class', '')}")
    if record.get("status") not in TERMINAL:
        errors.append(f"invalid terminal status: {record.get('status', '')}")
    applicability = record.get("applicability", "")
    if applicability == "not_applicable":
        if record.get("status") != "not_applicable_reviewed":
            errors.append("not_applicable route must close as not_applicable_reviewed")
        if len(str(record.get("applicability_rationale", "")).strip()) < 20:
            errors.append("not_applicable rationale is too short")
        return errors
    if applicability != "applicable":
        errors.append(f"invalid applicability: {applicability}")
        return errors
    for field in ("validation_artifact", "outcome_artifact"):
        value = str(record.get(field, ""))
        if not value or not path_at(root, value).is_file():
            errors.append(f"missing {field}")
    membership_value = str(record.get("query_membership_artifact", ""))
    membership = path_at(root, membership_value) if membership_value else Path()
    if not membership_value or not membership.is_file():
        errors.append("missing query_membership_artifact")
    else:
        expected = str(record.get("query_membership_sha256", ""))
        if not expected or sha256(membership) != expected:
            errors.append("query_membership_sha256 mismatch")
        try:
            if membership_count(membership) != int(float(record.get("n_query", 0) or 0)):
                errors.append("query membership count differs from n_query")
        except (ValueError, TypeError):
            errors.append("invalid query membership or n_query")
    route = record.get("route_class")
    if route == "targeted_rctd_review":
        try:
            high = int(float(record.get("rctd_high_n", 0) or 0))
            moderate = int(float(record.get("rctd_moderate_n", 0) or 0))
            low = int(float(record.get("rctd_low_n", 0) or 0))
            fine = int(float(record.get("rctd_fine_return_n", 0) or 0))
            broad = int(float(record.get("rctd_broad_return_n", 0) or 0))
            rerouted = int(float(record.get("n_rerouted", 0) or 0))
            n_query = int(float(record.get("n_query", 0) or 0))
            if min(high, moderate, low, fine, broad, rerouted) < 0:
                errors.append("RCTD counts must be nonnegative")
            if high + moderate + low != n_query:
                errors.append("RCTD tiers do not partition n_query")
            if fine > high or fine + broad > high + moderate:
                errors.append("RCTD return confidence exceeds its tier evidence")
            if fine and not truth(record.get("independent_fine_evidence")):
                errors.append("RCTD fine returns require independent evidence")
            if fine + broad + rerouted != n_query or low > rerouted:
                errors.append("RCTD returns plus QC reroute do not partition n_query")
            if truth(record.get("fine_anchor_eligible")):
                errors.append("RCTD-assisted output cannot become a fine anchor")
        except (ValueError, TypeError):
            errors.append("invalid RCTD counts")
    if route == "residual_qc_atlas_review":
        if record.get("source_state") != "residual_qc_holdout_after_all_cohorts":
            errors.append("Atlas source must be the final residual QC holdout")
        if record.get("calibration_origin") != "query_like_heldout_current_query_anchors":
            errors.append("Atlas calibration origin is invalid")
        if not truth(record.get("depth_matched_validation")) or not truth(record.get("observed_density_spatial_prior")):
            errors.append("Atlas requires depth-matched validation and observed-density spatial evidence")
        manifest = str(record.get("calibration_manifest", ""))
        if not manifest or not path_at(root, manifest).is_file():
            errors.append("missing calibration_manifest")
        try:
            n_query = int(float(record.get("n_query", 0) or 0))
            broad = int(float(record.get("n_defined_broad_only", 0) or 0))
            qc = int(float(record.get("n_qc_retained", 0) or 0))
            fine = int(float(record.get("n_defined_fine", 0) or 0))
            if broad + qc != n_query or fine != 0:
                errors.append("Atlas broad-only returns and QC rejects do not partition n_query")
            if truth(record.get("fine_anchor_eligible")):
                errors.append("Atlas rescue cannot become a fine anchor")
        except (ValueError, TypeError):
            errors.append("invalid Atlas outcome counts")
        if record.get("atlas_confidence_enum") != "high|moderate_only|low_reject":
            errors.append("Atlas confidence enum must be high|moderate_only|low_reject")
        partitions: dict[str, set[str]] = {}
        for prefix, expected_n in (("accepted_membership", record.get("n_defined_broad_only", 0)), ("rejected_membership", record.get("n_qc_retained", 0))):
            value = str(record.get(f"{prefix}_artifact", ""))
            path = path_at(root, value) if value else Path()
            if not value or not path.is_file():
                errors.append(f"missing {prefix}_artifact")
                continue
            expected_hash = str(record.get(f"{prefix}_sha256", ""))
            if not expected_hash or sha256(path) != expected_hash:
                errors.append(f"{prefix}_sha256 mismatch")
            try:
                partitions[prefix] = membership_set(path)
                observed_n = len(partitions[prefix])
                registry_n = int(float(record.get(f"{prefix}_n_observations", expected_n) or 0))
                if observed_n != registry_n or observed_n != int(float(expected_n or 0)):
                    errors.append(f"{prefix} count mismatch")
            except (ValueError, TypeError):
                errors.append(f"invalid {prefix} membership/count")
        try:
            query = membership_set(path_at(root, str(record.get("query_membership_artifact", ""))))
            accepted = partitions.get("accepted_membership", set())
            rejected = partitions.get("rejected_membership", set())
            if accepted & rejected or accepted | rejected != query:
                errors.append("Atlas accepted/rejected memberships do not exactly partition the query")
        except (OSError, ValueError):
            errors.append("unable to validate Atlas membership partition")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--record", required=True, type=Path)
    parser.add_argument("--biological-profile", type=Path, help="retained for provenance/profile binding")
    args = parser.parse_args()
    root = args.project_root.resolve()
    registry = root / "state/route_attempt_registry.tsv"
    record = json.loads(args.record.read_text(encoding="utf-8"))
    required = ["route_attempt_id", "sample_id", "route_class", "applicability", "status"]
    missing = [field for field in required if not str(record.get(field, "")).strip()]
    if missing:
        raise SystemExit(f"missing route fields: {missing}")
    if record["route_class"] not in ROUTE_CLASSES:
        raise SystemExit(f"invalid route_class: {record['route_class']}")
    now = datetime.now(timezone.utc).isoformat()
    record.setdefault("created_at", now)
    if record.get("status") in TERMINAL:
        record.setdefault("closed_at", now)
        errors = validate_terminal(root, record)
        if errors:
            raise SystemExit("terminal assisted route is invalid: " + "; ".join(errors))
    lock_path = registry.with_suffix(registry.suffix + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("w") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        if not registry.exists():
            raise SystemExit("route_attempt_registry.tsv is missing; initialize or migrate the project")
        with registry.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames or []
            rows = list(reader)
        old = [row for row in rows if row.get("route_attempt_id") == record["route_attempt_id"]]
        if old and old[0].get("status") in TERMINAL:
            raise SystemExit("terminal route_attempt_id is immutable; create a superseding ID")
        rows = [row for row in rows if row.get("route_attempt_id") != record["route_attempt_id"]]
        rows.append({field: record.get(field, "") for field in fields})
        fd, tmp_name = tempfile.mkstemp(prefix=registry.name + ".", suffix=".tmp", dir=registry.parent)
        try:
            with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
                handle.flush()
                os.fsync(handle.fileno())
            Path(tmp_name).replace(registry)
        finally:
            Path(tmp_name).unlink(missing_ok=True)
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
    print(record["route_attempt_id"], record["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
