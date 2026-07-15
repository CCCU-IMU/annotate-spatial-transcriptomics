#!/usr/bin/env python3
"""Legacy v1.1-v1.4 three-view migration helper; never use for a v1.5 release."""

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


QC_RESCUE = {"qc_anchor_recluster_broad_rescue", "depth_matched_atlas_anchor_mapping_fullfeature_rescue"}
INTERFACE_RETURN = {"interface_rctd_broad_return", "interface_calibrated_broad_return"}
VIEW_FIELDS = [
    "analysis_scope", "source_key", "parent_decision_id", "pool_snapshot_id", "generation", "route_attempt_id",
    "state_tags", "spatial_tags", "qc_tags", "candidate_lineages",
    "strict_state", "strict_broad_label", "strict_fine_label",
    "inclusive_state", "inclusive_broad_label", "inclusive_fine_label",
    "display_state", "display_broad_label", "display_fine_label", "display_policy",
]


def opener(path: Path, mode: str):
    return gzip.open(path, mode, newline="", encoding="utf-8") if path.suffix == ".gz" else path.open(mode.replace("t", ""), newline="", encoding="utf-8")


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--cell-ledger", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--view-version", default="v001")
    parser.add_argument(
        "--legacy-migration-only", action="store_true",
        help="required acknowledgement that these compatibility fields will not be published",
    )
    parser.add_argument(
        "--mode", choices=["derive", "preserve"], default="derive",
        help="derive views from route policy, or preserve already adjudicated view fields and only refresh census/registry",
    )
    parser.add_argument(
        "--registry-membership-path", type=Path,
        help="authoritative membership file to hash/register; defaults to --out",
    )
    args = parser.parse_args()
    if not args.legacy_migration_only:
        raise SystemExit(
            "strict/inclusive/display generation is retired. Use build_final_annotation.py; "
            "pass --legacy-migration-only only to inspect an old project before migration."
        )
    route_registry = args.project_root / "state/route_attempt_registry.tsv"
    qc_routes, interface_routes = set(QC_RESCUE), set(INTERFACE_RETURN)
    if route_registry.exists():
        with route_registry.open(newline="", encoding="utf-8") as handle:
            for route in csv.DictReader(handle, delimiter="\t"):
                if route.get("status") != "validated": continue
                if route.get("route_class") == "qc_atlas_review": qc_routes.add(route.get("route_attempt_id", ""))
                if route.get("route_class") == "interface_deconvolution_review": interface_routes.add(route.get("route_attempt_id", ""))
    with opener(args.cell_ledger, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        input_fields = reader.fieldnames or []
        if args.mode == "preserve":
            required_existing = {
                "analysis_scope",
                "strict_state", "strict_broad_label", "strict_fine_label",
                "inclusive_state", "inclusive_broad_label", "inclusive_fine_label",
                "display_state", "display_broad_label", "display_fine_label", "display_policy",
            }
            missing_existing = sorted(required_existing - set(input_fields))
            if missing_existing:
                raise SystemExit(
                    "--mode preserve requires existing adjudicated view fields: "
                    + ",".join(missing_existing)
                )
        fields = input_fields + [x for x in VIEW_FIELDS if x not in input_fields]
        args.out.parent.mkdir(parents=True, exist_ok=True)
        tmp = args.out.with_name(args.out.stem + ".tmp" + args.out.suffix)
        counts = {view: {} for view in ["strict", "inclusive", "display"]}
        n_analysis = 0
        with opener(tmp, "wt") as output:
            writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for row in reader:
                state, broad, fine, route = row.get("state", ""), row.get("broad_label", ""), row.get("fine_label", ""), row.get("route", "")
                scope = row.get("analysis_scope") or ("excluded_initial_qc" if state == "excluded_initial_qc" else "analysis_set")
                row["analysis_scope"] = scope
                row["source_key"] = row.get("source_key") or f"{row.get('source_run_id','')}|{row.get('source_cluster','')}"
                row["generation"] = row.get("generation") or row.get("iteration", "1")
                for field in ["parent_decision_id", "pool_snapshot_id", "route_attempt_id", "state_tags", "spatial_tags", "qc_tags", "candidate_lineages"]:
                    row.setdefault(field, "")
                if args.mode == "preserve":
                    strict = (
                        row.get("strict_state", ""), row.get("strict_broad_label", ""),
                        row.get("strict_fine_label", ""),
                    )
                    inclusive = (
                        row.get("inclusive_state", ""), row.get("inclusive_broad_label", ""),
                        row.get("inclusive_fine_label", ""),
                    )
                    display = (
                        row.get("display_state", ""), row.get("display_broad_label", ""),
                        row.get("display_fine_label", ""),
                    )
                    policy = row.get("display_policy", "")
                    if not all(value[0] for value in [strict, inclusive, display]):
                        raise SystemExit(
                            f"empty preserved view state for cell_id={row.get('cell_id','')}"
                        )
                elif scope == "excluded_initial_qc":
                    strict = inclusive = display = ("excluded_initial_qc", "", "")
                    policy = "excluded_from_analysis_set_retained_in_full_ledger"
                elif route in qc_routes:
                    strict = ("qc_holdout", "", "")
                    inclusive = display = (state, broad, "")
                    policy = "calibrated_qc_broad_rescue_inclusive_display_only"
                    row["fine_anchor_eligible"] = "false"
                elif route in interface_routes and state == "defined_broad_only":
                    strict = ("interface_review", "", "")
                    inclusive = display = (state, broad, "")
                    policy = "calibrated_interface_broad_return_inclusive_display_only"
                    row["fine_anchor_eligible"] = "false"
                else:
                    strict = inclusive = display = (state, broad, fine)
                    policy = "same_all_views_reviewed_direct_or_retained_state"
                for view, value in [("strict", strict), ("inclusive", inclusive), ("display", display)]:
                    row[f"{view}_state"], row[f"{view}_broad_label"], row[f"{view}_fine_label"] = value
                    if scope == "analysis_set":
                        key = value
                        counts[view][key] = counts[view].get(key, 0) + 1
                row["display_policy"] = policy
                if scope == "analysis_set": n_analysis += 1
                writer.writerow({field: row.get(field, "") for field in fields})
    tmp.replace(args.out)
    tables = args.project_root / "tables"; tables.mkdir(parents=True, exist_ok=True)
    for view, values in counts.items():
        path = tables / f"{view}_annotation_census.tsv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t"); writer.writerow(["state", "broad_label", "fine_label", "n_observations"])
            for (state, broad, fine), n in sorted(values.items(), key=lambda x: -x[1]): writer.writerow([state, broad, fine, n])
    policy_path = args.project_root / "provenance/annotation_view_policy.json"; policy_path.parent.mkdir(parents=True, exist_ok=True)
    policy = {
        "status":"PASS", "mode": args.mode, "analysis_set_n":n_analysis,
        "qc_rescue_routes":sorted(qc_routes), "interface_return_routes":sorted(interface_routes),
        "rule": (
            "preserve pre-adjudicated strict/inclusive/display fields; refresh census and registry only"
            if args.mode == "preserve" else
            "strict excludes calibrated broad rescue; inclusive/display accept reviewed broad-only rescue; retained uncertainty remains explicit"
        ),
    }
    policy_path.write_text(json.dumps(policy, ensure_ascii=False, indent=2) + "\n")
    registry = args.project_root / "state/annotation_view_registry.tsv"
    lock_path = registry.with_suffix(registry.suffix + ".lock")
    with lock_path.open("w") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        with registry.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t"); registry_fields = reader.fieldnames or []; old = [x for x in reader if x.get("view") not in {"strict", "inclusive", "display"}]
        membership_path = (args.registry_membership_path or args.out).resolve()
        if not membership_path.exists():
            raise SystemExit(f"registry membership path does not exist: {membership_path}")
        now = datetime.now(timezone.utc).isoformat(); digest = sha256(membership_path)
        for view in ["strict", "inclusive", "display"]:
            old.append({"view_id":f"{args.sample}_{view}_{args.view_version}","sample_id":args.sample,"view":view,"membership_path":str(membership_path),"membership_sha256":digest,"n_observations":n_analysis,"policy":policy["rule"],"marker_deg_eligible":"true" if view in {"strict", "inclusive"} else "false","status":"validated","artifact":str((tables/f"{view}_annotation_census.tsv").resolve()),"created_at":now})
        fd, name = tempfile.mkstemp(prefix=registry.name + ".", suffix=".tmp", dir=registry.parent)
        try:
            with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=registry_fields, delimiter="\t"); writer.writeheader(); writer.writerows(old); handle.flush(); os.fsync(handle.fileno())
            Path(name).replace(registry)
        finally:
            Path(name).unlink(missing_ok=True)
        fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
    print(json.dumps(policy, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
