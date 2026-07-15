#!/usr/bin/env python3
"""Fail-closed stage detector for an autonomous annotation project.

This script does not make biological decisions.  It prevents an Agent from
mistaking a successful intermediate analysis for a completed release and emits
the next deterministic control action.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path

from master_quality_lib import validate_bound_record, validate_master_approval


TERMINAL_RUN_STATES = {"validated_done", "failed_preserved", "cancelled_preserved"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def json_status(path: Path, key: str = "status") -> str:
    if not path.exists():
        return "MISSING"
    try:
        return str(json.loads(path.read_text(encoding="utf-8")).get(key, "MISSING"))
    except (OSError, json.JSONDecodeError):
        return "INVALID"


def newest_mtime(paths: list[Path]) -> float:
    return max((path.stat().st_mtime for path in paths if path.exists()), default=0.0)


def stale(target: Path, dependencies: list[Path]) -> bool:
    return not target.exists() or target.stat().st_mtime < newest_mtime(dependencies)


def newest_match(root: Path, patterns: list[str]) -> Path | None:
    matches = [path for pattern in patterns for path in root.glob(pattern) if path.is_file()]
    return max(matches, key=lambda path: path.stat().st_mtime) if matches else None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def confirmation_valid(path: Path, root: Path) -> bool:
    if not path.is_file():
        return False
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    if record.get("status") != "CONFIRMED":
        return False
    for key, hash_key in (
        ("cell_ledger", "cell_ledger_sha256"),
        ("cluster_ledger", "cluster_ledger_sha256"),
        ("completion_gate", "completion_gate_sha256"),
        ("release_taxonomy_audit", "release_taxonomy_audit_sha256"),
        ("master_quality_approval", "master_quality_approval_sha256"),
        ("confirmation_review_report", "confirmation_review_report_sha256"),
        ("confirmation_review_manifest", "confirmation_review_manifest_sha256"),
        ("annotation_support_registry", "annotation_support_registry_sha256"),
    ):
        target = root / str(record.get(key, ""))
        if not target.is_file() or sha256(target) != record.get(hash_key):
            return False
    return True


def read_confirmation(path: Path) -> dict:
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def parse_time(value: str) -> datetime | None:
    try:
        return datetime.fromisoformat(value.replace("Z", "+00:00"))
    except (TypeError, ValueError):
        return None


def release_only_run(row: dict[str, str]) -> bool:
    """Return whether a post-confirmation run can only build release evidence.

    These runs may update the run registry, but cannot alter biological ledgers.
    They therefore must be monitored and audited without invalidating the
    user's hash-bound annotation confirmation.
    """
    stage = row.get("stage", "").strip().lower()
    return stage.startswith(("final_", "release_", "report_")) or stage in {
        "marker_asset_generation",
        "annotation_map_generation",
        "spatial_gene_map_generation",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()

    actions: list[dict[str, str]] = []

    def add(phase: str, command: str, reason: str) -> None:
        actions.append({"priority": str(len(actions) + 1), "phase": phase, "command": command, "reason": reason})

    project = root / "config/project.json"
    context = root / "config/biological_context.json"
    context_validation = root / "provenance/biological_context_validation.json"
    discovery = root / "input_discovery"
    cluster_ledger = root / "state/cluster_decision_ledger.tsv"
    cell_ledger = root / "state/cell_ledger.tsv.gz"
    run_registry = root / "state/run_registry.tsv"
    route_registry = root / "state/route_attempt_registry.tsv"
    pool_registry = root / "state/pool_registry.tsv"
    view_registry = root / "state/annotation_view_registry.tsv"
    state_validation = root / "provenance/state_validation.json"
    iteration_plan = root / "provenance/iteration_plan.json"
    queue = root / "state/next_action_queue.tsv"
    completion = root / "provenance/completion_gate.json"
    taxonomy_audit = root / "provenance/release_taxonomy_audit.json"
    confirmation_request = root / "provenance/final_annotation_confirmation_request.json"
    review_asset_manifest = root / "review/confirmation/assets/review_asset_manifest.json"
    master_review_request = root / "provenance/master_quality_review_request.json"
    master_approval = root / "state/master_quality_approval.json"
    confirmation_review = root / "review/confirmation/index.html"
    confirmation_review_manifest = root / "provenance/confirmation_review_manifest.json"
    confirmation = root / "state/final_annotation_confirmation.json"
    report = root / "report/index.html"
    release_manifest = root / "provenance/release_manifest.tsv"
    checksums = root / "provenance/checksums.sha256"
    release_audit = root / "provenance/release_audit.json"
    final_census = root / "tables/final_annotation_census.tsv"
    confirmation_record = read_confirmation(confirmation)
    user_confirmed = confirmation_valid(confirmation, root)

    if not project.exists():
        add("initialize", "discover_inputs.py, then init_annotation_project.py", "project configuration is missing")
    if not discovery.exists() or not any(discovery.rglob("*")):
        add("discover", "discover_inputs.py", "frozen input discovery inventory is missing")
    if not context.exists():
        add("context", "create config/biological_context.json", "species, tissue, stage, platform and priority lineages are required")
    elif json_status(context_validation) != "PASS":
        add("context", "validate_biological_context.py", "biological context has not passed validation")

    if not cluster_ledger.exists() or not cell_ledger.exists():
        add("broad_annotation", "select clustering and write the initial broad/pool ledgers", "cell- and cluster-level annotation state is incomplete")

    runs = read_tsv(run_registry)
    unfinished = [row for row in runs if row.get("status", "") not in TERMINAL_RUN_STATES]
    if unfinished:
        add("run_control", "monitor, inspect logs, validate artifacts, then update run_registry.tsv", f"{len(unfinished)} submitted/running/nonterminal run(s) remain")

    confirmed_at = parse_time(str(confirmation_record.get("confirmed_at", ""))) if user_confirmed else None
    post_confirmation_runs = []
    if confirmed_at is not None:
        for row in runs:
            started = parse_time(row.get("started_at", ""))
            if started is not None and started > confirmed_at:
                post_confirmation_runs.append(row)
    unexpected_post_confirmation = [row for row in post_confirmation_runs if not release_only_run(row)]
    if unexpected_post_confirmation:
        add(
            "confirmation_invalidation",
            "review the post-confirmation biological runs, rerun biological gates, and request a new user confirmation",
            f"{len(unexpected_post_confirmation)} non-release run(s) started after the confirmed annotation snapshot",
        )

    final_views = [row for row in read_tsv(view_registry) if row.get("view") == "final" and row.get("status") == "validated"]
    final_membership = None
    if final_views:
        value = final_views[-1].get("membership_path", "")
        final_membership = Path(value) if value else None
        if final_membership is not None and not final_membership.is_absolute():
            final_membership = root / final_membership
    if cell_ledger.exists() and (
        not final_views or final_membership is None or not final_membership.exists()
        or stale(final_census, [cell_ledger])
    ):
        add(
            "final_annotation",
            "build_final_annotation.py PROJECT_ROOT --cell-ledger state/cell_ledger.tsv.gz --out state/cell_ledger.tsv.gz --sample SAMPLE --version DECISION_VERSION",
            "single final annotation is missing or stale; materialize moderate-or-higher broad and high-confidence fine labels before state/taxonomy/completion audits",
        )

    # Once the cell ledger, cluster ledger and completion gate have been
    # hash-bound by explicit confirmation, release-only jobs may append to the
    # run registry without making the biological snapshot stale. Unfinished
    # jobs are still caught by run_control, and the final manifest/audit still
    # depend on the current run registry.
    state_dependencies = [cluster_ledger, cell_ledger, pool_registry, route_registry, view_registry, final_census]
    if not user_confirmed or unexpected_post_confirmation:
        state_dependencies.append(run_registry)
    if cluster_ledger.exists() and cell_ledger.exists() and (stale(state_validation, state_dependencies) or json_status(state_validation) != "PASS"):
        add("state_validation", "validate_state.py", "state validation is missing, failed or stale after a writeback")

    planning_dependencies = [cluster_ledger, pool_registry, route_registry, view_registry, state_validation]
    if cluster_ledger.exists() and (stale(iteration_plan, planning_dependencies) or json_status(iteration_plan) not in {"READY_FOR_COMPLETION_AUDIT", "ITERATION_REQUIRED"}):
        add("iteration_planning", "plan_next_iteration.py", "the next-action plan is missing, invalid or stale")
    queued = read_tsv(queue)
    if queued:
        add("iterative_annotation", "execute every row in state/next_action_queue.tsv and re-plan after each atomic writeback", f"{len(queued)} route action(s) remain")

    taxonomy_record = read_confirmation(taxonomy_audit)
    if cell_ledger.exists() and final_views and (
        stale(taxonomy_audit, [cell_ledger]) or taxonomy_record.get("pass") is not True
    ):
        add(
            "release_taxonomy_audit",
            "audit_release_taxonomy.py state/cell_ledger.tsv.gz --profile ACTIVE_BIOLOGICAL_PROFILE --broad-column final_broad_label --status-column final_state --pool-column target_pool --out provenance/release_taxonomy_audit.json",
            "biological broad classes and retained anatomical/QC/technical states have not passed the current taxonomy/pool audit",
        )

    completion_dependencies = planning_dependencies + [iteration_plan, queue, taxonomy_audit]
    if not user_confirmed or unexpected_post_confirmation:
        completion_dependencies.append(run_registry)
    if cluster_ledger.exists() and not queued and (stale(completion, completion_dependencies) or json_status(completion) != "PASS"):
        add("completion_gate", "check_completion_gate.py", "biological completion gate is missing, blocked or stale")

    if json_status(completion) == "PASS" and not user_confirmed:
        review_dependencies = [completion, cluster_ledger, cell_ledger, taxonomy_audit, root / "state/annotation_support_registry.tsv"]
        master_request_record = read_confirmation(master_review_request)
        master_request_ok = (
            master_request_record.get("status") == "AWAITING_MASTER_QUALITY_APPROVAL"
            and not validate_bound_record(root, master_request_record)
        )
        master_ok, master_errors, master_record = validate_master_approval(root, master_approval)
        if stale(review_asset_manifest, review_dependencies) or json_status(review_asset_manifest) != "PASS":
            add(
                "confirmation_review_assets",
                "build_confirmation_review_assets.R",
                "build the frozen broad spatial and canonical-marker evidence only after all annotation and completion gates pass",
            )
        elif not master_request_ok:
            add(
                "master_quality_review_request",
                "request_master_quality_review.py",
                "freeze the fully annotated sample for main-Agent biological quality review",
            )
        elif not master_ok:
            action = "return the sample worker to targeted annotation iteration" if master_record.get("status") == "RETURN_FOR_ITERATION" else "main conversation Agent reviews the frozen annotation and runs record_master_quality_approval.py"
            add(
                "master_quality_approval",
                action,
                "quality approval is required after all annotation routes/final writeback, never after broad annotation alone" + (": " + "; ".join(master_errors) if master_errors else ""),
            )
        elif stale(confirmation_review, review_dependencies + [review_asset_manifest, master_approval]) or json_status(confirmation_review_manifest) != "PASS":
            add(
                "confirmation_review",
                "build_confirmation_review.py",
                "build the user-facing lightweight HTML only after the main-Agent annotation-quality approval",
            )
        elif stale(confirmation_request, review_dependencies + [master_approval, confirmation_review, confirmation_review_manifest]):
            add(
                "user_confirmation",
                "request_final_annotation_confirmation.py",
                "freeze and present the lightweight review plus final census before spending compute on full release assets",
            )
        else:
            add(
                "user_confirmation",
                "present provenance/final_annotation_confirmation_request.json to the user; after explicit approval run record_final_annotation_confirmation.py",
                "final annotation has not received explicit user approval",
            )

    final_metadata = newest_match(root, ["tables/final_cell_metadata.tsv", "tables/final_cell_metadata.tsv.gz", "tables/final_cell_metadata_v*.tsv", "tables/final_cell_metadata_v*.tsv.gz"])
    has_final_fine = False
    if final_census.exists():
        try:
            has_final_fine = any(row.get("fine_label", "") and int(float(row.get("n_observations", 0) or 0)) > 0 for row in read_tsv(final_census))
        except (TypeError, ValueError):
            has_final_fine = True
    required_assets = [
        final_metadata or root / "tables/final_cell_metadata.tsv.gz",
        root / "tables/broad_DEG_one_vs_rest_all.tsv",
        root / "figures/marker_dotplots/marker_dotplot_asset_index.tsv",
        root / "figures/final_broad_UMAP.png",
        root / "provenance/release_sessionInfo.txt",
    ]
    if has_final_fine:
        required_assets.extend([root / "tables/subtype_DEG_one_vs_rest_all.tsv", root / "figures/final_subtype_UMAP.png"])
    modality = ""
    if project.exists():
        try:
            modality = str(json.loads(project.read_text(encoding="utf-8")).get("modality", ""))
        except json.JSONDecodeError:
            pass
    if modality == "spatial":
        required_assets.extend([root / "figures/final_broad_spatial.png", root / "tables/spatial_node_asset_index.tsv"])
        if has_final_fine:
            required_assets.append(root / "figures/final_subtype_spatial.png")
    missing_assets = [str(path.relative_to(root)) for path in required_assets if not path.exists()]
    if json_status(completion) == "PASS" and user_confirmed and missing_assets:
        add("final_assets", "prepare final metadata, all-accepted broad DEG/dotplots and high-confidence subtype assets", f"{len(missing_assets)} required final asset(s) are missing")

    asset_dependencies = [path for path in required_assets if path.exists()] + [completion]
    if json_status(completion) == "PASS" and user_confirmed and not missing_assets and (not report.exists() or stale(report, asset_dependencies + [confirmation])):
        add("report", "build_report.py", "final report is missing or older than its evidence assets")

    release_dependencies = [report, completion, cluster_ledger, cell_ledger, run_registry] + [path for path in required_assets if path.exists()]
    if user_confirmed and report.exists() and (stale(release_manifest, release_dependencies + [confirmation]) or stale(checksums, release_dependencies + [confirmation])):
        add("release_manifest", "build_release_manifest.py", "release manifest/checksums are missing or stale")
    audit_dependencies = [report, completion, release_manifest, checksums]
    if user_confirmed and release_manifest.exists() and checksums.exists() and (stale(release_audit, audit_dependencies) or json_status(release_audit) != "PASS"):
        add("release_audit", "audit_release.py --profile full", "release audit is missing, failed or stale")

    if actions:
        status = "CONTINUE"
        phase = actions[0]["phase"]
        terminal = False
    else:
        status = "COMPLETE"
        phase = "complete"
        terminal = True

    result = {
        "status": status,
        "phase": phase,
        "terminal": terminal,
        "next_actions": actions,
        "rule": "Do not stop or call the project final while status is CONTINUE.",
    }
    output = root / "provenance/autopilot_status.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    queue_out = root / "state/autopilot_next_actions.tsv"
    queue_out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["priority", "phase", "command", "reason"]
    with queue_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(actions)

    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if terminal else 2


if __name__ == "__main__":
    raise SystemExit(main())
