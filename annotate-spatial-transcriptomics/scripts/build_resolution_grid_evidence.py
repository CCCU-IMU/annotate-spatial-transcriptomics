#!/usr/bin/env python3
"""Build deterministic biological resolution evidence from a full-grid score run."""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

from lineage_controller_lib import (
    canonical_cluster_challenger,
    candidate_can_release,
    catalog_candidates,
    dominant_generic_remainder_group,
    effective_broad_writeback_strategy,
    group_candidate_detected,
    group_candidate_score,
    group_identity_core_fraction,
    group_identity_core_direct_fraction,
    group_orthogonal_support_count,
    group_pseudobulk_contrast,
    group_release_supported_fraction,
    local_split_worthy_group_program,
    number,
    read_tsv,
    sha256,
    write_manifest,
    write_tsv,
)


def clamp(value: float) -> float:
    return max(0.0, min(1.0, value))


def candidate_detected(
    row: dict[str, str], candidate: dict | None = None
) -> bool:
    """Compatibility wrapper around the identity-core group detector."""
    return group_candidate_detected(row, candidate)


def candidate_rank(
    row: dict[str, str], candidate: dict | None = None
) -> tuple[float, ...]:
    return (
        float(candidate_detected(row, candidate)),
        group_candidate_score(row, candidate),
        group_identity_core_fraction(row),
        number(row.get("marker_deg_log2fc_mean")),
    )


def pseudobulk_score(row: dict[str, str]) -> float:
    n = max(1.0, number(row.get("n_observations"), 1.0))
    positive = max(0.0, number(row.get("positive_marker_pseudobulk_sum"))) / n
    anti = max(0.0, number(row.get("anti_marker_pseudobulk_sum"))) / n
    if positive <= 0:
        return 0.0
    contrast = positive / (positive + anti + 1e-12)
    detection = clamp(number(row.get("positive_marker_detection_fraction")))
    return clamp(math.sqrt(contrast * detection))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scoring-output", required=True, type=Path)
    parser.add_argument("--catalog", required=True, type=Path)
    parser.add_argument(
        "--selection-purpose", required=True,
        choices=["whole_tissue_cohort_partition", "cohort_identity_resolution"],
    )
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    scoring = args.scoring_output.resolve()
    manifest_path = scoring / "observation_scoring_manifest.json"
    evidence_path = scoring / "tables/cluster_candidate_multichannel_evidence.tsv.gz"
    if not manifest_path.is_file() or not evidence_path.is_file():
        raise SystemExit("full-grid scoring output is incomplete")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS":
        raise SystemExit("full-grid observation scoring did not pass")
    scorer = manifest.get("scorer", {})
    scorer_path = Path(str(scorer.get("path", "")))
    if (
        not scorer_path.is_file()
        or scorer.get("sha256") != sha256(scorer_path)
    ):
        raise SystemExit("full-grid scorer path/hash is absent or stale")
    catalog_doc = json.loads(args.catalog.read_text(encoding="utf-8"))
    candidate_map = catalog_candidates(catalog_doc)
    candidate_ids = set(candidate_map)
    if set(manifest.get("candidate_universe", [])) != candidate_ids:
        raise SystemExit("full-grid scoring candidate universe differs from catalog")

    rows = read_tsv(evidence_path)
    if not rows or "resolution" not in rows[0]:
        raise SystemExit("cluster multichannel evidence lacks full-grid resolution")
    by_resolution: dict[float, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_resolution[number(row.get("resolution"))].append(row)
    if len(by_resolution) < 3:
        raise SystemExit("resolution evidence requires at least three partitions")

    complete: dict[float, bool] = {}
    best: dict[float, dict[str, dict[str, str]]] = {}
    detected: dict[float, set[str]] = {}
    for resolution, resolution_rows in by_resolution.items():
        by_cluster: dict[tuple[str, str], set[str]] = defaultdict(set)
        best[resolution] = {}
        for row in resolution_rows:
            group = (row.get("source_boundary", ""), row.get("source_cluster", ""))
            candidate = row.get("candidate_id", "")
            by_cluster[group].add(candidate)
            previous = best[resolution].get(candidate)
            candidate_meta = candidate_map.get(str(candidate), {})
            if previous is None or candidate_rank(
                row, candidate_meta
            ) > candidate_rank(previous, candidate_meta):
                best[resolution][candidate] = row
        complete[resolution] = bool(by_cluster) and all(
            candidates == candidate_ids for candidates in by_cluster.values()
        )
        detected[resolution] = {
            candidate
            for candidate, row in best[resolution].items()
            if candidate_detected(row, candidate_map.get(candidate, {}))
        }
    positive_universe = set().union(*detected.values())

    output_rows: list[dict[str, object]] = []
    for resolution in sorted(by_resolution):
        resolution_rows = by_resolution[resolution]
        cluster_sizes = {
            (row.get("source_boundary", ""), row.get("source_cluster", "")):
            int(number(row.get("n_observations")))
            for row in resolution_rows
        }
        total_n = sum(cluster_sizes.values())
        small_cutoff = max(100, int(math.ceil(0.001 * total_n)))
        technical_fragmentation = (
            sum(size for size in cluster_sizes.values() if size < small_cutoff) / total_n
            if total_n else 1.0
        )
        candidates = [
            best[resolution][candidate]
            for candidate in sorted(detected[resolution])
        ]
        divisor = len(candidates) or 1
        rows_by_group: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for row in resolution_rows:
            rows_by_group[(
                str(row.get("source_boundary", "")),
                str(row.get("source_cluster", "")),
            )].append(row)
        state_driven_n = 0
        for group, group_rows in rows_by_group.items():
            identity_scores = []
            nonidentity_scores = []
            for row in group_rows:
                candidate = candidate_map.get(str(row.get("candidate_id", "")), {})
                if not group_candidate_detected(row, candidate):
                    continue
                score = group_candidate_score(row, candidate)
                if candidate_can_release(candidate):
                    identity_scores.append(score)
                elif str(candidate.get("candidate_role", "")) in {
                    "state", "exploratory"
                }:
                    nonidentity_scores.append(score)
            identity_best = max(identity_scores, default=0.0)
            nonidentity_best = max(nonidentity_scores, default=0.0)
            if nonidentity_best >= identity_best + 0.10:
                state_driven_n += cluster_sizes[group]
        state_overfragmentation = (
            state_driven_n / total_n if total_n else 1.0
        )
        directly_resolved_n = 0
        mixed_n = 0
        unresolved_n = 0
        weighted_margin = 0.0
        weighted_winner_core = 0.0
        for group, group_rows in rows_by_group.items():
            broad_by_label: dict[str, tuple[float, dict[str, str], dict]] = {}
            for row in group_rows:
                candidate_id = str(row.get("candidate_id", ""))
                candidate = candidate_map.get(candidate_id, {})
                if not candidate_can_release(candidate):
                    continue
                label = str(candidate.get("release_broad_label", ""))
                if not label:
                    continue
                value = (group_candidate_score(row, candidate), row, candidate)
                if label not in broad_by_label or value[0] > broad_by_label[label][0]:
                    broad_by_label[label] = value
            positive = [
                value for value in broad_by_label.values()
                if group_candidate_detected(value[1], value[2])
            ]
            specific_positive = [
                value for value in positive
                if str(value[2].get("candidate_id", "")) != "stromal_mesenchymal"
            ]
            specific_split_worthy = [
                value for value in specific_positive
                if local_split_worthy_group_program(value[1], value[2])
            ]
            generic_positive = [
                value for value in positive
                if str(value[2].get("candidate_id", "")) == "stromal_mesenchymal"
            ]
            ranked = sorted(
                (
                    specific_split_worthy
                    or generic_positive
                    or specific_positive
                    or positive
                ),
                key=lambda value: (
                    -value[0],
                    str(value[1].get("candidate_id", "")),
                ),
            )
            winner = ranked[0] if ranked else None
            competitors = [value for value in positive if value is not winner]
            group_n = cluster_sizes[group]
            if winner is None:
                unresolved_n += group_n
                continue
            winner_core = group_release_supported_fraction(winner[1], winner[2])
            runner_core = max(
                (
                    group_release_supported_fraction(value[1], value[2])
                    for value in competitors
                ),
                default=0.0,
            )
            winner_margin = max(0.0, winner_core - runner_core)
            split_worthy_competitors = [
                value for value in competitors
                if local_split_worthy_group_program(value[1], value[2])
            ]
            whole_strategy_ok = bool(
                effective_broad_writeback_strategy(winner[2])
                != "candidate_local_component_never_parent_expansion"
            )
            strict_whole_pass = (
                whole_strategy_ok
                and winner_core >= 0.40
                and winner_margin >= 0.20
                and number(winner[1].get("hard_contradiction_fraction")) <= 0.05
                and group_orthogonal_support_count(winner[1], winner[2]) >= 3
                and not split_worthy_competitors
            )
            dominant_whole_pass = (
                whole_strategy_ok
                and str(winner[2].get("candidate_id", ""))
                != "stromal_mesenchymal"
                and not canonical_cluster_challenger(winner[1], winner[2])
                and number(winner[1].get("observation_seed_fraction")) >= 0.70
                and group_identity_core_direct_fraction(winner[1]) >= 0.80
                and group_identity_core_fraction(winner[1]) >= 0.80
                and number(winner[1].get("hard_contradiction_fraction")) <= 0.20
                and group_orthogonal_support_count(winner[1], winner[2]) >= 3
                and not split_worthy_competitors
            )
            generic_whole_pass = (
                whole_strategy_ok
                and dominant_generic_remainder_group(winner[1], winner[2])
                and not split_worthy_competitors
                and not specific_split_worthy
            )
            whole_pass = (
                strict_whole_pass
                or dominant_whole_pass
                or generic_whole_pass
            )
            if whole_pass:
                directly_resolved_n += group_n
            elif specific_split_worthy:
                mixed_n += group_n
            else:
                unresolved_n += group_n
            weighted_margin += group_n * winner_margin
            weighted_winner_core += group_n * winner_core
        catalog_recall = (
            len(detected[resolution]) / len(positive_universe)
            if positive_universe else 1.0
        )
        embedded_program_separation = (
            weighted_margin / total_n if total_n else 0.0
        )
        deg_antideg_coherence = sum(
            clamp((
                number(row.get("marker_deg_log2fc_mean"))
                - 0.25 * max(
                    0.0, number(row.get("anti_marker_deg_log2fc_mean"))
                )
                + 0.25
            ) / 1.5)
            for row in candidates
        ) / divisor
        pseudobulk_coherence = sum(
            group_pseudobulk_contrast(row)
            for row in candidates
        ) / divisor
        spatial_morphology_coherence = sum(
            clamp(number(row.get("spatial_group_connectivity_fraction")))
            for row in candidates
        ) / divisor
        adjacent_resolution_stability = sum(
            clamp(number(row.get("cross_resolution_stable_fraction")))
            for row in candidates
        ) / divisor
        observation_support_coherence = sum(
            group_identity_core_fraction(row)
            for row in candidates
        ) / divisor
        output_rows.append({
            "resolution": resolution,
            "selection_purpose": args.selection_purpose,
            "complete_catalog_scanned": str(complete[resolution]).lower(),
            "zero_census_audited": str(complete[resolution]).lower(),
            "catalog_recall": catalog_recall,
            "embedded_program_separation": embedded_program_separation,
            "deg_antideg_coherence": deg_antideg_coherence,
            "pseudobulk_coherence": pseudobulk_coherence,
            "spatial_morphology_coherence": spatial_morphology_coherence,
            "adjacent_resolution_stability": adjacent_resolution_stability,
            "observation_support_coherence": observation_support_coherence,
            "directly_resolved_observation_fraction": (
                directly_resolved_n / total_n if total_n else 0.0
            ),
            "mixed_observation_fraction": mixed_n / total_n if total_n else 1.0,
            "unresolved_observation_fraction": (
                unresolved_n / total_n if total_n else 1.0
            ),
            "mean_identity_margin": weighted_margin / total_n if total_n else 0.0,
            "mean_winner_identity_core_fraction": (
                weighted_winner_core / total_n if total_n else 0.0
            ),
            "technical_fragmentation": technical_fragmentation,
            "state_overfragmentation": state_overfragmentation,
            "complexity": len(cluster_sizes),
            "positive_candidate_n": len(detected[resolution]),
            "positive_candidate_universe_n": len(positive_universe),
        })

    args.out.mkdir(parents=True, exist_ok=True)
    output_path = args.out / "resolution_grid_evidence.tsv"
    write_tsv(output_path, output_rows)
    result = {
        "status": "PASS",
        "schema_version": "2.2",
        "builder_version": "2.2.0",
        "selection_purpose": args.selection_purpose,
        "builder": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "scoring_manifest": {
            "path": str(manifest_path),
            "sha256": sha256(manifest_path),
        },
        "scorer": scorer,
        "cluster_multichannel_evidence": {
            "path": str(evidence_path),
            "sha256": sha256(evidence_path),
        },
        "candidate_catalog": {
            "path": str(args.catalog.resolve()),
            "sha256": sha256(args.catalog),
        },
        "resolution_grid_evidence": {
            "path": str(output_path.resolve()),
            "sha256": sha256(output_path),
        },
        "resolutions": sorted(by_resolution),
        "complete_catalog_all_resolutions": all(complete.values()),
    }
    write_manifest(args.out / "resolution_grid_evidence_manifest.json", result)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
