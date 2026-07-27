#!/usr/bin/env python3
"""Match stable nontechnical DEG/coexpression programs across resolutions."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from lineage_controller_lib import number, read_tsv, sha256, truth, write_manifest, write_tsv


PROGRAM_GENE_SETS = {
    "stress": {
        "FOS", "FOSB", "JUN", "JUNB", "JUND", "DUSP1", "ATF3",
        "HSPA1A", "HSPA1B", "HSP90AA1", "DNAJB1",
    },
    "hypoxia": {
        "HIF1A", "VEGFA", "BNIP3", "BNIP3L", "CA9", "EGLN3",
        "LDHA", "PGK1", "ENO1", "SLC2A1",
    },
    "cell_cycle": {
        "MKI67", "TOP2A", "UBE2C", "CENPF", "PCNA", "CCNB1",
        "CCNB2", "CDC20", "CDK1", "TYMS", "MCM2", "MCM3",
        "MCM4", "MCM5", "MCM6", "MCM7",
    },
    "pure_ecm": {
        "DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "COL5A1",
        "COL5A2", "COL6A1", "COL6A2", "COL6A3", "BGN", "SPARC",
        "FN1", "FBLN1", "FBLN2", "FBLN5", "ELN",
    },
}


def genes(row: dict[str, str]) -> set[str]:
    return {
        value.strip().upper()
        for value in str(row.get("genes", "")).replace(",", ";").split(";")
        if value.strip()
    }


def excluded_program_class(gene_set: set[str]) -> tuple[str, float]:
    if not gene_set:
        return "empty_program", 1.0
    prefix_counts = {
        "ribosomal": sum(
            gene.startswith(("RPL", "RPS", "MRPL", "MRPS"))
            for gene in gene_set
        ),
        "mitochondrial": sum(gene.startswith("MT-") for gene in gene_set),
    }
    fractions = {
        name: count / len(gene_set) for name, count in prefix_counts.items()
    }
    fractions.update({
        name: len(gene_set & members) / len(gene_set)
        for name, members in PROGRAM_GENE_SETS.items()
    })
    name, fraction = sorted(
        fractions.items(), key=lambda item: (-item[1], item[0])
    )[0]
    return (name, fraction) if fraction >= 0.40 else ("", fraction)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--programs", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--minimum-jaccard", type=float, default=0.25)
    args = parser.parse_args()
    eligible_rows = [
        row for row in read_tsv(args.programs)
        if row.get("candidate_status") == "unmodeled_program_seed"
        and truth(row.get("spatially_coherent"))
        and number(row.get("catalog_marker_overlap_fraction"), 1) < 0.30
        and int(number(row.get("coexpressed_gene_count"), 0)) >= 2
    ]
    excluded_rows: list[dict[str, object]] = []
    rows = []
    for row in eligible_rows:
        program_genes = genes(row)
        excluded_class, excluded_fraction = excluded_program_class(program_genes)
        if (
            not excluded_class
            and (
                int(number(row.get("best_catalog_overlap_family_count"), 0)) >= 2
                or int(number(row.get("best_catalog_overlap_gene_count"), 0)) >= 3
            )
        ):
            excluded_class = "modeled_catalog_program"
            excluded_fraction = number(
                row.get("best_catalog_overlap_gene_count"), 0
            ) / max(1, len(program_genes))
        if excluded_class:
            excluded_rows.append({
                "program_id": row.get("program_id", ""),
                "resolution": row.get("resolution", ""),
                "excluded_program_class": excluded_class,
                "excluded_gene_fraction": excluded_fraction,
                "best_catalog_candidate_id": row.get(
                    "best_catalog_candidate_id", ""
                ),
                "best_catalog_overlap_gene_count": row.get(
                    "best_catalog_overlap_gene_count", ""
                ),
                "best_catalog_overlap_family_count": row.get(
                    "best_catalog_overlap_family_count", ""
                ),
                "genes": ";".join(sorted(program_genes)),
            })
        else:
            rows.append(row)
    parent = list(range(len(rows)))

    def root(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        a, b = root(left), root(right)
        if a != b:
            parent[max(a, b)] = min(a, b)

    gene_sets = [genes(row) for row in rows]
    for left in range(len(rows)):
        for right in range(left + 1, len(rows)):
            if rows[left].get("resolution") == rows[right].get("resolution"):
                continue
            union_size = len(gene_sets[left] | gene_sets[right])
            jaccard = (
                len(gene_sets[left] & gene_sets[right]) / union_size if union_size else 0
            )
            if jaccard >= args.minimum_jaccard:
                union(left, right)

    groups: dict[int, list[int]] = defaultdict(list)
    for index in range(len(rows)):
        groups[root(index)].append(index)
    output = []
    for group_index, members in sorted(groups.items()):
        resolutions = sorted({str(rows[index]["resolution"]) for index in members})
        if len(resolutions) < 2:
            continue
        gene_count: dict[str, int] = defaultdict(int)
        for index in members:
            for gene in gene_sets[index]:
                gene_count[gene] += 1
        recurrent = sorted(
            (gene for gene, count in gene_count.items() if count >= 2),
            key=lambda gene: (-gene_count[gene], gene),
        )
        if len(recurrent) < 2:
            continue
        output.append({
            "program_id": f"unmodeled_lineage_candidate_{len(output) + 1:03d}",
            "resolutions": ";".join(resolutions),
            "source_program_ids": ";".join(
                sorted(str(rows[index]["program_id"]) for index in members)
            ),
            "genes": ";".join(recurrent),
            "coexpressed_gene_count": len(recurrent),
            "spatially_coherent": "true",
            "excluded_program_classes": "",
            "catalog_match": "",
            "status": "Unmodeled lineage candidate",
        })
    args.out.mkdir(parents=True, exist_ok=True)
    table = args.out / "unmodeled_lineage_candidates.tsv"
    write_tsv(
        table,
        output,
        fields=(
            list(output[0]) if output else [
                "program_id", "resolutions", "source_program_ids", "genes",
                "coexpressed_gene_count", "spatially_coherent",
                "excluded_program_classes", "catalog_match", "status",
            ]
        ),
    )
    excluded_table = args.out / "excluded_state_or_technical_programs.tsv"
    write_tsv(
        excluded_table, excluded_rows,
        fields=[
            "program_id", "resolution", "excluded_program_class",
            "excluded_gene_fraction", "best_catalog_candidate_id",
            "best_catalog_overlap_gene_count",
            "best_catalog_overlap_family_count", "genes",
        ],
    )
    manifest = {
        "status": "PASS",
        "schema_version": "2.2",
        "input_programs": {
            "path": str(args.programs.resolve()),
            "sha256": sha256(args.programs),
        },
        "unmodeled_lineage_candidate_n": len(output),
        "accepted_program_n": len(output),
        "accepted_programs": {
            "path": str(table.resolve()),
            "sha256": sha256(table),
        },
        "excluded_state_or_technical_programs": {
            "path": str(excluded_table.resolve()),
            "sha256": sha256(excluded_table),
            "n_programs": len(excluded_rows),
        },
        "minimum_repeated_resolutions": 2,
        "formal_labels_written": False,
    }
    write_manifest(args.out / "unmodeled_discovery_manifest.json", manifest)
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
