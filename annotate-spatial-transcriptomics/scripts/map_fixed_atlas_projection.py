#!/usr/bin/env python3
"""Project a frozen query matrix through a validated reference-only Atlas transform."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from sklearn.preprocessing import normalize


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True, type=Path)
    ap.add_argument("--cells", required=True, type=Path)
    ap.add_argument("--features", required=True, type=Path)
    ap.add_argument("--transform", required=True, type=Path)
    ap.add_argument("--prototypes", required=True, type=Path)
    ap.add_argument("--reference-heldout", required=True, type=Path)
    ap.add_argument(
        "--atlas-bundle-manifest", type=Path,
        default=(
            Path(__file__).resolve().parents[1]
            / "references/atlases/sheep_ovary_GSE233801_split_wall_v2.json"
        ),
    )
    ap.add_argument("--query-source", required=True, type=Path)
    ap.add_argument("--query-source-sha256-record", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    bundle = json.loads(args.atlas_bundle_manifest.read_text(encoding="utf-8"))
    canonical_bundle = (
        Path(__file__).resolve().parents[1]
        / "references/atlases/sheep_ovary_GSE233801_split_wall_v2.json"
    )
    if (
        args.atlas_bundle_manifest.resolve() != canonical_bundle.resolve()
        or sha256(args.atlas_bundle_manifest) != sha256(canonical_bundle)
        or bundle.get("bundle_id") != "sheep_ovary_GSE233801_split_wall_v2"
        or bundle.get("immutable") is not True
    ):
        raise SystemExit("projection requires the fixed GSE233801 Atlas bundle")
    assets = bundle.get("asset_hashes", {})
    for path, name in (
        (args.transform, "feature_transform.joblib"),
        (args.prototypes, "reference_prototypes.npz"),
        (args.features, "fixed_features.tsv"),
        (args.reference_heldout, "reference_heldout_predictions.tsv.gz"),
    ):
        if sha256(path) != str(assets.get(name, "")):
            raise SystemExit(f"fixed Atlas asset differs from bundle: {name}")

    if args.matrix.suffix == ".gz":
        with gzip.open(args.matrix, "rb") as handle:
            matrix_value = mmread(handle).tocsr().astype(np.float32)
    else:
        matrix_value = mmread(args.matrix).tocsr().astype(np.float32)
    cells = pd.read_csv(args.cells, sep="\t", dtype=str)["cell_id"].astype(str)
    features = pd.read_csv(args.features, sep="\t", dtype=str)["gene"].astype(str).tolist()
    model = joblib.load(args.transform)
    prototype_payload = np.load(args.prototypes)
    labels = prototype_payload["labels"].astype(str)
    prototypes = prototype_payload["prototypes"]
    if features != list(model["features"]):
        raise SystemExit("exported features differ from the fixed Atlas transform")
    if matrix_value.shape != (len(cells), len(features)) or not cells.is_unique:
        raise SystemExit("fixed-feature matrix/cell boundary is invalid")
    if matrix_value.data.size and (
        np.min(matrix_value.data) < 0 or not np.isfinite(matrix_value.data).all()
    ):
        raise SystemExit("query matrix contains negative or non-finite values")

    total = np.asarray(matrix_value.sum(axis=1)).ravel()
    normalized = sparse.diags(1e4 / np.maximum(total, 1)) @ matrix_value
    normalized.data = np.log1p(normalized.data)
    embedding = normalize(model["svd"].transform(normalized))
    similarity = embedding @ prototypes.T
    order = np.argsort(similarity, axis=1)
    top = order[:, -1]
    second = order[:, -2]
    mapped = pd.DataFrame({
        "cell_id": cells,
        "predicted_label": labels[top],
        "confidence": (similarity[np.arange(len(top)), top] + 1.0) / 2.0,
        "margin": similarity[np.arange(len(top)), top] - similarity[np.arange(len(top)), second],
        "top_similarity": similarity[np.arange(len(top)), top],
    })
    heldout = pd.read_csv(args.reference_heldout, sep="\t")
    required = {"predicted_label", "top_similarity"}
    if not required.issubset(heldout.columns):
        raise SystemExit("reference held-out predictions lack OOD fields")
    ood_threshold = heldout.groupby("predicted_label")["top_similarity"].quantile(0.01).to_dict()
    mapped["out_of_distribution"] = [
        similarity_value < ood_threshold.get(label, -1.0)
        for label, similarity_value in zip(mapped["predicted_label"], mapped["top_similarity"])
    ]
    mapped["ontology_conflict"] = False
    mapped["calibration_origin"] = "fresh_current_query_fixed_reference_projection"
    mapping_path = args.out / "all_cell_atlas_raw_mapping.tsv.gz"
    mapped.to_csv(mapping_path, sep="\t", index=False, compression="gzip")
    query_source_sha256 = sha256(args.query_source)
    if args.query_source_sha256_record:
        parts = args.query_source_sha256_record.read_text(encoding="utf-8").strip().split()
        if len(parts) < 2 or Path(parts[-1].lstrip("*")).resolve() != args.query_source.resolve():
            raise SystemExit("query-source sha256sum record is invalid or mismatched")
        query_source_sha256 = parts[0].lower()
    manifest = {
        "status": "PASS",
        "schema_version": "2.2",
        "mode": "fixed_reference_fresh_query_projection",
        "atlas_bundle": {
            "path": str(args.atlas_bundle_manifest.resolve()),
            "sha256": sha256(args.atlas_bundle_manifest),
            "bundle_id": bundle["bundle_id"],
            "reference_id": bundle["reference_id"],
        },
        "n_query": len(mapped),
        "n_features": len(features),
        "query_matrix": {"path": str(args.matrix.resolve()), "sha256": sha256(args.matrix)},
        "query_source": {"path": str(args.query_source.resolve()), "sha256": query_source_sha256},
        "feature_transform": {"path": str(args.transform.resolve()), "sha256": sha256(args.transform)},
        "reference_prototypes": {"path": str(args.prototypes.resolve()), "sha256": sha256(args.prototypes)},
        "mapping": {"path": str(mapping_path.resolve()), "sha256": sha256(mapping_path)},
        "query_reference_joint_retraining": False,
        "historical_query_labels_read": False,
        "tier_status": "raw_only_until_disjoint_current_query_anchor_calibration",
    }
    (args.out / "mapping_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
