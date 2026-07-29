from __future__ import annotations

import csv
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT / "annotate-spatial-transcriptomics/scripts/bind_atlas_routing_mapping.py"
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_table(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class AtlasSchemaBindingTests(unittest.TestCase):
    def test_optional_routing_fields_are_added_before_union_is_routed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            base = {
                "predicted_label": "Granulosa", "mapping_tier": "high",
                "out_of_distribution": "false", "ontology_conflict": "false",
            }
            target = root / "target.tsv"
            heldout = root / "heldout.tsv"
            combined = root / "combined.tsv"
            write_table(target, [{"cell_id": "q1", **base}])
            write_table(heldout, [{"cell_id": "h1", **base}])
            write_table(combined, [
                {"cell_id": "q1", **base}, {"cell_id": "h1", **base},
            ])
            calibration = root / "calibration.json"
            calibration.write_text(json.dumps({
                "artifacts": {"query_mapping": {
                    "path": str(target.resolve()), "sha256": sha(target),
                }},
            }), encoding="utf-8")
            out = root / "bound.json"
            result = subprocess.run([
                sys.executable, str(SCRIPT),
                "--calibration-manifest", str(calibration),
                "--heldout-mapping", str(heldout),
                "--combined-mapping", str(combined), "--out", str(out),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
            document = json.loads(out.read_text())
            normalized = Path(document["artifacts"]["query_mapping"]["path"])
            with gzip.open(normalized, "rt", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["fine_anchor_eligible"] for row in rows}, {"false"})
            self.assertEqual({row["review_required"] for row in rows}, {"false"})
            self.assertEqual({row["cell_id"] for row in rows}, {"q1", "h1"})


if __name__ == "__main__":
    unittest.main()
