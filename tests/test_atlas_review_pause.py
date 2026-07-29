from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "annotate-spatial-transcriptomics/scripts/validate_global_atlas_v2.py"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def artifact(path: Path) -> dict[str, str]:
    return {"path": str(path.resolve()), "sha256": sha(path)}


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


class AtlasReviewPauseTests(unittest.TestCase):
    def test_nonempty_queue_without_decisions_is_controlled_pause(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            placeholders = []
            for name in ("ledger", "mapping", "calibration", "profile"):
                path = root / f"{name}.txt"
                path.write_text(name + "\n", encoding="utf-8")
                placeholders.append(path)
            routing = root / "routing.tsv"
            write_tsv(routing, [{
                "cell_id": "c1", "atlas_state_route": "defined_label_disagreement",
                "writeback_status": "proposal_only_requires_atomic_commit",
                "fine_anchor_eligible": "false", "review_required": "true",
                "review_id": "atlas_review_1", "primary_broad": "Granulosa",
                "proposed_broad_label": "Granulosa",
            }])
            queue = root / "queue.tsv"
            write_tsv(queue, [{"review_id": "atlas_review_1"}])
            manifest = root / "routing.json"
            manifest.write_text(json.dumps({
                "schema_version": "2.0",
                "authoritative_router": "route_global_atlas_v2.py",
                "ledger_writeback_performed": False,
                "fine_anchor_eligible": False,
                "n_analysis_set": 1, "status": "REVIEW_REQUIRED",
                "cell_ledger": artifact(placeholders[0]),
                "atlas_mapping": artifact(placeholders[1]),
                "calibration_manifest": artifact(placeholders[2]),
                "workflow_profile": artifact(placeholders[3]),
                "artifacts": {
                    "routing": artifact(routing), "review_queue": artifact(queue),
                },
            }), encoding="utf-8")
            out = root / "validation.json"
            result = subprocess.run([
                sys.executable, str(SCRIPT), "--routing-manifest", str(manifest),
                "--out", str(out),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertEqual(json.loads(out.read_text())["status"], "REVIEW_REQUIRED")


if __name__ == "__main__":
    unittest.main()
