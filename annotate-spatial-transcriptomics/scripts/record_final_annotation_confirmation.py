#!/usr/bin/env python3
"""Record explicit user approval of a frozen annotation snapshot."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root", type=Path)
    parser.add_argument("--confirmed-by", required=True)
    parser.add_argument("--user-message", required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    request_path = root / "provenance/final_annotation_confirmation_request.json"
    if not request_path.is_file():
        raise SystemExit("create the frozen confirmation request before recording approval")
    request = json.loads(request_path.read_text(encoding="utf-8"))
    if request.get("status") != "AWAITING_USER_CONFIRMATION":
        raise SystemExit("confirmation request has an invalid status")
    for key, hash_key in (
        ("cell_ledger", "cell_ledger_sha256"),
        ("cluster_ledger", "cluster_ledger_sha256"),
        ("completion_gate", "completion_gate_sha256"),
        ("release_taxonomy_audit", "release_taxonomy_audit_sha256"),
    ):
        path = root / request[key]
        if not path.is_file() or sha256(path) != request[hash_key]:
            raise SystemExit(f"annotation snapshot changed after the confirmation request: {key}")
    confirmation = {
        "status": "CONFIRMED",
        "confirmed_at": datetime.now(timezone.utc).isoformat(),
        "confirmed_by": args.confirmed_by,
        "user_message": args.user_message,
        "sample_id": request.get("sample_id"),
        "decision_version": request.get("decision_version"),
        "confirmation_request": "provenance/final_annotation_confirmation_request.json",
        "confirmation_request_sha256": sha256(request_path),
        "cell_ledger": request["cell_ledger"],
        "cell_ledger_sha256": request["cell_ledger_sha256"],
        "cluster_ledger": request["cluster_ledger"],
        "cluster_ledger_sha256": request["cluster_ledger_sha256"],
        "completion_gate": request["completion_gate"],
        "completion_gate_sha256": request["completion_gate_sha256"],
        "release_taxonomy_audit": request["release_taxonomy_audit"],
        "release_taxonomy_audit_sha256": request["release_taxonomy_audit_sha256"],
        "release_scope": "Final assets and HTML may now be generated for exactly this frozen annotation snapshot.",
    }
    out = root / "state/final_annotation_confirmation.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(confirmation, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(confirmation, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
