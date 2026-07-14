#!/usr/bin/env python3
"""Bind the user's actual approval to the frozen cohort request."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("cohort_root", type=Path)
    parser.add_argument("--user-message", required=True)
    args = parser.parse_args()
    if len(args.user_message.strip()) < 2:
        raise SystemExit("the actual user confirmation message is required")
    root = args.cohort_root.resolve()
    request_path = root / "provenance/cohort_confirmation_request.json"
    request = json.loads(request_path.read_text(encoding="utf-8"))
    if request.get("status") != "COHORT_CONFIRMATION_PENDING":
        raise SystemExit("cohort confirmation request is not pending")
    confirmation = {"status": "CONFIRMED", "request_path": str(request_path.resolve()), "request_sha256": sha(request_path), "user_message": args.user_message, "confirmed_at": datetime.now(timezone.utc).isoformat()}
    out = root / "state/cohort_confirmation.json"
    out.write_text(json.dumps(confirmation, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(confirmation, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
