# v1.9.1 — Prelabel ledger preservation fix

This patch release fixes the reviewed-mapping writer so the v1.9 prelabel evidence artifact, SHA-256, frozen flag, winner, runner-up and winning margin remain intact in `cluster_decision_ledger.tsv`.

It prevents a state commit from silently downgrading the hallucination-audit boundary introduced in v1.9.0.
