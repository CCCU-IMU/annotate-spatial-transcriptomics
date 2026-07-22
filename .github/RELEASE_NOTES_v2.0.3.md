# annotate-spatial-transcriptomics v2.0.3

This patch aligns terminal controller state and final reports without changing the v2.0.2 biological writeback rules.

- Terminal ledgers may bind cells to active direct returns, active Atlas routes or typed residual-QC freeze decisions; zero-count Atlas partitions remain valid and auditable.
- Completion, autopilot, report and release audit share biological-context alias handling and recognize canonical v2 final-state/report asset names.
- Residual-QC validation auto-detects `final_state` and `qc_reason` while retaining explicit column overrides.
- Every released broad/fine spatial highlight now displays its validated route, quantitative support, competing evidence, spatial review, evidence source and current-data top DEGs.
- All-cell marker spatial assets require the exact analysis-set scope and can fail early on an observation-count mismatch.
- No new public completion gate was added.
