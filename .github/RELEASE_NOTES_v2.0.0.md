# v2.0.0 — Contract-bound evidence and release integrity

This major release makes annotation decisions and release outputs fail closed against one immutable project contract.

Highlights:

- separates upstream BANKSY candidate grids from fresh whole-tissue and query-reclustering grids;
- requires a complete full-feature cluster-by-lineage-by-marker-family evidence matrix before broad naming;
- replaces overlapping Atlas adjudicators with one calibrated, state-aware all-cell router;
- enforces the `Vascular-associated` release hierarchy and rejects technical/QC semantics as biological fine labels;
- adds typed residual-QC and read-only result-directory audits;
- binds profiles, candidate catalogs, input ancestry, calibration, Atlas mapping and release ledgers by SHA256;
- provides explicit v1.10 migration without grandfathering historical labels or completion status.

Projects created under v1.x remain readable but do not automatically satisfy v2 release gates. Migrate, regenerate the new evidence artifacts, and rerun completion validation before publication.
