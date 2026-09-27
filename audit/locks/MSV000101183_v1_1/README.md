# MSV000101183 StageBridge lock v1.1

This directory is an independent prospective-style lock for the public
validation cohort `MSV000101183`. It supersedes the unopened
`MSV000101183_v1` technical protocol while preserving that timestamped record,
and it does not modify or supersede the terminal `PXD067770` result.

At lock preparation, the public outcome vocabulary and aggregate condition
counts were known, but no sample-level row from `metadata/meta.csv` had been
opened. The only repository metadata content inspected was its header:
`No.,Raw files,Groups`.

The primary confirmatory endpoint is PDAC versus healthy control. IPMC has only
four participants and is therefore excluded from the PASS/FAIL rule; it is
retained as a descriptive endpoint and in a prespecified invasive-cancer pooled
secondary endpoint.

Key files:

- `OUTCOME_DICTIONARY.tsv`: exact, predeclared source-label aliases and mappings.
- `AGGREGATE_CONDITION_TABLE.tsv`: expected public group totals.
- `DECISION_RULE.json`: machine-readable PASS/FAIL/INCONCLUSIVE rule.
- `TECHNICAL_QC_RULE.json`: label-blind gate that must pass before outcomes open.
- `SUPERSESSION.md`: exact scope and predecessor commitment for version 1.1.
- `LOCKED_PROTOCOL.md`: complete scientific and operational protocol.
- `PUBLIC_EVIDENCE.md`: public sources inspected before lock.
- `frozen/lock_manifest_sha256.csv`: immutable-file commitment.
- `frozen/LOCK_ID.txt`: SHA-256 commitment over the manifest rows.
- `public_timestamp/`: independent RFC 3161 and OpenTimestamps evidence.

Version 1.1 maps the source-study one-decimal transition annotations to exact
theoretical ions, extracts at 20 ppm, excludes the chemically inconsistent
TGOLN2 annotation before data access, pins the execution environment, and
freezes absolute-panel assay coverage before labels are opened. These are
technical safeguards; the outcome mapping, endpoints, comparators, and decision
thresholds remain those of the unopened predecessor.

No result in this directory may overwrite or reinterpret the earlier
`PXD067770` one-shot outcome.
