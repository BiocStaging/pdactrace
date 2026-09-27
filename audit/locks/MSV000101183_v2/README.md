# MSV000101183 StageBridge outcome-blind lock v2

This namespace is a separately disclosed technical amendment created after the
terminal label-blind failure of `MSV000101183_v1_1` and before any sample-level
outcome row was opened. It does not overwrite or reinterpret v1.1.

The amendment fixes one instrument-resolution mismatch. Validation MS2 spectra
were acquired as low-resolution `ITMS` CID scans, whereas v1.1 had required a
20-ppm product-ion XIC. Version 2 restores the 0.5-Da product half-window that
was already frozen in the pre-data v1 protocol. Exact theoretical precursor
matching remains at 20 ppm. No biological endpoint, label mapping, comparator,
score, direction, or PASS/FAIL threshold changes.

The broad product window was checked only on the five label-blind pooled-QC
runs before this lock. Broad-window biological abundance was not extracted, and
the sealed metadata bytes were not parsed. The resulting workflow is an
outcome-blind locked validation, not a claim of preprocessing fixed before RAW
access.

Key files:

- `AMENDMENT.md`: full version history, failure, diagnosis, and allowed change.
- `PRODUCT_WINDOW_RULE.json`: executable mass-window specification.
- `DECISION_RULE.json`: unchanged biological decision rule.
- `TECHNICAL_QC_RULE.json`: unchanged pre-outcome technical gate.
- `LOCKED_PROTOCOL.md`: complete v2 workflow and interpretation ceiling.
- `PREOUTCOME_ACCESS_LOG.tsv`: exact data-access boundary before v2.
- `frozen/LOCK_ID.txt`: SHA-256 commitment over the immutable manifest.

The public assay dictionary contains none of the absolute StageBridge top-20.
Consequently, this cohort tests ranking and direct orientation within an
independently selected 20-target panel; it cannot validate the original
absolute StageBridge signature.
