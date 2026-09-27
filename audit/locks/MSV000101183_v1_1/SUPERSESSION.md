# Supersession statement

`MSV000101183_v1_1` supersedes `MSV000101183_v1` before any validation RAW
file, abundance matrix, sealed sample-level metadata file, or outcome row was
downloaded into either namespace.

The predecessor remains immutable and independently timestamped. It is not
deleted, rewritten, or treated as a failed scientific analysis. Version 1.1
changes only the pre-outcome technical execution layer:

1. reported one-decimal product ions are mapped to unique theoretical b/y ion
   exact masses and extracted at 20 ppm rather than a 0.5 Da half-window;
2. the inconsistent TGOLN2 precursor annotation is excluded before data access;
3. R, Python, package, rawrr assembly version, and assembly SHA-256 are pinned;
4. pooled-QC and per-run technical gates can terminate the workflow as
   `INCONCLUSIVE_TECHNICAL_QC` before outcomes are opened; and
5. executable lock-specific tests are included.

The outcome dictionary, aggregate counts, biological endpoints, discovery
scores, direct orientation, comparators, and PASS/FAIL criteria are unchanged.

Predecessor LOCK_ID:

```text
06392ef304074c3b83ea7ea93da5e0f28e016274e137af22f323757d98e82b4e
```
