# Outcome-blind technical amendment

## Immutable predecessors

- Pre-data v1 LOCK_ID:
  `06392ef304074c3b83ea7ea93da5e0f28e016274e137af22f323757d98e82b4e`
- Pre-data v1.1 LOCK_ID:
  `710ad1251eb2cd7511bb53169430303f15ab3fd7208a32bc6448f18e5abc611d`
- v1.1 terminal status: `INCONCLUSIVE_TECHNICAL_QC`
- v1.1 technical-terminal RFC 3161 time: 2026-07-18 00:47:03 UTC

Version 1.1 remains terminal. Its thresholds are not relaxed and its result is
not relabeled as PASS, FAIL, or a biological validation.

## Information accessed before v2

The 45 public validation RAW files were downloaded and processed without sample
labels under v1.1. The metadata file remains sealed as 742 uninterpreted bytes
with SHA-256
`88294963f5cdcbaa6e5c9625956b0a270a95068434376ce665a9a6ecce0a45ae`.
No `Groups` value or sample-to-condition assignment was opened.

Version 1.1 mapped 20 precursors consistently but failed its product-ion
technical gate: 11 eligible assays, minimum 9 detected assays per run, median
pooled-QC CV 1.0059, and zero eligible assays with CV at most 0.30.

## Root cause and fixed change

All acquired target scan filters are low-resolution ion-trap CID scans beginning
with `ITMS + c ESI Full ms2`. A 20-ppm fragment window is therefore inconsistent
with the product analyzer. The original pre-data v1 protocol had used the
source-study product m/z with a 0.5-Da half-window.

Before v2, that already prespecified v1 window was rerun on the five pooled-QC
files only. All 20 acquired targets were detected in 5/5 pooled runs, median CV
was 0.1369749, and 19/20 targets had CV at most 0.30. No broad-window biological
abundance was extracted before this lock.

Version 2 therefore changes only product-ion extraction:

1. precursor: unique theoretical peptide precursor, 20 ppm;
2. product center: public source-study reported product m/z;
3. product half-window: 0.5 Da, converted to ppm for rawrr;
4. integration and baseline: unchanged from v1 and v1.1.

The 20 acquired precursor targets are fixed from scan metadata shared by all 45
runs. TGOLN2 is not among them. The outcome dictionary, exact count gate,
candidate scores, comparators, direct orientation, feature endpoint, patient
endpoint, and decision thresholds are byte-identical to v1.1.

## Evidential status

Version 2 is outcome-blind but technically data-informed. It must be reported as
such. A v2 PASS would support transport within an independently selected
targeted panel under a blinded amendment; it would not erase the v1.1 technical
failure or constitute fully prospective preprocessing.
