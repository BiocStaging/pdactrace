# MSV000101183 StageBridge outcome-blind protocol v2

## Status

This is a locked reanalysis of a retrospective public cohort after label-blind
technical development and before outcome access. It is not a prospective trial
and not a preprocessing-before-RAW claim. `AMENDMENT.md` records the terminal
v1.1 failure and the sole allowed technical change.

At v2 lock time, public aggregate counts are Healthy n=20, IPMC n=4, and PDAC
n=16. The sealed metadata hash is known, but no sample-level outcome row has
been opened. Broad-window pooled-QC behavior is known; broad-window biological
abundance and all group comparisons are unknown.

## Scientific question and ceiling

The primary question is whether the frozen StageBridge PDAC ranking transports
within an independently selected serum PRM panel for PDAC versus healthy
control. IPMC is descriptive because n=4. The cohort has no pancreatitis control
and cannot establish inflammatory specificity or pan-cancer generalization.

The valid source dictionary and the acquired 20-target subset contain zero
proteins from the absolute StageBridge top-20. All conclusions are therefore
panel-conditional. The study does not validate the original absolute signature.

## Frozen scientific objects

- Candidate universe: 1,356 PDAC candidates.
- StageBridge and seven comparator scores: unchanged.
- Tissue direction: frozen signed-Stouffer discovery direction.
- Destination relation: direct concordant only.
- Outcome dictionary and aggregate count gate: byte-identical to v1.1.
- Biological decision rule: byte-identical to v1.1.
- Acquired assay subset: 20 precursor targets present in every one of 45
  label-blind RAW files and uniquely mapped to the public peptide dictionary.

No inverse orientation, outcome-driven feature deletion, endpoint substitution,
weight fitting, or label-based normalization is allowed.

## Locked extraction

1. Verify the v1.1 lock, terminal disclosure, RAW manifest, sealed metadata hash,
   and absence of any sample-level metadata output.
2. For each acquired assay, match its theoretical precursor within 20 ppm.
3. Extract the public reported product m/z with a 0.5-Da half-window from its
   exact `ITMS` scan filter.
4. Apply a 10th-percentile baseline, positive-part correction, trapezoidal
   integration, and log2 transformation of positive areas.
5. Mark an assay eligible when finite in at least 70% of 40 biological runs and
   detected in at least three of five pooled-QC runs.
6. Apply the unchanged `TECHNICAL_QC_RULE.json`. Failure keeps outcomes sealed.
7. Select assay-conditional top 5, 10, and 20 panels from frozen scores only;
   separately report absolute-panel coverage.
8. Freeze abundance, eligibility, panels, raw hashes, source sealed metadata,
   extraction identity, and technical QC in a prelabel manifest and obtain an
   independent external timestamp.

No manual peak editing, retention-time tuning, group-specific filtering,
imputation, normalization, or batch correction is permitted.

## One-shot outcome opening

Only after the v2 lock and prelabel disclosures independently verify may the
sealed metadata be parsed once. Labels are normalized only through
`OUTCOME_DICTIONARY.tsv`. Raw filenames must map exactly to `1.raw` through
`40.raw`, and counts must equal Healthy 20, IPMC 4, PDAC 16. Any alias, row, file,
or count mismatch is terminal `INCONCLUSIVE` without statistical analysis.

## Feature endpoint

Among eligible acquired candidates, a feature is positive when limma robust
empirical-Bayes BH-FDR for PDAC versus Healthy is below 0.10 and the effect has
the frozen direct tissue sign. Primary ranking performance is average precision
on the all-method common universe. AUROC and top-20 recall are secondary.

## Patient endpoint

The primary patient endpoint is PDAC versus Healthy using the assay-conditional
StageBridge top-10. Rank weights are `1/log2(rank + 1)` and no coefficient is
trained. AUROC and average precision use 2,000 stratified bootstrap resamples
with seed 20260717. Invasive versus Healthy is supportive; IPMC versus Healthy
is descriptive.

## Decision rule

The analysis is evaluable only with at least five common-universe positive
features and at least five StageBridge panel proteins. It is PASS only if:

1. StageBridge common-universe AP is at least prevalence + 0.05;
2. StageBridge AP is at least best prespecified comparator AP + 0.05; and
3. the lower 95% bootstrap bound for top-10 PDAC-versus-Healthy AUROC exceeds
   0.50.

An evaluable result missing any condition is FAIL. IPMC cannot rescue or overturn
the primary result.
