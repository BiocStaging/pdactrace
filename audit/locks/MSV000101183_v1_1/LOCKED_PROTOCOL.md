# MSV000101183 StageBridge prospective lock v1.1

## Status and separation

This is a prospectively locked reanalysis of a retrospective public cohort, not
a prospective clinical trial. It is a new namespace and does not revise,
replace, rescue, or pool the terminal `PXD067770` result.

It supersedes the unopened `MSV000101183_v1` execution layer. The predecessor
remains immutable and independently timestamped. `SUPERSESSION.md` fixes the
scope of this change before validation RAW data or sample-level outcomes are
accessed.

At lock time, the outcome vocabulary and aggregate counts were publicly known:
Healthy n=20, IPMC n=4, and PDAC n=16. No sample-level row from the public
`metadata/meta.csv` had been opened. The metadata header, repository file
inventory, numeric run names, and aggregate file counts had been inspected.

## Scientific question

The primary question is whether the previously frozen StageBridge PDAC ranking
transports to an independent, targeted serum PRM validation cohort for PDAC
versus healthy control. The cohort was not used to train StageBridge.

The study does not contain chronic pancreatitis or another inflammatory control.
It therefore cannot establish inflammatory specificity. It is a same-cancer,
same-specimen, independent-cohort test, not pan-cancer validation.

Because the assay targets were selected by the source study, all conclusions are
conditional on the deposited targeted panel. The result must be described as a
locked targeted-panel audit rather than an untargeted proteome-wide validation.

## Public outcome information allowed before lock

The following public information was intentionally known and is frozen:

1. Clinical meaning of HC, IPMC, and PDAC.
2. Aggregate counts Healthy 20, IPMC 4, and PDAC 16.
3. All validation-cohort PDAC and IPMC cases were reported as pathologically
   confirmed.
4. The source article's validation cohort was reported as independent and not
   used for source-study feature selection or model training.

The original article's reported classifier performance is not an input to any
StageBridge threshold, score, preprocessing choice, or decision criterion.

## Frozen discovery objects

- Candidate universe: 1,356 existing PDAC StageBridge candidates.
- StageBridge score and seven comparator scores: unchanged from the prior sealed
  discovery analysis.
- Tissue direction: frozen signed-Stouffer direction across discovery RNA
  `early_vs_normal` effects.
- Destination relation: direct concordant only.
- Weights: unchanged default StageBridge weights.
- Absolute panels: top 5, 10, and 20 for every method before assay inspection.
- Source-study PRM dictionary: 94 public target assays deposited with the paired
  discovery-cohort repository, fixed by gene symbol, accession, peptide,
  precursor m/z, and quantitative product m/z. Ninety-three annotations map
  uniquely to theoretical precursor and b/y product ions. TGOLN2 is excluded
  before validation data access because its reported precursor is inconsistent
  with its peptide sequence.

No inverse orientation, validation-driven weight change, feature deletion,
label-based target selection, or endpoint substitution is allowed.

## Label-blind preprocessing after external timestamp

Only after the lock disclosure verifies against an independent RFC 3161
timestamp may the workflow download or inspect quantitative files.

1. Download numbered biological RAW files `1.raw` through `40.raw`, pooled-QC
   files `acc_pooled1.raw` through `acc_pooled5.raw`, and the metadata file.
2. Store metadata as sealed bytes; compute its hash without parsing rows.
3. Assign blinded sample IDs from LOCK_ID and raw filename.
4. Match each of the 93 valid theoretical precursors within 20 ppm and extract
   its uniquely mapped theoretical product-ion XIC at 20 ppm. Apply
   10th-percentile baseline subtraction, trapezoidal integration over the
   scheduled trace, and log2 transformation of positive areas.
5. Exclude pooled-QC files from biological abundance. They may be used only for
   prespecified label-blind detection and reproducibility gates.
6. Mark an assay eligible only if abundance is finite in at least 70% of the 40
   biological samples and in at least three of five pooled-QC files.
7. Require all label-blind technical gates in `TECHNICAL_QC_RULE.json`: exact
   run counts; at least 20 assays mapped in at least 90% of biological runs; at
   least 15 eligible assays; at least 20 mapped and 15 detected assays in every
   run; median eligible-assay pooled-QC CV at most 0.30; and at least 80% of
   eligible assays with pooled-QC CV at most 0.30. Failure is
   `INCONCLUSIVE_TECHNICAL_QC` and forbids outcome access.
8. Select each method's assay-conditional top 5, 10, and 20 using frozen scores
   and directions only. Separately tabulate availability of the absolute top 5,
   10, and 20 panels first against the public assay dictionary and then against
   label-blind assay eligibility; this coverage is reported and cannot alter the
   decision rule. The public dictionary contains zero proteins from the absolute
   StageBridge top-10, so this cohort is explicitly a panel-conditional transport
   audit, not a validation of the original absolute top-10 signature.
9. Commit abundance, eligibility, conditional panels, absolute-panel coverage,
   download hashes, extraction record, technical-QC result, and software
   identity to a prelabel manifest and obtain a second external timestamp.

No group-specific imputation, filtering, normalization, peak selection, batch
correction, or manual chromatogram editing is permitted. The analysis runner
replaces missing abundance by the feature minimum minus 0.5 observed standard
deviations, then applies feature-wise robust scaling and clips values to [-8, 8].
R, Python, package versions, rawrr assembly version, and the rawrr executable
SHA-256 must equal `frozen/execution_environment.json` at every executable stage.

## One-shot outcome opening and count gate

The metadata file may be opened once only after both external timestamps and all
prelabel hashes verify. Raw group values are normalized by lowercasing, replacing
non-alphanumeric runs with one space, and trimming. The result must match an
entry in `OUTCOME_DICTIONARY.tsv` exactly.

The structural raw-file field is reduced to its basename; a bare integer is
normalized by appending `.raw`. Only `1.raw` through `40.raw` are accepted. This
predeclared filename normalization does not use the outcome field.

The observed condition counts must equal the frozen aggregate table exactly:
Healthy 20, IPMC 4, and PDAC 16. Unknown aliases, duplicate raw filenames,
missing biological files, or count disagreement produce a terminal
`INCONCLUSIVE` with reason `OUTCOME_SCHEMA_OR_COUNT_MISMATCH`. The mapping may not
be expanded or corrected after opening.

## Feature endpoint

The primary biological contrast is PDAC versus Healthy. Among assay-eligible
frozen candidates, a feature is positive only if:

1. limma robust empirical-Bayes BH-FDR is below 0.10; and
2. the observed PDAC-minus-Healthy effect has the frozen direct tissue sign.

Primary ranking metric: average precision on the strict all-method common
universe. Secondary ranking metrics: AUROC and top-20 recall. Full method-specific
coverage and non-common metrics are reported but do not replace the common-
universe decision metric.

IPMC versus Healthy is descriptive because n=4. A pooled invasive endpoint
(PDAC plus IPMC versus Healthy) is prespecified as supportive and cannot rescue
or overturn the primary result.

## Patient endpoint

The primary patient endpoint is PDAC versus Healthy using each method's frozen
direct-orientation assay-conditional top-10 signature. Feature rank weights are
`1/log2(rank + 1)` and no coefficient is trained.

AUROC and average precision are reported with 2,000 stratified bootstrap
resamples using seed 20260717. Pooled invasive versus Healthy and IPMC versus
Healthy are secondary and descriptive, respectively. No threshold, calibration
model, covariate model, or post hoc panel is fitted.

## Decision rule

After the exact count gate, the result is evaluable only when:

1. At least five positive features exist on the all-method common universe.
2. At least five proteins are available in the StageBridge top-10 signature.

If either condition fails, the result is `INCONCLUSIVE`. Otherwise, the result is
`PASS` only when every condition below holds:

1. StageBridge common-universe average precision is at least 0.05 above endpoint
   prevalence.
2. StageBridge common-universe average precision is at least 0.05 above the best
   prespecified comparator.
3. The lower 95% bootstrap AUROC bound for the StageBridge top-10 PDAC-versus-
   Healthy score exceeds 0.50.

An evaluable result missing any PASS condition is `FAIL`. IPMC performance is
not required for PASS because its prespecified aggregate size is four.

## Technical failure policy

A failed label-blind assay or pooled-QC gate is a terminal technical
inconclusive for this lock and outcomes remain sealed. It may not be repaired by
relaxing extraction tolerances, changing targets, or weakening thresholds.

A network interruption, dependency failure, disk error, or process crash is not
a scientific result. It is recorded as `ONE_SHOT_TECHNICAL_FAILURE.json`. Any
recovery must use byte-identical metadata, abundance, lock, prelabel artifacts,
and runner; it may not repeat mapping choices or change rules. A schema/count
mismatch is not technical and must remain terminal `INCONCLUSIVE`.

## Interpretation ceiling

A PASS supports independent targeted-panel transport for PDAC versus healthy
serum. It does not establish clinical utility, early-detection performance,
inflammatory specificity, IPMC performance, or cross-cancer generalization. A
matched benign inflammatory cohort and a separately locked second cancer remain
necessary for a strong Briefings in Bioinformatics generalization claim.
