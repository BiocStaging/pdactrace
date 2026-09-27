# MSV000101183 v1.1 input schema

## Label-blind abundance

`inbox/abundance.csv` is a protein-by-sample matrix.

- First column: `feature`, using the frozen HGNC gene symbols in
  `frozen/prm_target_assays_exact.csv`.
- Exactly 93 assay rows are present. TGOLN2 is absent because its public
  precursor annotation is chemically inconsistent with the peptide sequence.
- Remaining columns: exactly 40 blinded sample IDs.
- Sample IDs are `SB` plus the first 12 uppercase hexadecimal characters of
  `SHA256(LOCK_ID + "\t" + raw_filename)`.
- Values are log2 baseline-subtracted product-ion peak areas.
- Missing values are empty CSV cells.
- The five `acc_pooled` files are assay controls and may not appear as patient
  columns.

The abundance matrix is generated without reading `Groups` from the sealed
metadata. Assay eligibility is fixed before outcome access: finite abundance in
at least 70% of the 40 biological samples and detectable in at least three of
the five pooled-QC files.

`preprocessing/work/technical_qc_decision.json` must be `PASS` under
`TECHNICAL_QC_RULE.json`. Otherwise the metadata remains sealed and the workflow
ends as `INCONCLUSIVE_TECHNICAL_QC`.

`frozen/assay_dictionary_panel_coverage.csv` records absolute-panel coverage
against the 93 valid public assay annotations. The later
`absolute_panel_coverage_prelabel.csv` records coverage after label-blind
eligibility. Neither table is allowed to rerank a method or alter PASS/FAIL.

## One-shot metadata

The public repository file has the locked header:

```text
No.,Raw files,Groups
```

It remains sealed until the one-shot wrapper runs. The wrapper writes
`inbox/metadata.csv` with exactly two columns:

```text
sample_id,condition
```

`condition` must be one of `Healthy`, `IPMC`, or `PDAC`, obtained only through
an exact normalized match in `OUTCOME_DICTIONARY.tsv`. Broad regular-expression
or substring rescue is forbidden.

For the structural `Raw files` field only, directory components are removed and
a bare integer such as `1` is normalized to `1.raw`; `1.raw` through `40.raw`
are the only accepted biological identifiers. This rule does not inspect or
transform the outcome field.

## Prelabel freeze

Before metadata is opened, label-blind preprocessing writes:

- `preprocessing/prelabel_freeze/assay_eligibility.csv`
- `preprocessing/prelabel_freeze/frozen_panels_prelabel.csv`
- `preprocessing/prelabel_freeze/absolute_panel_coverage_prelabel.csv`
- `preprocessing/prelabel_freeze/prelabel_manifest_sha256.csv`
- `preprocessing/prelabel_freeze/PRELABEL_ID.txt`

The prelabel manifest binds the abundance matrix, assay eligibility, conditional
panels, absolute-panel coverage, technical-QC result, raw download manifest,
extraction record, and active lock ID. A second external timestamp is required
before one-shot outcome access.

## Terminal outputs

The immutable runner writes under `results/` only. `decision.json` contains one
of `PASS`, `FAIL`, or `INCONCLUSIVE`. An exact-count or label-dictionary mismatch
is a terminal `INCONCLUSIVE`, not a technical exception and not a reason to
change the mapping.
