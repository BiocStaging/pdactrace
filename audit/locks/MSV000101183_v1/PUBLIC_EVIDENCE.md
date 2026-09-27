# Public evidence checked before lock

Checked UTC: 2026-07-17T14:01:02Z

## Cohort and outcomes

The open-access article reports an independent validation cohort of 40 serum
participants: healthy control (HC), n=20; IPMN-associated invasive carcinoma
(IPMC), n=4; and pancreatic ductal adenocarcinoma without IPMN (PDAC), n=16.
These counts appear in Table 1 and again in the Figure 6 description.

The article defines IPMC as IPMN-associated invasive carcinoma and distinguishes
PDAC without intraductal papillary mucinous neoplasm from IPMC. It also states
that all validation-cohort PDAC and IPMC diagnoses were pathologically
confirmed and that the validation cohort was not used for feature selection or
model training.

- Article: https://doi.org/10.1002/cam4.71941
- Full text: https://onlinelibrary.wiley.com/doi/full/10.1002/cam4.71941
- Count locator: Table 1, Validation Cohort; Figure 6 legend
- Definition locator: Sections 1 and 2.1; Table 1 footnotes

## Repository and assay

The article identifies `MSV000101183` as the public validation-set PRM dataset.
The public MassIVE page describes 40 human serum samples and PRM LC-MS/MS,
marks the dataset Public, and provides dataset DOI `10.25345/C5ST7F99B`.

The label-free repository inventory contained 93 files totaling 1.21 GB:
40 numbered biological RAW files, 40 corresponding mzML files, five pooled-QC
RAW files, five corresponding pooled-QC mzML files, one metadata CSV, one
parameter XML, and one repository summary TSV. Biological filenames are numeric
(`1` through `40`) and do not disclose the condition.

- Dataset page: https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000101183
- Dataset DOI: https://doi.org/10.25345/C5ST7F99B
- MassIVE task: `6aab5886577844e098bea31f835c4954`
- Public FTP root advertised by MassIVE: `ftp://massive-ftp.ucsd.edu/v12/MSV000101183/`
- Repository summary at check: 93 files, 1.21 GB, 80,301 spectra

MassIVE reports zero standardized experimental-design conditions because the
submission uses a small `metadata/meta.csv` rather than a repository-standard
condition table. Therefore the per-condition totals are committed from the
peer-reviewed article, while the repository independently confirms the total
40 biological samples and exposes the sample-to-group file for the later
one-shot check.

The MassIVE free-text description generically mentions patients with IPMN,
whereas the article's exact validation table contains HC, IPMC, and PDAC only.
This discrepancy was known before lock. Noninvasive IPMN is deliberately absent
from the outcome dictionary; if it appears in the sample-level file, the exact
count gate returns terminal `INCONCLUSIVE` rather than silently pooling it.

## Metadata access boundary

Before lock, only the first header row of `metadata/meta.csv` was requested:

```text
No.,Raw files,Groups
```

No sample-level metadata row, raw group value, sample-to-group assignment,
quantitative PRM result, or validation performance was inspected. Search-engine
excerpts from the publication necessarily exposed the already published
aggregate counts and the original authors' reported validation performance;
those published performance values are not used in the StageBridge decision
rule or preprocessing.
