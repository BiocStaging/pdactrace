# MSV000101183 v2 input schema

`inbox/abundance.csv` contains exactly 20 acquired assay rows and 40 blinded
biological sample columns. The first column is `feature`; sample IDs are inherited
from the v1.1 LOCK_ID-derived run map. Values are log2 product-ion peak areas and
missing values are empty cells. Pooled-QC runs never appear as patient columns.

The acquired assay dictionary is `frozen/acquired_target_assays.csv`. Each row
must map to one precursor schedule shared by all 45 RAW files. The product XIC
uses `PRODUCT_WINDOW_RULE.json`.

Before outcome access, `preprocessing/prelabel_freeze/` contains:

- `assay_eligibility.csv`
- `frozen_panels_prelabel.csv`
- `absolute_panel_coverage_prelabel.csv`
- `prelabel_manifest_sha256.csv`
- `PRELABEL_ID.txt`

The prelabel manifest also binds the v1.1 RAW manifest, sealed metadata bytes,
v2 abundance, extraction record, run-level QC, and technical decision.

The one-shot wrapper writes `inbox/metadata.csv` with `sample_id,condition` only
after both timestamps verify. Conditions must be exact dictionary mappings and
must match the frozen 20/4/16 aggregate table.
