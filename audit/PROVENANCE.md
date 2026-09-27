# Provenance of the frozen PDAC score layer and one-shot audit records

## Where these files came from

The frozen score table, the absolute panels and the two one-shot audit records
in this repository were produced in a separate working repository that carried a
wider scope (PDAC + LUAD + CRC + HCC) under the working name `stagebridge`.
Every file that tree tracks is preserved unmodified at:

| | |
|---|---|
| Tag | `v0.1.0-pancancer-archive` |
| Tag object | `e268771c54b681c1984b50633df71ea8a942ac75` |
| Commit | `bc077e6d22ad8b583e2e2124beb78ed577e352eb` |
| Frozen | 2026-07-20 |

Only the PDAC subset was carried into `pdactrace`. The LUAD, CRC and HCC
analyses, and the three non-PDAC one-shot audits, remain in the archive and are
out of scope here.

## Why the archive must not be flattened into this repository

Each lock carries a `frozen/lock_manifest_sha256.csv` whose `path` column is
**relative to the root of the archived tree** (`DESCRIPTION`, `NAMESPACE`,
`R/*.R`, `prospective/locks/...`). Those manifests are covered by external
RFC 3161 timestamps (DigiCert) and OpenTimestamps attestations.

Consequently:

- Lock verification resolves **only inside the archive**, where the recorded
  paths still exist. It is not expected to resolve against this repository.
- The manifests must never be regenerated. Reissuing a lock destroys the
  property the lock exists to establish — that the candidate set, scores,
  directions and decision rule were fixed before any outcome label was opened.
- The copies under `audit/locks/` are a **disclosure record**, not a verifiable
  artefact in place. To verify, check out the tag above.

## How to verify, and one known discrepancy

Verification was last run on 2026-09-27 against the archive as tagged. All nine
`lock_manifest_sha256.csv` files resolve, with **one caveat that anyone
repeating the check will hit immediately**.

Five of the nine manifests hash a root-level `NAMESPACE` of 278 bytes
(`db4c580afbbf64ab9685526a4e841403f589531a9c5eadd152b126cc3c52f14d`), but the
`NAMESPACE` at the tagged commit is 331 bytes. Exports were added after those
five locks were frozen, roxygen rewrote the file, and the tag captured the
later state. Run naively, those five manifests therefore report exactly one
mismatched path each, and that path is always `NAMESPACE`.

The 278-byte file itself survives in the archive, tracked by git, inside a
later lock's frozen snapshot:

```
prospective/locks/PXD046295_hcc_discovery_precoverage_v1/frozen/snapshots/project/NAMESPACE
```

Substituting it for the root `NAMESPACE` makes all nine manifests verify with
zero mismatches. Every other hashed path — decision rules, candidate sets,
frozen scores, directions, `R/*.R` — matches as recorded, in all nine
manifests, without substitution.

This is a drift in a generated file, not in any frozen decision artefact, and
it is recorded here rather than repaired: editing the archive to make the check
pass cleanly would be the one thing the protocol forbids.

## What the archive still contains

On 2026-09-27 the archive was pruned from 5.3 GB to 72 MB by deleting 280 files
of third-party bulk data — vendor mass-spectrometry `.raw` files, DIA-NN report
bodies, `proteinGroups.txt`, GEO series matrices — together with build
artefacts. All of it is re-obtainable from PRIDE, MassIVE and GEO, and none of
it was redistributable in the first place; the archive's own `.gitignore`
excluded it for that reason.

Nothing tracked by git was removed, and the pruning did not touch four
untracked files that a manifest hashes
(`prospective/locks/MSV000101183_v1_1/preprocessing/...`, 20 KB in total). The
nine-manifest verification above was re-run after the pruning with the same
result. A `git bundle` of all refs, including the tag, is held separately.

## Files copied, with content hashes at the time of transfer

Hashes are SHA-256 of the file as copied into this repository.

```
8aae1dbd91c760bf48da840c77508f2affc0a32be6f30702fb5179ca512b3afc  data-raw/stagebridge/pdac_effects.csv
02a343bafd42ca975c607d67904a9b4d2107258d517875f94decd1736683bd8d  data-raw/stagebridge/frozen/pdac_frozen_score_table.csv
7471734c83db3e4cb8df3cc2fca8a77feeaa90d150209ffb949f410e83c5c982  data-raw/stagebridge/frozen/pdac_absolute_panels.csv
2d4781b0779d85e5f15f699575322030848fec20733d2d56e271c4dc3371406b  data-raw/stagebridge/pdac_detection_holdout.csv
64aa34361d4e42d9aa91cd992f952824dcd7eae589c83b875ce0cfa746c5ec90  data-raw/stagebridge/pdac_validation_labels.csv
```

`pdac_frozen_score_table.csv` holds 1,356 features with the frozen score,
`frozen_direction`, `target_relation` and `relation_status`. The 1,356 features
are the Early-only RNA x protein intersection of the atlas in `data/`, so the
frozen layer is a derived view over `pdactrace_reference`, not an independent
dataset.

## The two PDAC one-shot audits, and their outcome

Both were prespecified, externally timestamped, opened once, and terminated
without computing a performance statistic.

| Lock | Dataset | Endpoint | Terminal state |
|---|---|---|---|
| `MSV000101183` (v1, v1_1, v2) | MassIVE MSV000101183 serum | PDAC vs healthy control | `INCONCLUSIVE_TECHNICAL_QC` — prespecified label-blind technical QC failed |
| `PXD067770` | PRIDE PXD067770 serum | PDAC vs control | `INCONCLUSIVE` — mapping/coverage gate mismatch |

Neither audit reached an AUC or a p-value. This is the designed behaviour of a
terminal inconclusive state, not a missing analysis, and it must be reported
alongside any use of the frozen scores rather than omitted.

The recurring cause in both cases was assay coverage: the frozen candidate set
could not be detected at the required depth in undepleted public serum
proteomes. This is consistent with the published behaviour of the plasma and
serum matrix, where only a small fraction of tumour tissue proteins is
observable and secreted proteins are detected far more readily than
intracellular ones. The two terminations are therefore treated as evidence
about the feasibility of public-deposit transfer, and as the documented
justification for collecting a purpose-designed cohort, rather than as failed
attempts to be set aside.
