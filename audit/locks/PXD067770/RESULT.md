# PXD067770 locked validation result

- Lock ID: `fe595c20d577b0ddea62c1a622a3be65d110e319e89d21e79a72700855927b3a`
- Decision: **INCONCLUSIVE**
- Decision origin: terminal one-shot technical adjudication
- Outcome access started: 2026-07-17 12:30:08 UTC
- Statistical runner executed: no
- Performance metrics available: no
- Rerun or rescue analysis in this locked result: forbidden

## Reason

The source disease labels mapped to the three allowed conditions, but the
resulting counts did not equal the externally timestamped expectation of 12
PDAC-condition, 36 IPMN-condition and 12 Healthy samples. The locked endpoints
could not be instantiated without changing the mapping or accessing additional
outcome information.

The locked protocol assigns `INCONCLUSIVE` when evaluable group information is
insufficient. The mismatch occurred before the immutable statistical runner was
invoked, so no feature- or patient-level performance estimate exists and no
claim about StageBridge transport performance is made.

Any corrected mapping or alternative endpoint requires a new public lock and a
separately versioned analysis. It must not overwrite this result.
