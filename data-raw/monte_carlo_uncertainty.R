#!/usr/bin/env Rscript
# Monte Carlo uncertainty summary for the bundled atlas.
#
# Writes data-raw/audit_score_mc_v1.csv, which build_audit_scores_v2.R
# attaches to data/pdactrace_reference.rda as the audit_score_median,
# audit_score_lo95/hi95, audit_rank_* and audit_confidence_class columns.
#
# The perturbation loop is the package's own .audit_mc_table(), so the
# stored summary is exactly propagate_uncertainty(genes, n_mc = 500,
# seed = 42) and uses the same 3-axis formula as compute_audit_score().
# Only template rho, cohort agreement and I2 are perturbed; weights and
# the leakage gate are not.
#
# Run from the package root, then Rscript data-raw/build_audit_scores_v2.R
suppressPackageStartupMessages(library(data.table))
PKG <- rprojroot::find_package_root_file()
pkgload::load_all(PKG, quiet = TRUE)

load(file.path(PKG, "data", "pdactrace_reference.rda"))
mc <- .audit_mc_table(pdactrace_reference, n_mc = 500L, seed = 42)
fwrite(mc, file.path(PKG, "data-raw", "audit_score_mc_v1.csv"))
print(table(mc$confidence_class))
