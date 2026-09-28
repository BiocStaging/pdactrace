#!/usr/bin/env Rscript
# Recompute the seven audit feature columns of data/pdactrace_reference.rda
# (audit_score_layer, _direction, _early, _serum, _rescue, the two gate
# multipliers and the flags) from the atlas's current columns, using the same
# code that scores a user evidence table (.audit_compute_features()).
#
# Before 0.99.25 these columns came from a v0.3.0 CSV built on the 8-template
# labels, so 286 genes whose best template became Monotonic or Late in the
# 12-template catalog still carried Early-onset credit.
#
# Pipeline, from the package root:
#   Rscript data-raw/attach_audit_columns.R      # features
#   Rscript data-raw/monte_carlo_uncertainty.R   # MC summary from the features
#   Rscript data-raw/build_audit_scores_v2.R     # 3-axis score, class, MC columns
suppressPackageStartupMessages(library(data.table))
PKG <- rprojroot::find_package_root_file()
pkgload::load_all(PKG, quiet = TRUE)

ref_path <- file.path(PKG, "data", "pdactrace_reference.rda")
load(ref_path)
feat <- .audit_compute_features(pdactrace_reference)
stopifnot(identical(feat$gene_symbol, pdactrace_reference$gene_symbol))

col_map <- c(score_layer = "audit_score_layer", score_direction = "audit_score_direction",
             score_early = "audit_score_early", score_serum = "audit_score_serum",
             score_rescue = "audit_score_rescue", positive_score = "audit_positive_score",
             leakage_mult = "audit_leakage_mult", het_mult = "audit_heterogeneity_mult",
             audit_score_raw = "audit_score_raw", audit_score = "audit_score",
             is_hk = "audit_is_housekeeping", is_plasma_hi = "audit_is_plasma_high_abundance",
             rescue_eligible = "audit_rescue_eligible")
changed <- 0L
for (src in names(col_map)) {
  dst <- col_map[[src]]
  value <- feat[[src]]
  if (is.integer(pdactrace_reference[[dst]])) value <- as.integer(value)
  changed <- changed + sum(value != pdactrace_reference[[dst]], na.rm = TRUE)
  set(pdactrace_reference, j = dst, value = value)
}
save(pdactrace_reference, file = ref_path, compress = "xz")
message(sprintf("Recomputed audit feature columns (%d cell changes) in %s", changed, ref_path))
