#!/usr/bin/env Rscript
# Fill serum_log2fc_PDAC_vs_HC in the shipped atlas with the pooled phase77
# log2FC wherever the phase42 value is missing (the rule build_reference.R
# now applies). Only that column changes: it lets strict TRACE-D see all 38
# serum-detected genes instead of 19.
#
# Run from the package root: Rscript data-raw/attach_phase77_serum_log2fc.R
suppressPackageStartupMessages(library(data.table))
PKG <- rprojroot::find_package_root_file()
load(file.path(PKG, "data", "pdactrace_reference.rda"))
p77 <- fread(cmd = paste("xz -dc", shQuote(file.path(PKG, "inst", "extdata",
             "phase77_strict_RNAprotConvergent_serum.csv.xz"))))
i <- match(p77$gene, pdactrace_reference$gene_symbol)
fill <- !is.na(i) & is.na(pdactrace_reference$serum_log2fc_PDAC_vs_HC[i])
set(pdactrace_reference, i[fill], "serum_log2fc_PDAC_vs_HC", p77$serum_logFC[fill])
save(pdactrace_reference, file = file.path(PKG, "data", "pdactrace_reference.rda"), compress = "xz")
cat(sprintf("filled %d genes; %d genes now carry a serum log2FC\n", sum(fill),
            sum(!is.na(pdactrace_reference$serum_log2fc_PDAC_vs_HC))))
