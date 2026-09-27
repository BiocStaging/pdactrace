#!/usr/bin/env Rscript
source(file.path(dirname(grep("^--file=", commandArgs(FALSE), value = TRUE) |> sub("^--file=", "", x = _)), "common.R"))

if (!requireNamespace("data.table", quietly = TRUE)) stop("Install data.table to prepare the PDAC data.")
cache <- sb_path("data-raw", "cache", "pdac")
out_dir <- sb_dir("data", "processed")

candidate_table <- data.table::fread(file.path(cache, "gene_list_1356.csv"))
candidates <- sort(unique(candidate_table$gene))

rna <- data.table::fread(file.path(cache, "per_cohort_betas_long.csv"))
rna <- rna[
  gene_symbol %in% candidates &
    cohort %in% c("TCGA", "CPTAC", "GSE224564", "GSE79668") &
    contrast %in% c("Normal_vs_Early", "Mid_vs_Early", "Late_vs_Early")
]
data.table::setorder(rna, gene_symbol, cohort, contrast, padj)
rna <- rna[!duplicated(rna[, .(gene_symbol, cohort, contrast)])]
rna[, `:=`(
  feature = gene_symbol,
  modality = "rna",
  estimate = ifelse(contrast == "Normal_vs_Early", -beta, beta),
  std_error = se,
  contrast_std = c(
    Normal_vs_Early = "early_vs_normal",
    Mid_vs_Early = "mid_vs_early",
    Late_vs_Early = "late_vs_early"
  )[contrast]
)]
rna_effects <- rna[, .(feature, cohort, modality, contrast = contrast_std, estimate, std_error, padj)]

protein <- data.table::fread(file.path(cache, "combined_protein_lrt.csv"))
protein <- protein[gene %in% candidates]
protein_effects <- data.table::rbindlist(list(
  protein[, .(feature = gene, cohort = "CPTAC_PROT", modality = "protein", contrast = "mid_vs_early", estimate = cpt_mM - cpt_mE)],
  protein[, .(feature = gene, cohort = "CPTAC_PROT", modality = "protein", contrast = "late_vs_early", estimate = cpt_mL - cpt_mE)],
  protein[, .(feature = gene, cohort = "KU_PROT", modality = "protein", contrast = "mid_vs_early", estimate = ku_mM - ku_mE)],
  protein[, .(feature = gene, cohort = "KU_PROT", modality = "protein", contrast = "late_vs_early", estimate = ku_mL - ku_mE)]
), fill = TRUE)
protein_effects[, `:=`(std_error = NA_real_, padj = NA_real_)]
effects <- data.table::rbindlist(list(rna_effects, protein_effects), use.names = TRUE)
effects <- effects[is.finite(estimate)]

holdout_cohorts <- c("PXD046438", "PXD053603", "PXD065581", "PXD066048")
detection <- data.table::rbindlist(lapply(holdout_cohorts, function(cohort_id) {
  matrix_file <- file.path(cache, paste0(cohort_id, "_log2_norm.tsv"))
  matrix_data <- data.table::fread(matrix_file)
  feature_column <- names(matrix_data)[[1]]
  measured <- as.character(matrix_data[[feature_column]])
  values <- as.matrix(matrix_data[, -1, with = FALSE])
  storage.mode(values) <- "double"
  detected_genes <- measured[rowSums(is.finite(values)) > 0]
  data.table::data.table(feature = candidates, cohort = cohort_id, detected = candidates %in% detected_genes)
}))

cell_types <- data.table::fread(file.path(cache, "celltype_mean_1356.csv"))
tumor_stroma <- do.call(pmax, c(cell_types[, .(Ductal, iCAF, myCAF)], list(na.rm = TRUE)))
other_origin <- do.call(pmax, c(cell_types[, .(Acinar, B_cell, Endothelial, Macrophage, Mast, NK_cell, Plasma, T_cell)], list(na.rm = TRUE)))
source <- data.frame(
  feature = cell_types$gene,
  source_score = (tumor_stroma + 1e-6) / (tumor_stroma + other_origin + 2e-6),
  dominant_source = ifelse(tumor_stroma >= other_origin, "tumor_or_stroma", "non_tumor_host"),
  stringsAsFactors = FALSE
)

validation <- data.table::fread(file.path(cache, "two_cohort_signed_stouffer.csv"))
validation <- data.frame(
  feature = validation$gene,
  validation_tested = validation$n_cohorts >= 1,
  validation_replicated = validation$n_cohorts == 2,
  validation_positive = validation$n_cohorts == 2 & validation$meta_padj < 0.10,
  validation_exploratory = validation$n_cohorts >= 1 & validation$meta_padj < 0.10,
  serum_meta_z = validation$meta_z,
  serum_meta_padj = validation$meta_padj,
  stringsAsFactors = FALSE
)

pdac <- list(
  effects = as.data.frame(effects),
  detection = as.data.frame(detection),
  source = source,
  validation = validation,
  universe = candidates,
  split = list(
    rank_input_serum_cohorts = holdout_cohorts,
    validation_serum_cohorts = c("PXD039273", "PXD048034")
  ),
  provenance = list(
    transcriptomics = "PDAC_biomarker_thesis_data_v1: per_cohort_betas_long.csv",
    tissue_proteomics = "PDAC_biomarker_thesis_data_v1: combined_protein_lrt.csv",
    validation = "joint limma per cohort followed by signed Stouffer meta-analysis"
  )
)
saveRDS(pdac, file.path(out_dir, "pdac_stagebridge.rds"))
sb_write_csv(effects, "data", "processed", "pdac_effects.csv")
sb_write_csv(detection, "data", "processed", "pdac_detection_holdout.csv")
sb_write_csv(validation, "data", "processed", "pdac_validation_labels.csv")
message(sprintf(
  "PDAC prepared: %d candidates, %d effect cells, %d independent detection cohorts, %d replicated serum positives.",
  length(candidates), nrow(effects), length(holdout_cohorts), sum(validation$validation_positive, na.rm = TRUE)
))

