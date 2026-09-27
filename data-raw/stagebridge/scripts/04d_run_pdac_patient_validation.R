#!/usr/bin/env Rscript
source(file.path(dirname(grep("^--file=", commandArgs(FALSE), value = TRUE) |> sub("^--file=", "", x = _)), "common.R"))

required <- c("data.table", "limma")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))

library(data.table)

cache_dir <- sb_path("data-raw", "cache", "pdac")
pdac <- readRDS(sb_path("data", "processed", "pdac_stagebridge.rds"))
method_scores <- read.csv(sb_path("results", "tables", "pdac_comparator_scores.csv"), check.names = FALSE)
score_methods <- setdiff(names(method_scores), "feature")

direction_effects <- as.data.table(pdac$effects)[
  modality == "rna" & contrast == "early_vs_normal" & is.finite(std_error) & std_error > 0
]
frozen_direction <- direction_effects[, {
  z <- estimate / std_error
  signed_stouffer <- sum(z) / sqrt(.N)
  list(
    frozen_direction = if (is.finite(signed_stouffer) && signed_stouffer != 0) sign(signed_stouffer) else NA_real_,
    direction_coherence = abs(mean(sign(estimate))),
    direction_cohorts = .N,
    direction_z = signed_stouffer
  )
}, by = feature]
method_scores <- merge(method_scores, frozen_direction, by = "feature", all.x = TRUE, sort = FALSE)

read_matrix <- function(filename) {
  x <- fread(file.path(cache_dir, filename))
  feature <- trimws(as.character(x[[1L]]))
  matrix <- as.matrix(x[, -1L, with = FALSE])
  storage.mode(matrix) <- "double"
  rownames(matrix) <- feature
  keep <- nzchar(feature) & !duplicated(feature)
  matrix[keep, , drop = FALSE]
}

matrix39 <- read_matrix("PXD039273_log2_norm.tsv")
metadata39 <- fread(file.path(cache_dir, "PXD039273_clinical_stage.csv"))
metadata39 <- metadata39[match(colnames(matrix39), sample_id)]
if (anyNA(metadata39$group)) stop("PXD039273 metadata did not match every matrix column.")
metadata39[, condition := fifelse(
  group == "PDAC", "PDAC", fifelse(group == "Pancreatitis", "Pancreatitis", "Control")
)]

matrix48 <- read_matrix("PXD048034_log2_norm.tsv")
metadata48 <- data.table(
  sample_id = colnames(matrix48),
  condition = ifelse(grepl("_Cancer_", colnames(matrix48)), "PDAC", "Control")
)
if (!identical(c(sum(metadata48$condition == "PDAC"), sum(metadata48$condition == "Control")), c(40L, 40L))) {
  stop("Expected 40 PDAC and 40 controls in PXD048034.")
}

cohorts <- list(
  PXD039273 = list(expression = matrix39, metadata = metadata39[, .(sample_id, condition)]),
  PXD048034 = list(expression = matrix48, metadata = metadata48)
)

robust_scale <- function(x) {
  observed <- x[is.finite(x)]
  if (!length(observed)) return(rep(0, length(x)))
  spread <- stats::sd(observed)
  if (!is.finite(spread) || spread <= 0) spread <- 1
  x[!is.finite(x)] <- min(observed) - 0.5 * spread
  center <- stats::median(x)
  scale <- stats::mad(x, center = center, constant = 1)
  if (!is.finite(scale) || scale <= 0) scale <- stats::sd(x)
  if (!is.finite(scale) || scale <= 0) return(rep(0, length(x)))
  pmax(-8, pmin(8, (x - center) / scale))
}

# Freeze cohort-specific panels from discovery scores, tissue direction and
# label-free assay completeness before evaluating any patient labels.
panel_rows <- list()
signature_rows <- list()
for (cohort_name in names(cohorts)) {
  cohort <- cohorts[[cohort_name]]
  expression <- cohort$expression
  assay_eligible <- rowMeans(is.finite(expression)) >= 0.70
  scaled <- t(apply(expression[assay_eligible, , drop = FALSE], 1L, robust_scale))
  rownames(scaled) <- rownames(expression)[assay_eligible]
  colnames(scaled) <- colnames(expression)
  eligible <- method_scores[
    method_scores$feature %in% rownames(scaled) & is.finite(method_scores$frozen_direction),
  ]
  cohorts[[cohort_name]]$scaled_expression <- scaled
  for (method_name in score_methods) {
    ordering <- order(-eligible[[method_name]], eligible$feature, na.last = NA)
    ordered <- eligible[ordering, ]
    for (k in c(5L, 10L, 20L)) {
      selected <- head(ordered, min(k, nrow(ordered)))
      rank_weight <- 1 / log2(seq_len(nrow(selected)) + 1)
      for (target_relation in c("concordant", "inverse")) {
        relation_sign <- if (target_relation == "inverse") -1 else 1
        effective_direction <- selected$frozen_direction * relation_sign
        panel_rows[[length(panel_rows) + 1L]] <- data.frame(
          cohort = cohort_name,
          method = method_name,
          panel_size = k,
          target_relation = target_relation,
          panel_order = seq_len(nrow(selected)),
          feature = selected$feature,
          frozen_tissue_direction = selected$frozen_direction,
          effective_blood_direction = effective_direction,
          direction_coherence = selected$direction_coherence,
          rank_weight = rank_weight,
          discovery_score = selected[[method_name]],
          stringsAsFactors = FALSE
        )
        oriented <- scaled[selected$feature, , drop = FALSE] *
          (effective_direction * rank_weight)
        signature <- colSums(oriented) / sum(rank_weight)
        signature_rows[[length(signature_rows) + 1L]] <- data.frame(
          cohort = cohort_name,
          sample_id = names(signature),
          method = method_name,
          panel_size = k,
          target_relation = target_relation,
          signature_score = unname(signature),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
panels <- rbindlist(panel_rows)
signature_scores <- rbindlist(signature_rows)
sb_write_csv(panels, "results", "tables", "pdac_patient_frozen_panels.csv")

# Patient labels enter only below this line.
metadata_all <- rbindlist(lapply(names(cohorts), function(cohort_name) {
  out <- copy(cohorts[[cohort_name]]$metadata)
  out[, cohort := cohort_name]
  out
}), fill = TRUE)
signature_scores <- merge(
  signature_scores,
  metadata_all[, .(cohort, sample_id, condition)],
  by = c("cohort", "sample_id"),
  all.x = TRUE,
  sort = FALSE
)

bootstrap_metrics <- function(truth, score, n_boot = 1000L, seed = 20270719L) {
  positive <- which(truth)
  negative <- which(!truth)
  set.seed(seed)
  values <- replicate(n_boot, {
    index <- c(sample(positive, length(positive), replace = TRUE),
               sample(negative, length(negative), replace = TRUE))
    c(
      auroc = sb_binary_auc(truth[index], score[index]),
      auprc = sb_average_precision(truth[index], score[index])
    )
  })
  data.frame(
    auroc = sb_binary_auc(truth, score),
    auroc_low = unname(stats::quantile(values["auroc", ], 0.025)),
    auroc_high = unname(stats::quantile(values["auroc", ], 0.975)),
    auprc = sb_average_precision(truth, score),
    auprc_low = unname(stats::quantile(values["auprc", ], 0.025)),
    auprc_high = unname(stats::quantile(values["auprc", ], 0.975))
  )
}

patient_endpoints <- list(
  PXD039273_control = list(cohort = "PXD039273", negative = "Control"),
  PXD039273_pancreatitis = list(cohort = "PXD039273", negative = "Pancreatitis"),
  PXD048034_control = list(cohort = "PXD048034", negative = "Control")
)
patient_metric_rows <- list()
for (endpoint_name in names(patient_endpoints)) {
  endpoint <- patient_endpoints[[endpoint_name]]
  x <- signature_scores[
    signature_scores$cohort == endpoint$cohort &
      signature_scores$condition %in% c("PDAC", endpoint$negative)
  ]
  grouped <- split(x, by = c("method", "panel_size", "target_relation"))
  patient_metric_rows[[endpoint_name]] <- rbindlist(lapply(grouped, function(z) {
    truth <- z$condition == "PDAC"
    metrics <- bootstrap_metrics(truth, z$signature_score)
    data.frame(
      endpoint = endpoint_name,
      cohort = endpoint$cohort,
      reference = endpoint$negative,
      method = z$method[[1L]],
      panel_size = z$panel_size[[1L]],
      target_relation = z$target_relation[[1L]],
      n_case = sum(truth),
      n_reference = sum(!truth),
      metrics,
      stringsAsFactors = FALSE
    )
  }))
}
patient_metrics <- rbindlist(patient_metric_rows)
setorder(patient_metrics, endpoint, panel_size, target_relation, -auroc)
sb_write_csv(signature_scores, "results", "tables", "pdac_patient_scores.csv")
sb_write_csv(patient_metrics, "results", "tables", "pdac_patient_metrics.csv")

fit_cohort_differential <- function(cohort_name) {
  cohort <- cohorts[[cohort_name]]
  keep_samples <- cohort$metadata$condition %in% c("Control", "PDAC")
  metadata <- cohort$metadata[keep_samples]
  expression <- cohort$expression[, metadata$sample_id, drop = FALSE]
  keep_features <- rowSums(is.finite(expression[, metadata$condition == "Control", drop = FALSE])) >= 10L &
    rowSums(is.finite(expression[, metadata$condition == "PDAC", drop = FALSE])) >= 10L
  expression <- expression[keep_features, , drop = FALSE]
  condition <- factor(metadata$condition, levels = c("Control", "PDAC"))
  design <- stats::model.matrix(~ 0 + condition)
  colnames(design) <- levels(condition)
  contrast <- limma::makeContrasts(PDAC_vs_Control = PDAC - Control, levels = design)
  fit <- suppressWarnings(limma::lmFit(expression, design))
  fit <- suppressWarnings(limma::eBayes(limma::contrasts.fit(fit, contrast), robust = TRUE))
  out <- limma::topTable(fit, coef = 1L, number = Inf, sort.by = "none")
  out$feature <- rownames(out)
  rownames(out) <- NULL
  out$cohort <- cohort_name
  out$external_differential_fdr10 <- out$adj.P.Val < 0.10
  out
}
differential <- rbindlist(lapply(names(cohorts), fit_cohort_differential))
setorder(differential, cohort, P.Value)
sb_write_csv(differential, "results", "tables", "pdac_patient_differential.csv")

feature_metric_rows <- list()
overlap_rows <- list()
for (cohort_name in names(cohorts)) {
  de <- differential[differential$cohort == cohort_name]
  overlap <- merge(de, method_scores, by = "feature", all = FALSE, sort = FALSE)
  truth <- overlap$external_differential_fdr10
  metrics <- sb_evaluate_methods(
    truth,
    as.data.frame(overlap)[, c("feature", score_methods), drop = FALSE],
    top_k = min(20L, nrow(overlap))
  )
  metrics$endpoint <- paste0(cohort_name, " FDR < 0.10")
  metrics$n_assayed <- nrow(de)
  metrics$n_overlap <- nrow(overlap)
  metrics$n_positive <- sum(truth)
  metrics$prevalence <- mean(truth)
  feature_metric_rows[[cohort_name]] <- metrics
  overlap_rows[[cohort_name]] <- overlap
}

de39 <- differential[differential$cohort == "PXD039273"]
de48 <- differential[differential$cohort == "PXD048034"]
replicated <- merge(de39, de48, by = "feature", suffixes = c("_39", "_48"))
z39 <- stats::qnorm(pmax(replicated$P.Value_39, 1e-300) / 2, lower.tail = FALSE) * sign(replicated$logFC_39)
z48 <- stats::qnorm(pmax(replicated$P.Value_48, 1e-300) / 2, lower.tail = FALSE) * sign(replicated$logFC_48)
replicated$meta_z <- (z39 + z48) / sqrt(2)
replicated$meta_p <- 2 * stats::pnorm(-abs(replicated$meta_z))
replicated$meta_fdr <- p.adjust(replicated$meta_p, method = "BH")
replicated$external_differential_fdr10 <- replicated$meta_fdr < 0.10 &
  sign(replicated$logFC_39) == sign(replicated$logFC_48)
replicated_scores <- merge(replicated, method_scores, by = "feature", all = FALSE, sort = FALSE)
replicated_truth <- replicated_scores$external_differential_fdr10
replicated_metrics <- sb_evaluate_methods(
  replicated_truth,
  as.data.frame(replicated_scores)[, c("feature", score_methods), drop = FALSE],
  top_k = min(20L, nrow(replicated_scores))
)
replicated_metrics$endpoint <- "two-cohort signed Stouffer FDR < 0.10 with direction agreement"
replicated_metrics$n_assayed <- nrow(replicated)
replicated_metrics$n_overlap <- nrow(replicated_scores)
replicated_metrics$n_positive <- sum(replicated_truth)
replicated_metrics$prevalence <- mean(replicated_truth)
feature_metrics <- rbindlist(c(feature_metric_rows, list(replicated = replicated_metrics)), fill = TRUE)
sb_write_csv(rbindlist(overlap_rows, fill = TRUE), "results", "tables", "pdac_patient_overlap.csv")
sb_write_csv(replicated_scores, "results", "tables", "pdac_patient_replicated.csv")
sb_write_csv(feature_metrics, "results", "tables", "pdac_patient_feature_metrics.csv")

saveRDS(
  list(
    cohorts = cohorts,
    panels = as.data.frame(panels),
    signature_scores = as.data.frame(signature_scores),
    patient_metrics = as.data.frame(patient_metrics),
    differential = as.data.frame(differential),
    replicated = replicated_scores,
    feature_metrics = as.data.frame(feature_metrics)
  ),
  sb_path("results", "pdac_patient_validation.rds")
)

primary_rows <- patient_metrics[
  method == "StageBridge" & panel_size == 10L & target_relation == "concordant"
]
message("PDAC patient-level relation audit complete:")
for (i in seq_len(nrow(primary_rows))) {
  message(sprintf(
    "  %s AUROC %.3f (95%% CI %.3f-%.3f)",
    primary_rows$endpoint[[i]], primary_rows$auroc[[i]],
    primary_rows$auroc_low[[i]], primary_rows$auroc_high[[i]]
  ))
}
