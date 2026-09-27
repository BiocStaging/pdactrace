# General, dataset-agnostic prospective one-shot runner.
#
# This factors the gate stack + statistical endpoint + PASS/FAIL/INCONCLUSIVE
# adjudication that the per-accession audit scripts (e.g. scripts/38_open_pxd055218_crc_once.R)
# each implement by hand into a single exported function, so the "protocol" is an
# executable, testable runner rather than a discipline enacted by bespoke scripts.
# `sb_one_shot()` intentionally never mutates a locked record; it is the pure decision
# core that a locked wrapper calls once after outcome access.

sb_panel_score <- function(abundance, panel, selected) {
  wide <- reshape(
    abundance[abundance$gene %in% selected, c("gene", "sample", "log2_abundance")],
    idvar = "gene", timevar = "sample", direction = "wide"
  )
  rownames(wide) <- wide$gene
  mat <- as.matrix(wide[selected, setdiff(names(wide), "gene"), drop = FALSE])
  colnames(mat) <- sub("^log2_abundance\\.", "", colnames(mat))
  storage.mode(mat) <- "double"
  for (i in seq_len(nrow(mat))) {
    obs <- mat[i, is.finite(mat[i, ])]
    if (!length(obs)) stop("Selected panel feature has no finite abundance: ", selected[[i]])
    mat[i, !is.finite(mat[i, ])] <- stats::median(obs)
  }
  row_sd <- apply(mat, 1L, stats::sd)
  if (any(!is.finite(row_sd) | row_sd <= 0)) stop("Selected panel contains a zero-variance feature.")
  z <- sweep(sweep(mat, 1L, rowMeans(mat), "-"), 1L, row_sd, "/")
  sel <- panel[match(selected, panel$feature), ]
  w <- sel$frozen_direction * sel$rank_weight
  w <- w / sum(abs(w))
  colSums(z * w)
}

sb_permutation_p <- function(truth, score, replicates, seed) {
  observed <- sb_binary_auc(truth, score)
  set.seed(seed)
  null <- replicate(replicates, sb_binary_auc(sample(truth, replace = FALSE), score))
  (1 + sum(null >= observed, na.rm = TRUE)) / (replicates + 1)
}

sb_auc_interval <- function(truth, score, replicates, seed) {
  truth <- as.logical(truth)
  pos <- which(truth); neg <- which(!truth)
  set.seed(seed)
  vals <- replicate(replicates, {
    idx <- c(sample(pos, length(pos), replace = TRUE), sample(neg, length(neg), replace = TRUE))
    sb_binary_auc(truth[idx], score[idx])
  })
  stats::quantile(vals, c(0.025, 0.975), na.rm = TRUE, names = FALSE, type = 8)
}

sb_feature_replication <- function(abundance, frozen, selected, groups, min_detect = 3L) {
  eligible <- abundance[abundance$gene %in% frozen$feature, ]
  det <- table(eligible$gene, eligible$group)
  testable <- rownames(det)[
    det[, groups$case] >= min_detect &
      det[, groups$primary_reference] >= min_detect
  ]
  eligible <- eligible[eligible$gene %in% testable, ]
  if (!nrow(eligible)) return(data.frame())
  wide <- reshape(
    eligible[, c("gene", "sample", "log2_abundance")],
    idvar = "gene", timevar = "sample", direction = "wide"
  )
  rownames(wide) <- wide$gene
  mat <- as.matrix(wide[, setdiff(names(wide), "gene"), drop = FALSE])
  colnames(mat) <- sub("^log2_abundance\\.", "", colnames(mat))
  storage.mode(mat) <- "double"
  smeta <- unique(abundance[, c("sample", "group")])
  grp <- factor(smeta$group[match(colnames(mat), smeta$sample)])
  keep_two <- grp %in% c(groups$case, groups$primary_reference)
  mat <- mat[, keep_two, drop = FALSE]
  grp <- droplevels(grp[keep_two])
  grp <- stats::relevel(grp, ref = groups$primary_reference)
  # simple imputation for limma stability
  for (i in seq_len(nrow(mat))) {
    obs <- mat[i, is.finite(mat[i, ])]
    if (length(obs)) mat[i, !is.finite(mat[i, ])] <- stats::median(obs)
  }
  keep_row <- apply(mat, 1L, function(r) all(is.finite(r)) && stats::sd(r) > 0)
  mat <- mat[keep_row, , drop = FALSE]
  if (!nrow(mat)) return(data.frame())
  design <- stats::model.matrix(~grp)
  fit <- limma::eBayes(limma::lmFit(mat, design), robust = TRUE)
  coef_col <- ncol(design)
  data.frame(
    feature = rownames(mat),
    estimate = fit$coefficients[, coef_col],
    p_value = fit$p.value[, coef_col],
    frozen_direction = frozen$frozen_direction[match(rownames(mat), frozen$feature)],
    stringsAsFactors = FALSE
  )
}

#' Run one prospectively locked destination-cohort audit to a terminal decision
#'
#' Applies the frozen gate stack (ontology/group-count, broad coverage, top-10 and
#' top-20 panel coverage, complete panel) and, only if every gate passes, the locked
#' statistical endpoint (rank-weighted panel AUROC with a permutation test and a
#' feature-replication count). Returns a terminal decision whose `status` is exactly
#' one of `PASS`, `FAIL` or `INCONCLUSIVE`. `INCONCLUSIVE` ("the test never ran")
#' is returned whenever a gate halts execution and is operationally distinct from a
#' biological `FAIL` (the endpoint ran and did not meet the criteria).
#'
#' @param abundance Long data frame: `gene`, `sample`, `group`, `log2_abundance`.
#' @param frozen Broad frozen score table with `feature` and `frozen_direction`.
#' @param panel Ordered top-k panel with `feature`, `frozen_direction`, `rank_weight`.
#' @param groups List with `case`, `primary_reference`, `secondary_reference` group labels.
#' @param expected_counts Named integer vector of expected samples per group.
#' @param rule Named list of frozen thresholds (see [sb_default_one_shot_rule()]).
#' @return A terminal decision list with `status`, `gates`, `statistical_runner_entered`
#'   and, when the runner ran, the endpoint statistics.
#' @export
sb_one_shot <- function(abundance, frozen, panel, groups, expected_counts,
                        rule = sb_default_one_shot_rule()) {
  gate_result <- function(status, gates, extra = list()) {
    c(list(status = status, statistical_runner_entered = FALSE,
           failed_gates = names(gates)[!unlist(gates)], gates = gates), extra)
  }
  tryCatch({
    smeta <- unique(abundance[, c("sample", "group")])
    obs <- table(factor(smeta$group, levels = names(expected_counts)))
    n_expected <- sum(expected_counts)
    ontology <- nrow(smeta) == n_expected &&
      !anyNA(smeta$group) &&
      all(smeta$group %in% names(expected_counts)) &&
      identical(as.integer(obs[names(expected_counts)]), as.integer(expected_counts))

    covered <- intersect(frozen$feature, unique(abundance$gene))
    top10 <- sum(panel$feature[seq_len(min(10L, nrow(panel)))] %in% covered)
    top20 <- sum(panel$feature %in% covered)

    det <- table(abundance$gene[abundance$gene %in% panel$feature],
                 abundance$group[abundance$gene %in% panel$feature])
    per_group_ok <- function(g) if (g %in% rownames(det)) det[g, ] else stats::setNames(rep(0L, ncol(det)), colnames(det))
    eligible_panel <- vapply(panel$feature, function(g) {
      row <- if (g %in% rownames(det)) det[g, ] else return(FALSE)
      all(row[names(expected_counts)] >= rule$minimum_per_group_detections_for_panel_gene)
    }, logical(1))
    selected <- utils::head(panel$feature[eligible_panel], rule$panel_size)

    gates <- list(
      ontology = ontology,
      broad_coverage = length(covered) >= rule$minimum_broad_candidate_coverage,
      top10_coverage = top10 >= rule$minimum_top10_coverage,
      top20_coverage = top20 >= rule$minimum_top20_coverage,
      complete_panel = length(selected) == rule$panel_size
    )
    audit <- list(observed_counts = as.list(as.integer(obs[names(expected_counts)])),
                  broad_coverage = length(covered), top10_coverage = top10,
                  top20_coverage = top20, selected_panel = selected)

    if (!all(unlist(gates))) {
      return(gate_result("INCONCLUSIVE", gates, c(audit,
        list(interpretation = "The locked statistical endpoint was not entered; this is not a performance result."))))
    }

    score <- sb_panel_score(abundance, panel, selected)
    grp <- smeta$group[match(names(score), smeta$sample)]
    prim <- grp %in% c(groups$case, groups$primary_reference)
    truth <- grp[prim] == groups$case
    s <- score[prim]
    auc <- sb_binary_auc(truth, s)
    pperm <- sb_permutation_p(truth, s, rule$permutation_replicates, rule$permutation_seed)
    ci <- sb_auc_interval(truth, s, rule$bootstrap_replicates, rule$bootstrap_seed)

    repl <- sb_feature_replication(abundance, frozen, selected, groups)
    repl_sel <- repl[repl$feature %in% selected, , drop = FALSE]
    replications <- if (nrow(repl_sel)) sum(
      sign(repl_sel$estimate) == repl_sel$frozen_direction &
        repl_sel$p_value < rule$feature_replication_nominal_p, na.rm = TRUE) else 0L

    status <- if (auc >= rule$minimum_pass_auc &&
                  pperm <= rule$maximum_pass_permutation_p &&
                  replications >= rule$minimum_direction_consistent_nominal_feature_replications)
      "PASS" else "FAIL"

    c(list(status = status, statistical_runner_entered = TRUE, gates = gates,
           primary = list(endpoint = paste0(groups$case, "_vs_", groups$primary_reference),
                          auc = auc, bootstrap_95_ci = unname(ci), permutation_p = pperm,
                          n_case = sum(truth), n_reference = sum(!truth)),
           direction_consistent_nominal_replications = replications,
           interpretation = if (status == "PASS")
             "The locked endpoint is evaluable and meets the frozen transfer criteria."
           else "The locked endpoint is evaluable but does not meet the frozen transfer criteria."),
      audit)
  }, error = function(e) {
    list(status = "INCONCLUSIVE", statistical_runner_entered = FALSE,
         error = conditionMessage(e),
         interpretation = "A locked technical or schema failure prevented the endpoint; this is not a performance result.")
  })
}

#' Default frozen thresholds for [sb_one_shot()]
#'
#' Mirrors the PXD055218 CRC decision rule (`prospective/locks/PXD055218_crc_v1/DECISION_RULE.json`).
#' @export
sb_default_one_shot_rule <- function() {
  list(
    minimum_broad_candidate_coverage = 200L,
    minimum_top10_coverage = 5L,
    minimum_top20_coverage = 10L,
    minimum_per_group_detections_for_panel_gene = 14L,
    panel_size = 10L,
    permutation_replicates = 2000L,
    permutation_seed = 20260718L,
    bootstrap_replicates = 1000L,
    bootstrap_seed = 20260719L,
    minimum_pass_auc = 0.70,
    maximum_pass_permutation_p = 0.05,
    minimum_direction_consistent_nominal_feature_replications = 2L,
    feature_replication_nominal_p = 0.05
  )
}
