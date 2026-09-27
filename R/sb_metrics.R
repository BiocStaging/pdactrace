#' Binary area under the ROC curve
#' @param truth Logical or binary vector.
#' @param score Numeric score; larger values predict positives.
#' @export
sb_binary_auc <- function(truth, score) {
  keep <- !is.na(truth) & is.finite(score)
  truth <- as.logical(truth[keep])
  score <- score[keep]
  n_pos <- sum(truth)
  n_neg <- sum(!truth)
  if (!n_pos || !n_neg) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[truth]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

#' Average precision
#' @param truth Logical or binary vector.
#' @param score Numeric score; larger values predict positives.
#' @export
sb_average_precision <- function(truth, score) {
  keep <- !is.na(truth) & is.finite(score)
  truth <- as.logical(truth[keep])
  score <- score[keep]
  if (!sum(truth)) return(NA_real_)
  ordering <- order(score, decreasing = TRUE)
  y <- truth[ordering]
  ordered_score <- score[ordering]
  score_group <- cumsum(c(TRUE, diff(ordered_score) != 0))
  positives_at_threshold <- as.numeric(rowsum(as.integer(y), score_group, reorder = FALSE))
  observations_at_threshold <- tabulate(score_group)
  precision <- cumsum(positives_at_threshold) / cumsum(observations_at_threshold)
  sum(precision * positives_at_threshold) / sum(y)
}

#' Recall among the top-ranked candidates
#' @param truth Logical or binary vector.
#' @param score Numeric score.
#' @param k Number of top candidates.
#' @export
sb_topk_recall <- function(truth, score, k = 50L) {
  keep <- !is.na(truth) & is.finite(score)
  truth <- as.logical(truth[keep])
  score <- score[keep]
  if (!sum(truth)) return(NA_real_)
  selected <- order(score, decreasing = TRUE)[seq_len(min(k, length(score)))]
  sum(truth[selected]) / sum(truth)
}

sb_evaluate_methods <- function(truth, score_table, top_k = 50L) {
  methods <- setdiff(names(score_table), "feature")
  rows <- lapply(methods, function(method) {
    score <- score_table[[method]]
    keep <- !is.na(truth) & is.finite(score)
    evaluated_truth <- as.logical(truth[keep])
    data.frame(
      method = method,
      auroc = sb_binary_auc(truth, score),
      auprc = sb_average_precision(truth, score),
      top_k_recall = sb_topk_recall(truth, score, top_k),
      n_evaluated = sum(keep),
      n_positive_evaluated = sum(evaluated_truth),
      prevalence_evaluated = if (any(keep)) mean(evaluated_truth) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

sb_common_method_rows <- function(truth, score_table) {
  methods <- setdiff(names(score_table), "feature")
  if (!length(methods)) stop("A score table must contain at least one method.", call. = FALSE)
  score_matrix <- as.matrix(score_table[methods])
  !is.na(truth) & rowSums(is.finite(score_matrix)) == length(methods)
}

sb_evaluate_common_methods <- function(truth, score_table, top_k = 50L) {
  keep <- sb_common_method_rows(truth, score_table)
  if (!any(keep)) stop("No rows have finite scores for every method.", call. = FALSE)
  out <- sb_evaluate_methods(
    truth[keep],
    score_table[keep, , drop = FALSE],
    top_k = min(top_k, sum(keep))
  )
  out$evaluation_scope <- "strict_all_methods_common"
  out$n_source_universe <- nrow(score_table)
  out$n_common <- sum(keep)
  out$n_positive_common <- sum(as.logical(truth[keep]))
  out$prevalence_common <- mean(as.logical(truth[keep]))
  out
}
