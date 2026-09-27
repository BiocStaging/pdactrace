#' Calculate depth-adjusted blood observability
#'
#' A detection in a shallow cohort contributes more information than the same
#' detection in a cohort where nearly every candidate is observed.
#'
#' @param detection Data frame with feature, cohort, and detected columns.
#' @return Data frame with one observability score per feature.
#' @export
sb_observability <- function(detection) {
  if (is.null(detection)) return(NULL)
  if (!is.data.frame(detection)) stop("detection must be a data frame.", call. = FALSE)
  sb_required(detection, c("feature", "cohort", "detected"), "detection")
  x <- as.data.frame(detection, stringsAsFactors = FALSE)
  x$feature <- as.character(x$feature)
  x$cohort <- as.character(x$cohort)
  x$detected <- as.logical(x$detected)
  if (anyNA(x$detected)) stop("detection$detected must be TRUE or FALSE.", call. = FALSE)
  key <- paste(x$feature, x$cohort, sep = "\r")
  if (anyDuplicated(key)) stop("detection contains duplicate feature/cohort rows.", call. = FALSE)

  cohort_rate <- tapply(x$detected, x$cohort, mean)
  information <- -log(pmax(cohort_rate, 0.01))
  breadth_weight <- 1 / sqrt(pmax(cohort_rate, 0.01))
  x$information <- information[x$cohort]
  x$breadth_weight <- breadth_weight[x$cohort]

  rows <- lapply(split(x, x$feature), function(z) {
    breadth <- sum(z$detected * z$breadth_weight) / sum(z$breadth_weight)
    max_information <- sum(z$information)
    info <- if (max_information > 0) sum(z$detected * z$information) / max_information else breadth
    data.frame(
      feature = z$feature[[1]],
      observability = sb_clamp(0.5 * breadth + 0.5 * info),
      detection_breadth = mean(z$detected),
      n_detection_cohorts = nrow(z),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

