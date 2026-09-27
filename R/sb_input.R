#' Validate a StageBridge effect table
#'
#' @param effects Data frame with one row per feature, cohort, modality, and contrast.
#' @return A normalized data frame with numeric estimates and standard errors.
#' @export
sb_validate_effects <- function(effects) {
  if (!is.data.frame(effects)) stop("effects must be a data frame.", call. = FALSE)
  required <- c("feature", "cohort", "modality", "contrast", "estimate")
  sb_required(effects, required, "effects")

  out <- as.data.frame(effects, stringsAsFactors = FALSE)
  for (column in c("feature", "cohort", "modality", "contrast")) {
    out[[column]] <- trimws(as.character(out[[column]]))
    if (any(!nzchar(out[[column]]))) stop(sprintf("effects$%s contains empty values.", column), call. = FALSE)
  }
  out$estimate <- as.numeric(out$estimate)
  if (!"std_error" %in% names(out)) out$std_error <- NA_real_
  out$std_error <- as.numeric(out$std_error)
  out$std_error[!is.finite(out$std_error) | out$std_error <= 0] <- NA_real_
  out <- out[is.finite(out$estimate), , drop = FALSE]
  if (!nrow(out)) stop("effects contains no finite estimates.", call. = FALSE)

  key <- do.call(paste, c(out[c("feature", "cohort", "modality", "contrast")], sep = "\r"))
  if (anyDuplicated(key)) {
    stop("effects contains duplicate feature/cohort/modality/contrast rows.", call. = FALSE)
  }
  rownames(out) <- NULL
  out
}

sb_prepare_effects <- function(effects) {
  out <- sb_validate_effects(effects)
  layer <- interaction(out$cohort, out$modality, out$contrast, drop = TRUE, lex.order = TRUE)
  out$robust_effect <- stats::ave(out$estimate, layer, FUN = sb_robust_z)
  raw_z <- out$estimate / out$std_error
  out$signal_z <- ifelse(is.finite(raw_z), raw_z, out$robust_effect)
  out$signal_z <- sb_clamp(out$signal_z, -8, 8)
  out$confidence <- sb_clamp(abs(out$signal_z) / 3)
  out
}
