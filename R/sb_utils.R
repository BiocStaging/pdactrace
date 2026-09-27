sb_clamp <- function(x, lower = 0, upper = 1) {
  pmin(upper, pmax(lower, x))
}

sb_required <- function(x, columns, object_name) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    stop(
      sprintf("%s is missing required columns: %s", object_name, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
}

sb_robust_scale <- function(x) {
  s <- stats::mad(x, center = stats::median(x, na.rm = TRUE), constant = 1.4826, na.rm = TRUE)
  if (!is.finite(s) || s <= .Machine$double.eps) {
    s <- stats::sd(x, na.rm = TRUE)
  }
  if (!is.finite(s) || s <= .Machine$double.eps) 1 else s
}

sb_robust_z <- function(x) {
  center <- stats::median(x, na.rm = TRUE)
  (x - center) / sb_robust_scale(x)
}

sb_named_score <- function(x, features, value_name) {
  if (is.null(x)) return(rep(NA_real_, length(features)))
  if (is.numeric(x) && !is.null(names(x))) {
    return(sb_clamp(as.numeric(x[features])))
  }
  if (is.data.frame(x)) {
    sb_required(x, c("feature", value_name), deparse(substitute(x)))
    if (anyDuplicated(x$feature)) stop("Feature-level score tables must contain one row per feature.", call. = FALSE)
    return(sb_clamp(x[[value_name]][match(features, x$feature)]))
  }
  stop(sprintf("Expected a named numeric vector or data frame for %s.", value_name), call. = FALSE)
}

sb_geometric_score <- function(components, weights) {
  keep <- is.finite(components) & names(components) %in% names(weights)
  if (!any(keep)) return(NA_real_)
  w <- weights[names(components)[keep]]
  w <- w / sum(w)
  exp(sum(w * log(0.02 + 0.98 * sb_clamp(components[keep]))))
}

