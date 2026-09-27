sb_random_effects <- function(estimate, std_error) {
  keep <- is.finite(estimate) & is.finite(std_error) & std_error > 0
  estimate <- estimate[keep]
  std_error <- std_error[keep]
  k <- length(estimate)
  if (!k) return(c(estimate = NA_real_, std_error = NA_real_, i2 = NA_real_, tau2 = NA_real_))
  if (k == 1L) return(c(estimate = estimate, std_error = std_error, i2 = NA_real_, tau2 = 0))

  w <- 1 / std_error^2
  fixed <- sum(w * estimate) / sum(w)
  q <- sum(w * (estimate - fixed)^2)
  c_term <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (q - (k - 1)) / c_term)
  wr <- 1 / (std_error^2 + tau2)
  pooled <- sum(wr * estimate) / sum(wr)
  pooled_se <- sqrt(1 / sum(wr))
  i2 <- if (q > 0) max(0, (q - (k - 1)) / q) * 100 else 0
  c(estimate = pooled, std_error = pooled_se, i2 = i2, tau2 = tau2)
}

sb_layer_summary <- function(x) {
  meta <- sb_random_effects(x$estimate, x$std_error)
  precision <- x$confidence + 0.05
  direction <- abs(stats::weighted.mean(sign(x$estimate), precision, na.rm = TRUE))
  replication <- if (nrow(x) == 1L) {
    0.5
  } else if (is.finite(meta[["i2"]])) {
    0.5 * direction + 0.5 * (1 - meta[["i2"]] / 100)
  } else {
    direction
  }
  pooled_effect <- if (is.finite(meta[["estimate"]])) {
    stats::weighted.mean(x$robust_effect, 1 / x$std_error^2, na.rm = TRUE)
  } else {
    stats::median(x$robust_effect, na.rm = TRUE)
  }
  data.frame(
    feature = x$feature[[1]],
    modality = x$modality[[1]],
    contrast = x$contrast[[1]],
    pooled_effect = pooled_effect,
    evidence = 1 - exp(-mean(abs(x$signal_z), na.rm = TRUE) / 2),
    replication = sb_clamp(replication),
    i2 = unname(meta[["i2"]]),
    n_cohorts = nrow(x),
    stringsAsFactors = FALSE
  )
}

sb_cross_modal_concordance <- function(x) {
  modalities <- unique(x$modality)
  if (length(modalities) < 2L) return(NA_real_)
  pairs <- utils::combn(modalities, 2, simplify = FALSE)
  values <- vapply(pairs, function(pair) {
    a <- x[x$modality == pair[[1]], c("contrast", "pooled_effect")]
    b <- x[x$modality == pair[[2]], c("contrast", "pooled_effect")]
    joined <- merge(a, b, by = "contrast", suffixes = c("_a", "_b"))
    joined <- joined[is.finite(joined$pooled_effect_a) & is.finite(joined$pooled_effect_b), ]
    if (!nrow(joined)) return(NA_real_)
    if (nrow(joined) == 1L) {
      return(if (sign(joined$pooled_effect_a) == sign(joined$pooled_effect_b)) 1 else 0)
    }
    denom <- sqrt(sum(joined$pooled_effect_a^2) * sum(joined$pooled_effect_b^2))
    if (!is.finite(denom) || denom == 0) return(0.5)
    sb_clamp((sum(joined$pooled_effect_a * joined$pooled_effect_b) / denom + 1) / 2)
  }, numeric(1))
  if (all(!is.finite(values))) NA_real_ else mean(values, na.rm = TRUE)
}

