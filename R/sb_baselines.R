#' Score candidates with StageBridge and benchmark comparators
#'
#' @param effects Stage-resolved effect table.
#' @param detection Optional blood detection table.
#' @param source_score Optional source score used by StageBridge only.
#' @return Wide feature-by-method score table.
#' @examples
#' sim <- sb_simulate(n_features = 20L, seed = 1L)
#' head(sb_method_scores(sim$effects, sim$detection), 3)
#' @export
sb_method_scores <- function(effects, detection = NULL, source_score = NULL) {
  x <- sb_prepare_effects(effects)
  x$layer <- interaction(x$cohort, x$modality, x$contrast, drop = TRUE, lex.order = TRUE)
  features <- sort(unique(x$feature))
  layers <- sort(unique(x$layer))
  matrix_z <- matrix(NA_real_, nrow = length(features), ncol = length(layers), dimnames = list(features, layers))
  matrix_z[cbind(match(x$feature, features), match(x$layer, layers))] <- x$signal_z

  mean_abs <- rowMeans(abs(matrix_z), na.rm = TRUE)
  vote <- rowMeans(abs(matrix_z) >= 1.96, na.rm = TRUE)
  stouffer <- apply(matrix_z, 1, function(z) abs(sum(z, na.rm = TRUE)) / sqrt(sum(is.finite(z))))
  intersection <- apply(matrix_z, 1, function(z) {
    z <- z[is.finite(z)]
    if (!length(z)) return(NA_real_)
    mean(abs(z) >= 1.96) * stats::median(abs(z))
  })

  percentile <- apply(abs(matrix_z), 2, function(z) {
    keep <- is.finite(z)
    out <- rep(NA_real_, length(z))
    out[keep] <- rank(-z[keep], ties.method = "average") / sum(keep)
    out
  })
  if (is.null(dim(percentile))) percentile <- matrix(percentile, ncol = 1)
  rank_product <- apply(percentile, 1, function(z) {
    z <- z[is.finite(z)]
    if (!length(z)) return(NA_real_)
    1 - exp(mean(log(pmax(z, 1e-8))))
  })

  svd_matrix <- matrix_z
  svd_matrix[!is.finite(svd_matrix)] <- 0
  sv <- tryCatch(base::svd(svd_matrix, nu = 1, nv = 1), error = function(e) NULL)
  svd_loading <- if (is.null(sv)) rep(NA_real_, length(features)) else abs(sv$u[, 1] * sv$d[[1]])

  bridge <- stagebridge(x, detection = detection, source_score = source_score, compute_pareto = FALSE)
  bridge_score <- bridge$stagebridge_score[match(features, bridge$feature)]
  obs <- sb_observability(detection)
  obs_score <- if (is.null(obs)) rep(1, length(features)) else obs$observability[match(features, obs$feature)]
  observed_rank_product <- rank_product * (0.25 + 0.75 * ifelse(is.finite(obs_score), obs_score, 0))

  out <- data.frame(
    feature = features,
    mean_abs_z = mean_abs,
    vote_count = vote,
    signed_stouffer = stouffer,
    intersection = intersection,
    rank_product = rank_product,
    svd_loading = svd_loading,
    observed_rank_product = observed_rank_product,
    StageBridge = bridge_score,
    stringsAsFactors = FALSE
  )
  attr(out, "stagebridge_components") <- bridge
  out
}

#' Score shared features with matched-sample multiblock PLS-DA
#'
#' Fits the dense multiblock PLS-DA model implemented by mixOmics and combines
#' block-specific loading magnitudes into a shared-feature priority score. This
#' comparator requires the same observations, in the same order, in every
#' block; it is therefore not applicable to unpaired cohort-level effect tables.
#'
#' @param blocks Named list of sample-by-feature numeric matrices.
#' @param outcome Factor or vector of class labels, one per sample.
#' @param ncomp Number of latent components.
#' @param design Off-diagonal connection strength between omics blocks.
#' @return A feature-by-score data frame. The fitted model is attached as the
#'   `model` attribute.
#' @examples
#' if (requireNamespace("mixOmics", quietly = TRUE)) {
#'   outcome <- factor(rep(c("Normal", "Tumor"), each = 12))
#'   dn <- list(paste0("S", seq_len(24)), paste0("F", seq_len(8)))
#'   rna <- matrix(withr::with_seed(7, stats::rnorm(24 * 8)), 24, dimnames = dn)
#'   protein <- rna + withr::with_seed(8, stats::rnorm(24 * 8, sd = 0.5))
#'   head(sb_multiblock_plsda_scores(list(RNA = rna, Protein = protein), outcome), 3)
#' }
#' @export
sb_multiblock_plsda_scores <- function(blocks, outcome, ncomp = 2L, design = 0.1) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Install mixOmics to run the multiblock PLS-DA comparator.", call. = FALSE)
  }
  if (!is.list(blocks) || length(blocks) < 2L || is.null(names(blocks)) || any(!nzchar(names(blocks)))) {
    stop("blocks must be a named list containing at least two matrices.", call. = FALSE)
  }
  blocks <- lapply(blocks, function(x) {
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    x
  })
  sample_counts <- vapply(blocks, nrow, integer(1))
  if (length(unique(sample_counts)) != 1L || sample_counts[[1L]] != length(outcome)) {
    stop("Every block and outcome must contain the same number of samples.", call. = FALSE)
  }
  reference_rows <- rownames(blocks[[1L]])
  if (is.null(reference_rows) || any(vapply(blocks, function(x) {
    is.null(rownames(x)) || !identical(rownames(x), reference_rows)
  }, logical(1)))) {
    stop("Every block must use identical, non-null sample row names.", call. = FALSE)
  }
  shared_features <- Reduce(intersect, lapply(blocks, colnames))
  if (length(shared_features) < 2L) {
    stop("At least two named features must be shared across blocks.", call. = FALSE)
  }
  blocks <- lapply(blocks, function(x) x[, shared_features, drop = FALSE])
  if (any(!vapply(blocks, function(x) all(is.finite(x)), logical(1)))) {
    stop("blocks must contain finite values; impute or filter missing values before fitting.", call. = FALSE)
  }
  if (length(unique(outcome)) < 2L) stop("outcome must contain at least two classes.", call. = FALSE)
  if (!is.numeric(ncomp) || length(ncomp) != 1L || ncomp < 1L) {
    stop("ncomp must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(design) || length(design) != 1L || !is.finite(design) || design < 0 || design > 1) {
    stop("design must be a finite value between 0 and 1.", call. = FALSE)
  }

  design_matrix <- matrix(
    design,
    nrow = length(blocks), ncol = length(blocks),
    dimnames = list(names(blocks), names(blocks))
  )
  diag(design_matrix) <- 0
  fit <- suppressMessages(mixOmics::block.plsda(
    X = blocks,
    Y = factor(outcome),
    ncomp = as.integer(ncomp),
    design = design_matrix,
    scale = TRUE,
    near.zero.var = FALSE
  ))

  block_scores <- lapply(names(blocks), function(block_name) {
    loading <- as.matrix(fit$loadings[[block_name]])
    loading <- loading[shared_features, seq_len(min(ncol(loading), as.integer(ncomp))), drop = FALSE]
    strength <- sqrt(rowSums(loading^2))
    stats::setNames(rank(strength, ties.method = "average") / length(strength), shared_features)
  })
  names(block_scores) <- names(blocks)
  score_matrix <- do.call(cbind, block_scores)
  combined <- exp(rowMeans(log(pmax(score_matrix, 1 / length(shared_features)))))
  out <- data.frame(feature = shared_features, MBPLSDA = combined, stringsAsFactors = FALSE)
  for (block_name in names(blocks)) {
    out[[paste0("MBPLSDA_", block_name)]] <- score_matrix[, block_name]
  }
  attr(out, "model") <- fit
  out[order(out$MBPLSDA, decreasing = TRUE), , drop = FALSE]
}
