sb_default_weights <- function() {
  c(
    target = 0.55,
    anchor = 0.10,
    evidence = 0.05,
    replication = 0.08,
    concordance = 0.07,
    observability = 0.05,
    coverage = 0.05,
    source = 0.05
  )
}

sb_without_component <- function(component, weights = sb_default_weights()) {
  if (length(component) != 1L || is.na(component) || !component %in% names(weights)) {
    stop("Ablation component must name exactly one available weight.", call. = FALSE)
  }
  out <- weights
  out[[component]] <- 0
  out
}

sb_target_weight_profile <- function(target_weight) {
  if (length(target_weight) != 1L || !is.finite(target_weight) ||
      target_weight <= 0 || target_weight >= 1) {
    stop("Target weight must be one finite value strictly between zero and one.", call. = FALSE)
  }
  defaults <- sb_default_weights()
  other <- defaults[names(defaults) != "target"]
  c(target = target_weight, other * (1 - target_weight) / sum(other))
}

sb_weight_profiles <- function() {
  profiles <- list(
    balanced = stats::setNames(rep(1 / length(sb_default_weights()), length(sb_default_weights())), names(sb_default_weights())),
    target_0.20 = sb_target_weight_profile(0.20),
    target_0.35 = sb_target_weight_profile(0.35),
    default = sb_default_weights(),
    target_0.70 = sb_target_weight_profile(0.70),
    target_0.85 = sb_target_weight_profile(0.85)
  )
  profiles
}

sb_rescore_stagebridge <- function(component_table, weights = sb_default_weights()) {
  component_names <- names(sb_default_weights())
  if (!all(component_names %in% names(component_table))) {
    stop("Component table is missing one or more StageBridge components.", call. = FALSE)
  }
  vapply(seq_len(nrow(component_table)), function(i) {
    components <- unlist(component_table[i, component_names, drop = FALSE], use.names = TRUE)
    sb_geometric_score(components, weights)
  }, numeric(1))
}

#' Rank stage-resolved multi-omic biomarker candidates
#'
#' @param effects Stage-resolved effect table accepted by [sb_validate_effects()].
#' @param detection Optional blood detection table accepted by [sb_observability()].
#' @param source_score Optional named numeric vector or data frame with feature and source_score.
#' @param target_modality Optional destination modality. If omitted, plasma,
#'   serum, or blood is detected automatically.
#' @param target_relation Predeclared relationship between tissue and destination
#'   effects: concordant, inverse, agnostic, or exclude.
#' @param weights Named component weights.
#' @param compute_pareto Whether to calculate exact non-dominated fronts.
#' @return Ranked data frame and auditable component scores. `pareto_tier` is
#'   populated only when `compute_pareto = TRUE`; `score_percentile_tier` is an
#'   explicitly separate score-only top-5%/top-20% tier.
#' @examples
#' sim <- sb_simulate(n_features = 40L, seed = 1L)
#' ranked <- stagebridge(sim$effects, detection = sim$detection)
#' head(ranked[, c("feature", "stagebridge_score", "pareto_front")], 3)
#' @export
stagebridge <- function(
    effects,
    detection = NULL,
    source_score = NULL,
    target_modality = NULL,
    target_relation = c("concordant", "inverse", "agnostic", "exclude"),
    weights = c(
      target = 0.55,
      anchor = 0.10,
      evidence = 0.05,
      replication = 0.08,
      concordance = 0.07,
      observability = 0.05,
      coverage = 0.05,
      source = 0.05
    ),
    compute_pareto = TRUE) {
  target_relation <- match.arg(target_relation)
  x <- sb_prepare_effects(effects)
  if (is.null(target_modality)) {
    target_candidates <- unique(x$modality[grepl("plasma|serum|blood", x$modality, ignore.case = TRUE)])
    target_modality <- if (length(target_candidates)) target_candidates[[1]] else NA_character_
  }
  group <- interaction(x$feature, x$modality, x$contrast, drop = TRUE, lex.order = TRUE)
  layer_rows <- lapply(split(x, group), sb_layer_summary)
  layers <- do.call(rbind, layer_rows)
  rownames(layers) <- NULL
  layers_by_feature <- split(layers, layers$feature)

  all_modalities <- unique(x$modality)
  all_contrasts <- unique(x$contrast)
  features <- sort(unique(x$feature))
  obs <- sb_observability(detection)
  obs_value <- if (is.null(obs)) rep(NA_real_, length(features)) else obs$observability[match(features, obs$feature)]
  source_value <- sb_named_score(source_score, features, "source_score")

  rows <- lapply(seq_along(features), function(i) {
    feature_layers <- layers_by_feature[[features[[i]]]]
    modality_coverage <- length(unique(feature_layers$modality)) / length(all_modalities)
    contrast_coverage <- length(unique(feature_layers$contrast)) / length(all_contrasts)
    coverage <- sqrt(modality_coverage * contrast_coverage)
    anchor_layers <- grepl("normal|reference", feature_layers$contrast, ignore.case = TRUE)
    target_layers <- if (is.na(target_modality)) rep(FALSE, nrow(feature_layers)) else feature_layers$modality == target_modality
    primary_target <- target_layers & grepl("^(tumor|cancer|disease)_vs_", feature_layers$contrast, ignore.case = TRUE)
    if (any(primary_target)) target_layers <- primary_target
    target_score <- if (any(target_layers)) {
      target_evidence <- feature_layers$evidence[target_layers]
      target_direction <- abs(mean(sign(feature_layers$pooled_effect[target_layers]), na.rm = TRUE))
      min(target_evidence, na.rm = TRUE) * target_direction
    } else {
      NA_real_
    }
    tissue_layers <- if (is.na(target_modality)) feature_layers else
      feature_layers[feature_layers$modality != target_modality, , drop = FALSE]
    concordance_score <- if (is.na(target_modality) || !any(target_layers)) {
      sb_cross_modal_concordance(tissue_layers)
    } else if (target_relation == "exclude") {
      sb_cross_modal_concordance(tissue_layers)
    } else {
      direct_layers <- feature_layers
      inverse_layers <- feature_layers
      inverse_layers$pooled_effect[inverse_layers$modality == target_modality] <-
        -inverse_layers$pooled_effect[inverse_layers$modality == target_modality]
      direct_score <- sb_cross_modal_concordance(direct_layers)
      inverse_score <- sb_cross_modal_concordance(inverse_layers)
      if (target_relation == "concordant") direct_score else if (target_relation == "inverse") inverse_score else
        if (all(!is.finite(c(direct_score, inverse_score)))) NA_real_ else
          max(c(direct_score, inverse_score), na.rm = TRUE)
    }
    components <- c(
      target = target_score,
      anchor = if (any(anchor_layers)) mean(feature_layers$evidence[anchor_layers], na.rm = TRUE) else NA_real_,
      evidence = mean(feature_layers$evidence, na.rm = TRUE),
      replication = mean(feature_layers$replication, na.rm = TRUE),
      concordance = concordance_score,
      observability = obs_value[[i]],
      coverage = coverage,
      source = source_value[[i]]
    )
    data.frame(
      feature = features[[i]],
      target = components[["target"]],
      anchor = components[["anchor"]],
      evidence = components[["evidence"]],
      replication = components[["replication"]],
      concordance = components[["concordance"]],
      observability = components[["observability"]],
      coverage = components[["coverage"]],
      source = components[["source"]],
      mean_i2 = if (all(!is.finite(feature_layers$i2))) NA_real_ else mean(feature_layers$i2, na.rm = TRUE),
      n_layers = nrow(feature_layers),
      stagebridge_score = sb_geometric_score(components, weights),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  score_order <- order(-out$stagebridge_score, out$feature)
  out$score_rank <- integer(nrow(out))
  out$score_rank[score_order] <- seq_len(nrow(out))
  pareto_components <- out[c("target", "anchor", "evidence", "replication", "concordance", "observability", "coverage")]
  out$pareto_front <- if (compute_pareto) sb_pareto_fronts(pareto_components) else NA_integer_
  ordering <- if (compute_pareto) {
    order(out$pareto_front, -out$stagebridge_score, out$feature)
  } else {
    order(-out$stagebridge_score, out$feature)
  }
  out <- out[ordering, , drop = FALSE]
  out$stagebridge_rank <- seq_len(nrow(out))
  out$rank_basis <- if (compute_pareto) "pareto_then_score" else "score_only"
  out$pareto_tier <- if (compute_pareto) {
    ifelse(out$pareto_front == 1L, "A", ifelse(out$pareto_front <= 3L, "B", "C"))
  } else {
    NA_character_
  }
  out$score_percentile_tier <- cut(
    out$score_rank / nrow(out),
    c(0, 0.05, 0.20, 1),
    labels = c("A", "B", "C"),
    include.lowest = TRUE
  )
  rownames(out) <- NULL
  attr(out, "layer_summary") <- layers
  out
}

#' Leave-one-cohort-out rank stability
#'
#' @param effects StageBridge effect table.
#' @param detection Optional detection table.
#' @param source_score Optional source score.
#' @param top_k Top-ranked set used for selection frequency.
#' @return Stability summary by feature.
#' @examples
#' sim <- sb_simulate(n_features = 40L, seed = 1L)
#' head(sb_loco_stability(sim$effects, sim$detection, top_k = 10L), 3)
#' @export
sb_loco_stability <- function(effects, detection = NULL, source_score = NULL, top_k = 50L) {
  x <- sb_validate_effects(effects)
  cohorts <- sort(unique(x$cohort))
  if (length(cohorts) < 2L) stop("LOCO stability requires at least two cohorts.", call. = FALSE)
  rankings <- lapply(cohorts, function(held_out) {
    fit <- stagebridge(x[x$cohort != held_out, , drop = FALSE], detection, source_score, compute_pareto = FALSE)
    data.frame(feature = fit$feature, held_out = held_out, rank = fit$stagebridge_rank, stringsAsFactors = FALSE)
  })
  long <- do.call(rbind, rankings)
  rows <- lapply(split(long, long$feature), function(z) {
    data.frame(
      feature = z$feature[[1]],
      median_loco_rank = stats::median(z$rank),
      worst_loco_rank = max(z$rank),
      loco_rank_iqr = stats::IQR(z$rank),
      top_k_frequency = mean(z$rank <= top_k),
      n_loco = nrow(z),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$median_loco_rank, out$worst_loco_rank), ]
  rownames(out) <- NULL
  attr(out, "rankings") <- long
  out
}
