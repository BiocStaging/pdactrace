#' Simulate a multi-cohort, multi-omic biomarker benchmark
#'
#' @param n_features Number of simulated features.
#' @param positive_fraction Fraction of true translatable markers.
#' @param scenario One of generic_signal, translation, heterogeneity, dropout,
#'   neutral_observability, mixed_relation, partial_destination, or
#'   anti_observability. Translation scenarios include coherent tissue-only and
#'   inflammation-associated blood decoys; the final three deliberately violate
#'   a default StageBridge assumption.
#' @param seed Random seed.
#' @return List with effects, detection, and feature-level truth.
#' @examples
#' sim <- sb_simulate(n_features = 40L, seed = 1L)
#' vapply(sim, nrow, integer(1))
#' @export
sb_simulate <- function(
    n_features = 1200L,
    positive_fraction = 0.06,
    scenario = c(
      "translation", "generic_signal", "heterogeneity", "dropout",
      "neutral_observability", "mixed_relation", "partial_destination",
      "anti_observability"
    ),
    seed = 1L) {
  scenario <- match.arg(scenario)
  withr::local_seed(seed)
  features <- sprintf("G%05d", seq_len(n_features))
  n_positive <- max(2L, round(n_features * positive_fraction))
  truth <- rep(FALSE, n_features)
  truth[sample.int(n_features, n_positive)] <- TRUE
  inverse_truth <- rep(FALSE, n_features)
  partial_destination_truth <- rep(FALSE, n_features)
  if (scenario == "mixed_relation") {
    inverse_truth[sample(which(truth), floor(n_positive / 2))] <- TRUE
  }
  if (scenario == "partial_destination") {
    partial_destination_truth[sample(which(truth), floor(n_positive / 2))] <- TRUE
  }
  translation_task <- scenario != "generic_signal"
  available <- which(!truth)
  n_decoy <- if (translation_task) round(n_features * 0.08) else 0L
  tissue_decoy <- rep(FALSE, n_features)
  inflammation_decoy <- rep(FALSE, n_features)
  if (n_decoy) {
    tissue_decoy[sample(available, n_decoy)] <- TRUE
    available <- which(!truth & !tissue_decoy)
    inflammation_decoy[sample(available, n_decoy)] <- TRUE
  }
  available <- which(!truth & !tissue_decoy & !inflammation_decoy)
  confounder <- rep(FALSE, n_features)
  confounder[sample(available, min(length(available), round(n_features * 0.06)))] <- TRUE
  direction <- sample(c(-1, 1), n_features, replace = TRUE)
  tissue_signal <- truth | tissue_decoy
  strength <- ifelse(tissue_signal, stats::rlnorm(n_features, log(1.15), 0.25), 0)
  templates <- rbind(c(1, 0.75, 0.55), c(0.7, 1, 0.65), c(0.8, 0.6, 1))
  template_id <- sample(seq_len(nrow(templates)), n_features, replace = TRUE)
  contrasts <- c("early_vs_reference", "mid_vs_early", "late_vs_early")
  modalities <- c("rna", "protein")
  cohorts <- list(rna = paste0("R", 1:4), protein = paste0("P", 1:2))
  rows <- list()
  cursor <- 1L

  for (modality in modalities) {
    modality_scale <- if (modality == "rna") 1 else 0.85
    for (cohort in cohorts[[modality]]) {
      cohort_shift <- stats::rnorm(1, 0, if (scenario == "heterogeneity") 0.35 else 0.12)
      cohort_flip <- if (scenario == "heterogeneity" && stats::runif(1) < 0.18) -1 else 1
      for (contrast_index in seq_along(contrasts)) {
        template_value <- templates[cbind(template_id, contrast_index)]
        mu <- tissue_signal * direction * strength * template_value * modality_scale * cohort_flip
        if (modality == "rna" && cohort == "R1") {
          mu <- mu + confounder * direction * stats::rlnorm(n_features, log(1.6), 0.2)
        }
        se <- stats::runif(n_features, 0.28, if (scenario == "heterogeneity") 0.72 else 0.55)
        estimate <- stats::rnorm(n_features, mu + cohort_shift, se)
        missing_rate <- if (scenario == "dropout" && modality == "protein") 0.55 else 0.08
        observed <- stats::runif(n_features) > missing_rate
        rows[[cursor]] <- data.frame(
          feature = features[observed],
          cohort = cohort,
          modality = modality,
          contrast = contrasts[[contrast_index]],
          estimate = estimate[observed],
          std_error = se[observed],
          stringsAsFactors = FALSE
        )
        cursor <- cursor + 1L
      }
    }
  }
  effects <- do.call(rbind, rows)

  if (translation_task) {
    target_contrasts <- c("tumor_vs_normal", "tumor_vs_inflammation")
    for (cohort in c("B_DISC_1", "B_DISC_2")) {
      cohort_scale <- stats::rnorm(1, 1, if (scenario == "heterogeneity") 0.25 else 0.08)
      for (contrast in target_contrasts) {
        target_scale <- if (contrast == "tumor_vs_normal") 1 else 0.78
        relation_sign <- ifelse(inverse_truth, -1, 1)
        truth_signal <- truth * direction * relation_sign * target_scale * cohort_scale
        if (scenario == "partial_destination" && contrast == "tumor_vs_inflammation") {
          truth_signal[partial_destination_truth] <- 0
        }
        mu <- truth_signal
        if (contrast == "tumor_vs_normal") {
          mu <- mu + inflammation_decoy * direction * 1.05 * cohort_scale
        } else {
          mu <- mu - inflammation_decoy * direction * 0.10 * cohort_scale
        }
        se <- stats::runif(n_features, 0.35, if (scenario == "heterogeneity") 0.78 else 0.62)
        estimate <- stats::rnorm(n_features, mu, se)
        missing_rate <- if (scenario == "dropout") 0.35 else 0.06
        observed <- stats::runif(n_features) > missing_rate
        rows[[cursor]] <- data.frame(
          feature = features[observed],
          cohort = cohort,
          modality = "plasma",
          contrast = contrast,
          estimate = estimate[observed],
          std_error = se[observed],
          stringsAsFactors = FALSE
        )
        cursor <- cursor + 1L
      }
    }
    effects <- do.call(rbind, rows)
  }

  plasma_cohorts <- paste0("B", 1:4)
  depth_offset <- c(-0.8, -0.2, 0.35, 0.8)
  feature_observability <- stats::rnorm(n_features, -0.5, 0.9)
  truth_shift <- if (scenario == "neutral_observability") {
    0
  } else if (scenario == "anti_observability") {
    -1.2
  } else {
    1.6
  }
  detection <- do.call(rbind, lapply(seq_along(plasma_cohorts), function(i) {
    tissue_decoy_shift <- if (scenario == "neutral_observability") {
      0
    } else if (scenario == "anti_observability") {
      0.8
    } else {
      -0.8
    }
    inflammation_shift <- if (scenario == "neutral_observability") 0 else 1.2
    probability <- stats::plogis(
      feature_observability + depth_offset[[i]] +
        truth * truth_shift + tissue_decoy * tissue_decoy_shift +
        inflammation_decoy * inflammation_shift
    )
    data.frame(
      feature = features,
      cohort = plasma_cohorts[[i]],
      detected = stats::runif(n_features) < probability,
      stringsAsFactors = FALSE
    )
  }))
  feature_data <- data.frame(
    feature = features,
    truth = truth,
    tissue_decoy = tissue_decoy,
    inflammation_decoy = inflammation_decoy,
    one_cohort_confounder = confounder,
    inverse_truth = inverse_truth,
    partial_destination_truth = partial_destination_truth,
    stringsAsFactors = FALSE
  )
  list(effects = effects, detection = detection, truth = feature_data)
}
