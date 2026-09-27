test_that("effect validation rejects duplicate cells", {
  x <- data.frame(
    feature = c("A", "A"), cohort = "C1", modality = "rna",
    contrast = "early", estimate = c(1, 2)
  )
  expect_error(sb_validate_effects(x), "duplicate")
})

test_that("coherent replicated evidence outranks a one-cohort signal", {
  x <- expand.grid(
    feature = c("coherent", "single"),
    cohort = c("C1", "C2", "C3"),
    modality = c("rna", "protein"),
    contrast = c("early", "late"),
    stringsAsFactors = FALSE
  )
  x$estimate <- ifelse(x$feature == "coherent", ifelse(x$contrast == "early", 1.2, 0.8), 0)
  x$estimate[x$feature == "single" & x$cohort == "C1" & x$modality == "rna"] <- 3
  x$std_error <- 0.35
  fit <- stagebridge(x, compute_pareto = FALSE)
  expect_lt(fit$stagebridge_rank[fit$feature == "coherent"], fit$stagebridge_rank[fit$feature == "single"])
})

test_that("ranking metrics have known values", {
  truth <- c(TRUE, FALSE, TRUE, FALSE)
  perfect <- c(4, 2, 3, 1)
  expect_equal(sb_binary_auc(truth, perfect), 1)
  expect_equal(sb_average_precision(truth, perfect), 1)
  expect_equal(sb_topk_recall(truth, perfect, 1), 0.5)
})

test_that("average precision is invariant to ordering within score ties", {
  truth <- c(TRUE, FALSE, TRUE, FALSE)
  score <- c(1, 1, 0, 0)
  permutation <- c(2, 1, 4, 3)
  expect_equal(sb_average_precision(truth, score), 0.5)
  expect_equal(
    sb_average_precision(truth, score),
    sb_average_precision(truth[permutation], score[permutation])
  )
})

test_that("method evaluation discloses method-specific score coverage", {
  truth <- c(TRUE, FALSE, TRUE, FALSE)
  scores <- data.frame(
    feature = LETTERS[1:4],
    complete = c(4, 3, 2, 1),
    sparse = c(2, NA, 1, NA)
  )
  metrics <- sb_evaluate_methods(truth, scores, top_k = 1)
  expect_equal(metrics$n_evaluated, c(4, 2))
  expect_equal(metrics$n_positive_evaluated, c(2, 2))
  expect_equal(metrics$prevalence_evaluated, c(0.5, 1))
})

test_that("common-universe evaluation uses identical rows for every method", {
  truth <- c(TRUE, FALSE, TRUE, FALSE)
  scores <- data.frame(
    feature = LETTERS[1:4],
    complete = c(4, 3, 2, 1),
    sparse = c(2, NA, 1, NA)
  )
  metrics <- sb_evaluate_common_methods(truth, scores, top_k = 1)
  expect_true(all(metrics$n_evaluated == 2L))
  expect_true(all(metrics$n_positive_evaluated == 2L))
  expect_true(all(metrics$n_common == 2L))
  expect_true(all(metrics$n_positive_common == 2L))
  expect_true(all(metrics$evaluation_scope == "strict_all_methods_common"))
})

test_that("a shallow-cohort detection carries more information", {
  detection <- rbind(
    data.frame(feature = c("A", "B", "C", "D"), cohort = "shallow", detected = c(TRUE, FALSE, FALSE, FALSE)),
    data.frame(feature = c("A", "B", "C", "D"), cohort = "deep", detected = c(FALSE, TRUE, TRUE, TRUE))
  )
  obs <- sb_observability(detection)
  expect_gt(obs$observability[obs$feature == "A"], obs$observability[obs$feature == "B"])
})

test_that("translation simulations contain target effects and disjoint decoys", {
  simulated <- sb_simulate(n_features = 200, scenario = "translation", seed = 42)
  expect_true("plasma" %in% simulated$effects$modality)
  expect_true(any(simulated$truth$tissue_decoy))
  expect_true(any(simulated$truth$inflammation_decoy))
  expect_false(any(simulated$truth$truth & simulated$truth$tissue_decoy))
  expect_false(any(simulated$truth$truth & simulated$truth$inflammation_decoy))
})

test_that("destination relationship can be predeclared", {
  x <- expand.grid(
    feature = c("direct", "inverse"),
    cohort = c("C1", "C2"),
    modality = c("rna", "protein"),
    contrast = "early_vs_normal",
    stringsAsFactors = FALSE
  )
  x$estimate <- 1
  x$std_error <- 0.3
  plasma <- expand.grid(
    feature = c("direct", "inverse"),
    cohort = c("B1", "B2"),
    modality = "plasma",
    contrast = c("tumor_vs_normal", "tumor_vs_inflammation"),
    stringsAsFactors = FALSE
  )
  plasma$estimate <- ifelse(plasma$feature == "direct", 1, -1)
  plasma$std_error <- 0.3
  direct_fit <- stagebridge(rbind(x, plasma), target_relation = "concordant", compute_pareto = FALSE)
  inverse_fit <- stagebridge(rbind(x, plasma), target_relation = "inverse", compute_pareto = FALSE)
  agnostic_fit <- stagebridge(rbind(x, plasma), target_relation = "agnostic", compute_pareto = FALSE)
  expect_equal(
    direct_fit$stagebridge_score[direct_fit$feature == "direct"],
    inverse_fit$stagebridge_score[inverse_fit$feature == "inverse"],
    tolerance = 1e-10
  )
  expect_equal(
    agnostic_fit$stagebridge_score[agnostic_fit$feature == "direct"],
    agnostic_fit$stagebridge_score[agnostic_fit$feature == "inverse"],
    tolerance = 1e-10
  )
})

test_that("component ablations change exactly one default weight", {
  defaults <- sb_default_weights()
  expect_equal(eval(formals(stagebridge)$weights), defaults)
  for (component in names(defaults)) {
    ablated <- sb_without_component(component)
    expect_identical(names(ablated), names(defaults))
    expect_equal(ablated[[component]], 0)
    expect_equal(ablated[names(defaults) != component], defaults[names(defaults) != component])
  }
  expect_error(sb_without_component("not_a_component"), "exactly one")
})

test_that("fixed weight profiles are normalized and default rescoring is exact", {
  profiles <- sb_weight_profiles()
  expect_equal(unname(vapply(profiles, sum, numeric(1))), rep(1, length(profiles)), tolerance = 1e-12)
  expect_equal(profiles$default, sb_default_weights())
  expect_equal(
    unname(vapply(profiles[c("target_0.20", "target_0.35", "default", "target_0.70", "target_0.85")], `[[`, numeric(1), "target")),
    c(0.20, 0.35, 0.55, 0.70, 0.85)
  )

  simulated <- sb_simulate(n_features = 80, scenario = "translation", seed = 19)
  scores <- sb_method_scores(simulated$effects, simulated$detection)
  components <- attr(scores, "stagebridge_components")
  rescored <- sb_rescore_stagebridge(components)
  expect_equal(rescored, components$stagebridge_score, tolerance = 1e-12)
})

test_that("Pareto and score-percentile tiers are not conflated", {
  x <- expand.grid(
    feature = paste0("F", seq_len(40)),
    cohort = c("C1", "C2"),
    modality = c("rna", "protein"),
    contrast = "early_vs_normal",
    stringsAsFactors = FALSE
  )
  x$estimate <- rep(seq_len(40), each = 4)
  x$std_error <- 1

  score_only <- stagebridge(x, compute_pareto = FALSE)
  expect_true(all(is.na(score_only$pareto_front)))
  expect_true(all(is.na(score_only$pareto_tier)))
  expect_true(all(score_only$rank_basis == "score_only"))
  expect_equal(score_only$stagebridge_rank, score_only$score_rank)
  expect_equal(unname(as.integer(table(score_only$score_percentile_tier))), c(2L, 6L, 32L))

  pareto <- stagebridge(x, compute_pareto = TRUE)
  expect_true(all(is.finite(pareto$pareto_front)))
  expect_false(any(is.na(pareto$pareto_tier)))
  expect_true(all(pareto$rank_basis == "pareto_then_score"))
  expect_equal(nrow(pareto), sum(table(pareto$score_percentile_tier)))
})

test_that("out-of-model simulations encode the declared assumption violations", {
  mixed <- sb_simulate(n_features = 200, scenario = "mixed_relation", seed = 21)
  partial <- sb_simulate(n_features = 200, scenario = "partial_destination", seed = 22)
  anti <- sb_simulate(n_features = 200, scenario = "anti_observability", seed = 23)
  expect_true(any(mixed$truth$truth & mixed$truth$inverse_truth))
  expect_true(any(partial$truth$truth & partial$truth$partial_destination_truth))
  truth_detection <- aggregate(detected ~ feature, anti$detection, mean)
  truth_detection$truth <- anti$truth$truth[match(truth_detection$feature, anti$truth$feature)]
  expect_lt(mean(truth_detection$detected[truth_detection$truth]), mean(truth_detection$detected[!truth_detection$truth]))
})

test_that("multiblock PLS-DA requires aligned matched-sample blocks", {
  skip_if_not_installed("mixOmics")
  set.seed(7)
  outcome <- factor(rep(c("Normal", "Tumor"), each = 12))
  signal <- ifelse(outcome == "Tumor", 2, -2)
  block_a <- matrix(stats::rnorm(24 * 8), nrow = 24, dimnames = list(paste0("S", 1:24), paste0("F", 1:8)))
  block_b <- matrix(stats::rnorm(24 * 8), nrow = 24, dimnames = dimnames(block_a))
  block_a[, "F1"] <- block_a[, "F1"] + signal
  block_b[, "F1"] <- block_b[, "F1"] + signal
  scores <- sb_multiblock_plsda_scores(list(RNA = block_a, Protein = block_b), outcome)
  expect_equal(nrow(scores), 8L)
  expect_true(all(c("feature", "MBPLSDA", "MBPLSDA_RNA", "MBPLSDA_Protein") %in% names(scores)))
  expect_gt(scores$MBPLSDA[scores$feature == "F1"], stats::median(scores$MBPLSDA))
  rownames(block_b)[1] <- "misaligned"
  expect_error(
    sb_multiblock_plsda_scores(list(RNA = block_a, Protein = block_b), outcome),
    "identical"
  )
})

test_that("prospective lock verification detects changed files", {
  skip_if_not_installed("digest")
  root <- tempfile("stagebridge-lock-")
  dir.create(root)
  path <- file.path(root, "frozen.txt")
  writeLines("frozen", path)
  manifest <- sb_lock_manifest(path, root)
  expect_equal(sb_verify_lock_manifest(manifest, root)$status, "PASS")
  writeLines("changed", path)
  expect_true(sb_verify_lock_manifest(manifest, root)$status %in% c("SIZE_MISMATCH", "HASH_MISMATCH"))
})
