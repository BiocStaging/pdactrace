# Positive-control unit tests: the general locked one-shot runner must reach all three
# terminal states on controlled synthetic inputs and be calibrated (PASS only when a
# real planted signal exists). This is what shows "zero PASS" on real data reflects the
# data, not an over-conservative gate.

sb_test_rule <- function() {
  rule <- sb_default_one_shot_rule()
  rule$permutation_replicates <- 300L
  rule$bootstrap_replicates <- 200L
  rule
}

sb_test_fixture <- function(delta, genes = sprintf("G%03d", 1:320), drop_counts = NULL,
                            seed = 20260719L) {
  set.seed(seed)
  frozen <- data.frame(
    feature = sprintf("G%03d", 1:320),
    frozen_direction = rep(c(1L, -1L), length.out = 320L),
    stringsAsFactors = FALSE
  )
  panel <- data.frame(
    feature = frozen$feature[1:20],
    frozen_direction = frozen$frozen_direction[1:20],
    rank_weight = 1 / log2(2:21),
    stringsAsFactors = FALSE
  )
  signal <- panel$feature[1:10]
  samples <- c(sprintf("CP_%02d", 1:20), sprintf("CN_%02d", 1:20), sprintf("HC_%02d", 1:20))
  grp <- sub("_.*$", "", samples)
  if (!is.null(drop_counts)) {
    keep <- unlist(lapply(names(drop_counts), function(g) which(grp == g)[seq_len(drop_counts[[g]])]))
    samples <- samples[keep]; grp <- grp[keep]
  }
  grid <- expand.grid(gene = genes, sample = samples, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid$group <- grp[match(grid$sample, samples)]
  grid$log2_abundance <- stats::rnorm(nrow(grid), 8, 1)
  if (delta != 0) {
    hit <- grid$group == "CP" & grid$gene %in% signal
    dir <- frozen$frozen_direction[match(grid$gene, frozen$feature)]
    grid$log2_abundance[hit] <- grid$log2_abundance[hit] + delta * dir[hit]
  }
  list(abundance = grid, frozen = frozen, panel = panel,
       groups = list(case = "CP", primary_reference = "CN", secondary_reference = "HC"),
       counts = c(CP = 20L, CN = 20L, HC = 20L))
}

run_fixture <- function(delta, ...) {
  f <- sb_test_fixture(delta, ...)
  sb_one_shot(f$abundance, f$frozen, f$panel, f$groups, f$counts, sb_test_rule())
}

test_that("a strongly evaluable cohort reaches terminal PASS through the runner", {
  d <- run_fixture(1.2)
  expect_true(d$statistical_runner_entered)
  expect_identical(d$status, "PASS")
  expect_gt(d$primary$auc, 0.7)
})

test_that("an evaluable cohort with no signal reaches terminal biological FAIL, not INCONCLUSIVE", {
  d <- run_fixture(0)
  expect_true(d$statistical_runner_entered)
  expect_identical(d$status, "FAIL")
})

test_that("a broken group-count schema halts at the ontology gate (INCONCLUSIVE)", {
  d <- run_fixture(1.2, drop_counts = c(HC = 12L))
  expect_false(d$statistical_runner_entered)
  expect_identical(d$status, "INCONCLUSIVE")
  expect_true("ontology" %in% d$failed_gates)
})

test_that("insufficient candidate coverage halts at the coverage gate (INCONCLUSIVE)", {
  d <- run_fixture(1.2, genes = sprintf("G%03d", 1:40))
  expect_false(d$statistical_runner_entered)
  expect_identical(d$status, "INCONCLUSIVE")
  expect_true("broad_coverage" %in% d$failed_gates)
})

test_that("the runner is calibrated: AUROC increases monotonically with the planted signal", {
  aucs <- vapply(c(0, 0.4, 0.8, 1.4), function(x) run_fixture(x)$primary$auc, numeric(1))
  expect_true(all(diff(aucs) > 0))
  expect_lt(aucs[1], 0.7)   # no signal cannot pass
  expect_gt(aucs[4], 0.9)   # strong signal separates
})

test_that("golden equivalence: sb_one_shot() reproduces the PXD055218 CRC recorded terminal decision", {
  # The CRC one-shot (scripts/38_open_pxd055218_crc_once.R) recorded a terminal INCONCLUSIVE
  # with failed gates [ontology, broad_coverage, top10_coverage, top20_coverage, complete_panel]:
  # HC samples were labelled `Normal` (ontology), broad coverage was 108/320, and only 3 of the
  # top-20 panel genes were quantifiable (3/10, 3/20). Feeding those recorded intermediates into
  # sb_one_shot() must yield the identical terminal decision, binding the general runner to the
  # bespoke per-accession audit logic.
  set.seed(1)
  frozen <- data.frame(feature = sprintf("G%03d", 1:320),
                       frozen_direction = rep(c(1L, -1L), 160), stringsAsFactors = FALSE)
  panel <- data.frame(feature = frozen$feature[1:20], frozen_direction = frozen$frozen_direction[1:20],
                      rank_weight = 1 / log2(2:21), stringsAsFactors = FALSE)
  present <- c(frozen$feature[1:3], frozen$feature[21:125])   # 3 panel genes + 105 others = 108 covered
  samples <- c(sprintf("CP_%02d", 1:20), sprintf("CN_%02d", 1:20), sprintf("Normal_%02d", 1:20))
  grp <- sub("_.*$", "", samples)
  grid <- expand.grid(gene = present, sample = samples, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid$group <- grp[match(grid$sample, samples)]
  grid$log2_abundance <- stats::rnorm(nrow(grid), 8, 1)
  d <- sb_one_shot(grid, frozen, panel,
                   list(case = "CP", primary_reference = "CN", secondary_reference = "HC"),
                   c(CP = 20L, CN = 20L, HC = 20L), sb_default_one_shot_rule())
  expect_identical(d$status, "INCONCLUSIVE")
  expect_false(d$statistical_runner_entered)
  expect_setequal(d$failed_gates,
                  c("ontology", "broad_coverage", "top10_coverage", "top20_coverage", "complete_panel"))
  expect_equal(d$broad_coverage, 108)
  expect_equal(d$top20_coverage, 3)
})
