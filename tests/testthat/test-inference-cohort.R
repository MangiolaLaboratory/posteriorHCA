library(testthat)

skip_if_not_installed("edgeR")
skip_if_not_installed("cli")

suppressPackageStartupMessages({
  library(cli)
  library(edgeR)
})

pkg_root <- testthat::test_path("..", "..")
if (!dir.exists(file.path(pkg_root, "R"))) {
  pkg_root <- "."
}
sys.source(file.path(pkg_root, "R", "cohort.R"), envir = environment())
sys.source(file.path(pkg_root, "R", "gene-id.R"), envir = environment())
sys.source(file.path(pkg_root, "R", "inference-cohort.R"), envir = environment())
utils_file <- file.path(pkg_root, "R", "utlis.R")
if (file.exists(utils_file)) {
  sys.source(utils_file, envir = environment())
}

toy_counts <- function() {
  genes <- paste0("g", 1:10)
  set.seed(42)
  user <- matrix(
    rnbinom(10 * 6, mu = 50, size = 10),
    nrow = 10,
    dimnames = list(genes, paste0("s", 1:6))
  )
  ref <- setNames(as.numeric(rnbinom(10, mu = 200, size = 10)), genes)
  list(user = user, ref = ref)
}

test_that("draw_dirichlet_weights sums to n and are positive", {
  w <- draw_dirichlet_weights(5L)
  expect_equal(length(w), 5L)
  expect_equal(sum(w), 5, tolerance = 1e-6)
  expect_true(all(w > 0))
})

test_that("bootstrap_cohort_logmu returns n_boot draws", {
  toy <- toy_counts()
  aligned <- suppressMessages(scale_to_hca_reference(toy$user, toy$ref))
  draws <- bootstrap_cohort_logmu(
    counts = aligned,
    gene = "g1",
    n_boot = 50L,
    seed = 123
  )
  expect_equal(length(draws), 50L)
  expect_true(all(is.finite(draws)))
})

test_that("bootstrap_cohort_logmu increases when counts increase", {
  genes <- paste0("g", 1:10)
  set.seed(42)
  counts_low <- matrix(rnbinom(10 * 5, mu = 10, size = 10), nrow = 10, dimnames = list(genes, paste0("s", 1:5)))
  counts_high <- matrix(rnbinom(10 * 5, mu = 500, size = 10), nrow = 10, dimnames = list(genes, paste0("s", 1:5)))
  offset <- rep(0, 5)

  draws_low <- bootstrap_cohort_logmu(counts_low, offset = offset, dispersion = 0.1, gene = "g1", n_boot = 50L, seed = 1)
  draws_high <- bootstrap_cohort_logmu(counts_high, offset = offset, dispersion = 0.1, gene = "g1", n_boot = 50L, seed = 1)

  expect_gt(median(draws_high), median(draws_low))
})

test_that("bootstrap_cohort_logmu supports group subsetting", {
  toy <- toy_counts()
  aligned <- suppressMessages(scale_to_hca_reference(toy$user, toy$ref))
  group <- c("A", "A", "A", "B", "B", "B", "reference")

  draws_a <- bootstrap_cohort_logmu(
    aligned,
    dispersion = 0.1,
    group = group,
    gene = "g1",
    n_boot = 30L,
    seed = 42
  )
  expect_equal(length(draws_a), 30L)
  expect_true(all(is.finite(draws_a)))
})

test_that("bootstrap_cohort_logmu_batch returns one row per cohort", {
  toy <- toy_counts()
  aligned <- suppressMessages(scale_to_hca_reference(toy$user, toy$ref))
  group <- c("A", "A", "A", "B", "B", "B", "reference")

  boot_est <- bootstrap_cohort_logmu_batch(
    aligned,
    group = group,
    genes = "g1",
    dispersion = 0.1,
    n_boot = 30L,
    seed = 42L
  )

  expect_equal(nrow(boot_est), 2L)
  expect_setequal(boot_est$group, c("A", "B"))
  expect_equal(boot_est$method, c("bootstrap", "bootstrap"))
  expect_true(all(c("log_mu", "se", "boot_q025", "boot_q975") %in% names(boot_est)))
})

test_that("bootstrap_cohort_logmu_batch works with cohort_est from estimate_cohort_logmu", {
  toy <- toy_counts()
  aligned <- suppressMessages(scale_to_hca_reference(toy$user, toy$ref))
  meta <- aligned_fields(aligned)$sample_metadata
  meta$cohort <- c("A", "A", "A", "B", "B", "B", "reference")
  est <- estimate_cohort_logmu(
    aligned,
    metadata = meta,
    formula = ~ 0 + cohort,
    genes = "g1"
  )

  boot_est <- bootstrap_cohort_logmu_batch(
    aligned,
    group = c("A", "A", "A", "B", "B", "B", "reference"),
    cohort_est = est,
    n_boot = 20L,
    seed = 1L
  )

  expect_equal(nrow(boot_est), 2L)
  expect_true(all(boot_est$gene == "g1"))
})

test_that("bootstrap_cohort_logmu_batch uses cohort_est dispersion when provided", {
  toy <- toy_counts()
  aligned <- suppressMessages(scale_to_hca_reference(toy$user, toy$ref))
  meta <- aligned_fields(aligned)$sample_metadata
  meta$cohort <- c("A", "A", "A", "B", "B", "B", "reference")
  est <- estimate_cohort_logmu(
    aligned,
    metadata = meta,
    formula = ~ 0 + cohort,
    genes = "g1"
  )

  boot_est <- bootstrap_cohort_logmu_batch(
    aligned,
    group = c("A", "A", "A", "B", "B", "B", "reference"),
    cohort_est = est,
    genes = "g1",
    n_boot = 20L,
    seed = 1L
  )

  expect_true(all(is.finite(boot_est$dispersion)))
})

test_that("welch_t_test_cohort_hca accepts bootstrap_cohort_logmu_batch output", {
  hca_draws <- list(
    draws = rnorm(500, mean = 2.0, sd = 0.2),
    cell_type = "monocytic",
    gene_ensg = "ENSG00000169252",
    gene_symbol = "ADRB2"
  )

  cohort_est <- data.frame(
    gene = rep("ENSG00000169252", 3),
    gene_symbol = rep("ADRB2", 3),
    cell_type = rep("monocytic", 3),
    group = c("disease", "control", "suppressed"),
    log_mu = c(4.0, 2.02, 0.5),
    se = c(0.1, 0.15, 0.1),
    n = c(10, 10, 10),
    stringsAsFactors = FALSE
  )

  res <- welch_t_test_cohort_hca(cohort_est, hca_draws, exclude_groups = character(0))

  expect_equal(nrow(res), 3L)
  expect_equal(res$direction[res$cohort == "disease"], "above_hca")
  expect_equal(res$direction[res$cohort == "control"], "consistent_with_hca")
  expect_equal(res$direction[res$cohort == "suppressed"], "below_hca")
  expect_lt(res$p_value[res$cohort == "disease"], 0.001)
  expect_equal(res$empirical_rank[res$cohort == "disease"], 1)
})

test_that("welch_t_test_cohort_hca excludes reference by default", {
  hca_draws <- list(draws = rnorm(100, 2, 0.1), gene_ensg = "ENSG00000169252")
  cohort_est <- data.frame(
    gene = c("ENSG00000169252", "ENSG00000169252"),
    group = c("reference", "SAVI"),
    log_mu = c(2, 3),
    se = c(0.1, 0.1),
    n = c(1, 5),
    stringsAsFactors = FALSE
  )

  res <- welch_t_test_cohort_hca(cohort_est, hca_draws)
  expect_equal(nrow(res), 1L)
  expect_equal(res$cohort, "SAVI")
})

test_that("welch_t_test_cohort_hca errors on metadata mismatch", {
  hca_draws <- list(
    draws = rnorm(100, 2, 0.1),
    cell_type = "monocytic",
    gene_ensg = "ENSG00000169252"
  )
  cohort_est <- data.frame(
    gene = "ENSG00000073756",
    cell_type = "monocytic",
    group = "A",
    log_mu = 2,
    se = 0.1,
    n = 5,
    stringsAsFactors = FALSE
  )
  expect_error(
    welch_t_test_cohort_hca(cohort_est, hca_draws),
    "Mismatched gene_ensg"
  )
})

test_that("welch_t_test_cohort_hca accepts cohort subset", {
  hca_draws <- list(draws = rnorm(200, 2, 0.1), gene_ensg = "ENSG00000169252")
  cohort_est <- data.frame(
    gene = rep("ENSG00000169252", 2),
    group = c("A", "B"),
    log_mu = c(2.5, 3.5),
    se = c(0.1, 0.1),
    n = c(5, 5),
    stringsAsFactors = FALSE
  )

  res <- welch_t_test_cohort_hca(cohort_est, hca_draws, cohorts = "B")
  expect_equal(nrow(res), 1L)
  expect_equal(res$cohort, "B")
})

test_that("welch_test_means returns expected statistics", {
  out <- welch_test_means(3, 0.1, 2, 0.1, n1 = 10, n2 = 500)
  expect_equal(out$delta, 1)
  expect_equal(out$se_diff, sqrt(0.02))
  expect_lt(out$p_value, 0.001)
  expect_true(is.finite(out$df))
})

test_that("summarize_posterior_draws computes rank", {
  draws <- rnorm(100, mean = 2, sd = 0.2)
  out <- summarize_posterior_draws(draws, value = 2.5)
  expect_equal(out$n, 100L)
  expect_true(out$empirical_rank > 0.9)
})

test_that("manual workflow matches welch_t_test_cohort_hca", {
  set.seed(1)
  hca_draws <- list(
    draws = rnorm(200, mean = 2.0, sd = 0.2),
    cell_type = "monocytic",
    gene_ensg = "ENSG00000169252",
    gene_symbol = "ADRB2"
  )
  cohort_est <- data.frame(
    gene = "ENSG00000169252",
    gene_symbol = "ADRB2",
    cell_type = "monocytic",
    group = "SAVI",
    log_mu = 3.0,
    se = 0.12,
    n = 5,
    stringsAsFactors = FALSE
  )

  auto <- welch_t_test_cohort_hca(cohort_est, hca_draws)
  cohort <- cohort_estimate_at(cohort_est, group = "SAVI")
  baseline <- summarize_posterior_draws(hca_draws, value = cohort$mu)
  manual <- welch_test_means(
    cohort$mu, cohort$se, baseline$mean, baseline$sd,
    n1 = cohort$n, n2 = baseline$n
  )

  expect_equal(auto$cohort_log_mu, manual$mu1)
  expect_equal(auto$hca_mean, manual$mu2)
  expect_equal(auto$t_stat, manual$t_stat)
  expect_equal(auto$p_value, manual$p_value)
  expect_equal(auto$empirical_rank, baseline$empirical_rank)
})
