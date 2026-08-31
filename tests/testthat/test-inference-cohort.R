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
    aligned$counts,
    aligned$offset,
    dispersion = 0.1,
    group = group,
    gene = "g1",
    n_boot = 30L,
    seed = 42
  )
  expect_equal(length(draws_a), 30L)
  expect_true(all(is.finite(draws_a)))
})

test_that("test_cohort_vs_hca correctly identifies difference and rank", {
  hca_draws <- rnorm(500, mean = 2.0, sd = 0.2)

  # Cohort significantly above HCA
  cohort_high <- data.frame(gene = "g1", group = "disease", log_mu = 4.0, se = 0.1, n = 10)
  res_high <- test_cohort_vs_hca(cohort_high, hca_draws)

  expect_equal(res_high$direction, "above_hca")
  expect_lt(res_high$p_value, 0.001)
  expect_gt(res_high$t_stat, 5)
  expect_equal(res_high$empirical_rank, 1)

  # Cohort consistent with HCA
  cohort_mid <- data.frame(gene = "g1", group = "control", log_mu = 2.02, se = 0.15, n = 10)
  res_mid <- test_cohort_vs_hca(cohort_mid, hca_draws)

  expect_equal(res_mid$direction, "consistent_with_hca")
  expect_gt(res_mid$p_value, 0.05)

  # Cohort significantly below HCA
  cohort_low <- data.frame(gene = "g1", group = "suppressed", log_mu = 0.5, se = 0.1, n = 10)
  res_low <- test_cohort_vs_hca(cohort_low, hca_draws)

  expect_equal(res_low$direction, "below_hca")
  expect_lt(res_low$p_value, 0.001)
  expect_lt(res_low$t_stat, -5)
  expect_equal(res_low$empirical_rank, 0)
})

test_that("test_cohort_vs_hca accepts bootstrap vector directly", {
  hca_draws <- rnorm(200, mean = 1.5, sd = 0.1)
  cohort_boot <- rnorm(100, mean = 3.0, sd = 0.2)

  res <- test_cohort_vs_hca(cohort_boot, hca_draws, gene = "ACTB", cohort_name = "test_boot")
  expect_equal(res$gene, "ACTB")
  expect_equal(res$cohort, "test_boot")
  expect_equal(res$direction, "above_hca")
  expect_lt(res$p_value, 0.001)
})
