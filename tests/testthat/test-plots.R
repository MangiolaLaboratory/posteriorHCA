library(testthat)

suppressPackageStartupMessages({
  library(cli)
  library(ggplot2)
})

pkg_root <- testthat::test_path("..", "..")
if (!dir.exists(file.path(pkg_root, "R"))) {
  pkg_root <- "."
}
gene_id_file <- file.path(pkg_root, "R", "gene-id.R")
if (file.exists(gene_id_file)) {
  sys.source(gene_id_file, envir = environment())
}
sys.source(file.path(pkg_root, "R", "plots.R"), envir = environment())

test_that("plot_hca_draws returns ggplot from expr_draws-like input", {
  set.seed(1)
  draws_obj <- list(
    draws = rnorm(200, mean = 1.5, sd = 0.1),
    quantity = "linpred",
    gene_symbol = "ADRB2",
    gene_ensg = "ENSG00000169252",
    cell_type = "monocytic"
  )

  p <- plot_hca_draws(draws_obj)
  expect_s3_class(p, "ggplot")
})

test_that("plot_hca_draws builds for predict quantity", {
  set.seed(4)
  draws_obj <- list(
    draws = rnbinom(200, mu = 20, size = 5),
    quantity = "predict"
  )

  p <- plot_hca_draws(draws_obj)
  expect_s3_class(p, "ggplot")
})

test_that("plot_cohort_vs_hca overlays cohort test results", {
  set.seed(2)
  draws_obj <- list(
    draws = rnorm(200, mean = 0, sd = 0.2),
    quantity = "linpred",
    gene_symbol = "GENE1"
  )
  test_tbl <- data.frame(
    cohort = c("case_a", "case_b"),
    cohort_log_mu = c(0.5, -0.4),
    cohort_se = c(0.08, 0.07),
    direction = c("above_hca", "below_hca"),
    p_value = c(0.01, 0.02),
    empirical_rank = c(0.95, 0.05),
    stringsAsFactors = FALSE
  )

  p <- plot_cohort_vs_hca(draws_obj, test_tbl)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_gt(length(built$data), 1L)
})

test_that("plot_cohort_vs_hca errors for non-linpred draws", {
  draws_obj <- list(draws = rnbinom(100, mu = 10, size = 5), quantity = "predict")
  test_tbl <- data.frame(
    cohort = "case",
    cohort_log_mu = 2,
    cohort_se = 0.1,
    stringsAsFactors = FALSE
  )

  expect_error(
    plot_cohort_vs_hca(draws_obj, test_tbl),
    "requires `quantity = \"linpred\"`"
  )
})

test_that("plot_cohort_vs_hca staggers cohort marker heights", {
  set.seed(6)
  draws_obj <- list(draws = rnorm(200, mean = 0, sd = 0.2), quantity = "linpred")
  test_tbl <- data.frame(
    cohort = c("case_a", "case_b", "case_c"),
    cohort_log_mu = c(0.1, 0.12, 0.11),
    cohort_se = c(0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  )

  cohort_df <- normalize_test_results_plot_df(test_tbl)
  staggered <- stagger_cohort_y_positions(cohort_df, draws_obj$draws)
  expect_equal(length(unique(staggered$y)), 3L)
  expect_true(all(diff(staggered$y) < 0))
})

test_that("plot_cohort_vs_hca accepts bootstrap cohort_est", {
  set.seed(6)
  draws_obj <- list(draws = rnorm(200, mean = 0, sd = 0.2), quantity = "linpred")
  boot_est <- data.frame(
    group = "case_a",
    log_mu = 0.5,
    se = 0.08,
    method = "bootstrap",
    empirical_rank = 0.9,
    stringsAsFactors = FALSE
  )

  p <- plot_cohort_vs_hca(draws_obj, cohort_est = boot_est)
  expect_s3_class(p, "ggplot")
})

test_that("plot_cohort_vs_hca validates test_results columns", {
  draws_obj <- list(draws = rnorm(100), quantity = "linpred")
  expect_error(
    plot_cohort_vs_hca(draws_obj, data.frame(x = 1)),
    "must be from"
  )
})
