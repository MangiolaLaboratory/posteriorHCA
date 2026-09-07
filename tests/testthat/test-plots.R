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

test_that("plot_hca_draws returns ggplot from expression_draws-like input", {
  set.seed(1)
  draws_obj <- list(
    draws = rnorm(200, mean = 1.5, sd = 0.1),
    quantity = "linpred",
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

test_that("plot_cohort_vs_hca overlays cohort estimates", {
  set.seed(2)
  draws_obj <- list(
    draws = rnorm(200, mean = 0, sd = 0.2),
    quantity = "linpred",
    gene_ensg = "ENSG00000000001"
  )
  cohort_est <- data.frame(
    group = c("case_a", "case_b"),
    log_mu = c(0.5, -0.4),
    se = c(0.08, 0.07),
    direction = c("above_hca", "below_hca"),
    p_value = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  p <- plot_cohort_vs_hca(draws_obj, cohort_est = cohort_est)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_gt(length(built$data), 1L)
})

test_that("plot_cohort_vs_hca errors for non-linpred draws", {
  draws_obj <- list(draws = rnbinom(100, mu = 10, size = 5), quantity = "predict")
  cohort_est <- data.frame(
    group = "case",
    log_mu = 2,
    se = 0.1,
    stringsAsFactors = FALSE
  )

  expect_error(
    plot_cohort_vs_hca(draws_obj, cohort_est = cohort_est),
    "requires `quantity = \"linpred\"`"
  )
})

test_that("plot_cohort_vs_hca validates cohort_est columns", {
  draws_obj <- list(draws = rnorm(100), quantity = "linpred")
  expect_error(
    plot_cohort_vs_hca(draws_obj, cohort_est = data.frame(x = 1)),
    "must contain `group` and `log_mu`"
  )
})
