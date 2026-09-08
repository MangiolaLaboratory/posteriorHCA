library(testthat)

suppressPackageStartupMessages({
  library(cli)
})

pkg_root <- testthat::test_path("..", "..")
if (!dir.exists(file.path(pkg_root, "R"))) {
  pkg_root <- "."
}
sys.source(file.path(pkg_root, "R", "inference-welch.R"), envir = environment())

test_that("welch_test_means returns expected statistics", {
  out <- welch_test_means(3, 0.1, 2, 0.1, n1 = 10, n2 = 500)
  expect_equal(out$delta, 1)
  expect_true(out$p_value < 0.05)
  expect_true(is.finite(out$df))
})

test_that("welch_test_means works without n1/n2", {
  out <- welch_test_means(3, 0.1, 2, 0.1)
  expect_equal(out$df, Inf)
  expect_true(out$p_value < 0.05)
})

test_that("summarize_posterior_draws computes log_mu, se, n, and rank", {
  draws <- c(1, 2, 3, 4, 5)
  out <- summarize_posterior_draws(draws, value = 2.5)
  expect_equal(out$log_mu, 3)
  expect_equal(out$se, sd(draws))
  expect_equal(out$n, 5L)
  expect_equal(out$empirical_rank, mean(draws <= 2.5))
})

test_that("summarize_posterior_draws accepts expression_draws-like lists", {
  out <- summarize_posterior_draws(list(draws = c(1, 2, 3, 4)))
  expect_equal(out$log_mu, 2.5)
  expect_equal(out$n, 4L)
})

test_that("compare_cohort_to_hca matches summarize + welch_test_means", {
  cohort_estimates <- data.frame(
    gene = c("g1", "g1"),
    group = c("A", "B"),
    n = c(4L, 5L),
    log_mu = c(3.0, 2.5),
    se = c(0.2, 0.15),
    stringsAsFactors = FALSE
  )
  hca_draws <- list(draws = c(1.8, 2.0, 2.1, 1.9, 2.2))

  wrapper <- compare_cohort_to_hca(cohort_estimates, hca_draws)
  expect_equal(nrow(wrapper), 2L)

  hca_summary <- summarize_posterior_draws(hca_draws)
  core_a <- welch_test_means(
    mu1 = cohort_estimates$log_mu[[1]],
    se1 = cohort_estimates$se[[1]],
    mu2 = hca_summary$log_mu,
    se2 = hca_summary$se,
    n1 = cohort_estimates$n[[1]],
    n2 = hca_summary$n
  )
  expect_equal(wrapper$p_value[[1]], core_a$p_value)
  expect_equal(wrapper$t_stat[[1]], core_a$t_stat)
  expect_equal(wrapper$delta[[1]], core_a$delta)
  expect_equal(wrapper$df[[1]], core_a$df)
  expect_equal(wrapper$hca_log_mu[[1]], hca_summary$log_mu)
  expect_equal(wrapper$hca_se[[1]], hca_summary$se)
  expect_equal(wrapper$group, c("A", "B"))
})
