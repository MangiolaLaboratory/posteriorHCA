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

test_that("summarize_posterior_draws computes mean, sd, n, and rank", {
  draws <- c(1, 2, 3, 4, 5)
  out <- summarize_posterior_draws(draws, value = 2.5)
  expect_equal(out$mean, 3)
  expect_equal(out$sd, sd(draws))
  expect_equal(out$n, 5L)
  expect_equal(out$empirical_rank, mean(draws <= 2.5))
})

test_that("summarize_posterior_draws accepts expression_draws-like lists", {
  out <- summarize_posterior_draws(list(draws = c(1, 2, 3, 4)))
  expect_equal(out$mean, 2.5)
  expect_equal(out$n, 4L)
})
