library(testthat)

skip_if_not_installed("cli")

suppressPackageStartupMessages({
  library(cli)
})

pkg_root <- testthat::test_path("..", "..")
if (!dir.exists(file.path(pkg_root, "R"))) {
  pkg_root <- "."
}
sys.source(file.path(pkg_root, "R", "inference-sample.R"), envir = environment())

test_that("test_sample_predictive calculates expected statistics for single sample", {
  set.seed(42)
  # Healthy predicts counts around 10
  predictive_draws <- rnbinom(1000, mu = 10, size = 5)

  # Observed count = 10 (typical)
  res_typical <- test_sample_predictive(
    observed_count = 10,
    offset = 0,
    predictive_draws = predictive_draws,
    sample_id = "patient_1",
    gene = "GENE1"
  )

  expect_equal(res_typical$sample_id, "patient_1")
  expect_equal(res_typical$gene, "GENE1")
  expect_equal(res_typical$observed_count, 10)
  expect_equal(res_typical$direction, "consistent_with_personalised_healthy")
  expect_gt(res_typical$predictive_tail_p, 0.05)
  expect_lt(abs(res_typical$predictive_z), 2.0)

  # Observed count = 100 (extreme high)
  res_high <- test_sample_predictive(
    observed_count = 100,
    offset = 0,
    predictive_draws = predictive_draws,
    sample_id = "patient_high"
  )
  expect_equal(res_high$direction, "above_personalised_healthy")
  expect_lt(res_high$predictive_tail_p, 0.01)
  expect_gt(res_high$predictive_z, 2.5)

  # Observed count = 0 (extreme low if distribution is high)
  predictive_high <- rnbinom(1000, mu = 100, size = 10)
  res_low <- test_sample_predictive(
    observed_count = 0,
    offset = 0,
    predictive_draws = predictive_high,
    sample_id = "patient_low"
  )
  expect_equal(res_low$direction, "below_personalised_healthy")
  expect_lt(res_low$predictive_tail_p, 0.01)
  expect_lt(res_low$predictive_z, -2.5)
})

test_that("test_sample_predictive handles vector of samples and matrix of draws", {
  set.seed(42)
  draws_mat <- cbind(
    s1 = rnbinom(500, mu = 20, size = 5),
    s2 = rnbinom(500, mu = 50, size = 5)
  )
  obs <- c(s1 = 22, s2 = 200)
  offsets <- c(s1 = 0, s2 = 0.5)

  res <- test_sample_predictive(
    observed_count = obs,
    offset = offsets,
    predictive_draws = draws_mat,
    gene = "MYC"
  )

  expect_equal(nrow(res), 2L)
  expect_equal(res$sample_id, c("s1", "s2"))
  expect_equal(res$direction[[1]], "consistent_with_personalised_healthy")
  expect_equal(res$direction[[2]], "above_personalised_healthy")
})

test_that("test_sample_predictive incorporates epred and linpred tests", {
  set.seed(42)
  pred_draws <- rnbinom(500, mu = 20, size = 5)
  epred_draws <- rgamma(500, shape = 20, rate = 1)
  linpred_draws <- rnorm(500, mean = 3.0, sd = 0.15) # log(20) ~ 3.0

  res <- test_sample_predictive(
    observed_count = 80,
    offset = 0.2,
    predictive_draws = pred_draws,
    epred_draws = epred_draws,
    linpred_draws = linpred_draws,
    nb_shape = 5,
    sample_id = "test_sample",
    gene = "IFNG"
  )

  expect_true(!is.na(res$healthy_epred_median))
  expect_true(!is.na(res$log2fc_vs_healthy_epred))
  expect_gt(res$log2fc_vs_healthy_epred, 1.0)

  expect_true(!is.na(res$patient_log_mu))
  expect_true(!is.na(res$welch_style_z))
  expect_gt(res$welch_style_z, 2.0)
  expect_lt(res$welch_style_p, 0.05)
})
