library(testthat)
library(posteriorHCA)
library(mvtnorm)

# Test core compression functions that don't require model fits

test_that("sample_from_mclust requires valid compressed_rds", {
  expect_error(
    sample_from_mclust("nonexistent_file.rds"),
    "cannot open"
  )
})

test_that("sample_from_mclust handles n_draws parameter", {
  # Create a minimal compressed object for testing
  compressed <- list(
    param_names = c("alpha", "beta"),
    n_components = 2,
    weights = c(0.5, 0.5),
    means = matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2),
    covariances = array(c(1, 0, 0, 1, 1, 0, 0, 1), dim = c(2, 2, 2)),
    n_draws = 1000,
    method = "mclust",
    model_name = "VVV",
    loglik = -100,
    bic = -200
  )
  
  temp_file <- tempfile(fileext = ".rds")
  saveRDS(compressed, temp_file)
  on.exit(unlink(temp_file))
  
  # Test with default n_draws
  samples_default <- sample_from_mclust(temp_file)
  expect_equal(nrow(samples_default), 1000)
  expect_equal(ncol(samples_default), 2)
  expect_equal(colnames(samples_default), c("alpha", "beta"))
  
  # Test with custom n_draws
  samples_custom <- sample_from_mclust(temp_file, n_draws = 500)
  expect_equal(nrow(samples_custom), 500)
  expect_equal(ncol(samples_custom), 2)
})

test_that("sample_from_mclust handles EEE covariance structure", {
  # Create compressed object with EEE structure (single shared covariance)
  compressed <- list(
    param_names = c("alpha", "beta"),
    n_components = 2,
    weights = c(0.5, 0.5),
    means = matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2),
    covariances = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2),  # Single matrix
    n_draws = 100,
    method = "mclust",
    model_name = "EEE",
    loglik = -100,
    bic = -200
  )
  
  temp_file <- tempfile(fileext = ".rds")
  saveRDS(compressed, temp_file)
  on.exit(unlink(temp_file))
  
  samples <- sample_from_mclust(temp_file, n_draws = 100)
  expect_equal(nrow(samples), 100)
  expect_equal(ncol(samples), 2)
})

test_that("sample_from_mclust handles VVV covariance structure", {
  # Create compressed object with VVV structure (array of matrices)
  compressed <- list(
    param_names = c("alpha", "beta"),
    n_components = 2,
    weights = c(0.5, 0.5),
    means = matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2),
    covariances = array(c(1, 0, 0, 1, 1, 0, 0, 1), dim = c(2, 2, 2)),  # Array
    n_draws = 100,
    method = "mclust",
    model_name = "VVV",
    loglik = -100,
    bic = -200
  )
  
  temp_file <- tempfile(fileext = ".rds")
  saveRDS(compressed, temp_file)
  on.exit(unlink(temp_file))
  
  samples <- sample_from_mclust(temp_file, n_draws = 100)
  expect_equal(nrow(samples), 100)
  expect_equal(ncol(samples), 2)
})

test_that("sample_from_mclust handles unknown covariance structure", {
  compressed <- list(
    param_names = c("alpha", "beta"),
    n_components = 1,
    weights = 1,
    means = matrix(c(0, 0), nrow = 2, ncol = 1),
    covariances = "invalid",  # Invalid structure
    n_draws = 100,
    method = "mclust",
    model_name = "UNKNOWN",
    loglik = -100,
    bic = -200
  )
  
  temp_file <- tempfile(fileext = ".rds")
  saveRDS(compressed, temp_file)
  on.exit(unlink(temp_file))
  
  expect_error(
    sample_from_mclust(temp_file, n_draws = 100),
    "Unknown covariance structure"
  )
})

test_that("evaluate_posterior_density_mclust requires valid compressed_rds", {
  expect_error(
    evaluate_posterior_density_mclust("nonexistent_file.rds", matrix(1)),
    "cannot open"
  )
})

test_that("evaluate_posterior_density_mclust handles matrix input", {
  compressed <- list(
    param_names = c("alpha", "beta"),
    n_components = 2,
    weights = c(0.5, 0.5),
    means = matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2),
    covariances = array(c(1, 0, 0, 1, 1, 0, 0, 1), dim = c(2, 2, 2)),
    n_draws = 1000,
    method = "mclust",
    model_name = "VVV",
    loglik = -100,
    bic = -200
  )
  
  temp_file <- tempfile(fileext = ".rds")
  saveRDS(compressed, temp_file)
  on.exit(unlink(temp_file))
  
  # Test with matrix input
  test_points <- matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2)
  densities <- evaluate_posterior_density_mclust(temp_file, test_points)
  
  expect_type(densities, "double")
  expect_length(densities, 2)
  expect_true(all(densities >= 0))
})

test_that("evaluate_posterior_density_mclust handles vector input", {
  compressed <- list(
    param_names = c("alpha", "beta"),
    n_components = 2,
    weights = c(0.5, 0.5),
    means = matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2),
    covariances = array(c(1, 0, 0, 1, 1, 0, 0, 1), dim = c(2, 2, 2)),
    n_draws = 1000,
    method = "mclust",
    model_name = "VVV",
    loglik = -100,
    bic = -200
  )
  
  temp_file <- tempfile(fileext = ".rds")
  saveRDS(compressed, temp_file)
  on.exit(unlink(temp_file))
  
  # Test with vector input (should be converted to matrix)
  test_point <- c(0, 0)
  densities <- evaluate_posterior_density_mclust(temp_file, test_point)
  
  expect_type(densities, "double")
  expect_length(densities, 1)
  expect_true(densities >= 0)
})

test_that("evaluate_posterior_density_mclust handles EEE covariance structure", {
  compressed <- list(
    param_names = c("alpha", "beta"),
    n_components = 2,
    weights = c(0.5, 0.5),
    means = matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2),
    covariances = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2),  # Single matrix
    n_draws = 1000,
    method = "mclust",
    model_name = "EEE",
    loglik = -100,
    bic = -200
  )
  
  temp_file <- tempfile(fileext = ".rds")
  saveRDS(compressed, temp_file)
  on.exit(unlink(temp_file))
  
  test_points <- matrix(c(0, 0), nrow = 1, ncol = 2)
  densities <- evaluate_posterior_density_mclust(temp_file, test_points)
  
  expect_type(densities, "double")
  expect_length(densities, 1)
  expect_true(densities >= 0)
})

test_that("evaluate_posterior_density_mclust handles unknown covariance structure", {
  compressed <- list(
    param_names = c("alpha", "beta"),
    n_components = 1,
    weights = 1,
    means = matrix(c(0, 0), nrow = 2, ncol = 1),
    covariances = "invalid",  # Invalid structure
    n_draws = 1000,
    method = "mclust",
    model_name = "UNKNOWN",
    loglik = -100,
    bic = -200
  )
  
  temp_file <- tempfile(fileext = ".rds")
  saveRDS(compressed, temp_file)
  on.exit(unlink(temp_file))
  
  test_points <- matrix(c(0, 0), nrow = 1, ncol = 2)
  expect_error(
    evaluate_posterior_density_mclust(temp_file, test_points),
    "Unknown covariance structure"
  )
})

test_that("sample_from_mclust produces samples with correct dimensions", {
  compressed <- list(
    param_names = c("alpha", "beta", "gamma"),
    n_components = 3,
    weights = c(0.33, 0.33, 0.34),
    means = matrix(rnorm(9), nrow = 3, ncol = 3),
    covariances = array(rnorm(27), dim = c(3, 3, 3)),
    n_draws = 2000,
    method = "mclust",
    model_name = "VVV",
    loglik = -100,
    bic = -200
  )
  
  # Make covariances positive definite
  for (i in 1:3) {
    cov_mat <- compressed$covariances[, , i]
    compressed$covariances[, , i] <- cov_mat %*% t(cov_mat) + diag(0.1, 3)
  }
  
  temp_file <- tempfile(fileext = ".rds")
  saveRDS(compressed, temp_file)
  on.exit(unlink(temp_file))
  
  samples <- sample_from_mclust(temp_file, n_draws = 500)
  
  expect_equal(nrow(samples), 500)
  expect_equal(ncol(samples), 3)
  expect_equal(colnames(samples), c("alpha", "beta", "gamma"))
  expect_true(all(is.finite(samples)))
})

test_that("evaluate_posterior_density_mclust produces valid densities", {
  compressed <- list(
    param_names = c("alpha", "beta"),
    n_components = 2,
    weights = c(0.5, 0.5),
    means = matrix(c(0, 0, 2, 2), nrow = 2, ncol = 2),
    covariances = array(c(1, 0, 0, 1, 1, 0, 0, 1), dim = c(2, 2, 2)),
    n_draws = 1000,
    method = "mclust",
    model_name = "VVV",
    loglik = -100,
    bic = -200
  )
  
  temp_file <- tempfile(fileext = ".rds")
  saveRDS(compressed, temp_file)
  on.exit(unlink(temp_file))
  
  # Test multiple points
  test_points <- matrix(c(0, 0, 1, 1, 2, 2), nrow = 3, ncol = 2, byrow = TRUE)
  densities <- evaluate_posterior_density_mclust(temp_file, test_points)
  
  expect_type(densities, "double")
  expect_length(densities, 3)
  expect_true(all(densities >= 0))
  expect_true(all(is.finite(densities)))
  
  # Density at component means should be relatively high
  # (exact value depends on weights and covariances)
  expect_true(densities[1] > 0)  # At first component mean
  expect_true(densities[3] > 0)  # At second component mean
})

test_that("sample_from_mclust and evaluate_posterior_density_mclust are consistent", {
  compressed <- list(
    param_names = c("alpha", "beta"),
    n_components = 2,
    weights = c(0.5, 0.5),
    means = matrix(c(0, 0, 2, 2), nrow = 2, ncol = 2),
    covariances = array(c(1, 0, 0, 1, 1, 0, 0, 1), dim = c(2, 2, 2)),
    n_draws = 1000,
    method = "mclust",
    model_name = "VVV",
    loglik = -100,
    bic = -200
  )
  
  temp_file <- tempfile(fileext = ".rds")
  saveRDS(compressed, temp_file)
  on.exit(unlink(temp_file))
  
  # Generate samples
  samples <- sample_from_mclust(temp_file, n_draws = 100)
  
  # Evaluate density at sample points
  densities <- evaluate_posterior_density_mclust(temp_file, samples)
  
  expect_length(densities, 100)
  expect_true(all(densities >= 0))
  expect_true(all(is.finite(densities)))
  
  # Most samples should have non-zero density
  expect_true(mean(densities > 0) > 0.9)
})

