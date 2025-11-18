library(testthat)
library(posteriorHCA)

# Skip tests if brms is not available
skip_if_not_installed("brms")
skip_if_not_installed("cmdstanr")

# Load packages if available
if (requireNamespace("brms", quietly = TRUE)) {
  library(brms)
}
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  library(cmdstanr)
  # Set brms to use cmdstanr as default backend
  options(brms.backend = "cmdstanr")
}

test_that("compress_brmsfit requires brmsfit object", {
  expect_error(
    compress_brmsfit("not_a_brmsfit"),
    "Input must be a brmsfit object"
  )
  
  expect_error(
    compress_brmsfit(NULL),
    "Input must be a brmsfit object"
  )
})

test_that("compress_brmsfit requires cmdstanr backend", {
  # Create a mock brmsfit without cmdstanr
  mock_fit <- structure(
    list(),
    class = "brmsfit"
  )
  
  expect_error(
    compress_brmsfit(mock_fit),
    "brms fit must use cmdstanr backend"
  )
})

test_that("compress_brmsfit returns list with compressed and structure", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  
  # Fit model using kidney data
  data("kidney", package = "brms")
  
  fit <- brm(
    formula = time | cens(censored) ~ age * sex + disease + (1 + age|patient),
    data = kidney,
    family = lognormal(),
    prior = c(
      set_prior("normal(0,5)", class = "b"),
      set_prior("cauchy(0,2)", class = "sd"),
      set_prior("lkj(2)", class = "cor")
    ),
    warmup = 1000,
    iter = 2000,
    chains = 4,
    backend = "cmdstanr",
    control = list(adapt_delta = 0.95),
    refresh = 0,
    silent = 2
  )
  
  # Verify backend is cmdstanr
  if (is.null(fit$fit) || !inherits(fit$fit, "CmdStanMCMC")) {
    skip("brms fit did not use cmdstanr backend (may have fallen back to rstan)")
  }
  
  # Compress
  result <- compress_brmsfit(fit, n_components = 3, remove_csvs = FALSE)
  
  # Check structure
  expect_type(result, "list")
  expect_named(result, c("compressed", "structure"))
  expect_s3_class(result$structure, "brmsfit")
  expect_type(result$compressed, "list")
  expect_named(result$compressed, c("param_names", "n_components", "weights", 
                                    "means", "covariances", "n_draws", "method", 
                                    "model_name", "loglik", "bic"))
  
  # Check compressed object structure
  expect_type(result$compressed$param_names, "character")
  expect_equal(result$compressed$n_components, 3)
  expect_equal(result$compressed$method, "mclust")
  expect_true(result$compressed$n_draws > 0)
})

test_that("compress_brmsfit handles variables parameter", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  
  # Fit model
  data("kidney", package = "brms")
  fit <- brm(
    formula = time | cens(censored) ~ age * sex + disease + (1 + age|patient),
    data = kidney,
    family = lognormal(),
    warmup = 1000,
    iter = 2000,
    chains = 4,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )
  
  # Verify backend is cmdstanr
  if (is.null(fit$fit) || !inherits(fit$fit, "CmdStanMCMC")) {
    skip("brms fit did not use cmdstanr backend (may have fallen back to rstan)")
  }
  
  # Get some variable names
  all_vars <- posterior::variables(fit)
  selected_vars <- all_vars[seq_len(min(5, length(all_vars)))]
  
  # Test that variables parameter works
  result <- compress_brmsfit(fit, variables = selected_vars, n_components = 2, remove_csvs = FALSE)
  expect_true(all(result$compressed$param_names %in% selected_vars))
  expect_equal(length(result$compressed$param_names), length(selected_vars))
})

test_that("compress_brmsfit handles n_components parameter", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  
  # Fit model
  data("kidney", package = "brms")
  fit <- brm(
    formula = time | cens(censored) ~ age * sex + disease + (1 + age|patient),
    data = kidney,
    family = lognormal(),
    warmup = 1000,
    iter = 2000,
    chains = 4,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )
  
  # Verify backend is cmdstanr
  if (is.null(fit$fit) || !inherits(fit$fit, "CmdStanMCMC")) {
    skip("brms fit did not use cmdstanr backend (may have fallen back to rstan)")
  }
  
  # Test different numbers of components
  result_3 <- compress_brmsfit(fit, n_components = 3, remove_csvs = FALSE)
  result_5 <- compress_brmsfit(fit, n_components = 5, remove_csvs = FALSE)
  expect_equal(result_3$compressed$n_components, 3)
  expect_equal(result_5$compressed$n_components, 5)
})

test_that("reconstruct_brms_from_mclust requires correct input format", {
  expect_error(
    reconstruct_brms_from_mclust("not_a_list"),
    "x must be a list with 'compressed' and 'structure' components"
  )
  
  expect_error(
    reconstruct_brms_from_mclust(list(wrong = "name")),
    "x must be a list with 'compressed' and 'structure' components"
  )
  
  expect_error(
    reconstruct_brms_from_mclust(list(compressed = "test")),
    "x must be a list with 'compressed' and 'structure' components"
  )
})

test_that("reconstruct_brms_from_mclust requires brmsfit in structure", {
  # Create invalid structure
  invalid_result <- list(
    compressed = list(
      param_names = c("alpha", "beta"),
      n_components = 3,
      weights = c(0.3, 0.4, 0.3),
      means = matrix(rnorm(6), nrow = 2, ncol = 3),
      covariances = array(rnorm(18), dim = c(2, 2, 3)),
      n_draws = 1000,
      method = "mclust",
      model_name = "VVV",
      loglik = -100,
      bic = -200
    ),
    structure = "not_a_brmsfit"
  )
  
  expect_error(
    reconstruct_brms_from_mclust(invalid_result),
    "x\\$structure must be a brmsfit object"
  )
})

test_that("reconstruct_brms_from_mclust returns brmsfit object", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  
  # Fit model
  data("kidney", package = "brms")
  fit <- brm(
    formula = time | cens(censored) ~ age * sex + disease + (1 + age|patient),
    data = kidney,
    family = lognormal(),
    warmup = 1000,
    iter = 2000,
    chains = 4,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )
  
  # Verify backend is cmdstanr
  if (is.null(fit$fit) || !inherits(fit$fit, "CmdStanMCMC")) {
    skip("brms fit did not use cmdstanr backend (may have fallen back to rstan)")
  }
  
  # Compress
  result <- compress_brmsfit(fit, n_components = 3, remove_csvs = FALSE)
  
  # Reconstruct
  fit_recon <- reconstruct_brms_from_mclust(result)
  expect_s3_class(fit_recon, "brmsfit")
  
  # Check that it has draws
  draws <- posterior::as_draws_array(fit_recon)
  expect_true(!is.null(draws))
  expect_true(posterior::ndraws(draws) > 0)
})

test_that("compress_brmsfit and reconstruct_brms_from_mclust work together", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  
  # Fit model
  data("kidney", package = "brms")
  fit <- brm(
    formula = time | cens(censored) ~ age * sex + disease + (1 + age|patient),
    data = kidney,
    family = lognormal(),
    warmup = 1000,
    iter = 2000,
    chains = 4,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )
  
  # Verify backend is cmdstanr
  if (is.null(fit$fit) || !inherits(fit$fit, "CmdStanMCMC")) {
    skip("brms fit did not use cmdstanr backend (may have fallen back to rstan)")
  }
  
  # Integration test
  # 1. Compress
  result <- compress_brmsfit(fit, n_components = 3, remove_csvs = FALSE)
  expect_type(result, "list")
  expect_named(result, c("compressed", "structure"))
  
  # 2. Reconstruct
  fit_recon <- reconstruct_brms_from_mclust(result)
  expect_s3_class(fit_recon, "brmsfit")
  
  # 3. Check that reconstructed fit has regenerated draws
  draws <- posterior::as_draws_array(fit_recon)
  expect_true(!is.null(draws))
  expect_true(posterior::ndraws(draws) > 0)
  
  # 4. Check parameter names match
  param_names_orig <- posterior::variables(fit)
  param_names_recon <- posterior::variables(fit_recon)
  expect_true(length(intersect(param_names_orig, param_names_recon)) > 0)
})

test_that("compress_brmsfit handles remove_csvs parameter", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  
  # Fit model
  data("kidney", package = "brms")
  fit <- brm(
    formula = time | cens(censored) ~ age * sex + disease + (1 + age|patient),
    data = kidney,
    family = lognormal(),
    warmup = 1000,
    iter = 2000,
    chains = 4,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )
  
  # Verify backend is cmdstanr
  if (is.null(fit$fit) || !inherits(fit$fit, "CmdStanMCMC")) {
    skip("brms fit did not use cmdstanr backend (may have fallen back to rstan)")
  }
  
  # Test that CSV files are preserved by default
  csv_files_before <- fit$fit$output_files()
  result <- compress_brmsfit(fit, remove_csvs = FALSE)
  csv_files_after <- fit$fit$output_files()
  expect_true(all(file.exists(csv_files_after)))
})

test_that("reconstruct_brms_from_mclust handles n_draws parameter", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  
  # Fit model
  data("kidney", package = "brms")
  fit <- brm(
    formula = time | cens(censored) ~ age * sex + disease + (1 + age|patient),
    data = kidney,
    family = lognormal(),
    warmup = 1000,
    iter = 2000,
    chains = 4,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )
  
  # Verify backend is cmdstanr
  if (is.null(fit$fit) || !inherits(fit$fit, "CmdStanMCMC")) {
    skip("brms fit did not use cmdstanr backend (may have fallen back to rstan)")
  }
  
  # Compress
  result <- compress_brmsfit(fit, n_components = 3, remove_csvs = FALSE)
  
  # Test custom number of draws
  fit_recon_1000 <- reconstruct_brms_from_mclust(result, n_draws = 1000)
  draws_1000 <- posterior::as_draws_array(fit_recon_1000)
  expect_equal(posterior::ndraws(draws_1000), 1000)
})

test_that("compressed object can be saved and loaded", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  
  # Fit model
  data("kidney", package = "brms")
  fit <- brm(
    formula = time | cens(censored) ~ age * sex + disease + (1 + age|patient),
    data = kidney,
    family = lognormal(),
    warmup = 1000,
    iter = 2000,
    chains = 4,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )
  
  # Verify backend is cmdstanr
  if (is.null(fit$fit) || !inherits(fit$fit, "CmdStanMCMC")) {
    skip("brms fit did not use cmdstanr backend (may have fallen back to rstan)")
  }
  
  # Compress
  result <- compress_brmsfit(fit, n_components = 3, remove_csvs = FALSE)
  
  # Save
  temp_compressed <- tempfile(fileext = ".rds")
  temp_structure <- tempfile(fileext = ".rds")
  saveRDS(result$compressed, temp_compressed, compress = "xz")
  saveRDS(result$structure, temp_structure)
  on.exit(unlink(c(temp_compressed, temp_structure)))
  
  # Load
  result_loaded <- list(
    compressed = readRDS(temp_compressed),
    structure = readRDS(temp_structure)
  )
  
  # Reconstruct from loaded objects
  fit_recon <- reconstruct_brms_from_mclust(result_loaded)
  expect_s3_class(fit_recon, "brmsfit")
  
  # Check it works
  draws <- posterior::as_draws_array(fit_recon)
  expect_true(posterior::ndraws(draws) > 0)
})

