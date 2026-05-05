library(testthat)
library(posteriorHCA)
library(dplyr)

# Skip tests if sccomp is not available
skip_if_not_installed("sccomp")
skip_if_not_installed("cmdstanr")

# Load sccomp if available
if (requireNamespace("sccomp", quietly = TRUE)) {
  library(sccomp)
}

test_that("compress_sccomp requires object with fit attribute", {
  # Test with object without fit attribute
  obj_no_fit <- structure(list(), class = "sccomp")
  
  expect_error(
    compress_sccomp(obj_no_fit),
    "sccomp object must have cmdstanr fit in attr\\(x, 'fit'\\)"
  )
  
  # Test with NULL fit
  obj_null_fit <- structure(list(), class = "sccomp")
  attr(obj_null_fit, "fit") <- NULL
  
  expect_error(
    compress_sccomp(obj_null_fit),
    "sccomp object must have cmdstanr fit in attr\\(x, 'fit'\\)"
  )
})

test_that("compress_sccomp requires cmdstanr fit in attribute", {
  # Create object with non-cmdstanr fit
  obj_wrong_fit <- structure(list(), class = "sccomp")
  attr(obj_wrong_fit, "fit") <- "not_cmdstanr"
  
  expect_error(
    compress_sccomp(obj_wrong_fit),
    "sccomp fit must use cmdstanr backend"
  )
})

test_that("compress_sccomp returns correct structure", {
  skip_on_cran()
  skip_if_not_installed("sccomp")
  skip_if_not_installed("cmdstanr")
  
  # Generate data and fit sccomp model
  set.seed(123)
  n_samples <- 30  # Smaller for faster testing
  n_cell_types <- 4
  
  sample_data <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    condition = rep(c("A", "B"), each = n_samples/2),
    group = rep(1:3, each = n_samples/3)
  )
  
  counts <- matrix(
    rnbinom(n_samples * n_cell_types, size = 10, mu = 100),
    nrow = n_samples,
    ncol = n_cell_types
  )
  colnames(counts) <- paste0("cell_type_", 1:n_cell_types)
  
  sccomp_data <- counts %>%
    as.data.frame() %>%
    mutate(sample = sample_data$sample_id) %>%
    left_join(sample_data, by = c("sample" = "sample_id"))
  
  # Fit model
  fit_sccomp <- sccomp_model(
    formula = ~ condition + (1 | group),
    data = sccomp_data,
    sample = sample,
    cell_type = starts_with("cell_type_"),
    backend = "cmdstanr",
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    seed = 123,
    refresh = 0,
    silent = 2
  )
  
  # Compress
  result <- compress_sccomp(fit_sccomp, n_components = 3, remove_csvs = FALSE)
  
  # Check structure
  expect_type(result, "list")
  expect_named(result, c("compressed", "structure"))
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

test_that("compress_sccomp removes fit attribute from structure", {
  skip_on_cran()
  skip_if_not_installed("sccomp")
  skip_if_not_installed("cmdstanr")
  
  # Generate data and fit
  set.seed(123)
  n_samples <- 30
  n_cell_types <- 4
  
  sample_data <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    condition = rep(c("A", "B"), each = n_samples/2),
    group = rep(1:3, each = n_samples/3)
  )
  
  counts <- matrix(
    rnbinom(n_samples * n_cell_types, size = 10, mu = 100),
    nrow = n_samples,
    ncol = n_cell_types
  )
  colnames(counts) <- paste0("cell_type_", 1:n_cell_types)
  
  sccomp_data <- counts %>%
    as.data.frame() %>%
    mutate(sample = sample_data$sample_id) %>%
    left_join(sample_data, by = c("sample" = "sample_id"))
  
  fit_sccomp <- sccomp_model(
    formula = ~ condition + (1 | group),
    data = sccomp_data,
    sample = sample,
    cell_type = starts_with("cell_type_"),
    backend = "cmdstanr",
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    seed = 123,
    refresh = 0,
    silent = 2
  )
  
  # Test that structure doesn't have fit attribute
  result <- compress_sccomp(fit_sccomp, n_components = 2, remove_csvs = FALSE)
  expect_null(attr(result$structure, "fit"))
  expect_false(is.null(attr(fit_sccomp, "fit")))  # Original should still have it
})

test_that("compress_sccomp handles variables parameter", {
  skip_on_cran()
  skip_if_not_installed("sccomp")
  skip_if_not_installed("cmdstanr")
  
  # Generate data and fit
  set.seed(123)
  n_samples <- 30
  n_cell_types <- 4
  
  sample_data <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    condition = rep(c("A", "B"), each = n_samples/2),
    group = rep(1:3, each = n_samples/3)
  )
  
  counts <- matrix(
    rnbinom(n_samples * n_cell_types, size = 10, mu = 100),
    nrow = n_samples,
    ncol = n_cell_types
  )
  colnames(counts) <- paste0("cell_type_", 1:n_cell_types)
  
  sccomp_data <- counts %>%
    as.data.frame() %>%
    mutate(sample = sample_data$sample_id) %>%
    left_join(sample_data, by = c("sample" = "sample_id"))
  
  fit_sccomp <- sccomp_model(
    formula = ~ condition + (1 | group),
    data = sccomp_data,
    sample = sample,
    cell_type = starts_with("cell_type_"),
    backend = "cmdstanr",
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    seed = 123,
    refresh = 0,
    silent = 2
  )
  
  # Get fit object to extract variable names
  fit_obj <- attr(fit_sccomp, "fit")
  draws <- fit_obj$draws(format = "matrix")
  all_vars <- colnames(draws)
  selected_vars <- all_vars[seq_len(min(5, length(all_vars)))]
  
  # Test that variables parameter works
  result <- compress_sccomp(fit_sccomp, variables = selected_vars, 
                            n_components = 2, remove_csvs = FALSE)
  expect_true(all(result$compressed$param_names %in% selected_vars))
  expect_equal(length(result$compressed$param_names), length(selected_vars))
})

test_that("compress_sccomp handles n_components parameter", {
  skip_on_cran()
  skip_if_not_installed("sccomp")
  skip_if_not_installed("cmdstanr")
  
  # Generate data and fit
  set.seed(123)
  n_samples <- 30
  n_cell_types <- 4
  
  sample_data <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    condition = rep(c("A", "B"), each = n_samples/2),
    group = rep(1:3, each = n_samples/3)
  )
  
  counts <- matrix(
    rnbinom(n_samples * n_cell_types, size = 10, mu = 100),
    nrow = n_samples,
    ncol = n_cell_types
  )
  colnames(counts) <- paste0("cell_type_", 1:n_cell_types)
  
  sccomp_data <- counts %>%
    as.data.frame() %>%
    mutate(sample = sample_data$sample_id) %>%
    left_join(sample_data, by = c("sample" = "sample_id"))
  
  fit_sccomp <- sccomp_model(
    formula = ~ condition + (1 | group),
    data = sccomp_data,
    sample = sample,
    cell_type = starts_with("cell_type_"),
    backend = "cmdstanr",
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    seed = 123,
    refresh = 0,
    silent = 2
  )
  
  # Test different numbers of components
  result_3 <- compress_sccomp(fit_sccomp, n_components = 3, remove_csvs = FALSE)
  result_5 <- compress_sccomp(fit_sccomp, n_components = 5, remove_csvs = FALSE)
  expect_equal(result_3$compressed$n_components, 3)
  expect_equal(result_5$compressed$n_components, 5)
})

test_that("reconstruct_sccomp_from_mclust requires correct input format", {
  expect_error(
    reconstruct_sccomp_from_mclust("not_a_list"),
    "x must be a list with 'compressed' and 'structure' components"
  )
  
  expect_error(
    reconstruct_sccomp_from_mclust(list(wrong = "name")),
    "x must be a list with 'compressed' and 'structure' components"
  )
  
  expect_error(
    reconstruct_sccomp_from_mclust(list(compressed = "test")),
    "x must be a list with 'compressed' and 'structure' components"
  )
})

test_that("reconstruct_sccomp_from_mclust stores fit_compressed attribute", {
  skip_on_cran()
  skip_if_not_installed("sccomp")
  skip_if_not_installed("cmdstanr")
  
  # Generate data and fit
  set.seed(123)
  n_samples <- 30
  n_cell_types <- 4
  
  sample_data <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    condition = rep(c("A", "B"), each = n_samples/2),
    group = rep(1:3, each = n_samples/3)
  )
  
  counts <- matrix(
    rnbinom(n_samples * n_cell_types, size = 10, mu = 100),
    nrow = n_samples,
    ncol = n_cell_types
  )
  colnames(counts) <- paste0("cell_type_", 1:n_cell_types)
  
  sccomp_data <- counts %>%
    as.data.frame() %>%
    mutate(sample = sample_data$sample_id) %>%
    left_join(sample_data, by = c("sample" = "sample_id"))
  
  fit_sccomp <- sccomp_model(
    formula = ~ condition + (1 | group),
    data = sccomp_data,
    sample = sample,
    cell_type = starts_with("cell_type_"),
    backend = "cmdstanr",
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    seed = 123,
    refresh = 0,
    silent = 2
  )
  
  # Compress
  result <- compress_sccomp(fit_sccomp, n_components = 3, remove_csvs = FALSE)
  
  # Reconstruct
  sccomp_recon <- reconstruct_sccomp_from_mclust(result)
  
  # Test that fit_compressed is stored correctly
  fit_compressed <- attr(sccomp_recon, "fit_compressed")
  expect_type(fit_compressed, "list")
  expect_named(fit_compressed, c("compressed", "regenerated_samples", 
                                 "n_draws", "param_names"))
  expect_null(attr(sccomp_recon, "fit"))  # Original fit should be removed
})

test_that("reconstruct_sccomp_from_mclust removes original fit attribute", {
  skip_on_cran()
  skip_if_not_installed("sccomp")
  skip_if_not_installed("cmdstanr")
  
  # Generate data and fit
  set.seed(123)
  n_samples <- 30
  n_cell_types <- 4
  
  sample_data <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    condition = rep(c("A", "B"), each = n_samples/2),
    group = rep(1:3, each = n_samples/3)
  )
  
  counts <- matrix(
    rnbinom(n_samples * n_cell_types, size = 10, mu = 100),
    nrow = n_samples,
    ncol = n_cell_types
  )
  colnames(counts) <- paste0("cell_type_", 1:n_cell_types)
  
  sccomp_data <- counts %>%
    as.data.frame() %>%
    mutate(sample = sample_data$sample_id) %>%
    left_join(sample_data, by = c("sample" = "sample_id"))
  
  fit_sccomp <- sccomp_model(
    formula = ~ condition + (1 | group),
    data = sccomp_data,
    sample = sample,
    cell_type = starts_with("cell_type_"),
    backend = "cmdstanr",
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    seed = 123,
    refresh = 0,
    silent = 2
  )
  
  # Compress and reconstruct
  result <- compress_sccomp(fit_sccomp, n_components = 2, remove_csvs = FALSE)
  sccomp_recon <- reconstruct_sccomp_from_mclust(result)
  
  # Test that original fit is removed
  expect_null(attr(sccomp_recon, "fit"))
  expect_false(is.null(attr(fit_sccomp, "fit")))  # Original should still have it
})

test_that("reconstruct_sccomp_from_mclust adds metadata attributes", {
  skip_on_cran()
  skip_if_not_installed("sccomp")
  skip_if_not_installed("cmdstanr")
  
  # Generate data and fit
  set.seed(123)
  n_samples <- 30
  n_cell_types <- 4
  
  sample_data <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    condition = rep(c("A", "B"), each = n_samples/2),
    group = rep(1:3, each = n_samples/3)
  )
  
  counts <- matrix(
    rnbinom(n_samples * n_cell_types, size = 10, mu = 100),
    nrow = n_samples,
    ncol = n_cell_types
  )
  colnames(counts) <- paste0("cell_type_", 1:n_cell_types)
  
  sccomp_data <- counts %>%
    as.data.frame() %>%
    mutate(sample = sample_data$sample_id) %>%
    left_join(sample_data, by = c("sample" = "sample_id"))
  
  fit_sccomp <- sccomp_model(
    formula = ~ condition + (1 | group),
    data = sccomp_data,
    sample = sample,
    cell_type = starts_with("cell_type_"),
    backend = "cmdstanr",
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    seed = 123,
    refresh = 0,
    silent = 2
  )
  
  # Compress and reconstruct
  result <- compress_sccomp(fit_sccomp, n_components = 2, remove_csvs = FALSE)
  sccomp_recon <- reconstruct_sccomp_from_mclust(result)
  
  # Test metadata attributes
  expect_equal(attr(sccomp_recon, "compression_method"), "mclust")
  expect_true(attr(sccomp_recon, "reconstructed"))
})

test_that("compress_sccomp and reconstruct_sccomp_from_mclust work together", {
  skip_on_cran()
  skip_if_not_installed("sccomp")
  skip_if_not_installed("cmdstanr")
  
  # Generate data and fit
  set.seed(123)
  n_samples <- 30
  n_cell_types <- 4
  
  sample_data <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    condition = rep(c("A", "B"), each = n_samples/2),
    group = rep(1:3, each = n_samples/3)
  )
  
  counts <- matrix(
    rnbinom(n_samples * n_cell_types, size = 10, mu = 100),
    nrow = n_samples,
    ncol = n_cell_types
  )
  colnames(counts) <- paste0("cell_type_", 1:n_cell_types)
  
  sccomp_data <- counts %>%
    as.data.frame() %>%
    mutate(sample = sample_data$sample_id) %>%
    left_join(sample_data, by = c("sample" = "sample_id"))
  
  fit_sccomp <- sccomp_model(
    formula = ~ condition + (1 | group),
    data = sccomp_data,
    sample = sample,
    cell_type = starts_with("cell_type_"),
    backend = "cmdstanr",
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    seed = 123,
    refresh = 0,
    silent = 2
  )
  
  # Integration test
  # 1. Compress
  result <- compress_sccomp(fit_sccomp, n_components = 3, remove_csvs = FALSE)
  expect_type(result, "list")
  expect_named(result, c("compressed", "structure"))
  
  # 2. Reconstruct
  sccomp_recon <- reconstruct_sccomp_from_mclust(result)
  expect_s3_class(sccomp_recon, "sccomp")
  
  # 3. Check fit_compressed attribute
  fit_compressed <- attr(sccomp_recon, "fit_compressed")
  expect_type(fit_compressed, "list")
  expect_true(nrow(fit_compressed$regenerated_samples) > 0)
})

test_that("reconstruct_sccomp_from_mclust handles n_draws parameter", {
  skip_on_cran()
  skip_if_not_installed("sccomp")
  skip_if_not_installed("cmdstanr")
  
  # Generate data and fit
  set.seed(123)
  n_samples <- 30
  n_cell_types <- 4
  
  sample_data <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    condition = rep(c("A", "B"), each = n_samples/2),
    group = rep(1:3, each = n_samples/3)
  )
  
  counts <- matrix(
    rnbinom(n_samples * n_cell_types, size = 10, mu = 100),
    nrow = n_samples,
    ncol = n_cell_types
  )
  colnames(counts) <- paste0("cell_type_", 1:n_cell_types)
  
  sccomp_data <- counts %>%
    as.data.frame() %>%
    mutate(sample = sample_data$sample_id) %>%
    left_join(sample_data, by = c("sample" = "sample_id"))
  
  fit_sccomp <- sccomp_model(
    formula = ~ condition + (1 | group),
    data = sccomp_data,
    sample = sample,
    cell_type = starts_with("cell_type_"),
    backend = "cmdstanr",
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    seed = 123,
    refresh = 0,
    silent = 2
  )
  
  # Compress
  result <- compress_sccomp(fit_sccomp, n_components = 2, remove_csvs = FALSE)
  
  # Test custom number of draws
  sccomp_recon_1000 <- reconstruct_sccomp_from_mclust(result, n_draws = 1000)
  
  fit_compressed <- attr(sccomp_recon_1000, "fit_compressed")
  expect_equal(nrow(fit_compressed$regenerated_samples), 1000)
})

test_that("compressed sccomp object can be saved and loaded", {
  skip_on_cran()
  skip_if_not_installed("sccomp")
  skip_if_not_installed("cmdstanr")
  
  # Generate data and fit
  set.seed(123)
  n_samples <- 30
  n_cell_types <- 4
  
  sample_data <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    condition = rep(c("A", "B"), each = n_samples/2),
    group = rep(1:3, each = n_samples/3)
  )
  
  counts <- matrix(
    rnbinom(n_samples * n_cell_types, size = 10, mu = 100),
    nrow = n_samples,
    ncol = n_cell_types
  )
  colnames(counts) <- paste0("cell_type_", 1:n_cell_types)
  
  sccomp_data <- counts %>%
    as.data.frame() %>%
    mutate(sample = sample_data$sample_id) %>%
    left_join(sample_data, by = c("sample" = "sample_id"))
  
  fit_sccomp <- sccomp_model(
    formula = ~ condition + (1 | group),
    data = sccomp_data,
    sample = sample,
    cell_type = starts_with("cell_type_"),
    backend = "cmdstanr",
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    seed = 123,
    refresh = 0,
    silent = 2
  )
  
  # Compress
  result <- compress_sccomp(fit_sccomp, n_components = 2, remove_csvs = FALSE)
  
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
  sccomp_recon <- reconstruct_sccomp_from_mclust(result_loaded)
  expect_s3_class(sccomp_recon, "sccomp")
  expect_false(is.null(attr(sccomp_recon, "fit_compressed")))
})

test_that("fit_compressed contains expected components", {
  skip_on_cran()
  skip_if_not_installed("sccomp")
  skip_if_not_installed("cmdstanr")
  
  # Generate data and fit
  set.seed(123)
  n_samples <- 30
  n_cell_types <- 4
  
  sample_data <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    condition = rep(c("A", "B"), each = n_samples/2),
    group = rep(1:3, each = n_samples/3)
  )
  
  counts <- matrix(
    rnbinom(n_samples * n_cell_types, size = 10, mu = 100),
    nrow = n_samples,
    ncol = n_cell_types
  )
  colnames(counts) <- paste0("cell_type_", 1:n_cell_types)
  
  sccomp_data <- counts %>%
    as.data.frame() %>%
    mutate(sample = sample_data$sample_id) %>%
    left_join(sample_data, by = c("sample" = "sample_id"))
  
  fit_sccomp <- sccomp_model(
    formula = ~ condition + (1 | group),
    data = sccomp_data,
    sample = sample,
    cell_type = starts_with("cell_type_"),
    backend = "cmdstanr",
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    seed = 123,
    refresh = 0,
    silent = 2
  )
  
  # Compress and reconstruct
  result <- compress_sccomp(fit_sccomp, n_components = 2, remove_csvs = FALSE)
  sccomp_recon <- reconstruct_sccomp_from_mclust(result)
  
  # Test structure of fit_compressed
  fit_compressed <- attr(sccomp_recon, "fit_compressed")
  expect_named(fit_compressed, c("compressed", "regenerated_samples", 
                                "n_draws", "param_names"))
  expect_type(fit_compressed$compressed, "list")
  expect_true(is.matrix(fit_compressed$regenerated_samples))
  expect_type(fit_compressed$n_draws, "integer")
  expect_type(fit_compressed$param_names, "character")
})

test_that("compress_sccomp handles remove_csvs parameter", {
  skip_on_cran()
  skip_if_not_installed("sccomp")
  skip_if_not_installed("cmdstanr")
  
  # Generate data and fit
  set.seed(123)
  n_samples <- 30
  n_cell_types <- 4
  
  sample_data <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    condition = rep(c("A", "B"), each = n_samples/2),
    group = rep(1:3, each = n_samples/3)
  )
  
  counts <- matrix(
    rnbinom(n_samples * n_cell_types, size = 10, mu = 100),
    nrow = n_samples,
    ncol = n_cell_types
  )
  colnames(counts) <- paste0("cell_type_", 1:n_cell_types)
  
  sccomp_data <- counts %>%
    as.data.frame() %>%
    mutate(sample = sample_data$sample_id) %>%
    left_join(sample_data, by = c("sample" = "sample_id"))
  
  fit_sccomp <- sccomp_model(
    formula = ~ condition + (1 | group),
    data = sccomp_data,
    sample = sample,
    cell_type = starts_with("cell_type_"),
    backend = "cmdstanr",
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    seed = 123,
    refresh = 0,
    silent = 2
  )
  
  # Test that CSV files are preserved by default
  fit_obj <- attr(fit_sccomp, "fit")
  csv_files_before <- fit_obj$output_files()
  result <- compress_sccomp(fit_sccomp, remove_csvs = FALSE)
  csv_files_after <- fit_obj$output_files()
  expect_true(all(file.exists(csv_files_after)))
})

