library(testthat)

suppressPackageStartupMessages({
  library(cli)
  library(tidyr)
})

expr_query_file <- file.path(testthat::test_path("..", ".."), "R", "expr-query.R")
if (!file.exists(expr_query_file)) {
  expr_query_file <- file.path("R", "expr-query.R")
}
gene_id_file <- file.path(dirname(expr_query_file), "gene-id.R")
if (file.exists(gene_id_file)) {
  sys.source(gene_id_file, envir = environment())
}
sys.source(expr_query_file, envir = environment())

fake_expression_fit <- function() {
  dat <- data.frame(
    age_decade = factor(c("4", "5", "7", "4")),
    sex = factor(c("female", "male", "female", "male"), levels = c("female", "male")),
    disease_groups_altered = factor(
      c("Normal", "Normal", "COVID-19 related_blood", "Normal")
    ),
    ethnicity_groups = factor(c("European", "East Asian", "European", "African")),
    assay_groups_altered = factor(
      c("10x Genomics 3", "10x Genomics 5", "10x Genomics 3", "10x Genomics 3"),
      levels = c("10x Genomics 3", "10x Genomics 5")
    ),
    tissue_groups = factor(c("blood", "liver", "blood", "blood")),
    dataset_id_altered = factor(c("ds1", "ds2", "ds1", "ds3")),
    offset = c(0, 0.1, 0, 0),
    counts = c(10, 20, 30, 40)
  )
  list(
    data = dat,
    parnames = c(
      "Intercept",
      "assay_groups_altered10x Genomics 5",
      "sexmale"
    )
  )
}

test_that("expression_model_levels includes the reference category from data", {
  fit <- fake_expression_fit()
  assay_lv <- expression_model_levels(fit, "assay_groups")
  expect_true("10x Genomics 3" %in% assay_lv)
  expect_true("10x Genomics 5" %in% assay_lv)
})

test_that("build_newdata_grid expands missing covariates", {
  fit <- fake_expression_fit()
  grid <- suppressMessages(build_newdata_grid(
    fit,
    disease_groups = "Normal",
    tissue_groups = "blood",
    assay_groups = "10x Genomics 3"
  ))
  expect_true(nrow(grid) > 1L)
  expect_true(all(grid$disease_groups_altered == "Normal"))
  expect_true(all(grid$tissue_groups == "blood"))
  expect_true(all(grid$assay_groups_altered == "10x Genomics 3"))
  expect_true(all(grid$dataset_id_altered == "__new_study__"))
})

test_that("build_newdata_grid can fix every covariate", {
  fit <- fake_expression_fit()
  grid <- suppressMessages(build_newdata_grid(
    fit,
    age_decade = "4",
    sex = "female",
    disease_groups = "Normal",
    ethnicity_groups = "European",
    assay_groups = "10x Genomics 3",
    tissue_groups = "blood"
  ))
  expect_equal(nrow(grid), 1L)
  expect_equal(attr(grid, "fixed"), names(.expression_meta_map))
})

test_that("build_newdata_grid errors on unknown levels", {
  fit <- fake_expression_fit()
  expect_error(
    build_newdata_grid(fit, assay_groups = "not-a-real-assay"),
    "Unknown assay_groups"
  )
})

test_that("build_newdata_grid errors when a missing covariate is requested", {
  fit <- fake_expression_fit()
  fit$data$tissue_groups <- NULL
  expect_error(
    build_newdata_grid(fit, tissue_groups = "blood"),
    "no `tissue_groups` column"
  )
})

test_that("expression_draws rejects a bad quantity before calling brms", {
  fit <- fake_expression_fit()
  grid <- suppressMessages(build_newdata_grid(
    fit,
    sex = "female",
    age_decade = "4",
    disease_groups = "Normal",
    ethnicity_groups = "European",
    assay_groups = "10x Genomics 3",
    tissue_groups = "blood"
  ))
  expect_error(
    expression_draws(fit, grid, quantity = "log_mu"),
    "arg"
  )
})
