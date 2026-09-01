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

fake_expr_fit <- function() {
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
  fit <- list(
    data = dat,
    parnames = c(
      "Intercept",
      "assay_groups_altered10x Genomics 5",
      "sexmale"
    )
  )
  fit
}

test_that("levels_from_coef_names strips the column prefix", {
  nms <- c("Intercept", "assay_groups_altered10x Genomics 5", "sexmale")
  expect_equal(
    levels_from_coef_names(nms, "assay_groups_altered"),
    "10x Genomics 5"
  )
  expect_equal(
    levels_from_coef_names(nms, "sex"),
    "male"
  )
})

test_that("expr_model_levels includes the reference category from data", {
  fit <- fake_expr_fit()
  assay_lv <- expr_model_levels(fit, "assay_groups")
  expect_true("10x Genomics 3" %in% assay_lv)
  expect_true("10x Genomics 5" %in% assay_lv)
})

test_that("a single assay value is checked and kept fixed", {
  fit <- fake_expr_fit()
  grid <- suppressMessages(build_newdata_grid(
    fit,
    assay_groups = "10x Genomics 3",
    sex = "female",
    age_decade = "4",
    disease_groups = "Normal",
    ethnicity_groups = "European",
    tissue_groups = "blood"
  ))
  expect_equal(nrow(grid), 1L)
  expect_equal(as.character(grid$assay_groups_altered), "10x Genomics 3")
  expect_equal(attr(grid, "fixed"), c(
    "age_decade", "sex", "disease_groups",
    "ethnicity_groups", "assay_groups", "tissue_groups"
  ))
  expect_equal(attr(grid, "marginalised"), character(0))
})

test_that("several assay values expand only those levels", {
  fit <- fake_expr_fit()
  grid <- suppressMessages(build_newdata_grid(
    fit,
    assay_groups = c("10x Genomics 3", "10x Genomics 5"),
    sex = "female",
    age_decade = "4",
    disease_groups = "Normal",
    ethnicity_groups = "European",
    tissue_groups = "blood"
  ))
  expect_equal(nrow(grid), 2L)
  expect_setequal(
    as.character(grid$assay_groups_altered),
    c("10x Genomics 3", "10x Genomics 5")
  )
  expect_true("assay_groups" %in% attr(grid, "marginalised"))
})

test_that("NA assay expands every assay in the model", {
  fit <- fake_expr_fit()
  grid <- suppressMessages(build_newdata_grid(
    fit,
    assay_groups = NA,
    sex = "female",
    age_decade = "4",
    disease_groups = "Normal",
    ethnicity_groups = "European",
    tissue_groups = "blood"
  ))
  expect_equal(nrow(grid), 2L)
  expect_setequal(
    as.character(unique(grid$assay_groups_altered)),
    expr_model_levels(fit, "assay_groups")
  )
})

test_that("an unknown assay value errors with the allowed list", {
  fit <- fake_expr_fit()
  expect_error(
    build_newdata_grid(fit, assay_groups = "Smart-seq2"),
    "Unknown assay_groups"
  )
})

test_that("NA on several covariates takes the Cartesian product", {
  fit <- fake_expr_fit()
  grid <- suppressMessages(build_newdata_grid(
    fit,
    sex = NA,
    age_decade = c("4", "5"),
    disease_groups = "Normal",
    ethnicity_groups = "European",
    assay_groups = "10x Genomics 3",
    tissue_groups = "blood"
  ))
  # sex: female, male; age: 4, 5
  expect_equal(nrow(grid), 4L)
  expect_setequal(as.character(grid$sex), c("female", "male"))
  expect_setequal(as.character(grid$age_decade), c("4", "5"))
})

test_that("dataset_id is not expanded over training studies", {
  fit <- fake_expr_fit()
  grid <- suppressMessages(build_newdata_grid(
    fit,
    sex = "female",
    age_decade = "4",
    disease_groups = "Normal",
    ethnicity_groups = "European",
    assay_groups = "10x Genomics 3",
    tissue_groups = "blood"
  ))
  expect_equal(unique(as.character(grid$dataset_id_altered)), "__new_study__")
})

test_that("factor levels on the grid keep unused training levels", {
  fit <- fake_expr_fit()
  grid <- suppressMessages(build_newdata_grid(
    fit,
    assay_groups = "10x Genomics 3",
    sex = "female",
    age_decade = "4",
    disease_groups = "Normal",
    ethnicity_groups = "European",
    tissue_groups = "blood"
  ))
  expect_true(is.factor(grid$assay_groups_altered))
  expect_equal(
    levels(grid$assay_groups_altered),
    c("10x Genomics 3", "10x Genomics 5")
  )
})

test_that("offset is 0 and counts dummy is present", {
  fit <- fake_expr_fit()
  grid <- suppressMessages(build_newdata_grid(
    fit,
    sex = "female",
    age_decade = "4",
    disease_groups = "Normal",
    ethnicity_groups = "European",
    assay_groups = "10x Genomics 3",
    tissue_groups = "blood"
  ))
  expect_equal(grid$offset, 0)
  expect_equal(grid$counts, 1)
})

test_that("marginalize_draw_matrix mean, pool, and sample have expected length", {
  y <- matrix(
    c(1, 2, 3, 10, 20, 30),
    nrow = 3,
    ncol = 2
  )
  expect_equal(marginalize_draw_matrix(y, "mean"), c(5.5, 11, 16.5))
  expect_equal(length(marginalize_draw_matrix(y, "pool")), 6L)

  set.seed(1)
  sampled <- marginalize_draw_matrix(y, "sample")
  expect_equal(length(sampled), 3L)
  expect_true(all(sampled %in% as.numeric(y)))
})

test_that("marginalize_draw_matrix is a no-op for a single profile", {
  y <- matrix(1:4, ncol = 1)
  expect_equal(marginalize_draw_matrix(y, "mean"), 1:4)
  expect_equal(marginalize_draw_matrix(y, "pool"), 1:4)
})

test_that("expr_draws rejects a bad quantity before calling brms", {
  fit <- fake_expr_fit()
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
    expr_draws(fit, grid, quantity = "log_mu"),
    "arg"
  )
})
