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

test_that("build_newdata_grid still accepts character vectors", {
  fit <- fake_expression_fit()
  assays <- c("10x Genomics 5", "10x Genomics 3")
  grid <- suppressMessages(build_newdata_grid(
    fit,
    age_decade = "7",
    sex = "female",
    disease_groups = "Normal",
    ethnicity_groups = "European",
    assay_groups = assays,
    tissue_groups = "blood"
  ))
  expect_equal(
    as.character(grid$assay_groups_altered),
    c("10x Genomics 5", "10x Genomics 3")
  )
})

test_that("build_newdata_grid filters assay levels with str_detect", {
  fit <- fake_expression_fit()
  grid <- suppressMessages(build_newdata_grid(
    fit,
    age_decade = "7",
    sex = "female",
    disease_groups = "Normal",
    ethnicity_groups = "European",
    assay_groups = str_detect("10x"),
    tissue_groups = "blood"
  ))

  expect_equal(nrow(grid), 2L)
  expect_setequal(
    as.character(grid$assay_groups_altered),
    c("10x Genomics 3", "10x Genomics 5")
  )
  expect_true(all(grid$sex == "female"))
  expect_true(all(grid$age_decade == "7"))
  expect_true("assay_groups" %in% attr(grid, "marginalised"))
})

test_that("str_detect ignores case only when the pattern asks it to", {
  fit <- fake_expression_fit()
  ignore_case <- suppressMessages(build_newdata_grid(
    fit,
    sex = "female",
    age_decade = "7",
    disease_groups = "Normal",
    ethnicity_groups = "European",
    assay_groups = str_detect("(?i)10X"),
    tissue_groups = "blood"
  ))
  expect_setequal(
    as.character(ignore_case$assay_groups_altered),
    c("10x Genomics 3", "10x Genomics 5")
  )

  expect_error(
    build_newdata_grid(fit, assay_groups = str_detect("10X")),
    "No `assay_groups` levels matched"
  )
})

test_that("build_newdata_grid filters numeric-looking levels with between", {
  fit <- fake_expression_fit()
  grid <- suppressMessages(build_newdata_grid(
    fit,
    age_decade = between(6, 8),
    sex = "female",
    disease_groups = "Normal",
    ethnicity_groups = "European",
    assay_groups = "10x Genomics 3",
    tissue_groups = "blood"
  ))
  expect_equal(as.character(grid$age_decade), "7")
})

test_that("build_newdata_grid accepts a caller function and an explicit column", {
  fit <- fake_expression_fit()
  matches_10x <- function(x) grepl("10X", x, ignore.case = TRUE)
  from_function <- suppressMessages(build_newdata_grid(
    fit,
    sex = "female",
    age_decade = "7",
    disease_groups = "Normal",
    ethnicity_groups = "European",
    assay_groups = matches_10x(),
    tissue_groups = "blood"
  ))
  from_expr <- suppressMessages(build_newdata_grid(
    fit,
    sex = "female",
    age_decade = "7",
    disease_groups = "Normal",
    ethnicity_groups = "European",
    assay_groups = grepl("10X", assay_groups, ignore.case = TRUE),
    tissue_groups = "blood"
  ))
  expect_setequal(
    as.character(from_function$assay_groups_altered),
    c("10x Genomics 3", "10x Genomics 5")
  )
  expect_equal(
    as.character(from_function$assay_groups_altered),
    as.character(from_expr$assay_groups_altered)
  )

  qualified <- suppressMessages(build_newdata_grid(
    fit,
    sex = "female",
    age_decade = "7",
    disease_groups = "Normal",
    ethnicity_groups = "European",
    assay_groups = stringr::str_detect("10x"),
    tissue_groups = "blood"
  ))
  expect_setequal(
    as.character(qualified$assay_groups_altered),
    c("10x Genomics 3", "10x Genomics 5")
  )
})

test_that("expression_baseline_draws forwards filter expressions", {
  fit <- fake_expression_fit()
  wrap_env <- environment(expression_baseline_draws)
  original_draws <- wrap_env$expression_draws
  on.exit(assign("expression_draws", original_draws, envir = wrap_env), add = TRUE)
  assign(
    "expression_draws",
    function(fit, newdata, ...) {
      list(grid = newdata, n = nrow(newdata))
    },
    envir = wrap_env
  )

  out <- suppressMessages(expression_baseline_draws(
    fit = fit,
    age_decade = "7",
    sex = "female",
    disease_groups = "Normal",
    ethnicity_groups = "European",
    assay_groups = str_detect("10x"),
    tissue_groups = "blood"
  ))
  expect_setequal(
    as.character(out$grid$assay_groups_altered),
    c("10x Genomics 3", "10x Genomics 5")
  )
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

test_that("expression_baseline_draws requires fit or cell_type + gene_ensg", {
  expect_error(
    expression_baseline_draws(cell_type = "monocytic"),
    "Provide `fit`"
  )
  expect_error(
    expression_baseline_draws(gene_ensg = "ENSG00000169252"),
    "Provide `fit`"
  )
})

test_that("expression_baseline_draws matches build_newdata_grid + expression_draws", {
  fit <- fake_expression_fit()
  wrap_env <- environment(expression_baseline_draws)
  original_draws <- wrap_env$expression_draws
  on.exit(assign("expression_draws", original_draws, envir = wrap_env), add = TRUE)

  stub_draws <- function(
    fit,
    newdata,
    quantity = c("linpred", "predict", "epred"),
    marginalise = c("mean", "pool", "sample"),
    ndraws = NULL,
    transform = FALSE,
    re_formula = NULL,
    allow_new_levels = TRUE,
    sample_new_levels = "gaussian",
    seed = NULL
  ) {
    quantity <- match.arg(quantity)
    marginalise <- match.arg(marginalise)
    list(
      draws = rep(as.numeric(if (is.null(seed)) 0 else seed), nrow(newdata)),
      grid = newdata,
      quantity = quantity,
      marginalise = marginalise,
      n_grid = nrow(newdata),
      cell_type = NA_character_,
      gene_ensg = NA_character_
    )
  }
  assign("expression_draws", stub_draws, envir = wrap_env)

  wrapper_result <- suppressMessages(expression_baseline_draws(
    fit = fit,
    disease_groups = "Normal",
    tissue_groups = "blood",
    assay_groups = "10x Genomics 3",
    seed = 42
  ))

  grid <- suppressMessages(build_newdata_grid(
    fit,
    disease_groups = "Normal",
    tissue_groups = "blood",
    assay_groups = "10x Genomics 3"
  ))
  core_result <- stub_draws(
    fit,
    newdata = grid,
    quantity = "linpred",
    marginalise = "mean",
    seed = 42
  )

  expect_equal(wrapper_result$draws, core_result$draws)
  expect_equal(wrapper_result$grid, core_result$grid)
  expect_equal(wrapper_result$n_grid, core_result$n_grid)
  expect_equal(wrapper_result$quantity, "linpred")
})
