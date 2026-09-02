library(testthat)

skip_if_not_installed("sccomp")
skip_if_not_installed("dittoSeq")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

pkg_root <- testthat::test_path("..", "..")
if (!dir.exists(file.path(pkg_root, "R"))) {
  pkg_root <- "."
}
sys.source(file.path(pkg_root, "R", "utlis.R"), envir = environment())
sys.source(file.path(pkg_root, "R", "composition-query.R"), envir = environment())
sys.source(file.path(pkg_root, "R", "method.R"), envir = environment())

local_fit <- function() {
  local_path <- "/home/a1237163/lab/Mangiola_ImmuneAtlas/taskforce_shared_folder/sccomp_on_cellNexus_1_0_12/estimates_age_decade___L3___disease_FALSE___immune_only_TRUE.rds"
  skip_if_not(file.exists(local_path), "Local healthy sccomp model not available")
  structure(
    list(
      fit = readRDS(local_path),
      release = "cellNexus_1_0_12",
      model = "estimates_age_decade___L3___disease_FALSE___immune_only_TRUE.rds",
      path = local_path
    ),
    class = c("posteriorHCA_sccomp_fit", "list")
  )
}

test_that("composition_posterior_test returns a list", {
  fit <- local_fit()
  result <- composition_posterior_test(
    proportions = data.frame(
      sample_id = c("ID1", "ID2"),
      cell_type = c("b memory", "cd14 mono"),
      proportion = c(0.3, 0.2)
    ),
    sex = "male",
    age_decade = "4",
    ethnicity_groups = "European",
    assay_groups = "10x Genomics 3",
    tissue_groups = "blood",
    fit = fit
  )
  expect_type(result, "list")
  expect_named(result, c("result_table", "plot"))
})

test_that("composition_posterior_test errors for invalid input", {
  fit <- local_fit()
  expect_error(
    composition_posterior_test(
      proportions = data.frame(sample_id = c("ID1", "ID2")),
      fit = fit
    )
  )
  expect_error(
    composition_posterior_test(sex = "invalid_sex", fit = fit)
  )
  expect_error(
    composition_posterior_test(age_decade = "200", fit = fit)
  )
})
