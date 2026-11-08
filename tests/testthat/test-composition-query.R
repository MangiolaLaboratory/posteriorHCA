library(testthat)

suppressPackageStartupMessages({
  library(cli)
  library(tidyr)
  library(dplyr)
  library(rlang)
  library(stringr)
})

pkg_root <- testthat::test_path("..", "..")
if (!dir.exists(file.path(pkg_root, "R"))) {
  pkg_root <- "."
}
sys.source(file.path(pkg_root, "R", "utlis.R"), envir = environment())
sys.source(file.path(pkg_root, "R", "metadata-filter.R"), envir = environment())
sys.source(file.path(pkg_root, "R", "composition-query.R"), envir = environment())

fake_sccomp_fit <- function() {
  cd <- data.frame(
    age_decade = c("4", "5", "7", "4"),
    sex = c("male", "female", "male", "female"),
    ethnicity_groups_imputed = c("European", "European", "African", "East Asian"),
    assay_groups___altered = c(
      "10x Genomics 3",
      "10x Genomics 5",
      "Smart seq",
      "10x Genomics 3"
    ),
    tissue_groups = c("blood", "blood", "liver", "blood"),
    stringsAsFactors = FALSE
  )
  fit <- list()
  attr(fit, "count_data") <- cd
  fit
}

test_that("get_sccomp_ready builds the healthy CellNexus path", {
  skip_if_not_installed("sccomp")
  skip_on_cran()
  url <- nectar_object_url(
    "sccomp_est",
    "cellNexus_1_0_12",
    "estimates_age_decade___L3___disease_FALSE___immune_only_TRUE.rds"
  )
  expect_match(
    url,
    "sccomp_est/cellNexus_1_0_12/estimates_age_decade___L3___disease_FALSE___immune_only_TRUE.rds"
  )
})

test_that("build_sccomp_newdata passes NA metadata through", {
  nd <- build_sccomp_newdata(
    fake_sccomp_fit(),
    sample_id = "s1",
    age_decade = NA,
    sex = "male",
    ethnicity_groups = NA,
    assay_groups = NA,
    tissue_groups = NA
  )

  expect_equal(nd$sample_id, "s1")
  expect_true(is.na(nd$age_decade))
  expect_equal(nd$sex, "male")
  expect_true(is.na(nd$dataset_id___altered))
  expect_equal(nrow(nd), 1L)
})

test_that("build_sccomp_newdata accepts str_detect(\"10x\")", {
  nd <- build_sccomp_newdata(
    fake_sccomp_fit(),
    sample_id = "s1",
    age_decade = "7",
    sex = "male",
    ethnicity_groups = "African",
    assay_groups = str_detect("10x"),
    tissue_groups = "liver"
  )
  expect_equal(nrow(nd), 2L)
  expect_setequal(
    nd$assay_groups___altered,
    c("10x Genomics 3", "10x Genomics 5")
  )
  expect_equal(nd$sample_id, c("s1__1", "s1__2"))
})

test_that("build_sccomp_newdata filters age with between(6, 8)", {
  nd <- build_sccomp_newdata(
    fake_sccomp_fit(),
    sample_id = "s1",
    age_decade = between(6, 8),
    sex = "male",
    ethnicity_groups = "African",
    assay_groups = "Smart seq",
    tissue_groups = "liver"
  )
  expect_equal(nrow(nd), 1L)
  expect_equal(nd$age_decade, "7")
})

test_that("load_sccomp_fit and composition_draws return draws", {
  skip_if_not_installed("sccomp")
  local_path <- "/home/a1237163/lab/Mangiola_ImmuneAtlas/taskforce_shared_folder/sccomp_on_cellNexus_1_0_12/estimates_age_decade___L3___disease_FALSE___immune_only_TRUE.rds"
  skip_if_not(file.exists(local_path), "Local healthy sccomp model not available")

  fit <- structure(
    list(
      fit = readRDS(local_path),
      release = "cellNexus_1_0_12",
      model = "estimates_age_decade___L3___disease_FALSE___immune_only_TRUE.rds",
      path = local_path
    ),
    class = c("posteriorHCA_sccomp_fit", "list")
  )

  out <- composition_draws(
    fit,
    sex = "male",
    age_decade = "4",
    ethnicity_groups = "European",
    assay_groups = "10x Genomics 3",
    tissue_groups = "blood"
  )

  expect_true("draws" %in% names(out))
  expect_true(all(c("cell_type", "proportion") %in% names(out$draws)))
  expect_gt(nrow(out$draws), 0L)
})
