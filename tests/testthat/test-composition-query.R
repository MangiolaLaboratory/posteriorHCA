library(testthat)

skip_if_not_installed("sccomp")

pkg_root <- testthat::test_path("..", "..")
if (!dir.exists(file.path(pkg_root, "R"))) {
  pkg_root <- "."
}
sys.source(file.path(pkg_root, "R", "utlis.R"), envir = environment())
sys.source(file.path(pkg_root, "R", "composition-query.R"), envir = environment())

test_that("get_sccomp_ready builds the healthy CellNexus path", {
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
  local_path <- "/home/a1237163/lab/Mangiola_ImmuneAtlas/taskforce_shared_folder/sccomp_on_cellNexus_1_0_12/estimates_age_decade___L3___disease_FALSE___immune_only_TRUE.rds"
  skip_if_not(file.exists(local_path), "Local healthy sccomp model not available")

  fit <- readRDS(local_path)
  nd <- build_sccomp_newdata(
    fit,
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
})

test_that("load_sccomp_fit and composition_draws return draws", {
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
