library(testthat)
library(posteriorHCA)

test_that("get_sub_formula extracts correct terms", {
  full_formula <- as.formula("~ 1 + age_bin + sex + disease_groups + ethnicity_groups")
  factor_names <- c("age_bin", "sex")

  suppressWarnings({
    sub_formula <- get_sub_formula(full_formula, factor_names)
  })

  expect_true(grepl("age_bin", sub_formula))
  expect_true(grepl("sex", sub_formula))
  expect_false(grepl("disease_groups", sub_formula))  # Should not be included
})

test_that("get_sub_formula handles empty factor names", {
  full_formula <- as.formula("~ 1 + age_bin + sex")
  factor_names <- c()

  suppressWarnings({
    sub_formula <- get_sub_formula(full_formula, factor_names)
  })
  expect_equal(sub_formula, "~ 1")
})

test_that("get_sub_formula handles random effect formula", {
  full_formula <- as.formula("~ 1 + age_bin + sex + disease_groups + ethnicity_groups | tissue_groups")
  factor_names <- c()

  suppressWarnings({
    sub_formula <- get_sub_formula(full_formula, factor_names)
  })
  expect_equal(sub_formula, "~ 1")
})

test_that("get_sub_formula does not modify original formula", {
  full_formula <- as.formula("~ 1 + age_bin + sex + age_bin:sex")
  factor_names <- c("age_bin", "sex")

  suppressWarnings({
    get_sub_formula(full_formula, factor_names)
  })

  # Ensure original formula is unchanged
  expect_equal(full_formula, as.formula("~ 1 + age_bin + sex + age_bin:sex"))
})

test_that("get_sub_formula extract random effect formula with proper format", {
  full_formula <- as.formula("~ 1 + age_bin + sex + age_bin:sex + (1 + age_bin + sex + disease_groups + ethnicity_groups | tissue_groups)")
  factor_names <- c("age_bin", "sex", "disease_groups", "ethnicity_groups", "tissue_groups")

  suppressWarnings({
    sub_formula <- get_sub_formula(full_formula, factor_names)
  })

  expect_equal(sub_formula, "~ 1 + age_bin + sex + ( 1 + age_bin + sex + disease_groups + ethnicity_groups | tissue_groups ) + age_bin:sex")
})
