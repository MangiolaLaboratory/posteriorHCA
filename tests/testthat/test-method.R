library(testthat)
library(posteriorHCA)

test_that("posterior_test returns a list", {
  result <- posterior_test(
    proportions = data.frame(sample_id = c('ID1', 'ID2'),cell_type = c("b memory", "cd14 mono"), proportion = c(0.3, 0.2)),
    sex = "male",
    age_bin = "Senior",
    disease_groups = "Normal",
    ethnicity_groups = "European",
    tissue_groups = "blood"
  )
  expect_type(result, "list")
  expect_named(result, c("result_table", "plot"))
})

test_that("posterior_test errors for invalid input", {
  expect_error(posterior_test(proportions = data.frame(sample_id = c('ID1', 'ID2'))))
  expect_error(posterior_test(sex = "invalid_sex"))
  expect_error(posterior_test(age_bin = "200")) # Invalid age category
})
