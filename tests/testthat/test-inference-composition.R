library(testthat)

pkg_root <- testthat::test_path("..", "..")
if (!dir.exists(file.path(pkg_root, "R"))) {
  pkg_root <- "."
}
sys.source(file.path(pkg_root, "R", "inference-composition.R"), envir = environment())

fake_draws <- function() {
  data.frame(
    cell_type = rep(c("b memory", "cd14 mono"), each = 100),
    proportion = c(rnorm(100, 0.05, 0.01), rnorm(100, 0.12, 0.02)),
    stringsAsFactors = FALSE
  )
}

test_that("empirical_proportion_confidence is between 0 and 1", {
  set.seed(1)
  draws <- rnorm(200, 0.1, 0.02)
  expect_gte(empirical_proportion_confidence(0.1, draws), 0)
  expect_lte(empirical_proportion_confidence(0.1, draws), 1)
  expect_lt(empirical_proportion_confidence(0.5, draws), 0.05)
})

test_that("normalize_proportions handles 2- and 3-column input", {
  two <- data.frame(cell = "b memory", p = 0.2)
  out2 <- normalize_proportions(two)
  expect_equal(names(out2), c("sample_id", "cell_type", "proportion"))
  expect_equal(out2$sample_id, "sample_1")

  three <- data.frame(id = "A", cell = "b memory", p = 0.2)
  out3 <- normalize_proportions(three)
  expect_equal(out3$sample_id, "A")
})

test_that("composition_test returns one row per observation", {
  set.seed(2)
  draws <- fake_draws()
  props <- data.frame(
    sample_id = c("S1", "S2"),
    cell_type = c("b memory", "cd14 mono"),
    proportion = c(0.3, 0.2)
  )
  res <- composition_test(props, draws)
  expect_equal(nrow(res), 2L)
  expect_true(all(c(
    "sample_id", "cell_type", "proportion_observed",
    "empirical_confidence", "mean", "lower", "upper"
  ) %in% names(res)))
})

test_that("composition_test_sample subsets to one sample", {
  set.seed(3)
  draws <- fake_draws()
  props <- data.frame(
    sample_id = c("S1", "S1", "S2"),
    cell_type = c("b memory", "cd14 mono", "b memory"),
    proportion = c(0.2, 0.1, 0.3)
  )
  res <- composition_test_sample(props, draws, sample_id = "S1")
  expect_equal(unique(res$sample_id), "S1")
  expect_equal(nrow(res), 2L)
})

test_that("composition_test_cohort mirrors composition_test", {
  set.seed(4)
  draws <- fake_draws()
  props <- data.frame(
    sample_id = c("S1", "S2"),
    cell_type = c("b memory", "cd14 mono"),
    proportion = c(0.2, 0.15)
  )
  expect_equal(
    composition_test_cohort(props, draws),
    composition_test(props, draws)
  )
})
