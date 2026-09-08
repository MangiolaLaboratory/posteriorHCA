# Internal bootstrap / mglm helpers (not part of the public cohort API).
# Sourced directly so tests do not depend on exports.
library(testthat)

skip_if_not_installed("edgeR")
skip_if_not_installed("cli")

suppressPackageStartupMessages({
  library(cli)
  library(edgeR)
})

pkg_root <- testthat::test_path("..", "..")
if (!dir.exists(file.path(pkg_root, "R"))) {
  pkg_root <- "."
}
sys.source(file.path(pkg_root, "R", "cohort.R"), envir = environment())
sys.source(file.path(pkg_root, "R", "gene-id.R"), envir = environment())
sys.source(file.path(pkg_root, "R", "inference-cohort.R"), envir = environment())
utils_file <- file.path(pkg_root, "R", "utlis.R")
if (file.exists(utils_file)) {
  sys.source(utils_file, envir = environment())
}

test_that("draw_dirichlet_weights sums to n and are positive", {
  w <- draw_dirichlet_weights(5L)
  expect_equal(length(w), 5L)
  expect_equal(sum(w), 5, tolerance = 1e-6)
  expect_true(all(w > 0))
})

test_that("estimate_dispersion_nb returns named dispersions", {
  genes <- paste0("g", 1:20)
  set.seed(3)
  counts <- matrix(
    rnbinom(20 * 6, mu = 30, size = 5),
    nrow = 20,
    dimnames = list(genes, paste0("s", 1:6))
  )
  offset <- rep(0, 6)
  disp <- estimate_dispersion_nb(counts, offset = offset)
  expect_equal(length(disp), 20L)
  expect_equal(names(disp), genes)
  expect_true(all(is.finite(disp) & disp > 0))
})

test_that("bootstrap_logmu_mglm returns n_boot draws", {
  genes <- paste0("g", 1:10)
  set.seed(7)
  counts <- matrix(
    rnbinom(10 * 5, mu = 40, size = 8),
    nrow = 10,
    dimnames = list(genes, paste0("s", 1:5))
  )
  offset <- rep(0, 5)
  draws <- bootstrap_logmu_mglm(
    y = counts["g1", , drop = FALSE],
    offset = offset,
    dispersion = 0.15,
    n_boot = 40L,
    seed = 99
  )
  expect_equal(length(draws), 40L)
  expect_true(all(is.finite(draws)))
})

test_that("bootstrap_logmu_mglm increases when counts increase", {
  genes <- paste0("g", 1:10)
  set.seed(42)
  counts_low <- matrix(
    rnbinom(10 * 5, mu = 10, size = 10),
    nrow = 10,
    dimnames = list(genes, paste0("s", 1:5))
  )
  counts_high <- matrix(
    rnbinom(10 * 5, mu = 500, size = 10),
    nrow = 10,
    dimnames = list(genes, paste0("s", 1:5))
  )
  offset <- rep(0, 5)

  draws_low <- bootstrap_logmu_mglm(
    counts_low["g1", , drop = FALSE],
    offset = offset,
    dispersion = 0.1,
    n_boot = 50L,
    seed = 1
  )
  draws_high <- bootstrap_logmu_mglm(
    counts_high["g1", , drop = FALSE],
    offset = offset,
    dispersion = 0.1,
    n_boot = 50L,
    seed = 1
  )

  expect_gt(median(draws_high), median(draws_low))
})
