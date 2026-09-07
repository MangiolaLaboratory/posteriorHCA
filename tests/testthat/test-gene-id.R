library(testthat)

suppressPackageStartupMessages({
  library(cli)
})

pkg_root <- testthat::test_path("..", "..")
if (!dir.exists(file.path(pkg_root, "R"))) {
  pkg_root <- "."
}
sys.source(file.path(pkg_root, "R", "gene-id.R"), envir = environment())

test_that("is_ensembl_gene_id recognises ENSG ids", {
  expect_true(is_ensembl_gene_id("ENSG00000169252"))
  expect_false(is_ensembl_gene_id("ADRB2"))
})

test_that("strip_ensembl_version removes version suffix", {
  expect_equal(strip_ensembl_version("ENSG00000169252.1"), "ENSG00000169252")
})

test_that("expression metadata helpers propagate fit annotations", {
  fit <- new_expression_fit(
    fit = list(data = data.frame(x = 1)),
    cell_type = "monocytic",
    gene_ensg = "ENSG00000169252"
  )
  meta <- expression_metadata(fit)
  expect_equal(meta$cell_type, "monocytic")
  expect_equal(meta$gene_ensg, "ENSG00000169252")
  expect_true(is_expression_fit(fit))
})
