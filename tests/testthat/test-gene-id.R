library(testthat)

suppressPackageStartupMessages({
  library(cli)
})

pkg_root <- testthat::test_path("..", "..")
if (!dir.exists(file.path(pkg_root, "R"))) {
  pkg_root <- "."
}
cohort_file <- file.path(pkg_root, "R", "cohort.R")
if (file.exists(cohort_file)) {
  sys.source(cohort_file, envir = environment())
}
sys.source(file.path(pkg_root, "R", "gene-id.R"), envir = environment())
utils_file <- file.path(pkg_root, "R", "utlis.R")
if (file.exists(utils_file)) {
  sys.source(utils_file, envir = environment())
}

test_that("is_ensembl_gene_id recognises ENSG ids", {
  expect_true(is_ensembl_gene_id("ENSG00000169252"))
  expect_false(is_ensembl_gene_id("ADRB2"))
})

test_that("strip_ensembl_version removes version suffix", {
  expect_equal(strip_ensembl_version("ENSG00000169252.1"), "ENSG00000169252")
})

test_that("resolve_gene keeps ENSG ids", {
  out <- resolve_gene("ENSG00000169252", id_type = "ensembl", strict = FALSE)
  expect_equal(unname(out), "ENSG00000169252")
})

test_that("resolve_gene maps symbols when orgdb is available", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("AnnotationDbi")
  out <- resolve_gene("ADRB2", id_type = "symbol", strict = TRUE)
  expect_equal(unname(out), "ENSG00000169252")
})

test_that("harmonise_gene_ids renames matrix rownames", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("AnnotationDbi")
  mat <- matrix(1:6, nrow = 2, dimnames = list(c("ADRB2", "PTGS2"), c("s1", "s2", "s3")))
  out <- harmonise_gene_ids(mat, id_type = "symbol", strict = TRUE)
  expect_true(is.matrix(out))
  expect_true(all(is_ensembl_gene_id(rownames(out))))
  expect_equal(ncol(out), 3L)
})

test_that("harmonise_gene_ids returns SummarizedExperiment for SE input", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("SummarizedExperiment")
  mat <- matrix(
    1:6,
    nrow = 2,
    dimnames = list(c("ADRB2", "PTGS2"), c("s1", "s2", "s3"))
  )
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
  out <- harmonise_gene_ids(se, id_type = "symbol", strict = TRUE)
  expect_true(inherits(out, "SummarizedExperiment"))
  expect_true(all(is_ensembl_gene_id(rownames(out))))
  expect_equal(colnames(out), colnames(mat))
})

test_that("expr metadata helpers propagate fit annotations", {
  fit <- new_expr_fit(
    fit = list(data = data.frame(x = 1)),
    cell_type = "monocytic",
    gene_ensg = "ENSG00000169252",
    gene_symbol = "ADRB2"
  )
  meta <- expr_metadata(fit)
  expect_equal(meta$cell_type, "monocytic")
  expect_equal(meta$gene_ensg, "ENSG00000169252")
  expect_equal(meta$gene_symbol, "ADRB2")

  draws <- list(
    draws = 1:5,
    cell_type = "monocytic",
    gene_ensg = "ENSG00000169252"
  )
  cohort <- data.frame(
    gene = "ENSG00000169252",
    cell_type = "monocytic",
    group = "A",
    log_mu = 1,
    se = 0.1,
    n = 3
  )
  expect_silent(
    validate_expr_metadata_match(
      cohort_est_metadata(cohort),
      expr_metadata(draws)
    )
  )
})
