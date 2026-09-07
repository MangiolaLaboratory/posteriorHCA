library(testthat)

skip_if_not_installed("edgeR")
skip_if_not_installed("cli")

suppressPackageStartupMessages({
  library(cli)
  if (requireNamespace("readr", quietly = TRUE)) {
    library(readr)
  }
  if (requireNamespace("httr", quietly = TRUE)) {
    library(httr)
  }
})

pkg_root <- testthat::test_path("..", "..")
if (!dir.exists(file.path(pkg_root, "R"))) {
  pkg_root <- "."
}
sys.source(file.path(pkg_root, "R", "cohort.R"), envir = environment())
gene_id_file <- file.path(pkg_root, "R", "gene-id.R")
if (file.exists(gene_id_file)) {
  sys.source(gene_id_file, envir = environment())
}
utils_file <- file.path(pkg_root, "R", "utlis.R")
if (file.exists(utils_file)) {
  sys.source(utils_file, envir = environment())
}

cd8_naive_rds_path <- function() {
  path <- file.path(
    tools::R_user_dir("posteriorHCA", "cache"),
    "V1",
    "reference_samples",
    "cd8.naive.rds"
  )
  if (isTRUE(file.exists(path) && isTRUE(file.size(path) > 0))) {
    path
  } else {
    NULL
  }
}

toy_counts <- function() {
  genes <- paste0("g", 1:20)
  set.seed(1)
  user <- matrix(
    rnbinom(20 * 4, mu = 50, size = 10),
    nrow = 20,
    dimnames = list(genes, paste0("s", 1:4))
  )
  ref <- setNames(as.numeric(rnbinom(20, mu = 200, size = 10)), genes)
  list(user = user, ref = ref)
}

hca_offset <- function(x) {
  if (inherits(x, "Seurat") || inherits(x, "SummarizedExperiment")) {
    meta <- extract_sample_metadata(x)
    return(setNames(meta$hca_offset, rownames(meta)))
  }
  attr(x, "hca_offset")
}

hca_multiplier <- function(x) {
  if (inherits(x, "Seurat") || inherits(x, "SummarizedExperiment")) {
    meta <- extract_sample_metadata(x)
    return(setNames(meta$hca_multiplier, rownames(meta)))
  }
  attr(x, "hca_multiplier")
}

sample_role <- function(x) {
  if (inherits(x, "Seurat") || inherits(x, "SummarizedExperiment")) {
    meta <- extract_sample_metadata(x)
    return(setNames(as.character(meta$sample_role), rownames(meta)))
  }
  attr(x, "sample_role")
}

reference_name_of <- function(x) {
  if (inherits(x, "Seurat") || inherits(x, "SummarizedExperiment")) {
    meta <- extract_sample_metadata(x)
    return(meta$hca_reference_name[[1]])
  }
  attr(x, "reference_name")
}

test_that("merge_with_reference_sample appends the reference column", {
  toy <- toy_counts()
  combined <- merge_with_reference_sample(
    toy$user,
    toy$ref,
    reference_name = "hca_ref"
  )
  expect_true(is.matrix(combined))
  expect_equal(ncol(combined), 5L)
  expect_equal(colnames(combined)[[5]], "hca_ref")
  expect_equal(nrow(combined), 20L)
  expect_equal(attr(combined, "reference_name"), "hca_ref")
  expect_setequal(attr(combined, "shared_features"), paste0("g", 1:20))
  expect_equal(unname(combined[, "hca_ref"]), unname(toy$ref[rownames(combined)]))
})

test_that("calculate_tmm_offset puts the reference at offset 0", {
  toy <- toy_counts()
  combined <- merge_with_reference_sample(toy$user, toy$ref, reference_name = "hca_ref")
  scaling <- calculate_tmm_offset(combined, reference_name = "hca_ref")
  expect_equal(unname(scaling$offset[["hca_ref"]]), 0)
  expect_equal(unname(scaling$multiplier[["hca_ref"]]), 1)
  expect_equal(names(scaling$offset), colnames(combined))
  expect_true(all(is.finite(scaling$offset)))
})

test_that("scale_to_hca_reference matches the core helpers on a matrix", {
  toy <- toy_counts()
  combined <- merge_with_reference_sample(toy$user, toy$ref, reference_name = "hca_ref")
  scaling <- calculate_tmm_offset(combined, reference_name = "hca_ref")
  aligned <- suppressMessages(
    scale_to_hca_reference(toy$user, toy$ref, reference_name = "hca_ref")
  )
  expect_equal(as.matrix(aligned), combined, ignore_attr = TRUE)
  expect_equal(unname(hca_offset(aligned)), unname(scaling$offset))
  expect_equal(unname(hca_multiplier(aligned)), unname(scaling$multiplier))
})

test_that("scale_to_hca_reference puts the reference at offset 0", {
  toy <- toy_counts()
  aligned <- suppressMessages(scale_to_hca_reference(toy$user, toy$ref))
  expect_true(is.matrix(aligned))
  expect_equal(reference_name_of(aligned), "hca_reference")
  expect_equal(unname(hca_offset(aligned)[["hca_reference"]]), 0)
  expect_equal(unname(hca_multiplier(aligned)[["hca_reference"]]), 1)
  expect_equal(ncol(aligned), 5L)
  expect_equal(sample_role(aligned)[["hca_reference"]], "reference")
  expect_equal(unique(sample_role(aligned)[colnames(toy$user)]), "user")
})

test_that("scale_to_hca_reference keeps only shared genes", {
  toy <- toy_counts()
  extra_user <- rbind(toy$user, extra = c(1, 2, 3, 4))
  ref <- c(toy$ref, other = 9)
  aligned <- suppressMessages(scale_to_hca_reference(extra_user, ref))
  expect_setequal(attr(aligned, "shared_features"), paste0("g", 1:20))
  expect_equal(nrow(aligned), 20L)
})

test_that("scale_to_hca_reference errors on a name clash", {
  toy <- toy_counts()
  expect_error(
    scale_to_hca_reference(toy$user, toy$ref, reference_name = "s1"),
    "already in `counts`"
  )
})

test_that("scale_to_hca_reference errors when too few genes are shared", {
  counts <- matrix(1:4, nrow = 2, dimnames = list(c("a", "b"), c("s1", "s2")))
  ref <- c(z = 10)
  expect_error(scale_to_hca_reference(counts, ref), "shared genes")
})

test_that("a legacy sample-id RDS cannot be merged", {
  expect_error(
    parse_reference_object(
      "e7cfc2caa50e7dffa8b15e540da31358___b memory",
      cell_type = "b.memory"
    ),
    "sample id"
  )
})

test_that("parse_reference_object reads a named count vector", {
  parsed <- parse_reference_object(c(g1 = 10, g2 = 20))
  expect_equal(parsed$counts[["g1"]], 10)
  expect_equal(parsed$sample_id, "hca_reference")
})

test_that("a Nectar download list with a one-column RDS is merged", {
  toy <- toy_counts()
  sid <- "e7cfc2caa50e7dffa8b15e540da31358___b memory"
  ref_mat <- matrix(toy$ref, ncol = 1, dimnames = list(names(toy$ref), sid))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(ref_mat, path)

  download <- list(status = "success", path = path, cell_type = "b.memory")
  aligned <- suppressMessages(scale_to_hca_reference(toy$user, download))
  expect_equal(reference_name_of(aligned), sid)
  expect_equal(ncol(aligned), 5L)
  expect_equal(unname(hca_offset(aligned)[[sid]]), 0)
  expect_equal(sample_role(aligned)[[sid]], "reference")
  expect_equal(sum(sample_role(aligned) == "user"), 4L)
})

test_that("a sample-id-only reference errors", {
  toy <- toy_counts()
  parsed <- list(
    sample_id = "e7cfc2caa50e7dffa8b15e540da31358___b memory",
    counts = NULL
  )
  expect_error(
    scale_to_hca_reference(toy$user, parsed),
    "no count vector"
  )
})

test_that("SummarizedExperiment input returns SummarizedExperiment", {
  skip_if_not_installed("SummarizedExperiment")
  toy <- toy_counts()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = toy$user),
    colData = S4Vectors::DataFrame(condition = c("A", "A", "B", "B"))
  )
  aligned <- suppressMessages(scale_to_hca_reference(se, toy$ref))
  expect_s4_class(aligned, "SummarizedExperiment")
  expect_equal(ncol(aligned), 5L)
  expect_equal(unname(hca_offset(aligned)[["hca_reference"]]), 0)
})

test_that("Seurat input returns Seurat", {
  skip_if_not_installed("Seurat")
  toy <- toy_counts()
  so <- suppressWarnings(Seurat::CreateSeuratObject(counts = toy$user))
  aligned <- suppressMessages(scale_to_hca_reference(so, toy$ref))
  expect_true(inherits(aligned, "Seurat"))
  expect_equal(unname(hca_offset(aligned)[["hca_reference"]]), 0)
  expect_equal(sum(sample_role(aligned) == "user"), ncol(toy$user))
})

test_that("Seurat colnames are preserved after scale_to_hca_reference", {
  skip_if_not_installed("Seurat")
  genes <- paste0("g", 1:10)
  sample_names <- c(
    "C1_17. Disease-associated monocytes",
    "P1-STING-ht_17. Disease-associated monocytes"
  )
  mat <- matrix(1:20, nrow = 10, dimnames = list(genes, sample_names))
  so <- suppressWarnings(Seurat::CreateSeuratObject(counts = mat))
  so$Category <- factor(c("CTRL", "SAVI"), levels = c("CTRL", "SAVI", "SAVI_treated"))
  ref <- setNames(rep(100, 10), genes)
  aligned <- suppressMessages(
    scale_to_hca_reference(so, ref, reference_name = "hca_reference")
  )
  expect_equal(colnames(aligned), c(sample_names, "hca_reference"))
  expect_s3_class(aligned$Category, "factor")
  expect_equal(
    as.character(aligned$Category[colnames(aligned) %in% sample_names]),
    c("CTRL", "SAVI")
  )
})

test_that("estimate_logmu_ql returns one row per gene x group on user libraries", {
  toy <- toy_counts()
  combined <- merge_with_reference_sample(toy$user, toy$ref, reference_name = "hca_ref")
  scaling <- calculate_tmm_offset(combined, reference_name = "hca_ref")
  user_offset <- scaling$offset[colnames(toy$user)]

  sample_metadata <- data.frame(
    cohort = factor(c("A", "A", "B", "B")),
    row.names = colnames(toy$user),
    stringsAsFactors = FALSE
  )
  design_matrix <- model.matrix(~ 0 + cohort, data = sample_metadata)
  colnames(design_matrix) <- sub("^cohort", "", colnames(design_matrix))

  expression_estimates <- estimate_logmu_ql(toy$user, user_offset, design_matrix)
  expect_equal(nrow(expression_estimates), nrow(toy$user) * ncol(design_matrix))
  expect_setequal(expression_estimates$group, c("A", "B"))
  expect_true(all(c("gene", "group", "n", "log_mu", "mu", "se", "df", "dispersion") %in% names(expression_estimates)))
  expect_false("cell_type" %in% names(expression_estimates))
  expect_true(all(is.finite(expression_estimates$log_mu)))
  expect_true(all(expression_estimates$se > 0))
})

test_that("estimate_logmu_ql errors when design rows do not match samples", {
  toy <- toy_counts()
  design_matrix <- model.matrix(~ 0 + factor(c("A", "B")))
  expect_error(
    estimate_logmu_ql(toy$user, rep(0, 4), design_matrix),
    "must match"
  )
})

test_that("a gene-specific count increase raises log_mu", {
  genes <- paste0("g", 1:30)
  set.seed(2)
  counts <- matrix(
    rnbinom(30 * 6, mu = 40, size = 8),
    nrow = 30,
    dimnames = list(genes, paste0("s", 1:6))
  )
  counts["g1", 1:3] <- rnbinom(3, mu = 8, size = 8)
  counts["g1", 4:6] <- rnbinom(3, mu = 120, size = 8)
  ref <- setNames(as.numeric(rnbinom(30, mu = 40, size = 8)), genes)

  combined <- merge_with_reference_sample(counts, ref, reference_name = "hca_ref")
  scaling <- calculate_tmm_offset(combined, reference_name = "hca_ref")
  user_offset <- scaling$offset[colnames(counts)]

  sample_metadata <- data.frame(
    cohort = factor(c(rep("low", 3), rep("high", 3))),
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )
  design_matrix <- model.matrix(~ 0 + cohort, data = sample_metadata)
  colnames(design_matrix) <- sub("^cohort", "", colnames(design_matrix))

  expression_estimates <- estimate_logmu_ql(counts, user_offset, design_matrix)
  expression_estimates <- expression_estimates[
    expression_estimates$gene == "g1",
    ,
    drop = FALSE
  ]
  expect_gt(
    expression_estimates$log_mu[expression_estimates$group == "high"],
    expression_estimates$log_mu[expression_estimates$group == "low"]
  )
})

test_that("cd8.naive Nectar SE is parsed, merged, and scaled with user counts", {
  skip_if_not_installed("SummarizedExperiment")
  path <- cd8_naive_rds_path()
  if (is.null(path)) {
    skip("Cached cd8.naive.rds is not available")
  }

  ref <- readRDS(path)
  expect_s4_class(ref, "SummarizedExperiment")
  expect_equal(ncol(ref), 1L)

  sid <- colnames(ref)[[1]]
  parsed <- parse_reference_object(ref, cell_type = "cd8.naive", path = path)
  expect_equal(parsed$sample_id, sid)
  expect_equal(length(parsed$counts), nrow(ref))

  genes <- rownames(ref)
  set.seed(3)
  n_user <- 3L
  user <- matrix(
    rnbinom(length(genes) * n_user, mu = 40, size = 8),
    nrow = length(genes),
    dimnames = list(genes, paste0("user", seq_len(n_user)))
  )

  aligned_se <- suppressMessages(scale_to_hca_reference(user, ref))
  expect_equal(reference_name_of(aligned_se), sid)
  expect_equal(unname(hca_offset(aligned_se)[[sid]]), 0)
  expect_equal(ncol(aligned_se), n_user + 1L)
  expect_equal(sum(sample_role(aligned_se) == "user"), n_user)

  download <- list(status = "success", path = path, cell_type = "cd8.naive")
  aligned_dl <- suppressMessages(scale_to_hca_reference(user, download))
  expect_equal(reference_name_of(aligned_dl), sid)
  expect_equal(unname(hca_offset(aligned_dl)[[sid]]), 0)
})

test_that("reference merge uses counts, not counts_scaled", {
  skip_if_not_installed("SummarizedExperiment")
  path <- cd8_naive_rds_path()
  if (is.null(path)) {
    skip("Cached cd8.naive.rds is not available")
  }

  ref <- readRDS(path)
  SummarizedExperiment::assay(ref, "counts_scaled") <-
    SummarizedExperiment::assay(ref, "counts_scaled") * 10
  parsed <- parse_reference_object(ref)
  expect_equal(
    unname(parsed$counts),
    unname(as.numeric(SummarizedExperiment::assay(ref, "counts")[, 1L]))
  )
})
