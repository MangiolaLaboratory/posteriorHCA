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

scaling_field <- function(x, field) {
  if (inherits(x, "Seurat") || inherits(x, "SummarizedExperiment")) {
    meta <- extract_sample_metadata(x)
    return(setNames(meta[[field]], rownames(meta)))
  }
  attr(x, field)
}

hca_multiplier <- function(x) {
  scaling_field(x, "hca_multiplier")
}

log_effective_library_size <- function(x) {
  scaling_field(x, "log_effective_library_size")
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

hca_log_E_of <- function(x) {
  if (inherits(x, "Seurat") || inherits(x, "SummarizedExperiment")) {
    meta <- extract_sample_metadata(x)
    return(meta$hca_log_effective_library_size[[1]])
  }
  attr(x, "hca_log_effective_library_size")
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

test_that("calculate_tmm_scaling retains log(E) and HCA log(E_H)", {
  toy <- toy_counts()
  combined <- merge_with_reference_sample(toy$user, toy$ref, reference_name = "hca_ref")
  scaling <- calculate_tmm_scaling(combined, reference_name = "hca_ref")

  expect_equal(unname(scaling$multiplier[["hca_ref"]]), 1)
  expect_equal(
    unname(scaling$hca_effective_library_size),
    unname(scaling$effective_size[["hca_ref"]])
  )
  expect_equal(
    unname(scaling$hca_log_effective_library_size),
    unname(log(scaling$effective_size[["hca_ref"]]))
  )
  expect_equal(
    unname(scaling$log_effective_library_size),
    unname(log(scaling$effective_size))
  )
  expect_false(isTRUE(all.equal(unname(scaling$hca_log_effective_library_size), 0)))
  expect_null(scaling$offset)
  expect_equal(names(scaling$effective_size), colnames(combined))
  expect_true(all(is.finite(scaling$log_effective_library_size)))
})

test_that("scale_to_hca_reference matches the core helpers on a matrix", {
  toy <- toy_counts()
  combined <- merge_with_reference_sample(toy$user, toy$ref, reference_name = "hca_ref")
  scaling <- calculate_tmm_scaling(combined, reference_name = "hca_ref")
  aligned <- suppressMessages(
    scale_to_hca_reference(toy$user, toy$ref, reference_name = "hca_ref")
  )
  expect_equal(as.matrix(aligned), combined, ignore_attr = TRUE)
  expect_equal(
    unname(log_effective_library_size(aligned)),
    unname(scaling$log_effective_library_size)
  )
  expect_equal(unname(hca_multiplier(aligned)), unname(scaling$multiplier))
  expect_equal(hca_log_E_of(aligned), scaling$hca_log_effective_library_size)
})

test_that("scale_to_hca_reference keeps reference log(E_H), not offset 0", {
  toy <- toy_counts()
  aligned <- suppressMessages(scale_to_hca_reference(toy$user, toy$ref))
  expect_true(is.matrix(aligned))
  expect_equal(reference_name_of(aligned), "hca_reference")
  expect_equal(unname(hca_multiplier(aligned)[["hca_reference"]]), 1)
  expect_equal(
    unname(log_effective_library_size(aligned)[["hca_reference"]]),
    hca_log_E_of(aligned)
  )
  expect_false(isTRUE(all.equal(hca_log_E_of(aligned), 0)))
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
  expect_equal(
    unname(log_effective_library_size(aligned)[[sid]]),
    hca_log_E_of(aligned)
  )
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
  expect_equal(
    unname(log_effective_library_size(aligned)[["hca_reference"]]),
    hca_log_E_of(aligned)
  )
  expect_true("log_effective_library_size" %in% names(SummarizedExperiment::colData(aligned)))
  expect_false("hca_offset" %in% names(SummarizedExperiment::colData(aligned)))
})

test_that("Seurat input returns Seurat", {
  skip_if_not_installed("Seurat")
  toy <- toy_counts()
  so <- suppressWarnings(Seurat::CreateSeuratObject(counts = toy$user))
  aligned <- suppressMessages(scale_to_hca_reference(so, toy$ref))
  expect_true(inherits(aligned, "Seurat"))
  expect_equal(
    unname(log_effective_library_size(aligned)[["hca_reference"]]),
    hca_log_E_of(aligned)
  )
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

test_that("estimate_ql returns coefficients for a one-hot group design", {
  toy <- toy_counts()
  combined <- merge_with_reference_sample(toy$user, toy$ref, reference_name = "hca_ref")
  scaling <- calculate_tmm_scaling(combined, reference_name = "hca_ref")
  user_offset <- scaling$log_effective_library_size[colnames(toy$user)]

  sample_metadata <- data.frame(
    Category = factor(c("A", "A", "B", "B")),
    row.names = colnames(toy$user),
    stringsAsFactors = FALSE
  )

  expression_estimates <- suppressMessages(estimate_ql(
    toy$user,
    user_offset,
    metadata = sample_metadata,
    formula = ~ 0 + Category,
    contrast = NULL
  ))
  expect_equal(nrow(expression_estimates), nrow(toy$user) * 2L)
  expect_setequal(expression_estimates$contrast, c("CategoryA", "CategoryB"))
  expect_true(all(c("gene", "contrast", "estimate", "se", "df", "dispersion") %in% names(expression_estimates)))
  expect_false("log_mu" %in% names(expression_estimates))
  expect_false("group" %in% names(expression_estimates))
  expect_true(all(is.finite(expression_estimates$estimate)))
  expect_true(all(is.finite(expression_estimates$se) & expression_estimates$se > 0))
})

test_that("edgeR receives log(E_user), not HCA-centred offsets", {
  toy <- toy_counts()
  combined <- merge_with_reference_sample(toy$user, toy$ref, reference_name = "hca_ref")
  scaling <- calculate_tmm_scaling(combined, reference_name = "hca_ref")
  user_offset <- scaling$log_effective_library_size[colnames(toy$user)]
  design_matrix <- model.matrix(
    ~ 0 + factor(c("A", "A", "B", "B")),
    data = data.frame(row.names = colnames(toy$user))
  )

  ql <- fit_nb_ql(toy$user, offset = user_offset, design = design_matrix)
  expect_equal(colnames(toy$user), setdiff(colnames(combined), "hca_ref"))
  expect_false("hca_ref" %in% colnames(toy$user))
  expect_equal(dim(ql$offset), dim(toy$user))
  expect_equal(
    unname(ql$offset[1, ]),
    unname(as.numeric(user_offset))
  )
  expect_equal(
    unname(ql$offset[1, ]),
    unname(log(scaling$effective_size[colnames(toy$user)]))
  )
  centred <- log(
    scaling$effective_size[colnames(toy$user)] /
      scaling$hca_effective_library_size
  )
  expect_false(isTRUE(all.equal(unname(ql$offset[1, ]), unname(as.numeric(centred)))))
  expect_true(all(abs(ql$offset - matrix(ql$offset[1, ], nrow = nrow(ql$offset), ncol = ncol(ql$offset), byrow = TRUE)) < 1e-12))
})

test_that("old centred-offset and new log(E) parameterisations match", {
  toy <- toy_counts()
  combined <- merge_with_reference_sample(toy$user, toy$ref, reference_name = "hca_ref")
  scaling <- calculate_tmm_scaling(combined, reference_name = "hca_ref")
  user_samples <- colnames(toy$user)
  E_user <- scaling$effective_size[user_samples]
  E_H <- scaling$hca_effective_library_size
  log_E_H <- scaling$hca_log_effective_library_size

  sample_metadata <- data.frame(
    Category = factor(c("A", "A", "B", "B")),
    row.names = user_samples,
    stringsAsFactors = FALSE
  )
  design_matrix <- model.matrix(~ 0 + Category, data = sample_metadata)

  offset_old <- log(E_user / E_H)
  offset_new <- log(E_user)

  fit_old <- fit_nb_ql(toy$user, offset = offset_old, design = design_matrix)
  fit_new <- fit_nb_ql(toy$user, offset = offset_new, design = design_matrix)

  expect_equal(fit_old$fit$fitted.values, fit_new$fit$fitted.values, tolerance = 1e-6)
  expect_equal(
    fit_old$fit$coefficients,
    fit_new$fit$coefficients + log_E_H,
    tolerance = 1e-6
  )

  est_new <- suppressMessages(estimate_ql(
    toy$user,
    offset_new,
    metadata = sample_metadata,
    formula = ~ 0 + Category
  ))
  log_mu_hca <- est_new$estimate + log_E_H
  log_mu_old <- as.numeric(fit_old$fit$coefficients)
  expect_equal(log_mu_hca, log_mu_old, tolerance = 1e-6)
  expect_equal(exp(log_mu_hca), exp(log_mu_old), tolerance = 1e-6)

  coef_old <- fit_old$fit$coefficients
  coef_new <- fit_new$fit$coefficients
  expect_equal(
    coef_old[, "CategoryA"] - coef_old[, "CategoryB"],
    coef_new[, "CategoryA"] - coef_new[, "CategoryB"],
    tolerance = 1e-6
  )
})

test_that("HCA-scale estimates are invariant to global effective-size rescaling", {
  toy <- toy_counts()
  combined <- merge_with_reference_sample(toy$user, toy$ref, reference_name = "hca_ref")
  scaling <- calculate_tmm_scaling(combined, reference_name = "hca_ref")
  user_samples <- colnames(toy$user)

  sample_metadata <- data.frame(
    Category = factor(c("A", "A", "B", "B")),
    row.names = user_samples,
    stringsAsFactors = FALSE
  )

  c_factor <- 3.5
  offset_base <- log(scaling$effective_size[user_samples])
  offset_scaled <- offset_base + log(c_factor)
  log_E_H <- scaling$hca_log_effective_library_size
  log_E_H_scaled <- log_E_H + log(c_factor)

  est_base <- suppressMessages(estimate_ql(
    toy$user, offset_base, metadata = sample_metadata, formula = ~ 0 + Category
  ))
  est_scaled <- suppressMessages(estimate_ql(
    toy$user, offset_scaled, metadata = sample_metadata, formula = ~ 0 + Category
  ))
  expect_equal(
    est_base$estimate + log_E_H,
    est_scaled$estimate + log_E_H_scaled,
    tolerance = 1e-6
  )
})

test_that("formula+contrast estimate equals c'beta and SE equals sqrt(c'Vc)", {
  toy <- toy_counts()
  combined <- merge_with_reference_sample(toy$user, toy$ref, reference_name = "hca_ref")
  scaling <- calculate_tmm_scaling(combined, reference_name = "hca_ref")
  user_offset <- scaling$log_effective_library_size[colnames(toy$user)]
  meta <- data.frame(
    Category = factor(c("A", "A", "B", "B"), levels = c("A", "B")),
    row.names = colnames(toy$user)
  )

  est <- suppressMessages(estimate_ql(
    toy$user,
    user_offset,
    metadata = meta,
    formula = ~ 0 + Category,
    contrast = "CategoryA"
  ))
  expect_equal(unique(est$contrast), "CategoryA")
  expect_false("log_mu" %in% names(est))

  fit <- attr(est, "fit")
  design <- attr(est, "design")
  C <- attr(est, "contrast")
  g <- 1L
  V <- vcov_ql_gene(fit, design, g)
  expect_equal(
    est$estimate[est$gene == rownames(toy$user)[[g]]],
    as.numeric(crossprod(C[, 1], fit$coefficients[g, ])),
    tolerance = 1e-10
  )
  expect_equal(
    est$se[est$gene == rownames(toy$user)[[g]]],
    sqrt(as.numeric(crossprod(C[, 1], V %*% C[, 1]))),
    tolerance = 1e-10
  )
})

test_that("Category + Experiment estimands are c'beta with SE from same fit", {
  toy <- toy_counts()
  combined <- merge_with_reference_sample(toy$user, toy$ref, reference_name = "hca_ref")
  scaling <- calculate_tmm_scaling(combined, reference_name = "hca_ref")
  user_offset <- scaling$log_effective_library_size[colnames(toy$user)]
  meta <- data.frame(
    Category = factor(c("Control", "Control", "SAVI", "SAVI"), levels = c("Control", "SAVI")),
    Experiment = factor(c("E1", "E2", "E1", "E2")),
    row.names = colnames(toy$user)
  )

  est <- suppressMessages(estimate_ql(
    toy$user,
    user_offset,
    metadata = meta,
    formula = ~ 0 + Category + Experiment,
    contrast = c("CategorySAVI", "CategorySAVI + ExperimentE2")
  ))
  expect_true(all(is.finite(est$estimate) & is.finite(est$se) & est$se > 0))
  expect_false("log_mu" %in% names(est))

  fit <- attr(est, "fit")
  C <- attr(est, "contrast")
  expect_equal(est$estimate, as.numeric(fit$coefficients %*% C), tolerance = 1e-10)

  # Relative contrast: no automatic HCA shift in core output.
  rel <- suppressMessages(estimate_ql(
    toy$user,
    user_offset,
    metadata = meta,
    formula = ~ 0 + Category + Experiment,
    contrast = "CategorySAVI - CategoryControl"
  ))
  expect_equal(
    rel$estimate,
    as.numeric(fit$coefficients %*% attr(rel, "contrast")),
    tolerance = 1e-8
  )
  expect_false(isTRUE(all.equal(
    rel$estimate,
    rel$estimate + scaling$hca_log_effective_library_size
  )))
})

test_that("general SE matches one-hot closed form for ~ 0 + Category and ~ 1", {
  toy <- toy_counts()
  combined <- merge_with_reference_sample(toy$user, toy$ref, reference_name = "hca_ref")
  scaling <- calculate_tmm_scaling(combined, reference_name = "hca_ref")
  user_offset <- scaling$log_effective_library_size[colnames(toy$user)]
  meta <- data.frame(
    Category = factor(c("A", "A", "B", "B")),
    row.names = colnames(toy$user)
  )

  est <- suppressMessages(estimate_ql(
    toy$user, user_offset, metadata = meta, formula = ~ 0 + Category, contrast = "CategoryA"
  ))
  design <- attr(est, "design")
  fit <- attr(est, "fit")
  phi <- ql_fit_nb_dispersion(fit, nrow(toy$user))
  se_closed <- se_group_mean_closed_form(fit, design, phi)[, "CategoryA"]
  expect_equal(est$se, se_closed, tolerance = 1e-10, ignore_attr = TRUE)

  meta1 <- data.frame(row.names = colnames(toy$user))
  est1 <- suppressMessages(estimate_ql(
    toy$user, user_offset, metadata = meta1, formula = ~ 1, contrast = "(Intercept)"
  ))
  design1 <- attr(est1, "design")
  fit1 <- attr(est1, "fit")
  phi1 <- ql_fit_nb_dispersion(fit1, nrow(toy$user))
  expect_equal(
    est1$se,
    se_group_mean_closed_form(fit1, design1, phi1)[, 1],
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
})

test_that("algebraically equivalent contrasts agree (Bioconductor-style)", {
  # (B - A) - (C - A) == B - C as linear algebra on coefficients.
  set.seed(21)
  genes <- paste0("g", 1:25)
  counts <- matrix(
    rnbinom(25 * 9, mu = 45, size = 6),
    nrow = 25,
    dimnames = list(genes, paste0("s", 1:9))
  )
  meta <- data.frame(
    Category = factor(
      rep(c("A", "B", "C"), each = 3),
      levels = c("A", "B", "C")
    ),
    row.names = colnames(counts)
  )
  offset <- log(colSums(counts))

  est_bc <- suppressMessages(estimate_ql(
    counts, offset, metadata = meta, formula = ~ 0 + Category,
    contrast = "CategoryB - CategoryC"
  ))
  est_equiv <- suppressMessages(estimate_ql(
    counts, offset, metadata = meta, formula = ~ 0 + Category,
    contrast = "(CategoryB - CategoryA) - (CategoryC - CategoryA)"
  ))

  expect_equal(est_bc$estimate, est_equiv$estimate, tolerance = 1e-10)
  expect_equal(est_bc$se, est_equiv$se, tolerance = 1e-10)
  # Contrast vectors themselves reduce to the same coefficients.
  expect_equal(
    as.numeric(attr(est_bc, "contrast")),
    as.numeric(attr(est_equiv, "contrast")),
    tolerance = 1e-12
  )
})

test_that("contrast sign reversal flips estimate and preserves SE", {
  set.seed(22)
  genes <- paste0("g", 1:20)
  counts <- matrix(
    rnbinom(20 * 6, mu = 35, size = 7),
    nrow = 20,
    dimnames = list(genes, paste0("s", 1:6))
  )
  meta <- data.frame(
    Category = factor(rep(c("B", "C"), each = 3), levels = c("B", "C")),
    row.names = colnames(counts)
  )
  offset <- log(colSums(counts))

  est_fwd <- suppressMessages(estimate_ql(
    counts, offset, metadata = meta, formula = ~ 0 + Category,
    contrast = "CategoryB - CategoryC"
  ))
  est_rev <- suppressMessages(estimate_ql(
    counts, offset, metadata = meta, formula = ~ 0 + Category,
    contrast = "CategoryC - CategoryB"
  ))

  expect_equal(est_rev$estimate, -est_fwd$estimate, tolerance = 1e-10)
  expect_equal(est_rev$se, est_fwd$se, tolerance = 1e-10)
})

test_that("glmQLFTest logFC matches natural-log estimate / log(2)", {
  # edgeR DE tables report log2 fold-change:
  #   logFC = (coefficients %*% contrast) / log(2)
  # Our estimate is the natural-log contrast c' beta.
  # glmQLFTest F is LR / df.test / s2.post (LR-based), not the Wald
  # (estimate / SE)^2; we only assert the estimate/logFC link here.
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")

  set.seed(23)
  genes <- paste0("g", 1:40)
  counts <- matrix(
    rnbinom(40 * 8, mu = 40, size = 5),
    nrow = 40,
    dimnames = list(genes, paste0("s", 1:8))
  )
  meta <- data.frame(
    Category = factor(
      rep(c("Control", "SAVI"), each = 4),
      levels = c("Control", "SAVI")
    ),
    row.names = colnames(counts)
  )
  offset <- log(colSums(counts))

  est <- suppressMessages(estimate_ql(
    counts, offset, metadata = meta, formula = ~ 0 + Category,
    contrast = "CategorySAVI - CategoryControl"
  ))
  fit <- attr(est, "fit")
  C <- attr(est, "contrast")
  qlf <- edgeR::glmQLFTest(fit, contrast = C)

  expect_equal(
    est$estimate,
    qlf$table$logFC * log(2),
    tolerance = 1e-10,
    ignore_attr = TRUE
  )

  # Sanity: Wald and LR-based F are related but not identical definitions.
  se <- est$se
  t2 <- (est$estimate / se)^2
  expect_true(cor(t2, qlf$table$F) > 0.99)
  expect_false(isTRUE(all.equal(t2, qlf$table$F, tolerance = 1e-8)))
})

test_that("~ 1 + Experiment contrasts share one fit and return finite SE", {
  set.seed(24)
  genes <- paste0("g", 1:20)
  counts <- matrix(
    rnbinom(20 * 6, mu = 30, size = 8),
    nrow = 20,
    dimnames = list(genes, paste0("s", 1:6))
  )
  meta <- data.frame(
    Experiment = factor(c("E1", "E1", "E2", "E2", "E3", "E3")),
    row.names = colnames(counts)
  )
  offset <- log(colSums(counts))

  est <- suppressMessages(estimate_ql(
    counts,
    offset,
    metadata = meta,
    formula = ~ 1 + Experiment,
    contrast = c("(Intercept)", "(Intercept) + ExperimentE2")
  ))
  est_int <- est[est$contrast == "(Intercept)", , drop = FALSE]
  est_e2 <- est[est$contrast == "(Intercept) + ExperimentE2", , drop = FALSE]

  expect_true(all(is.finite(est_int$estimate) & is.finite(est_int$se) & est_int$se > 0))
  expect_true(all(is.finite(est_e2$estimate) & is.finite(est_e2$se) & est_e2$se > 0))

  fit <- attr(est, "fit")
  C <- attr(est, "contrast")
  expect_equal(
    est$estimate,
    as.numeric(fit$coefficients %*% C),
    tolerance = 1e-10
  )
})

test_that("~ 1 + Experiment coefficient SEs differ from naive shortcut", {
  genes <- paste0("g", 1:30)
  set.seed(11)
  counts <- matrix(
    rnbinom(30 * 6, mu = 35, size = 8),
    nrow = 30,
    dimnames = list(genes, paste0("s", 1:6))
  )
  meta <- data.frame(
    Experiment = factor(c("E1", "E1", "E2", "E2", "E3", "E3")),
    row.names = colnames(counts)
  )
  offset <- log(colSums(counts))

  est <- suppressMessages(estimate_ql(
    counts, offset, metadata = meta, formula = ~ 1 + Experiment
  ))
  expect_true(all(c("(Intercept)", "ExperimentE2", "ExperimentE3") %in% est$contrast))
  expect_true(all(is.finite(est$se) & est$se > 0))

  fit <- attr(est, "fit")
  design <- attr(est, "design")
  phi <- ql_fit_nb_dispersion(fit, nrow(counts))
  naive <- se_group_mean_closed_form(fit, design, phi)
  expect_false(isTRUE(all.equal(se_ql_coefficients(fit, design), naive)))
})

test_that("core does not add log(E_H) to arbitrary coefficients", {
  toy <- toy_counts()
  combined <- merge_with_reference_sample(toy$user, toy$ref, reference_name = "hca_ref")
  scaling <- calculate_tmm_scaling(combined, reference_name = "hca_ref")
  user_offset <- scaling$log_effective_library_size[colnames(toy$user)]
  meta <- data.frame(
    Category = factor(c("A", "A", "B", "B")),
    Experiment = factor(c("e1", "e2", "e1", "e2")),
    row.names = colnames(toy$user)
  )
  est <- suppressMessages(estimate_ql(
    toy$user, user_offset, metadata = meta, formula = ~ 0 + Category + Experiment
  ))
  fit <- attr(est, "fit")
  expect_equal(est$estimate, as.numeric(fit$coefficients), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(
    est$estimate,
    as.numeric(fit$coefficients) + scaling$hca_log_effective_library_size
  )))
})

test_that("estimate_ql errors on mismatched metadata or offset", {
  toy <- toy_counts()
  meta <- data.frame(Category = factor(c("A", "B")), row.names = c("x", "y"))
  expect_error(
    estimate_ql(toy$user, rep(0, 4), metadata = meta, formula = ~ 0 + Category),
    "metadata|match"
  )

  meta_ok <- data.frame(
    Category = factor(c("A", "A", "B", "B")),
    row.names = colnames(toy$user)
  )
  expect_error(
    estimate_ql(toy$user, offset = c(0, 0), metadata = meta_ok, formula = ~ 0 + Category),
    "offset"
  )
  expect_error(
    estimate_ql(
      toy$user,
      offset = matrix(0, 2, 2),
      metadata = meta_ok,
      formula = ~ 0 + Category
    ),
    "offset"
  )
})

test_that("welch_test_means accepts generic user-provided mean and SE", {
  out <- posteriorHCA::welch_test_means(
    mu1 = 1.2, se1 = 0.3, mu2 = 0.8, se2 = 0.25, n1 = 5, n2 = 40
  )
  expect_true(is.finite(out$p_value))
  expect_equal(out$mu1, 1.2)
  expect_equal(out$se1, 0.3)
})

test_that("a gene-specific count increase raises the coefficient estimate", {
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
  scaling <- calculate_tmm_scaling(combined, reference_name = "hca_ref")
  user_offset <- scaling$log_effective_library_size[colnames(counts)]

  sample_metadata <- data.frame(
    cohort = factor(c(rep("low", 3), rep("high", 3))),
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )

  expression_estimates <- suppressMessages(estimate_ql(
    counts,
    user_offset,
    metadata = sample_metadata,
    formula = ~ 0 + cohort
  ))
  expression_estimates <- expression_estimates[
    expression_estimates$gene == "g1",
    ,
    drop = FALSE
  ]
  expect_gt(
    expression_estimates$estimate[expression_estimates$contrast == "cohorthigh"],
    expression_estimates$estimate[expression_estimates$contrast == "cohortlow"]
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
  expect_equal(
    unname(log_effective_library_size(aligned_se)[[sid]]),
    hca_log_E_of(aligned_se)
  )
  expect_equal(ncol(aligned_se), n_user + 1L)
  expect_equal(sum(sample_role(aligned_se) == "user"), n_user)

  download <- list(status = "success", path = path, cell_type = "cd8.naive")
  aligned_dl <- suppressMessages(scale_to_hca_reference(user, download))
  expect_equal(reference_name_of(aligned_dl), sid)
  expect_equal(
    unname(log_effective_library_size(aligned_dl)[[sid]]),
    hca_log_E_of(aligned_dl)
  )
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

test_that("estimate_cohort_logmu returns HCA-scale group means for ~ 0 + Category", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  toy <- toy_counts()
  aligned <- suppressMessages(
    scale_to_hca_reference(toy$user, toy$ref, reference_name = "hca_ref")
  )

  sample_metadata <- data.frame(
    sample_role = unname(attr(aligned, "sample_role")),
    log_effective_library_size = unname(
      as.numeric(attr(aligned, "log_effective_library_size"))
    ),
    hca_log_effective_library_size = attr(
      aligned,
      "hca_log_effective_library_size"
    ),
    Category = factor(c("A", "A", "B", "B", NA_character_), levels = c("A", "B")),
    row.names = colnames(aligned),
    stringsAsFactors = FALSE
  )
  scaled_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = aligned),
    colData = S4Vectors::DataFrame(sample_metadata)
  )

  wrapper_result <- estimate_cohort_logmu(
    scaled_se,
    formula = ~ 0 + Category,
    gene_ensg = "g1"
  )

  user_samples <- sample_metadata$sample_role == "user"
  user_counts <- aligned[, user_samples, drop = FALSE]
  user_metadata <- droplevels(sample_metadata[user_samples, , drop = FALSE])
  user_offset <- setNames(
    as.numeric(user_metadata$log_effective_library_size),
    rownames(user_metadata)
  )
  design_matrix <- stats::model.matrix(~ 0 + Category, data = user_metadata)
  core_result <- suppressMessages(estimate_ql(
    user_counts,
    user_offset,
    metadata = user_metadata,
    formula = ~ 0 + Category
  ))
  core_result <- core_result[core_result$gene == "g1", , drop = FALSE]
  hca_log_E <- user_metadata$hca_log_effective_library_size[[1]]

  expect_equal(nrow(wrapper_result), nrow(core_result))
  expect_equal(wrapper_result$group, core_result$contrast)
  expect_equal(wrapper_result$log_mu, core_result$estimate + hca_log_E)
  expect_equal(wrapper_result$se, core_result$se)
  expect_equal(wrapper_result$n, c(2L, 2L))
})

test_that("estimate_cohort_logmu accepts ~ 1 for a single group", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  toy <- toy_counts()
  aligned <- suppressMessages(
    scale_to_hca_reference(toy$user[, 1:2, drop = FALSE], toy$ref, reference_name = "hca_ref")
  )
  sample_metadata <- data.frame(
    sample_role = unname(attr(aligned, "sample_role")),
    log_effective_library_size = unname(
      as.numeric(attr(aligned, "log_effective_library_size"))
    ),
    hca_log_effective_library_size = attr(
      aligned,
      "hca_log_effective_library_size"
    ),
    row.names = colnames(aligned),
    stringsAsFactors = FALSE
  )
  scaled_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = aligned),
    colData = S4Vectors::DataFrame(sample_metadata)
  )
  out <- estimate_cohort_logmu(scaled_se, formula = ~ 1, gene_ensg = "g1")
  expect_equal(unique(out$group), "all")
  expect_true(all(is.finite(out$log_mu)))
})

test_that("estimate_cohort_logmu rejects designs outside the group-mean workflow", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  toy <- toy_counts()
  aligned <- suppressMessages(
    scale_to_hca_reference(toy$user, toy$ref, reference_name = "hca_ref")
  )
  sample_metadata <- data.frame(
    sample_role = unname(attr(aligned, "sample_role")),
    log_effective_library_size = unname(
      as.numeric(attr(aligned, "log_effective_library_size"))
    ),
    hca_log_effective_library_size = attr(
      aligned,
      "hca_log_effective_library_size"
    ),
    Category = factor(c("A", "A", "B", "B", NA_character_), levels = c("A", "B")),
    Experiment = factor(c("e1", "e2", "e1", "e2", NA_character_)),
    batch = factor(c("b1", "b1", "b2", "b2", NA_character_)),
    row.names = colnames(aligned),
    stringsAsFactors = FALSE
  )
  scaled_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = aligned),
    colData = S4Vectors::DataFrame(sample_metadata)
  )

  expect_error(
    estimate_cohort_logmu(scaled_se, formula = ~ Category),
    "convenience workflow"
  )
  expect_error(
    estimate_cohort_logmu(scaled_se, formula = ~ 0 + Category + Experiment),
    "lower-level core API"
  )
  expect_error(
    estimate_cohort_logmu(scaled_se, formula = ~ 1 + batch),
    "convenience workflow"
  )
})
