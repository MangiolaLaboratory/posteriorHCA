# User counts on the atlas scale
#
#   Core:
#     merge_with_reference_sample() / calculate_tmm_scaling()
#     fit_nb_ql() / estimate_ql() — thin edgeR helper (formula + contrast)
#   Wrappers (opinionated posteriorHCA workflow):
#     scale_to_hca_reference()   — load + merge + TMM
#     estimate_cohort_logmu()    — one-hot group means + HCA-scale log_mu
#
#   edgeR fits with documented offset = log(effective library size).
#   The core estimator returns estimate + se for the user contrast; it does
#   not add log(E_HCA) or interpret biological meaning of the formula.

#' Extract a gene-by-sample count matrix
#' @keywords internal
#' @noRd
extract_counts <- function(x, assay = NULL) {
  if (inherits(x, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      cli::cli_abort("Package Seurat is required for Seurat input.")
    }
    if (is.null(assay)) {
      assay <- Seurat::DefaultAssay(x)
    }
    return(as.matrix(Seurat::GetAssayData(x, assay = assay, layer = "counts")))
  }

  if (inherits(x, "SummarizedExperiment")) {
    assay_names <- SummarizedExperiment::assayNames(x)
    if (is.null(assay)) {
      assay <- if ("counts" %in% assay_names) "counts" else assay_names[[1]]
    }
    return(as.matrix(SummarizedExperiment::assay(x, assay)))
  }

  as.matrix(x)
}

#' Sample metadata from Seurat / SummarizedExperiment, else empty
#' @keywords internal
#' @noRd
extract_sample_metadata <- function(x) {
  if (inherits(x, "Seurat")) {
    return(as.data.frame(x[[]]))
  }
  if (inherits(x, "SummarizedExperiment")) {
    return(as.data.frame(SummarizedExperiment::colData(x)))
  }
  data.frame()
}

#' Coerce a one-library reference to a named gene count vector
#' @keywords internal
#' @noRd
reference_counts_vector <- function(reference) {
  if (is.numeric(reference) && is.null(dim(reference))) {
    if (is.null(names(reference))) {
      cli::cli_abort("`reference` must be a named numeric vector.")
    }
    return(setNames(as.numeric(reference), names(reference)))
  }

  mat <- extract_counts(reference, assay = "counts")
  if (ncol(mat) != 1L) {
    cli::cli_abort("`reference` must be one library; found {ncol(mat)} columns.")
  }
  setNames(as.numeric(mat[, 1]), rownames(mat))
}

#' Parse a downloaded atlas reference RDS into a named count vector
#' @keywords internal
#' @noRd
parse_reference_object <- function(obj, cell_type = NA_character_, path = NULL) {
  if (is.character(obj) && length(obj) == 1L) {
    cli::cli_abort(c(
      "Reference file is a sample id (`{obj}`), not counts.",
      "i" = "Nectar should store one gene-by-1 SummarizedExperiment per cell type."
    ))
  }

  if (is.numeric(obj) && is.null(dim(obj)) && !is.null(names(obj))) {
    return(list(
      cell_type = cell_type,
      sample_id = "hca_reference",
      counts = setNames(as.numeric(obj), names(obj)),
      path = path
    ))
  }

  mat <- extract_counts(obj, assay = "counts")
  if (ncol(mat) != 1L) {
    cli::cli_abort("Atlas reference must be one sample; found {ncol(mat)} columns.")
  }
  sample_id <- colnames(mat)[[1]]
  if (is.null(sample_id) || !nzchar(sample_id)) {
    sample_id <- "hca_reference"
  }
  list(
    cell_type = cell_type,
    sample_id = sample_id,
    counts = setNames(as.numeric(mat[, 1]), rownames(mat)),
    path = path
  )
}

#' Resolve a reference argument to a parsed atlas library
#' @keywords internal
#' @noRd
resolve_reference <- function(reference, version = "latest") {
  if (is.character(reference) && length(reference) == 1L) {
    return(load_reference_sample(reference, version = version))
  }
  if (is.list(reference) && !is.data.frame(reference) &&
      all(c("status", "path") %in% names(reference))) {
    return(load_reference_sample(reference, version = version))
  }
  if (is.list(reference) && !is.data.frame(reference) && "counts" %in% names(reference)) {
    if (is.null(reference$counts) || !length(reference$counts)) {
      cli::cli_abort("Atlas reference has no count vector to merge.")
    }
    if (is.null(reference$sample_id) || !nzchar(as.character(reference$sample_id))) {
      reference$sample_id <- "hca_reference"
    }
    return(reference)
  }
  parse_reference_object(reference)
}

#' Load the atlas reference library for a cell type
#'
#' @inheritParams get_brms_ready
#' @param cell_type Cell type name, or a download list from
#'   [get_reference_sample_ready()].
#' @return A list with `cell_type`, `sample_id`, `counts` (named numeric
#'   vector), and `path`.
#' @export
#' @importFrom cli cli_abort cli_alert_info
load_reference_sample <- function(
  cell_type,
  version = "latest",
  cache_directory = get_default_cache_dir(),
  use_cache = TRUE
) {
  if (is.list(cell_type) && all(c("status", "path") %in% names(cell_type))) {
    res <- cell_type
  } else {
    res <- get_reference_sample_ready(
      cell_type = cell_type,
      version = version,
      cache_directory = cache_directory,
      use_cache = use_cache
    )
  }

  if (!identical(res$status, "success")) {
    extra <- if (!is.null(res$error)) res$error else res$status
    cli_abort("Failed to retrieve reference sample: {extra}")
  }

  parsed <- parse_reference_object(
    readRDS(res$path),
    cell_type = if (!is.null(res$cell_type)) res$cell_type else NA_character_,
    path = res$path
  )
  cli_alert_info(
    "Loaded atlas reference `{parsed$sample_id}` ({length(parsed$counts)} gene{?s})."
  )
  parsed
}

#' Merge user counts with a one-library reference sample
#'
#' Keeps genes shared with the reference and appends the reference as one
#' column. Does not scale or compute offsets.
#'
#' @param counts Gene-by-sample numeric matrix.
#' @param reference Named numeric gene vector, or a one-column matrix.
#' @param reference_name Column name for the appended reference library.
#' @return Gene-by-sample matrix with the reference column last.
#' @seealso [calculate_tmm_scaling()], [scale_to_hca_reference()]
#' @export
#' @importFrom cli cli_abort
merge_with_reference_sample <- function(
  counts,
  reference,
  reference_name = "hca_reference"
) {
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  reference_name <- as.character(reference_name[[1]])
  if (reference_name %in% colnames(counts)) {
    cli_abort(
      "Sample name `{reference_name}` is already in `counts`. Choose another `reference_name`."
    )
  }

  ref_vec <- reference_counts_vector(reference)
  shared <- intersect(rownames(counts), names(ref_vec))
  if (length(shared) < 2L) {
    cli_abort(
      "Need at least 2 shared genes between `counts` and the reference; found {length(shared)}."
    )
  }

  combined <- cbind(
    counts[shared, , drop = FALSE],
    matrix(ref_vec[shared], ncol = 1L, dimnames = list(shared, reference_name))
  )
  storage.mode(combined) <- "double"
  attr(combined, "shared_features") <- shared
  attr(combined, "reference_name") <- reference_name
  combined
}

#' TMM / TMMwsp effective library sizes with an HCA reference
#'
#' Runs [edgeR::calcNormFactors()] jointly on user libraries plus the HCA
#' reference (`refColumn = reference_name`). Returns effective library sizes on
#' that shared scale. The edgeR fitting offset for user samples is
#' `log(effective_size)`. HCA-scale reporting uses
#' `hca_log_effective_library_size` after model fitting; it is not folded into
#' the edgeR offset.
#'
#' @param counts Gene-by-sample count matrix that already includes the
#'   reference column.
#' @param reference_name Name of the reference column.
#' @param method Passed to [edgeR::calcNormFactors()]. Default `"TMMwsp"`.
#' @return A list with `norm_factors`, `library_size`, `effective_size`,
#'   `log_effective_library_size`, `hca_effective_library_size`,
#'   `hca_log_effective_library_size`, `multiplier`, `reference_name`, and
#'   `method`. `multiplier` is `E_H / E_j` (descriptive only).
#' @seealso [merge_with_reference_sample()], [scale_to_hca_reference()]
#' @export
#' @importFrom cli cli_abort
calculate_tmm_scaling <- function(counts, reference_name, method = "TMMwsp") {
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  reference_name <- as.character(reference_name)
  ref_col <- match(reference_name, colnames(counts))
  if (is.na(ref_col)) {
    cli_abort("Reference column `{reference_name}` was not found in `counts`.")
  }

  norm_factors <- edgeR::calcNormFactors(
    counts,
    refColumn = ref_col,
    method = method
  )
  library_size <- colSums(counts)
  effective_size <- library_size * norm_factors

  # multiplier <- as.numeric(effective_size[[reference_name]] / effective_size)
  # offset <- log(1 / multiplier)

  log_effective_library_size <- log(effective_size)
  hca_effective_library_size <- unname(as.numeric(effective_size[[reference_name]]))
  hca_log_effective_library_size <- log(hca_effective_library_size)
  multiplier <- as.numeric(hca_effective_library_size / effective_size)
  names(norm_factors) <- colnames(counts)
  names(library_size) <- colnames(counts)
  names(effective_size) <- colnames(counts)
  names(log_effective_library_size) <- colnames(counts)
  names(multiplier) <- colnames(counts)

  list(
    norm_factors = norm_factors,
    library_size = library_size,
    effective_size = effective_size,
    log_effective_library_size = log_effective_library_size,
    hca_effective_library_size = hca_effective_library_size,
    hca_log_effective_library_size = hca_log_effective_library_size,
    multiplier = multiplier,
    reference_name = reference_name,
    method = method
  )
}

#' TMM-align user counts to an atlas reference sample
#'
#' Convenience wrapper around [load_reference_sample()] (when `reference` is
#' a cell-type string), [merge_with_reference_sample()], and
#' [calculate_tmm_scaling()] for Seurat / SummarizedExperiment / matrix input.
#'
#' Returns an object of the same class with sample-level
#' `effective_library_size`, `log_effective_library_size`,
#' `hca_effective_library_size`, `hca_log_effective_library_size`,
#' `hca_multiplier`, and `sample_role` (or as matrix attributes). For a
#' matrix-only pipeline, prefer calling the core helpers directly.
#'
#' @param counts User libraries (matrix, SE/SCE, or Seurat).
#' @param reference Cell type string, Nectar download list, or one-library
#'   count object.
#' @param reference_name Column name for the bound atlas library.
#' @param method Passed to [edgeR::calcNormFactors()].
#' @param assay Assay name for SE / Seurat input.
#' @param version Nectar version pin when `reference` is a cell type.
#' @return Same class as `counts`, with the atlas library appended.
#' @seealso [merge_with_reference_sample()], [calculate_tmm_scaling()],
#'   [load_reference_sample()], [estimate_cohort_logmu()]
#' @export
#' @importFrom cli cli_abort cli_alert_info
scale_to_hca_reference <- function(
  counts,
  reference,
  reference_name = NULL,
  method = "TMMwsp",
  assay = NULL,
  version = "latest"
) {
  user <- extract_counts(counts, assay = assay)
  ref <- resolve_reference(reference, version = version)
  if (is.null(reference_name) || !nzchar(reference_name)) {
    reference_name <- ref$sample_id
  }
  if (is.null(reference_name) || !nzchar(reference_name)) {
    reference_name <- "hca_reference"
  }

  combined <- merge_with_reference_sample(
    user,
    reference = ref$counts,
    reference_name = reference_name
  )
  scaling <- calculate_tmm_scaling(
    combined,
    reference_name = reference_name,
    method = method
  )

  sample_role <- ifelse(colnames(combined) == reference_name, "reference", "user")
  names(sample_role) <- colnames(combined)
  cli_alert_info(
    "Aligned {sum(sample_role == 'user')} user sample{?s} to `{reference_name}` on {nrow(combined)} shared gene{?s}."
  )

  sample_metadata <- data.frame(
    sample_id = colnames(combined),
    sample_role = unname(sample_role),
    effective_library_size = unname(as.numeric(scaling$effective_size)),
    log_effective_library_size = unname(as.numeric(scaling$log_effective_library_size)),
    hca_effective_library_size = scaling$hca_effective_library_size,
    hca_log_effective_library_size = scaling$hca_log_effective_library_size,
    hca_multiplier = unname(scaling$multiplier),
    hca_reference_name = reference_name,
    hca_cell_type = if (!is.null(ref$cell_type)) ref$cell_type else NA_character_,
    row.names = colnames(combined),
    stringsAsFactors = FALSE
  )

  user_meta <- extract_sample_metadata(counts)
  if (nrow(user_meta) && !is.null(rownames(user_meta))) {
    extra <- setdiff(names(user_meta), names(sample_metadata))
    if (length(extra)) {
      idx <- match(sample_metadata$sample_id, rownames(user_meta))
      for (nm in extra) {
        sample_metadata[[nm]] <- user_meta[[nm]][idx]
      }
    }
  }

  if (inherits(counts, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      cli_abort("Package Seurat is required for Seurat input.")
    }
    sparse <- if (requireNamespace("Matrix", quietly = TRUE)) {
      Matrix::Matrix(combined, sparse = TRUE)
    } else {
      combined
    }
    out <- Seurat::CreateSeuratObject(counts = sparse, meta.data = sample_metadata)
    if (!identical(colnames(out), colnames(combined))) {
      out <- Seurat::RenameCells(out, new.names = colnames(combined))
    }
    return(out)
  }

  if (inherits(counts, "SummarizedExperiment")) {
    out <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = combined),
      colData = S4Vectors::DataFrame(sample_metadata)
    )
    return(out)
  }

  attr(combined, "effective_library_size") <- scaling$effective_size
  attr(combined, "log_effective_library_size") <- scaling$log_effective_library_size
  attr(combined, "hca_effective_library_size") <- scaling$hca_effective_library_size
  attr(combined, "hca_log_effective_library_size") <- scaling$hca_log_effective_library_size
  attr(combined, "hca_multiplier") <- scaling$multiplier
  attr(combined, "sample_role") <- sample_role
  attr(combined, "reference_name") <- reference_name
  attr(combined, "cell_type") <- sample_metadata$hca_cell_type[[1]]
  combined
}

#' Fit edgeR QL with prior.count = 0 and explicit log-library-size offset
#'
#' Passes `prior.count = 0` on Bioconductor release edgeR so coefficients are
#' unshrunk. Bioconductor devel `glmQLFit()` already hardcodes
#' `prior.count = 0` when calling `glmFit()` and also forwards `...`, so
#' passing it again errors with a duplicate formal argument.
#' @keywords internal
#' @noRd
fit_nb_ql <- function(counts, offset, design, robust = TRUE) {
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  design <- as.matrix(design)

  if (is.matrix(offset)) {
    if (!identical(dim(offset), dim(counts))) {
      cli::cli_abort("Matrix `offset` must have the same dimensions as `counts`.")
    }
    offset_mat <- offset
    storage.mode(offset_mat) <- "double"
  } else {
    offset <- as.numeric(offset)
    if (!is.null(names(offset)) && !is.null(colnames(counts))) {
      offset <- offset[colnames(counts)]
      if (anyNA(offset)) {
        cli::cli_abort("`offset` is missing values for one or more samples in `counts`.")
      }
    }
    if (length(offset) != ncol(counts)) {
      cli::cli_abort(
        "`offset` length must match the number of columns in `counts`."
      )
    }
    offset_mat <- matrix(
      offset,
      nrow = nrow(counts),
      ncol = ncol(counts),
      byrow = TRUE
    )
  }

  # Explicit edgeR inputs: raw counts + design + log(effective library size).
  # TMM factors are already inside `offset` from calculate_tmm_scaling();
  # do not recreate a DGEList or recompute norm factors here.
  disp <- edgeR::estimateDisp(
    y = counts,
    design = design,
    offset = offset_mat,
    robust = robust
  )
  ql_args <- list(
    y = counts,
    design = design,
    offset = offset_mat,
    robust = robust
  )
  # Release glmQLFit forwards ... to glmFit (default prior.count = 0.125).
  # Devel hardcodes prior.count = 0 and also forwards ..., so passing it
  # again fails with a duplicate formal.
  if (utils::packageVersion("edgeR") < "4.99.0") {
    ql_args$prior.count <- 0
  }
  fit <- do.call(edgeR::glmQLFit, ql_args)

  dispersion <- if (!is.null(disp$trended.dispersion)) {
    disp$trended.dispersion
  } else if (!is.null(disp$tagwise.dispersion)) {
    disp$tagwise.dispersion
  } else {
    rep_len(disp$common.dispersion, nrow(counts))
  }
  names(dispersion) <- rownames(counts)

  list(fit = fit, dispersion = dispersion, design = design, offset = offset_mat)
}

#' Check whether a design is a one-hot group-mean design
#' @keywords internal
#' @noRd
is_group_mean_design <- function(design) {
  design <- as.matrix(design)
  if (!nrow(design) || !ncol(design) || any(!is.finite(design))) {
    return(FALSE)
  }
  is_indicator <- abs(design) < .Machine$double.eps |
    abs(design - 1) < .Machine$double.eps
  all(is_indicator) &&
    all(abs(rowSums(design) - 1) < .Machine$double.eps) &&
    all(colSums(design) > .Machine$double.eps)
}

#' Validate a one-hot group-mean design for the posteriorHCA convenience workflow
#'
#' Used by [estimate_cohort_logmu()] only. Not a requirement of the core
#' modelling API [estimate_ql()], which accepts arbitrary edgeR-compatible
#' designs.
#'
#' @keywords internal
#' @noRd
assert_group_mean_design <- function(design) {
  design <- as.matrix(design)
  if (!nrow(design) || !ncol(design)) {
    cli::cli_abort("`design` must be a non-empty matrix.")
  }
  if (any(!is.finite(design))) {
    cli::cli_abort("`design` must contain only finite values.")
  }
  if (!is_group_mean_design(design)) {
    cli::cli_abort(c(
      "This convenience workflow requires a one-hot group-mean design.",
      "i" = "Use `~ 0 + group` (one mean per group) or `~ 1` (single group).",
      "i" = "Use the lower-level core API [estimate_ql()] for arbitrary edgeR-compatible designs."
    ))
  }
  invisible(design)
}

#' Basic technical checks for an edgeR design matrix
#' @keywords internal
#' @noRd
assert_design_matrix <- function(design, n_samples) {
  design <- as.matrix(design)
  if (!nrow(design) || !ncol(design)) {
    cli::cli_abort("`design` must be a non-empty matrix.")
  }
  if (any(!is.finite(design))) {
    cli::cli_abort("`design` must contain only finite values.")
  }
  if (nrow(design) != n_samples) {
    cli::cli_abort(
      "`design` rows ({nrow(design)}) must match `ncol(counts)` ({n_samples})."
    )
  }
  invisible(design)
}

#' Gene-specific QL coefficient covariance `s2.post * (X'WX)^{-1}`
#'
#' Uses the NB GLM working weights from the same `glmQLFit()` object:
#' `w = prior_weight * mu / (1 + phi * mu)`, with `phi` equal to the NB
#' dispersion used for the final fitted means (stored `fit$dispersion`
#' divided by `fit$average.ql.dispersion` when the latter is present).
#'
#' @return A `p x p` matrix, or `NULL` if the covariance is not estimable.
#' @keywords internal
#' @noRd
vcov_ql_gene <- function(fit, design, gene) {
  design <- as.matrix(design)
  mu <- as.matrix(fit$fitted.values)
  n_gene <- nrow(mu)
  g <- resolve_ql_gene_index(fit, gene, n_gene)

  if (is.null(fit$dispersion)) {
    cli::cli_abort(
      "QL fit is missing NB `dispersion` needed for coefficient covariance."
    )
  }
  phi <- ql_fit_nb_dispersion(fit, n_gene)
  if (length(phi) != n_gene) {
    cli::cli_abort(
      "NB dispersion length ({length(phi)}) does not match number of genes ({n_gene})."
    )
  }

  s2 <- as.numeric(fit$s2.post)[[g]]
  phi_g <- phi[[g]]
  if (!is.finite(s2) || s2 < 0 || !is.finite(phi_g) || phi_g < 0) {
    return(NULL)
  }

  prior_w <- fit$weights
  if (is.null(prior_w)) {
    w_g <- mu[g, ] / (1 + phi_g * mu[g, ])
  } else {
    w_g <- as.matrix(prior_w)[g, ] * mu[g, ] / (1 + phi_g * mu[g, ])
  }
  if (any(!is.finite(w_g)) || any(w_g < 0) || !any(w_g > 0)) {
    return(NULL)
  }

  xtwx <- crossprod(design, design * w_g)
  tryCatch(
    s2 * chol2inv(chol(xtwx)),
    error = function(e) NULL
  )
}

#' @keywords internal
#' @noRd
resolve_ql_gene_index <- function(fit, gene, n_gene = nrow(fit$fitted.values)) {
  if (is.numeric(gene) && length(gene) == 1L && is.finite(gene)) {
    g <- as.integer(gene)
    if (g < 1L || g > n_gene) {
      cli::cli_abort("`gene` index {g} is out of range (1..{n_gene}).")
    }
    return(g)
  }
  gene <- as.character(gene)
  if (length(gene) != 1L) {
    cli::cli_abort("`gene` must be a single id or row index.")
  }
  ids <- rownames(fit$coefficients)
  if (is.null(ids)) {
    ids <- rownames(fit$fitted.values)
  }
  g <- match(gene, ids)
  if (is.na(g)) {
    cli::cli_abort("Gene `{gene}` was not found in the QL fit.")
  }
  g
}

#' Coefficient SE from the QL fit via per-gene weighted information
#'
#' For a one-hot design this reduces to the closed-form
#' `sqrt(s2 / sum(w in group))`.
#'
#' @keywords internal
#' @noRd
se_ql_coefficients <- function(fit, design) {
  design <- as.matrix(design)
  n_gene <- nrow(fit$fitted.values)
  n_coef <- ncol(design)
  se_mat <- matrix(NA_real_, nrow = n_gene, ncol = n_coef)
  colnames(se_mat) <- colnames(design)

  for (g in seq_len(n_gene)) {
    vcov_g <- vcov_ql_gene(fit, design, g)
    if (is.null(vcov_g)) {
      next
    }
    diag_v <- diag(vcov_g)
    se_mat[g, ] <- ifelse(is.finite(diag_v) & diag_v >= 0, sqrt(diag_v), NA_real_)
  }

  se_mat
}

#' Contrast SE `sqrt(c' V c)` for one gene and one or more contrasts
#' @keywords internal
#' @noRd
se_ql_contrasts <- function(fit, design, contrast_mat) {
  design <- as.matrix(design)
  contrast_mat <- as.matrix(contrast_mat)
  n_gene <- nrow(fit$fitted.values)
  n_contrast <- ncol(contrast_mat)
  se_mat <- matrix(NA_real_, nrow = n_gene, ncol = n_contrast)
  colnames(se_mat) <- colnames(contrast_mat)

  for (g in seq_len(n_gene)) {
    vcov_g <- vcov_ql_gene(fit, design, g)
    if (is.null(vcov_g)) {
      next
    }
    for (j in seq_len(n_contrast)) {
      cvec <- contrast_mat[, j]
      var_j <- as.numeric(crossprod(cvec, vcov_g %*% cvec))
      se_mat[g, j] <- if (is.finite(var_j) && var_j >= 0) sqrt(var_j) else NA_real_
    }
  }
  se_mat
}

#' Build a limma/edgeR contrast matrix from character, numeric, or NULL
#'
#' `contrast = NULL` returns the identity (one contrast per coefficient).
#' Character contrasts use [limma::makeContrasts()]. Interaction `:` in
#' design column names is rewritten to `.` for parsing (tidybulk-style).
#'
#' @keywords internal
#' @noRd
build_contrast_matrix <- function(contrast, design) {
  design <- as.matrix(design)
  coef_names <- colnames(design)
  if (is.null(coef_names)) {
    coef_names <- paste0("X", seq_len(ncol(design)))
    colnames(design) <- coef_names
  }
  n_coef <- ncol(design)

  if (is.null(contrast)) {
    C <- diag(n_coef)
    dimnames(C) <- list(coef_names, coef_names)
    return(C)
  }

  if (is.character(contrast)) {
    if (!length(contrast) || any(!nzchar(contrast))) {
      cli::cli_abort("`contrast` character values must be non-empty.")
    }
    # limma requires syntactically valid names; `:` -> `.` (tidybulk pattern).
    design_levels <- design
    colnames(design_levels) <- gsub(":", ".", colnames(design_levels), fixed = TRUE)
    contrast_chr <- gsub(":", ".", contrast, fixed = TRUE)
    C <- tryCatch(
      limma::makeContrasts(contrasts = contrast_chr, levels = design_levels),
      error = function(e) {
        cli::cli_abort(c(
          "Could not construct contrast from `{paste(contrast, collapse = '; ')}`.",
          "i" = "Design columns: {paste(coef_names, collapse = ', ')}",
          "i" = "For interactions, use `.` in place of `:` in the contrast string.",
          "x" = conditionMessage(e)
        ))
      }
    )
    # limma may rename "(Intercept)" -> "Intercept" in rownames; keep position order.
    if (nrow(C) != n_coef) {
      cli::cli_abort("Contrast matrix rows ({nrow(C)}) do not match design columns ({n_coef}).")
    }
    rownames(C) <- coef_names
    if (is.null(colnames(C))) {
      colnames(C) <- contrast
    }
    return(C)
  }

  if (is.numeric(contrast) && is.null(dim(contrast))) {
    if (length(contrast) != n_coef) {
      cli::cli_abort(
        "Numeric `contrast` length ({length(contrast)}) must equal number of design columns ({n_coef})."
      )
    }
    C <- matrix(as.numeric(contrast), ncol = 1L, dimnames = list(coef_names, "contrast"))
    return(C)
  }

  if (is.matrix(contrast) || is.data.frame(contrast)) {
    C <- as.matrix(contrast)
    storage.mode(C) <- "double"
    if (nrow(C) != n_coef) {
      cli::cli_abort(
        "Contrast matrix rows ({nrow(C)}) must equal number of design columns ({n_coef})."
      )
    }
    if (!ncol(C)) {
      cli::cli_abort("`contrast` matrix must have at least one column.")
    }
    if (any(!is.finite(C))) {
      cli::cli_abort("`contrast` matrix must contain only finite values.")
    }
    rownames(C) <- coef_names
    if (is.null(colnames(C))) {
      colnames(C) <- paste0("contrast", seq_len(ncol(C)))
    }
    return(C)
  }

  cli::cli_abort(
    "`contrast` must be NULL, a character string/vector, a numeric vector, or a numeric matrix."
  )
}

#' Align sample metadata rows to count columns
#' @keywords internal
#' @noRd
align_sample_metadata <- function(metadata, counts) {
  metadata <- as.data.frame(metadata)
  n_samples <- ncol(counts)
  if (nrow(metadata) == n_samples) {
    if (!is.null(colnames(counts)) && !is.null(rownames(metadata))) {
      if (!identical(rownames(metadata), colnames(counts))) {
        if (all(colnames(counts) %in% rownames(metadata))) {
          metadata <- metadata[colnames(counts), , drop = FALSE]
        }
      }
    }
    return(metadata)
  }
  if (!is.null(colnames(counts)) && !is.null(rownames(metadata)) &&
      all(colnames(counts) %in% rownames(metadata))) {
    return(metadata[colnames(counts), , drop = FALSE])
  }
  cli::cli_abort(c(
    "`metadata` rows ({nrow(metadata)}) must match `ncol(counts)` ({n_samples}).",
    "i" = "Provide one metadata row per count column, ideally with matching row/column names."
  ))
}

#' One-hot closed-form SE (special case of [se_ql_coefficients()])
#'
#' `SE_j = sqrt(s2 / sum(w_i for samples with design_ij = 1))`.
#' Used in tests to confirm the general information-matrix SE reduces to
#' this formula on one-hot designs when given the same `phi` and `s2`.
#'
#' @keywords internal
#' @noRd
se_group_mean_closed_form <- function(fit, design, phi) {
  design <- as.matrix(design)
  mu_hat <- as.matrix(fit$fitted.values)
  phi <- as.numeric(phi)
  if (length(phi) == 1L) {
    phi <- rep(phi, nrow(mu_hat))
  }
  w <- mu_hat / (1 + phi * mu_hat)
  if (!is.null(fit$weights)) {
    w <- w * as.matrix(fit$weights)
  }
  s2 <- as.numeric(fit$s2.post)
  n_gene <- nrow(mu_hat)
  n_coef <- ncol(design)
  se_mat <- matrix(NA_real_, nrow = n_gene, ncol = n_coef)
  colnames(se_mat) <- colnames(design)
  for (j in seq_len(n_coef)) {
    idx <- which(abs(design[, j]) > .Machine$double.eps)
    cohort_weight <- rowSums(w[, idx, drop = FALSE])
    se_mat[, j] <- ifelse(cohort_weight > 0, sqrt(s2 / cohort_weight), Inf)
  }
  se_mat
}

#' NB dispersion used for working weights of a QL fit
#' @keywords internal
#' @noRd
ql_fit_nb_dispersion <- function(fit, n_gene) {
  phi <- as.numeric(fit$dispersion)
  if (!is.null(fit$average.ql.dispersion)) {
    ave_ql <- as.numeric(fit$average.ql.dispersion)[[1L]]
    if (is.finite(ave_ql) && ave_ql > 0) {
      phi <- phi / ave_ql
    }
  }
  if (length(phi) == 1L) {
    phi <- rep(phi, n_gene)
  }
  phi
}

#' Optional edgeR QL estimator for a user-defined linear estimand
#'
#' Fits a negative-binomial quasi-likelihood model that can account for
#' batch/covariates via `formula`, then extracts a user-defined linear
#' predictor and its SE. This is an **estimation** helper for comparison with
#' the HCA posterior, not a differential-expression testing workflow.
#'
#' Builds `design = model.matrix(formula, data = metadata)`, fits
#' [edgeR::estimateDisp()] / [edgeR::glmQLFit()] (`prior.count = 0` on
#' Bioconductor release; devel already forces this) with
#' offset `log(effective library size)`, optionally builds contrasts with
#' [limma::makeContrasts()], and returns
#'
#' ```
#' estimate = c' beta
#' SE       = sqrt(c' Var(beta) c)
#' ```
#'
#' with `Var(beta) = s2.post * (X'WX)^{-1}` from the same QL fit.
#'
#' The `formula` determines which effects enter the model (e.g. Category,
#' Experiment, batch). The `contrast` selects which fitted quantity to
#' report (an estimand selector), for example `"CategorySAVI"` or
#' `"CategorySAVI + ExperimentB"`. posteriorHCA does **not** compute
#' marginal means, equal-weight Experiment averages, or other emmeans-style
#' summaries; the user owns that definition.
#'
#' Relative contrasts such as `"CategorySAVI - CategoryControl"` remain
#' allowed and return generic `estimate`/`se` without automatic HCA
#' scaling. Absolute estimands may be shifted by the caller as
#' `estimate + log(E_HCA)`.
#'
#' The returned `estimate` is on the same **natural-log** scale as
#' [edgeR::glmQLFit()] coefficients (`log(mu) = X beta + offset`). It is
#' **not** edgeR DE-table `logFC` (log2). [edgeR::glmQLFTest()] is not used
#' in this estimation path.
#'
#' If `contrast = NULL`, every design coefficient is returned (identity
#' contrasts). There is no silent default such as `coef = 2`.
#'
#' @param counts Gene-by-sample numeric count matrix (user libraries).
#' @param offset Numeric vector (length `ncol(counts)`) or matrix of
#'   `log(effective library size)` values.
#' @param metadata Sample-level data frame used by [stats::model.matrix()];
#'   rows should match `colnames(counts)`.
#' @param formula Model formula for [stats::model.matrix()], e.g.
#'   `~ 0 + Category`, `~ 0 + Category + Experiment`, or `~ 1 + Experiment`.
#' @param contrast Character contrast(s) for [limma::makeContrasts()], a
#'   numeric contrast vector/matrix, or `NULL` to return all coefficients.
#'   Examples: `"CategorySAVI"`, `"CategorySAVI + ExperimentB"`,
#'   `"(Intercept)"`, `"CategorySAVI - CategoryControl"`.
#' @param robust Passed to edgeR.
#' @return Data frame with `gene`, `contrast`, `estimate`, `se`, `df`,
#'   `dispersion`. `estimate` is the requested linear predictor on the
#'   edgeR natural-log model scale and is not necessarily an absolute
#'   `log_mu`. Attributes `fit`, `design`, and `contrast` store the QL
#'   fit, design matrix, and contrast matrix. Capture attributes before
#'   subsetting the data frame.
#' @seealso [estimate_cohort_logmu()], [calculate_tmm_scaling()],
#'   [welch_test_means()]
#' @export
#' @importFrom cli cli_abort cli_inform
estimate_ql <- function(
  counts,
  offset,
  metadata,
  formula,
  contrast = NULL,
  robust = TRUE
) {
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  metadata <- align_sample_metadata(metadata, counts)

  if (missing(formula) || is.null(formula)) {
    cli::cli_abort("`formula` is required (e.g. `~ 0 + Category` or `~ 1 + Experiment`).")
  }

  design <- stats::model.matrix(formula, data = metadata)
  design <- assert_design_matrix(design, ncol(counts))
  cli::cli_inform("Design columns: {paste(colnames(design), collapse = ', ')}")

  contrast_mat <- build_contrast_matrix(contrast, design)

  ql <- fit_nb_ql(counts, offset = offset, design = design, robust = robust)
  fit <- ql$fit
  dispersion <- ql$dispersion

  gene_ids <- rownames(counts)
  if (is.null(gene_ids)) {
    gene_ids <- paste0("gene", seq_len(nrow(counts)))
  }

  beta <- as.matrix(fit$coefficients)
  n_gene <- nrow(beta)
  n_contrast <- ncol(contrast_mat)
  contrast_labels <- colnames(contrast_mat)
  if (is.null(contrast_labels)) {
    contrast_labels <- paste0("contrast", seq_len(n_contrast))
  }

  est_mat <- beta %*% contrast_mat
  se_mat <- se_ql_contrasts(fit, design, contrast_mat)

  out <- data.frame(
    gene = rep(gene_ids, times = n_contrast),
    contrast = rep(contrast_labels, each = n_gene),
    estimate = as.numeric(est_mat),
    se = as.numeric(se_mat),
    df = rep(as.numeric(fit$df.residual.adj), times = n_contrast),
    dispersion = rep(as.numeric(dispersion), times = n_contrast),
    stringsAsFactors = FALSE
  )
  attr(out, "fit") <- fit
  attr(out, "design") <- design
  attr(out, "contrast") <- contrast_mat
  out
}

#' Estimate HCA-scale cohort log(μ) for a one-hot group-mean design
#'
#' Opinionated convenience wrapper around [scale_to_hca_reference()] inputs
#' and [estimate_ql()]. Assumes a one-hot group-mean design such as
#' `~ 0 + Category` (or `~ 1` for a single group), treats each coefficient
#' as an absolute group log abundance, and returns HCA-scale
#' `log_mu = estimate + hca_log_effective_library_size`.
#'
#' For arbitrary formulas/contrasts, call [estimate_ql()] directly and
#' decide yourself whether the estimand should receive the HCA scale shift
#' before [welch_test_means()].
#'
#' @param data Scaled object from [scale_to_hca_reference()] (or equivalent
#'   container with counts, `sample_role`, `log_effective_library_size`, and
#'   `hca_log_effective_library_size`).
#' @param formula One-hot group-mean formula on user sample metadata
#'   (`~ 0 + Category` or `~ 1`).
#' @param gene_ensg Optional character vector of Ensembl gene ids to keep.
#' @param assay Assay name for Seurat / SummarizedExperiment input.
#' @param robust Passed to [estimate_ql()].
#' @return Data frame with `gene`, `group`, `n`, `log_mu`, `mu`, `se`, `df`,
#'   `dispersion` on the matched HCA reference effective-library-size scale.
#' @seealso [estimate_ql()], [scale_to_hca_reference()],
#'   [expression_baseline_draws()], [compare_cohort_to_hca()]
#' @export
#' @importFrom cli cli_abort
estimate_cohort_logmu <- function(
  data,
  formula,
  gene_ensg = NULL,
  assay = NULL,
  robust = TRUE
) {
  counts <- extract_counts(data, assay = assay)
  sample_metadata <- extract_sample_metadata(data)

  if (!nrow(sample_metadata) && !is.null(attr(data, "sample_role"))) {
    sample_metadata <- data.frame(
      sample_role = unname(attr(data, "sample_role")),
      log_effective_library_size = unname(
        as.numeric(attr(data, "log_effective_library_size"))
      ),
      hca_log_effective_library_size = attr(
        data,
        "hca_log_effective_library_size"
      ),
      row.names = colnames(counts),
      stringsAsFactors = FALSE
    )
  }

  if (!nrow(sample_metadata)) {
    cli_abort(c(
      "`data` has no sample metadata.",
      "i" = "Pass a Seurat / SummarizedExperiment from [scale_to_hca_reference()], or a matrix with `sample_role` and scaling attributes."
    ))
  }

  required_cols <- c("log_effective_library_size", "hca_log_effective_library_size")
  missing_cols <- setdiff(required_cols, names(sample_metadata))
  if (length(missing_cols)) {
    cli_abort(
      "`data` must contain {missing_cols} column{?s} (or matrix attribute{?s})."
    )
  }

  if ("sample_role" %in% names(sample_metadata)) {
    user_samples <- sample_metadata$sample_role == "user"
  } else {
    user_samples <- rep(TRUE, nrow(sample_metadata))
  }
  if (!any(user_samples)) {
    cli_abort("No user samples found (`sample_role == \"user\"`).")
  }

  user_counts <- counts[, user_samples, drop = FALSE]
  user_metadata <- droplevels(sample_metadata[user_samples, , drop = FALSE])
  user_offset <- setNames(
    as.numeric(user_metadata$log_effective_library_size),
    rownames(user_metadata)
  )
  hca_log_E <- unique(as.numeric(user_metadata$hca_log_effective_library_size))
  if (length(hca_log_E) != 1L || !is.finite(hca_log_E)) {
    cli_abort("`hca_log_effective_library_size` must be a single finite value.")
  }

  design_check <- stats::model.matrix(formula, data = user_metadata)
  if (nrow(design_check) != ncol(user_counts)) {
    cli_abort(
      "`formula` produced {nrow(design_check)} design row{?s} for {ncol(user_counts)} user sample{?s}."
    )
  }
  design_check <- assert_group_mean_design(design_check)

  coef_estimates <- estimate_ql(
    counts = user_counts,
    offset = user_offset,
    metadata = user_metadata,
    formula = formula,
    contrast = NULL,
    robust = robust
  )

  group_labels <- coef_estimates$contrast
  group_labels[group_labels == "(Intercept)"] <- "all"
  group_n <- vapply(seq_len(ncol(design_check)), function(j) {
    sum(abs(design_check[, j]) > .Machine$double.eps)
  }, integer(1))
  n_gene <- nrow(user_counts)
  log_mu <- coef_estimates$estimate + hca_log_E

  cohort_estimates <- data.frame(
    gene = coef_estimates$gene,
    group = group_labels,
    n = rep(as.integer(group_n), each = n_gene),
    log_mu = log_mu,
    mu = exp(log_mu),
    se = coef_estimates$se,
    df = coef_estimates$df,
    dispersion = coef_estimates$dispersion,
    stringsAsFactors = FALSE
  )

  if (!is.null(gene_ensg)) {
    gene_ensg <- as.character(gene_ensg)
    cohort_estimates <- cohort_estimates[
      cohort_estimates$gene %in% gene_ensg,
      ,
      drop = FALSE
    ]
  }

  cohort_estimates
}
