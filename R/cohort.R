# User counts on the atlas scale
#
#   Core:
#     merge_with_reference_sample() / calculate_tmm_offset()
#     estimate_logmu_ql()
#   Wrappers (compose cores only):
#     scale_to_hca_reference()   — load + merge + TMM
#     estimate_cohort_logmu()    — extract + model.matrix + estimate_logmu_ql

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
#' @seealso [calculate_tmm_offset()], [scale_to_hca_reference()]
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

#' TMM normalisation factors, multipliers, and log offsets
#'
#' Runs [edgeR::calcNormFactors()] with `refColumn = reference_name`. Returns
#' `multiplier_j = effectiveSize_ref / effectiveSize_j` and
#' `offset_j = log(1 / multiplier_j)` so the reference has multiplier 1 and
#' offset 0.
#'
#' @param counts Gene-by-sample count matrix that already includes the
#'   reference column.
#' @param reference_name Name of the reference column.
#' @param method Passed to [edgeR::calcNormFactors()]. Default `"TMMwsp"`.
#' @return A list with `norm_factors`, `library_size`, `effective_size`,
#'   `multiplier`, `offset`, `reference_name`, and `method`.
#' @seealso [merge_with_reference_sample()], [scale_to_hca_reference()]
#' @export
#' @importFrom cli cli_abort
calculate_tmm_offset <- function(counts, reference_name, method = "TMMwsp") {
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
  multiplier <- as.numeric(effective_size[[reference_name]] / effective_size)
  offset <- log(1 / multiplier)
  names(norm_factors) <- colnames(counts)
  names(library_size) <- colnames(counts)
  names(effective_size) <- colnames(counts)
  names(multiplier) <- colnames(counts)
  names(offset) <- colnames(counts)

  list(
    norm_factors = norm_factors,
    library_size = library_size,
    effective_size = effective_size,
    multiplier = multiplier,
    offset = offset,
    reference_name = reference_name,
    method = method
  )
}

#' TMM-align user counts to an atlas reference sample
#'
#' Convenience wrapper around [load_reference_sample()] (when `reference` is
#' a cell-type string), [merge_with_reference_sample()], and
#' [calculate_tmm_offset()] for Seurat / SummarizedExperiment / matrix input.
#'
#' Returns an object of the same class with `hca_offset`, `hca_multiplier`,
#' and `sample_role` on sample metadata (or as matrix attributes). For a
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
#' @seealso [merge_with_reference_sample()], [calculate_tmm_offset()],
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
  scaling <- calculate_tmm_offset(
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
    hca_offset = unname(scaling$offset),
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

  attr(combined, "hca_offset") <- scaling$offset
  attr(combined, "hca_multiplier") <- scaling$multiplier
  attr(combined, "sample_role") <- sample_role
  attr(combined, "reference_name") <- reference_name
  attr(combined, "cell_type") <- sample_metadata$hca_cell_type[[1]]
  combined
}

#' Fit edgeR QL with prior.count = 0
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
    offset_mat <- matrix(offset, nrow = nrow(counts), ncol = ncol(counts), byrow = TRUE)
  }

  dge <- edgeR::DGEList(counts = counts)
  dge$offset <- offset_mat
  
  # stopifnot(
  #   isTRUE(all.equal(
  #     as.matrix(edgeR::getOffset(dge)),
  #     as.matrix(offset_mat)
  #   ))
  # )
  
  dge <- edgeR::estimateDisp(dge, design, robust = robust)
  fit <- edgeR::glmQLFit(dge, design, robust = robust, prior.count = 0)

  disp <- if (!is.null(dge$trended.dispersion)) {
    dge$trended.dispersion
  } else if (!is.null(dge$tagwise.dispersion)) {
    dge$tagwise.dispersion
  } else {
    rep_len(dge$common.dispersion, nrow(dge))
  }
  names(disp) <- rownames(dge)

  list(fit = fit, dispersion = disp, design = design)
}

#' Cohort log(μ) and SE for all genes from one edgeR QL fit
#'
#' Fits [edgeR::estimateDisp()] and [edgeR::glmQLFit()] once with
#' `prior.count = 0`. Pass the **user** count matrix and corresponding
#' offsets; do not include the HCA reference as a design group.
#'
#' @param counts Gene-by-sample numeric count matrix.
#' @param offset Numeric vector (length `ncol(counts)`) or matrix.
#' @param design Numeric design matrix from [stats::model.matrix()].
#' @param robust Passed to edgeR.
#' @return Data frame with `gene`, `group`, `n`, `log_mu`, `mu`, `se`,
#'   `df`, `dispersion`.
#' @seealso [estimate_cohort_logmu()]
#' @export
#' @importFrom cli cli_abort
estimate_logmu_ql <- function(counts, offset, design, robust = TRUE) {
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  design <- as.matrix(design)
  if (nrow(design) != ncol(counts)) {
    cli_abort(
      "`design` rows ({nrow(design)}) must match `ncol(counts)` ({ncol(counts)})."
    )
  }

  ql <- fit_nb_ql(counts, offset = offset, design = design, robust = robust)
  fit <- ql$fit
  dispersion <- ql$dispersion

  gene_ids <- rownames(counts)
  group_labels <- colnames(design)
  if (is.null(group_labels)) {
    group_labels <- paste0("group", seq_len(ncol(design)))
  }
  group_labels[group_labels == "(Intercept)"] <- "all"

  group_n <- vapply(seq_len(ncol(design)), function(j) {
    sum(abs(design[, j]) > .Machine$double.eps)
  }, integer(1))

  mu_hat <- fit$fitted.values
  w <- mu_hat / (1 + as.numeric(dispersion) * mu_hat)
  s2 <- fit$s2.post
  coef_mat <- fit$coefficients
  n_gene <- nrow(coef_mat)
  n_group <- ncol(coef_mat)

  se_mat <- matrix(NA_real_, nrow = n_gene, ncol = n_group)
  for (j in seq_len(n_group)) {
    idx <- which(abs(design[, j]) > .Machine$double.eps)
    cohort_weight <- rowSums(w[, idx, drop = FALSE])
    se_mat[, j] <- ifelse(cohort_weight > 0, sqrt(s2 / cohort_weight), Inf)
  }

  data.frame(
    gene = rep(gene_ids, times = n_group),
    group = rep(group_labels, each = n_gene),
    n = rep(as.integer(group_n), each = n_gene),
    log_mu = as.numeric(coef_mat),
    mu = exp(as.numeric(coef_mat)),
    se = as.numeric(se_mat),
    df = rep(as.numeric(fit$df.residual.adj), times = n_group),
    dispersion = rep(as.numeric(dispersion), times = n_group),
    stringsAsFactors = FALSE
  )
}

#' Estimate cohort log(μ) from a scaled count container
#'
#' Convenience wrapper around [stats::model.matrix()] and
#' [estimate_logmu_ql()] for scaled Seurat / SummarizedExperiment / matrix
#' inputs from [scale_to_hca_reference()].
#'
#' Extracts user libraries and their `hca_offset` values, builds the design
#' matrix from `formula`, and delegates all edgeR QL estimation to
#' [estimate_logmu_ql()]. Does not call edgeR itself.
#'
#' @param data Scaled object from [scale_to_hca_reference()] (or equivalent
#'   container with counts, `sample_role`, and `hca_offset`).
#' @param formula Model formula evaluated on user sample metadata
#'   (for example `~ 0 + Category`).
#' @param gene_ensg Optional character vector of Ensembl gene ids to keep.
#' @param assay Assay name for Seurat / SummarizedExperiment input.
#' @param robust Passed to [estimate_logmu_ql()].
#' @return Data frame from [estimate_logmu_ql()], optionally filtered by
#'   `gene_ensg`.
#' @seealso [estimate_logmu_ql()], [scale_to_hca_reference()],
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
      hca_offset = unname(as.numeric(attr(data, "hca_offset"))),
      row.names = colnames(counts),
      stringsAsFactors = FALSE
    )
  }

  if (!nrow(sample_metadata)) {
    cli_abort(c(
      "`data` has no sample metadata.",
      "i" = "Pass a Seurat / SummarizedExperiment from [scale_to_hca_reference()], or a matrix with `sample_role` and `hca_offset` attributes."
    ))
  }

  if (!"hca_offset" %in% names(sample_metadata)) {
    cli_abort("`data` must contain an `hca_offset` column (or matrix attribute).")
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
    as.numeric(user_metadata$hca_offset),
    rownames(user_metadata)
  )

  design_matrix <- stats::model.matrix(formula, data = user_metadata)
  if (nrow(design_matrix) != ncol(user_counts)) {
    cli_abort(
      "`formula` produced {nrow(design_matrix)} design row{?s} for {ncol(user_counts)} user sample{?s}."
    )
  }

  cohort_estimates <- estimate_logmu_ql(
    counts = user_counts,
    offset = user_offset,
    design = design_matrix,
    robust = robust
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
