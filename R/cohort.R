# User counts on the atlas scale
#
# 1. Read a count matrix from a matrix, SummarizedExperiment/SCE, or Seurat object.
# 2. Resolve the atlas reference from Nectar (get_reference_sample_ready) or
#    from a count vector / one-column object.
# 3. TMMwsp with that library as refColumn, offset = log(1 / multiplier).
# 4. Estimate cohort log(mu) with mglmOneGroup (not glmQLFit coefficients).

# ---------------------------------------------------------------------------
# Input coercion
# ---------------------------------------------------------------------------

is_aligned_result <- function(x) {
  is.list(x) &&
    !is.data.frame(x) &&
    all(c("counts", "offset") %in% names(x))
}

is_nectar_download <- function(x) {
  is.list(x) &&
    !is.data.frame(x) &&
    all(c("status", "path") %in% names(x))
}

#' Count matrix from common single-cell / bulk containers
#' @keywords internal
#' @noRd
extract_count_matrix <- function(x, assay = NULL, arg_name = "counts") {
  if (is.null(x)) {
    cli::cli_abort("`{arg_name}` is missing.")
  }
  if (is_aligned_result(x)) {
    return(as.matrix(x$counts))
  }

  if (inherits(x, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      cli::cli_abort("Install SummarizedExperiment to use that input class.")
    }
    assay_names <- SummarizedExperiment::assayNames(x)
    if (!length(assay_names)) {
      cli::cli_abort("`{arg_name}` has no assays.")
    }
    if (is.null(assay)) {
      assay <- if ("counts" %in% assay_names) "counts" else assay_names[[1]]
    }
    if (!assay %in% assay_names) {
      cli::cli_abort("Assay `{assay}` not found. Available: {assay_names}.")
    }
    return(as.matrix(SummarizedExperiment::assay(x, assay)))
  }

  if (inherits(x, "Seurat")) {
    return(extract_seurat_counts(x, assay = assay, arg_name = arg_name))
  }

  if (inherits(x, "Matrix") || inherits(x, "delayedMatrix")) {
    return(as.matrix(x))
  }

  if (is.matrix(x) || is.data.frame(x)) {
    mat <- as.matrix(x)
    storage.mode(mat) <- "double"
    return(mat)
  }

  cli::cli_abort(
    c(
      "Don't know how to get counts from class `{paste(class(x), collapse = ', ')}`.",
      "i" = "Use a matrix, SummarizedExperiment, SingleCellExperiment, or Seurat object."
    )
  )
}

#' @keywords internal
#' @noRd
extract_seurat_counts <- function(x, assay = NULL, arg_name = "counts") {
  seurat_ns <- if (requireNamespace("SeuratObject", quietly = TRUE)) {
    "SeuratObject"
  } else if (requireNamespace("Seurat", quietly = TRUE)) {
    "Seurat"
  } else {
    cli::cli_abort("Install Seurat or SeuratObject to use a Seurat `{arg_name}`.")
  }

  if (is.null(assay)) {
    assay <- tryCatch(
      getExportedValue(seurat_ns, "DefaultAssay")(x),
      error = function(e) x@active.assay
    )
  }

  get_data <- getExportedValue(seurat_ns, "GetAssayData")
  mat <- tryCatch(
    get_data(x, assay = assay, layer = "counts"),
    error = function(e) get_data(x, assay = assay, slot = "counts")
  )
  as.matrix(mat)
}

#' Sample metadata from SE colData or Seurat meta.data
#' @keywords internal
#' @noRd
extract_sample_metadata <- function(x) {
  if (is_aligned_result(x)) {
    if (!is.null(x$sample_metadata)) {
      return(as.data.frame(x$sample_metadata))
    }
    x <- x$counts
  }
  if (inherits(x, "SummarizedExperiment")) {
    return(as.data.frame(SummarizedExperiment::colData(x)))
  }
  if (inherits(x, "Seurat")) {
    return(as.data.frame(x[[]]))
  }
  if (is.matrix(x) || is.data.frame(x) || inherits(x, "Matrix")) {
    nms <- colnames(x)
    if (is.null(nms)) {
      return(data.frame())
    }
    return(data.frame(row.names = nms, stringsAsFactors = FALSE))
  }
  data.frame()
}

#' Named gene-count vector (one library)
#' @keywords internal
#' @noRd
as_named_gene_counts <- function(x, arg_name = "reference") {
  if (is.null(x)) {
    cli::cli_abort("`{arg_name}` is missing.")
  }
  if (is.numeric(x) && is.null(dim(x))) {
    nms <- names(x)
    if (is.null(nms) || any(!nzchar(nms))) {
      cli::cli_abort("`{arg_name}` must be a named numeric vector (names = genes).")
    }
    out <- as.numeric(x)
    names(out) <- nms
    return(out)
  }

  mat <- extract_count_matrix(x, arg_name = arg_name)
  if (ncol(mat) != 1L) {
    cli::cli_abort(
      "`{arg_name}` must be one library; found {ncol(mat)} columns."
    )
  }
  nms <- rownames(mat)
  if (is.null(nms) || any(!nzchar(nms))) {
    cli::cli_abort("`{arg_name}` must have gene names as rownames.")
  }
  out <- as.numeric(mat[, 1])
  names(out) <- nms
  attr(out, "sample_id") <- colnames(mat)[[1]]
  out
}

# ---------------------------------------------------------------------------
# Nectar reference sample
# ---------------------------------------------------------------------------

#' Interpret a downloaded atlas reference RDS
#'
#' Published Nectar objects are one library per cell type: a one-column
#' `SummarizedExperiment` with a `counts` assay (and usually
#' `counts_scaled`). Scaling always uses `counts`. A length-1 character
#' (legacy sample id) cannot be merged with user counts.
#' @keywords internal
#' @noRd
parse_reference_object <- function(obj, cell_type = NA_character_, path = NULL) {
  path_msg <- if (!is.null(path)) paste0(" (", path, ")") else ""

  if (is.character(obj) && length(obj) == 1L) {
    cli::cli_abort(c(
      "Reference file{path_msg} is a sample id (`{obj}`), not counts.",
      "i" = "Nectar should store one gene-by-1 SummarizedExperiment per cell type."
    ))
  }

  if (is.numeric(obj) && is.null(dim(obj)) && !is.null(names(obj))) {
    out <- as.numeric(obj)
    names(out) <- names(obj)
    return(list(
      cell_type = cell_type,
      sample_id = "hca_reference",
      counts = out,
      path = path
    ))
  }

  mat <- extract_count_matrix(obj, assay = "counts", arg_name = "reference")
  if (ncol(mat) != 1L) {
    cli::cli_abort(
      "Atlas reference{path_msg} must be one sample; found {ncol(mat)} columns."
    )
  }
  sid <- colnames(mat)[[1]]
  if (is.null(sid) || !nzchar(sid)) {
    sid <- "hca_reference"
  }
  counts <- as.numeric(mat[, 1])
  names(counts) <- rownames(mat)
  list(
    cell_type = cell_type,
    sample_id = sid,
    counts = counts,
    path = path
  )
}

#' Load the atlas reference library for a cell type
#'
#' Calls [get_reference_sample_ready()] (or uses its return value), reads
#' the cached one-column `SummarizedExperiment`, and returns the `counts`
#' assay to merge with user data.
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
  if (is_nectar_download(cell_type)) {
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

  obj <- readRDS(res$path)
  parsed <- parse_reference_object(
    obj,
    cell_type = if (!is.null(res$cell_type)) res$cell_type else NA_character_,
    path = res$path
  )
  cli_alert_info(
    "Loaded atlas reference `{parsed$sample_id}` ({length(parsed$counts)} gene{?s})."
  )
  parsed
}

#' Load atlas reference counts from a cell type, Nectar download, or count object
#' @keywords internal
#' @noRd
as_atlas_reference <- function(reference, version = "latest") {
  if (is.null(reference)) {
    cli::cli_abort("`reference` is missing.")
  }
  if (is.character(reference) && length(reference) == 1L) {
    return(load_reference_sample(reference, version = version))
  }
  if (is_nectar_download(reference)) {
    return(load_reference_sample(reference, version = version))
  }
  if (is.list(reference) && !is.data.frame(reference) &&
      "counts" %in% names(reference) &&
      (is.null(reference$counts) ||
       (is.numeric(reference$counts) && is.null(dim(reference$counts))))) {
    if (!is.null(reference$counts) &&
        (is.null(reference$sample_id) || is.na(reference$sample_id) ||
         !nzchar(reference$sample_id))) {
      reference$sample_id <- "hca_reference"
    }
    return(reference)
  }

  if (inherits(reference, "SummarizedExperiment")) {
    return(parse_reference_object(reference))
  }

  vec <- as_named_gene_counts(reference, "reference")
  sid <- attr(vec, "sample_id")
  list(
    cell_type = NA_character_,
    sample_id = if (!is.null(sid) && !is.na(sid) && nzchar(sid)) sid else "hca_reference",
    counts = vec,
    path = NULL
  )
}

# ---------------------------------------------------------------------------
# Public methods
# ---------------------------------------------------------------------------

#' TMM-align user counts to an atlas reference sample
#'
#' The atlas reference is provided by this package: one sample per cell
#' type, downloaded from Nectar by [get_reference_sample_ready()] and
#' cached locally. That helper's return value records the file `path`.
#' User libraries are then bound to that sample and TMMwsp-scaled
#' together, with the atlas column as `refColumn`. The returned offset is
#' `log(1 / multiplier)` so the reference sits at offset 0.
#'
#' `counts` may be a matrix, `SummarizedExperiment` /
#' `SingleCellExperiment`, or `Seurat` object.
#'
#' `reference` is usually a cell type (`"cd8 naive"`) or the list returned
#' by [get_reference_sample_ready()]. A one-column `SummarizedExperiment`,
#' named count vector, or matrix is also accepted (tests / local files).
#'
#' @param counts User libraries (matrix, SE/SCE, or Seurat).
#' @param reference Cell type, Nectar download list from
#'   [get_reference_sample_ready()], or a one-library count object.
#' @param reference_name Column name for the bound atlas library.
#'   Default is the atlas sample id, or `"hca_reference"`.
#' @param method Passed to [edgeR::calcNormFactors()]. Default `"TMMwsp"`.
#' @param assay Assay name for SE / Seurat input. Default `"counts"` when present.
#' @param version Nectar version pin, used when `reference` is a cell type.
#' @return A list with `counts` (shared genes, user columns plus reference),
#'   `offset`, `multiplier`, `shared_features`, `reference_name`, and
#'   `sample_role`.
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
  user <- extract_count_matrix(counts, assay = assay, arg_name = "counts")
  if (is.null(rownames(user)) || is.null(colnames(user))) {
    cli_abort("`counts` must have gene rownames and sample colnames.")
  }

  ref <- as_atlas_reference(reference, version = version)
  if (is.null(ref$counts) || !length(ref$counts)) {
    cli_abort("Atlas reference has no count vector to merge.")
  }

  if (is.null(reference_name) || !nzchar(reference_name)) {
    reference_name <- ref$sample_id
  }
  if (is.null(reference_name) || is.na(reference_name) || !nzchar(reference_name)) {
    reference_name <- "hca_reference"
  }
  if (reference_name %in% colnames(user)) {
    cli_abort(
      "Sample name `{reference_name}` is already in `counts`. Choose another `reference_name`."
    )
  }

  shared <- intersect(rownames(user), names(ref$counts))
  if (length(shared) < 2L) {
    cli_abort(
      "Need at least 2 shared genes between `counts` and the reference; found {length(shared)}."
    )
  }

  combined <- cbind(
    user[shared, , drop = FALSE],
    matrix(
      ref$counts[shared],
      ncol = 1L,
      dimnames = list(shared, reference_name)
    )
  )
  storage.mode(combined) <- "double"

  ref_col <- match(reference_name, colnames(combined))
  if (is.na(ref_col)) {
    cli_abort("Reference column `{reference_name}` was not found after alignment.")
  }

  norm_factors <- edgeR::calcNormFactors(
    combined,
    refColumn = ref_col,
    method = method
  )
  library_size <- colSums(combined)
  effective_size <- library_size * norm_factors
  multiplier <- as.numeric(effective_size[[reference_name]] / effective_size)
  offset <- log(1 / multiplier)
  names(multiplier) <- colnames(combined)
  names(offset) <- colnames(combined)

  sample_role <- ifelse(
    colnames(combined) == reference_name,
    "reference",
    "user"
  )
  names(sample_role) <- colnames(combined)

  n_user <- sum(sample_role == "user")
  cli_alert_info(
    "Aligned {n_user} user sample{?s} to `{reference_name}` on {length(shared)} shared gene{?s}."
  )

  user_meta <- extract_sample_metadata(counts)
  sample_metadata <- data.frame(
    sample_id = colnames(combined),
    sample_role = unname(sample_role),
    row.names = colnames(combined),
    stringsAsFactors = FALSE
  )
  extra <- setdiff(names(user_meta), names(sample_metadata))
  if (length(extra) && NROW(user_meta)) {
    for (nm in extra) {
      sample_metadata[[nm]] <- NA
    }
    hit <- intersect(rownames(user_meta), rownames(sample_metadata))
    if (length(hit)) {
      sample_metadata[hit, extra] <- user_meta[hit, extra, drop = FALSE]
    }
  }

  list(
    counts = combined,
    offset = offset,
    multiplier = multiplier,
    shared_features = shared,
    reference_name = reference_name,
    sample_role = sample_role,
    sample_metadata = sample_metadata
  )
}

#' Negative-binomial dispersion vector from an edgeR DGEList
#' @keywords internal
#' @noRd
dispersion_from_dge <- function(dge) {
  disp <- if (!is.null(dge$trended.dispersion)) {
    dge$trended.dispersion
  } else if (!is.null(dge$tagwise.dispersion)) {
    dge$tagwise.dispersion
  } else {
    rep_len(dge$common.dispersion, nrow(dge))
  }
  names(disp) <- rownames(dge)
  disp
}

#' Resolve a group vector from a column name or a vector
#' @keywords internal
#' @noRd
resolve_group <- function(group, counts, n_lib, sample_role = NULL) {
  if (length(group) == 1L && is.character(group)) {
    meta <- extract_sample_metadata(counts)
    if (NROW(meta) && group %in% names(meta)) {
      group <- as.character(meta[[group]])
    }
  }

  if (!is.null(sample_role) && length(group) == sum(sample_role == "user") &&
      length(group) + 1L == n_lib) {
    full <- character(n_lib)
    full[sample_role == "user"] <- as.character(group)
    full[sample_role == "reference"] <- "reference"
    group <- full
  }

  if (length(group) != n_lib) {
    cli::cli_abort(
      "`group` length ({length(group)}) must match the number of libraries ({n_lib}), or be a metadata column name."
    )
  }
  group
}

#' Cohort log(mu) at the atlas offset-zero scale
#'
#' Point estimates come from [edgeR::mglmOneGroup()] with the TMM offset
#' applied directly. [edgeR::glmQLFit()] is used only for the QL Wald SE.
#'
#' `counts` may be a matrix, SE/SCE, Seurat object, or the list returned
#' by [scale_to_hca_reference()]. `group` may be a vector or the name of
#' a metadata column. If `counts` is an aligned list and `group` covers
#' only the user samples, `"reference"` is appended for the atlas library.
#'
#' @param counts Gene-by-sample counts, or an aligned list from
#'   [scale_to_hca_reference()].
#' @param offset Named or unnamed numeric vector, one value per sample.
#'   Taken from `counts$offset` when `counts` is an aligned list.
#' @param group Group label per sample, or a column name in SE / Seurat
#'   metadata.
#' @param genes Optional gene ids to report. Default is all rows. Dispersion
#'   is still estimated from the full `counts` matrix.
#' @param robust Passed to [edgeR::estimateDisp()] and [edgeR::glmQLFit()].
#' @param assay Assay name for SE / Seurat input.
#' @return A data frame with one row per gene x group: `gene`, `group`,
#'   `n`, `log_mu`, `mu`, `se`, `df`, `dispersion`.
#' @export
#' @importFrom cli cli_abort
estimate_cohort_logmu <- function(
  counts,
  offset = NULL,
  group,
  genes = NULL,
  robust = TRUE,
  assay = NULL
) {
  sample_role <- NULL
  source_obj <- counts
  if (is_aligned_result(counts)) {
    if (is.null(offset)) {
      offset <- counts$offset
    }
    sample_role <- counts$sample_role
    counts <- counts$counts
  } else {
    counts <- extract_count_matrix(source_obj, assay = assay, arg_name = "counts")
  }

  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    cli_abort("`counts` must have gene rownames and sample colnames.")
  }
  if (is.null(offset)) {
    cli_abort("`offset` is missing.")
  }
  if (length(offset) != ncol(counts)) {
    cli_abort("`offset` length ({length(offset)}) must match `ncol(counts)` ({ncol(counts)}).")
  }
  if (!is.null(names(offset)) && all(names(offset) %in% colnames(counts))) {
    offset <- offset[colnames(counts)]
  }

  group <- resolve_group(group, source_obj, ncol(counts), sample_role)
  group <- as.character(group)
  group[is.na(group) | !nzchar(group)] <- "reference"
  samples <- data.frame(
    sample_id = colnames(counts),
    group = factor(group),
    stringsAsFactors = FALSE
  )

  if (is.null(genes)) {
    genes <- rownames(counts)
  } else {
    missing <- setdiff(genes, rownames(counts))
    if (length(missing) > 0L) {
      cli_abort("Gene{?s} not in `counts`: {missing}.")
    }
  }

  dge <- edgeR::DGEList(counts = counts, samples = samples)
  dge$offset <- matrix(
    as.numeric(offset),
    nrow = nrow(dge),
    ncol = ncol(dge),
    byrow = TRUE
  )

  design <- stats::model.matrix(~ 0 + group, data = samples)
  colnames(design) <- sub("^group", "", colnames(design))

  dge <- edgeR::estimateDisp(dge, design, robust = robust)
  fit <- edgeR::glmQLFit(dge, design, robust = robust)
  dispersion <- dispersion_from_dge(dge)

  out <- lapply(genes, function(gene) {
    gi <- match(gene, rownames(fit))
    disp <- dispersion[[gene]]
    mu_all <- fit$fitted.values[gi, ]
    w <- mu_all / (1 + disp * mu_all)

    lapply(levels(samples$group), function(g) {
      j <- which(samples$group == g)
      log_mu <- edgeR::mglmOneGroup(
        matrix(counts[gi, j], nrow = 1L),
        dispersion = disp,
        offset = as.numeric(offset)[j]
      )
      cohort_weight <- sum(w[j])
      se <- if (is.finite(cohort_weight) && cohort_weight > 0) {
        sqrt(fit$s2.post[gi] / cohort_weight)
      } else {
        Inf
      }
      data.frame(
        gene = gene,
        group = g,
        n = length(j),
        log_mu = unname(as.numeric(log_mu)),
        mu = unname(exp(as.numeric(log_mu))),
        se = unname(as.numeric(se)),
        df = unname(fit$df.residual.adj[gi]),
        dispersion = unname(disp),
        stringsAsFactors = FALSE
      )
    })
  })

  do.call(rbind, unlist(out, recursive = FALSE))
}
