# User counts on the atlas scale
#
#   merge_with_reference_sample() — bind user matrix + one reference library
#   calculate_tmm_offset() — TMMwsp; offset = log(1 / multiplier), ref at 0
#   estimate_logmu_ql() — estimateDisp + glmQLFit(prior.count = 0)
#   scale_to_hca_reference() — Seurat / SummarizedExperiment convenience
#   bootstrap_logmu_mglm() — weighted mglmOneGroup (inference-cohort.R)


# ---------------------------------------------------------------------------
# Input coercion
# ---------------------------------------------------------------------------

.posteriorhca_aligned_attr <- "posteriorHCA_aligned"

is_aligned_list <- function(x) {
  is.list(x) &&
    !is.data.frame(x) &&
    !inherits(x, "Seurat") &&
    !inherits(x, "SummarizedExperiment") &&
    all(c("counts", "offset") %in% names(x))
}

is_aligned_result <- function(x) {
  is_aligned_list(x) || !is.null(attr(x, .posteriorhca_aligned_attr))
}

#' @keywords internal
#' @noRd
aligned_sample_metadata <- function(x) {
  if (inherits(x, "SummarizedExperiment")) {
    return(as.data.frame(SummarizedExperiment::colData(x)))
  }
  if (inherits(x, "Seurat")) {
    return(as.data.frame(x[[]]))
  }
  attr_meta <- attr(x, .posteriorhca_aligned_attr)
  if (!is.null(attr_meta$sample_metadata)) {
    return(as.data.frame(attr_meta$sample_metadata))
  }
  data.frame()
}

#' @keywords internal
#' @noRd
aligned_fields <- function(x) {
  if (is_aligned_list(x)) {
    return(x)
  }

  attr_meta <- attr(x, .posteriorhca_aligned_attr)
  sample_df <- aligned_sample_metadata(x)
  if (nrow(sample_df) && "hca_offset" %in% names(sample_df)) {
    sid <- rownames(sample_df)
    if (is.null(sid) || !length(sid)) {
      sid <- sample_df$sample_id
    }
    return(list(
      counts = extract_count_matrix(x, arg_name = "counts"),
      offset = setNames(sample_df$hca_offset, sid),
      multiplier = setNames(sample_df$hca_multiplier, sid),
      sample_role = setNames(sample_df$sample_role, sid),
      reference_name = sample_df$hca_reference_name[[1]],
      cell_type = sample_df$hca_cell_type[[1]],
      shared_features = if (!is.null(attr_meta)) attr_meta$shared_features else NULL,
      sample_metadata = sample_df
    ))
  }

  if (!is.null(attr_meta)) {
    c(list(counts = extract_count_matrix(x, arg_name = "counts")), attr_meta)
  } else {
    NULL
  }
}

#' @keywords internal
#' @noRd
attach_aligned_metadata <- function(x, fields) {
  meta <- fields
  meta$counts <- NULL
  attr(x, .posteriorhca_aligned_attr) <- meta
  x
}

is_nectar_download <- function(x) {
  is.list(x) &&
    !is.data.frame(x) &&
    all(c("status", "path") %in% names(x))
}

#' Detect the container class of a count input
#' @keywords internal
#' @noRd
count_input_class <- function(x) {
  if (inherits(x, "Seurat")) {
    return("Seurat")
  }
  if (inherits(x, "SummarizedExperiment")) {
    return("SummarizedExperiment")
  }
  if (is_aligned_list(x)) {
    return("aligned")
  }
  if (inherits(x, "Matrix") || is.matrix(x) || is.data.frame(x)) {
    return("matrix")
  }
  "unknown"
}

#' Sample identifiers for a count container
#' @keywords internal
#' @noRd
sample_ids_from_counts <- function(x, assay = NULL) {
  colnames(extract_count_matrix(x, assay = assay, arg_name = "counts"))
}

#' Rebuild a count container with a harmonised matrix
#' @keywords internal
#' @noRd
rebuild_count_container <- function(x, mat, assay = NULL) {
  cls <- count_input_class(x)
  switch(
    cls,
    Seurat = rebuild_seurat_counts(x, mat, assay = assay),
    SummarizedExperiment = rebuild_se_counts(x, mat, assay = assay),
    aligned = rebuild_aligned_counts(x, mat),
    matrix = {
      storage.mode(mat) <- "double"
      mat
    },
    cli::cli_abort(
      "Don't know how to rebuild class `{paste(class(x), collapse = ', ')}`."
    )
  )
}

#' Coerce count matrix to sparse dgCMatrix when Matrix is available
#' @keywords internal
#' @noRd
as_sparse_counts <- function(mat) {
  if (inherits(mat, "Matrix") || inherits(mat, "dgCMatrix")) {
    return(mat)
  }
  if (requireNamespace("Matrix", quietly = TRUE)) {
    return(Matrix::Matrix(mat, sparse = TRUE))
  }
  mat
}

#' Restore original cell/sample names on a rebuilt Seurat object
#' @keywords internal
#' @noRd
restore_seurat_cell_names <- function(obj, cell_names, seurat_ns) {
  desired <- as.character(cell_names)
  current <- colnames(obj)
  if (length(current) != length(desired)) {
    return(obj)
  }
  if (!identical(current, desired)) {
    rename_fn <- getExportedValue(seurat_ns, "RenameCells")
    obj <- rename_fn(obj, new.names = desired)
  }
  obj
}

#' Merge user metadata into alignment sample metadata without coercing types
#' @keywords internal
#' @noRd
merge_user_sample_metadata <- function(sample_metadata, user_meta) {
  if (!is.data.frame(user_meta) || nrow(user_meta) == 0L) {
    return(sample_metadata)
  }
  if (is.null(rownames(user_meta))) {
    return(sample_metadata)
  }

  extra <- setdiff(names(user_meta), names(sample_metadata))
  if (!length(extra)) {
    return(sample_metadata)
  }

  idx <- match(sample_metadata$sample_id, rownames(user_meta))
  for (nm in extra) {
    col <- user_meta[[nm]]
    merged <- col[idx]
    if (is.factor(col)) {
      merged <- factor(merged, levels = levels(col))
    }
    sample_metadata[[nm]] <- merged
  }

  sample_metadata
}

#' @keywords internal
#' @noRd
rebuild_seurat_counts <- function(obj, mat, assay = NULL) {
  seurat_ns <- if (requireNamespace("SeuratObject", quietly = TRUE)) {
    "SeuratObject"
  } else if (requireNamespace("Seurat", quietly = TRUE)) {
    "Seurat"
  } else {
    cli::cli_abort("Install Seurat or SeuratObject to rebuild a Seurat object.")
  }

  if (is.null(assay)) {
    assay <- tryCatch(
      getExportedValue(seurat_ns, "DefaultAssay")(obj),
      error = function(e) obj@active.assay
    )
  }

  cells <- intersect(colnames(mat), colnames(obj))
  if (!length(cells)) {
    cli_abort("No shared sample names between harmonised counts and Seurat object.")
  }
  mat <- mat[, cells, drop = FALSE]
  mat <- as_sparse_counts(mat)

  meta <- tryCatch(
    obj[[]][cells, , drop = FALSE],
    error = function(e) NULL
  )
  if (!is.null(meta) && nrow(meta) > 0L) {
    rownames(meta) <- cells
  }

  create_sobj <- getExportedValue(seurat_ns, "CreateSeuratObject")
  create_args <- list(counts = mat, meta.data = meta)
  if ("assay" %in% names(formals(create_sobj))) {
    create_args$assay <- assay
  }
  obj <- do.call(create_sobj, create_args)
  restore_seurat_cell_names(obj, cells, seurat_ns)
}

#' @keywords internal
#' @noRd
rebuild_se_counts <- function(se, mat, assay = NULL) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    cli::cli_abort("Install SummarizedExperiment to rebuild a SummarizedExperiment.")
  }

  assay_names <- SummarizedExperiment::assayNames(se)
  if (!length(assay_names)) {
    cli_abort("SummarizedExperiment has no assays.")
  }
  if (is.null(assay)) {
    assay <- if ("counts" %in% assay_names) "counts" else assay_names[[1]]
  }
  if (!assay %in% assay_names) {
    cli_abort("Assay `{assay}` not found. Available: {assay_names}.")
  }

  cells <- intersect(colnames(mat), colnames(se))
  if (!length(cells)) {
    cli_abort("No shared sample names between harmonised counts and SummarizedExperiment.")
  }
  mat <- mat[, cells, drop = FALSE]

  rd <- tryCatch(
    as.data.frame(SummarizedExperiment::rowData(se)),
    error = function(e) data.frame()
  )
  if (nrow(rd) && !is.null(rownames(rd))) {
    rd <- rd[rownames(mat), , drop = FALSE]
    rownames(rd) <- rownames(mat)
  } else {
    rd <- data.frame(row.names = rownames(mat))
  }

  cd <- SummarizedExperiment::colData(se)[cells, , drop = FALSE]
  SummarizedExperiment::SummarizedExperiment(
    assays = setNames(list(mat), assay),
    rowData = rd,
    colData = cd
  )
}

#' @keywords internal
#' @noRd
rebuild_aligned_counts <- function(aligned, mat) {
  if (is_aligned_list(aligned)) {
    aligned$counts <- mat
    if (!is.null(aligned$shared_features)) {
      aligned$shared_features <- intersect(aligned$shared_features, rownames(mat))
    }
    return(aligned)
  }
  attach_aligned_metadata(
    mat,
    c(aligned_fields(aligned), list(counts = mat))
  )
}

#' @keywords internal
#' @noRd
rebuild_scaled_container <- function(x, mat, sample_metadata, assay = NULL) {
  cls <- count_input_class(x)
  if (cls == "aligned") {
    cls <- "matrix"
  }
  meta_df <- as.data.frame(sample_metadata)
  rownames(meta_df) <- meta_df$sample_id

  out <- switch(
    cls,
    Seurat = rebuild_seurat_scaled(x, mat, meta_df, assay = assay),
    SummarizedExperiment = rebuild_se_scaled(x, mat, meta_df, assay = assay),
    matrix = {
      storage.mode(mat) <- "double"
      mat
    },
    cli::cli_abort(
      "Don't know how to rebuild class `{paste(class(x), collapse = ', ')}`."
    )
  )
  out
}

#' @keywords internal
#' @noRd
rebuild_seurat_scaled <- function(obj, mat, sample_metadata, assay = NULL) {
  seurat_ns <- if (requireNamespace("SeuratObject", quietly = TRUE)) {
    "SeuratObject"
  } else if (requireNamespace("Seurat", quietly = TRUE)) {
    "Seurat"
  } else {
    cli::cli_abort("Install Seurat or SeuratObject to rebuild a Seurat object.")
  }

  if (is.null(assay)) {
    assay <- tryCatch(
      getExportedValue(seurat_ns, "DefaultAssay")(obj),
      error = function(e) obj@active.assay
    )
  }

  mat <- as_sparse_counts(mat)
  meta <- sample_metadata
  rownames(meta) <- colnames(mat)

  create_sobj <- getExportedValue(seurat_ns, "CreateSeuratObject")
  create_args <- list(counts = mat, meta.data = meta)
  if ("assay" %in% names(formals(create_sobj))) {
    create_args$assay <- assay
  }
  obj <- do.call(create_sobj, create_args)
  restore_seurat_cell_names(obj, colnames(mat), seurat_ns)
}

#' @keywords internal
#' @noRd
rebuild_se_scaled <- function(se, mat, sample_metadata, assay = NULL) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    cli::cli_abort("Install SummarizedExperiment to rebuild a SummarizedExperiment.")
  }

  assay_names <- SummarizedExperiment::assayNames(se)
  if (!length(assay_names)) {
    cli_abort("SummarizedExperiment has no assays.")
  }
  if (is.null(assay)) {
    assay <- if ("counts" %in% assay_names) "counts" else assay_names[[1]]
  }
  if (!assay %in% assay_names) {
    cli_abort("Assay `{assay}` not found. Available: {assay_names}.")
  }

  rd <- tryCatch(
    as.data.frame(SummarizedExperiment::rowData(se)),
    error = function(e) data.frame()
  )
  if (nrow(rd) && !is.null(rownames(rd))) {
    rd <- rd[rownames(mat), , drop = FALSE]
    rownames(rd) <- rownames(mat)
  } else {
    rd <- data.frame(row.names = rownames(mat))
  }

  SummarizedExperiment::SummarizedExperiment(
    assays = setNames(list(mat), assay),
    rowData = rd,
    colData = S4Vectors::DataFrame(sample_metadata)
  )
}

#' Count matrix from common single-cell / bulk containers
#' @keywords internal
#' @noRd
extract_count_matrix <- function(x, assay = NULL, arg_name = "counts") {
  if (is.null(x)) {
    cli::cli_abort("`{arg_name}` is missing.")
  }
  if (is_aligned_list(x)) {
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
  if (is_aligned_list(x)) {
    if (!is.null(x$sample_metadata)) {
      return(as.data.frame(x$sample_metadata))
    }
    x <- x$counts
  }

  sample_df <- aligned_sample_metadata(x)
  if (nrow(sample_df)) {
    return(sample_df)
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
# Public methods — core matrix API + container wrappers
# ---------------------------------------------------------------------------

#' Merge user counts with a one-library reference sample
#'
#' Core matrix helper: keeps genes shared with the reference and appends the
#' reference as one column. Does not scale or compute offsets.
#'
#' @param counts Gene-by-sample numeric matrix (raw counts). Must have gene
#'   rownames and sample colnames.
#' @param reference Named numeric gene vector, or a one-column matrix / data
#'   frame of counts.
#' @param reference_name Column name for the appended reference library.
#'   Default `"hca_reference"`.
#' @return A gene-by-sample matrix of shared genes with the reference column
#'   last. Attribute `shared_features` stores the gene ids used.
#' @seealso [calculate_tmm_offset()], [scale_to_hca_reference()]
#' @export
#' @importFrom cli cli_abort
merge_with_reference_sample <- function(
  counts,
  reference,
  reference_name = "hca_reference"
) {
  if (is.null(counts) || !(is.matrix(counts) || is.data.frame(counts) ||
      inherits(counts, "Matrix"))) {
    cli_abort("`counts` must be a matrix (or coercible to one).")
  }
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    cli_abort("`counts` must have gene rownames and sample colnames.")
  }

  if (is.null(reference_name) || is.na(reference_name) || !nzchar(reference_name)) {
    reference_name <- "hca_reference"
  }
  reference_name <- as.character(reference_name[[1L]])
  if (reference_name %in% colnames(counts)) {
    cli_abort(
      "Sample name `{reference_name}` is already in `counts`. Choose another `reference_name`."
    )
  }

  ref_vec <- as_named_gene_counts(reference, "reference")
  shared <- intersect(rownames(counts), names(ref_vec))
  if (length(shared) < 2L) {
    cli_abort(
      "Need at least 2 shared genes between `counts` and the reference; found {length(shared)}."
    )
  }

  combined <- cbind(
    counts[shared, , drop = FALSE],
    matrix(
      ref_vec[shared],
      ncol = 1L,
      dimnames = list(shared, reference_name)
    )
  )
  storage.mode(combined) <- "double"
  attr(combined, "shared_features") <- shared
  attr(combined, "reference_name") <- reference_name
  combined
}

#' TMM normalisation factors, multipliers, and log offsets
#'
#' Core matrix helper. Runs [edgeR::calcNormFactors()] with
#' `refColumn = reference_name`, then returns:
#' \itemize{
#'   \item `multiplier_j = effectiveSize_ref / effectiveSize_j`
#'   \item `offset_j = log(1 / multiplier_j)` (natural log)
#' }
#' so the reference sample has multiplier 1 and offset 0. This is the same
#' offset convention as `tidybulk::scale_abundance()` followed by
#' `log(1 / multiplier)`, and what [edgeR::mglmOneGroup()] expects.
#'
#' @param counts Gene-by-sample count matrix that already includes the
#'   reference column (e.g. from [merge_with_reference_sample()]).
#' @param reference_name Name of the reference column used as
#'   `refColumn` in [edgeR::calcNormFactors()].
#' @param method Passed to [edgeR::calcNormFactors()]. Default `"TMMwsp"`.
#' @return A list with named numeric vectors `norm_factors`, `library_size`,
#'   `effective_size`, `multiplier`, `offset`, plus scalars `reference_name`
#'   and `method`.
#' @seealso [merge_with_reference_sample()], [scale_to_hca_reference()]
#' @export
#' @importFrom cli cli_abort
calculate_tmm_offset <- function(
  counts,
  reference_name,
  method = "TMMwsp"
) {
  if (is.null(counts) || !(is.matrix(counts) || is.data.frame(counts) ||
      inherits(counts, "Matrix"))) {
    cli_abort("`counts` must be a matrix (or coercible to one).")
  }
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  if (is.null(colnames(counts))) {
    cli_abort("`counts` must have sample colnames.")
  }
  if (missing(reference_name) || is.null(reference_name) ||
      length(reference_name) != 1L || is.na(reference_name) ||
      !nzchar(reference_name)) {
    cli_abort("`reference_name` must be a single non-empty sample name.")
  }
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
#' Wrapper over [merge_with_reference_sample()] and
#' [calculate_tmm_offset()]. Accepts a matrix,
#' `SummarizedExperiment` / `SingleCellExperiment`, or `Seurat` object and
#' returns the same class with the atlas library appended and alignment
#' metadata attached (`hca_offset`, `hca_multiplier`, `sample_role`, …).
#'
#' The atlas reference is usually a cell type string resolved via
#' [load_reference_sample()] / Nectar. A named count vector, one-column
#' matrix/SE, or [get_reference_sample_ready()] download list is also
#' accepted. The reference sits at offset 0.
#'
#' For a matrix-only pipeline without container I/O, call the two core
#' helpers directly (see examples in `examples/savi_adrb2_workflow.R`).
#'
#' @param counts User libraries (matrix, SE/SCE, or Seurat).
#' @param reference Cell type, Nectar download list from
#'   [get_reference_sample_ready()], or a one-library count object.
#' @param reference_name Column name for the bound atlas library.
#'   Default is the atlas sample id, or `"hca_reference"`.
#' @param method Passed to [edgeR::calcNormFactors()]. Default `"TMMwsp"`.
#' @param assay Assay name for SE / Seurat input. Default `"counts"` when present.
#' @param version Nectar version pin, used when `reference` is a cell type.
#' @return An object of the same class as `counts` with the atlas reference
#'   library appended and alignment fields attached.
#' @seealso [merge_with_reference_sample()], [calculate_tmm_offset()],
#'   [estimate_logmu_ql()]
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
  source_container <- counts
  user <- extract_count_matrix(counts, assay = assay, arg_name = "counts")

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

  combined <- merge_with_reference_sample(
    user,
    reference = ref$counts,
    reference_name = reference_name
  )
  shared <- attr(combined, "shared_features")
  scaling <- calculate_tmm_offset(
    combined,
    reference_name = reference_name,
    method = method
  )

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

  user_meta <- extract_sample_metadata(source_container)
  cell_type_val <- if (!is.null(ref$cell_type) && !is.na(ref$cell_type)) {
    ref$cell_type
  } else {
    NA_character_
  }
  sample_metadata <- data.frame(
    sample_id = colnames(combined),
    sample_role = unname(sample_role),
    hca_offset = unname(scaling$offset),
    hca_multiplier = unname(scaling$multiplier),
    hca_reference_name = reference_name,
    hca_cell_type = cell_type_val,
    row.names = colnames(combined),
    stringsAsFactors = FALSE
  )
  if (NROW(user_meta)) {
    sample_metadata <- merge_user_sample_metadata(sample_metadata, user_meta)
  }

  alignment <- list(
    offset = scaling$offset,
    multiplier = scaling$multiplier,
    shared_features = shared,
    reference_name = reference_name,
    sample_role = sample_role,
    sample_metadata = sample_metadata,
    cell_type = cell_type_val
  )

  out <- rebuild_scaled_container(
    source_container,
    combined,
    sample_metadata,
    assay = assay
  )
  attach_aligned_metadata(out, c(alignment, list(counts = combined)))
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

#' Compress offset to a gene-by-sample matrix
#' @keywords internal
#' @noRd
compress_offset_matrix <- function(offset, counts) {
  n_gene <- nrow(counts)
  n_lib <- ncol(counts)
  if (is.matrix(offset)) {
    storage.mode(offset) <- "double"
    return(offset)
  }
  offset <- as.numeric(offset)
  if (length(offset) == 1L) {
    offset <- rep(offset, n_lib)
  }
  if (!is.null(names(offset)) && !is.null(colnames(counts))) {
    offset <- offset[colnames(counts)]
  }
  matrix(offset, nrow = n_gene, ncol = n_lib, byrow = TRUE)
}

#' Fit edgeR QL with prior.count = 0 (absolute log(μ) on offset scale)
#' @keywords internal
#' @noRd
fit_nb_ql <- function(counts, offset, design, robust = TRUE) {
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  design <- as.matrix(design)
  offset_mat <- compress_offset_matrix(offset, counts)

  dge <- edgeR::DGEList(counts = counts)
  dge$offset <- offset_mat
  dge <- edgeR::estimateDisp(dge, design, robust = robust)
  # prior.count = 0 keeps coefficients as absolute log(μ) on the offset scale
  fit <- edgeR::glmQLFit(dge, design, robust = robust, prior.count = 0)

  list(
    fit = fit,
    dispersion = dispersion_from_dge(dge),
    design = design
  )
}

#' Cohort log(μ) and SE for all genes from one edgeR QL fit
#'
#' Fits [edgeR::estimateDisp()] and [edgeR::glmQLFit()] once with
#' `prior.count = 0`, then returns absolute group means (QL coefficients on
#' the supplied offset scale) and QL Wald SEs for every gene.
#'
#' Pass the **user** count matrix and the corresponding offsets. The HCA
#' reference library is used when computing TMM offsets, but should not be
#' included here as an artificial design group.
#'
#' Build the design with ordinary [stats::model.matrix()], for example
#' `model.matrix(~ 0 + Category, sample_metadata)`.
#'
#' @param counts Gene-by-sample numeric count matrix.
#' @param offset Numeric vector (length `ncol(counts)`) or matrix. Typically
#'   `log(1 / multiplier)` from [calculate_tmm_offset()] for the user samples.
#' @param design Numeric design matrix (`nrow` = number of samples). Column
#'   names become the `group` labels in the output.
#' @param robust Passed to [edgeR::estimateDisp()] and [edgeR::glmQLFit()].
#' @return A data frame with one row per gene × design column: `gene`,
#'   `group`, `n`, `log_mu`, `mu`, `se`, `df`, `dispersion`.
#' @seealso [calculate_tmm_offset()], [welch_test_means()]
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

  # QL Wald SE for cell-means columns: sqrt(s2.post / sum(w_j))
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
