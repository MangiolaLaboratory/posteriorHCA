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
#' `counts` may be a matrix, `SummarizedExperiment` / `SingleCellExperiment`,
#' or `Seurat` object. The returned object matches the input class. Sample
#' metadata columns `sample_role`, `hca_offset`, `hca_multiplier`,
#' `hca_reference_name`, and `hca_cell_type` describe the atlas alignment.
#'
#' @return An object of the same class as `counts` (matrix, `Seurat`, or
#'   `SummarizedExperiment`) with the atlas reference library appended.
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
#' @inheritParams scale_to_hca_reference
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

  user_meta <- extract_sample_metadata(source_container)
  cell_type_val <- if (!is.null(ref$cell_type) && !is.na(ref$cell_type)) {
    ref$cell_type
  } else {
    NA_character_
  }
  sample_metadata <- data.frame(
    sample_id = colnames(combined),
    sample_role = unname(sample_role),
    hca_offset = unname(offset),
    hca_multiplier = unname(multiplier),
    hca_reference_name = reference_name,
    hca_cell_type = cell_type_val,
    row.names = colnames(combined),
    stringsAsFactors = FALSE
  )
  extra <- setdiff(names(user_meta), names(sample_metadata))
  if (length(extra) && NROW(user_meta)) {
    sample_metadata <- merge_user_sample_metadata(sample_metadata, user_meta)
  }

  alignment <- list(
    offset = offset,
    multiplier = multiplier,
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

#' Subset an aligned count container to one cohort plus the reference sample
#' @keywords internal
#' @noRd
subset_aligned_by_group <- function(counts, resolved_group, cohort_label, assay = NULL) {
  meta <- extract_sample_metadata(counts)
  if (!nrow(meta)) {
    cli::cli_abort("Could not extract sample metadata from `counts`.")
  }

  sample_ids <- sample_ids_from_counts(counts, assay = assay)
  if (length(resolved_group) != length(sample_ids)) {
    cli::cli_abort(
      "`resolved_group` length ({length(resolved_group)}) must match sample count ({length(sample_ids)})."
    )
  }

  keep <- as.character(resolved_group) == cohort_label
  if ("sample_role" %in% names(meta)) {
    role <- as.character(meta$sample_role)
    if (!is.null(rownames(meta))) {
      role <- role[match(sample_ids, rownames(meta))]
    }
    keep <- keep | role == "reference"
  }
  cells <- sample_ids[keep]
  if (!length(cells)) {
    cli::cli_abort("No samples selected for cohort `{cohort_label}`.")
  }

  if (inherits(counts, "Seurat")) {
    return(counts[, cells])
  }
  if (inherits(counts, "SummarizedExperiment")) {
    return(counts[, cells])
  }

  mat <- extract_count_matrix(counts, assay = assay)
  mat <- mat[, cells, drop = FALSE]
  sample_metadata <- meta[cells, , drop = FALSE]

  if (is_aligned_result(counts)) {
    out <- rebuild_scaled_container(counts, mat, sample_metadata)
    fields <- aligned_fields(counts)
    if (!is.null(fields)) {
      fields$counts <- mat
      fields$sample_metadata <- sample_metadata
      if (!is.null(fields$offset)) {
        fields$offset <- fields$offset[cells]
      }
      if (!is.null(fields$multiplier)) {
        fields$multiplier <- fields$multiplier[cells]
      }
      if (!is.null(fields$sample_role)) {
        fields$sample_role <- fields$sample_role[cells]
      }
      out <- attach_aligned_metadata(out, fields)
    }
    return(out)
  }

  mat
}

subset_aligned_cohort <- function(counts, group_col, cohort_label) {
  meta <- extract_sample_metadata(counts)
  if (!nrow(meta)) {
    cli::cli_abort("Could not extract sample metadata from `counts`.")
  }
  if (!group_col %in% names(meta)) {
    cli::cli_abort("Metadata column `{group_col}` not found in `counts`.")
  }

  sample_ids <- rownames(meta)
  if (is.null(sample_ids) || !length(sample_ids)) {
    sample_ids <- sample_ids_from_counts(counts)
  }
  resolved_group <- as.character(meta[[group_col]])
  if (!is.null(rownames(meta))) {
    resolved_group <- resolved_group[match(sample_ids, rownames(meta))]
  }
  subset_aligned_by_group(counts, resolved_group, cohort_label)
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

#' Resolve a group vector from metadata, names, or values
#'
#' For `Seurat` / `SummarizedExperiment` inputs, `group` may be a metadata
#' column name. For matrix inputs, pass a named vector or list keyed by sample
#' id (column name).
#' @keywords internal
#' @noRd
resolve_group <- function(group, counts, n_lib, sample_role = NULL, assay = NULL) {
  if (missing(group) || is.null(group)) {
    cli::cli_abort("`group` is missing.")
  }

  sample_ids <- sample_ids_from_counts(counts, assay = assay)
  if (is.null(sample_ids) || !length(sample_ids)) {
    sample_ids <- seq_len(n_lib)
  }

  if (is.list(group) && !is.data.frame(group)) {
    group <- unlist(group, use.names = TRUE)
  }

  if (is.vector(group) && !is.null(names(group)) && any(nzchar(names(group)))) {
    mapped <- as.character(group)[match(sample_ids, names(group))]
    if (any(is.na(mapped))) {
      missing_ids <- sample_ids[is.na(mapped)]
      cli::cli_abort(
        "Named `group` is missing sample{?s}: {missing_ids}."
      )
    }
    group <- mapped
  } else if (length(group) == 1L && is.character(group)) {
    col_name <- group[[1]]
    meta <- extract_sample_metadata(counts)
    if (NROW(meta) && col_name %in% names(meta)) {
      col_vals <- as.character(meta[[col_name]])
      if (!is.null(rownames(meta))) {
        group <- col_vals[match(sample_ids, rownames(meta))]
      } else {
        group <- col_vals
      }
      if (!is.null(sample_role)) {
        group[is.na(group) & sample_role == "reference"] <- "reference"
      }
      if (any(is.na(group))) {
        cli::cli_abort(
          "Metadata column `{col_name}` could not be matched to all samples."
        )
      }
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
      c(
        "`group` length ({length(group)}) must match the number of libraries ({n_lib}).",
        "i" = "Pass a metadata column name, a named vector keyed by sample id, or a vector in sample order."
      )
    )
  }
  as.character(group)
}

#' Cohort log(mu) at the atlas offset-zero scale
#'
#' Point estimates come from [edgeR::mglmOneGroup()] with the TMM offset
#' applied directly. [edgeR::glmQLFit()] is used only for the QL Wald SE.
#'
#' `counts` may be a matrix, SE/SCE, Seurat object, or the list returned
#' by [scale_to_hca_reference()]. `group` may be:
#' \itemize{
#'   \item a metadata column name for Seurat / SummarizedExperiment inputs,
#'   \item a named vector or list keyed by sample id for matrix inputs, or
#'   \item a vector in sample order.
#' }
#' If `counts` is an aligned list and `group` covers only the user samples,
#' `"reference"` is appended for the atlas library.
#'
#' @param counts Gene-by-sample counts, or an aligned list from
#'   [scale_to_hca_reference()].
#' @param offset Named or unnamed numeric vector, one value per sample.
#'   Taken from `counts$offset` when `counts` is an aligned list.
#' @param group Group labels per sample. See details above.
#' @param genes Optional gene identifiers (symbols or ENSG). Default is all
#'   rows. Values are resolved to canonical ENSG ids. Dispersion is still
#'   estimated from the full `counts` matrix.
#' @param cell_type Optional cell type label stored in the output. Taken from
#'   an aligned list when available.
#' @param version Nectar version pin used when validating genes against the
#'   cell-type universe.
#' @param robust Passed to [edgeR::estimateDisp()] and [edgeR::glmQLFit()].
#' @param assay Assay name for SE / Seurat input.
#' @return A data frame with one row per gene x group: `gene` (ENSG),
#'   `gene_symbol`, `cell_type`, `group`, `n`, `log_mu`, `mu`, `se`, `df`,
#'   `dispersion`.
#' @export
#' @importFrom cli cli_abort
estimate_cohort_logmu <- function(
  counts,
  offset = NULL,
  group,
  genes = NULL,
  cell_type = NULL,
  version = "latest",
  robust = TRUE,
  assay = NULL
) {
  sample_role <- NULL
  source_obj <- counts
  if (is_aligned_result(counts)) {
    aligned <- aligned_fields(counts)
    if (is.null(offset)) {
      offset <- aligned$offset
    }
    if (is.null(cell_type) && !is.null(aligned$cell_type)) {
      cell_type <- aligned$cell_type
    }
    sample_role <- aligned$sample_role
    counts <- aligned$counts
  } else {
    counts <- extract_count_matrix(source_obj, assay = assay, arg_name = "counts")
  }

  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    cli_abort("`counts` must have gene rownames and sample colnames.")
  }
  if (is.null(offset)) {
    fields <- aligned_fields(counts)
    if (!is.null(fields) && !is.null(fields$offset)) {
      offset <- fields$offset
    }
  }
  if (is.null(offset) && inherits(source_obj, "SummarizedExperiment")) {
    cd <- as.data.frame(SummarizedExperiment::colData(source_obj))
    if ("hca_offset" %in% names(cd)) {
      offset <- cd$hca_offset
      names(offset) <- rownames(cd)
    }
  }
  if (is.null(offset) && inherits(source_obj, "Seurat")) {
    meta <- tryCatch(source_obj[[]], error = function(e) NULL)
    if (!is.null(meta) && "hca_offset" %in% names(meta)) {
      offset <- meta$hca_offset
      names(offset) <- rownames(meta)
    }
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

  group <- resolve_group(group, source_obj, ncol(counts), sample_role, assay = assay)
  group <- as.character(group)
  group[is.na(group) | !nzchar(group)] <- "reference"
  samples <- data.frame(
    sample_id = colnames(counts),
    group = factor(group),
    stringsAsFactors = FALSE
  )

  if (is.null(genes)) {
    gene_ids <- rownames(counts)
    gene_symbols <- rep(NA_character_, length(gene_ids))
    names(gene_symbols) <- gene_ids
  } else if (all(genes %in% rownames(counts))) {
    gene_ids <- as.character(genes)
    gene_symbols <- ifelse(is_ensembl_gene_id(gene_ids), NA_character_, gene_ids)
  } else {
    resolved <- resolve_gene(
      genes,
      cell_type = cell_type,
      version = version,
      strict = TRUE
    )
    gene_ids <- unname(resolved)
    gene_symbols <- ifelse(
      is_ensembl_gene_id(names(resolved)),
      NA_character_,
      names(resolved)
    )
    missing <- setdiff(gene_ids, rownames(counts))
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

  out <- lapply(seq_along(gene_ids), function(i) {
    gene <- gene_ids[[i]]
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
        gene_symbol = gene_symbols[[i]],
        cell_type = if (!is.null(cell_type) && !is.na(cell_type)) cell_type else NA_character_,
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
