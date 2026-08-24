# Shared helpers for the PAK2 / IPF data-preparation pipeline.
# Intentionally does not filter cells, spots, genes, or samples.

suppressPackageStartupMessages({
  library(Matrix)
  library(methods)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

na_chr <- function(n = 1L) rep(NA_character_, n)
na_int <- function(n = 1L) rep(NA_integer_, n)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

pak2_case_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
  if (length(file_arg) == 1L && nzchar(file_arg)) {
    return(normalizePath(file.path(dirname(file_arg), ".."), mustWork = TRUE))
  }
  wd <- normalizePath(getwd(), mustWork = TRUE)
  if (file.exists(file.path(wd, "scripts", "utils.R"))) return(wd)
  if (file.exists(file.path(wd, "utils.R"))) {
    return(normalizePath(file.path(wd, ".."), mustWork = TRUE))
  }
  stop("Cannot locate PAK2_case directory. Run from that folder or via Rscript scripts/*.R")
}

pak2_paths <- function(case_dir = pak2_case_dir()) {
  paths <- list(
    case_dir = case_dir,
    scripts = file.path(case_dir, "scripts"),
    data = file.path(case_dir, "data"),
    raw = file.path(case_dir, "data", "raw"),
    raw_gse285246 = file.path(case_dir, "data", "raw", "GSE285246"),
    raw_gse173896 = file.path(case_dir, "data", "raw", "GSE173896"),
    prepared = file.path(case_dir, "data", "raw", "prepared"),
    metadata = file.path(case_dir, "data", "metadata"),
    seurat_spatial = file.path(case_dir, "data", "seurat", "spatial"),
    seurat_scrna = file.path(case_dir, "data", "seurat", "scrna"),
    pb_spatial = file.path(case_dir, "data", "pseudobulk", "spatial"),
    pb_scrna = file.path(case_dir, "data", "pseudobulk", "scrna"),
    logs = file.path(case_dir, "logs")
  )
  invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))
  paths
}

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

pak2_log <- function(..., .log_file = NULL) {
  msg <- paste0(...)
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", stamp, "] ", msg)
  message(line)
  if (!is.null(.log_file)) {
    cat(line, "\n", file = .log_file, append = TRUE)
  }
  invisible(line)
}

start_log <- function(paths, name) {
  f <- file.path(paths$logs, paste0(name, ".log"))
  cat("", file = f, append = TRUE)
  pak2_log("Log file: ", f, .log_file = f)
  f
}

# ---------------------------------------------------------------------------
# Downloads
# ---------------------------------------------------------------------------

ftp_to_https <- function(url) {
  sub("^ftp://ftp\\.ncbi\\.nlm\\.nih\\.gov/", "https://ftp.ncbi.nlm.nih.gov/", url)
}

is_valid_file <- function(path, expected_size = NA_real_, min_bytes = 1) {
  if (!file.exists(path)) return(FALSE)
  sz <- file.info(path)$size
  if (is.na(sz) || sz < min_bytes) return(FALSE)
  if (!is.na(expected_size) && expected_size > 0) {
    return(as.numeric(sz) == as.numeric(expected_size))
  }
  TRUE
}

download_url <- function(url, destfile, expected_size = NA_real_, overwrite = FALSE,
                         log_file = NULL, tries = 12L) {
  url <- ftp_to_https(url)
  dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
  if (!overwrite && is_valid_file(destfile, expected_size = expected_size)) {
    pak2_log("Skip existing: ", destfile, " (", file.info(destfile)$size, " bytes)",
             .log_file = log_file)
    return(invisible(destfile))
  }
  tmp <- paste0(destfile, ".partial")
  if (file.exists(destfile) && (overwrite || !is_valid_file(destfile, expected_size))) {
    unlink(destfile)
  }
  pak2_log("Downloading ", url, " -> ", destfile, .log_file = log_file)
  cmd <- sprintf(
    "wget -c --tries=%d --timeout=60 --read-timeout=180 --retry-connrefused -O %s %s",
    as.integer(tries),
    shQuote(tmp),
    shQuote(url)
  )
  status <- system(cmd)
  if (status != 0L || !file.exists(tmp) || is.na(file.info(tmp)$size) || file.info(tmp)$size < 1) {
    stop("Download failed: ", url)
  }
  if (!is.na(expected_size) && expected_size > 0 && file.info(tmp)$size != expected_size) {
    stop(
      "Size mismatch for ", destfile, ": got ", file.info(tmp)$size,
      " expected ", expected_size
    )
  }
  if (file.exists(destfile)) unlink(destfile)
  file.rename(tmp, destfile)
  pak2_log("OK ", destfile, " (", file.info(destfile)$size, " bytes)", .log_file = log_file)
  invisible(destfile)
}

gunzip_copy <- function(src_gz, dest, overwrite = FALSE) {
  if (!overwrite && file.exists(dest) && file.info(dest)$size > 0) return(dest)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (!grepl("\\.gz$", src_gz, ignore.case = TRUE)) {
    if (!file.exists(dest)) file.copy(src_gz, dest, overwrite = TRUE)
    return(dest)
  }
  status <- system2("gzip", c("-dc", src_gz), stdout = dest)
  # gzip -dc writes to stdout; system2 with stdout=dest captures it.
  if (!file.exists(dest) || file.info(dest)$size < 1) {
    stop("gunzip failed: ", src_gz)
  }
  dest
}

link_or_copy <- function(from, to, overwrite = FALSE) {
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(to)) {
    if (!overwrite) return(to)
    unlink(to)
  }
  ok <- file.symlink(normalizePath(from, mustWork = TRUE), to)
  if (!ok || !file.exists(to)) {
    file.copy(from, to, overwrite = TRUE)
  }
  to
}

# ---------------------------------------------------------------------------
# SOFT / GEO parsing
# ---------------------------------------------------------------------------

parse_soft_samples <- function(soft_path) {
  lines <- readLines(soft_path, warn = FALSE)
  sample_starts <- which(grepl("^\\^SAMPLE = ", lines))
  if (!length(sample_starts)) {
    return(data.frame())
  }
  sample_ends <- c(sample_starts[-1] - 1L, length(lines))
  rows <- lapply(seq_along(sample_starts), function(i) {
    block <- lines[sample_starts[i]:sample_ends[i]]
    kv <- function(key) {
      hits <- sub(paste0("^", key, " = "), "", block[grepl(paste0("^", key, " = "), block)])
      if (!length(hits)) NA_character_ else hits
    }
    chars <- kv("!Sample_characteristics_ch1")
    char_map <- list()
    if (!all(is.na(chars))) {
      for (ch in chars) {
        if (grepl(":", ch, fixed = TRUE)) {
          k <- trimws(sub(":.*$", "", ch))
          v <- trimws(sub("^[^:]+:\\s*", "", ch))
          char_map[[k]] <- v
        }
      }
    }
    files <- kv("!Sample_supplementary_file_[0-9]+")
    # The regex above does not work in grepl that way; parse files separately.
    files <- sub(
      "^!Sample_supplementary_file_[0-9]+ = ",
      "",
      block[grepl("^!Sample_supplementary_file_", block)]
    )
    data.frame(
      gsm_accession = sub("^\\^SAMPLE = ", "", block[1]),
      sample_title = kv("!Sample_title")[1],
      source_name = kv("!Sample_source_name_ch1")[1],
      organism = kv("!Sample_organism_ch1")[1],
      platform_id = kv("!Sample_platform_id")[1],
      instrument_model = kv("!Sample_instrument_model")[1],
      library_strategy = kv("!Sample_library_strategy")[1],
      library_source = kv("!Sample_library_source")[1],
      series_id = kv("!Sample_series_id")[1],
      tissue_characteristic = char_map$tissue %||% NA_character_,
      disease_characteristic = char_map$disease %||% NA_character_,
      cell_type_characteristic = char_map[["cell type"]] %||% NA_character_,
      characteristics_raw = paste(chars, collapse = " | "),
      supplementary_files = paste(files, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

parse_geo_filelist <- function(path) {
  raw <- readLines(path, warn = FALSE)
  raw <- raw[!grepl("^#", raw) & nzchar(raw)]
  if (!length(raw) || !grepl("^Archive|^File", raw[1])) {
    stop("Unexpected GEO filelist format: ", path)
  }
  # Tab-separated: type, name, time, size, type2
  parts <- strsplit(raw, "\t", fixed = TRUE)
  keep <- vapply(parts, function(x) x[[1]] == "File" && length(x) >= 4L, logical(1))
  parts <- parts[keep]
  data.frame(
    file_name = vapply(parts, `[[`, character(1), 2L),
    expected_size = as.numeric(vapply(parts, `[[`, character(1), 4L)),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# 10x matrices
# ---------------------------------------------------------------------------

read_features_tsv <- function(path) {
  feat <- utils::read.delim(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c("gene_id", "gene_symbol", "feature_type")[seq_len(max(1L, ncol(utils::read.delim(path, header = FALSE, nrows = 1))))]
  )
  if (ncol(feat) == 1L) {
    names(feat) <- "gene_id"
    feat$gene_symbol <- feat$gene_id
    feat$feature_type <- NA_character_
  } else if (ncol(feat) == 2L) {
    names(feat) <- c("gene_id", "gene_symbol")
    feat$feature_type <- NA_character_
  } else {
    names(feat)[1:3] <- c("gene_id", "gene_symbol", "feature_type")
    feat <- feat[, 1:3, drop = FALSE]
  }
  feat
}

read_10x_counts <- function(matrix_dir, gene_column = 1L) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required to read 10x MTX directories")
  }
  counts <- Seurat::Read10X(
    data.dir = matrix_dir,
    gene.column = gene_column,
    unique.features = TRUE,
    strip.suffix = FALSE
  )
  if (is.list(counts) && !inherits(counts, "Matrix")) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }
  counts
}

prepare_10x_mtx_dir <- function(dest_dir, barcodes, features, matrix, overwrite = FALSE) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  link_or_copy(barcodes, file.path(dest_dir, "barcodes.tsv.gz"), overwrite = overwrite)
  link_or_copy(features, file.path(dest_dir, "features.tsv.gz"), overwrite = overwrite)
  link_or_copy(matrix, file.path(dest_dir, "matrix.mtx.gz"), overwrite = overwrite)
  dest_dir
}

prepare_visium_dir <- function(dest_dir, files, overwrite = FALSE) {
  # files: named list/vector with barcodes, features, matrix, positions, scalefactors,
  # hires, lowres (optional jpg extras ignored for Seurat)
  mtx_dir <- file.path(dest_dir, "filtered_feature_bc_matrix")
  img_dir <- file.path(dest_dir, "spatial")
  prepare_10x_mtx_dir(
    mtx_dir,
    barcodes = files[["barcodes"]],
    features = files[["features"]],
    matrix = files[["matrix"]],
    overwrite = overwrite
  )
  dir.create(img_dir, recursive = TRUE, showWarnings = FALSE)
  gunzip_copy(files[["positions"]], file.path(img_dir, "tissue_positions_list.csv"), overwrite)
  gunzip_copy(files[["scalefactors"]], file.path(img_dir, "scalefactors_json.json"), overwrite)
  gunzip_copy(files[["hires"]], file.path(img_dir, "tissue_hires_image.png"), overwrite)
  gunzip_copy(files[["lowres"]], file.path(img_dir, "tissue_lowres_image.png"), overwrite)
  dest_dir
}

read_tissue_positions <- function(path) {
  first <- readLines(path, n = 1L)
  has_header <- grepl("barcode", first, ignore.case = TRUE)
  pos <- utils::read.csv(path, header = has_header, stringsAsFactors = FALSE)
  if (!has_header) {
    names(pos) <- c(
      "barcode", "in_tissue", "array_row", "array_col",
      "pxl_row_in_fullres", "pxl_col_in_fullres"
    )[seq_len(ncol(pos))]
  }
  names(pos) <- gsub("\\.", "_", names(pos))
  pos
}

# ---------------------------------------------------------------------------
# Counts / Seurat helpers
# ---------------------------------------------------------------------------

get_raw_counts <- function(obj, assay = NULL) {
  if (is.null(assay)) assay <- Seurat::DefaultAssay(obj)
  SeuratObject::LayerData(obj, assay = assay, layer = "counts")
}

feature_table_from_10x <- function(features_path, rownames_used) {
  feat <- read_features_tsv(features_path)
  idx <- match(rownames_used, feat$gene_id)
  if (anyNA(idx)) {
    idx2 <- match(rownames_used, feat$gene_symbol)
    if (!anyNA(idx2) && anyNA(idx)) idx <- idx2
  }
  if (anyNA(idx)) {
    stop(
      "Could not map ", sum(is.na(idx)), " features from matrix rownames to ",
      features_path
    )
  }
  out <- feat[idx, , drop = FALSE]
  rownames(out) <- rownames_used
  out
}

standard_spot_or_cell_meta <- function(n, values) {
  stopifnot(is.list(values))
  as.data.frame(
    lapply(values, function(v) {
      if (length(v) == 1L) rep(v, n) else v
    }),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

# ---------------------------------------------------------------------------
# Pseudobulk aggregation of raw integer counts
# ---------------------------------------------------------------------------

aggregate_counts_by_group <- function(counts, group, group_order = NULL) {
  stopifnot(length(group) == ncol(counts))
  group <- as.character(group)
  if (anyNA(group) || any(!nzchar(group))) {
    stop("group contains NA/empty values; caller must split annotated vs unannotated")
  }
  if (is.null(group_order)) {
    group_order <- sort(unique(group))
  } else {
    group_order <- unique(as.character(group_order))
  }
  j <- match(group, group_order)
  if (anyNA(j)) stop("group values missing from group_order")
  dummy <- Matrix::sparseMatrix(
    i = seq_along(j),
    j = j,
    x = 1,
    dims = c(ncol(counts), length(group_order))
  )
  colnames(dummy) <- group_order
  pb <- counts %*% dummy
  pb <- methods::as(pb, "dgCMatrix")
  if (length(pb@x) && any(abs(pb@x - round(pb@x)) > 1e-6)) {
    stop("Pseudobulk aggregation produced non-integer values")
  }
  if (length(pb@x)) pb@x <- round(pb@x)
  pb
}

make_summarized_experiment <- function(counts, row_data, col_data, metadata = list()) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment is required")
  }
  counts <- methods::as(counts, "dgCMatrix")
  row_data <- S4Vectors::DataFrame(row_data, row.names = rownames(counts))
  col_data <- S4Vectors::DataFrame(col_data, row.names = colnames(counts))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = row_data,
    colData = col_data,
    metadata = metadata
  )
}

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

validate_count_sum <- function(source_counts, pb_counts, label) {
  s1 <- sum(source_counts)
  s2 <- sum(pb_counts)
  if (!isTRUE(all.equal(s1, s2, tolerance = 0, countEQ = TRUE)) && s1 != s2) {
    stop(label, ": count-sum mismatch source=", s1, " pseudobulk=", s2)
  }
  TRUE
}

n_detected_features <- function(counts) {
  as.integer(Matrix::colSums(counts > 0))
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  path
}

save_rds_atomic <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  saveRDS(object, tmp)
  if (file.exists(path)) unlink(path)
  file.rename(tmp, path)
  path
}

looks_like_barcode_column <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(FALSE)
  probe <- x[seq_len(min(200L, length(x)))]
  probe <- sub("_[0-9]+$", "", probe)
  mean(grepl("^[ACGT]{14,}-?[0-9]*$", probe)) > 0.8
}

normalise_barcode <- function(x) {
  x <- as.character(x)
  # Drop a trailing Seurat merge suffix such as "_3", keeping the 10x "-1".
  sub("_[0-9]+$", "", x)
}
