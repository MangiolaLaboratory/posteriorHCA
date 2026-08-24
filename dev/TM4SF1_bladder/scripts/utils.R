# Shared helpers for the GSE293189 TM4SF1 bladder-cancer case study.

suppressPackageStartupMessages({
  library(Matrix)
})

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

case_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(file_arg) == 1L && nzchar(file_arg)) {
    return(normalizePath(file.path(dirname(file_arg), ".."), mustWork = TRUE))
  }
  wd <- normalizePath(getwd(), mustWork = TRUE)
  if (file.exists(file.path(wd, "scripts", "utils.R"))) return(wd)
  if (file.exists(file.path(wd, "utils.R"))) {
    return(normalizePath(file.path(wd, ".."), mustWork = TRUE))
  }
  stop("Cannot locate TM4SF1_bladder case-study directory")
}

case_paths <- function(root = case_dir()) {
  p <- list(
    root = root,
    scripts = file.path(root, "scripts"),
    raw = file.path(root, "data", "raw"),
    raw_dge = file.path(root, "data", "raw", "dge"),
    raw_paper = file.path(root, "data", "raw", "paper"),
    author_code = file.path(root, "data", "raw", "Histological_Variant_Bladder_Cancer_Analysis"),
    metadata = file.path(root, "data", "metadata"),
    processed = file.path(root, "data", "processed"),
    logs = file.path(root, "logs")
  )
  invisible(lapply(p, dir.create, recursive = TRUE, showWarnings = FALSE))
  p
}

case_log <- function(..., file = NULL) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  message(line)
  if (!is.null(file)) cat(line, "\n", file = file, append = TRUE)
  invisible(line)
}

start_log <- function(paths, name) {
  out <- file.path(paths$logs, paste0(name, ".log"))
  cat("", file = out)
  case_log("Log file: ", out, file = out)
  out
}

is_valid_file <- function(path, expected_bytes = NA_real_, min_bytes = 1) {
  if (!file.exists(path)) return(FALSE)
  size <- file.info(path)$size
  if (is.na(size) || size < min_bytes) return(FALSE)
  is.na(expected_bytes) || expected_bytes <= 0 || size == expected_bytes
}

download_url <- function(url, dest, expected_bytes = NA_real_, log_file = NULL) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (is_valid_file(dest, expected_bytes)) {
    case_log("Skip existing: ", dest, file = log_file)
    return(invisible(dest))
  }
  partial <- paste0(dest, ".partial")
  args <- c(
    "-q", "-c", "--tries=12", "--timeout=60", "--read-timeout=180",
    "--retry-connrefused", "-O", partial, url
  )
  case_log("Downloading ", url, file = log_file)
  status <- system2("wget", args)
  if (status != 0L || !is_valid_file(partial, expected_bytes)) {
    stop("Download failed or size mismatch: ", url)
  }
  if (file.exists(dest)) unlink(dest)
  if (!file.rename(partial, dest)) stop("Could not finalize download: ", dest)
  case_log("Downloaded ", dest, " (", file.info(dest)$size, " bytes)", file = log_file)
  invisible(dest)
}

sha256_file <- function(path) {
  unname(tools::md5sum(path))
}

parse_soft_samples <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  starts <- which(grepl("^\\^SAMPLE = ", lines))
  ends <- c(starts[-1L] - 1L, length(lines))
  one <- function(i) {
    block <- lines[starts[[i]]:ends[[i]]]
    value <- function(prefix) {
      hit <- block[grepl(prefix, block)]
      if (!length(hit)) return(NA_character_)
      sub(paste0("^", prefix, " = "), "", hit[[1L]])
    }
    relations <- sub("^!Sample_relation = ", "", block[grepl("^!Sample_relation = ", block)])
    supplement <- sub(
      "^!Sample_supplementary_file_[0-9]+ = ", "",
      block[grepl("^!Sample_supplementary_file_", block)]
    )
    data.frame(
      GSM = sub("^\\^SAMPLE = ", "", block[[1L]]),
      geo_title = value("!Sample_title"),
      geo_library_id = sub("^Library name: ", "", value("!Sample_description")),
      biosample = sub(".*biosample/", "", relations[grepl("biosample/", relations)][1L]),
      experiment = sub(".*term=", "", relations[grepl("sra\\?term=", relations)][1L]),
      filename = basename(supplement[[1L]]),
      supplementary_url = supplement[[1L]],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, lapply(seq_along(starts), one))
}

paper_patient_map <- function() {
  data.frame(
    original_patient_id = c(
      "12049", "SG", "HG", "21032", "21222", "PS", "JM", "12041", "FG",
      "20847", "21226", "11734", "21262", "12050", "21217", "PG_M", "PG_S"
    ),
    patient_id = c(
      "UC01", "UC02", "UC03", "UC04", "VAR01", "VAR02", "VAR03", "VAR04",
      "VAR05", "VAR06", "VAR07", "VAR08", "VAR09", "VAR10", "VAR11",
      "PG_M", "PG_S"
    ),
    paper_id = c(
      "UC01", "UC02", "UC03", "UC04", "VAR01", "VAR02", "VAR03", "VAR04",
      "VAR05", "VAR06", "VAR07", "VAR08", "VAR09", "VAR10", "VAR11", NA, NA
    ),
    histology = c(
      "pure UC", "pure UC", "pure UC", "CIS", "micropapillary", "micropapillary",
      "pleomorphic giant cell-like + micropapillary", "nested", "nested",
      "lymphoepithelioma-like", "squamous differentiation", "plasmacytoid",
      "small cell + micropapillary", "pure squamous",
      "squamous + sarcomatoid differentiation", "paraganglioma", "paraganglioma"
    ),
    histology_group = c(rep("UC", 4), rep("HV", 11), "other", "other"),
    tumour_cohort_status = c(
      rep("published_retained", 3), "published_excluded_lt150",
      rep("published_retained", 9), rep("published_excluded_lt150", 2),
      rep("not_in_published_15_tumour_cohort", 2)
    ),
    stringsAsFactors = FALSE
  )
}

patient_from_library <- function(x) {
  keys <- c(
    "20847", "21032", "21217", "21222", "21226", "21262", "11734", "12041",
    "12049", "12050", "PG_M", "PG_S", "FG", "HG", "JM", "PS", "SG"
  )
  vapply(x, function(z) {
    if (grepl("PG_M", z, fixed = TRUE)) return("PG_M")
    if (grepl("PG_S", z, fixed = TRUE)) return("PG_S")
    hit <- keys[vapply(keys, grepl, logical(1), x = z, fixed = TRUE)]
    if (!length(hit)) NA_character_ else hit[[1L]]
  }, character(1))
}

patient_from_title <- function(x) sub("^Patient ([^ ]+).*$", "\\1", x)

sample_fields_from_geo <- function(title) {
  sample_type <- ifelse(grepl("Paired Normal", title), "paired_normal", "tumour")
  tumour_piece <- ifelse(
    sample_type == "paired_normal",
    sub(".*Paired Normal([0-9]*).*$", "N\\1", title),
    sub(".*Tumor([0-9]*).*$", "T\\1", title)
  )
  tumour_piece[tumour_piece %in% c("N", "T")] <- tumour_piece[tumour_piece %in% c("N", "T")]
  technical_rep <- rep(1L, length(title))
  has_rep <- grepl("Replica[0-9]+", title)
  technical_rep[has_rep] <- as.integer(
    sub(".*Replica([0-9]+).*$", "\\1", title[has_rep])
  )
  data.frame(sample_type, tumour_piece, technical_rep, stringsAsFactors = FALSE)
}

sample_fields_from_library <- function(library_id) {
  normal <- grepl(
    "20847.*_N|21032_N|21217_N|21222_N|21226_N|FG_N|JMx_N|PG_MN|PG_SN|12050N",
    library_id
  )
  sample_type <- ifelse(normal, "paired_normal", "tumour")
  piece <- rep(NA_character_, length(library_id))
  for (i in seq_along(library_id)) {
    z <- library_id[[i]]
    if (grepl("11734", z)) piece[[i]] <- "T"
    if (grepl("12041", z)) piece[[i]] <- "T"
    if (grepl("12049", z)) piece[[i]] <- "T"
    if (grepl("12050N", z)) piece[[i]] <- "N"
    if (grepl("12050T1", z)) piece[[i]] <- "T1"
    if (grepl("12050T2", z)) piece[[i]] <- "T2"
  }
  data.frame(sample_type, tumour_piece = piece, technical_rep = 1L, stringsAsFactors = FALSE)
}

read_dge_sparse <- function(path, cell_prefix = NULL) {
  tab <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  genes <- as.character(tab[[1L]])
  if (anyDuplicated(genes)) stop("Duplicated gene rows in ", path)
  mat <- as.matrix(tab[, -1L, drop = FALSE])
  storage.mode(mat) <- "integer"
  mat <- Matrix(mat, sparse = TRUE)
  rownames(mat) <- genes
  if (!is.null(cell_prefix)) {
    colnames(mat) <- paste0(cell_prefix, "__", colnames(mat))
  }
  mat
}

get_counts <- function(object) {
  SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
}

align_sparse_rows <- function(mat, target_rows) {
  if (identical(rownames(mat), target_rows)) return(mat)
  triplet <- summary(mat)
  out <- Matrix::sparseMatrix(
    i = match(rownames(mat)[triplet$i], target_rows),
    j = triplet$j,
    x = triplet$x,
    dims = c(length(target_rows), ncol(mat)),
    dimnames = list(target_rows, colnames(mat)),
    giveCsparse = TRUE
  )
  storage.mode(out@x) <- "double"
  out
}

cbind_sparse <- function(mats, target_rows = NULL) {
  if (is.null(target_rows)) {
    target_rows <- sort(unique(unlist(lapply(mats, rownames), use.names = FALSE)))
  }
  mats <- lapply(mats, align_sparse_rows, target_rows = target_rows)
  do.call(cbind, mats)
}

write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", row.names = FALSE, col.names = TRUE,
    quote = FALSE, na = ""
  )
}
