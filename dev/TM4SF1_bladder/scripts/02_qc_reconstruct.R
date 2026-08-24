#!/usr/bin/env Rscript
# Reconstruct counts from the 67 DGE matrices, apply sequential QC, and remove
# likely doublets with scDblFinder's multi-sample model and library-specific
# thresholds (the installed reproducible alternative to DoubletFinder 2.0.3).

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(S4Vectors)
  library(scDblFinder)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
source(file.path(.script_dir, "utils.R"))

paths <- case_paths()
logf <- start_log(paths, "02_qc_reconstruct")
set.seed(293189)

manifest <- utils::read.delim(
  file.path(paths$metadata, "GSE293189_sample_manifest.tsv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(manifest) != 67L) stop("Sample manifest must contain 67 libraries")

intermediate_dir <- file.path(paths$processed, "intermediate")
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
checkpoint_path <- file.path(intermediate_dir, "GSE293189_pre_doublet_checkpoint.rds")
if (file.exists(checkpoint_path)) {
  case_log("Loading pre-doublet checkpoint", file = logf)
  checkpoint <- readRDS(checkpoint_path)
  raw_counts <- checkpoint$raw_counts
  cell_meta <- checkpoint$cell_meta
  feature_map <- checkpoint$feature_map
  all_genes <- rownames(raw_counts)
} else {
raw_mats <- vector("list", nrow(manifest))
names(raw_mats) <- manifest$library_id
cell_meta <- vector("list", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, , drop = FALSE]
  case_log("Reading/QC ", i, "/", nrow(manifest), ": ", row$library_id, file = logf)
  mat <- read_dge_sparse(row$local_path, cell_prefix = row$library_id)
  raw_mats[[i]] <- mat

  n_count <- Matrix::colSums(mat)
  n_feature <- Matrix::colSums(mat > 0)
  mt <- grepl("^MT-", rownames(mat), ignore.case = FALSE)
  percent_mt <- if (any(mt)) 100 * Matrix::colSums(mat[mt, , drop = FALSE]) / pmax(n_count, 1) else 0

  pass_gene <- n_feature >= 300
  pass_umi <- pass_gene & n_count >= 500
  pass_mito <- pass_umi & percent_mt < 20

  md <- data.frame(
    cell_id = colnames(mat),
    barcode = sub("^[^_]+(?:_[^_]+)*__", "", colnames(mat), perl = TRUE),
    GSM = row$GSM,
    library_id = row$library_id,
    patient_id = row$patient_id,
    original_patient_id = row$original_patient_id,
    paper_id = row$paper_id,
    histology = row$histology,
    histology_group = row$histology_group,
    sample_type = row$sample_type,
    tumour_piece = row$tumour_piece,
    technical_rep = row$technical_rep,
    nCount_RNA_raw = as.numeric(n_count),
    nFeature_RNA_raw = as.integer(n_feature),
    percent_mt_raw = as.numeric(percent_mt),
    pass_min_300_genes = pass_gene,
    pass_min_500_umis = pass_umi,
    pass_max_20pct_mito = pass_mito,
    doublet_score = NA_real_,
    doublet_class = NA_character_,
    doublet_method = ifelse(pass_mito, "pending", "not_eligible_basic_qc"),
    pass_doublet_filter = FALSE,
    stringsAsFactors = FALSE
  )
  rownames(md) <- md$cell_id
  cell_meta[[i]] <- md
  rm(mat)
  invisible(gc(FALSE))
}

cell_meta <- do.call(rbind, cell_meta)
all_genes <- sort(unique(unlist(lapply(raw_mats, rownames), use.names = FALSE)))
case_log("Union feature space: ", length(all_genes), " gene symbols", file = logf)

valid_symbols <- intersect(all_genes, keys(org.Hs.eg.db, keytype = "SYMBOL"))
ensembl <- entrez <- rep(NA_character_, length(all_genes))
names(ensembl) <- names(entrez) <- all_genes
ensembl[valid_symbols] <- mapIds(
  org.Hs.eg.db, valid_symbols, "ENSEMBL", "SYMBOL", multiVals = "first"
)
entrez[valid_symbols] <- mapIds(
  org.Hs.eg.db, valid_symbols, "ENTREZID", "SYMBOL", multiVals = "first"
)
feature_map <- data.frame(
  gene_symbol = all_genes,
  ensembl_id = unname(ensembl),
  entrez_id = unname(entrez),
  is_mitochondrial = grepl("^MT-", all_genes),
  stringsAsFactors = FALSE,
  row.names = all_genes
)
if (feature_map["TM4SF1", "ensembl_id"] != "ENSG00000169908") {
  stop("TM4SF1 did not map to ENSG00000169908")
}
write_tsv(feature_map, file.path(paths$metadata, "GSE293189_feature_map.tsv"))

case_log("Aligning/saving deposited raw count matrix", file = logf)
raw_counts <- cbind_sparse(raw_mats, all_genes)
if (ncol(raw_counts) != nrow(cell_meta)) stop("Raw cell metadata/count mismatch")
cell_meta <- cell_meta[colnames(raw_counts), , drop = FALSE]
saveRDS(
  list(raw_counts = raw_counts, cell_meta = cell_meta, feature_map = feature_map),
  checkpoint_path, compress = FALSE
)
}

basic_cells <- which(cell_meta$pass_max_20pct_mito)
post_qc_per_library <- table(cell_meta$library_id[basic_cells])
eligible_libraries <- names(post_qc_per_library)[post_qc_per_library >= 50L]
eligible_cells <- basic_cells[cell_meta$library_id[basic_cells] %in% eligible_libraries]
ineligible_cells <- setdiff(basic_cells, eligible_cells)
case_log(
  "Calling doublets for ", length(eligible_cells), " basic-QC cells in ",
  length(eligible_libraries), " libraries using scDblFinder artificial-doublet scores and explicit expected-rate calls; ",
  length(ineligible_cells), " cells in <50-cell libraries retained without a call", file = logf
)
set.seed(293189)
dbl_sce <- SingleCellExperiment(
  assays = list(counts = raw_counts[, eligible_cells, drop = FALSE]),
  colData = DataFrame(library_id = cell_meta$library_id[eligible_cells])
)
score_checkpoint <- file.path(intermediate_dir, "GSE293189_scDblFinder_score_table.rds")
if (file.exists(score_checkpoint)) {
  dbl_table <- readRDS(score_checkpoint)
} else {
  dbl_table <- scDblFinder(
    dbl_sce,
    clusters = FALSE,
    samples = "library_id",
    multiSampleMode = "singleModel",
    score = "ratio",
    iter = 1,
    threshold = FALSE,
    returnType = "table",
    verbose = TRUE,
    BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)
  )
  saveRDS(dbl_table, score_checkpoint, compress = "gzip")
}
real_scores <- as.data.frame(dbl_table[dbl_table$type == "real", , drop = FALSE])
score_column <- intersect(c("score", "ratio", "weighted"), names(real_scores))[[1L]]
case_log(
  "Doublet table columns: ", paste(names(real_scores), collapse = ","),
  "; using ", score_column, file = logf
)
score_by_cell <- setNames(as.numeric(real_scores[[score_column]]), rownames(real_scores))
if (!all(colnames(dbl_sce) %in% names(score_by_cell))) {
  stop("scDblFinder score table lost real-cell identifiers")
}
cell_meta$doublet_score[eligible_cells] <- score_by_cell[colnames(dbl_sce)]
cell_meta$doublet_class[eligible_cells] <- "singlet"
for (lib in eligible_libraries) {
  ii <- eligible_cells[cell_meta$library_id[eligible_cells] == lib]
  n <- length(ii)
  expected_rate <- 0.008 * n / 1000
  n_expected <- as.integer(round(n * expected_rate))
  if (n_expected > 0L) {
    called <- head(ii[order(cell_meta$doublet_score[ii], decreasing = TRUE)], n_expected)
    cell_meta$doublet_class[called] <- "doublet"
  }
}
cell_meta$doublet_method[eligible_cells] <-
  paste0(
    "scDblFinder_1.24.10_", score_column,
    "_score_manual_per_library_expected_rate"
  )
cell_meta$doublet_class[ineligible_cells] <- "singlet"
cell_meta$doublet_method[ineligible_cells] <- "not_run_lt50_post_qc_cells_kept"
cell_meta$pass_doublet_filter <-
  cell_meta$pass_max_20pct_mito & cell_meta$doublet_class != "doublet"
cell_meta$pass_doublet_filter[is.na(cell_meta$pass_doublet_filter)] <- FALSE

raw_sce <- SingleCellExperiment(
  assays = list(counts = raw_counts),
  rowData = DataFrame(feature_map),
  colData = DataFrame(cell_meta)
)
metadata(raw_sce)$source <- "GEO GSE293189 per-library DGE matrices"
metadata(raw_sce)$counts_are_raw_integer_umis <- TRUE
metadata(raw_sce)$fastq_downloaded <- FALSE
metadata(raw_sce)$feature_identifier_note <- "GEO DGE rows are gene symbols; Ensembl IDs are mapped in rowData"
saveRDS(raw_sce, file.path(paths$processed, "GSE293189_raw_counts.rds"), compress = "gzip")

qc_counts <- raw_counts[, cell_meta$pass_doublet_filter, drop = FALSE]
qc_meta <- cell_meta[colnames(qc_counts), , drop = FALSE]
saveRDS(
  list(counts = qc_counts, cell_metadata = qc_meta, feature_map = feature_map),
  file.path(paths$processed, "GSE293189_qc_counts.rds"), compress = "gzip"
)
saveRDS(
  cell_meta,
  file.path(paths$processed, "GSE293189_all_cell_qc_metadata.rds"), compress = "gzip"
)

stage_counts <- function(df, group_vars) {
  key <- interaction(df[, group_vars, drop = FALSE], drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(df)), key)
  out <- lapply(groups, function(ii) {
    z <- df[ii, , drop = FALSE]
    cbind(
      z[1L, group_vars, drop = FALSE],
      data.frame(
        n_raw = nrow(z),
        n_after_min_300_genes = sum(z$pass_min_300_genes),
        n_after_min_500_umis = sum(z$pass_min_500_umis),
        n_after_max_20pct_mito = sum(z$pass_max_20pct_mito),
        n_after_doublet_filter = sum(z$pass_doublet_filter),
        n_doublets_removed = sum(z$doublet_class == "doublet", na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    )
  })
  do.call(rbind, out)
}

qc_library <- stage_counts(cell_meta, c("GSM", "library_id", "patient_id", "sample_type"))
qc_patient <- stage_counts(cell_meta, c("patient_id", "original_patient_id", "paper_id", "histology_group"))
qc_overall <- data.frame(
  level = "overall",
  n_raw = nrow(cell_meta),
  n_after_min_300_genes = sum(cell_meta$pass_min_300_genes),
  n_after_min_500_umis = sum(cell_meta$pass_min_500_umis),
  n_after_max_20pct_mito = sum(cell_meta$pass_max_20pct_mito),
  n_after_doublet_filter = sum(cell_meta$pass_doublet_filter),
  n_doublets_removed = sum(cell_meta$doublet_class == "doublet", na.rm = TRUE)
)
write_tsv(qc_library, file.path(paths$metadata, "qc_cell_counts_by_library.tsv"))
write_tsv(qc_patient, file.path(paths$metadata, "qc_cell_counts_by_patient.tsv"))
write_tsv(qc_overall, file.path(paths$metadata, "qc_cell_counts_overall.tsv"))

if (any(qc_counts@x != floor(qc_counts@x))) stop("QC count matrix is not integer-valued")
case_log(
  "QC complete: ", nrow(cell_meta), " raw -> ", ncol(qc_counts),
  " retained cells; ", nrow(feature_map), " union genes", file = logf
)
