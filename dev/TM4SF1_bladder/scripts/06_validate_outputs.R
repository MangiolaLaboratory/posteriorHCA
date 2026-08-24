#!/usr/bin/env Rscript
# Independent structural and accounting checks for every deliverable class.

suppressPackageStartupMessages({
  library(Matrix)
  library(SummarizedExperiment)
  library(SeuratObject)
})

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
source(file.path(.script_dir, "utils.R"))

paths <- case_paths()
logf <- start_log(paths, "06_validate_outputs")
checks <- list()
add_check <- function(name, observed, expected, pass = identical(observed, expected)) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    observed = paste(observed, collapse = ","),
    expected = paste(expected, collapse = ","),
    pass = isTRUE(pass),
    stringsAsFactors = FALSE
  )
}

sample_manifest <- utils::read.delim(
  file.path(paths$metadata, "GSE293189_sample_manifest.tsv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
add_check("sample_manifest_rows", nrow(sample_manifest), 67L)
add_check("unique_GSMs", length(unique(sample_manifest$GSM)), 67L)
add_check("unique_SRA_runs", length(unique(sample_manifest$Run)), 67L)
add_check(
  "all_DGE_files_exist",
  sum(file.exists(sample_manifest$local_path)), 67L
)

raw <- readRDS(file.path(paths$processed, "GSE293189_raw_counts.rds"))
raw_counts <- assay(raw, "counts")
add_check("raw_matrix_dimension", dim(raw_counts), c(36994L, 56192L))
add_check("raw_counts_integer_valued", all(raw_counts@x == floor(raw_counts@x)), TRUE)
add_check("raw_metadata_alignment", identical(colnames(raw_counts), rownames(colData(raw))), TRUE)
rm(raw, raw_counts)
invisible(gc(FALSE))

qc <- readRDS(file.path(paths$processed, "GSE293189_qc_counts.rds"))
add_check("QC_matrix_dimension", dim(qc$counts), c(36994L, 26402L))
add_check("QC_counts_integer_valued", all(qc$counts@x == floor(qc$counts@x)), TRUE)
add_check("QC_metadata_alignment", identical(colnames(qc$counts), rownames(qc$cell_metadata)), TRUE)
bladder_cells <- qc$cell_metadata$histology_group %in% c("HV", "UC")
expected_bladder_sum <- sum(qc$counts[, bladder_cells, drop = FALSE])
rm(qc)
invisible(gc(FALSE))

obj <- readRDS(file.path(paths$processed, "GSE293189_processed_seurat.rds"))
obj_counts <- LayerData(obj, assay = "RNA", layer = "counts")
add_check("processed_Seurat_dimension", dim(obj_counts), c(36994L, 21782L))
add_check("processed_raw_count_sum_preserved", sum(obj_counts), expected_bladder_sum)
add_check("Seurat_HVG_count", length(VariableFeatures(obj)), 2000L)
add_check("Seurat_PC_count", ncol(Embeddings(obj, "pca")), 100L)
add_check("Seurat_cluster_count", length(unique(obj$seurat_clusters)), 35L)
rm(obj, obj_counts)
invisible(gc(FALSE))

validate_se <- function(filename, expected_dim, expected_cells) {
  x <- readRDS(file.path(paths$processed, filename))
  m <- assay(x, "counts")
  add_check(paste0(filename, "_dimension"), dim(m), expected_dim)
  add_check(
    paste0(filename, "_integer_valued"),
    all(m@x == floor(m@x)), TRUE
  )
  add_check(
    paste0(filename, "_library_sizes"),
    all(as.numeric(colSums(m)) == colData(x)$library_size), TRUE
  )
  add_check(
    paste0(filename, "_cell_accounting"),
    sum(colData(x)$n_cells), expected_cells
  )
  rm(x, m)
  invisible(gc(FALSE))
}

validate_se("pseudobulk_tumour_epithelial_counts.rds", c(36994L, 12L), 7862L)
validate_se("pseudobulk_paired_normal_epithelial_counts.rds", c(36994L, 7L), 278L)
validate_se("pseudobulk_all_broad_celltypes_counts.rds", c(36994L, 80L), 21782L)

tm <- utils::read.delim(
  file.path(paths$metadata, "TM4SF1_patient_summary.tsv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
add_check("TM4SF1_symbol", unique(tm$gene_symbol), "TM4SF1")
add_check("TM4SF1_Ensembl", unique(tm$ensembl_id), "ENSG00000169908")
add_check("TM4SF1_CPM_finite", all(is.finite(tm$CPM)), TRUE)
add_check("TM4SF1_detection_range", all(tm$detection_rate >= 0 & tm$detection_rate <= 1), TRUE)

result <- do.call(rbind, checks)
write_tsv(result, file.path(paths$metadata, "validation_checks.tsv"))
if (!all(result$pass)) {
  stop(
    "Validation failed: ",
    paste(result$check[!result$pass], collapse = ", ")
  )
}
writeLines(
  paste("PASS", nrow(result), "checks", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  file.path(paths$metadata, "validation_passed.txt")
)
case_log("All ", nrow(result), " validation checks passed", file = logf)

