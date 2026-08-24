#!/usr/bin/env Rscript
# Sum raw integer counts to patient x sample type x broad cell type and produce
# tumour-epithelial, paired-normal epithelial, broad, and TM4SF1 summaries.

suppressPackageStartupMessages({
  library(Matrix)
  library(SummarizedExperiment)
  library(S4Vectors)
})

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
source(file.path(.script_dir, "utils.R"))

paths <- case_paths()
logf <- start_log(paths, "04_pseudobulk_tm4sf1")

# Load only the lean raw-count object here.  The processed Seurat object contains
# the same counts plus multiple normalized/scaled layers and embeddings, which
# needlessly raises peak memory for a purely count-summing operation.
qc <- readRDS(file.path(paths$processed, "GSE293189_qc_counts.rds"))
counts <- qc$counts
md <- utils::read.delim(
  gzfile(file.path(paths$metadata, "GSE293189_cell_metadata.tsv.gz")),
  stringsAsFactors = FALSE, check.names = FALSE
)
rownames(md) <- md$cell_id
md <- md[colnames(counts), , drop = FALSE]
if (anyNA(md$cell_id) || !identical(md$cell_id, colnames(counts))) {
  stop("Cell annotations do not align to the raw count columns")
}
feature_map <- qc$feature_map
rownames(feature_map) <- feature_map$gene_symbol
feature_map <- feature_map[rownames(counts), , drop = FALSE]

aggregate_counts <- function(counts, md, group_cols, prefix) {
  if (!ncol(counts) || !nrow(md)) stop("Cannot pseudobulk an empty subset")
  id_parts <- lapply(md[, group_cols, drop = FALSE], function(x) {
    x <- as.character(x)
    x[is.na(x) | !nzchar(x)] <- "NA"
    gsub("[^A-Za-z0-9.-]+", "-", x)
  })
  ids <- do.call(paste, c(id_parts, sep = "__"))
  levels_in_order <- unique(ids)
  f <- factor(ids, levels = levels_in_order)
  design <- sparseMatrix(
    i = seq_along(f), j = as.integer(f), x = 1,
    dims = c(length(f), nlevels(f)),
    dimnames = list(rownames(md), paste0(prefix, "__", levels(f)))
  )
  pb <- counts %*% design
  pb <- as(pb, "dgCMatrix")
  storage.mode(pb@x) <- "double"
  first <- match(levels(f), ids)
  cd <- md[first, group_cols, drop = FALSE]
  cd$pseudobulk_id <- colnames(pb)
  cd$n_cells <- as.integer(table(f)[levels(f)])
  cd$library_size <- as.numeric(Matrix::colSums(pb))
  rownames(cd) <- cd$pseudobulk_id
  list(counts = pb, coldata = cd, cell_group = f)
}

make_se <- function(pb, feature_map, metadata_list) {
  se <- SummarizedExperiment(
    assays = list(counts = pb$counts),
    rowData = DataFrame(feature_map),
    colData = DataFrame(pb$coldata)
  )
  metadata(se) <- metadata_list
  se
}

expected_retained <- c(paste0("UC0", 1:3), paste0("VAR0", 1:9))
is_tumour_epi <-
  md$sample_type == "tumour" & md$broad_cell_type == "epithelial" &
  md$patient_id %in% expected_retained
tumour_pb <- aggregate_counts(
  counts[, is_tumour_epi, drop = FALSE], md[is_tumour_epi, , drop = FALSE],
  c("patient_id", "original_patient_id", "paper_id", "histology", "histology_group", "sample_type", "broad_cell_type"),
  "tumour_epithelial"
)
tumour_pb$coldata$reconstructed_n150_threshold_met <- tumour_pb$coldata$n_cells >= 150L
tumour_pb$coldata$published_cohort_member <- tumour_pb$coldata$patient_id %in% expected_retained

source_xlsx <- file.path(paths$raw_paper, "41467_2025_59888_MOESM4_ESM.xlsx")
published <- as.data.frame(readxl::read_excel(source_xlsx, sheet = "Figure1B"))
names(published) <- c("patient_id", "published_tumour_epithelial_cells")
published$published_tumour_epithelial_cells <- as.integer(published$published_tumour_epithelial_cells)
if (sum(published$published_tumour_epithelial_cells) != 8553L) {
  stop("Paper source-data tumour cell counts no longer sum to 8,553")
}
pub_idx <- match(tumour_pb$coldata$patient_id, published$patient_id)
tumour_pb$coldata$published_tumour_epithelial_cells <-
  published$published_tumour_epithelial_cells[pub_idx]
tumour_pb$coldata$cell_count_delta_vs_published <-
  tumour_pb$coldata$n_cells - tumour_pb$coldata$published_tumour_epithelial_cells

tumour_se <- make_se(
  tumour_pb, feature_map,
  list(
    biological_unit = "patient x tumour x epithelial",
    aggregation = "sum of raw integer UMI counts across cells, technical replicates, and tumour pieces",
    cohort_policy = "all 12 published cohort patients retained; reconstructed >=150 flag in colData",
    published_source = "Nature Communications source data, Figure1B"
  )
)
saveRDS(
  tumour_se,
  file.path(paths$processed, "pseudobulk_tumour_epithelial_counts.rds"),
  compress = "gzip"
)

is_normal_epi <- md$sample_type == "paired_normal" & md$broad_cell_type == "epithelial"
normal_pb <- aggregate_counts(
  counts[, is_normal_epi, drop = FALSE], md[is_normal_epi, , drop = FALSE],
  c("patient_id", "original_patient_id", "paper_id", "histology", "histology_group", "sample_type", "broad_cell_type"),
  "paired_normal_epithelial"
)
normal_pb$coldata$sufficient_cells_ge20 <- normal_pb$coldata$n_cells >= 20L
normal_se <- make_se(
  normal_pb, feature_map,
  list(
    biological_unit = "patient x paired_normal x epithelial",
    aggregation = "sum of raw integer UMI counts across cells and libraries",
    tissue_note = "paired/adjacent normal bladder; not assumed fully healthy"
  )
)
saveRDS(
  normal_se,
  file.path(paths$processed, "pseudobulk_paired_normal_epithelial_counts.rds"),
  compress = "gzip"
)

is_annotated_broad <- md$broad_cell_type %in% c("epithelial", "endothelial", "stromal", "immune")
counts_broad <- counts[, is_annotated_broad, drop = FALSE]
md_broad <- md[is_annotated_broad, , drop = FALSE]
all_pb <- aggregate_counts(
  counts_broad, md_broad,
  c("patient_id", "original_patient_id", "paper_id", "histology", "histology_group", "sample_type", "broad_cell_type"),
  "broad"
)
all_se <- make_se(
  all_pb, feature_map,
  list(
    biological_unit = "patient x sample_type x harmonised broad_cell_type",
    aggregation = "sum of raw integer UMI counts across cells, technical replicates, and pieces",
    broad_cell_types = c("epithelial", "endothelial", "stromal", "immune"),
    scope = "annotated bladder cohort; paraganglioma cells remain in raw/QC objects but were outside paper bladder analysis"
  )
)
saveRDS(
  all_se,
  file.path(paths$processed, "pseudobulk_all_broad_celltypes_counts.rds"),
  compress = "gzip"
)

all_cd <- as.data.frame(colData(all_se))
all_cd$matrix <- "pseudobulk_all_broad_celltypes_counts.rds"
write_tsv(all_cd, file.path(paths$metadata, "pseudobulk_sample_metadata.tsv"))
write_tsv(as.data.frame(colData(tumour_se)), file.path(paths$metadata, "pseudobulk_tumour_epithelial_sample_metadata.tsv"))
write_tsv(as.data.frame(colData(normal_se)), file.path(paths$metadata, "pseudobulk_paired_normal_epithelial_sample_metadata.tsv"))

tm_gene <- "TM4SF1"
tm_ens <- feature_map[tm_gene, "ensembl_id"]
if (tm_ens != "ENSG00000169908") stop("TM4SF1 Ensembl identity mismatch")
tm_counts <- as.numeric(assay(all_se, "counts")[tm_gene, ])
tm_summary <- all_cd
tm_summary$gene_symbol <- tm_gene
tm_summary$ensembl_id <- tm_ens
tm_summary$raw_count <- tm_counts
tm_summary$CPM <- 1e6 * tm_summary$raw_count / pmax(tm_summary$library_size, 1)

tm_cell_counts <- as.numeric(counts_broad[tm_gene, ])
all_group_ids <- levels(all_pb$cell_group)
tm_summary$detection_rate <- vapply(
  all_group_ids,
  function(z) mean(tm_cell_counts[all_pb$cell_group == z] > 0),
  numeric(1)
)
tm_summary$n_detected_cells <- as.integer(round(tm_summary$detection_rate * tm_summary$n_cells))
tm_summary$paired_normal_tissue_note <- ifelse(
  tm_summary$sample_type == "paired_normal",
  "paired/adjacent normal bladder; not fully healthy",
  NA_character_
)
tm_summary$published_tumour_cohort_member <-
  tm_summary$sample_type == "tumour" & tm_summary$broad_cell_type == "epithelial" &
  tm_summary$patient_id %in% expected_retained
tm_summary$reconstructed_n150_threshold_met <- ifelse(
  tm_summary$published_tumour_cohort_member,
  tm_summary$n_cells >= 150L,
  NA
)
write_tsv(tm_summary, file.path(paths$metadata, "TM4SF1_patient_summary.tsv"))

group_key <- interaction(
  tm_summary$sample_type, tm_summary$broad_cell_type, tm_summary$histology_group,
  drop = TRUE, lex.order = TRUE
)
group_summary <- do.call(rbind, lapply(split(seq_len(nrow(tm_summary)), group_key), function(ii) {
  z <- tm_summary[ii, , drop = FALSE]
  data.frame(
    sample_type = z$sample_type[[1L]],
    broad_cell_type = z$broad_cell_type[[1L]],
    histology_group = z$histology_group[[1L]],
    n_patients = length(unique(z$patient_id)),
    total_raw_TM4SF1_count = sum(z$raw_count),
    median_patient_CPM = median(z$CPM),
    mean_patient_CPM = mean(z$CPM),
    median_cell_detection_rate = median(z$detection_rate),
    stringsAsFactors = FALSE
  )
}))
write_tsv(group_summary, file.path(paths$metadata, "TM4SF1_group_summary.tsv"))

tumour_cohort_tm <- tm_summary[
  tm_summary$sample_type == "tumour" &
    tm_summary$broad_cell_type == "epithelial" &
    tm_summary$patient_id %in% expected_retained,
  , drop = FALSE
]
write_tsv(
  tumour_cohort_tm,
  file.path(paths$metadata, "TM4SF1_published_tumour_cohort_by_patient.tsv")
)
tumour_cohort_group <- do.call(rbind, lapply(
  split(seq_len(nrow(tumour_cohort_tm)), tumour_cohort_tm$histology_group),
  function(ii) {
    z <- tumour_cohort_tm[ii, , drop = FALSE]
    data.frame(
      histology_group = z$histology_group[[1L]],
      n_published_patients = length(unique(z$patient_id)),
      n_reconstructed_patients_ge150 = sum(z$reconstructed_n150_threshold_met),
      total_cells = sum(z$n_cells),
      total_raw_TM4SF1_count = sum(z$raw_count),
      pooled_library_size = sum(z$library_size),
      pooled_CPM = 1e6 * sum(z$raw_count) / sum(z$library_size),
      median_patient_CPM = median(z$CPM),
      pooled_detection_rate = sum(z$n_detected_cells) / sum(z$n_cells),
      stringsAsFactors = FALSE
    )
  }
))
write_tsv(
  tumour_cohort_group,
  file.path(paths$metadata, "TM4SF1_published_tumour_cohort_group_summary.tsv")
)

population_summary <- aggregate(
  cbind(raw_count, library_size, n_cells, n_detected_cells) ~ sample_type + broad_cell_type,
  data = tm_summary, FUN = sum
)
population_summary$pooled_CPM <-
  1e6 * population_summary$raw_count / pmax(population_summary$library_size, 1)
population_summary$pooled_detection_rate <-
  population_summary$n_detected_cells / pmax(population_summary$n_cells, 1)
population_summary$share_of_TM4SF1_counts_within_sample_type <- ave(
  population_summary$raw_count, population_summary$sample_type,
  FUN = function(x) x / pmax(sum(x), 1)
)
write_tsv(population_summary, file.path(paths$metadata, "TM4SF1_population_source_summary.tsv"))

if (any(assay(all_se, "counts")@x != floor(assay(all_se, "counts")@x))) {
  stop("Pseudobulk counts are not integer-valued")
}
if (sum(assay(all_se, "counts")) != sum(counts_broad)) {
  stop("Broad pseudobulk UMI sum does not match annotated bladder single-cell counts")
}
if (sum(colData(tumour_se)$n_cells) != sum(is_tumour_epi)) {
  stop("Tumour epithelial pseudobulk cell accounting failed")
}

comparison <- merge(
  published,
  as.data.frame(colData(tumour_se))[, c("patient_id", "n_cells", "reconstructed_n150_threshold_met")],
  by = "patient_id", all = TRUE
)
comparison$delta <- comparison$n_cells - comparison$published_tumour_epithelial_cells
comparison$percent_delta <- 100 * comparison$delta / comparison$published_tumour_epithelial_cells
write_tsv(comparison, file.path(paths$metadata, "published_vs_reconstructed_tumour_cells.tsv"))

case_log(
  "Pseudobulk complete: ", ncol(tumour_se), " published tumour patients, ",
  ncol(normal_se), " paired-normal epithelial patients, ", ncol(all_se),
  " patient/sample-type/broad-cell-type profiles", file = logf
)
