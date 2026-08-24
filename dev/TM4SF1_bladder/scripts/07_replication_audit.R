#!/usr/bin/env Rscript
# Audit the published-vs-reconstructed tumour cohort without mutating canonical
# objects. All outputs are written to results/replication_audit/.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
  library(readxl)
})

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
source(file.path(.script_dir, "utils.R"))

paths <- case_paths()
audit <- file.path(paths$root, "results", "replication_audit")
dir.create(audit, recursive = TRUE, showWarnings = FALSE)
logf <- start_log(paths, "07_replication_audit")
set.seed(293189)

read_meta <- function(name) utils::read.delim(
  file.path(paths$metadata, name), stringsAsFactors = FALSE,
  check.names = FALSE, na.strings = c("", "NA")
)
write_gz_tsv <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

manifest <- read_meta("GSE293189_sample_manifest.tsv")
patients <- read_meta("GSE293189_patient_manifest.tsv")
all_md <- utils::read.delim(
  gzfile(file.path(paths$metadata, "GSE293189_cell_metadata.tsv.gz")),
  stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA")
)
obj <- readRDS(file.path(paths$processed, "GSE293189_processed_seurat.rds"))
counts <- LayerData(obj, assay = "RNA", layer = "counts")
norm <- LayerData(obj, assay = "RNA", layer = "data")

source_xlsx <- file.path(paths$raw_paper, "41467_2025_59888_MOESM4_ESM.xlsx")
published <- as.data.frame(read_excel(source_xlsx, sheet = "Figure1B"))
names(published) <- c("patient_id", "published_tumour_epithelial_cells")
published$published_tumour_epithelial_cells <-
  as.integer(published$published_tumour_epithelial_cells)

# Patient-level flow from every deposited barcode through current reconstruction.
patient_ids <- patients$patient_id[patients$histology_group %in% c("HV", "UC")]
flow_one <- function(pid) {
  z <- all_md[all_md$patient_id == pid, , drop = FALSE]
  data.frame(
    patient_id = pid,
    original_patient_id = unique(z$original_patient_id)[1L],
    histology_group = unique(z$histology_group)[1L],
    n_raw = nrow(z),
    n_after_gene = sum(z$pass_min_300_genes, na.rm = TRUE),
    n_after_umi = sum(z$pass_min_500_umis, na.rm = TRUE),
    n_after_mito = sum(z$pass_max_20pct_mito, na.rm = TRUE),
    n_after_doublet = sum(z$pass_doublet_filter, na.rm = TRUE),
    n_processed = sum(z$analysis_inclusion == "processed_paper_cohort", na.rm = TRUE),
    n_broad_epithelial = sum(z$broad_cell_type == "epithelial", na.rm = TRUE),
    n_candidate_tumour_epithelial = sum(
      z$sample_type == "tumour" & z$broad_cell_type == "epithelial", na.rm = TRUE
    ),
    n_final_tumour_epithelial = sum(z$tumour_epithelial_identity, na.rm = TRUE),
    reconstructed_ge150 = sum(z$tumour_epithelial_identity, na.rm = TRUE) >= 150L,
    stringsAsFactors = FALSE
  )
}
cell_flow <- do.call(rbind, lapply(patient_ids, flow_one))
cell_flow <- merge(cell_flow, published, by = "patient_id", all.x = TRUE, sort = FALSE)
cell_flow <- cell_flow[match(patient_ids, cell_flow$patient_id), ]
cell_flow$delta_tumour_epithelial <-
  cell_flow$n_final_tumour_epithelial - cell_flow$published_tumour_epithelial_cells
write_tsv(cell_flow, file.path(audit, "all_patient_cell_flow.tsv"))

uc_flow_wide <- cell_flow[cell_flow$patient_id == "UC01", , drop = FALSE]
uc_flow <- data.frame(
  patient_id = "UC01",
  stage = c(
    "deposited_barcodes", "pass_gene", "pass_umi", "pass_mito",
    "after_doublet", "retained_bladder_cells", "broad_epithelial",
    "candidate_tumour_epithelial", "final_tumour_epithelial"
  ),
  n_cells = as.integer(c(
    uc_flow_wide$n_raw, uc_flow_wide$n_after_gene, uc_flow_wide$n_after_umi,
    uc_flow_wide$n_after_mito, uc_flow_wide$n_after_doublet,
    uc_flow_wide$n_processed, uc_flow_wide$n_broad_epithelial,
    uc_flow_wide$n_candidate_tumour_epithelial,
    uc_flow_wide$n_final_tumour_epithelial
  )),
  stringsAsFactors = FALSE
)
write_tsv(uc_flow, file.path(audit, "uc01_cell_flow.tsv"))

# Library-level audit includes both the DGE-filename assignment and the
# conflicting GEO-title assignment, plus the original author-code identifier.
lib_qc <- read_meta("qc_cell_counts_by_library.tsv")
uc_rows <- unique(c(
  which(manifest$patient_id == "UC01"),
  which(manifest$geo_original_patient_id == "12049"),
  which(manifest$geo_library_id == "PA_U67-12049_Pool_1_2_3_S20_L001")
))
uc_lib <- manifest[uc_rows, c(
  "GSM", "library_id", "patient_id", "original_patient_id", "sample_type",
  "filename", "geo_title", "geo_library_id", "geo_original_patient_id",
  "mapping_status", "Run", "local_path"
)]
uc_lib <- merge(
  uc_lib, lib_qc[, c(
    "GSM", "n_raw", "n_after_min_300_genes", "n_after_min_500_umis",
    "n_after_max_20pct_mito", "n_after_doublet_filter", "n_doublets_removed"
  )], by = "GSM", all.x = TRUE, sort = FALSE
)
uc_lib$uc01_claim_basis <- ifelse(
  uc_lib$patient_id == "UC01", "DGE_filename_12049",
  "GEO_title_or_description_12049"
)
uc_lib$downloaded <- file.exists(uc_lib$local_path)
uc_lib$assigned_to_current_UC01 <- uc_lib$patient_id == "UC01"
missing_12923 <- as.list(rep(NA, ncol(uc_lib)))
names(missing_12923) <- names(uc_lib)
missing_12923$GSM <- "not_deposited"
missing_12923$library_id <- "12923"
missing_12923$patient_id <- "UC01"
missing_12923$original_patient_id <- "12923"
missing_12923$sample_type <- "tumour"
missing_12923$mapping_status <- "original_public_code_and_downstream_scripts"
missing_12923$uc01_claim_basis <- "authors_pre2026_code_12923"
missing_12923$downloaded <- FALSE
missing_12923$assigned_to_current_UC01 <- FALSE
uc_lib <- rbind(uc_lib, as.data.frame(missing_12923, stringsAsFactors = FALSE))
write_tsv(uc_lib, file.path(audit, "uc01_library_audit.tsv"))

# The paper source workbook contains cell identifiers for all 8,553 reported
# tumour cells in its TM4SF1 worksheet, plus a separate 458-cell Plasma group.
source_cells <- as.data.frame(read_excel(source_xlsx, sheet = "Figure5B"))
source_cells$base_barcode <- ifelse(
  grepl("^[ACGT]{12}_", source_cells$Barcode),
  sub("_.*$", "", source_cells$Barcode),
  NA_character_
)
dge_barcodes <- lapply(manifest$local_path, function(f) {
  strsplit(readLines(gzfile(f), n = 1L), "\t", fixed = TRUE)[[1L]][-1L]
})
names(dge_barcodes) <- manifest$library_id
all_deposited_barcodes <- unique(unlist(dge_barcodes, use.names = FALSE))

source_match_summary <- do.call(rbind, lapply(
  setdiff(unique(source_cells$Name), "Plasma"),
  function(nm) {
    z <- unique(stats::na.omit(source_cells$base_barcode[source_cells$Name == nm]))
    overlap <- vapply(dge_barcodes, function(b) sum(z %in% b), integer(1))
    best <- which.max(overlap)
    assigned_idx <- which(manifest$patient_id == nm)
    data.frame(
      patient_id = nm,
      source_rows = sum(source_cells$Name == nm),
      source_simple_unique_barcodes = length(z),
      simple_barcodes_present_anywhere = sum(z %in% all_deposited_barcodes),
      best_single_library_overlap = overlap[[best]],
      best_library = names(overlap)[[best]],
      best_library_assigned_patient = manifest$patient_id[[best]],
      overlap_with_assigned_patient_libraries = if (length(assigned_idx)) {
        sum(z %in% unique(unlist(dge_barcodes[assigned_idx], use.names = FALSE)))
      } else NA_integer_,
      barcode_interpretation_note = if (any(is.na(source_cells$base_barcode[source_cells$Name == nm]))) {
        "source identifiers include pooled-library names rather than recoverable 12-bp barcodes"
      } else {
        "12-bp source barcodes can be compared to deposited DGE headers"
      },
      stringsAsFactors = FALSE
    )
  }
))
write_tsv(source_match_summary, file.path(audit, "published_source_barcode_match.tsv"))

uc_source <- source_cells[source_cells$Name == "UC01", , drop = FALSE]
current_uc_libs <- which(manifest$patient_id == "UC01")
current_uc_barcodes <- unique(unlist(dge_barcodes[current_uc_libs], use.names = FALSE))
uc_barcode_audit <- data.frame(
  published_barcode = uc_source$Barcode,
  base_barcode = uc_source$base_barcode,
  published_TM4SF1_log_normalized = uc_source$TM4SF1,
  found_in_current_12049_DGE = uc_source$base_barcode %in% current_uc_barcodes,
  found_in_any_deposited_DGE = uc_source$base_barcode %in% all_deposited_barcodes,
  matching_libraries = vapply(uc_source$base_barcode, function(b) {
    paste(names(dge_barcodes)[vapply(dge_barcodes, function(z) b %in% z, logical(1))], collapse = ",")
  }, character(1)),
  stringsAsFactors = FALSE
)
write_tsv(uc_barcode_audit, file.path(audit, "uc01_published_barcode_audit.tsv"))

# UC01 cluster and marker audit using current deposited 12049 cells.
uc_cells <- rownames(obj@meta.data)[obj$patient_id == "UC01"]
uc_md <- obj@meta.data[uc_cells, , drop = FALSE]
markers <- intersect(c(
  "EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "KRT20", "KRT5", "KRT14",
  "UPK1A", "UPK1B", "UPK2", "UPK3A", "PTPRC", "DCN", "ACTA2", "SELE",
  "PECAM1", "VWF"
), rownames(obj))
uc_raw <- counts[markers, uc_cells, drop = FALSE]
uc_norm <- norm[markers, uc_cells, drop = FALSE]

cluster_ids <- sort(unique(as.character(uc_md$seurat_clusters)))
cluster_annotation <- do.call(rbind, lapply(cluster_ids, function(cl) {
  cells <- rownames(uc_md)[as.character(uc_md$seurat_clusters) == cl]
  z <- uc_md[cells, , drop = FALSE]
  data.frame(
    seurat_cluster = cl,
    n_cells = nrow(z),
    current_broad_cell_type = names(sort(table(z$broad_cell_type), decreasing = TRUE))[1L],
    annotation_confidence = names(sort(table(z$annotation_confidence), decreasing = TRUE))[1L],
    median_score_epithelial = median(z$score_epithelial),
    median_score_immune = median(z$score_immune),
    median_score_stromal = median(z$score_stromal),
    median_score_endothelial = median(z$score_endothelial),
    median_nCount_RNA = median(z$nCount_RNA),
    median_nFeature_RNA = median(z$nFeature_RNA),
    median_percent_mt = median(z$percent_mt_raw),
    n_current_tumour_epithelial = sum(z$tumour_epithelial_identity),
    stringsAsFactors = FALSE
  )
}))
write_tsv(cluster_annotation, file.path(audit, "uc01_cluster_annotation.tsv"))

marker_summary <- do.call(rbind, lapply(cluster_ids, function(cl) {
  cells <- rownames(uc_md)[as.character(uc_md$seurat_clusters) == cl]
  data.frame(
    seurat_cluster = cl,
    n_cells = length(cells),
    current_broad_cell_type = unique(uc_md[cells, "broad_cell_type"])[1L],
    gene = markers,
    pct_detected = as.numeric(Matrix::rowMeans(uc_raw[, cells, drop = FALSE] > 0)),
    mean_raw_UMI = as.numeric(Matrix::rowMeans(uc_raw[, cells, drop = FALSE])),
    mean_log_normalized = as.numeric(Matrix::rowMeans(uc_norm[, cells, drop = FALSE])),
    stringsAsFactors = FALSE
  )
}))
write_tsv(marker_summary, file.path(audit, "uc01_marker_summary.tsv"))

cell_audit <- uc_md[, c(
  "GSM", "library_id", "patient_id", "sample_type", "seurat_clusters",
  "broad_cell_type", "annotation_confidence", "score_epithelial", "score_immune",
  "score_stromal", "score_endothelial", "nCount_RNA", "nFeature_RNA",
  "percent_mt_raw", "doublet_score", "doublet_class",
  "tumour_epithelial_identity"
), drop = FALSE]
cell_audit$cell_id <- rownames(cell_audit)
for (g in markers) cell_audit[[paste0(g, "_raw")]] <- as.numeric(uc_raw[g, rownames(cell_audit)])
write_gz_tsv(cell_audit, file.path(audit, "uc01_qc_pass_cell_audit.tsv.gz"))

uc_sub <- subset(obj, cells = uc_cells)
Idents(uc_sub) <- "seurat_clusters"
eligible_clusters <- names(table(Idents(uc_sub)))[table(Idents(uc_sub)) >= 5L]
uc_marker_obj <- subset(uc_sub, idents = eligible_clusters)
uc_cluster_markers <- FindAllMarkers(
  uc_marker_obj, assay = "RNA", only.pos = TRUE, min.pct = 0.10,
  logfc.threshold = 0.25, test.use = "wilcox", verbose = FALSE
)
write_tsv(uc_cluster_markers, file.path(audit, "uc01_cluster_markers.tsv"))

umap <- Embeddings(obj, "umap")[uc_cells, , drop = FALSE]
plot_df <- cbind(uc_md, UMAP_1 = umap[, 1L], UMAP_2 = umap[, 2L])
p_type <- ggplot(plot_df, aes(UMAP_1, UMAP_2, colour = broad_cell_type)) +
  geom_point(size = 0.65, alpha = 0.85) +
  coord_equal() + theme_bw(base_size = 10) +
  labs(title = "UC01/current 12049: broad annotation", colour = "Broad type")
p_cluster <- ggplot(plot_df, aes(UMAP_1, UMAP_2, colour = factor(seurat_clusters))) +
  geom_point(size = 0.65, alpha = 0.85) +
  coord_equal() + theme_bw(base_size = 10) +
  labs(title = "UC01/current 12049: global Seurat cluster", colour = "Cluster")
ggsave(
  file.path(audit, "uc01_umap_annotations.pdf"),
  p_type + p_cluster, width = 12, height = 5.8
)

feature_plot <- FeaturePlot(
  uc_sub, features = markers, reduction = "umap", ncol = 3,
  order = TRUE, keep.scale = "all", raster = TRUE
)
ggsave(
  file.path(audit, "uc01_marker_featureplots.pdf"),
  feature_plot, width = 13, height = 18, limitsize = FALSE
)

uc_annotation_summary <- as.data.frame(addmargins(table(
  sample_type = uc_md$sample_type,
  broad_cell_type = uc_md$broad_cell_type
)))
names(uc_annotation_summary)[3L] <- "n_cells"
write_tsv(uc_annotation_summary, file.path(audit, "uc01_tumour_normal_annotation_summary.tsv"))

infercnv_audit <- data.frame(
  item = c(
    "paper_settings", "public_code_preparation", "public_runtime_output",
    "current_reconstruction", "UC01_interpretability"
  ),
  finding = c(
    "InferCNV 1.4.0; non-tumour references; genes in >=5 cells; cutoff 0.1; HMM",
    "Private merged object and Final_ID labels; code downsamples 500 per active identity",
    "No InferCNV result, gene-order file, HMM calls, or epithelial keep/drop list released",
    "Input prepared only; exact old InferCNV stack and hg19 gene-order input unavailable",
    "Published 12923/UC01 cells are absent from deposited DGE matrices, so their CNV cannot be rerun"
  ),
  reproducible = c(TRUE, FALSE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)
write_tsv(infercnv_audit, file.path(audit, "infercnv_audit.tsv"))

comparison <- cell_flow[, c(
  "patient_id", "histology_group", "n_processed", "n_broad_epithelial",
  "n_final_tumour_epithelial", "published_tumour_epithelial_cells",
  "delta_tumour_epithelial", "reconstructed_ge150"
)]
write_tsv(comparison, file.path(audit, "published_vs_reconstructed_counts.tsv"))

case_log(
  "Audit tables/UC01 plots complete: source UC01 overlap with current 12049 DGE = ",
  sum(uc_barcode_audit$found_in_current_12049_DGE), "/", nrow(uc_barcode_audit),
  "; current broad epithelial = ", sum(uc_md$broad_cell_type == "epithelial"),
  file = logf
)

