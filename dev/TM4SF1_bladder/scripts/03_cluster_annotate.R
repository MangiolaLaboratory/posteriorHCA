#!/usr/bin/env Rscript
# Reproduce the paper-level Seurat workflow and assign broad cell classes from
# multi-gene module scores plus cluster-level marker detection evidence.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
})

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
source(file.path(.script_dir, "utils.R"))

paths <- case_paths()
logf <- start_log(paths, "03_cluster_annotate")
set.seed(293189)

qc <- readRDS(file.path(paths$processed, "GSE293189_qc_counts.rds"))
keep <- qc$cell_metadata$histology_group %in% c("HV", "UC")
counts <- qc$counts[, keep, drop = FALSE]
md <- qc$cell_metadata[keep, , drop = FALSE]
if (ncol(counts) != nrow(md)) stop("Processed cohort count/metadata mismatch")
case_log(
  "Paper-level cohort after QC: ", ncol(counts),
  " cells (paraganglioma subjects retained only in raw/QC objects)", file = logf
)

obj <- CreateSeuratObject(
  counts = counts,
  project = "GSE293189_bladder",
  min.cells = 0,
  min.features = 0,
  meta.data = md
)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
if (length(VariableFeatures(obj)) != 2000L) stop("Did not obtain 2,000 HVGs")
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 100, verbose = FALSE, seed.use = 293189)
obj <- FindNeighbors(obj, dims = 1:75, k.param = 30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, random.seed = 293189, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:75, seed.use = 293189, verbose = FALSE)

panels <- list(
  epithelial = c(
    "EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "KRT20", "KRT4", "KRT13",
    "KRT5", "KRT14", "TP63", "UPK1A", "UPK1B", "UPK2", "UPK3A"
  ),
  immune = c(
    "PTPRC", "CD3D", "CD3E", "TRBC1", "CD79A", "MS4A1", "CD68", "LYZ",
    "LST1", "FCER1G", "NKG7", "TYROBP"
  ),
  stromal = c(
    "DCN", "COL1A1", "COL1A2", "COL3A1", "LUM", "C7", "ACTA2", "TAGLN",
    "MYL9", "MYH11", "DES", "RGS5"
  ),
  endothelial = c(
    "SELE", "PECAM1", "VWF", "KDR", "EMCN", "ENG", "RAMP2", "PLVAP",
    "CLDN5", "ESAM", "EPAS1", "RGCC"
  )
)
panels <- lapply(panels, intersect, y = rownames(obj))
if (any(lengths(panels) < 4L)) stop("Insufficient marker genes for a broad panel")

for (nm in names(panels)) {
  set.seed(293189)
  raw_name <- paste0("module_", nm)
  obj <- AddModuleScore(
    obj, features = list(panels[[nm]]), assay = "RNA", ctrl = 5,
    name = raw_name, seed = 293189
  )
  names(obj@meta.data)[names(obj@meta.data) == paste0(raw_name, "1")] <- paste0("score_", nm)
}

score_cols <- paste0("score_", names(panels))
marker_counts <- LayerData(
  obj, assay = "RNA", layer = "counts",
  features = unique(unlist(panels, use.names = FALSE))
)
clusters <- sort(unique(as.character(obj$seurat_clusters)))
cluster_rows <- lapply(clusters, function(cl) {
  cells <- rownames(obj@meta.data)[as.character(obj$seurat_clusters) == cl]
  means <- colMeans(obj@meta.data[cells, score_cols, drop = FALSE])
  names(means) <- sub("^score_", "", names(means))
  ord <- order(means, decreasing = TRUE)
  top <- names(means)[ord[[1L]]]
  second <- names(means)[ord[[2L]]]
  pct <- Matrix::rowMeans(marker_counts[panels[[top]], cells, drop = FALSE] > 0)
  evidence <- names(sort(pct, decreasing = TRUE))[pct[order(pct, decreasing = TRUE)] >= 0.05]
  if (!length(evidence)) evidence <- names(sort(pct, decreasing = TRUE))[seq_len(min(3L, length(pct)))]
  n_evidence <- sum(pct >= 0.10)
  margin <- unname(means[[top]] - means[[second]])
  confidence <- if (margin >= 0.15 && n_evidence >= 2L) {
    "high"
  } else if (margin >= 0.05 && n_evidence >= 2L) {
    "moderate"
  } else {
    "low"
  }
  data.frame(
    seurat_cluster = cl,
    n_cells = length(cells),
    broad_cell_type = top,
    annotation_confidence = confidence,
    top_vs_second_score_margin = margin,
    annotation_evidence = paste(head(evidence, 8L), collapse = ","),
    mean_score_epithelial = means[["epithelial"]],
    mean_score_immune = means[["immune"]],
    mean_score_stromal = means[["stromal"]],
    mean_score_endothelial = means[["endothelial"]],
    stringsAsFactors = FALSE
  )
})
cluster_annotation <- do.call(rbind, cluster_rows)
write_tsv(cluster_annotation, file.path(paths$metadata, "broad_cluster_annotation.tsv"))

idx <- match(as.character(obj$seurat_clusters), cluster_annotation$seurat_cluster)
obj$broad_cell_type <- cluster_annotation$broad_cell_type[idx]
obj$annotation_confidence <- cluster_annotation$annotation_confidence[idx]
obj$annotation_evidence <- cluster_annotation$annotation_evidence[idx]
obj$annotation_method <-
  "cluster_top_multi_gene_AddModuleScore_with_marker_detection_evidence"
obj$tumour_epithelial_identity <-
  obj$sample_type == "tumour" & obj$broad_cell_type == "epithelial"
obj$tumour_identity_support <- ifelse(
  obj$tumour_epithelial_identity,
  paste0("tumour_source;epithelial_module;", obj$annotation_evidence),
  NA_character_
)
obj$infercnv_status <-
  "input_prepared_not_run_infercnv_package_and_hg19_gene_order_reference_unavailable"

tumour <- obj@meta.data[
  obj$tumour_epithelial_identity & obj$histology_group %in% c("HV", "UC") & !is.na(obj$paper_id),
  , drop = FALSE
]
tumour_n <- table(tumour$patient_id)
included_patients <- names(tumour_n)[tumour_n >= 150L]
obj$tumour_analysis_included <- obj$tumour_epithelial_identity & obj$patient_id %in% included_patients
obj$tumour_analysis_exclusion_reason <- ifelse(
  !obj$tumour_epithelial_identity,
  "not_tumour_epithelial",
  ifelse(obj$patient_id %in% included_patients, NA_character_, "fewer_than_150_tumour_epithelial_cells")
)

expected_retained <- c(paste0("VAR0", 1:9), paste0("UC0", 1:3))
patient_tumour <- data.frame(
  patient_id = names(tumour_n),
  n_tumour_epithelial = as.integer(tumour_n),
  reconstructed_retained = names(tumour_n) %in% included_patients,
  published_expected_retained = names(tumour_n) %in% expected_retained,
  stringsAsFactors = FALSE
)
patient_info <- unique(obj@meta.data[, c(
  "patient_id", "original_patient_id", "paper_id", "histology", "histology_group"
)])
patient_tumour <- merge(patient_info, patient_tumour, by = "patient_id", all.x = TRUE, sort = FALSE)
patient_tumour$n_tumour_epithelial[is.na(patient_tumour$n_tumour_epithelial)] <- 0L
patient_tumour$reconstructed_retained[is.na(patient_tumour$reconstructed_retained)] <- FALSE
patient_tumour$published_expected_retained[is.na(patient_tumour$published_expected_retained)] <-
  patient_tumour$patient_id[is.na(patient_tumour$published_expected_retained)] %in% expected_retained
write_tsv(patient_tumour, file.path(paths$metadata, "tumour_epithelial_counts_by_patient.tsv"))

# Prepare the public-code-compatible InferCNV count/annotation inputs. The HMM
# run is deliberately not claimed without infercnv and the study's hg19 gene-order reference.
infer_cells <- rownames(obj@meta.data)[
  obj$tumour_epithelial_identity | obj$broad_cell_type %in% c("immune", "stromal", "endothelial")
]
infer_counts <- LayerData(obj, assay = "RNA", layer = "counts")[, infer_cells, drop = FALSE]
infer_counts <- infer_counts[Matrix::rowSums(infer_counts > 0) >= 5L, , drop = FALSE]
infer_group <- ifelse(
  obj@meta.data[infer_cells, "tumour_epithelial_identity"],
  obj@meta.data[infer_cells, "patient_id"],
  obj@meta.data[infer_cells, "broad_cell_type"]
)
saveRDS(
  list(
    counts = infer_counts,
    annotation = data.frame(cell_id = infer_cells, group = infer_group),
    reference_groups = c("immune", "stromal", "endothelial"),
    author_parameters = list(min_cells_per_gene = 5L, cutoff = 0.1, HMM = TRUE),
    status = "prepared_not_run"
  ),
  file.path(paths$processed, "GSE293189_infercnv_input.rds"), compress = "gzip"
)

saveRDS(obj, file.path(paths$processed, "GSE293189_processed_seurat.rds"), compress = "gzip")

final_md <- obj@meta.data
umap <- Embeddings(obj, "umap")
final_md$UMAP_1 <- umap[rownames(final_md), 1]
final_md$UMAP_2 <- umap[rownames(final_md), 2]
final_md$cell_id <- rownames(final_md)
final_md[] <- lapply(final_md, function(x) if (is.factor(x)) as.character(x) else x)

all_md <- readRDS(file.path(paths$processed, "GSE293189_all_cell_qc_metadata.rds"))
add_cols <- setdiff(names(final_md), names(all_md))
all_md[add_cols] <- NA
matched <- match(rownames(final_md), rownames(all_md))
all_md[matched, names(final_md)] <- final_md[, names(final_md), drop = FALSE]
all_md$analysis_inclusion <- ifelse(
  !all_md$pass_doublet_filter,
  "failed_qc_or_doublet",
  ifelse(all_md$histology_group == "other", "excluded_non_paper_paraganglioma", "processed_paper_cohort")
)
all_md$cell_id <- rownames(all_md)
con <- gzfile(file.path(paths$metadata, "GSE293189_cell_metadata.tsv.gz"), "wt")
utils::write.table(all_md, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
close(con)

tumour_md <- final_md[final_md$tumour_epithelial_identity, , drop = FALSE]
con <- gzfile(file.path(paths$metadata, "tumour_epithelial_cell_metadata.tsv.gz"), "wt")
utils::write.table(tumour_md, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
close(con)

case_log(
  "Clustering/annotation complete: ", length(clusters), " clusters; ",
  sum(obj$broad_cell_type == "epithelial"), " epithelial cells; ",
  sum(obj$tumour_analysis_included), " tumour epithelial cells in >=150-cell patients; retained ",
  paste(sort(included_patients), collapse = ","), file = logf
)
