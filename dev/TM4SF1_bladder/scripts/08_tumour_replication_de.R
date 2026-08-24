#!/usr/bin/env Rscript
# Closest defensible tumour reclustering and HV-vs-UC replication analyses.
# The public DGE reconstruction cannot supply 150 UC01 epithelial cells, so the
# requested 12-patient/150-cell analysis is not mislabeled as reproduced.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
  library(presto)
  library(ggplot2)
  library(ggrepel)
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
logf <- start_log(paths, "08_tumour_replication_de")
fixed_seed <- 293189L
set.seed(fixed_seed)

obj <- readRDS(file.path(paths$processed, "GSE293189_processed_seurat.rds"))
expected_patients <- c(paste0("VAR0", 1:9), paste0("UC0", 1:3))
strict_patients <- c(paste0("VAR0", 1:9), "UC02", "UC03")
tumour_cells <- rownames(obj@meta.data)[
  obj$tumour_epithelial_identity & obj$patient_id %in% expected_patients
]
tumour_counts <- LayerData(obj, assay = "RNA", layer = "counts")[, tumour_cells, drop = FALSE]
tumour_md <- obj@meta.data[tumour_cells, , drop = FALSE]

# Recluster current 12-ID public-DGE reconstruction with the paper Methods
# parameters. This remains an audit object because current UC01 has only 90 cells.
tumour_obj <- CreateSeuratObject(
  tumour_counts, project = "GSE293189_tumour_replication_audit",
  min.cells = 0, min.features = 0, meta.data = tumour_md
)
tumour_obj <- NormalizeData(tumour_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
tumour_obj <- FindVariableFeatures(tumour_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
tumour_obj <- ScaleData(tumour_obj, features = VariableFeatures(tumour_obj), verbose = FALSE)
tumour_obj <- RunPCA(tumour_obj, features = VariableFeatures(tumour_obj), npcs = 100, seed.use = fixed_seed, verbose = FALSE)
tumour_obj <- FindNeighbors(tumour_obj, dims = 1:75, k.param = 30, verbose = FALSE)
tumour_obj <- FindClusters(tumour_obj, resolution = 0.5, random.seed = fixed_seed, verbose = FALSE)
tumour_obj <- RunUMAP(tumour_obj, dims = 1:75, seed.use = fixed_seed, verbose = FALSE)

signature_genes <- intersect(c("MUC16", "MUC4", "KRT24", "WISP2"), rownames(tumour_obj))
tumour_obj <- AddModuleScore(
  tumour_obj, features = list(signature_genes), ctrl = 5,
  name = "paper_cluster13_signature", seed = fixed_seed
)
clusters <- sort(unique(as.character(tumour_obj$seurat_clusters)))
tumour_norm <- LayerData(tumour_obj, assay = "RNA", layer = "data")
cluster_summary <- do.call(rbind, lapply(clusters, function(cl) {
  cells <- rownames(tumour_obj@meta.data)[as.character(tumour_obj$seurat_clusters) == cl]
  z <- tumour_obj@meta.data[cells, , drop = FALSE]
  means <- Matrix::rowMeans(tumour_norm[signature_genes, cells, drop = FALSE])
  data.frame(
    tumour_cluster = cl,
    n_cells = length(cells),
    n_patients = length(unique(z$patient_id)),
    n_HV_patients = length(unique(z$patient_id[z$histology_group == "HV"])),
    n_UC_patients = length(unique(z$patient_id[z$histology_group == "UC"])),
    mean_cluster13_module_score = mean(z$paper_cluster13_signature1),
    mean_MUC16 = unname(means["MUC16"]),
    mean_MUC4 = unname(means["MUC4"]),
    mean_KRT24 = unname(means["KRT24"]),
    mean_WISP2 = unname(means["WISP2"]),
    top_patient = names(sort(table(z$patient_id), decreasing = TRUE))[1L],
    top_patient_fraction = max(table(z$patient_id)) / length(cells),
    stringsAsFactors = FALSE
  )
}))
candidate_pool <- cluster_summary$n_cells >= 50L & cluster_summary$n_patients >= 3L
candidate_cluster <- cluster_summary$tumour_cluster[
  which.max(replace(cluster_summary$mean_cluster13_module_score, !candidate_pool, -Inf))
]
cluster_summary$audit_cluster13_candidate <- cluster_summary$tumour_cluster == candidate_cluster
write_tsv(cluster_summary, file.path(audit, "tumour_recluster_summary.tsv"))

cluster_composition <- as.data.frame(table(
  tumour_cluster = as.character(tumour_obj$seurat_clusters),
  patient_id = tumour_obj$patient_id
), stringsAsFactors = FALSE)
cluster_composition <- cluster_composition[cluster_composition$Freq > 0, , drop = FALSE]
names(cluster_composition)[3L] <- "n_cells"
cluster_composition$within_cluster_fraction <- ave(
  cluster_composition$n_cells, cluster_composition$tumour_cluster,
  FUN = function(x) x / sum(x)
)
write_tsv(
  cluster_composition,
  file.path(audit, "tumour_recluster_composition.tsv")
)

tumour_obj$audit_cluster13_candidate <-
  as.character(tumour_obj$seurat_clusters) == candidate_cluster
Idents(tumour_obj) <- "seurat_clusters"
candidate_markers <- FindMarkers(
  tumour_obj, ident.1 = candidate_cluster, only.pos = TRUE,
  min.pct = 0.10, logfc.threshold = 0.25, test.use = "wilcox",
  verbose = FALSE
)
candidate_markers$gene <- rownames(candidate_markers)
write_tsv(
  candidate_markers,
  file.path(audit, "tumour_cluster13_candidate_markers.tsv")
)
saveRDS(tumour_obj, file.path(audit, "tumour_epithelial_reclustered_audit.rds"), compress = "gzip")

signature_summary <- do.call(rbind, lapply(
  split(rownames(tumour_obj@meta.data), tumour_obj$patient_id),
  function(cells) {
    data.frame(
      patient_id = tumour_obj$patient_id[cells][1L],
      histology_group = tumour_obj$histology_group[cells][1L],
      gene = signature_genes,
      n_cells = length(cells),
      pct_detected = as.numeric(Matrix::rowMeans(
        LayerData(tumour_obj, assay = "RNA", layer = "counts")[signature_genes, cells, drop = FALSE] > 0
      )),
      mean_log_normalized = as.numeric(Matrix::rowMeans(tumour_norm[signature_genes, cells, drop = FALSE])),
      stringsAsFactors = FALSE
    )
  }
))
write_tsv(signature_summary, file.path(audit, "MUC16_MUC4_KRT24_WISP2_expression.tsv"))

# Requested UMAPs and signature views.
p_patient <- DimPlot(tumour_obj, reduction = "umap", group.by = "patient_id", raster = TRUE) +
  ggtitle("Tumour epithelial audit UMAP: patient") + theme_bw(base_size = 10)
p_histology <- DimPlot(tumour_obj, reduction = "umap", group.by = "histology_group", raster = TRUE) +
  ggtitle("Tumour epithelial audit UMAP: HV versus UC") + theme_bw(base_size = 10)
p_clusters <- DimPlot(tumour_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = TRUE) +
  ggtitle("Tumour epithelial audit reclustering") + theme_bw(base_size = 10)
p_candidate <- DimPlot(tumour_obj, reduction = "umap", group.by = "audit_cluster13_candidate", raster = TRUE) +
  ggtitle(paste0("Closest multi-patient signature cluster: ", candidate_cluster)) + theme_bw(base_size = 10)
ggsave(file.path(audit, "tumour_epithelial_umap_patient.pdf"), p_patient, width = 9, height = 7)
ggsave(file.path(audit, "tumour_epithelial_umap_histology.pdf"), p_histology, width = 8, height = 7)
ggsave(
  file.path(audit, "tumour_reclustering_cluster13_candidate.pdf"),
  p_clusters + p_candidate, width = 14, height = 6.5
)
signature_plot <- FeaturePlot(
  tumour_obj, features = c(signature_genes, "TM4SF1"), reduction = "umap",
  ncol = 3, keep.scale = "all", order = TRUE, raster = TRUE
)
ggsave(
  file.path(audit, "tumour_signature_expression.pdf"), signature_plot,
  width = 13, height = 9, limitsize = FALSE
)

# Use the already LogNormalized full-cohort layer for all direct DEG analyses.
full_norm <- LayerData(obj, assay = "RNA", layer = "data")
feature_names <- rownames(full_norm)

sample_by_patient <- function(cells, n_per_patient, seed) {
  set.seed(seed)
  groups <- split(cells, obj@meta.data[cells, "patient_id"])
  if (any(lengths(groups) < n_per_patient)) {
    stop("At least one patient has fewer than ", n_per_patient, " tumour cells")
  }
  unlist(lapply(groups, sample, size = n_per_patient, replace = FALSE), use.names = FALSE)
}

run_de <- function(cells) {
  x <- full_norm[, cells, drop = FALSE]
  group <- as.character(obj@meta.data[cells, "histology_group"])
  w <- presto::wilcoxauc(x, group)
  w <- w[w$group == "HV", , drop = FALSE]
  w <- w[match(feature_names, w$feature), , drop = FALSE]
  hv <- group == "HV"
  mean_hv <- log2((Matrix::rowSums(expm1(x[, hv, drop = FALSE])) + 1) / sum(hv))
  mean_uc <- log2((Matrix::rowSums(expm1(x[, !hv, drop = FALSE])) + 1) / sum(!hv))
  out <- data.frame(
    gene = feature_names,
    p_value = w$pval,
    FDR_BH = p.adjust(w$pval, method = "BH"),
    p_value_adj_Seurat_Bonferroni = p.adjust(w$pval, method = "bonferroni", n = length(feature_names)),
    avg_log2FC = as.numeric(mean_hv - mean_uc),
    pct_HV = w$pct_in / 100,
    pct_UC = w$pct_out / 100,
    AUC_HV = w$auc,
    stringsAsFactors = FALSE
  )
  positive <- which(out$avg_log2FC > 0 & !is.na(out$p_value))
  ord_p <- positive[order(out$p_value[positive], -out$avg_log2FC[positive], out$gene[positive])]
  ord_fc <- positive[order(-out$avg_log2FC[positive], out$p_value[positive], out$gene[positive])]
  out$HV_enriched_rank_by_p <- NA_integer_
  out$HV_enriched_rank_by_logFC <- NA_integer_
  out$HV_enriched_rank_by_p[ord_p] <- seq_along(ord_p)
  out$HV_enriched_rank_by_logFC[ord_fc] <- seq_along(ord_fc)
  out
}

strict_cells <- tumour_cells[obj$patient_id[tumour_cells] %in% strict_patients]
fixed_cells_150 <- sample_by_patient(strict_cells, 150L, fixed_seed)
fixed_de <- run_de(fixed_cells_150)
fixed_de$analysis <- "closest_public_DGE_9HV_2UC_150_cells_per_patient_seed293189"
fixed_de$n_HV_patients <- 9L
fixed_de$n_UC_patients <- 2L
fixed_de$n_cells_per_patient <- 150L
write_tsv(fixed_de, file.path(audit, "hv_uc_150cells_DE.tsv"))

top30 <- fixed_de[
  !is.na(fixed_de$HV_enriched_rank_by_p) & fixed_de$HV_enriched_rank_by_p <= 30L,
  , drop = FALSE
]
top30 <- top30[order(top30$HV_enriched_rank_by_p), ]
write_tsv(top30, file.path(audit, "top30_HV_enriched_genes.tsv"))

# A 12-patient analysis is possible only at 90 cells per patient in the current
# reconstruction; it is explicitly secondary and not the requested replication.
fixed_cells_90 <- sample_by_patient(tumour_cells, 90L, fixed_seed)
fixed_de_90 <- run_de(fixed_cells_90)
write_tsv(fixed_de_90, file.path(audit, "hv_uc_all12_90cells_DE.tsv"))

# Extract the paper's released fixed-run DEG table for an external truth check.
published_de <- as.data.frame(read_excel(
  file.path(paths$raw_paper, "41467_2025_59888_MOESM4_ESM.xlsx"),
  sheet = "Figure5A"
))
published_de$p_val <- as.numeric(published_de$p_val)
published_de$p_val_adj <- as.numeric(published_de$p_val_adj)
published_variant <- published_de[published_de$cluster == "Variant", , drop = FALSE]
positive_pub <- which(published_variant$avg_log2FC > 0 & !is.na(published_variant$p_val))
ord_pub <- positive_pub[order(
  published_variant$p_val[positive_pub],
  -published_variant$avg_log2FC[positive_pub],
  published_variant$gene[positive_pub]
)]
published_variant$HV_enriched_rank_by_p <- NA_integer_
published_variant$HV_enriched_rank_by_p[ord_pub] <- seq_along(ord_pub)
write_tsv(published_variant, file.path(audit, "published_source_hv_uc_DE.tsv"))
published_top30 <- published_variant[
  !is.na(published_variant$HV_enriched_rank_by_p) &
    published_variant$HV_enriched_rank_by_p <= 30L,
  , drop = FALSE
]
published_top30 <- published_top30[order(published_top30$HV_enriched_rank_by_p), ]
write_tsv(
  published_top30,
  file.path(audit, "published_source_top30_HV_enriched_genes.tsv")
)

tm_rows <- rbind(
  transform(
    published_variant[published_variant$gene == "TM4SF1", c(
      "gene", "p_val", "p_val_adj", "avg_log2FC", "pct.1", "pct.2",
      "HV_enriched_rank_by_p"
    )],
    analysis = "paper_released_source_DE_9HV_3UC_150_cells_per_patient_unknown_seed",
    FDR_BH = NA_real_, n_HV_patients = 9L, n_UC_patients = 3L,
    n_cells_per_patient = 150L
  ),
  data.frame(
    gene = "TM4SF1",
    p_val = fixed_de$p_value[fixed_de$gene == "TM4SF1"],
    p_val_adj = fixed_de$p_value_adj_Seurat_Bonferroni[fixed_de$gene == "TM4SF1"],
    avg_log2FC = fixed_de$avg_log2FC[fixed_de$gene == "TM4SF1"],
    pct.1 = fixed_de$pct_HV[fixed_de$gene == "TM4SF1"],
    pct.2 = fixed_de$pct_UC[fixed_de$gene == "TM4SF1"],
    HV_enriched_rank_by_p = fixed_de$HV_enriched_rank_by_p[fixed_de$gene == "TM4SF1"],
    analysis = "closest_public_DGE_9HV_2UC_150_cells_per_patient_seed293189",
    FDR_BH = fixed_de$FDR_BH[fixed_de$gene == "TM4SF1"],
    n_HV_patients = 9L, n_UC_patients = 2L, n_cells_per_patient = 150L
  ),
  data.frame(
    gene = "TM4SF1",
    p_val = fixed_de_90$p_value[fixed_de_90$gene == "TM4SF1"],
    p_val_adj = fixed_de_90$p_value_adj_Seurat_Bonferroni[fixed_de_90$gene == "TM4SF1"],
    avg_log2FC = fixed_de_90$avg_log2FC[fixed_de_90$gene == "TM4SF1"],
    pct.1 = fixed_de_90$pct_HV[fixed_de_90$gene == "TM4SF1"],
    pct.2 = fixed_de_90$pct_UC[fixed_de_90$gene == "TM4SF1"],
    HV_enriched_rank_by_p = fixed_de_90$HV_enriched_rank_by_p[fixed_de_90$gene == "TM4SF1"],
    analysis = "secondary_public_DGE_9HV_3UC_90_cells_per_patient_seed293189",
    FDR_BH = fixed_de_90$FDR_BH[fixed_de_90$gene == "TM4SF1"],
    n_HV_patients = 9L, n_UC_patients = 3L, n_cells_per_patient = 90L
  )
)
write_tsv(tm_rows, file.path(audit, "TM4SF1_original_DE_result.tsv"))

# Rank sensitivity for the two analyses feasible from deposited count matrices.
tm_seed_one <- function(cells, n_per_patient, seed, analysis) {
  sampled <- sample_by_patient(cells, n_per_patient, seed)
  x <- full_norm[, sampled, drop = FALSE]
  group <- as.character(obj@meta.data[sampled, "histology_group"])
  w <- presto::wilcoxauc(x, group)
  w <- w[w$group == "HV", , drop = FALSE]
  positive <- which(w$logFC > 0 & !is.na(w$pval))
  ord <- positive[order(w$pval[positive], -w$logFC[positive], w$feature[positive])]
  tm <- which(w$feature == "TM4SF1")
  data.frame(
    analysis = analysis,
    seed = seed,
    n_cells_per_patient = n_per_patient,
    n_HV_patients = length(unique(obj$patient_id[sampled][group == "HV"])),
    n_UC_patients = length(unique(obj$patient_id[sampled][group == "UC"])),
    TM4SF1_logFC_mean_log_normalized = w$logFC[tm],
    TM4SF1_p_value = w$pval[tm],
    TM4SF1_FDR_BH = w$padj[tm],
    TM4SF1_rank_among_HV_enriched_by_p = match(tm, ord),
    n_HV_enriched_tested = length(ord),
    stringsAsFactors = FALSE
  )
}

sensitivity <- vector("list", 200L)
for (seed in 1:100) {
  sensitivity[[seed]] <- tm_seed_one(
    strict_cells, 150L, seed,
    "closest_public_DGE_9HV_2UC_150_cells_per_patient"
  )
  sensitivity[[100L + seed]] <- tm_seed_one(
    tumour_cells, 90L, seed,
    "secondary_public_DGE_9HV_3UC_90_cells_per_patient"
  )
  if (seed %% 20L == 0L) case_log("DE sensitivity seed ", seed, "/100", file = logf)
}
sensitivity <- do.call(rbind, sensitivity)
write_tsv(sensitivity, file.path(audit, "TM4SF1_downsampling_sensitivity.tsv"))

sensitivity_summary <- do.call(rbind, lapply(
  split(sensitivity, sensitivity$analysis),
  function(z) data.frame(
    analysis = z$analysis[1L],
    n_runs = nrow(z),
    median_rank = median(z$TM4SF1_rank_among_HV_enriched_by_p, na.rm = TRUE),
    rank_min = min(z$TM4SF1_rank_among_HV_enriched_by_p, na.rm = TRUE),
    rank_Q1 = unname(quantile(z$TM4SF1_rank_among_HV_enriched_by_p, 0.25, na.rm = TRUE)),
    rank_Q3 = unname(quantile(z$TM4SF1_rank_among_HV_enriched_by_p, 0.75, na.rm = TRUE)),
    rank_max = max(z$TM4SF1_rank_among_HV_enriched_by_p, na.rm = TRUE),
    fraction_top1 = mean(z$TM4SF1_rank_among_HV_enriched_by_p <= 1, na.rm = TRUE),
    fraction_top5 = mean(z$TM4SF1_rank_among_HV_enriched_by_p <= 5, na.rm = TRUE),
    fraction_top10 = mean(z$TM4SF1_rank_among_HV_enriched_by_p <= 10, na.rm = TRUE),
    median_logFC = median(z$TM4SF1_logFC_mean_log_normalized, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
))
write_tsv(
  sensitivity_summary,
  file.path(audit, "TM4SF1_downsampling_sensitivity_summary.tsv")
)

# The released workbook permits an exact 12-patient TM4SF1-only resampling,
# though not genome-wide ranks because it releases only TM4SF1 at cell level.
source_tm <- as.data.frame(read_excel(
  file.path(paths$raw_paper, "41467_2025_59888_MOESM4_ESM.xlsx"), sheet = "Figure5B"
))
source_tm <- source_tm[source_tm$Name %in% expected_patients, , drop = FALSE]
source_tm$histology_group <- ifelse(grepl("^VAR", source_tm$Name), "HV", "UC")
source_sensitivity <- do.call(rbind, lapply(1:100, function(seed) {
  set.seed(seed)
  ii <- unlist(lapply(split(seq_len(nrow(source_tm)), source_tm$Name), sample, 150L))
  z <- source_tm[ii, , drop = FALSE]
  wt <- wilcox.test(TM4SF1 ~ histology_group, data = z, alternative = "two.sided", exact = FALSE)
  data.frame(
    seed = seed,
    mean_log_normalized_HV = mean(z$TM4SF1[z$histology_group == "HV"]),
    mean_log_normalized_UC = mean(z$TM4SF1[z$histology_group == "UC"]),
    mean_difference_HV_minus_UC =
      mean(z$TM4SF1[z$histology_group == "HV"]) - mean(z$TM4SF1[z$histology_group == "UC"]),
    p_value = wt$p.value,
    note = "exact released TM4SF1 values; genome-wide rank unavailable",
    stringsAsFactors = FALSE
  )
}))
write_tsv(
  source_sensitivity,
  file.path(audit, "TM4SF1_published_source_expression_sensitivity.tsv")
)

# Paper-facing expression and DEG plots from the feasible fixed analysis.
plot_cells <- fixed_cells_150
plot_df <- data.frame(
  cell_id = plot_cells,
  patient_id = obj$patient_id[plot_cells],
  histology_group = obj$histology_group[plot_cells],
  TM4SF1 = as.numeric(full_norm["TM4SF1", plot_cells]),
  stringsAsFactors = FALSE
)
p_tm_patient <- ggplot(plot_df, aes(patient_id, TM4SF1, fill = histology_group)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, colour = "black") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "TM4SF1 by patient: closest 9 HV + 2 UC / 150-cell analysis", x = NULL)
p_tm_group <- ggplot(plot_df, aes(histology_group, TM4SF1, fill = histology_group)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.13, outlier.shape = NA) +
  theme_bw(base_size = 11) + guides(fill = "none") +
  labs(title = "TM4SF1: HV versus UC", x = NULL)
ggsave(file.path(audit, "TM4SF1_expression_by_patient.pdf"), p_tm_patient, width = 10, height = 6)
ggsave(file.path(audit, "TM4SF1_expression_HV_vs_UC.pdf"), p_tm_group, width = 5.5, height = 5.5)

volcano <- fixed_de
positive_p <- volcano$p_value[is.finite(volcano$p_value) & volcano$p_value > 0]
p_floor <- if (length(positive_p)) min(positive_p) / 10 else .Machine$double.xmin
volcano$minus_log10_p <- -log10(pmax(volcano$p_value, p_floor, na.rm = TRUE))
volcano$highlight <- ifelse(volcano$gene == "TM4SF1", "TM4SF1", "Other")
label_genes <- unique(c("TM4SF1", head(top30$gene, 10L)))
p_volcano <- ggplot(volcano, aes(avg_log2FC, minus_log10_p)) +
  geom_point(aes(colour = highlight), size = 0.65, alpha = 0.55) +
  geom_text_repel(
    data = volcano[volcano$gene %in% label_genes, ], aes(label = gene),
    size = 3, max.overlaps = Inf
  ) +
  scale_colour_manual(values = c(Other = "grey65", TM4SF1 = "#D73027")) +
  theme_bw(base_size = 10) + guides(colour = "none") +
  labs(
    title = "Closest public-DGE HV-vs-UC replication (9 HV + 2 UC)",
    x = "Average log2 fold-change (HV - UC)", y = "-log10(two-sided Wilcoxon p)"
  )
ggsave(file.path(audit, "hv_uc_150cells_volcano.pdf"), p_volcano, width = 8.5, height = 7)

top30$gene <- factor(top30$gene, levels = rev(top30$gene))
p_top <- ggplot(top30, aes(avg_log2FC, gene)) +
  geom_col(fill = "#B2182B") + theme_bw(base_size = 9) +
  labs(
    title = "Top 30 HV-enriched genes by Wilcoxon p-value",
    subtitle = "Closest public-DGE 9 HV + 2 UC / 150 cells per patient",
    x = "Average log2 fold-change", y = NULL
  )
ggsave(file.path(audit, "top30_HV_enriched_genes.pdf"), p_top, width = 8, height = 9)

case_log(
  "Tumour audit/DE complete: candidate shared-signature cluster ", candidate_cluster,
  "; fixed 9HV+2UC TM4SF1 rank = ",
  fixed_de$HV_enriched_rank_by_p[fixed_de$gene == "TM4SF1"],
  file = logf
)
