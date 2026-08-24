# Build SummarizedExperiment objects for Code Ocean GeoMx (all segments).
#
# Outputs under data/processed/:
#   GeoMx_all_segments_SE.rds
#       gene × ROI (93 segments), colData = ROI metadata
#   GeoMx_patient_compartment_pseudobulk_SE.rds
#       gene × (patient × compartment), assays = expression (sum) + counts (rounded)
#
# Run:
#   bash -lc 'use_R_env && Rscript /home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC/data/build_geomx_summarized_experiments.R'

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(S4Vectors)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})

case_dir <- "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC"
data_dir <- file.path(case_dir, "data")
proc_dir <- file.path(data_dir, "processed")

src_rds <- file.path(proc_dir, "TIMES_CodeOcean_GeoMx_TC_processed.rds")
biomarker_csv <- file.path(data_dir, "times_biomarkers.csv")
stopifnot(file.exists(src_rds), file.exists(biomarker_csv))

message("Loading Code Ocean extract: ", src_rds)
co <- readRDS(src_rds)
biomarkers <- read_csv(biomarker_csv, show_col_types = FALSE)

expr <- as.matrix(co$all_expression)
mode(expr) <- "numeric"
meta <- as.data.frame(co$all_metadata, stringsAsFactors = FALSE)

# Align columns to metadata order
stopifnot(all(colnames(expr) %in% meta$author_column))
meta <- meta[match(colnames(expr), meta$author_column), , drop = FALSE]
stopifnot(identical(colnames(expr), meta$author_column))

# ---------------------------------------------------------------------------
# rowData
# ---------------------------------------------------------------------------

genes <- rownames(expr)
row_df <- data.frame(
  gene = genes,
  .feature = biomarkers$.feature[match(genes, biomarkers$gene)],
  is_times_biomarker = genes %in% biomarkers$gene,
  row.names = genes,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ---------------------------------------------------------------------------
# 1) ROI / all-segments SE
# ---------------------------------------------------------------------------

roi_coldata <- DataFrame(
  meta |>
    transmute(
      author_column,
      author_index = as.integer(author_index),
      sampleid_raw,
      label,
      patient_id,
      recurrence,
      relapse,
      death,
      age = as.numeric(age),
      gender,
      compartment,
      cell_context,
      platform,
      included_in_author_nk_analysis = as.logical(included_in_author_nk_analysis)
    ),
  row.names = meta$author_column
)

se_roi <- SummarizedExperiment(
  assays = list(
    expression = expr,
    counts = round(expr)
  ),
  rowData = DataFrame(row_df),
  colData = roi_coldata,
  metadata = list(
    level = "GeoMx ROI / geometric segment",
    source_file = co$source_file,
    platform = "GeoMx_DSP",
    value_scale = paste(
      "Author-processed Code Ocean transcript-read-like values;",
      "assay 'expression' is continuous; assay 'counts' is round(expression).",
      "Not OMIX raw DSP counts."
    ),
    n_genes = nrow(expr),
    n_rois = ncol(expr),
    created_by = "build_geomx_summarized_experiments.R",
    created_at = as.character(Sys.time())
  )
)

out_roi <- file.path(proc_dir, "GeoMx_all_segments_SE.rds")
saveRDS(se_roi, out_roi)
message("Wrote: ", out_roi, " [", paste(dim(se_roi), collapse = " x "), "]")

# ---------------------------------------------------------------------------
# 2) Patient × compartment pseudobulk SE (all segments aggregated)
# ---------------------------------------------------------------------------

pb_keys <- meta |>
  distinct(patient_id, recurrence, age, gender, compartment, platform) |>
  arrange(patient_id, compartment)

pb_ids <- paste(pb_keys$patient_id, pb_keys$compartment, sep = "_")

pb_expr <- matrix(
  NA_real_,
  nrow = nrow(expr),
  ncol = nrow(pb_keys),
  dimnames = list(genes, pb_ids)
)
pb_n_rois <- integer(nrow(pb_keys))
pb_n_nk_rois <- integer(nrow(pb_keys))
pb_n_rois_any_value <- integer(nrow(pb_keys))
names(pb_n_rois) <- pb_ids
names(pb_n_nk_rois) <- pb_ids
names(pb_n_rois_any_value) <- pb_ids

n_non_missing <- matrix(
  0L,
  nrow = nrow(expr),
  ncol = nrow(pb_keys),
  dimnames = list(genes, pb_ids)
)

for (i in seq_len(nrow(pb_keys))) {
  key <- pb_keys[i, ]
  id <- pb_ids[i]
  cols <- meta$author_column[
    meta$patient_id == key$patient_id & meta$compartment == key$compartment
  ]
  mat <- expr[, cols, drop = FALSE]
  nnm <- rowSums(!is.na(mat))
  s <- rowSums(mat, na.rm = TRUE)
  s[nnm == 0] <- NA_real_
  pb_expr[, id] <- s
  n_non_missing[, id] <- as.integer(nnm)
  pb_n_rois[id] <- length(cols)
  pb_n_nk_rois[id] <- sum(
    meta$included_in_author_nk_analysis[
      meta$author_column %in% cols
    ],
    na.rm = TRUE
  )
  pb_n_rois_any_value[id] <- sum(colSums(!is.na(mat)) > 0)
}

pb_coldata <- DataFrame(
  pb_keys |>
    mutate(
      pb_id = pb_ids,
      age = as.numeric(age),
      n_rois = as.integer(pb_n_rois[pb_ids]),
      n_author_nk_rois = as.integer(pb_n_nk_rois[pb_ids]),
      n_rois_with_any_value = as.integer(pb_n_rois_any_value[pb_ids]),
      aggregation = "sum of non-NA ROI expression within patient x compartment",
      cell_context = "all_segments_collapsed"
    ),
  row.names = pb_ids
)

se_pb <- SummarizedExperiment(
  assays = list(
    expression = pb_expr,
    counts = round(pb_expr),
    n_non_missing = n_non_missing
  ),
  rowData = DataFrame(row_df),
  colData = pb_coldata,
  metadata = list(
    level = "patient x compartment pseudobulk (all GeoMx segments)",
    source_file = co$source_file,
    platform = "GeoMx_DSP",
    aggregation = paste(
      "For each gene, sum non-NA ROI values within patient_id x compartment;",
      "all cell_context labels collapsed. Assay n_non_missing = ROI count used."
    ),
    value_scale = paste(
      "Author-processed Code Ocean transcript-read-like values;",
      "assay 'expression' is continuous sum; assay 'counts' is round(expression).",
      "Not OMIX raw DSP counts."
    ),
    n_genes = nrow(pb_expr),
    n_pseudobulk = ncol(pb_expr),
    created_by = "build_geomx_summarized_experiments.R",
    created_at = as.character(Sys.time())
  )
)

out_pb <- file.path(proc_dir, "GeoMx_patient_compartment_pseudobulk_SE.rds")
saveRDS(se_pb, out_pb)
message("Wrote: ", out_pb, " [", paste(dim(se_pb), collapse = " x "), "]")

# ---------------------------------------------------------------------------
# Quick validation summary
# ---------------------------------------------------------------------------

summary_tbl <- tibble(
  object = c("GeoMx_all_segments_SE", "GeoMx_patient_compartment_pseudobulk_SE"),
  path = c(out_roi, out_pb),
  n_genes = c(nrow(se_roi), nrow(se_pb)),
  n_cols = c(ncol(se_roi), ncol(se_pb)),
  assays = c(
    paste(assayNames(se_roi), collapse = ", "),
    paste(assayNames(se_pb), collapse = ", ")
  ),
  n_patients = c(
    length(unique(colData(se_roi)$patient_id)),
    length(unique(colData(se_pb)$patient_id))
  ),
  compartments = c(
    paste(sort(unique(colData(se_roi)$compartment)), collapse = "/"),
    paste(sort(unique(colData(se_pb)$compartment)), collapse = "/")
  )
)

out_summary <- file.path(proc_dir, "GeoMx_SE_build_summary.csv")
write_csv(summary_tbl, out_summary)
message("Wrote: ", out_summary)
print(summary_tbl)
message("colData (pseudobulk):")
print(as.data.frame(colData(se_pb)))
message("Done.")
