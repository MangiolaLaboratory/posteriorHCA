# Extract processed TIMES GeoMx data from the authors' Code Ocean capsule.
#
# Run from the posteriorHCA project root:
#   conda activate R_env
#   Rscript dev/Jie_HCC/data/extract_codeocean_processed_data.R
#
# The Code Ocean capsule stores processed GeoMx NK-enriched ROI measurements in
# RData workspaces. This script extracts the complete GeoMx gene expression
# matrix across AS, IF, and TC ROIs, then writes reusable matrix, segment-level,
# and patient-level files under dev/Jie_HCC/data/processed/. TC/NK and TIMES
# biomarker files are also written as validation/convenience subsets.

message("Extracting Code Ocean processed GeoMx expression data")

script_path <- sub(
  "^--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
)
project_root <- Sys.getenv("POSTERIORHCA_ROOT", unset = NA_character_)
if (is.na(project_root) && length(script_path)) {
  project_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), mustWork = FALSE)
}
if (is.na(project_root) || !dir.exists(project_root)) {
  project_root <- normalizePath(getwd(), mustWork = FALSE)
}

case_dir <- file.path(project_root, "dev", "Jie_HCC")
data_dir <- file.path(case_dir, "data")
processed_dir <- file.path(data_dir, "processed")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

codeocean_input <- file.path(
  case_dir,
  "codeocean",
  "data",
  "Input_1_Spatial_transcriptomics_analysis",
  "transcriptome_differential_analysis_input_1.RData"
)

biomarker_path <- file.path(data_dir, "times_biomarkers.csv")
if (!file.exists(codeocean_input)) {
  stop("Missing Code Ocean spatial transcriptomics input: ", codeocean_input)
}
if (!file.exists(biomarker_path)) {
  stop("Missing biomarker registry: ", biomarker_path)
}

biomarkers <- read.csv(biomarker_path, stringsAsFactors = FALSE, check.names = FALSE)
times_genes <- biomarkers$gene

env <- new.env(parent = emptyenv())
loaded <- load(codeocean_input, envir = env)

required_objects <- c(
  "data_gene_df_ed",
  "data_meta_df_ed",
  "NKatStroma_indx",
  "NKatBorder_indx",
  "NKatTumor_indx",
  "Relp_NKatTumor_indx",
  "Nonrelp_NKatTumor_indx"
)
missing_objects <- setdiff(required_objects, loaded)
if (length(missing_objects)) {
  stop("Code Ocean RData is missing required objects: ", paste(missing_objects, collapse = ", "))
}

gene_matrix <- get("data_gene_df_ed", env)
metadata <- get("data_meta_df_ed", env)

missing_genes <- setdiff(times_genes, rownames(gene_matrix))
if (length(missing_genes)) {
  stop("TIMES genes missing from Code Ocean matrix: ", paste(missing_genes, collapse = ", "))
}

all_columns <- colnames(gene_matrix)
missing_columns <- setdiff(all_columns, colnames(metadata))
if (length(missing_columns)) {
  stop("Expression columns missing from Code Ocean metadata: ", paste(missing_columns, collapse = ", "))
}

classify_compartment <- function(label) {
  ifelse(
    grepl("stroma", label, ignore.case = TRUE), "AS",
    ifelse(
      grepl("border", label, ignore.case = TRUE), "IF",
      ifelse(grepl("tumor|tumour", label, ignore.case = TRUE), "TC", NA_character_)
    )
  )
}

classify_cell_context <- function(label) {
  ifelse(
    grepl("cold tumor", label, ignore.case = TRUE), "cold_tumor",
    ifelse(
      grepl("cd57", label, ignore.case = TRUE) & grepl("cd3", label, ignore.case = TRUE),
      "CD57_rich_CD3_rich_overlap",
      ifelse(
        grepl("cd57", label, ignore.case = TRUE), "CD57_rich",
        ifelse(grepl("cd3", label, ignore.case = TRUE), "CD3_rich", "unspecified")
      )
    )
  )
}

all_meta <- as.data.frame(t(metadata[, all_columns, drop = FALSE]), stringsAsFactors = FALSE)
all_meta$author_column <- rownames(all_meta)
all_meta$author_index <- as.integer(sub("^\\.\\.\\.", "", all_meta$author_column))
rownames(all_meta) <- NULL

all_meta$compartment <- classify_compartment(all_meta$label)
all_meta$platform <- "GeoMx_DSP"
all_meta$cell_context <- classify_cell_context(all_meta$label)
all_meta$recurrence <- ifelse(all_meta$relapse == "Y", "REC", "non-REC")
all_meta$patient_id <- all_meta$sampleid
all_meta$included_in_author_nk_analysis <- all_meta$author_index %in% c(
  get("NKatStroma_indx", env),
  get("NKatBorder_indx", env),
  get("NKatTumor_indx", env)
)

tc_indices <- get("NKatTumor_indx", env)
tc_columns <- paste0("...", tc_indices)
tc_meta <- all_meta[all_meta$author_column %in% tc_columns, , drop = FALSE]

all_genes <- rownames(gene_matrix)
all_expression <- as.matrix(gene_matrix[all_genes, all_columns, drop = FALSE])
mode(all_expression) <- "numeric"
tc_expression <- all_expression[, tc_columns, drop = FALSE]
mode(tc_expression) <- "numeric"

all_expression_df <- data.frame(
  gene = rownames(all_expression),
  all_expression,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

tc_expression_df <- data.frame(
  gene = rownames(tc_expression),
  tc_expression,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

gene_feature <- biomarkers$.feature[match(all_genes, biomarkers$gene)]
is_times_biomarker <- all_genes %in% times_genes

all_segment_long <- data.frame(
  gene = rep(all_genes, times = length(all_columns)),
  .feature = rep(gene_feature, times = length(all_columns)),
  is_times_biomarker = rep(is_times_biomarker, times = length(all_columns)),
  author_column = rep(all_columns, each = length(all_genes)),
  value = as.vector(all_expression),
  stringsAsFactors = FALSE
)

all_segment_long <- merge(all_segment_long, all_meta, by = "author_column", all.x = TRUE, sort = FALSE)
all_segment_long <- all_segment_long[order(all_segment_long$gene, all_segment_long$author_index), ]

tc_segment_long <- data.frame(
  gene = rep(all_genes, times = length(tc_columns)),
  .feature = rep(gene_feature, times = length(tc_columns)),
  is_times_biomarker = rep(is_times_biomarker, times = length(tc_columns)),
  author_column = rep(tc_columns, each = length(all_genes)),
  value = as.vector(tc_expression),
  stringsAsFactors = FALSE
)

tc_segment_long <- merge(tc_segment_long, tc_meta, by = "author_column", all.x = TRUE, sort = FALSE)
tc_segment_long <- tc_segment_long[order(tc_segment_long$gene, tc_segment_long$author_index), ]

all_patient_keys <- unique(all_meta[, c(
  "patient_id", "recurrence", "age", "gender", "compartment", "platform", "cell_context"
)])

tc_patient_keys <- unique(tc_meta[, c(
  "patient_id", "recurrence", "age", "gender", "compartment", "platform", "cell_context"
)])

summarise_matrix_by_keys <- function(expression_matrix, metadata_df, keys_df) {
  do.call(rbind, lapply(seq_len(nrow(keys_df)), function(i) {
    key <- keys_df[i, , drop = FALSE]
    columns <- metadata_df$author_column[
      metadata_df$patient_id == key$patient_id &
        metadata_df$compartment == key$compartment &
        metadata_df$cell_context == key$cell_context
    ]
    values <- expression_matrix[, columns, drop = FALSE]
  n_non_missing <- rowSums(!is.na(values))
  value_sum <- rowSums(values, na.rm = TRUE)
  value_sum[n_non_missing == 0] <- NA_real_
  value_mean <- rowMeans(values, na.rm = TRUE)
  value_mean[n_non_missing == 0] <- NA_real_
  value_median <- apply(values, 1, function(x) {
    if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
  })

  data.frame(
    gene = all_genes,
    .feature = gene_feature,
    is_times_biomarker = is_times_biomarker,
    key[rep(1, length(all_genes)), , drop = FALSE],
    n_segments_total = length(columns),
    n_non_missing = n_non_missing,
    value_sum = value_sum,
    value_mean = value_mean,
    value_median = value_median,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  }))
}

all_patient_summary <- summarise_matrix_by_keys(all_expression, all_meta, all_patient_keys)
rownames(all_patient_summary) <- NULL
all_patient_summary$n_segments_total <- as.integer(all_patient_summary$n_segments_total)
all_patient_summary$n_non_missing <- as.integer(all_patient_summary$n_non_missing)

tc_patient_summary <- summarise_matrix_by_keys(tc_expression, tc_meta, tc_patient_keys)
rownames(tc_patient_summary) <- NULL
tc_patient_summary$n_segments_total <- as.integer(tc_patient_summary$n_segments_total)
tc_patient_summary$n_non_missing <- as.integer(tc_patient_summary$n_non_missing)

times_segment_long <- tc_segment_long[tc_segment_long$is_times_biomarker, ]
times_patient_summary <- tc_patient_summary[tc_patient_summary$is_times_biomarker, ]

validation <- data.frame(
  metric = c(
    "source_file",
    "genes_expected",
    "genes_found",
    "genes_total",
    "all_roi_segments_total",
    "all_roi_segments_with_any_gene_value",
    "all_roi_matrix_dimensions",
    "compartment_counts",
    "tc_segments_total",
    "tc_segments_with_any_gene_value",
    "tc_segment_matrix_dimensions",
    "tc_segments_with_any_times_value",
    "patients_total",
    "patients_with_any_gene_value",
    "patients_with_any_times_value",
    "recurrence_levels",
    "value_scale_note"
  ),
  value = c(
    codeocean_input,
    paste(times_genes, collapse = ", "),
    paste(intersect(times_genes, rownames(gene_matrix)), collapse = ", "),
    length(all_genes),
    length(all_columns),
    sum(tapply(!is.na(all_segment_long$value), all_segment_long$author_column, any)),
    paste(dim(all_expression), collapse = " x "),
    paste(names(table(all_meta$compartment)), table(all_meta$compartment), sep = "=", collapse = "; "),
    length(tc_columns),
    sum(tapply(!is.na(tc_segment_long$value), tc_segment_long$author_column, any)),
    paste(dim(tc_expression), collapse = " x "),
    sum(tapply(!is.na(times_segment_long$value), times_segment_long$author_column, any)),
    length(unique(tc_meta$patient_id)),
    length(unique(tc_segment_long$patient_id[!is.na(tc_segment_long$value)])),
    length(unique(times_segment_long$patient_id[!is.na(times_segment_long$value)])),
    paste(names(table(tc_meta$recurrence)), table(tc_meta$recurrence), sep = "=", collapse = "; "),
    "Processed Code Ocean GeoMx ROI values; treated as author processed transcript reads, not raw Visium counts."
  ),
  stringsAsFactors = FALSE
)

all_matrix_csv <- file.path(processed_dir, "codeocean_GeoMx_all_gene_expression_matrix.csv")
all_metadata_csv <- file.path(processed_dir, "codeocean_GeoMx_all_roi_metadata.csv")
all_segment_csv <- file.path(processed_dir, "codeocean_GeoMx_all_segments_long_all_genes.csv")
all_patient_csv <- file.path(processed_dir, "codeocean_GeoMx_all_patient_compartment_summary_all_genes.csv")
matrix_csv <- file.path(processed_dir, "codeocean_GeoMx_TC_gene_expression_matrix.csv")
metadata_csv <- file.path(processed_dir, "codeocean_GeoMx_TC_roi_metadata.csv")
segment_all_csv <- file.path(processed_dir, "codeocean_GeoMx_TC_segments_long_all_genes.csv")
patient_all_csv <- file.path(processed_dir, "codeocean_GeoMx_TC_patient_summary_all_genes.csv")
segment_csv <- file.path(processed_dir, "codeocean_TIMES_GeoMx_TC_segments.csv")
patient_csv <- file.path(processed_dir, "codeocean_TIMES_GeoMx_TC_patient_summary.csv")
validation_csv <- file.path(processed_dir, "codeocean_TIMES_GeoMx_TC_validation.csv")
rds_path <- file.path(processed_dir, "TIMES_CodeOcean_GeoMx_TC_processed.rds")

write.csv(all_expression_df, all_matrix_csv, row.names = FALSE)
write.csv(all_meta, all_metadata_csv, row.names = FALSE)
write.csv(all_segment_long, all_segment_csv, row.names = FALSE)
write.csv(all_patient_summary, all_patient_csv, row.names = FALSE)
write.csv(tc_expression_df, matrix_csv, row.names = FALSE)
write.csv(tc_meta, metadata_csv, row.names = FALSE)
write.csv(tc_segment_long, segment_all_csv, row.names = FALSE)
write.csv(tc_patient_summary, patient_all_csv, row.names = FALSE)
write.csv(times_segment_long, segment_csv, row.names = FALSE)
write.csv(times_patient_summary, patient_csv, row.names = FALSE)
write.csv(validation, validation_csv, row.names = FALSE)
saveRDS(
  list(
    all_expression = all_expression,
    all_metadata = all_meta,
    all_segment_long = all_segment_long,
    all_patient_summary = all_patient_summary,
    tc_expression = tc_expression,
    tc_metadata = tc_meta,
    tc_segment_long = tc_segment_long,
    tc_patient_summary = tc_patient_summary,
    times_segment_long = times_segment_long,
    times_patient_summary = times_patient_summary,
    validation = validation,
    source_file = codeocean_input
  ),
  rds_path
)

message("Wrote: ", all_matrix_csv)
message("Wrote: ", all_metadata_csv)
message("Wrote: ", all_segment_csv)
message("Wrote: ", all_patient_csv)
message("Wrote: ", matrix_csv)
message("Wrote: ", metadata_csv)
message("Wrote: ", segment_all_csv)
message("Wrote: ", patient_all_csv)
message("Wrote: ", segment_csv)
message("Wrote: ", patient_csv)
message("Wrote: ", validation_csv)
message("Wrote: ", rds_path)
message("Full all-compartment expression matrix: ", paste(dim(all_expression), collapse = " x "))
message("Full TC expression matrix: ", paste(dim(tc_expression), collapse = " x "))
message("Patients with any TIMES TC value: ", validation$value[validation$metric == "patients_with_any_times_value"])
