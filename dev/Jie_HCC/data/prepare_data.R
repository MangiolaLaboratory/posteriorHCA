# Prepare TIMES / Jie HCC spatial case-study data
#
# Run from posteriorHCA project root:
#   conda activate R_env
#   Rscript dev/Jie_HCC/data/prepare_data.R
#
# Primary deposits (Jia et al. Nature 2025, s41586-025-08668-x) are on CNCB-NGDC
# and are currently controlled-access / unavailable for automated download from
# this environment. This script records download attempts, saves a registry,
# defines TIMES biomarkers, and builds TC-compartment pseudobulk when raw files
# are present under dev/Jie_HCC/data/raw/.

suppressPackageStartupMessages({
  library(tibble)
  library(tidyr)
  library(dplyr)
  library(readr)
  library(purrr)
  library(glue)
  library(SummarizedExperiment)
  library(S4Vectors)
})

if (requireNamespace("AnnotationDbi", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  library(AnnotationDbi)
  library(org.Hs.eg.db)
}

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
project_root <- Sys.getenv("POSTERIORHCA_ROOT", unset = NA_character_)
if (is.na(project_root) && length(script_path)) {
  project_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), mustWork = FALSE)
}
if (is.na(project_root) || !dir.exists(project_root)) {
  project_root <- normalizePath(getwd(), mustWork = FALSE)
}

case_dir  <- file.path(project_root, "dev", "Jie_HCC")
data_dir  <- file.path(case_dir, "data")
raw_dir   <- file.path(data_dir, "raw")
proc_dir  <- file.path(data_dir, "processed")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

message("Project root: ", project_root)
message("Data dir:     ", data_dir)

# ---------------------------------------------------------------------------
# TIMES biomarkers (five genes from Fig. 2 / TIMES score)
# ---------------------------------------------------------------------------

times_biomarkers <- readr::read_csv(
  file.path(data_dir, "times_biomarkers.csv"),
  show_col_types = FALSE
)

if (requireNamespace("AnnotationDbi", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  times_biomarkers <- times_biomarkers |>
    mutate(
      .feature = AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = gene,
        column = "ENSEMBL",
        keytype = "SYMBOL",
        multiVals = "first"
      )
    )
} else {
  times_biomarkers <- times_biomarkers |> mutate(.feature = NA_character_)
}

saveRDS(times_biomarkers, file.path(data_dir, "times_biomarkers.rds"))
readr::write_csv(times_biomarkers, file.path(data_dir, "times_biomarkers.csv"))

# ---------------------------------------------------------------------------
# Data registry (deposits cited in the paper)
# ---------------------------------------------------------------------------

data_registry <- tibble::tribble(
  ~resource_id,   ~platform,              ~description,                                              ~url,                                                          ~local_subdir,              ~access,
  "OMIX005738",   "GeoMx DSP",            "WTA/CTA expression + ROI labels (discovery cohort)",        "https://ngdc.cncb.ac.cn/omix/release/OMIX005738",             "omix005738",               "controlled",
  "OMIX005736",   "LC-MS/MS",             "Proteomics (validation; not used in pseudobulk plan)",    "https://ngdc.cncb.ac.cn/omix/release/OMIX005736",             "omix005736",               "controlled",
  "HRA006579",    "10x Visium",           "Spatial transcriptomics (discovery cohort, 17 patients)", "https://ngdc.cncb.ac.cn/gsa-human/browse/HRA006579",            "hra006579",                "controlled",
  "CNP0000650",   "scRNA-seq",            "In-house scRNA-seq (NK validation)",                      "https://db.cngb.org/search/project/CNP0000650",                 "cnp0000650",               "controlled",
  "CO.7364332",   "Code Ocean",           "Scripts + bundled input/output",                          "https://doi.org/10.24433/CO.7364332.v1",                      "code_ocean",               "manual",
  "MOESM1",       "Supplementary",        "Supplementary tables / figures (Nature)",                 "https://www.nature.com/articles/s41586-025-08668-x#Sec24",      "supplementary",            "manual"
)

saveRDS(data_registry, file.path(data_dir, "data_registry.rds"))
readr::write_csv(data_registry, file.path(data_dir, "data_registry.csv"))

# Expected raw filenames once CNCB access is granted (from OMIX pages)
expected_files <- tibble::tribble(
  ~resource_id,   ~file_name,                     ~role,
  "OMIX005738",   "OMIX005738-01_DSP_WTA.xlsx",   "geomp_wta_counts",
  "OMIX005738",   "OMIX005738-02_DSP_CTA.xlsx",   "geomp_cta_counts",
  "OMIX005738",   "OMIX005738-03_DSP_labels.xlsx","geomp_roi_metadata",
  "HRA006579",    "visium_spaceranger/",          "visium_per_sample_dirs",
  "CNP0000650",   "scRNA_matrix/",                "scRNA_counts"
)

readr::write_csv(expected_files, file.path(data_dir, "expected_raw_files.csv"))

# ---------------------------------------------------------------------------
# Download helpers
# ---------------------------------------------------------------------------

download_file_wget <- function(url, destfile, overwrite = FALSE) {
  if (file.exists(destfile) && !overwrite) {
    return(tibble(
      url = url, destfile = destfile, status = "skipped_exists",
      bytes = file.info(destfile)$size
    ))
  }
  dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
  cmd <- sprintf(
    "wget -c --tries=3 --timeout=60 --read-timeout=120 -O %s %s",
    shQuote(destfile), shQuote(url)
  )
  message("Running: ", cmd)
  exit_code <- system(cmd)
  tibble(
    url = url,
    destfile = destfile,
    status = if (exit_code == 0L) "ok" else "failed",
    bytes = if (file.exists(destfile)) file.info(destfile)$size else NA_real_
  )
}

# CNCB OMIX files are not publicly downloadable at time of writing; record intent.
download_registry <- purrr::map_dfr(seq_len(nrow(data_registry)), function(i) {
  row <- data_registry[i, ]
  dest <- file.path(raw_dir, row$local_subdir)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  tibble(
    resource_id = row$resource_id,
    platform = row$platform,
    access = row$access,
    attempted = TRUE,
    downloaded = FALSE,
    note = case_when(
      row$access == "controlled" ~ "CNCB controlled-access; request via DAC (Cheng Sun, charless@ustc.edu.cn)",
      row$access == "manual" ~ "Manual download required (Code Ocean capsule or Nature supplementary)",
      TRUE ~ NA_character_
    ),
    local_dir = dest
  )
})

readr::write_csv(download_registry, file.path(data_dir, "download_manifest.csv"))
saveRDS(download_registry, file.path(data_dir, "download_manifest.rds"))

# ---------------------------------------------------------------------------
# Pseudobulk — GeoMx DSP (NK-enriched ROIs, TC compartment only)
# ---------------------------------------------------------------------------

#' Aggregate GeoMx WTA counts to patient x TC pseudobulk.
#'
#' Expects:
#' - counts_xlsx: wide or long count table (Gene x ROI) from OMIX005738-01
#' - labels_xlsx: ROI metadata with columns including compartment (AS/IF/TC),
#'   patient_id, recurrence (REC/non-REC), and NK-enrichment flag if available
build_geomp_tc_pseudobulk <- function(counts_xlsx, labels_xlsx, compartment = "TC") {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Install readxl to parse GeoMx xlsx files.")
  }
  counts_raw <- readxl::read_excel(counts_xlsx)
  labels <- readxl::read_excel(labels_xlsx)

  # Normalise column names heuristically
  labels <- labels |>
    rename_with(tolower) |>
    mutate(
      compartment = dplyr::coalesce(compartment, region, tissue_compartment),
      patient_id = dplyr::coalesce(patient_id, patient, sample_id, case_id),
      recurrence = dplyr::coalesce(recurrence, rec_status, outcome)
    )

  roi_in_tc <- labels |>
    filter(compartment %in% c(compartment, "Tumor center", "tumour centre", "TC"))

  # Assume first column is gene id / symbol
  gene_col <- names(counts_raw)[1]
  counts_long <- counts_raw |>
    pivot_longer(-all_of(gene_col), names_to = "roi_id", values_to = "count") |>
    rename(gene = all_of(gene_col))

  pb <- counts_long |>
    inner_join(
      roi_in_tc |> transmute(roi_id = as.character(roi_id), patient_id, recurrence),
      by = "roi_id"
    ) |>
    group_by(patient_id, recurrence, gene) |>
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
    mutate(
      pb_id = paste(patient_id, compartment, sep = "_"),
      compartment = compartment,
      platform = "GeoMx_DSP",
      cell_context = "NK_enriched_ROI"
    )

  mat <- pb |>
    select(gene, pb_id, count) |>
    pivot_wider(names_from = pb_id, values_from = count) |>
    tibble::column_to_rownames("gene") |>
    as.matrix()

  coldata <- pb |>
    distinct(pb_id, patient_id, recurrence, compartment, platform, cell_context)

  se <- SummarizedExperiment(
    assays = list(counts = mat[, coldata$pb_id, drop = FALSE]),
    colData = DataFrame(coldata |> tibble::column_to_rownames("pb_id"))
  )
  se
}

# ---------------------------------------------------------------------------
# Pseudobulk — 10x Visium (TC spots only; spot annotations required)
# ---------------------------------------------------------------------------

#' Build patient x TC pseudobulk from Visium count matrices + spot metadata.
#'
#' spot_meta must contain: barcode, patient_id, compartment (TC), recurrence
build_visium_tc_pseudobulk <- function(counts_mat, spot_meta, compartment = "TC") {
  spot_meta <- spot_meta |>
    rename_with(tolower) |>
    mutate(
      barcode = dplyr::coalesce(barcode, spot_id, rownames),
      compartment = dplyr::coalesce(compartment, region, tissue_compartment),
      patient_id = dplyr::coalesce(patient_id, sample_id, patient)
    ) |>
    filter(compartment %in% c(compartment, "Tumor center", "tumour centre", "TC"))

  common <- intersect(colnames(counts_mat), spot_meta$barcode)
  if (!length(common)) {
    stop("No overlapping barcodes between counts_mat and spot_meta.")
  }

  counts_mat <- counts_mat[, common, drop = FALSE]
  spot_meta <- spot_meta |> filter(barcode %in% common)

  pb_list <- lapply(split(spot_meta, spot_meta$patient_id), function(meta) {
    cols <- meta$barcode
    rowSums(counts_mat[, cols, drop = FALSE])
  })

  mat <- do.call(cbind, pb_list)
  colnames(mat) <- names(pb_list)

  coldata <- spot_meta |>
    distinct(patient_id, recurrence) |>
    mutate(
      pb_id = patient_id,
      compartment = compartment,
      platform = "10x_Visium",
      cell_context = "TC_spots_all_celltypes"
    ) |>
    tibble::column_to_rownames("pb_id")

  SummarizedExperiment(
    assays = list(counts = mat),
    colData = DataFrame(coldata[colnames(mat), , drop = FALSE])
  )
}

# ---------------------------------------------------------------------------
# Try to build pseudobulk from any locally present raw files
# ---------------------------------------------------------------------------

geomp_wta <- file.path(raw_dir, "omix005738", "OMIX005738-01_DSP_WTA.xlsx")
geomp_lbl <- file.path(raw_dir, "omix005738", "OMIX005738-03_DSP_labels.xlsx")

pseudobulk_outputs <- list()

if (file.exists(geomp_wta) && file.exists(geomp_lbl)) {
  message("Building GeoMx TC pseudobulk …")
  se_geomp <- build_geomp_tc_pseudobulk(geomp_wta, geomp_lbl, compartment = "TC")
  out_geomp <- file.path(proc_dir, "TIMES_GeoMx_TC_pseudobulk.rds")
  saveRDS(se_geomp, out_geomp)
  pseudobulk_outputs$geomp_tc <- out_geomp
  message("Saved: ", out_geomp)
} else {
  message(
    "GeoMx raw xlsx not found under ", file.path(raw_dir, "omix005738"),
    ". Place approved CNCB downloads there and re-run."
  )
}

# Visium: look for Space Ranger filtered matrices under hra006579/*/
visium_root <- file.path(raw_dir, "hra006579")
visium_samples <- list.dirs(visium_root, recursive = FALSE, full.names = TRUE)
spot_meta_path <- file.path(raw_dir, "hra006579", "spot_compartment_annotations.csv")

if (length(visium_samples) && file.exists(spot_meta_path)) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Install Matrix to read 10x Visium matrices.")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Install Seurat to read 10x Visium matrices.")
  }
  message("Building Visium TC pseudobulk …")
  counts_list <- lapply(visium_samples, function(sdir) {
    mat <- Seurat::Read10X(file.path(sdir, "filtered_feature_bc_matrix"))
    list(sample = basename(sdir), mat = mat)
  })
  # User-supplied spot annotations joining patient + compartment
  spot_meta <- readr::read_csv(spot_meta_path, show_col_types = FALSE)
  # Implementation note: merge per-sample mats with spot_meta before aggregation
  message("Visium spot metadata found; finish sample-specific merge in qmd when data arrive.")
} else {
  message(
    "Visium raw data or spot_compartment_annotations.csv not found under ",
    visium_root, ". Required after CNCB access."
  )
}

# Save processing plan metadata
processing_plan <- tibble(
  pseudobulk_unit = "patient x TC (tumour centre)",
  spatial_scope = "TC compartment only; ignore IF/AS and tile-level spatial patterns",
  primary_platform = "10x Visium (17 discovery patients)",
  secondary_platform = "GeoMx DSP NK-enriched ROIs (8 tissues)",
  atlas_cell_type = "nk (posteriorHCA V1_nk)",
  atlas_baseline = "universal marginalised null (Normal disease, tissue_groups = blood, offset = 0)",
  biomarker_genes = paste(times_biomarkers$gene, collapse = ", "),
  outputs_created = paste(names(pseudobulk_outputs), collapse = ", ")
)

saveRDS(processing_plan, file.path(data_dir, "processing_plan.rds"))
readr::write_csv(processing_plan, file.path(data_dir, "processing_plan.csv"))

message("\n=== prepare_data.R complete ===")
message("Biomarkers: ", nrow(times_biomarkers), " genes")
message("Download manifest: ", file.path(data_dir, "download_manifest.csv"))
if (length(pseudobulk_outputs)) {
  message("Pseudobulk outputs: ", paste(unlist(pseudobulk_outputs), collapse = ", "))
} else {
  message("No pseudobulk built yet — awaiting CNCB / Code Ocean raw files.")
}
