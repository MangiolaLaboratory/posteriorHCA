#!/usr/bin/env Rscript
# Reconstruct one Seurat object per Visium section from GEO processed files.
# No QC, normalisation, integration, or clustering.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
source(file.path(.script_dir, "utils.R"))

paths <- pak2_paths()
logf <- start_log(paths, "02_build_spatial_seurat")
spatial_meta <- utils::read.csv(file.path(paths$metadata, "spatial_metadata.csv"), stringsAsFactors = FALSE)
availability <- readRDS(file.path(paths$metadata, "annotation_availability.rds"))
xlsx_contents <- readRDS(file.path(paths$metadata, "supplement_xlsx_parsed.rds"))

# Optional author FF/DF labels: require barcode + region columns.
region_map <- NULL
if (isTRUE(availability$ff_df_spot_annotations) && length(xlsx_contents)) {
  for (src in names(xlsx_contents)) {
    for (sh in names(xlsx_contents[[src]])) {
      dat <- xlsx_contents[[src]][[sh]]$data
      if (is.null(dat) || !nrow(dat)) next
      bcol <- names(dat)[vapply(dat, looks_like_barcode_column, logical(1))]
      if (!length(bcol)) next
      rcol <- names(dat)[grepl("ff|df|fibroblastic|dense.?fibr|patholog|region|annot", names(dat), ignore.case = TRUE)]
      if (!length(rcol)) next
      tmp <- data.frame(
        barcode = as.character(dat[[bcol[[1]]]]),
        pathological_region = as.character(dat[[rcol[[1]]]]),
        section_hint = if ("section_id" %in% names(dat)) as.character(dat$section_id) else NA_character_,
        source = paste(src, sh, sep = ":"),
        stringsAsFactors = FALSE
      )
      tmp <- tmp[!is.na(tmp$barcode) & nzchar(tmp$barcode), ]
      region_map <- rbind(region_map, tmp)
    }
  }
}

build_one_section <- function(row) {
  section <- row$section_id
  out_rds <- file.path(paths$seurat_spatial, paste0(section, ".rds"))
  raw_dir <- file.path(paths$raw_gse285246, "spatial", section)
  files <- list(
    barcodes = list.files(raw_dir, pattern = "barcodes\\.tsv\\.gz$", full.names = TRUE),
    features = list.files(raw_dir, pattern = "features\\.tsv\\.gz$", full.names = TRUE),
    matrix = list.files(raw_dir, pattern = "matrix\\.mtx\\.gz$", full.names = TRUE),
    positions = list.files(raw_dir, pattern = "tissue_positions_list\\.csv\\.gz$", full.names = TRUE),
    scalefactors = list.files(raw_dir, pattern = "scalefactors_json\\.json\\.gz$", full.names = TRUE),
    hires = list.files(raw_dir, pattern = "tissue_hires_image\\.png\\.gz$", full.names = TRUE),
    lowres = list.files(raw_dir, pattern = "tissue_lowres_image\\.png\\.gz$", full.names = TRUE)
  )
  if (any(lengths(files) != 1L)) {
    stop("Missing spatial files for ", section, ": ", paste(names(files)[lengths(files) != 1L], collapse = ", "))
  }
  files <- lapply(files, `[[`, 1L)
  prep <- file.path(paths$prepared, "spatial", section)
  prepare_visium_dir(prep, files)

  counts <- read_10x_counts(file.path(prep, "filtered_feature_bc_matrix"), gene_column = 1L)
  feat <- feature_table_from_10x(files$features, rownames(counts))
  obj <- CreateSeuratObject(
    counts = counts,
    assay = "Spatial",
    project = section,
    min.cells = 0,
    min.features = 0
  )
  image_ok <- FALSE
  img <- tryCatch({
    Seurat::Read10X_Image(
      image.dir = file.path(prep, "spatial"),
      image.name = "tissue_lowres_image.png",
      assay = "Spatial",
      slice = section,
      filter.matrix = FALSE,
      image.type = "VisiumV1"
    )
  }, error = function(e) {
    pak2_log("Read10X_Image VisiumV1 failed for ", section, ": ", conditionMessage(e), .log_file = logf)
    tryCatch(
      Seurat::Read10X_Image(
        image.dir = file.path(prep, "spatial"),
        image.name = "tissue_lowres_image.png",
        assay = "Spatial",
        slice = section,
        filter.matrix = FALSE,
        image.type = "VisiumV2"
      ),
      error = function(e2) {
        pak2_log("Read10X_Image VisiumV2 failed for ", section, ": ", conditionMessage(e2), .log_file = logf)
        NULL
      }
    )
  })
  if (!is.null(img)) {
    common <- intersect(Cells(obj), Cells(img))
    if (!length(common)) {
      pak2_log("No overlapping barcodes between matrix and image for ", section, .log_file = logf)
    } else {
      img <- img[Cells(obj)]
      obj[[section]] <- img
      image_ok <- TRUE
    }
  }

  pos <- read_tissue_positions(file.path(prep, "spatial", "tissue_positions_list.csv"))
  pos <- pos[match(colnames(obj), pos$barcode), ]
  n <- ncol(obj)
  meta <- standard_spot_or_cell_meta(n, list(
    donor_id = row$donor_id,
    sample_id = row$sample_id,
    section_id = row$section_id,
    condition = row$condition,
    disease = row$disease,
    source_GSE = row$source_GSE,
    geo_accession = row$geo_accession,
    gsm_accession = row$gsm_accession,
    tissue = row$tissue,
    tissue_group = row$tissue_group,
    sex = row$sex,
    age = row$age,
    smoking_status = row$smoking_status,
    smoking_cessation = row$smoking_cessation,
    pack_years = row$pack_years,
    bmi = row$bmi,
    library_chemistry = row$library_chemistry,
    genome_assembly = row$genome_assembly,
    sequencer = row$sequencer,
    processing_pipeline = row$processing_pipeline,
    pathological_region = NA_character_,
    pathological_region_source = NA_character_,
    in_tissue = pos$in_tissue,
    array_row = pos$array_row,
    array_col = pos$array_col,
    pxl_row_in_fullres = pos$pxl_row_in_fullres,
    pxl_col_in_fullres = pos$pxl_col_in_fullres,
    spatial_image_attached = image_ok
  ))
  rownames(meta) <- colnames(obj)
  obj <- AddMetaData(obj, metadata = meta)

  if (!is.null(region_map)) {
    rm <- region_map
    if (!all(is.na(rm$section_hint))) {
      rm_sec <- rm[!is.na(rm$section_hint) & rm$section_hint %in% c(section, row$sample_id, row$donor_id), ]
      if (nrow(rm_sec)) rm <- rm_sec
    }
    m <- match(colnames(obj), rm$barcode)
    if (all(is.na(m))) {
      m <- match(normalise_barcode(colnames(obj)), normalise_barcode(rm$barcode))
    }
    if (any(!is.na(m))) {
      obj$pathological_region <- rm$pathological_region[m]
      obj$pathological_region_source <- rm$source[m]
      pak2_log(
        section, ": attached pathological_region for ",
        sum(!is.na(m)), "/", ncol(obj), " spots from ", unique(na.omit(rm$source))[1],
        .log_file = logf
      )
    }
  }

  raw <- get_raw_counts(obj, assay = "Spatial")
  obj$n_features_detected <- as.integer(Matrix::colSums(raw > 0))
  obj$library_size_spot <- as.numeric(Matrix::colSums(raw))
  obj@misc$feature_table <- feat
  obj@misc$no_qc_note <- "Unfiltered GEO processed counts; no cells/spots/genes removed."
  save_rds_atomic(obj, out_rds)

  list(
    section_id = section,
    donor_id = row$donor_id,
    n_spots = ncol(obj),
    n_features = nrow(obj),
    total_UMI = as.numeric(sum(raw)),
    image_attached = image_ok,
    n_pathological_region_labelled = sum(!is.na(obj$pathological_region)),
    rds = out_rds
  )
}

reports <- lapply(seq_len(nrow(spatial_meta)), function(i) {
  row <- spatial_meta[i, ]
  pak2_log("Building spatial Seurat: ", row$section_id, .log_file = logf)
  build_one_section(row)
})
report_df <- do.call(rbind, lapply(reports, function(x) data.frame(x, stringsAsFactors = FALSE)))
write.csv(report_df, file.path(paths$logs, "spatial_seurat_validation.csv"), row.names = FALSE)

objs <- lapply(report_df$rds, readRDS)
names(objs) <- report_df$section_id
save_rds_atomic(objs, file.path(paths$seurat_spatial, "spatial_seurat_list.rds"))

pak2_log("Spatial sections: ", paste(sprintf("%s n_spots=%s total_UMI=%s", report_df$section_id, report_df$n_spots, report_df$total_UMI), collapse = "; "), .log_file = logf)
pak2_log("02_build_spatial_seurat.R complete.", .log_file = logf)
