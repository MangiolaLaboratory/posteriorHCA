#!/usr/bin/env Rscript
# Spatial raw-count pseudobulk as SummarizedExperiment.
# Section-level always; region-level only if verified author FF/DF labels exist.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(SummarizedExperiment)
  library(S4Vectors)
})

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
source(file.path(.script_dir, "utils.R"))

paths <- pak2_paths()
logf <- start_log(paths, "04_pseudobulk_spatial")
spatial_meta <- utils::read.csv(file.path(paths$metadata, "spatial_metadata.csv"), stringsAsFactors = FALSE)

obj_files <- file.path(paths$seurat_spatial, paste0(spatial_meta$section_id, ".rds"))
missing <- obj_files[!file.exists(obj_files)]
if (length(missing)) stop("Missing spatial Seurat objects:\n", paste(missing, collapse = "\n"))

objs <- lapply(obj_files, readRDS)
names(objs) <- spatial_meta$section_id

feat0 <- objs[[1]]@misc$feature_table
gene_ids <- rownames(feat0)
for (nm in names(objs)) {
  if (!identical(rownames(objs[[nm]]), gene_ids)) {
    stop("Spatial feature rownames differ for ", nm, "; refusing silent gene collapse")
  }
}

section_counts <- lapply(names(objs), function(nm) {
  raw <- get_raw_counts(objs[[nm]], assay = "Spatial")
  Matrix::rowSums(raw)
})
section_mat <- do.call(cbind, section_counts)
colnames(section_mat) <- names(objs)
section_mat <- methods::as(section_mat, "dgCMatrix")
if (length(section_mat@x)) section_mat@x <- round(section_mat@x)

col_df <- do.call(rbind, lapply(names(objs), function(nm) {
  obj <- objs[[nm]]
  raw <- get_raw_counts(obj, assay = "Spatial")
  md <- obj[[]]
  data.frame(
    donor_id = unique(md$donor_id),
    section_id = nm,
    sample_id = unique(md$sample_id),
    condition = unique(md$condition),
    disease = unique(md$disease),
    n_spots = ncol(obj),
    library_size = as.numeric(sum(raw)),
    n_features_detected = as.integer(sum(Matrix::rowSums(raw) > 0)),
    geo_accession = unique(md$geo_accession),
    gsm_accession = unique(md$gsm_accession),
    source_GSE = unique(md$source_GSE),
    tissue = unique(md$tissue),
    tissue_group = unique(md$tissue_group),
    stringsAsFactors = FALSE
  )
}))
rownames(col_df) <- col_df$section_id

row_df <- feat0
rownames(row_df) <- gene_ids

# Integrity: each section column sums to that object's total UMI; no spots dropped.
for (nm in names(objs)) {
  src <- sum(get_raw_counts(objs[[nm]], assay = "Spatial"))
  pb <- sum(section_mat[, nm])
  if (src != pb) stop("Section ", nm, " count-sum mismatch: ", src, " vs ", pb)
  if (col_df[nm, "n_spots"] != ncol(objs[[nm]])) {
    stop("Spot count changed during aggregation for ", nm)
  }
}

se_section <- make_summarized_experiment(
  counts = section_mat,
  row_data = row_df,
  col_data = col_df,
  metadata = list(
    aggregation = "section_id",
    note = "Raw integer UMI sums over all spots in each Visium section. No QC filtering."
  )
)
save_rds_atomic(se_section, file.path(paths$pb_spatial, "spatial_by_section_SE.rds"))
pak2_log("Wrote spatial_by_section_SE.rds with ", ncol(se_section), " sections.", .log_file = logf)

labelled <- vapply(objs, function(o) sum(!is.na(o$pathological_region) & nzchar(o$pathological_region)), integer(1))
if (sum(labelled) == 0L) {
  pak2_log("No verified pathological_region labels; not writing spatial_by_region_SE.rds.", .log_file = logf)
} else {
  pak2_log("Building region-level pseudobulk from author labels.", .log_file = logf)
  region_cols <- list()
  region_cd <- list()
  unlabelled_report <- list()
  for (nm in names(objs)) {
    obj <- objs[[nm]]
    raw <- get_raw_counts(obj, assay = "Spatial")
    lab <- as.character(obj$pathological_region)
    keep <- !is.na(lab) & nzchar(lab)
    n_drop <- sum(!keep)
    unlabelled_report[[nm]] <- data.frame(section_id = nm, n_spots = ncol(obj), n_labelled = sum(keep), n_unlabelled = n_drop)
    if (!any(keep)) next
    pb <- aggregate_counts_by_group(raw[, keep, drop = FALSE], lab[keep])
    regions <- colnames(pb)
    ids <- paste(nm, regions, sep = "__")
    colnames(pb) <- ids
    region_cols[[nm]] <- pb
    md <- obj[[]]
    region_cd[[nm]] <- data.frame(
      pseudobulk_id = ids,
      donor_id = unique(md$donor_id),
      section_id = nm,
      sample_id = unique(md$sample_id),
      pathological_region = regions,
      condition = unique(md$condition),
      disease = unique(md$disease),
      n_spots = as.integer(table(lab[keep])[regions]),
      library_size = as.numeric(Matrix::colSums(pb)),
      geo_accession = unique(md$geo_accession),
      gsm_accession = unique(md$gsm_accession),
      source_GSE = unique(md$source_GSE),
      tissue = unique(md$tissue),
      tissue_group = unique(md$tissue_group),
      n_unlabelled_spots_in_section = n_drop,
      stringsAsFactors = FALSE
    )
    if (sum(raw[, keep, drop = FALSE]) != sum(pb)) {
      stop("Region aggregation count mismatch for ", nm)
    }
  }
  region_mat <- do.call(cbind, region_cols)
  region_cd <- do.call(rbind, region_cd)
  rownames(region_cd) <- region_cd$pseudobulk_id
  se_region <- make_summarized_experiment(
    counts = region_mat,
    row_data = row_df,
    col_data = region_cd,
    metadata = list(
      aggregation = "donor_id x section_id x pathological_region",
      unlabelled = do.call(rbind, unlabelled_report),
      note = "Only author-provided region labels. Unlabelled spots are not assigned FF/DF."
    )
  )
  save_rds_atomic(se_region, file.path(paths$pb_spatial, "spatial_by_region_SE.rds"))
  utils::write.csv(do.call(rbind, unlabelled_report), file.path(paths$logs, "spatial_region_unlabelled.csv"), row.names = FALSE)
}

utils::write.csv(col_df, file.path(paths$logs, "spatial_section_pseudobulk_coldata.csv"), row.names = FALSE)
pak2_log("04_pseudobulk_spatial.R complete.", .log_file = logf)
