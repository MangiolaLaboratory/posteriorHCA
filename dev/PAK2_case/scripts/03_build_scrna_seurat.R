#!/usr/bin/env Rscript
# Reconstruct one Seurat object per scRNA donor from GEO processed MTX files.
# No QC, normalisation, integration, clustering, or de novo annotation.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(readxl)
})

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
source(file.path(.script_dir, "utils.R"))

paths <- pak2_paths()
logf <- start_log(paths, "03_build_scrna_seurat")
scrna_meta <- utils::read.csv(file.path(paths$metadata, "scrna_metadata.csv"), stringsAsFactors = FALSE)
xlsx_contents <- readRDS(file.path(paths$metadata, "supplement_xlsx_parsed.rds"))
availability <- readRDS(file.path(paths$metadata, "annotation_availability.rds"))

pick_column <- function(dat, patterns) {
  cols <- names(dat)
  for (p in patterns) {
    hit <- cols[grepl(p, cols, ignore.case = TRUE)]
    if (length(hit)) return(hit[[1]])
  }
  NA_character_
}

# Collect barcode-level author annotations from scanned tables.
cell_annot <- NULL
if (length(xlsx_contents)) {
  for (src in names(xlsx_contents)) {
    for (sh in names(xlsx_contents[[src]])) {
      dat <- xlsx_contents[[src]][[sh]]$data
      if (is.null(dat) || !nrow(dat)) next
      bcol <- names(dat)[vapply(dat, looks_like_barcode_column, logical(1))]
      if (!length(bcol)) {
        bcol <- names(dat)[grepl("barcode|cell_id|cellid", names(dat), ignore.case = TRUE)]
      }
      if (!length(bcol)) next
      ct <- pick_column(dat, c("cell_type_original", "^celltype$", "cell_type", "cell.annotation"))
      broad <- pick_column(dat, c("cell_type_broad", "major_cell", "compartment", "lineage"))
      donor_col <- pick_column(dat, c("donor_id", "orig.ident", "orig_ident", "sample", "patient", "donor"))
      tmp <- data.frame(
        barcode_raw = as.character(dat[[bcol[[1]]]]),
        cell_type_original = if (!is.na(ct)) as.character(dat[[ct]]) else NA_character_,
        cell_type_broad = if (!is.na(broad)) as.character(dat[[broad]]) else NA_character_,
        donor_hint = if (!is.na(donor_col)) as.character(dat[[donor_col]]) else NA_character_,
        annotation_source = paste(src, sh, sep = ":"),
        stringsAsFactors = FALSE
      )
      tmp <- tmp[!is.na(tmp$barcode_raw) & nzchar(tmp$barcode_raw), ]
      if (!nrow(tmp)) next
      if (all(is.na(tmp$cell_type_original))) next
      cell_annot <- rbind(cell_annot, tmp)
    }
  }
}
if (!is.null(cell_annot)) {
  write.csv(cell_annot, file.path(paths$metadata, "author_cell_annotations_raw.csv"), row.names = FALSE)
  pak2_log("Loaded ", nrow(cell_annot), " author annotation rows from supplements.", .log_file = logf)
} else {
  pak2_log("No barcode-level scRNA cell-type annotations found in supplements.", .log_file = logf)
}

copd_meta_path <- file.path(paths$raw_gse173896, "GSE173896_COPD_meta_new.xlsx")
copd_meta <- NULL
if (file.exists(copd_meta_path)) {
  copd_meta <- as.data.frame(readxl::read_excel(copd_meta_path, sheet = 1, col_types = "text"), stringsAsFactors = FALSE)
  names(copd_meta)[1] <- "barcode"
  pak2_log("Loaded COPD study metadata with ", nrow(copd_meta), " rows; Class is a disease group, not a cell type.", .log_file = logf)
}

raw_dir_for <- function(row) {
  if (row$source_GSE == "GSE285246") {
    file.path(paths$raw_gse285246, "scrna", row$donor_id)
  } else {
    file.path(paths$raw_gse173896, "scrna", row$donor_id)
  }
}

build_one_donor <- function(row) {
  donor <- row$donor_id
  out_rds <- file.path(paths$seurat_scrna, paste0(donor, ".rds"))
  raw_dir <- raw_dir_for(row)
  barcodes <- list.files(raw_dir, pattern = "barcodes\\.tsv\\.gz$", full.names = TRUE)
  features <- list.files(raw_dir, pattern = "features\\.tsv\\.gz$", full.names = TRUE)
  matrix <- list.files(raw_dir, pattern = "matrix\\.mtx\\.gz$", full.names = TRUE)
  if (length(barcodes) != 1L || length(features) != 1L || length(matrix) != 1L) {
    stop("Missing scRNA MTX files for ", donor, " in ", raw_dir)
  }
  prep <- file.path(paths$prepared, "scrna", donor)
  prepare_10x_mtx_dir(prep, barcodes[[1]], features[[1]], matrix[[1]])
  counts <- read_10x_counts(prep, gene_column = 1L)
  feat <- feature_table_from_10x(features[[1]], rownames(counts))
  obj <- CreateSeuratObject(
    counts = counts,
    assay = "RNA",
    project = donor,
    min.cells = 0,
    min.features = 0
  )
  n <- ncol(obj)
  meta <- standard_spot_or_cell_meta(n, list(
    donor_id = row$donor_id,
    sample_id = row$sample_id,
    section_id = NA_character_,
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
    cell_type_original = NA_character_,
    cell_type_broad = NA_character_,
    annotation_source = NA_character_
  ))
  rownames(meta) <- colnames(obj)
  obj <- AddMetaData(obj, metadata = meta)

  n_matched <- 0L
  if (!is.null(cell_annot)) {
    ca <- cell_annot
    if (!all(is.na(ca$donor_hint))) {
      ca_d <- ca[!is.na(ca$donor_hint) & grepl(paste0("(^|_)", donor, "(_|$)"), ca$donor_hint)]
      if (nrow(ca_d)) ca <- ca_d
    }
    cells <- colnames(obj)
    m <- match(cells, ca$barcode_raw)
    if (all(is.na(m))) m <- match(normalise_barcode(cells), normalise_barcode(ca$barcode_raw))
    if (all(is.na(m))) {
      m <- match(paste0(donor, "_", cells), ca$barcode_raw)
    }
    if (any(!is.na(m))) {
      obj$cell_type_original <- ca$cell_type_original[m]
      obj$cell_type_broad <- ca$cell_type_broad[m]
      obj$annotation_source <- ca$annotation_source[m]
      n_matched <- sum(!is.na(m) & !is.na(obj$cell_type_original) & nzchar(obj$cell_type_original))
      pak2_log(donor, ": matched author cell-type annotations for ", n_matched, "/", n, " cells", .log_file = logf)
    }
  }

  if (!is.null(copd_meta) && donor %in% copd_meta$orig.ident) {
    cm <- copd_meta[copd_meta$orig.ident == donor, , drop = FALSE]
    m <- match(normalise_barcode(colnames(obj)), normalise_barcode(cm$barcode))
    obj$copd_study_class <- cm$Class[m]
    if ("seurat_clusters" %in% names(cm)) obj$copd_study_seurat_clusters <- cm$seurat_clusters[m]
    if ("DoubletFinder" %in% names(cm)) obj$copd_study_DoubletFinder <- cm$DoubletFinder[m]
    pak2_log(
      donor, ": attached COPD-study metadata for ", sum(!is.na(m)), "/", n,
      " cells (Class is disease group, not cell type).",
      .log_file = logf
    )
  }

  raw <- get_raw_counts(obj, assay = "RNA")
  obj$n_features_detected <- as.integer(Matrix::colSums(raw > 0))
  obj$library_size_cell <- as.numeric(Matrix::colSums(raw))
  obj@misc$feature_table <- feat
  obj@misc$no_qc_note <- "Unfiltered GEO processed counts; no cells/genes removed."
  save_rds_atomic(obj, out_rds)

  list(
    donor_id = donor,
    condition = row$condition,
    geo_accession = row$geo_accession,
    n_cells = ncol(obj),
    n_features = nrow(obj),
    total_UMI = as.numeric(sum(raw)),
    n_annotated = as.integer(n_matched),
    n_unannotated = as.integer(n - n_matched),
    feature_id_class = if (grepl("^ENSG", rownames(obj)[1])) "Ensembl" else "other",
    rds = out_rds,
    features_path = features[[1]]
  )
}

reports <- lapply(seq_len(nrow(scrna_meta)), function(i) {
  row <- scrna_meta[i, ]
  pak2_log("Building scRNA Seurat: ", row$donor_id, .log_file = logf)
  build_one_donor(row)
})
report_df <- do.call(rbind, lapply(reports, function(x) {
  data.frame(
    donor_id = x$donor_id,
    condition = x$condition,
    geo_accession = x$geo_accession,
    n_cells = x$n_cells,
    n_features = x$n_features,
    total_UMI = x$total_UMI,
    n_annotated = x$n_annotated,
    n_unannotated = x$n_unannotated,
    feature_id_class = x$feature_id_class,
    rds = x$rds,
    features_path = x$features_path,
    stringsAsFactors = FALSE
  )
}))
write.csv(report_df, file.path(paths$logs, "scrna_seurat_validation.csv"), row.names = FALSE)

objs <- lapply(report_df$rds, readRDS)
names(objs) <- report_df$donor_id
save_rds_atomic(objs, file.path(paths$seurat_scrna, "scrna_seurat_list.rds"))

# Count-level merge for convenience; preserve raw counts and donor IDs.
pak2_log("Merging six donor objects at count level...", .log_file = logf)
merged <- merge(
  x = objs[[1]],
  y = objs[-1],
  add.cell.ids = names(objs),
  project = "PAK2_scrna",
  merge.data = FALSE
)
if (requireNamespace("SeuratObject", quietly = TRUE) && isTRUE(is(merged[["RNA"]], "Assay5"))) {
  merged[["RNA"]] <- SeuratObject::JoinLayers(merged[["RNA"]])
}
save_rds_atomic(merged, file.path(paths$seurat_scrna, "scrna_seurat_merged_counts.rds"))

# Feature map across studies (do not subset original objects).
feat_tables <- lapply(seq_len(nrow(report_df)), function(i) {
  ft <- read_features_tsv(report_df$features_path[[i]])
  ft$donor_id <- report_df$donor_id[[i]]
  ft$geo_accession <- report_df$geo_accession[[i]]
  ft
})
all_feat <- do.call(rbind, feat_tables)
id_sets <- split(all_feat$gene_id, all_feat$geo_accession)
id_sets <- lapply(id_sets, unique)
universe <- unique(unlist(id_sets, use.names = FALSE))
feature_map <- data.frame(
  gene_id = universe,
  stringsAsFactors = FALSE
)
for (gse in names(id_sets)) {
  feature_map[[paste0("present_", gse)]] <- universe %in% id_sets[[gse]]
}
sym_map <- unique(all_feat[, c("gene_id", "gene_symbol", "feature_type")])
# If a gene_id maps to multiple symbols across files, keep them all concatenated.
sym_agg <- stats::aggregate(
  cbind(gene_symbol, feature_type) ~ gene_id,
  data = sym_map,
  FUN = function(x) paste(unique(x), collapse = "|")
)
feature_map <- merge(feature_map, sym_agg, by = "gene_id", all.x = TRUE, sort = FALSE)
dup_symbol <- unique(feature_map$gene_symbol[duplicated(feature_map$gene_symbol)])
feature_map$duplicate_symbol <- feature_map$gene_symbol %in% dup_symbol
write.csv(feature_map, file.path(paths$metadata, "feature_map_scrna.csv"), row.names = FALSE)

n_ens <- sum(grepl("^ENSG", feature_map$gene_id))
pak2_log(
  "Feature map: ", nrow(feature_map), " unique gene_id values; Ensembl-like: ", n_ens,
  "; duplicated symbols: ", length(dup_symbol),
  .log_file = logf
)
if ("present_GSE285246" %in% names(feature_map) && "present_GSE173896" %in% names(feature_map)) {
  pak2_log(
    "Shared gene_id: ", sum(feature_map$present_GSE285246 & feature_map$present_GSE173896),
    "; GSE285246-only: ", sum(feature_map$present_GSE285246 & !feature_map$present_GSE173896),
    "; GSE173896-only: ", sum(!feature_map$present_GSE285246 & feature_map$present_GSE173896),
    .log_file = logf
  )
}

pak2_log("03_build_scrna_seurat.R complete.", .log_file = logf)
