#!/usr/bin/env Rscript
# scRNA raw-count pseudobulk as SummarizedExperiment.
# Donor x cell type only if author annotations exist for the relevant cells.
# Otherwise write an explicitly unannotated donor-level object.

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
logf <- start_log(paths, "05_pseudobulk_scrna")
scrna_meta <- utils::read.csv(file.path(paths$metadata, "scrna_metadata.csv"), stringsAsFactors = FALSE)
feature_map <- utils::read.csv(file.path(paths$metadata, "feature_map_scrna.csv"), stringsAsFactors = FALSE)

obj_files <- file.path(paths$seurat_scrna, paste0(scrna_meta$donor_id, ".rds"))
missing <- obj_files[!file.exists(obj_files)]
if (length(missing)) stop("Missing scRNA Seurat objects:\n", paste(missing, collapse = "\n"))
objs <- lapply(obj_files, readRDS)
names(objs) <- scrna_meta$donor_id

annot_counts <- vapply(objs, function(o) {
  x <- o$cell_type_original
  sum(!is.na(x) & nzchar(x))
}, integer(1))
pak2_log(
  "Annotated cells per donor: ",
  paste(sprintf("%s=%s/%s", names(objs), annot_counts, vapply(objs, function(o) as.integer(ncol(o)), integer(1))), collapse = "; "),
  .log_file = logf
)

# Harmonised gene universe: Ensembl IDs; absent-from-matrix distinguished from zero.
universe <- feature_map$gene_id
row_df <- feature_map[, intersect(c("gene_id", "gene_symbol", "feature_type", "present_GSE285246", "present_GSE173896", "duplicate_symbol"), names(feature_map)), drop = FALSE]
rownames(row_df) <- universe

expand_to_universe <- function(vec_or_mat, gene_ids, universe) {
  if (is.null(dim(vec_or_mat))) {
    out <- matrix(NA_real_, nrow = length(universe), ncol = 1,
                  dimnames = list(universe, NULL))
    out[gene_ids, 1] <- as.numeric(vec_or_mat)
    methods::as(out, "dgCMatrix")
  } else {
    out <- Matrix::Matrix(
      NA_real_,
      nrow = length(universe),
      ncol = ncol(vec_or_mat),
      sparse = TRUE,
      dimnames = list(universe, colnames(vec_or_mat))
    )
    out[gene_ids, ] <- vec_or_mat
    methods::as(out, "dgCMatrix")
  }
}

# Always write donor-level unannotated aggregation as a data-integrity object.
donor_mat_list <- list()
donor_cd <- list()
for (nm in names(objs)) {
  obj <- objs[[nm]]
  raw <- get_raw_counts(obj, assay = "RNA")
  vec <- Matrix::rowSums(raw)
  if (sum(raw) != sum(vec)) stop("Donor rowSums mismatch for ", nm)
  donor_mat_list[[nm]] <- expand_to_universe(vec, rownames(raw), universe)
  md <- obj[[]]
  donor_cd[[nm]] <- data.frame(
    pseudobulk_id = nm,
    donor_id = unique(md$donor_id),
    sample_id = unique(md$sample_id),
    condition = unique(md$condition),
    disease = unique(md$disease),
    cell_type_original = NA_character_,
    cell_type_broad = NA_character_,
    n_cells = ncol(obj),
    library_size = as.numeric(sum(raw)),
    n_features_detected = as.integer(sum(vec > 0)),
    source_GSE = unique(md$source_GSE),
    geo_accession = unique(md$geo_accession),
    gsm_accession = unique(md$gsm_accession),
    tissue = unique(md$tissue),
    tissue_group = unique(md$tissue_group),
    sex = unique(md$sex),
    age = unique(md$age),
    smoking_status = unique(md$smoking_status),
    pack_years = unique(md$pack_years),
    bmi = unique(md$bmi),
    stringsAsFactors = FALSE
  )
}
donor_mat <- do.call(cbind, donor_mat_list)
colnames(donor_mat) <- names(donor_mat_list)
donor_cd <- do.call(rbind, donor_cd)
rownames(donor_cd) <- donor_cd$pseudobulk_id

# NA entries are genes absent from that donor's source matrix, not zeros.
se_donor <- make_summarized_experiment(
  counts = donor_mat,
  row_data = row_df,
  col_data = donor_cd,
  metadata = list(
    aggregation = "donor_id",
    suitable_for_posteriorHCA_celltype = FALSE,
    note = paste(
      "TEMPORARY data-integrity output: donor-level raw UMI sums.",
      "NOT suitable for posteriorHCA cell-type-specific comparison.",
      "Genes missing from a source feature matrix are stored as NA, not 0."
    )
  )
)
save_rds_atomic(se_donor, file.path(paths$pb_scrna, "scrna_by_donor_unannotated_SE.rds"))
pak2_log("Wrote scrna_by_donor_unannotated_SE.rds (NOT for cell-type posteriorHCA).", .log_file = logf)

ipf_ok <- all(annot_counts[scrna_meta$donor_id[scrna_meta$condition == "IPF"]] > 0)
ctrl_ok <- all(annot_counts[scrna_meta$donor_id[scrna_meta$condition == "control"]] > 0)

make_celltype_pb <- function(keep_donors, outfile, note) {
  mats <- list()
  cds <- list()
  leftover <- list()
  for (nm in keep_donors) {
    obj <- objs[[nm]]
    raw <- get_raw_counts(obj, assay = "RNA")
    lab <- as.character(obj$cell_type_original)
    keep <- !is.na(lab) & nzchar(lab)
    leftover[[nm]] <- data.frame(
      donor_id = nm,
      n_cells = ncol(obj),
      n_annotated = sum(keep),
      n_unannotated = sum(!keep),
      stringsAsFactors = FALSE
    )
    if (!any(keep)) next
    pb <- aggregate_counts_by_group(raw[, keep, drop = FALSE], lab[keep])
    if (sum(raw[, keep, drop = FALSE]) != sum(pb)) {
      stop("Cell-type aggregation count mismatch for ", nm)
    }
    ctypes <- colnames(pb)
    pb_ids <- paste(nm, ctypes, sep = "__")
    colnames(pb) <- pb_ids
    mats[[nm]] <- expand_to_universe(pb, rownames(raw), universe)
    md <- obj[[]]
    ntab <- as.integer(table(lab[keep])[ctypes])
    cds[[nm]] <- data.frame(
      pseudobulk_id = pb_ids,
      donor_id = nm,
      sample_id = unique(md$sample_id),
      condition = unique(md$condition),
      disease = unique(md$disease),
      cell_type_original = ctypes,
      cell_type_broad = vapply(ctypes, function(ct) {
        v <- unique(md$cell_type_broad[lab == ct & keep])
        v <- v[!is.na(v) & nzchar(v)]
        if (!length(v)) NA_character_ else paste(v, collapse = "|")
      }, character(1)),
      n_cells = ntab,
      library_size = as.numeric(Matrix::colSums(pb)),
      source_GSE = unique(md$source_GSE),
      geo_accession = unique(md$geo_accession),
      gsm_accession = unique(md$gsm_accession),
      tissue = unique(md$tissue),
      tissue_group = unique(md$tissue_group),
      sex = unique(md$sex),
      age = unique(md$age),
      smoking_status = unique(md$smoking_status),
    pack_years = unique(md$pack_years),
    bmi = unique(md$bmi),
      n_unannotated_cells_in_donor = sum(!keep),
      annotation_source = unique(na.omit(md$annotation_source[keep]))[1],
      stringsAsFactors = FALSE
    )
    if (sum(leftover[[nm]]$n_annotated) != sum(cds[[nm]]$n_cells)) {
      stop("n_cells across pseudobulks != annotated cells for ", nm)
    }
  }
  if (!length(mats)) {
    pak2_log("No annotated cells for requested donors; skipped ", outfile, .log_file = logf)
    return(invisible(NULL))
  }
  mat <- do.call(cbind, mats)
  cd <- do.call(rbind, cds)
  rownames(cd) <- cd$pseudobulk_id
  se <- make_summarized_experiment(
    counts = mat,
    row_data = row_df,
    col_data = cd,
    metadata = list(
      aggregation = "donor_id x cell_type_original",
      leftover_unannotated = do.call(rbind, leftover),
      note = note
    )
  )
  save_rds_atomic(se, outfile)
  utils::write.csv(do.call(rbind, leftover), paste0(outfile, ".unannotated_cells.csv"), row.names = FALSE)
  pak2_log("Wrote ", outfile, " with ", ncol(se), " pseudobulk columns.", .log_file = logf)
  invisible(se)
}

if (ipf_ok && ctrl_ok) {
  make_celltype_pb(
    keep_donors = names(objs),
    outfile = file.path(paths$pb_scrna, "scrna_by_donor_celltype_SE.rds"),
    note = paste(
      "Author-provided cell_type_original labels aggregated by donor x cell type.",
      "No minimum n_cells filter. Genes absent from a source matrix are NA, not 0."
    )
  )
} else {
  pak2_log(
    "Not writing scrna_by_donor_celltype_SE.rds: author cell-type labels are not available for all six donors (IPF complete=",
    ipf_ok, "; control complete=", ctrl_ok, ").",
    .log_file = logf
  )
  pak2_log(
    "Individual-cell Seurat objects are preserved for a later annotation stage.",
    .log_file = logf
  )
  # If one study has labels, keep them as a clearly named extra object.
  if (ctrl_ok && !ipf_ok) {
    make_celltype_pb(
      keep_donors = scrna_meta$donor_id[scrna_meta$condition == "control"],
      outfile = file.path(paths$pb_scrna, "scrna_GSE173896_control_by_donor_celltype_SE.rds"),
      note = paste(
        "CONTROL-ONLY extra object from GSE173896 source-study labels.",
        "This is NOT the Watanabe IPF paper fibroblast annotation and is NOT a combined IPF/control cell-type pseudobulk."
      )
    )
  }
  if (ipf_ok && !ctrl_ok) {
    make_celltype_pb(
      keep_donors = scrna_meta$donor_id[scrna_meta$condition == "IPF"],
      outfile = file.path(paths$pb_scrna, "scrna_GSE285246_IPF_by_donor_celltype_SE.rds"),
      note = "IPF-ONLY extra object. Control donors lack matching author cell-type labels."
    )
  }
}

# Per-donor integrity log
integrity <- do.call(rbind, lapply(names(objs), function(nm) {
  obj <- objs[[nm]]
  raw <- get_raw_counts(obj, assay = "RNA")
  data.frame(
    donor_id = nm,
    n_cells = ncol(obj),
    n_features = nrow(obj),
    total_UMI = as.numeric(sum(raw)),
    condition = unique(obj$condition),
    geo_accession = unique(obj$geo_accession),
    n_annotated = annot_counts[[nm]],
    n_unannotated = ncol(obj) - annot_counts[[nm]],
    donor_pseudobulk_UMI = as.numeric(sum(donor_mat[, nm], na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}))
# Donor PB UMI should equal source UMI (NAs are absent genes, not counts).
if (any(integrity$total_UMI != integrity$donor_pseudobulk_UMI)) {
  stop("Donor-level pseudobulk UMI does not match source totals")
}
utils::write.csv(integrity, file.path(paths$logs, "scrna_pseudobulk_validation.csv"), row.names = FALSE)
pak2_log("05_pseudobulk_scrna.R complete.", .log_file = logf)
