#!/usr/bin/env Rscript
# Build canonical donor/sample metadata and search for author annotations.
# Does not infer missing covariates.

suppressPackageStartupMessages({
  library(readxl)
})

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
source(file.path(.script_dir, "utils.R"))

paths <- pak2_paths()
logf <- start_log(paths, "01_build_metadata")
acc <- readRDS(file.path(paths$metadata, "target_accessions.rds"))

soft_285 <- gunzip_copy(
  file.path(paths$raw_gse285246, "GSE285246_family.soft.gz"),
  file.path(paths$raw_gse285246, "GSE285246_family.soft")
)
soft_173 <- gunzip_copy(
  file.path(paths$raw_gse173896, "GSE173896_family.soft.gz"),
  file.path(paths$raw_gse173896, "GSE173896_family.soft")
)
geo_285 <- parse_soft_samples(soft_285)
geo_173 <- parse_soft_samples(soft_173)
geo <- rbind(geo_285, geo_173)
write.csv(geo, file.path(paths$metadata, "geo_soft_sample_records.csv"), row.names = FALSE)

lookup_geo <- function(gsm) {
  hit <- geo[geo$gsm_accession == gsm, , drop = FALSE]
  if (!nrow(hit)) stop("GSM not found in SOFT: ", gsm)
  hit[1, ]
}

# ---------------------------------------------------------------------------
# Paper Table S1 and other supplements
# ---------------------------------------------------------------------------

inspect_xlsx <- function(path, max_rows_preview = 8L) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size < 1) {
    return(NULL)
  }
  sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) character())
  out <- lapply(sheets, function(sh) {
    dat <- tryCatch(
      readxl::read_excel(path, sheet = sh, col_types = "text"),
      error = function(e) NULL
    )
    if (is.null(dat)) return(NULL)
    list(
      sheet = sh,
      ncol = ncol(dat),
      nrow = nrow(dat),
      columns = names(dat),
      data = dat
    )
  })
  names(out) <- sheets
  out[!vapply(out, is.null, logical(1))]
}

annotation_hits <- list()
record_hit <- function(source, sheet, kind, n, columns, note) {
  annotation_hits[[length(annotation_hits) + 1L]] <<- data.frame(
    source = source,
    sheet = sheet %||% NA_character_,
    kind = kind,
    n_rows = n,
    columns = paste(columns, collapse = " | "),
    note = note,
    stringsAsFactors = FALSE
  )
}

scan_table_for_annotations <- function(tbl, source, sheet) {
  if (is.null(tbl) || !nrow(tbl) || !ncol(tbl)) return(invisible(NULL))
  cols <- names(tbl)
  cols_l <- tolower(cols)
  barcode_cols <- cols[vapply(tbl, looks_like_barcode_column, logical(1))]
  if (!length(barcode_cols)) {
    barcode_cols <- cols[grepl("barcode|spot_id|cell_id|cellid", cols_l)]
  }
  ff_cols <- cols[grepl("(^|_)ff($|_)|fibroblastic|dense.?fibr|patholog|region_annot", cols_l)]
  ct_cols <- cols[grepl("cell.?type|celltype|annotation", cols_l)]
  if (length(barcode_cols) && length(ff_cols)) {
    vals <- unique(as.character(tbl[[ff_cols[[1]]]]))
    if (any(grepl("\\bFF\\b|\\bDF\\b|fibroblastic|dense fib", vals, ignore.case = TRUE))) {
      record_hit(source, sheet, "possible_spatial_region_annotation", nrow(tbl), union(barcode_cols, ff_cols),
                 "Barcode-like column together with FF/DF-like values")
    }
  }
  if (length(barcode_cols) && length(ct_cols)) {
    vals <- unique(as.character(tbl[[ct_cols[[1]]]]))
    if (any(grepl("fibroblast|macrophage|epithelial|endothelial|pericyte|myofibroblast|AT1|AT2|basaloid|smooth muscle", vals, ignore.case = TRUE))) {
      record_hit(source, sheet, "possible_celltype_annotation", nrow(tbl), union(barcode_cols, ct_cols),
                 "Barcode-like column together with cell-type-like values")
    }
  }
  invisible(NULL)
}

xlsx_targets <- c(
  file.path(paths$raw_gse285246, "paper", "ERJ-00022-2025.Tables.xlsx"),
  file.path(paths$raw_gse285246, "paper", "figshare_28138391_Supplymental_Table_4_6_7.xlsx"),
  file.path(paths$raw_gse173896, "GSE173896_COPD_meta_new.xlsx")
)

xlsx_contents <- list()
for (xp in xlsx_targets) {
  pak2_log("Inspecting ", xp, .log_file = logf)
  parsed <- inspect_xlsx(xp)
  xlsx_contents[[basename(xp)]] <- parsed
  if (is.null(parsed)) {
    pak2_log("  missing or unreadable", .log_file = logf)
    next
  }
  for (nm in names(parsed)) {
    info <- parsed[[nm]]
    pak2_log(
      "  sheet '", info$sheet, "': ", info$nrow, " rows x ", info$ncol,
      " cols: ", paste(info$columns, collapse = ", "),
      .log_file = logf
    )
    scan_table_for_annotations(info$data, basename(xp), info$sheet)
  }
}
saveRDS(xlsx_contents, file.path(paths$metadata, "supplement_xlsx_parsed.rds"))

# Parse Watanabe et al. Supplementary Table 1 (clinicopathological details).
# The Excel sheet has a title row, then a header row starting with "Donor".
paper_donor_extra <- data.frame(
  donor_id = character(),
  sex = character(),
  age = character(),
  smoking_status_paper = character(),
  pack_years = character(),
  bmi = character(),
  stringsAsFactors = FALSE
)

parse_table_s1 <- function(dat) {
  if (is.null(dat) || !nrow(dat)) return(NULL)
  first <- as.character(dat[[1]])
  header_i <- which(grepl("^Donor$", first, ignore.case = TRUE))[1]
  if (is.na(header_i)) return(NULL)
  hdr <- as.character(unlist(dat[header_i, ], use.names = FALSE))
  hdr[is.na(hdr)] <- paste0("col", seq_along(hdr))[is.na(hdr)]
  body <- dat[seq_len(nrow(dat)) > header_i, , drop = FALSE]
  names(body) <- hdr
  body <- body[!is.na(body$Donor) & nzchar(body$Donor), , drop = FALSE]
  if (!nrow(body)) return(NULL)
  findc <- function(pats) {
    hit <- names(body)[grepl(pats, names(body), ignore.case = TRUE)]
    if (!length(hit)) NA_character_ else hit[[1]]
  }
  data.frame(
    donor_id = as.character(body$Donor),
    sex = if (!is.na(findc("gend|sex"))) as.character(body[[findc("gend|sex")]]) else NA_character_,
    age = if (!is.na(findc("^age"))) as.character(body[[findc("^age")]]) else NA_character_,
    smoking_status_paper = if (!is.na(findc("smok"))) as.character(body[[findc("smok")]]) else NA_character_,
    pack_years = if (!is.na(findc("pack"))) as.character(body[[findc("pack")]]) else NA_character_,
    bmi = if (!is.na(findc("^bmi"))) as.character(body[[findc("^bmi")]]) else NA_character_,
    stringsAsFactors = FALSE
  )
}

tables_xlsx <- xlsx_contents[["ERJ-00022-2025.Tables.xlsx"]]
if (!is.null(tables_xlsx)) {
  s1_name <- names(tables_xlsx)[grepl("table 1|s1", names(tables_xlsx), ignore.case = TRUE)][1]
  if (!is.na(s1_name) && nzchar(s1_name)) {
    s1_raw <- tables_xlsx[[s1_name]]$data
    write.csv(s1_raw, file.path(paths$metadata, "paper_table_S1_raw.csv"), row.names = FALSE)
    parsed <- parse_table_s1(s1_raw)
    if (!is.null(parsed) && nrow(parsed)) {
      paper_donor_extra <- parsed
      paper_donor_extra$smoking_cessation <- paper_donor_extra$smoking_status_paper
      paper_donor_extra$smoking_status_paper <- ifelse(
        grepl("^never$", paper_donor_extra$smoking_cessation, ignore.case = TRUE),
        "never",
        ifelse(
          grepl("year", paper_donor_extra$smoking_cessation, ignore.case = TRUE),
          "former",
          paper_donor_extra$smoking_cessation
        )
      )
      write.csv(paper_donor_extra, file.path(paths$metadata, "paper_table_S1_parsed.csv"), row.names = FALSE)
      pak2_log("Parsed Table S1 donors: ", paste(paper_donor_extra$donor_id, collapse = ", "), .log_file = logf)
    } else {
      pak2_log("Table S1 present but donor rows could not be parsed.", .log_file = logf)
    }
  }
}

extra_val <- function(extra, col) {
  if (!nrow(extra) || !col %in% names(extra)) NA_character_ else extra[[col]][[1]]
}

spatial_rows <- lapply(seq_len(nrow(acc$gse285246_spatial)), function(i) {
  r <- acc$gse285246_spatial[i, ]
  g <- lookup_geo(r$gsm_accession)
  extra <- paper_donor_extra[paper_donor_extra$donor_id == r$donor_id, , drop = FALSE]
  data.frame(
    donor_id = r$donor_id,
    sample_id = r$section_id,
    section_id = r$section_id,
    condition = "IPF",
    disease = "IPF",
    modality = "spatial",
    source_GSE = "GSE285246",
    geo_accession = "GSE285246",
    gsm_accession = r$gsm_accession,
    tissue = "lung",
    tissue_group = "respiratory system",
    sex = extra_val(extra, "sex"),
    age = extra_val(extra, "age"),
    smoking_status = extra_val(extra, "smoking_status_paper"),
    smoking_cessation = extra_val(extra, "smoking_cessation"),
    pack_years = extra_val(extra, "pack_years"),
    bmi = extra_val(extra, "bmi"),
    library_chemistry = "10x Visium Spatial Gene Expression",
    genome_assembly = "GRCh38",
    sequencer = g$instrument_model,
    processing_pipeline = "Space Ranger",
    geo_source_name = g$source_name,
    geo_sample_title = g$sample_title,
    stringsAsFactors = FALSE
  )
})
spatial_metadata <- do.call(rbind, spatial_rows)

scrna_ipf_rows <- lapply(seq_len(nrow(acc$gse285246_scrna)), function(i) {
  r <- acc$gse285246_scrna[i, ]
  g <- lookup_geo(r$gsm_accession)
  extra <- paper_donor_extra[paper_donor_extra$donor_id == r$donor_id, , drop = FALSE]
  data.frame(
    donor_id = r$donor_id,
    sample_id = r$sample_id,
    section_id = NA_character_,
    condition = "IPF",
    disease = "IPF",
    modality = "scRNA",
    source_GSE = "GSE285246",
    geo_accession = "GSE285246",
    gsm_accession = r$gsm_accession,
    tissue = "lung",
    tissue_group = "respiratory system",
    sex = extra_val(extra, "sex"),
    age = extra_val(extra, "age"),
    smoking_status = extra_val(extra, "smoking_status_paper"),
    smoking_cessation = extra_val(extra, "smoking_cessation"),
    pack_years = extra_val(extra, "pack_years"),
    bmi = extra_val(extra, "bmi"),
    library_chemistry = "10x Chromium 3' v3.1",
    genome_assembly = "GRCh38",
    sequencer = g$instrument_model,
    processing_pipeline = "Cell Ranger 6.1.2",
    geo_source_name = g$source_name,
    geo_sample_title = g$sample_title,
    stringsAsFactors = FALSE
  )
})
scrna_ipf <- do.call(rbind, scrna_ipf_rows)

scrna_ctrl_rows <- lapply(seq_len(nrow(acc$gse173896_scrna)), function(i) {
  r <- acc$gse173896_scrna[i, ]
  g <- lookup_geo(r$gsm_accession)
  extra <- paper_donor_extra[paper_donor_extra$donor_id == r$donor_id, , drop = FALSE]
  data.frame(
    donor_id = r$donor_id,
    sample_id = r$sample_id,
    section_id = NA_character_,
    condition = "control",
    disease = "healthy",
    modality = "scRNA",
    source_GSE = "GSE173896",
    geo_accession = "GSE173896",
    gsm_accession = r$gsm_accession,
    tissue = "lung",
    tissue_group = "respiratory system",
    sex = extra_val(extra, "sex"),
    age = extra_val(extra, "age"),
    smoking_status = {
      s <- extra_val(extra, "smoking_status_paper")
      if (!is.na(s) && nzchar(s)) s else if (grepl("never", paste(g$source_name, g$disease_characteristic), ignore.case = TRUE)) "never" else NA_character_
    },
    smoking_cessation = extra_val(extra, "smoking_cessation"),
    pack_years = extra_val(extra, "pack_years"),
    bmi = extra_val(extra, "bmi"),
    library_chemistry = "10x Chromium 3' v3.1",
    genome_assembly = "GRCh38",
    sequencer = g$instrument_model,
    processing_pipeline = "Cell Ranger 3.0.2",
    geo_source_name = g$source_name,
    geo_sample_title = g$sample_title,
    stringsAsFactors = FALSE
  )
})
scrna_ctrl <- do.call(rbind, scrna_ctrl_rows)
scrna_metadata <- rbind(scrna_ipf, scrna_ctrl)

sample_manifest <- rbind(spatial_metadata, scrna_metadata)

donor_metadata <- unique(sample_manifest[, c(
  "donor_id", "condition", "disease", "source_GSE", "geo_accession",
  "tissue", "tissue_group", "sex", "age", "smoking_status", "smoking_cessation", "pack_years", "bmi"
)])

write.csv(sample_manifest, file.path(paths$metadata, "sample_manifest.csv"), row.names = FALSE)
write.csv(donor_metadata, file.path(paths$metadata, "donor_metadata.csv"), row.names = FALSE)
write.csv(spatial_metadata, file.path(paths$metadata, "spatial_metadata.csv"), row.names = FALSE)
write.csv(scrna_metadata, file.path(paths$metadata, "scrna_metadata.csv"), row.names = FALSE)

if (length(annotation_hits)) {
  hits_df <- do.call(rbind, annotation_hits)
} else {
  hits_df <- data.frame(
    source = character(), sheet = character(), kind = character(),
    n_rows = integer(), columns = character(), note = character(),
    stringsAsFactors = FALSE
  )
}
write.csv(hits_df, file.path(paths$metadata, "annotation_search_hits.csv"), row.names = FALSE)

ff_available <- any(hits_df$kind == "possible_spatial_region_annotation")
ct_available <- any(hits_df$kind == "possible_celltype_annotation")

availability <- list(
  ff_df_spot_annotations = ff_available,
  scrna_celltype_annotations = ct_available,
  age = any(!is.na(donor_metadata$age)),
  sex = any(!is.na(donor_metadata$sex)),
  smoking_status = any(!is.na(donor_metadata$smoking_status)),
  notes = c(
    "GEO GSE285246 processed files contain Space Ranger MTX/image outputs only; no spot-level FF/DF labels.",
    "GEO GSE285246 scRNA files contain MTX only; no cell-type labels.",
    "GSE173896_COPD_meta_new.xlsx has barcode-level QC/cluster columns and a Class field (COPD/NC/Never), not cell-type names.",
    "ERJ Supplementary Table 1 provides donor age, sex, smoking and pack-years; Tables 4-8 are gene/pathway lists, not barcode annotations.",
    "disease='healthy' for GSE173896 donors is taken from Watanabe et al. Eur Respir J 2025 ('three healthy controls'), not from GEO's 'Never smokers' characteristic.",
    "condition and source_GSE are stored as separate variables."
  )
)
saveRDS(availability, file.path(paths$metadata, "annotation_availability.rds"))

pak2_log("FF/DF spot annotations found: ", ff_available, .log_file = logf)
pak2_log("scRNA cell-type annotations found: ", ct_available, .log_file = logf)
pak2_log("Age available: ", availability$age, "; sex: ", availability$sex, "; smoking: ", availability$smoking_status, .log_file = logf)
pak2_log("01_build_metadata.R complete.", .log_file = logf)
