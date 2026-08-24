#!/usr/bin/env Rscript
# Download processed GEO files for GSE285246 (all samples) and GSE173896
# (JK06, JK11, JK12 only), plus paper supplementary tables from PMC.

suppressPackageStartupMessages({
  library(utils)
})

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
source(file.path(.script_dir, "utils.R"))

paths <- pak2_paths()
logf <- start_log(paths, "00_download_geo")

# ---------------------------------------------------------------------------
# Accessions we keep
# ---------------------------------------------------------------------------

gse285246_spatial <- data.frame(
  gsm_accession = paste0("GSM8699", 165:173),
  section_id = c(
    "JKPF1_1", "JKPF1_2",
    "JKPF2_1", "JKPF2_2", "JKPF2_3",
    "JKPF3_1", "JKPF3_2", "JKPF3_3", "JKPF3_5"
  ),
  donor_id = c("JKPF1", "JKPF1", "JKPF2", "JKPF2", "JKPF2", "JKPF3", "JKPF3", "JKPF3", "JKPF3"),
  stringsAsFactors = FALSE
)

gse285246_scrna <- data.frame(
  gsm_accession = c("GSM8699174", "GSM8699175", "GSM8699176"),
  donor_id = c("JKPF1", "JKPF2", "JKPF3"),
  sample_id = c("JKPF1", "JKPF2", "JKPF3"),
  file_stub = c("IPF_JKPF1", "IPF_JKPF2", "IPF_JKPF3"),
  stringsAsFactors = FALSE
)

gse173896_scrna <- data.frame(
  gsm_accession = c("GSM5282546", "GSM5282547", "GSM5282548"),
  donor_id = c("JK06", "JK11", "JK12"),
  sample_id = c("JK06", "JK11", "JK12"),
  file_stub = c("JK06", "JK11", "JK12"),
  stringsAsFactors = FALSE
)

spatial_suffixes <- c(
  barcodes = "barcodes.tsv.gz",
  features = "features.tsv.gz",
  matrix = "matrix.mtx.gz",
  positions = "tissue_positions_list.csv.gz",
  scalefactors = "scalefactors_json.json.gz",
  hires = "tissue_hires_image.png.gz",
  lowres = "tissue_lowres_image.png.gz",
  aligned_fiducials = "aligned_fiducials.jpg.gz",
  detected_tissue = "detected_tissue_image.jpg.gz"
)

scrna_suffixes <- c(
  barcodes = "barcodes.tsv.gz",
  features = "features.tsv.gz",
  matrix = "matrix.mtx.gz"
)

gsm_ftp_https <- function(gsm, file_name) {
  prefix <- paste0(substr(gsm, 1, 7), "nnn")
  sprintf(
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/%s/%s/suppl/%s",
    prefix, gsm, file_name
  )
}

# ---------------------------------------------------------------------------
# Series-level metadata files
# ---------------------------------------------------------------------------

series_files <- list(
  list(
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285246/soft/GSE285246_family.soft.gz",
    dest = file.path(paths$raw_gse285246, "GSE285246_family.soft.gz")
  ),
  list(
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285246/suppl/filelist.txt",
    dest = file.path(paths$raw_gse285246, "filelist.txt")
  ),
  list(
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173896/soft/GSE173896_family.soft.gz",
    dest = file.path(paths$raw_gse173896, "GSE173896_family.soft.gz")
  ),
  list(
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173896/suppl/filelist.txt",
    dest = file.path(paths$raw_gse173896, "filelist.txt")
  ),
  list(
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173896/suppl/GSE173896_COPD_meta_new.xlsx",
    dest = file.path(paths$raw_gse173896, "GSE173896_COPD_meta_new.xlsx")
  ),
  list(
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12501429/supplementaryFiles",
    dest = file.path(paths$raw_gse285246, "paper", "PMC12501429_supplementaryFiles.zip")
  ),
  list(
    url = "https://ndownloader.figshare.com/files/51484487",
    dest = file.path(paths$raw_gse285246, "paper", "figshare_28138391_Supplymental_Table_4_6_7.xlsx")
  )
)

for (sf in series_files) {
  download_url(sf$url, sf$dest, log_file = logf)
}

# PMC HTML interstitial blocks direct wget of bin/ files; EuropePMC zip is the reliable source.
epmc_zip <- file.path(paths$raw_gse285246, "paper", "PMC12501429_supplementaryFiles.zip")
if (file.exists(epmc_zip) && file.info(epmc_zip)$size > 1e5) {
  unzip(
    epmc_zip,
    files = c("ERJ-00022-2025.Tables.xlsx", "ERJ-00022-2025.Methods.pdf"),
    exdir = file.path(paths$raw_gse285246, "paper"),
    overwrite = TRUE
  )
  pak2_log("Extracted ERJ supplementary tables/methods from EuropePMC zip.", .log_file = logf)
}

fl_285 <- parse_geo_filelist(file.path(paths$raw_gse285246, "filelist.txt"))
fl_173 <- parse_geo_filelist(file.path(paths$raw_gse173896, "filelist.txt"))
expected_size <- function(file_name, fl) {
  hit <- fl$expected_size[fl$file_name == file_name]
  if (!length(hit)) NA_real_ else as.numeric(hit[[1]])
}

manifest_rows <- list()
add_manifest <- function(...) {
  manifest_rows[[length(manifest_rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# GSE285246 spatial
# ---------------------------------------------------------------------------

for (i in seq_len(nrow(gse285246_spatial))) {
  gsm <- gse285246_spatial$gsm_accession[[i]]
  section <- gse285246_spatial$section_id[[i]]
  donor <- gse285246_spatial$donor_id[[i]]
  out_dir <- file.path(paths$raw_gse285246, "spatial", section)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (kind in names(spatial_suffixes)) {
    fname <- sprintf("%s_%s_%s", gsm, section, spatial_suffixes[[kind]])
    dest <- file.path(out_dir, fname)
    url <- gsm_ftp_https(gsm, fname)
    download_url(url, dest, expected_size = expected_size(fname, fl_285), log_file = logf)
    add_manifest(
      geo_accession = "GSE285246",
      gsm_accession = gsm,
      modality = "spatial",
      condition = "IPF",
      donor = donor,
      sample_id = section,
      section_id = section,
      file_role = kind,
      local_path = dest,
      source_url_or_geo_reference = url,
      expected_bytes = expected_size(fname, fl_285),
      downloaded_bytes = file.info(dest)$size
    )
  }
}

# ---------------------------------------------------------------------------
# GSE285246 scRNA
# ---------------------------------------------------------------------------

for (i in seq_len(nrow(gse285246_scrna))) {
  gsm <- gse285246_scrna$gsm_accession[[i]]
  donor <- gse285246_scrna$donor_id[[i]]
  stub <- gse285246_scrna$file_stub[[i]]
  out_dir <- file.path(paths$raw_gse285246, "scrna", donor)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (kind in names(scrna_suffixes)) {
    fname <- sprintf("%s_%s_%s", gsm, stub, scrna_suffixes[[kind]])
    dest <- file.path(out_dir, fname)
    url <- gsm_ftp_https(gsm, fname)
    download_url(url, dest, expected_size = expected_size(fname, fl_285), log_file = logf)
    add_manifest(
      geo_accession = "GSE285246",
      gsm_accession = gsm,
      modality = "scRNA",
      condition = "IPF",
      donor = donor,
      sample_id = donor,
      section_id = NA_character_,
      file_role = kind,
      local_path = dest,
      source_url_or_geo_reference = url,
      expected_bytes = expected_size(fname, fl_285),
      downloaded_bytes = file.info(dest)$size
    )
  }
}

# ---------------------------------------------------------------------------
# GSE173896 control scRNA only (JK06, JK11, JK12)
# ---------------------------------------------------------------------------

pak2_log("Skipping remaining GSE173896 donors; keeping JK06, JK11, JK12 only.", .log_file = logf)
pak2_log("Not downloading GSE173896_COPD.rds.gz (11.9 Gb) or GSE173896_RAW.tar.", .log_file = logf)

for (i in seq_len(nrow(gse173896_scrna))) {
  gsm <- gse173896_scrna$gsm_accession[[i]]
  donor <- gse173896_scrna$donor_id[[i]]
  stub <- gse173896_scrna$file_stub[[i]]
  out_dir <- file.path(paths$raw_gse173896, "scrna", donor)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (kind in names(scrna_suffixes)) {
    fname <- sprintf("%s_%s_%s", gsm, stub, scrna_suffixes[[kind]])
    dest <- file.path(out_dir, fname)
    url <- gsm_ftp_https(gsm, fname)
    download_url(url, dest, expected_size = expected_size(fname, fl_173), log_file = logf)
    add_manifest(
      geo_accession = "GSE173896",
      gsm_accession = gsm,
      modality = "scRNA",
      condition = "control",
      donor = donor,
      sample_id = donor,
      section_id = NA_character_,
      file_role = kind,
      local_path = dest,
      source_url_or_geo_reference = url,
      expected_bytes = expected_size(fname, fl_173),
      downloaded_bytes = file.info(dest)$size
    )
  }
}

# Series-level extras in the manifest
add_manifest(
  geo_accession = "GSE173896",
  gsm_accession = NA_character_,
  modality = "scRNA",
  condition = "control",
  donor = NA_character_,
  sample_id = NA_character_,
  section_id = NA_character_,
  file_role = "copd_study_cell_metadata_xlsx",
  local_path = file.path(paths$raw_gse173896, "GSE173896_COPD_meta_new.xlsx"),
  source_url_or_geo_reference = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173896/suppl/GSE173896_COPD_meta_new.xlsx",
  expected_bytes = NA_real_,
  downloaded_bytes = file.info(file.path(paths$raw_gse173896, "GSE173896_COPD_meta_new.xlsx"))$size
)
add_manifest(
  geo_accession = "GSE285246",
  gsm_accession = NA_character_,
  modality = "paper_supplement",
  condition = NA_character_,
  donor = NA_character_,
  sample_id = NA_character_,
  section_id = NA_character_,
  file_role = "erj_supplementary_tables_xlsx",
  local_path = file.path(paths$raw_gse285246, "paper", "ERJ-00022-2025.Tables.xlsx"),
  source_url_or_geo_reference = "https://pmc.ncbi.nlm.nih.gov/articles/instance/12501429/bin/ERJ-00022-2025.Tables.xlsx",
  expected_bytes = NA_real_,
  downloaded_bytes = file.info(file.path(paths$raw_gse285246, "paper", "ERJ-00022-2025.Tables.xlsx"))$size
)

manifest <- do.call(rbind, manifest_rows)
missing <- manifest[!file.exists(manifest$local_path) | is.na(manifest$downloaded_bytes) | manifest$downloaded_bytes < 1, ]
if (nrow(missing)) {
  stop("Missing essential files:\n", paste(missing$local_path, collapse = "\n"))
}

size_fail <- which(
  !is.na(manifest$expected_bytes) &
    manifest$expected_bytes > 0 &
    manifest$downloaded_bytes != manifest$expected_bytes
)
if (length(size_fail)) {
  stop("Size verification failed for:\n", paste(manifest$local_path[size_fail], collapse = "\n"))
}

write.csv(manifest, file.path(paths$metadata, "download_manifest.csv"), row.names = FALSE)
saveRDS(
  list(
    gse285246_spatial = gse285246_spatial,
    gse285246_scrna = gse285246_scrna,
    gse173896_scrna = gse173896_scrna
  ),
  file.path(paths$metadata, "target_accessions.rds")
)

pak2_log("Downloaded ", nrow(manifest), " files. Manifest: data/metadata/download_manifest.csv", .log_file = logf)
pak2_log("00_download_geo.R complete.", .log_file = logf)
