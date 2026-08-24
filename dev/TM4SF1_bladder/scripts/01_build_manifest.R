#!/usr/bin/env Rscript
# Construct the library/patient manifest and retain conflicting GEO/file mappings.

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
source(file.path(.script_dir, "utils.R"))

paths <- case_paths()
logf <- start_log(paths, "01_build_manifest")

soft <- parse_soft_samples(file.path(paths$raw, "GSE293189_family.soft.gz"))
if (nrow(soft) != 67L) stop("Expected 67 GEO samples, found ", nrow(soft))

soft$filename <- basename(soft$filename)
soft$library_id <- sub("_dge\\.txt\\.gz$", "", sub("^GSM[0-9]+_", "", soft$filename))
soft$geo_original_patient_id <- patient_from_title(soft$geo_title)
soft$original_patient_id <- patient_from_library(soft$library_id)
soft$title_filename_patient_match <- soft$geo_original_patient_id == soft$original_patient_id
soft$geo_library_filename_match <- soft$geo_library_id == soft$library_id
soft$mapping_status <- "consistent"
soft$mapping_status[!soft$title_filename_patient_match & soft$geo_library_filename_match] <-
  "GEO_title_vs_description_and_DGE_filename_conflict"
soft$mapping_status[soft$title_filename_patient_match & !soft$geo_library_filename_match] <-
  "GEO_description_vs_DGE_filename_conflict"
soft$mapping_status[!soft$title_filename_patient_match & !soft$geo_library_filename_match] <-
  "GEO_title_description_vs_DGE_filename_conflict"

geo_fields <- sample_fields_from_geo(soft$geo_title)
file_fields <- sample_fields_from_library(soft$library_id)
conflict <- !soft$geo_library_filename_match
soft$sample_type <- geo_fields$sample_type
soft$tumour_piece <- geo_fields$tumour_piece
soft$technical_rep <- geo_fields$technical_rep
soft$sample_type[conflict] <- file_fields$sample_type[conflict]
soft$tumour_piece[conflict] <- file_fields$tumour_piece[conflict]
soft$technical_rep[conflict] <- file_fields$technical_rep[conflict]

map <- paper_patient_map()
idx <- match(soft$original_patient_id, map$original_patient_id)
if (anyNA(idx)) stop("Unmapped filename-derived patient(s)")
for (nm in setdiff(names(map), "original_patient_id")) soft[[nm]] <- map[[nm]][idx]

soft$local_path <- file.path(paths$raw_dge, soft$filename)
if (any(!file.exists(soft$local_path))) stop("Manifest refers to missing DGE files")
soft$bytes <- file.info(soft$local_path)$size
soft$md5 <- unname(tools::md5sum(soft$local_path))

sra <- utils::read.csv(
  file.path(paths$raw, "PRJNA1243332_SraRunInfo.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
sra_cols <- c("Run", "Experiment", "LibraryName", "BioSample", "spots", "bases", "size_MB")
sra <- sra[, sra_cols]
names(sra)[names(sra) == "LibraryName"] <- "GSM"
names(sra)[names(sra) == "Experiment"] <- "experiment"
names(sra)[names(sra) == "BioSample"] <- "biosample"
soft <- merge(soft, sra, by = c("GSM", "experiment", "biosample"), all.x = TRUE, sort = FALSE)
if (anyNA(soft$Run)) stop("Could not attach SRA run metadata to every GSM")

keep <- c(
  "GSM", "library_id", "patient_id", "original_patient_id", "paper_id",
  "histology", "histology_group", "sample_type", "tumour_piece", "technical_rep",
  "filename", "geo_title", "geo_library_id", "geo_original_patient_id",
  "title_filename_patient_match", "geo_library_filename_match", "mapping_status",
  "tumour_cohort_status", "Run", "experiment", "biosample",
  "spots", "bases", "size_MB", "supplementary_url", "local_path", "bytes", "md5"
)
manifest <- soft[, keep]
manifest <- manifest[order(as.integer(sub("GSM", "", manifest$GSM))), ]

out <- file.path(paths$metadata, "GSE293189_sample_manifest.tsv")
write_tsv(manifest, out)

issues <- manifest[manifest$mapping_status != "consistent", c(
  "GSM", "geo_title", "geo_library_id", "filename", "original_patient_id",
  "geo_original_patient_id", "title_filename_patient_match",
  "geo_library_filename_match", "mapping_status"
)]
write_tsv(issues, file.path(paths$metadata, "metadata_mapping_inconsistencies.tsv"))

patient_manifest <- unique(manifest[, c(
  "patient_id", "original_patient_id", "paper_id", "histology", "histology_group",
  "tumour_cohort_status"
)])
patient_manifest$n_libraries <- as.integer(table(manifest$patient_id)[patient_manifest$patient_id])
patient_manifest$n_tumour_libraries <- vapply(
  patient_manifest$patient_id,
  function(z) sum(manifest$patient_id == z & manifest$sample_type == "tumour"),
  integer(1)
)
patient_manifest$n_paired_normal_libraries <- vapply(
  patient_manifest$patient_id,
  function(z) sum(manifest$patient_id == z & manifest$sample_type == "paired_normal"),
  integer(1)
)
write_tsv(patient_manifest, file.path(paths$metadata, "GSE293189_patient_manifest.tsv"))

case_log(
  "Wrote ", nrow(manifest), " libraries for ", length(unique(manifest$patient_id)),
  " patients/subjects; mapping conflicts: ", nrow(issues), file = logf
)
