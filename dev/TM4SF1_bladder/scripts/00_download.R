#!/usr/bin/env Rscript
# Download GEO metadata/DGE matrices, SRA run metadata, paper supplements,
# and the authors' public code. FASTQs are deliberately not downloaded.

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
source(file.path(.script_dir, "utils.R"))

paths <- case_paths()
logf <- start_log(paths, "00_download")

downloads <- data.frame(
  url = c(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE293nnn/GSE293189/soft/GSE293189_family.soft.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE293nnn/GSE293189/miniml/GSE293189_family.xml.tgz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE293nnn/GSE293189/suppl/filelist.txt",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE293nnn/GSE293189/suppl/GSE293189_RAW.tar",
    "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12174346/supplementaryFiles",
    "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12174346/fullTextXML"
  ),
  dest = file.path(paths$raw, c(
    "GSE293189_family.soft.gz", "GSE293189_family.xml.tgz",
    "GSE293189_filelist.txt", "GSE293189_RAW.tar",
    "PMC12174346_supplementaryFiles.zip", "PMC12174346_fullTextXML.xml"
  )),
  expected_bytes = c(8310, 8211, 5081, 92139520, NA, NA),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(downloads))) {
  download_url(
    downloads$url[[i]], downloads$dest[[i]],
    downloads$expected_bytes[[i]], log_file = logf
  )
}

if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required")
search_url <- paste0(
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?",
  "db=sra&term=PRJNA1243332%5BBioProject%5D&retmax=100&retmode=json"
)
search_path <- file.path(paths$raw, "PRJNA1243332_sra_esearch.json")
download_url(search_url, search_path, log_file = logf)
ids <- jsonlite::fromJSON(search_path)$esearchresult$idlist
if (length(ids) != 67L) stop("Expected 67 SRA records, found ", length(ids))
runinfo_url <- paste0(
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=",
  paste(ids, collapse = ","), "&rettype=runinfo&retmode=text"
)
runinfo_path <- file.path(paths$raw, "PRJNA1243332_SraRunInfo.csv")
download_url(runinfo_url, runinfo_path, log_file = logf)
if (nrow(utils::read.csv(runinfo_path, check.names = FALSE)) != 67L) {
  stop("SRA run-info table does not contain 67 records")
}

raw_tar <- file.path(paths$raw, "GSE293189_RAW.tar")
if (length(list.files(paths$raw_dge, pattern = "_dge\\.txt\\.gz$")) != 67L) {
  status <- system2("tar", c("-xf", raw_tar, "-C", paths$raw_dge))
  if (status != 0L) stop("Could not extract GEO raw archive")
}

supp_zip <- file.path(paths$raw, "PMC12174346_supplementaryFiles.zip")
if (!file.exists(file.path(paths$raw_paper, "41467_2025_59888_MOESM4_ESM.xlsx"))) {
  utils::unzip(supp_zip, exdir = paths$raw_paper, overwrite = TRUE)
}

author_url <- "https://github.com/angelussong/Histological_Variant_Bladder_Cancer_Analysis.git"
if (!dir.exists(file.path(paths$author_code, ".git"))) {
  status <- system2("git", c("clone", author_url, paths$author_code))
  if (status != 0L) stop("Could not clone authors' code")
}
author_commit <- system2(
  "git", c("-C", paths$author_code, "rev-parse", "HEAD"),
  stdout = TRUE
)

local_files <- downloads
local_files$bytes <- file.info(local_files$dest)$size
local_files$md5 <- unname(tools::md5sum(local_files$dest))
write_tsv(local_files, file.path(paths$metadata, "download_manifest.tsv"))
writeLines(author_commit, file.path(paths$metadata, "authors_code_commit.txt"))

case_log("Downloaded/extracted 67 DGE matrices; authors' commit ", author_commit, file = logf)

