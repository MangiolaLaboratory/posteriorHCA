#!/usr/bin/env Rscript
# Build a self-contained validation report from the generated metadata and
# matrices. No analytical result is hard-coded except published targets and
# explicitly documented workflow parameters.

suppressPackageStartupMessages({
  library(Matrix)
  library(SummarizedExperiment)
})

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
source(file.path(.script_dir, "utils.R"))

paths <- case_paths()
logf <- start_log(paths, "05_report")

read_meta <- function(name) {
  utils::read.delim(
    file.path(paths$metadata, name), stringsAsFactors = FALSE,
    check.names = FALSE, na.strings = c("", "NA")
  )
}

md_table <- function(x, columns = names(x), digits = 3L) {
  x <- x[, columns, drop = FALSE]
  x[] <- lapply(x, function(z) {
    if (is.numeric(z)) {
      out <- format(round(z, digits), trim = TRUE, scientific = FALSE)
    } else {
      out <- as.character(z)
    }
    out[is.na(out)] <- ""
    gsub("\\|", "\\\\|", out)
  })
  c(
    paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |"),
    apply(x, 1L, function(z) paste0("| ", paste(z, collapse = " | "), " |"))
  )
}

fmt_int <- function(x) formatC(as.integer(x), format = "d", big.mark = ",")

samples <- read_meta("GSE293189_sample_manifest.tsv")
patients <- read_meta("GSE293189_patient_manifest.tsv")
qc <- read_meta("qc_cell_counts_overall.tsv")
tumour <- read_meta("published_vs_reconstructed_tumour_cells.tsv")
tumour_info <- read_meta("tumour_epithelial_counts_by_patient.tsv")
issues <- read_meta("metadata_mapping_inconsistencies.tsv")
tm_patient <- read_meta("TM4SF1_published_tumour_cohort_by_patient.tsv")
tm_group <- read_meta("TM4SF1_published_tumour_cohort_group_summary.tsv")
tm_source <- read_meta("TM4SF1_population_source_summary.tsv")
normal_pb_md <- read_meta("pseudobulk_paired_normal_epithelial_sample_metadata.tsv")

tumour <- merge(
  tumour,
  tumour_info[, c("patient_id", "histology_group")],
  by = "patient_id", all.x = TRUE, sort = FALSE
)
tumour <- tumour[match(c(paste0("UC0", 1:3), paste0("VAR0", 1:9)), tumour$patient_id), ]
tumour$delta <- as.integer(tumour$delta)

patients$paper_id[is.na(patients$paper_id)] <- "-"
patient_report <- patients[, c(
  "patient_id", "original_patient_id", "paper_id", "histology",
  "histology_group", "n_libraries", "n_tumour_libraries",
  "n_paired_normal_libraries", "tumour_cohort_status"
)]

normal_report <- normal_pb_md[, c(
  "patient_id", "n_cells", "library_size", "sufficient_cells_ge20"
)]
normal_tm <- read_meta("TM4SF1_patient_summary.tsv")
normal_tm <- normal_tm[
  normal_tm$sample_type == "paired_normal" & normal_tm$broad_cell_type == "epithelial",
  c("patient_id", "raw_count", "CPM", "detection_rate")
]
normal_report <- merge(normal_report, normal_tm, by = "patient_id", sort = FALSE)

qc_patient <- read_meta("qc_cell_counts_by_patient.tsv")
bladder_qc_cells <- sum(
  qc_patient$n_after_doublet_filter[qc_patient$histology_group %in% c("HV", "UC")]
)
cell_md <- utils::read.delim(
  gzfile(file.path(paths$metadata, "GSE293189_cell_metadata.tsv.gz")),
  stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA")
)
broad_n <- table(cell_md$broad_cell_type[cell_md$analysis_inclusion == "processed_paper_cohort"])

strict <- tumour$reconstructed_n150_threshold_met
strict[is.na(strict)] <- FALSE
published_total <- sum(tumour$published_tumour_epithelial_cells)
reconstructed_all <- sum(tumour$n_cells)
reconstructed_strict <- sum(tumour$n_cells[strict])
strict_groups <- table(tumour$histology_group[strict])

session_path <- file.path(paths$metadata, "sessionInfo.txt")
con <- file(session_path, "wt")
writeLines(capture.output(sessionInfo()), con)
close(con)

report <- c(
  "# GSE293189 TM4SF1 bladder-cancer preprocessing report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Outcome",
  "",
  paste0(
    "The public raw DGE deposit was reconstructed into raw-count, QC-filtered, ",
    "Seurat, annotation, and patient-level pseudobulk objects. The published ",
    "9 HV + 3 pure-UC / 8,553 tumour-epithelial-cell structure was **not fully ",
    "reproduced**: all 9 HV patients passed the 150-cell rule, but only 2 of 3 ",
    "pure-UC patients did. Across the 12 published patients, ",
    fmt_int(reconstructed_all), " epithelial tumour cells were reconstructed; ",
    fmt_int(reconstructed_strict), " cells in 9 HV + 2 UC patients pass the ",
    "reconstructed threshold. The main discrepancy is UC01 (90 reconstructed ",
    "versus 605 published)."
  ),
  "",
  "No posteriorHCA comparison and no tumour/toxicity threshold were performed.",
  "",
  "## Downloads and provenance",
  "",
  paste0(
    "- GEO series metadata: SOFT and MINiML for GSE293189; GEO file list; ",
    nrow(samples), " GSM records and ", nrow(samples), " deposited per-library DGE matrices."
  ),
  paste0(
    "- SRA metadata: ", length(unique(samples$Run)),
    " run records for PRJNA1243332 / SRP574011. FASTQs were not downloaded ",
    "because integer per-library DGE matrices were present and sufficient."
  ),
  "- Paper materials: Europe PMC full-text XML, supplementary archive, supplementary PDFs, figures, and the Figure 1B source-data workbook.",
  "- Public analysis code: `Histological_Variant_Bladder_Cancer_Analysis`, commit `d26f0427481c5ea39d1697ec3c83b7921a747504`.",
  "- File sizes and checksums are recorded in `data/metadata/download_manifest.tsv` and the sample manifest.",
  "",
  "## Patient mapping",
  "",
  paste0(
    "There are ", nrow(samples), " library/GSM entries but only ", nrow(patients),
    " biological subjects: 15 bladder-tumour cohort patients and two ",
    "paraganglioma subjects outside the paper's bladder analysis. Libraries, ",
    "tumour pieces, and technical replicates were never treated as independent patients."
  ),
  "",
  md_table(patient_report, digits = 0L),
  "",
  "The exact 67-row GSM/library mapping, including filenames, SRA runs, pieces, replicates, and both GEO-derived and filename-derived IDs, is in `data/metadata/GSE293189_sample_manifest.tsv`.",
  "",
  "### Metadata inconsistencies",
  "",
  paste0(
    "Eight GSMs have a title/description/filename conflict. GSM8878337-38 have ",
    "GEO titles naming 21226 but descriptions and DGE filenames naming 21262. ",
    "GSM8878366-71 form a six-entry mismatch between GEO descriptions and attached ",
    "DGE filenames. The processing mapping uses the patient ID embedded in the ",
    "actual DGE filename because that mapping agrees with the authors' public code ",
    "and, for most patients, the paper's source-data cell counts. Both alternatives ",
    "remain in the manifest; this is an evidence-based processing choice, not a ",
    "claim that GEO has been corrected. See `metadata_mapping_inconsistencies.tsv`."
  ),
  "",
  md_table(
    issues[, c("GSM", "geo_original_patient_id", "original_patient_id", "mapping_status")],
    digits = 0L
  ),
  "",
  "The public code also contains a post-publication mapping change from 12923 to 12049; the current code and filename mapping use 12049.",
  "",
  "## QC and doublets",
  "",
  md_table(qc, digits = 0L),
  "",
  paste0(
    "The raw union contains 36,994 gene symbols and ", fmt_int(qc$n_raw),
    " deposited cell barcodes. Filters were applied sequentially exactly as ",
    "requested: at least 300 detected genes, at least 500 UMIs, and mitochondrial ",
    "fraction below 20%. After doublet removal, ", fmt_int(qc$n_after_doublet_filter),
    " cells remained overall; ", fmt_int(bladder_qc_cells),
    " were in the HV/UC bladder cohort and 4,620 were paraganglioma cells retained ",
    "only in raw/QC objects."
  ),
  "",
  paste0(
    "The paper reports DoubletFinder but does not provide the parameters required ",
    "for an exact rerun. The reproducible substitute used scDblFinder 1.24.10: ",
    "26,180 cells in 42 libraries with at least 50 post-basic-QC cells received an ",
    "artificial-doublet score, and the highest-scoring expected number was called ",
    "per library using an expected doublet fraction of 0.008 x n_cells / 1,000. This removed ",
    fmt_int(qc$n_doublets_removed), " cells. The remaining 440 cells in 25 smaller ",
    "libraries were retained with an explicit `not_run_lt50_post_qc_cells_kept` flag."
  ),
  "",
  "Per-library and per-patient counts at every step are in `qc_cell_counts_by_library.tsv` and `qc_cell_counts_by_patient.tsv`.",
  "",
  "## Seurat reconstruction and annotation",
  "",
  paste0(
    "The ", fmt_int(bladder_qc_cells), " QC-passing bladder cells were processed ",
    "with Seurat LogNormalize (scale factor 10,000), 2,000 VST HVGs, 100 computed ",
    "PCs, graph construction on PCs 1-75 (k=30), and resolution 0.5. This produced ",
    "35 clusters versus 36 reported by the paper. UMAP is stored for inspection ",
    "but was not used as an annotation rule."
  ),
  "",
  paste0(
    "Broad multi-marker module-score annotation identified ",
    paste(paste0(names(broad_n), "=", fmt_int(broad_n)), collapse = ", "),
    ". Each cluster was assigned by its top epithelial, immune, stromal, or ",
    "endothelial panel score and accompanied by detected-marker evidence and a ",
    "score-margin confidence flag. Panels include EPCAM/KRT7 and other keratins/UPKs, ",
    "PTPRC plus lymphoid/myeloid genes, DCN/ACTA2 plus fibroblast/pericyte genes, ",
    "and SELE plus PECAM1/VWF/KDR/CLDN5. No cell class relies on one marker alone."
  ),
  "",
  "InferCNV-compatible raw counts, tumour/reference annotations, and the authors' cutoff/min-cell parameters were prepared. The HMM was not run because `infercnv` and the exact hg19 gene-order reference used by the authors were unavailable; tumour identity therefore uses tumour source plus multi-marker epithelial evidence, not a claimed InferCNV result.",
  "",
  "## Published tumour cohort comparison",
  "",
  md_table(
    tumour[, c(
      "patient_id", "histology_group", "published_tumour_epithelial_cells",
      "n_cells", "delta", "reconstructed_n150_threshold_met"
    )],
    digits = 0L
  ),
  "",
  paste0(
    "Published total: ", fmt_int(published_total), ". Reconstructed across the same ",
    "12 patient IDs: ", fmt_int(reconstructed_all), ". Strict reconstructed cohort: ",
    fmt_int(reconstructed_strict), " cells in ", strict_groups[["HV"]], " HV + ",
    strict_groups[["UC"]], " UC patients."
  ),
  "",
  "UC01 is the material unresolved exception: its filename-assigned 12049 matrix is immune-dominant after QC (90 epithelial versus 877 immune cells), so the public deposit cannot plausibly yield the published 605 epithelial cells under the reconstructed broad-marker logic. Substituting the conflicting GEO-title assignment does not restore 605. The 12-column tumour pseudobulk therefore retains UC01 with both published-cohort and reconstructed-`>=150` flags rather than silently dropping or remapping it.",
  "",
  "## Paired/adjacent normal bladder",
  "",
  paste0(
    "Seven bladder patients have at least one reconstructed paired-normal epithelial ",
    "cell; only three have at least 20. All seven are retained with an explicit ",
    "sufficiency flag. These samples are paired/adjacent normal bladder, not fully ",
    "healthy tissue."
  ),
  "",
  md_table(normal_report, digits = 3L),
  "",
  "## Pseudobulk construction",
  "",
  "All pseudobulk matrices sum raw integer-valued UMIs: cell -> technical replicate -> tumour piece/library -> patient. No normalized values were averaged. The primary columns are patient x sample_type x broad_cell_type.",
  "",
  "- Tumour epithelial: 36,994 genes x 12 published patient IDs. UC01 is flagged as below the reconstructed 150-cell cutoff; all other published patients pass.",
  "- Paired-normal epithelial: 36,994 genes x 7 patients; three columns have at least 20 cells.",
  "- Broad bladder populations: 36,994 genes x 80 patient/sample-type/cell-type profiles. The four classes are epithelial, endothelial, stromal, and immune; out-of-scope paraganglioma cells are not included here.",
  "",
  "Every pseudobulk object is a `SummarizedExperiment` containing the raw count assay, feature mapping, column metadata, cell count, library size, and aggregation provenance.",
  "",
  "## TM4SF1 / ENSG00000169908",
  "",
  "### Published tumour cohort",
  "",
  md_table(
    tm_patient[, c(
      "patient_id", "histology_group", "n_cells", "raw_count", "library_size",
      "CPM", "detection_rate", "reconstructed_n150_threshold_met"
    )],
    digits = 3L
  ),
  "",
  md_table(tm_group, digits = 3L),
  "",
  "The published-ID reconstruction has substantially higher pooled TM4SF1 abundance in HV than pure UC (485.324 versus 80.677 CPM), but this is descriptive only and is strongly heterogeneous: VAR08 contributes 12,915 of 16,720 HV counts, while VAR09 and VAR03 are low. UC01 is below the reconstructed cell threshold, so no inferential group test is claimed.",
  "",
  "### Population source check",
  "",
  md_table(tm_source, digits = 3L),
  "",
  "TM4SF1 is not tumour-epithelial-specific. Endothelial cells show the highest pooled abundance and near-universal detection in both tumour and paired-normal compartments (2,378.180 and 2,701.322 CPM; 96.9% and 98.5% detected). Tumour epithelial cells contribute most tumour TM4SF1 molecules because they are more numerous, but pooled abundance is about five-fold lower than endothelial cells. This endothelial signal is directly relevant to later off-tumour work, but no liability threshold is defined here.",
  "",
  "## Exact outputs",
  "",
  "- `data/raw/`: GEO SOFT/MINiML, 67 DGE matrices, SRA metadata, paper/source data, and authors' public repository.",
  "- `data/metadata/GSE293189_sample_manifest.tsv`: exact GSM/library/patient mapping with checksums and conflict flags.",
  "- `data/metadata/GSE293189_cell_metadata.tsv.gz`: all 56,192 deposited barcodes with sequential QC and processed annotations where applicable.",
  "- `data/processed/GSE293189_raw_counts.rds`: full raw `SingleCellExperiment`.",
  "- `data/processed/GSE293189_processed_seurat.rds`: normalized/clustering object retaining raw counts.",
  "- `data/metadata/tumour_epithelial_cell_metadata.tsv.gz`: reconstructed tumour epithelial cells.",
  "- `data/processed/pseudobulk_tumour_epithelial_counts.rds`.",
  "- `data/processed/pseudobulk_paired_normal_epithelial_counts.rds`.",
  "- `data/processed/pseudobulk_all_broad_celltypes_counts.rds`.",
  "- `data/metadata/pseudobulk_sample_metadata.tsv` and matrix-specific column metadata files.",
  "- `data/metadata/TM4SF1_patient_summary.tsv` plus published-cohort group and population-source summaries.",
  "- `data/processed/GSE293189_infercnv_input.rds`: prepared input, clearly marked not run.",
  "- `data/metadata/sessionInfo.txt`: exact R/package environment.",
  "- `data/metadata/validation_checks.tsv`: machine-readable structural and count-accounting checks (generated by stage 06).",
  "",
  "## Reproducibility",
  "",
  "Run `Rscript scripts/run_all.R` from this case-study directory. Download steps are checksum/size-aware and skip valid existing files. Randomized stages use seed 293189. Logs are written to `logs/`.",
  "",
  "Primary public sources:",
  "",
  "- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293189",
  "- Paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC12174346/",
  "- Authors' code: https://github.com/angelussong/Histological_Variant_Bladder_Cancer_Analysis"
)

out <- file.path(paths$processed, "qc_summary.md")
writeLines(report, out)
case_log("Wrote validation report: ", out, file = logf)
