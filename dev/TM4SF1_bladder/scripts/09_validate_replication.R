#!/usr/bin/env Rscript
# Publication-replication checks layered on top of the canonical structural
# validation. Known unresolved publication targets are reported as failures but
# do not prevent the audit report from being written.

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(Matrix)
})

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
source(file.path(.script_dir, "utils.R"))

paths <- case_paths()
audit <- file.path(paths$root, "results", "replication_audit")
dir.create(audit, recursive = TRUE, showWarnings = FALSE)
logf <- start_log(paths, "09_validate_replication")

flow <- utils::read.delim(
  file.path(audit, "all_patient_cell_flow.tsv"), stringsAsFactors = FALSE,
  check.names = FALSE
)
barcode <- utils::read.delim(
  file.path(audit, "uc01_published_barcode_audit.tsv"), stringsAsFactors = FALSE,
  check.names = FALSE
)
structural <- utils::read.delim(
  file.path(paths$metadata, "validation_checks.tsv"), stringsAsFactors = FALSE,
  check.names = FALSE
)
tumour_pb <- readRDS(file.path(paths$processed, "pseudobulk_tumour_epithelial_counts.rds"))

expected <- c(paste0("VAR0", 1:9), paste0("UC0", 1:3))
retained <- flow$patient_id[flow$n_final_tumour_epithelial >= 150L]

checks <- list()
add <- function(check, category, observed, expected_value, pass, evidence) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = check,
    category = category,
    observed = paste(observed, collapse = ","),
    expected = paste(expected_value, collapse = ","),
    pass = isTRUE(pass),
    evidence = evidence,
    stringsAsFactors = FALSE
  )
}

add(
  "expected_patients_present", "cohort_logic",
  sort(intersect(flow$patient_id, expected)), sort(expected),
  all(expected %in% flow$patient_id),
  "all_patient_cell_flow.tsv contains one row for each expected published ID"
)
add(
  "VAR01_to_VAR09_retained", "cohort_logic",
  sort(intersect(retained, paste0("VAR0", 1:9))), paste0("VAR0", 1:9),
  all(paste0("VAR0", 1:9) %in% retained),
  "retention recomputed from final tumour-epithelial counts >=150"
)
for (pid in c("VAR10", "VAR11", "UC04")) {
  n <- flow$n_final_tumour_epithelial[flow$patient_id == pid]
  add(
    paste0(pid, "_excluded"), "cohort_logic", n, "<150", length(n) == 1L && n < 150L,
    "final tumour-epithelial count evaluated against the published threshold"
  )
}
for (pid in c("UC01", "UC02", "UC03")) {
  n <- flow$n_final_tumour_epithelial[flow$patient_id == pid]
  add(
    paste0(pid, "_retained"), "cohort_logic", n, ">=150", length(n) == 1L && n >= 150L,
    "final tumour-epithelial count evaluated against the published threshold"
  )
}
retained_counts <- flow$n_final_tumour_epithelial[flow$patient_id %in% retained]
add(
  "all_reconstructed_retained_patients_ge150", "cohort_logic",
  min(retained_counts), ">=150", all(retained_counts >= 150L),
  "retained set is recomputed from biological patient counts, not library IDs"
)
add(
  "published_9HV_3UC_structure", "publication_target",
  c(
    HV = sum(flow$histology_group == "HV" & flow$n_final_tumour_epithelial >= 150L),
    UC = sum(flow$histology_group == "UC" & flow$n_final_tumour_epithelial >= 150L)
  ),
  c(HV = 9L, UC = 3L),
  sum(flow$histology_group == "HV" & flow$n_final_tumour_epithelial >= 150L) == 9L &&
    sum(flow$histology_group == "UC" & flow$n_final_tumour_epithelial >= 150L) == 3L,
  "patient-level counts; UC01 is not promoted by lowering the threshold"
)
published_rows <- flow[flow$patient_id %in% expected, , drop = FALSE]
add(
  "published_tumour_cell_total", "publication_target",
  sum(published_rows$n_final_tumour_epithelial), 8553L,
  sum(published_rows$n_final_tumour_epithelial) == 8553L,
  "sum over the same 12 published patient IDs"
)
add(
  "UC01_published_barcodes_in_current_12049_DGE", "data_provenance",
  sum(barcode$found_in_current_12049_DGE), nrow(barcode),
  all(barcode$found_in_current_12049_DGE),
  "exact Figure5B source-data barcodes compared to the deposited 12049 DGE header"
)

pb_cd <- as.data.frame(colData(tumour_pb))
add(
  "no_library_treated_as_independent_patient", "pseudobulk_logic",
  c(columns = ncol(tumour_pb), unique_patients = length(unique(pb_cd$patient_id))),
  c(columns = 12L, unique_patients = 12L),
  ncol(tumour_pb) == length(unique(pb_cd$patient_id)) &&
    !any(c("GSM", "library_id", "tumour_piece", "technical_rep") %in% names(pb_cd)),
  "tumour pseudobulk colData has exactly one column per biological patient"
)
raw_check <- structural$pass[structural$check == "processed_raw_count_sum_preserved"]
add(
  "raw_count_accounting_preserved", "count_integrity",
  raw_check, TRUE, length(raw_check) == 1L && isTRUE(raw_check),
  "canonical structural validation compares processed count sum to QC bladder count sum"
)

out <- do.call(rbind, checks)
write_tsv(out, file.path(audit, "validation_checks_updated.tsv"))
replication_ok <- all(out$pass[out$category %in% c("cohort_logic", "publication_target")])
writeLines(
  c(
    if (replication_ok) "PUBLICATION COHORT REPRODUCED" else "PUBLICATION COHORT NOT REPRODUCED",
    paste("Structural canonical validation:", if (all(structural$pass)) "PASS" else "FAIL"),
    paste("Replication checks passed:", sum(out$pass), "of", nrow(out)),
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  ),
  file.path(audit, "replication_validation_status.txt")
)
case_log(
  "Replication validation: ", sum(out$pass), "/", nrow(out),
  " checks pass; publication cohort reproduced = ", replication_ok,
  file = logf
)

