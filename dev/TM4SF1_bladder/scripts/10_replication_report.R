#!/usr/bin/env Rscript
# Assemble the evidence-backed replication discrepancy report from audit tables.
# This stage is read-only with respect to canonical data/processed objects.

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
source(file.path(.script_dir, "utils.R"))

paths <- case_paths()
audit <- file.path(paths$root, "results", "replication_audit")
dir.create(audit, recursive = TRUE, showWarnings = FALSE)
logf <- start_log(paths, "10_replication_report")

read_tsv <- function(name) {
  utils::read.delim(file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE)
}

md_table <- function(x) {
  x[] <- lapply(x, function(v) {
    v <- ifelse(is.na(v), "", as.character(v))
    gsub("\\|", "\\\\|", v)
  })
  c(
    paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |"),
    apply(x, 1L, function(z) paste0("| ", paste(z, collapse = " | "), " |"))
  )
}

fmt_int <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)
fmt_num <- function(x, digits = 3L) formatC(as.numeric(x), digits = digits, format = "fg")

flow <- read_tsv("all_patient_cell_flow.tsv")
compare <- read_tsv("published_vs_reconstructed_counts.tsv")
uc_lib <- read_tsv("uc01_library_audit.tsv")
uc_clusters <- read_tsv("uc01_cluster_annotation.tsv")
barcode <- read_tsv("uc01_published_barcode_audit.tsv")
de <- read_tsv("TM4SF1_original_DE_result.tsv")
sensitivity <- read_tsv("TM4SF1_downsampling_sensitivity_summary.tsv")
recluster <- read_tsv("tumour_recluster_summary.tsv")
validation <- read_tsv("validation_checks_updated.tsv")
structural <- utils::read.delim(
  file.path(paths$metadata, "validation_checks.tsv"),
  stringsAsFactors = FALSE, check.names = FALSE
)

methods <- data.frame(
  component = c(
    "Raw matrix assembly", "Ambient RNA correction", "Cell QC", "Doublet removal",
    "Whole-bladder clustering", "Broad annotation", "Epithelial selection",
    "Tumour support by InferCNV", "Tumour clustering", "UC01 identity",
    "HV-versus-UC DE"
  ),
  paper_or_public_code = c(
    "Public script starts from private dge_mergedall_032123.rds and omits merge code",
    "decontX is run on the private merged object",
    "Paper thresholds: >=300 genes, >=500 UMIs, <20% mitochondrial reads",
    "Paper names DoubletFinder; public parameters and calls are absent",
    "Public scripts use 3,000 HVGs, 100 PCs, k=30, final resolution 0.4",
    "Private label transfer, external marker files, and manual Final_ID cluster overrides",
    "cells_epithelial is taken from an unreleased Epithelial object",
    "Paper reports InferCNV 1.4.0/HMM; public script only prepares input from private labels",
    "Paper Methods specify 2,000 HVGs, top 75 PCs, resolution 0.5",
    "Pre-2026 code uses 12923; commit d26f042 changes it to 12049; other scripts retain PureUC12923",
    "Released Figure5A gives the final DE table; exact seed and complete input object are absent"
  ),
  audit_implementation = c(
    "Read and union all 67 deposited raw DGE matrices",
    "Not applied to canonical raw counts; omission cannot create missing barcodes",
    "Applied sequentially and recorded per library and patient",
    "scDblFinder-style closest reproducible implementation; scores/calls retained",
    "Canonical reconstruction follows paper Methods: 2,000 HVGs, 75-PC graph, resolution 0.5",
    "Multi-marker module evidence for epithelial, immune, stromal, endothelial classes",
    "Marker-supported broad epithelial cells; no access to authors' private membership list",
    "Prepared a reproducible input only; no HMM result or keep/drop list is claimed",
    "Reclustered 7,862 cells using 2,000 HVGs, 75-PC graph, resolution 0.5",
    "Audited filename, GEO title/description, Git history, and Figure5B source barcodes",
    "Compared released result; ran fixed-seed and 100-seed public-DGE sensitivity analyses"
  ),
  exact_reproduction = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  implication = c(
    "The public DGE deposit is the only reproducible starting point",
    "May alter expression/annotation, not cell identities or the absent UC01 source barcodes",
    "UC01 loss is not caused by threshold misapplication",
    "Likely contributes to the 249-cell global post-QC difference, but not the 605-barcode mismatch",
    "Paper and public code disagree; both settings are documented",
    "Exact authors' labels cannot be recreated from public materials",
    "Exact 8,553-cell membership cannot be reconstructed from code alone",
    "Cannot adjudicate the missing UC01 sample",
    "Settings were reproduced, but the input cohort and Cluster-13-like state were not",
    "The deposited matrix assigned to UC01 is not the matrix behind the published UC01 cells",
    "TM4SF1 direction is robust in available data; exact published genome-wide test is not reproducible"
  ),
  stringsAsFactors = FALSE
)
write_tsv(methods, file.path(audit, "authors_code_method_audit.tsv"))

published_ids <- c(paste0("VAR0", 1:9), paste0("UC0", 1:3))
published_flow <- flow[match(published_ids, flow$patient_id), ]
retained <- published_flow$n_final_tumour_epithelial >= 150L
hv_n <- sum(retained & published_flow$histology_group == "HV")
uc_n <- sum(retained & published_flow$histology_group == "UC")
tumour_total <- sum(published_flow$n_final_tumour_epithelial)
retained_total <- sum(published_flow$n_final_tumour_epithelial[retained])
all_tumour_total <- sum(flow$n_final_tumour_epithelial)

bladder_totals <- colSums(flow[, c(
  "n_raw", "n_after_gene", "n_after_umi", "n_after_mito", "n_after_doublet",
  "n_processed", "n_broad_epithelial", "n_candidate_tumour_epithelial",
  "n_final_tumour_epithelial"
)])
global_flow <- data.frame(
  stage = c(
    "Raw bladder barcodes", ">=300 detected genes", ">=500 UMIs",
    "<20% mitochondrial reads", "After doublet filter", "Processed bladder cells",
    "Broad epithelial", "Candidate tumour epithelial", "Final tumour epithelial"
  ),
  cells = as.integer(bladder_totals), stringsAsFactors = FALSE
)
global_flow$cells <- fmt_int(global_flow$cells)

uc01 <- flow[flow$patient_id == "UC01", ]
uc01_flow <- data.frame(
  stage = global_flow$stage,
  cells = as.integer(uc01[, c(
    "n_raw", "n_after_gene", "n_after_umi", "n_after_mito", "n_after_doublet",
    "n_processed", "n_broad_epithelial", "n_candidate_tumour_epithelial",
    "n_final_tumour_epithelial"
  )]), stringsAsFactors = FALSE
)

cluster32 <- uc_clusters[uc_clusters$seurat_cluster == 32, ]
major_non_epi <- uc_clusters[uc_clusters$seurat_cluster %in% c(8, 10, 14, 21),
                             c("seurat_cluster", "n_cells", "current_broad_cell_type", "annotation_confidence")]
candidate <- recluster[recluster$audit_cluster13_candidate, ]

paper_de <- de[de$analysis == "paper_released_source_DE_9HV_3UC_150_cells_per_patient_unknown_seed", ]
closest_de <- de[de$analysis == "closest_public_DGE_9HV_2UC_150_cells_per_patient_seed293189", ]
secondary_de <- de[de$analysis == "secondary_public_DGE_9HV_3UC_90_cells_per_patient_seed293189", ]

lines <- c(
  "# GSE293189 replication discrepancy audit",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Outcome",
  "",
  paste0(
    "Across the 12 patient IDs in the published cohort, the public deposit yields **",
    fmt_int(tumour_total), " reconstructed tumour epithelial cells**, including 90 UC01 cells ",
    "that fail the 150-cell rule. Applying that rule leaves **", hv_n, " HV + ", uc_n,
    " UC patients / ", fmt_int(retained_total), " cells**, not the published 9 HV + 3 UC / ",
    "8,553-cell cohort. Across all 15 tumours, including the three expected exclusions, the ",
    "reconstruction contains ", fmt_int(all_tumour_total), " candidate tumour epithelial cells. ",
    "The discrepancy is localized primarily to UC01: ",
    "90 reconstructed cells versus 605 published cells. No canonical object or pseudobulk ",
    "matrix was changed, because neither relabelling non-epithelial cells nor lowering the ",
    "150-cell threshold would be scientifically justified."
  ),
  "",
  paste0(
    "The strongest provenance result is exact: **0 of 605** Figure5B UC01 source barcodes ",
    "occur in the currently filename-assigned 12049 DGE, and 0 occur in any library assigned ",
    "to current UC01. Only ", sum(barcode$found_in_any_deposited_DGE),
    " barcode strings occur anywhere among all 67 deposited matrices, each as an unrelated ",
    "single-library collision. In contrast, the recoverable source barcodes for UC02 and UC03 ",
    "map back to their assigned deposited libraries. This is evidence of a missing or different ",
    "UC01 input matrix, not a routine QC or annotation error."
  ),
  "",
  "## 1. Cell flow",
  "",
  "### Bladder cohort overall",
  "",
  md_table(global_flow),
  "",
  "The authors report 21,533 post-QC cells globally. The public reconstruction has 21,782 bladder cells after doublet filtering (249 more). Starting from the same 21,970 post-mitochondrial cells, matching 21,533 would require removal of 437 cells rather than 188. Missing DoubletFinder parameters can plausibly explain much of this global difference, but cannot explain why all 605 published UC01 barcodes are absent.",
  "",
  "### UC01/current 12049 library",
  "",
  md_table(uc01_flow),
  "",
  "The 1,142 UC01 QC-passing cells consist of 90 epithelial, 877 immune, 134 stromal, and 41 endothelial cells. Cluster 32 supplies 88 of the 90 epithelial cells and has moderate multi-marker confidence. Its epithelial support includes EPCAM detection in 73.9% of cells and KRT8/KRT18/KRT19 detection in 68.2%/71.6%/65.9%. Major excluded clusters show coherent alternative identities (for example PTPRC in 94.3% of cluster 10, ACTA2 in 100% of cluster 21, and VWF in 96.9% of cluster 14). There is therefore no hidden 515-cell epithelial population that can be recovered by a reasonable broad-marker reinterpretation.",
  "",
  "Full evidence: `uc01_cell_flow.tsv`, `uc01_cluster_annotation.tsv`, `uc01_marker_summary.tsv`, `uc01_umap_annotations.pdf`, and `uc01_marker_featureplots.pdf`.",
  "",
  "## 2. UC01 identity and metadata conflict",
  "",
  md_table(uc_lib[, c(
    "GSM", "library_id", "patient_id", "original_patient_id", "geo_title",
    "geo_original_patient_id", "uc01_claim_basis", "downloaded"
  )]),
  "",
  "The original public code mapped UC1 to patient 12923. Commit `d26f0427481c5ea39d1697ec3c83b7921a747504` (16 January 2026) changed only this assignment to 12049, while downstream scripts still refer to `PureUC12923`. Patient 12923 is not represented by a deposited DGE filename. GEO record GSM8878368 carries a 12049 filename but a 11734 title/description; GSM8878370 has a 12049 title/description but carries a 12050T1 DGE filename. These conflicts are preserved, not silently resolved.",
  "",
  "## 3. Published versus reconstructed tumour cohort",
  "",
  md_table(compare[, c(
    "patient_id", "histology_group", "n_final_tumour_epithelial",
    "published_tumour_epithelial_cells", "delta_tumour_epithelial", "reconstructed_ge150"
  )]),
  "",
  "All nine expected HV patients pass the 150-cell threshold, as do UC02 and UC03. UC01 does not. VAR10, VAR11, and UC04 remain correctly excluded. Validation records this as a publication-replication failure rather than changing the threshold or treating libraries as patients.",
  "",
  "## 4. Authors' code and method audit",
  "",
  md_table(methods[, c("component", "paper_or_public_code", "audit_implementation", "exact_reproduction", "implication")]),
  "",
  "The exact author object cannot be regenerated from the repository because the raw-matrix merge, `dge_mergedall_032123.rds`, reference label object, external marker lists, final epithelial object, manual keep/drop decisions, and InferCNV outputs are not public. DecontX could change gene counts and downstream label scores, but it cannot produce the missing published UC01 barcodes. Exact InferCNV support is likewise not reproducible: the public code prepares an input from private `Final_ID` labels but releases no gene-order file, HMM calls, result object, or epithelial membership decision.",
  "",
  "## 5. Tumour reclustering and the reported Cluster 13 state",
  "",
  paste0(
    "The 7,862 cells from the expected 12 patient IDs were reclustered with 2,000 HVGs, ",
    "PCA, the top 75 PCs, k=30, and resolution 0.5. The closest signature-scoring candidate ",
    "was cluster ", candidate$tumour_cluster, " (", fmt_int(candidate$n_cells), " cells), but ",
    candidate$top_patient, " contributed ", fmt_num(100 * candidate$top_patient_fraction, 3),
    "% of it; VAR05 contributed 111 cells and UC01 one cell. This patient-dominated composition ",
    "does not match a convincing cross-patient recovery of the authors' Cluster 13 state. ",
    "Candidate markers include KRT6A, KRT5, S100A2, FN1, VIM, and other basal/mesenchymal genes."
  ),
  "",
  "See `tumour_recluster_summary.tsv`, `tumour_recluster_composition.tsv`, `tumour_cluster13_candidate_markers.tsv`, `tumour_reclustering_cluster13_candidate.pdf`, and `tumour_signature_expression.pdf`.",
  "",
  "## 6. TM4SF1 differential expression and downsampling",
  "",
  md_table(de[, c(
    "analysis", "n_HV_patients", "n_UC_patients", "n_cells_per_patient",
    "avg_log2FC", "p_val", "p_val_adj", "FDR_BH", "pct.1", "pct.2",
    "HV_enriched_rank_by_p"
  )]),
  "",
  paste0(
    "The released paper source table places TM4SF1 first among HV-enriched genes by p value ",
    "(log2FC ", fmt_num(paper_de$avg_log2FC), ", p=", format(paper_de$p_val, scientific = TRUE, digits = 3),
    ", source adjusted p=", format(paper_de$p_val_adj, scientific = TRUE, digits = 3),
    "). In the closest feasible public-DGE 150-cell comparison (9 HV + 2 UC), it is also rank 1 ",
    "(log2FC ", fmt_num(closest_de$avg_log2FC), ", p=", format(closest_de$p_val, scientific = TRUE, digits = 3),
    ", BH FDR=", format(closest_de$FDR_BH, scientific = TRUE, digits = 3),
    "). Across 100 seeds, its median rank is 1 (range 1-2; 79% rank 1 and 100% top 5)."
  ),
  "",
  paste0(
    "A secondary 9 HV + 3 UC comparison is possible only by lowering every patient to 90 cells. ",
    "At the fixed seed TM4SF1 ranks ", secondary_de$HV_enriched_rank_by_p,
    "; over 100 seeds its median rank is 15 (range 5-32; 1% top 5). This is a sensitivity analysis, ",
    "not a reproduction of the requested 150-cell design. Using the released per-cell TM4SF1 values ",
    "for the exact 12 published patients, the HV-minus-UC signal remains positive in all 100 resamples, ",
    "but genome-wide rank sensitivity cannot be calculated because the source workbook releases only ",
    "per-cell TM4SF1 for this purpose, not the complete underlying count matrix."
  ),
  "",
  md_table(sensitivity),
  "",
  "## 7. Validation and disposition",
  "",
  paste0(
    "Canonical structural validation: **", if (all(structural$pass)) "PASS" else "FAIL",
    "** (", sum(structural$pass), "/", nrow(structural), "). Replication audit validation: **",
    sum(validation$pass), "/", nrow(validation), " pass**. The four failures are the expected ",
    "UC01/publication-target failures and are retained in `validation_checks_updated.tsv`."
  ),
  "",
  "No corrected reconstruction was emitted. The existing patient-level pseudobulk remains valid as a transparent public-DGE reconstruction, with UC01 explicitly flagged below 150 cells; it must not be represented as an exact reproduction of the paper's 9 HV + 3 UC / 8,553-cell cohort. Resolving the discrepancy requires the authors' original UC01/12923 count matrix or the exact private epithelial object and cell-membership list.",
  "",
  "## 8. Reproducible outputs",
  "",
  "- Audit driver: `scripts/run_replication_audit.R`",
  "- UC01 and provenance audit: `scripts/07_replication_audit.R`",
  "- Reclustering, DE, and sensitivity: `scripts/08_tumour_replication_de.R`",
  "- Validation: `scripts/09_validate_replication.R`",
  "- This report: `scripts/10_replication_report.R`",
  "- Machine-readable method comparison: `results/replication_audit/authors_code_method_audit.tsv`",
  "- Original-object checksums: `results/replication_audit/original_reconstruction/SHA256SUMS.txt`",
  "",
  "All audit tables, plots, and the audit-only reclustered object are under `results/replication_audit/`. Canonical data remain under `data/processed/`; no posteriorHCA integration or tumour/toxicity threshold was attempted."
)

writeLines(lines, file.path(audit, "replication_discrepancy_report.md"))
case_log(
  "Wrote replication_discrepancy_report.md; canonical objects unchanged",
  file = logf
)
