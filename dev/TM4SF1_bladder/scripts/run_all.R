#!/usr/bin/env Rscript
# End-to-end, restartable driver. Existing valid downloads and QC checkpoints
# are reused; delete a checkpoint deliberately if that stage must be recomputed.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else normalizePath("scripts")
rscript <- file.path(R.home("bin"), "Rscript")
stages <- c(
  "00_download.R",
  "01_build_manifest.R",
  "02_qc_reconstruct.R",
  "03_cluster_annotate.R",
  "04_pseudobulk_tm4sf1.R",
  "05_report.R",
  "06_validate_outputs.R"
)

for (stage in stages) {
  message("\n=== ", stage, " ===")
  status <- system2(rscript, file.path(script_dir, stage))
  if (status != 0L) stop("Workflow stage failed: ", stage, call. = FALSE)
}

message("\nWorkflow complete.")
