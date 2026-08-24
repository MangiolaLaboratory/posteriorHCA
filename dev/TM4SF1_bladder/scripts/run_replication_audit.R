#!/usr/bin/env Rscript
# Reproducible, non-destructive audit layered on the completed canonical workflow.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else normalizePath("scripts")
rscript <- file.path(R.home("bin"), "Rscript")
stages <- c(
  "07_replication_audit.R",
  "08_tumour_replication_de.R",
  "09_validate_replication.R",
  "10_replication_report.R"
)

for (stage in stages) {
  message("\n=== ", stage, " ===")
  status <- system2(rscript, file.path(script_dir, stage))
  if (status != 0L) stop("Audit stage failed: ", stage, call. = FALSE)
}

message("\nReplication audit complete; canonical processed objects were not modified.")
