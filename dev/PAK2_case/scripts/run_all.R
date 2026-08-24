#!/usr/bin/env Rscript
# Run the PAK2 data-preparation pipeline in order.

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
.script_dir <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
source(file.path(.script_dir, "utils.R"))

paths <- pak2_paths()
logf <- start_log(paths, "run_all")

scripts <- c(
  "00_download_geo.R",
  "01_build_metadata.R",
  "02_build_spatial_seurat.R",
  "03_build_scrna_seurat.R",
  "04_pseudobulk_spatial.R",
  "05_pseudobulk_scrna.R"
)

for (s in scripts) {
  pak2_log("===== ", s, " =====", .log_file = logf)
  status <- system2("Rscript", file.path(.script_dir, s))
  if (!identical(status, 0L)) stop("Pipeline failed at ", s, " (status ", status, ")")
}

sink(file.path(paths$logs, "sessionInfo.txt"))
print(utils::sessionInfo())
sink()
pak2_log("Pipeline finished successfully.", .log_file = logf)
