#!/usr/bin/env Rscript
# Run formula_bench.
#
#   module load R/4.5.3-gfbf-2025b
#   Rscript targets/formula_bench/run.R
#
#   JUST_SCRIPT=true Rscript targets/formula_bench/run.R
#   NAMES=timing_summary Rscript targets/formula_bench/run.R

suppressPackageStartupMessages({
  library(targets)
  library(cli)
  library(glue)
})

case_dir <- "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC"
script_path <- file.path(case_dir, "targets", "formula_bench", "tar_script.R")
root_path <- file.path(case_dir, "data", "processed", "formula_bench")

dir.create(file.path(root_path, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_path, "fits"), recursive = TRUE, showWarnings = FALSE)

setwd(root_path)
cli_alert_info(glue("Store root: {getwd()}"))
source(script_path)

if (identical(Sys.getenv("JUST_SCRIPT"), "true")) {
  cli_alert_success("_targets.R written (JUST_SCRIPT=true).")
  quit(save = "no", status = 0)
}

names_raw <- Sys.getenv("NAMES", unset = "")
if (nzchar(names_raw)) {
  tar_make(
    names = tidyselect::all_of(trimws(strsplit(names_raw, ",", fixed = TRUE)[[1]])),
    callr_function = NULL
  )
} else {
  tar_make(callr_function = NULL)
}

cli_alert_success(glue("Done → {root_path}/timing_summary.csv"))
