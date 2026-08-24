#!/usr/bin/env Rscript
# Initialise / refresh the strong RE-vs-interaction targets pipeline and run it.
#
#   module load R/4.5.3-gfbf-2025b
#   Rscript targets/run_re_vs_interaction_strong.R
#
# Optional:
#   JUST_SCRIPT=true  Rscript ...   # only write _targets.R, do not tar_make
#   NAMES=fit_SPON2_re  Rscript ... # tar_make selected names only

suppressPackageStartupMessages({
  library(targets)
  library(cli)
  library(glue)
})

case_dir <- "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC"
script_path <- file.path(case_dir, "targets", "re_vs_interaction_strong_tar_script.R")
root_path <- "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/re_vs_interaction_strong"
out_dir <- file.path(case_dir, "data", "processed", "re_vs_interaction_strong")

dir.create(root_path, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

setwd(root_path)
cli_alert_info(glue("Working directory: {getwd()}"))
cli_alert_info(glue("Sourcing tar script: {script_path}"))

source(script_path)

just_script <- identical(Sys.getenv("JUST_SCRIPT"), "true")
if (just_script) {
  cli_alert_success("_targets.R written. Skipping tar_make (JUST_SCRIPT=true).")
  cli_alert_info(glue("Next: cd {root_path} && R -e 'targets::tar_make()'"))
  quit(save = "no", status = 0)
}

names_raw <- Sys.getenv("NAMES", unset = "")
tar_names <- if (nzchar(names_raw)) {
  trimws(strsplit(names_raw, ",", fixed = TRUE)[[1]])
} else {
  NULL
}

cli_alert_info(
  if (is.null(tar_names)) {
    "Starting tar_make() for the full pipeline (12 fits) ..."
  } else {
    glue("Starting tar_make(names = {paste(tar_names, collapse = ', ')}) ...")
  }
)

if (is.null(tar_names)) {
  tar_make(callr_function = NULL)
} else {
  tar_make(names = tidyselect::all_of(tar_names), callr_function = NULL)
}

cli_alert_success(glue("Done. Fits / CSVs also under {out_dir}"))
cli_alert_info(glue("Targets store: {file.path(root_path, '_targets')}"))
