#!/usr/bin/env Rscript
# Batch runner for RE vs interaction pilot.
# Run on an allocated COMPUTE node (the random-slopes RE fit is too heavy for the
# login node), e.g.:
#   salloc -c 32 --mem=64G -t 12:00:00
#   module load R/4.5.3-gfbf-2025b
#   REFRESH_FITS=true Rscript R/HCC_re_vs_interaction_pilot.R
#
# Targeted stronger refit example (ZFP36 interaction + rich shape only):
#   GENES=ZFP36 MODELS=interaction SHAPE_VARIANTS=rich \
#   CHAINS=4 WARMUP=1500 ITER=3000 ADAPT_DELTA=0.99 MAX_TREEDEPTH=15 SEED=43 \
#   REFRESH_FITS=true Rscript R/HCC_re_vs_interaction_pilot.R
#
# Progress is streamed with timestamps; also tee to a log to watch it:
#   ... 2>&1 | tee pilot_run.log

suppressPackageStartupMessages({
  library(targets)
  library(brms)
  library(posterior)
  library(SummarizedExperiment)
})

case_dir <- "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC"
out_dir <- file.path(case_dir, "data", "processed", "re_vs_interaction_pilot")

source(file.path(case_dir, "R", "re_vs_interaction_helpers.R"))

parse_csv_env <- function(name) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(NULL)
  trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
}

env_int <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) as.integer(default) else as.integer(raw)
}

env_dbl <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) as.numeric(default) else as.numeric(raw)
}

refresh <- identical(Sys.getenv("REFRESH_FITS"), "true")
genes <- parse_csv_env("GENES")
models <- parse_csv_env("MODELS")
shape_variants <- parse_csv_env("SHAPE_VARIANTS")

chains <- env_int("CHAINS", 2L)
warmup <- env_int("WARMUP", 400L)
iter <- env_int("ITER", 600L)
adapt_delta <- env_dbl("ADAPT_DELTA", 0.95)
max_treedepth <- env_int("MAX_TREEDEPTH", 12L)
seed <- env_int("SEED", 42L)
fit_cores <- env_int("FIT_CORES", 32L)
n_draws <- env_int("N_DRAWS", 400L)

t0 <- Sys.time()
pilot_log(sprintf(
  paste0(
    "Batch runner start (refresh_fits = %s, genes = %s, models = %s, ",
    "shape_variants = %s, chains = %d, warmup = %d, iter = %d, seed = %d, ",
    "adapt_delta = %s, max_treedepth = %d, out_dir = %s)"
  ),
  refresh,
  if (is.null(genes)) "ALL" else paste(genes, collapse = ","),
  if (is.null(models)) "ALL" else paste(models, collapse = ","),
  if (is.null(shape_variants)) "ALL" else paste(shape_variants, collapse = ","),
  chains, warmup, iter, seed, adapt_delta, max_treedepth, out_dir
))

pilot <- run_re_vs_interaction_pilot(
  case_dir = case_dir,
  out_dir = out_dir,
  n_draws = n_draws,
  fit_cores = fit_cores,
  refresh_fits = refresh,
  genes = genes,
  models = models,
  shape_variants = shape_variants,
  chains = chains,
  warmup = warmup,
  iter = iter,
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,
  seed = seed
)

cat("\n=== Baseline summary (per gene x family x shape) ===\n")
print(pilot$baseline_summary)
cat("\n=== Shape effect (simple vs rich SD) ===\n")
print(pilot$shape_effect)

pilot_log(sprintf(
  "Batch runner finished in %s. Outputs in %s",
  fmt_dur(difftime(Sys.time(), t0, units = "secs")), out_dir
))
