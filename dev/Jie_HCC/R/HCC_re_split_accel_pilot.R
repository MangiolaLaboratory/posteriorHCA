#!/usr/bin/env Rscript
# One-gene RE pilot: split correlated tissue random-effects block vs production baseline.
#
# Both variants use random-effect counts + simple shape (revised production candidate).
# The only structural change is how tissue random slopes are parameterised:
#
#   baseline (V1_nk production counts formula):
#     (1 + age_decade * sex + ethnicity_groups | tissue_groups)
#
#   split (benchmark — two smaller, separately correlated blocks):
#     (1 + age_decade * sex | tissue_groups)
#   + (1 + ethnicity_groups | tissue_groups)
#
# Baseline reference fit: load the stored V1_nk production brms fit from the NK
# targets store (rich shape — for timing / Rhat comparison only).
#
# Optional acceleration knobs (manual testing):
#   ALGORITHM=pathfinder|sampling
#   USE_OPENCL=true|false
#
# Examples (compute node recommended):
#   module load R/4.5.3-gfbf-2025b
#   cd /home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC
#
#   # fit split RE only
#   MODEL=re_split Rscript R/HCC_re_split_accel_pilot.R
#
#   # fit baseline RE (simple shape refit) for apples-to-apples benchmark
#   MODEL=re_baseline Rscript R/HCC_re_split_accel_pilot.R
#
#   # fit both back-to-back and write timing summary
#   RUN_BOTH=true ALGORITHM=sampling CHAINS=2 Rscript R/HCC_re_split_accel_pilot.R
#
# Env overrides:
#   GENE=ZFP36  MODEL=re_split|re_baseline  RUN_BOTH=true|false
#   ALGORITHM=pathfinder  USE_OPENCL=true  USE_SPARSE=false
#   CHAINS=2 WARMUP=400 ITER=600 SEED=42 FIT_CORES=8

suppressPackageStartupMessages({
  library(targets)
  library(brms)
  library(cmdstanr)
  library(dplyr)
  library(stringr)
  library(readr)
  library(glue)
  library(tibble)
})

case_dir <- "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC"
model_store_root <- "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE"
target_cell_type <- "nk"

source(file.path(case_dir, "R", "re_vs_interaction_helpers.R"))
specs <- formula_specs()

gene <- Sys.getenv("GENE", unset = "ZFP36")
# re_baseline = production counts RE block | re_split = two smaller RE blocks
model_name <- Sys.getenv("MODEL", unset = "re_split")
run_both <- identical(Sys.getenv("RUN_BOTH", unset = "false"), "true")
algorithm <- Sys.getenv("ALGORITHM", unset = "pathfinder")
use_sparse <- identical(Sys.getenv("USE_SPARSE", unset = "false"), "true")
use_opencl <- identical(Sys.getenv("USE_OPENCL", unset = "true"), "true")
opencl_ids <- c(0L, 0L)

chains_default <- if (identical(Sys.getenv("ALGORITHM", unset = "pathfinder"), "pathfinder")) "50" else "2"
chains <- as.integer(Sys.getenv("CHAINS", unset = chains_default))
warmup <- as.integer(Sys.getenv("WARMUP", unset = "400"))
iter <- as.integer(Sys.getenv("ITER", unset = "600"))
seed <- as.integer(Sys.getenv("SEED", unset = "42"))
fit_cores <- as.integer(Sys.getenv("FIT_CORES", unset = "8"))
adapt_delta <- as.numeric(Sys.getenv("ADAPT_DELTA", unset = "0.95"))
max_treedepth <- as.integer(Sys.getenv("MAX_TREEDEPTH", unset = "12"))
pathfinder_history_size <- as.integer(Sys.getenv("PF_HISTORY_SIZE", unset = "100"))
pathfinder_max_lbfgs_iters <- as.integer(Sys.getenv("PF_MAX_LBFGS", unset = "100"))

out_dir <- file.path(case_dir, "data", "processed", "re_split_accel_pilot")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Formulas — RE + simple shape only
# ---------------------------------------------------------------------------
shape_formula <- specs$re$shape_simple

# Production V1_nk counts formula (single correlated tissue RE block).
counts_re_baseline <- specs$re$counts

# Split tissue RE: (1 + age*sex) and (1 + ethnicity) as separate correlated blocks.
counts_re_split <- paste(
  "counts ~ 1 + offset(offset) + age_decade * sex + disease_groups_altered",
  "+ ethnicity_groups + assay_groups_altered + (1 | dataset_id_altered)",
  "+ (1 + age_decade * sex | tissue_groups)",
  "+ (0 + ethnicity_groups | tissue_groups)"
)

formula_for_model <- function(name) {
  switch(
    name,
    re_baseline = list(counts = counts_re_baseline, label = "re_baseline"),
    re_split = list(counts = counts_re_split, label = "re_split"),
    stop("Unknown MODEL=", name, " (use re_baseline|re_split)")
  )
}

models_to_run <- if (run_both) c("re_baseline", "re_split") else model_name

message(glue(
  "=== RE split pilot ===\n",
  "  gene={gene}  models={paste(models_to_run, collapse=', ')}\n",
  "  algorithm={algorithm}  sparse={use_sparse}  opencl={use_opencl}\n",
  "  chains={chains} warmup={warmup} iter={iter} seed={seed} cores={fit_cores}\n",
  "  shape (both): {shape_formula}"
))

# ---------------------------------------------------------------------------
# Load gene data + production NK reference fit
# ---------------------------------------------------------------------------
store <- file.path(model_store_root, paste0("V1_", target_cell_type), "_targets")

model_meta <- tar_meta(starts_with("estimates_chunk_"), store = store) |>
  filter(.data$type == "branch") |>
  mutate(.feature = str_extract(.data$warnings, "(?<=Gene:___)ENSG\\d+(?=___)"))

biomarkers <- read_csv(
  file.path(case_dir, "data", "times_biomarkers.csv"),
  show_col_types = FALSE
) |>
  filter(.data$gene == !!gene, .data$.feature %in% model_meta$.feature)

if (!nrow(biomarkers)) {
  stop("Gene ", gene, " not found in biomarkers x V1_nk store metadata.")
}

feature_id <- biomarkers$.feature[[1]]
target_name <- model_meta |>
  filter(.data$.feature == feature_id) |>
  slice(1) |>
  pull(.data$name)

message(glue("Loading V1_nk production fit + data: {gene} ({feature_id})"))
fit_prod <- tar_read_raw(target_name[[1]], store = store)$brms_fit[[1]]
dat <- droplevels(fit_prod$data)
colnames(dat) <- str_replace_all(colnames(dat), "_+", "_")
message(glue("  n = {nrow(dat)}, tissues = {n_distinct(dat$tissue_groups)}"))
message(glue("  production store formula (reference): {fit_prod$formula$formula}"))

prod_ref_path <- file.path(out_dir, paste0("ref_production_V1_nk_", gene, ".rds"))
if (!file.exists(prod_ref_path)) {
  saveRDS(fit_prod, prod_ref_path)
  message(glue("Saved production reference -> {prod_ref_path}"))
}

# ---------------------------------------------------------------------------
# Fit helper
# ---------------------------------------------------------------------------
fit_re_model <- function(counts_formula, label) {
  i <- mean(log1p(dat$counts / exp(dat$offset)))
  prior <- eval(substitute(
    c(
      prior(student_t(3, i, 1.5), class = Intercept),
      prior(student_t(3, 0, 1), class = Intercept, dpar = shape),
      prior(student_t(3, 0, 5), class = b),
      prior(student_t(3, 0, 2), class = b, dpar = shape)
    ),
    env = list(i = i)
  ))

  form <- bf(
    as.formula(counts_formula),
    as.formula(shape_formula),
    sparse = use_sparse
  )

  out_tag <- paste(
    gene, label, "simple", algorithm,
    if (use_sparse) "sparse" else "dense",
    if (use_opencl) "opencl" else "cpu",
    sep = "_"
  )
  out_rds <- file.path(out_dir, paste0("fit_", out_tag, ".rds"))

  brm_args <- list(
    formula = form,
    data = dat,
    family = zero_inflated_negbinomial(),
    prior = prior,
    chains = chains,
    cores = min(fit_cores, chains),
    warmup = warmup,
    iter = iter,
    refresh = 50,
    backend = "cmdstanr",
    algorithm = algorithm,
    silent = 1,
    seed = seed,
    file = out_rds,
    file_refit = "always"
  )

  if (use_opencl) {
    brm_args$opencl <- opencl(opencl_ids)
    message(glue("OpenCL: platform/device = ({opencl_ids[1]}, {opencl_ids[2]})"))
  }
  if (identical(algorithm, "sampling")) {
    brm_args$control <- list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  }
  if (identical(algorithm, "pathfinder")) {
    brm_args$history_size <- pathfinder_history_size
    brm_args$max_lbfgs_iters <- pathfinder_max_lbfgs_iters
    brm_args$psis_resample <- FALSE
    message(glue(
      "Pathfinder (sccomp-style): num_paths=chains={chains}, ",
      "history_size={pathfinder_history_size}, max_lbfgs_iters={pathfinder_max_lbfgs_iters}"
    ))
  }

  message(glue("\n--- {label} ---\ncounts: {counts_formula}"))
  message(glue("Fitting -> {out_rds}"))
  t0 <- Sys.time()
  fit <- do.call(brm, brm_args)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(glue("Done in {round(elapsed, 1)}s"))

  tibble(
    gene = gene,
    model = label,
    algorithm = algorithm,
    sparse = use_sparse,
    opencl = use_opencl,
    chains = chains,
    warmup = warmup,
    iter = iter,
    seed = seed,
    elapsed_sec = elapsed,
    out_rds = out_rds,
    counts_formula = counts_formula,
    shape_formula = shape_formula
  )
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
timing_rows <- list()
fits <- list()

for (m in models_to_run) {
  spec <- formula_for_model(m)
  timing_rows[[m]] <- fit_re_model(spec$counts, spec$label)
  fits[[m]] <- readRDS(timing_rows[[m]]$out_rds)
  print(fits[[m]])
}

timing <- bind_rows(timing_rows)
readr::write_csv(timing, file.path(out_dir, paste0("timing_", gene, ".csv")))
message(glue("\nTiming summary -> {file.path(out_dir, paste0('timing_', gene, '.csv'))}"))
print(timing)

message(glue(
  "\nBenchmark notes:\n",
  "  - re_baseline: production NK counts RE (one correlated block) + simple shape\n",
  "  - re_split:    two tissue RE blocks + simple shape\n",
  "  - ref_production_V1_nk_*: stored NK fit (rich shape) for external reference"
))
