#!/usr/bin/env Rscript
# Systematic combo grid: RE form × sparse × algorithm × OpenCL
# Based on HCC_opencl_one_gene_pilot.R settings (ZFP36, V1_nk).
#
# For GPU/OpenCL cells, activate the custom conda env first (sets OCL_ICD_VENDORS):
#   source ~/.bashrc && use_R_env
#   cd /home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC
#   Rscript R/HCC_accel_combo_grid.R 2>&1 | tee data/processed/re_split_accel_pilot/combo_grid.log
#
# GPU-only retest (dense + OpenCL): R/HCC_accel_gpu_retest.R
#
# HMC cells use a short smoke schedule so the grid finishes; Pathfinder uses
# reduced paths vs full sccomp (20) for the same reason. Flip FULL=true for
# production-like settings (much slower).

suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(targets)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
  library(tidyr)
})

options(brms.backend = "cmdstanr")
options(brms.opencl = NULL)
cmdstanr::check_cmdstan_toolchain()

full <- identical(Sys.getenv("FULL", unset = "false"), "true")

case_dir <- "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC"
model_store_root <- "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE"
store <- file.path(model_store_root, "V1_nk", "_targets")
gene <- "ZFP36"
out_dir <- file.path(case_dir, "data", "processed", "re_split_accel_pilot")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
results_csv <- file.path(out_dir, "accel_combo_grid_results.csv")
fits_dir <- file.path(out_dir, "combo_fits")
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Formulas
# ---------------------------------------------------------------------------
counts_baseline <- paste(
  "counts ~ 1 + offset(offset) + age_decade * sex + disease_groups_altered",
  "+ ethnicity_groups + assay_groups_altered + (1 | dataset_id_altered)",
  "+ (1 + age_decade * sex + ethnicity_groups | tissue_groups)"
)
counts_split <- paste(
  "counts ~ 1 + offset(offset) + age_decade * sex + disease_groups_altered",
  "+ ethnicity_groups + assay_groups_altered + (1 | dataset_id_altered)",
  "+ (1 + age_decade * sex | tissue_groups)",
  "+ (0 + ethnicity_groups | tissue_groups)"
)
shape_formula <- "shape ~ 1 + disease_groups_altered + assay_groups_altered + (1 | tissue_groups)"

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
model_meta <- tar_meta(starts_with("estimates_chunk_"), store = store) |>
  filter(.data$type == "branch") |>
  mutate(.feature = str_extract(.data$warnings, "(?<=Gene:___)ENSG\\d+(?=___)"))

biomarkers <- read_csv(
  file.path(case_dir, "data", "times_biomarkers.csv"),
  show_col_types = FALSE
) |>
  filter(.data$gene == !!gene, .data$.feature %in% model_meta$.feature)

stopifnot(nrow(biomarkers) >= 1L)
feature_id <- biomarkers$.feature[[1]]
target_name <- model_meta |>
  filter(.data$.feature == feature_id) |>
  slice(1) |>
  pull(.data$name)

message("Loading data: ", gene, " (", feature_id, ")")
fit_prod <- tar_read_raw(target_name[[1]], store = store)$brms_fit[[1]]
dat <- droplevels(fit_prod$data)
colnames(dat) <- str_replace_all(colnames(dat), "_+", "_")
message("n = ", nrow(dat), ", tissues = ", n_distinct(dat$tissue_groups))

i <- mean(log1p(dat$counts / exp(dat$offset)))
prior_dense <- eval(substitute(
  c(
    prior(student_t(3, i, 1.5), class = Intercept),
    prior(student_t(3, 0, 1), class = Intercept, dpar = shape),
    prior(student_t(3, 0, 5), class = b),
    prior(student_t(3, 0, 2), class = b, dpar = shape)
  ),
  env = list(i = i)
))
prior_sparse <- eval(substitute(
  c(
    prior(student_t(3, i, 1.5), class = b, coef = Intercept),
    prior(student_t(3, 0, 1), class = Intercept, dpar = shape),
    prior(student_t(3, 0, 5), class = b),
    prior(student_t(3, 0, 2), class = b, dpar = shape)
  ),
  env = list(i = i)
))

# ---------------------------------------------------------------------------
# Grid: 2 × 2 × 2 × 2 = 16
# ---------------------------------------------------------------------------
grid <- tidyr::expand_grid(
  re = c("baseline", "split"),
  sparse = c(FALSE, TRUE),
  algorithm = c("pathfinder", "sampling"),
  opencl = c(FALSE, TRUE)
) |>
  mutate(
    cell_id = sprintf(
      "%02d_%s_%s_%s_%s",
      row_number(),
      re,
      if_else(sparse, "sparse", "dense"),
      algorithm,
      if_else(opencl, "opencl", "cpu")
    )
  )

message("=== accel combo grid: ", nrow(grid), " cells (FULL=", full, ") ===")
message(
  "brms=", as.character(packageVersion("brms")),
  " cmdstanr=", as.character(packageVersion("cmdstanr")),
  " cmdstan=", cmdstan_version()
)

# Pathfinder / HMC schedules
pf_chains <- if (full) 50L else 20L
pf_draws <- if (full) 80L else 40L
hmc_chains <- if (full) 2L else 1L
hmc_warmup <- if (full) 400L else 50L
hmc_iter <- if (full) 600L else 100L

run_one <- function(cell) {
  tag <- cell$cell_id
  message("\n========== ", tag, " ==========")
  counts <- if (identical(cell$re, "split")) counts_split else counts_baseline
  form <- bf(
    as.formula(counts),
    as.formula(shape_formula),
    sparse = isTRUE(cell$sparse)
  )
  pr <- if (isTRUE(cell$sparse)) prior_sparse else prior_dense
  out_rds <- file.path(fits_dir, paste0("fit_", tag, ".rds"))

  # Unique stanvar so OpenCL vs CPU / sparse vs dense do not silently reuse wrong exe
  sv <- stanvar(
    scode = paste0("// combo_", tag),
    block = "functions"
  )

  args <- list(
    formula = form,
    data = dat,
    family = zero_inflated_negbinomial(),
    prior = pr,
    stanvars = sv,
    backend = "cmdstanr",
    algorithm = cell$algorithm,
    seed = 20260720L,
    refresh = 0,
    silent = 2,
    file = out_rds,
    file_refit = "always"
  )

  if (isTRUE(cell$opencl)) {
    args$opencl <- brms::opencl(c(0L, 0L))
  } else {
    args$opencl <- NULL
  }

  if (identical(cell$algorithm, "pathfinder")) {
    args$chains <- pf_chains
    args$cores <- 1L
    args$history_size <- 100L
    args$max_lbfgs_iters <- 100L
    args$single_path_draws <- pf_draws
    args$psis_resample <- FALSE
    # threads forbidden with sparse
    if (!isTRUE(cell$sparse)) {
      args$threads <- threading(8)
    }
  } else {
    args$chains <- hmc_chains
    args$cores <- min(hmc_chains, 2L)
    args$warmup <- hmc_warmup
    args$iter <- hmc_iter
    if (!isTRUE(cell$sparse)) {
      args$threads <- threading(4)
    }
  }

  t0 <- Sys.time()
  err <- NA_character_
  status <- "ok"
  fit <- tryCatch(
    do.call(brm, args),
    error = function(e) {
      err <<- conditionMessage(e)
      status <<- "error"
      NULL
    }
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Truncate long Stan compile errors for the CSV
  if (!is.na(err) && nchar(err) > 500) {
    err <- paste0(substr(err, 1L, 500L), "...")
  }

  # Classify known failure modes for reporting
  fail_class <- NA_character_
  if (identical(status, "error")) {
    fail_class <- dplyr::case_when(
      grepl("old array syntax|int vX\\[", err) ~ "sparse_stan_syntax",
      grepl("threading and sparse", err, ignore.case = TRUE) ~ "sparse_threads_forbidden",
      grepl("OpenCL|opencl", err) ~ "opencl_runtime",
      grepl("path__", err) ~ "brms_pathfinder_parse",
      grepl("prior", err, ignore.case = TRUE) ~ "prior_mismatch",
      TRUE ~ "other"
    )
  }

  n_div <- NA_real_
  max_rhat <- NA_real_
  if (!is.null(fit) && identical(cell$algorithm, "sampling")) {
    sm <- tryCatch(summary(fit)$fixed, error = function(e) NULL)
    max_rhat <- tryCatch(
      max(posterior::rhat(fit), na.rm = TRUE),
      error = function(e) NA_real_
    )
  }

  row <- tibble(
    cell_id = tag,
    gene = gene,
    re = cell$re,
    sparse = cell$sparse,
    algorithm = cell$algorithm,
    opencl = cell$opencl,
    status = status,
    fail_class = fail_class,
    elapsed_sec = elapsed,
    error = err,
    out_rds = if (file.exists(out_rds)) out_rds else NA_character_,
    brms = as.character(packageVersion("brms")),
    cmdstan = cmdstan_version(),
    pf_chains = if (identical(cell$algorithm, "pathfinder")) pf_chains else NA_integer_,
    hmc_iter = if (identical(cell$algorithm, "sampling")) hmc_iter else NA_integer_,
    max_rhat = max_rhat,
    finished_at = as.character(Sys.time())
  )

  # Append immediately so a crash still leaves partial results
  if (file.exists(results_csv)) {
    write_csv(row, results_csv, append = TRUE)
  } else {
    write_csv(row, results_csv)
  }

  message(
    tag, " -> ", status,
    if (!is.na(fail_class)) paste0(" [", fail_class, "]") else "",
    " (", round(elapsed, 1), "s)"
  )
  if (!is.na(err)) message("  error: ", err)

  invisible(row)
}

results <- lapply(seq_len(nrow(grid)), function(i) run_one(grid[i, ]))
results <- bind_rows(results)

# Rewrite clean full table (append may have duplicated if re-run)
write_csv(results, results_csv)

message("\n=== SUMMARY ===")
print(
  results |>
    count(re, sparse, algorithm, opencl, status, fail_class) |>
    arrange(status, fail_class)
)
message("Results: ", results_csv)
message("OK cells: ", sum(results$status == "ok"), " / ", nrow(results))
