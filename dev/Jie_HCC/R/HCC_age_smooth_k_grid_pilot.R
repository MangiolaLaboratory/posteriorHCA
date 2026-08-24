#!/usr/bin/env Rscript
# Age-smooth k-grid pilot — sequential fallback.
#
# Preferred: Slurm targets pipeline (4 genes × k ∈ {3,4,5}, healthy Normal only):
#   module load R/4.5.3-gfbf-2025b
#   Rscript targets/run_age_smooth_k_grid.R
#
# This script runs the same grid locally / sequentially (useful for debugging).
# Data filter: disease_groups_altered == "Normal".
#
# Optional env:
#   GENES=ZFP36,VIM   KS=3,5   SKIP_EXISTING=true   SEED=20260720

suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(targets)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(SummarizedExperiment)
})

options(brms.backend = "cmdstanr")
options(brms.opencl = NULL)
cmdstanr::check_cmdstan_toolchain()

case_dir <- "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC"
source(file.path(case_dir, "R", "re_vs_interaction_helpers.R"))
source(file.path(case_dir, "R", "age_smooth_helpers.R"))

model_store_root <- "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE"
store <- file.path(model_store_root, "V1_nk", "_targets")
out_dir <- file.path(case_dir, "data", "processed", "age_smooth_k_grid_pilot")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
results_csv <- file.path(out_dir, "age_smooth_k_grid_results.csv")

genes_raw <- Sys.getenv("GENES", unset = "")
genes <- if (nzchar(genes_raw)) {
  trimws(strsplit(genes_raw, ",", fixed = TRUE)[[1]])
} else {
  pilot_biomarker_genes()
}

ks_raw <- Sys.getenv("KS", unset = "")
ks <- if (nzchar(ks_raw)) {
  as.integer(trimws(strsplit(ks_raw, ",", fixed = TRUE)[[1]]))
} else {
  age_smooth_ks_default()
}

skip_existing <- identical(Sys.getenv("SKIP_EXISTING", unset = "false"), "true")
seed <- as.integer(Sys.getenv("SEED", unset = "20260720"))

grid <- tidyr::expand_grid(gene = genes, k = ks) |>
  mutate(
    out_rds = file.path(out_dir, paste0("fit_age_smooth_hmc_", gene, "_k", k, ".rds"))
  )

message(
  "Age-smooth k-grid (sequential): ", nrow(grid), " cells. ",
  "Prefer: Rscript targets/run_age_smooth_k_grid.R"
)

model_meta <- tar_meta(starts_with("estimates_chunk_"), store = store) |>
  filter(.data$type == "branch") |>
  mutate(.feature = str_extract(.data$warnings, "(?<=Gene:___)ENSG\\d+(?=___)"))

biomarkers <- read_csv(
  file.path(case_dir, "data", "times_biomarkers.csv"),
  show_col_types = FALSE
) |>
  filter(.data$gene %in% genes, .data$.feature %in% model_meta$.feature)

stopifnot(nrow(biomarkers) == length(unique(genes)))

se <- tar_read(pseudobulk_sample, store = store)
cd <- as.data.frame(SummarizedExperiment::colData(se))

load_gene_dat <- function(gene) {
  feature_id <- biomarkers |>
    filter(.data$gene == !!gene) |>
    pull(.data$.feature)
  stopifnot(length(feature_id) == 1L)
  target_name <- model_meta |>
    filter(.data$.feature == feature_id) |>
    slice(1) |>
    pull(.data$name)
  message("Loading data: ", gene, " (", feature_id, ") -> ", target_name)
  fit_prod <- tar_read_raw(target_name[[1]], store = store)$brms_fit[[1]]
  dat <- prepare_age_smooth_data(fit_prod$data, cd)
  message(
    "  n = ", nrow(dat),
    ", tissue_sex = ", nlevels(dat$tissue_sex)
  )
  list(dat = dat, .feature = feature_id)
}

summarise_fit <- function(fit, gene, k, .feature, out_rds, elapsed_sec, status, err = NA_character_) {
  if (identical(status, "ok") && !is.null(fit)) {
    rh <- tryCatch(posterior::rhat(fit), error = function(e) NA_real_)
    rh <- rh[is.finite(rh)]
    np <- tryCatch(brms::nuts_params(fit), error = function(e) NULL)
    pct_div <- if (is.null(np)) NA_real_ else 100 * mean(np$Value[np$Parameter == "divergent__"])
    tibble(
      gene = gene, .feature = .feature, k = as.integer(k), status = status,
      error = err, elapsed_sec = elapsed_sec,
      max_rhat = if (length(rh)) max(rh) else NA_real_,
      n_rhat_gt_1_05 = if (length(rh)) sum(rh > 1.05) else NA_integer_,
      pct_divergent = pct_div,
      out_rds = out_rds
    )
  } else {
    tibble(
      gene = gene, .feature = .feature, k = as.integer(k), status = status,
      error = err, elapsed_sec = elapsed_sec,
      max_rhat = NA_real_, n_rhat_gt_1_05 = NA_integer_,
      pct_divergent = NA_real_, out_rds = out_rds
    )
  }
}

dat_cache <- list()
result_rows <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  row <- grid[i, ]
  message(sprintf("\n=== [%d/%d] %s  k = %d ===", i, nrow(grid), row$gene, row$k))

  if (isTRUE(skip_existing) && file.exists(row$out_rds)) {
    message("Skipping existing: ", row$out_rds)
    fit_existing <- tryCatch(readRDS(row$out_rds), error = function(e) NULL)
    .feature <- biomarkers$.feature[biomarkers$gene == row$gene][[1]]
    result_rows[[i]] <- summarise_fit(
      fit_existing, row$gene, row$k, .feature, row$out_rds,
      elapsed_sec = NA_real_, status = "skipped_existing"
    )
    next
  }

  if (is.null(dat_cache[[row$gene]])) {
    dat_cache[[row$gene]] <- load_gene_dat(row$gene)
  }
  packed <- dat_cache[[row$gene]]

  message("Fitting HMC (CPU, control = NULL) -> ", row$out_rds)
  t0 <- Sys.time()
  fit <- tryCatch(
    fit_age_smooth_gene(
      data = packed$dat,
      gene = row$gene,
      k = as.integer(row$k),
      seed = seed,
      control = NULL
    ),
    error = function(e) e
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (inherits(fit, "error")) {
    message("FAILED in ", round(elapsed, 1), "s: ", conditionMessage(fit))
    result_rows[[i]] <- summarise_fit(
      NULL, row$gene, row$k, packed$.feature, row$out_rds,
      elapsed_sec = elapsed, status = "error", err = conditionMessage(fit)
    )
  } else {
    saveRDS(fit, row$out_rds, compress = "xz")
    message("Done in ", round(elapsed, 1), "s. Saved: ", row$out_rds)
    result_rows[[i]] <- summarise_fit(
      fit, row$gene, row$k, packed$.feature, row$out_rds,
      elapsed_sec = elapsed, status = "ok"
    )
  }

  bind_rows(result_rows) |> write_csv(results_csv)
}

results <- bind_rows(result_rows)
write_csv(results, results_csv)
message("\n=== Grid complete ===\nResults: ", results_csv)
print(results |> select(gene, k, status, elapsed_sec, max_rhat, pct_divergent))
