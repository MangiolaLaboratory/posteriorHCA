#!/usr/bin/env Rscript
# Re-test dense + OpenCL combos under use_R_env (OCL_ICD_VENDORS + conda R_env).
# Sparse cells omitted — they fail at Stan compile before OpenCL matters.
#
#   source ~/.bashrc && use_R_env
#   cd /home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC
#   Rscript R/HCC_accel_gpu_retest.R 2>&1 | tee data/processed/re_split_accel_pilot/gpu_retest.log

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

message("OCL_ICD_VENDORS=", Sys.getenv("OCL_ICD_VENDORS"))
message("CONDA_PREFIX=", Sys.getenv("CONDA_PREFIX"))
message("brms=", as.character(packageVersion("brms")),
        " cmdstan=", cmdstan_version())

case_dir <- "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC"
model_store_root <- "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE"
store <- file.path(model_store_root, "V1_nk", "_targets")
gene <- "ZFP36"
out_dir <- file.path(case_dir, "data", "processed", "re_split_accel_pilot")
fits_dir <- file.path(out_dir, "combo_fits_gpu")
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
results_csv <- file.path(out_dir, "accel_gpu_retest_results.csv")
if (file.exists(results_csv)) file.remove(results_csv)

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

model_meta <- tar_meta(starts_with("estimates_chunk_"), store = store) |>
  filter(.data$type == "branch") |>
  mutate(.feature = str_extract(.data$warnings, "(?<=Gene:___)ENSG\\d+(?=___)"))
biomarkers <- read_csv(
  file.path(case_dir, "data", "times_biomarkers.csv"),
  show_col_types = FALSE
) |>
  filter(.data$gene == !!gene, .data$.feature %in% model_meta$.feature)
feature_id <- biomarkers$.feature[[1]]
target_name <- model_meta |>
  filter(.data$.feature == feature_id) |>
  slice(1) |>
  pull(.data$name)

message("Loading data: ", gene)
fit_prod <- tar_read_raw(target_name[[1]], store = store)$brms_fit[[1]]
dat <- droplevels(fit_prod$data)
colnames(dat) <- str_replace_all(colnames(dat), "_+", "_")

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
# normal priors — optional sensitivity (OpenCL student_t failed before)
prior_n <- c(
  prior(normal(0, 5), class = Intercept),
  prior(normal(0, 2), class = Intercept, dpar = shape),
  prior(normal(0, 5), class = b),
  prior(normal(0, 2), class = b, dpar = shape),
  prior(normal(0, 2), class = sd)
)

grid <- tidyr::expand_grid(
  re = c("baseline", "split"),
  algorithm = c("pathfinder", "sampling"),
  prior_set = c("student_t", "normal")
) |>
  mutate(
    cell_id = sprintf(
      "gpu_%s_dense_%s_%s",
      re, algorithm, prior_set
    )
  )

message("=== GPU retest under use_R_env: ", nrow(grid), " cells ===")

run_one <- function(cell) {
  tag <- cell$cell_id
  message("\n========== ", tag, " ==========")
  counts <- if (identical(cell$re, "split")) counts_split else counts_baseline
  form <- bf(as.formula(counts), as.formula(shape_formula), sparse = FALSE)
  pr <- if (identical(cell$prior_set, "normal")) prior_n else prior_dense
  out_rds <- file.path(fits_dir, paste0("fit_", tag, ".rds"))
  sv <- stanvar(scode = paste0("// ", tag, "_use_R_env"), block = "functions")

  args <- list(
    formula = form,
    data = dat,
    family = zero_inflated_negbinomial(),
    prior = pr,
    stanvars = sv,
    backend = "cmdstanr",
    algorithm = cell$algorithm,
    opencl = brms::opencl(c(0L, 0L)),
    seed = 20260720L,
    refresh = 0,
    silent = 2,
    file = out_rds,
    file_refit = "always",
    cores = 1L
  )
  if (identical(cell$algorithm, "pathfinder")) {
    args$chains <- 20L
    args$history_size <- 100L
    args$max_lbfgs_iters <- 100L
    args$single_path_draws <- 40L
    args$psis_resample <- FALSE
    args$threads <- threading(8)
  } else {
    args$chains <- 1L
    args$warmup <- 50L
    args$iter <- 100L
    args$threads <- threading(4)
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
  if (!is.na(err) && nchar(err) > 400) err <- paste0(substr(err, 1, 400), "...")

  fail_class <- NA_character_
  if (identical(status, "error")) {
    fail_class <- dplyr::case_when(
      grepl("CL_PLATFORM_NOT_FOUND|opencl_context", err) ~ "opencl_platform_unavailable",
      grepl("OpenCL|opencl|lpdf\\(OpenCL\\)", err) ~ "opencl_runtime",
      grepl("Fitting failed", err) ~ "opencl_or_fit_failed",
      TRUE ~ "other"
    )
  }

  row <- tibble(
    cell_id = tag,
    gene = gene,
    re = cell$re,
    sparse = FALSE,
    algorithm = cell$algorithm,
    opencl = TRUE,
    prior_set = cell$prior_set,
    status = status,
    fail_class = fail_class,
    elapsed_sec = elapsed,
    error = err,
    ocl_icd_vendors = Sys.getenv("OCL_ICD_VENDORS"),
    conda_prefix = Sys.getenv("CONDA_PREFIX"),
    brms = as.character(packageVersion("brms")),
    cmdstan = cmdstan_version(),
    finished_at = as.character(Sys.time())
  )
  if (file.exists(results_csv)) write_csv(row, results_csv, append = TRUE) else write_csv(row, results_csv)
  message(tag, " -> ", status, " (", round(elapsed, 1), "s)")
  if (!is.na(err)) message("  ", err)
  invisible(row)
}

results <- bind_rows(lapply(seq_len(nrow(grid)), function(i) run_one(grid[i, ])))
write_csv(results, results_csv)
message("\n=== GPU RETEST SUMMARY ===")
print(as.data.frame(results %>% select(cell_id, status, fail_class, elapsed_sec)))
message("Wrote ", results_csv)
