#!/usr/bin/env Rscript
# Minimal accel pilot: one gene, load V1_nk data, fit split-RE + simple shape.
# Speed knobs under test: split RE, Pathfinder (sccomp-style), threads, OpenCL.
#
#   module load R/4.5.3-gfbf-2025b
#   cd /home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC
#   Rscript R/HCC_opencl_one_gene_pilot.R

library(brms)
library(cmdstanr)
library(targets)
library(dplyr)
library(stringr)
library(readr)

options(brms.backend = "cmdstanr")
# Ensure a previous OpenCL session option does not stick unless fit_gpu sets it
options(brms.opencl = NULL)
cmdstanr::check_cmdstan_toolchain()

# ---------------------------------------------------------------------------
# Paths / gene
# ---------------------------------------------------------------------------
case_dir <- "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC"
model_store_root <- "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE"
store <- file.path(model_store_root, "V1_nk", "_targets")
gene <- "ZFP36"
out_dir <- file.path(case_dir, "data", "processed", "re_split_accel_pilot")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Formulas (RE split + simple shape)
# ---------------------------------------------------------------------------
counts_formula <- paste(
  "counts ~ 1 + offset(offset) + age_decade * sex + disease_groups_altered",
  "+ ethnicity_groups + assay_groups_altered + (1 | dataset_id_altered)",
  "+ (1 + age_decade * sex | tissue_groups)",
  "+ (0 + ethnicity_groups | tissue_groups)"
)
shape_formula <- "shape ~ 1 + disease_groups_altered + assay_groups_altered + (1 | tissue_groups)"

form <- bf(
  as.formula(counts_formula),
  as.formula(shape_formula),
  # sparse = TRUE broken with CmdStan >= 2.33 (brms old array syntax)
  sparse = FALSE
)

form_sparse <- bf(
  as.formula(counts_formula),
  as.formula(shape_formula),
  # sparse = TRUE broken with CmdStan >= 2.33 (brms old array syntax)
  sparse = TRUE
)

# ---------------------------------------------------------------------------
# Load one gene's data from V1_nk production store
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

message("Loading data: ", gene, " (", feature_id, ") -> ", target_name)
fit_prod <- tar_read_raw(target_name[[1]], store = store)$brms_fit[[1]]
dat <- droplevels(fit_prod$data)
colnames(dat) <- str_replace_all(colnames(dat), "_+", "_")
message("n = ", nrow(dat), ", tissues = ", n_distinct(dat$tissue_groups))

# ---------------------------------------------------------------------------
# Priors
# prior        = production-style student_t (dense design)
# prior_sparse = same idea, but mean Intercept is class b under sparse=TRUE
# prior_n      = normal (kept for OpenCL / sensitivity experiments)
# ---------------------------------------------------------------------------
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

# sparse=TRUE does not center X → no class Intercept for the mean formula
prior_sparse <- eval(substitute(
  c(
    prior(student_t(3, i, 1.5), class = b, coef = Intercept),
    prior(student_t(3, 0, 1), class = Intercept, dpar = shape),
    prior(student_t(3, 0, 5), class = b),
    prior(student_t(3, 0, 2), class = b, dpar = shape)
  ),
  env = list(i = i)
))

prior_n <- c(
  prior(normal(0, 5), class = Intercept),
  prior(normal(0, 2), class = Intercept, dpar = shape),
  prior(normal(0, 5), class = b),
  prior(normal(0, 2), class = b, dpar = shape),
  prior(normal(0, 2), class = sd)
)

# Force a *new* Stan executable for CPU (do not reuse an OpenCL-compiled one).
force_cpu_stanvars <- stanvar(
  scode = "// force_cpu_recompile_20260721",
  block = "functions"
)

out_rds_pathfinder <- file.path(
  out_dir,
  paste0("fit_pathfinder_cpu_", gene, "_re_split_simple.rds")
)
out_rds_gpu <- file.path(
  out_dir,
  paste0("fit_pathfinder_opencl_", gene, "_re_split_simple.rds")
)

# ---------------------------------------------------------------------------
# Fit 1 — Pathfinder (sccomp-style), CPU; production student_t prior
# ---------------------------------------------------------------------------
fit_pathfinder <- brm(
  formula = form,
  data = dat,
  family = zero_inflated_negbinomial(),
  prior = prior,
  stanvars = force_cpu_stanvars,
  backend = "cmdstanr",
  algorithm = "pathfinder",
  opencl = NULL,
  # brms: num_paths = chains  (sccomp uses num_paths = 50)
  chains = 50,
  cores = 1,
  threads = threading(8),
  history_size = 100,
  max_lbfgs_iters = 100,
  single_path_draws = 80,
  psis_resample = FALSE,
  seed = 20260720,
  refresh = 50,
  silent = 0
)

fit_pathfinder_sparse <- brm(
  formula = form_sparse,
  data = dat,
  family = zero_inflated_negbinomial(),
  prior = prior_sparse,
  stanvars = force_cpu_stanvars,
  backend = "cmdstanr",
  algorithm = "pathfinder",
  opencl = NULL,
  # brms: num_paths = chains  (sccomp uses num_paths = 50)
  chains = 50,
  cores = 1,
  # ! Cannot use threading and sparse matrices at the same time.
  # threads = threading(8),
  history_size = 100,
  max_lbfgs_iters = 100,
  single_path_draws = 80,
  psis_resample = FALSE,
  seed = 20260720,
  refresh = 50,
  silent = 0
)

print(fit_pathfinder)
message("Saved: ", out_rds_pathfinder)

fit_pathfinder <- brm(
  formula = form,
  data = dat,
  family = zero_inflated_negbinomial(),
  prior = prior,
  stanvars = force_cpu_stanvars,
  backend = "cmdstanr",
  algorithm = "pathfinder",
  opencl = NULL,
  # brms: num_paths = chains  (sccomp uses num_paths = 50)
  chains = 50,
  cores = 1,
  threads = threading(8),
  history_size = 100,
  max_lbfgs_iters = 100,
  single_path_draws = 80,
  psis_resample = FALSE,
  seed = 20260720,
  refresh = 50,
  silent = 0
)

# ---------------------------------------------------------------------------
# Fit 2 — Pathfinder + OpenCL (GPU); same formula / prior
# Note: may fail on student_t_lpdf(OpenCL) with this hierarchical ZINB model.
# ---------------------------------------------------------------------------
fit_gpu <- brm(
  formula = form,
  data = dat,
  family = zero_inflated_negbinomial(),
  prior = prior,
  backend = "cmdstanr",
  algorithm = "pathfinder",
  opencl = brms::opencl(c(0L, 0L)),
  chains = 50,
  cores = 1,
  threads = threading(8),
  history_size = 100,
  max_lbfgs_iters = 100,
  single_path_draws = 80,
  psis_resample = FALSE,
  seed = 20260720,
  refresh = 50,
  silent = 0
)

fit_gpu <- brm(
  formula = form,
  data = dat,
  family = zero_inflated_negbinomial(),
  prior = prior_n,
  backend = "cmdstanr",
  algorithm = "pathfinder",
  opencl = brms::opencl(c(0L, 0L)),
  chains = 50,
  cores = 1,
  threads = threading(8),
  history_size = 100,
  max_lbfgs_iters = 100,
  single_path_draws = 80,
  psis_resample = FALSE,
  seed = 20260720,
  refresh = 50,
  silent = 0
)

print(fit_gpu)
message("Saved: ", out_rds_gpu)
