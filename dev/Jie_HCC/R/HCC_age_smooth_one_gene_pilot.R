#!/usr/bin/env Rscript
# One-gene feasibility pilot: continuous-age thin-plate smooth by tissue×sex.
#
# counts ~ offset + sex + s(age_days_scaled, by = tissue_sex, ...) +
#          ethnicity + assay + (1|dataset) + uncorrelated tissue REs
#
#   module load R/4.5.3-gfbf-2025b
#   cd /home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC
#   Rscript R/HCC_age_smooth_one_gene_pilot.R

library(brms)
library(cmdstanr)
library(targets)
library(dplyr)
library(stringr)
library(readr)
library(SummarizedExperiment)

options(brms.backend = "cmdstanr")
options(brms.opencl = NULL)
cmdstanr::check_cmdstan_toolchain()

# ---------------------------------------------------------------------------
# Paths / gene
# ---------------------------------------------------------------------------
case_dir <- "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC"
model_store_root <- "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE"
store <- file.path(model_store_root, "V1_nk", "_targets")
gene <- "ZFP36"
out_dir <- file.path(case_dir, "data", "processed", "age_smooth_pilot")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load one gene's data from V1_nk (fit$data + SE colData for age / assay)
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
  pull(.data$name)

message("Loading data: ", gene, " (", feature_id, ") -> ", target_name)
fit_prod <- tar_read_raw(target_name[[1]], store = store)$brms_fit[[1]]
dat <- droplevels(fit_prod$data)
colnames(dat) <- str_replace_all(colnames(dat), "_+", "_")

se <- tar_read(pseudobulk_sample, store = store)
cd <- as.data.frame(SummarizedExperiment::colData(se))
stopifnot(nrow(cd) == nrow(dat))

# Prefer SE age_days_scaled if present; otherwise scale age_days.
if ("age_days_scaled" %in% names(cd)) {
  dat$age_days_scaled <- as.numeric(cd$age_days_scaled)
} else {
  dat$age_days <- as.numeric(as.character(cd$age_days))
  dat$age_days_scaled <- as.numeric(scale(dat$age_days))
}

# User formula uses assay_groups / dataset_id (not *_altered)
dat$assay_groups <- factor(cd$assay_groups)
dat$dataset_id <- factor(cd$dataset_id)

# tissue × sex interaction for smooth `by`
dat$tissue_sex <- interaction(dat$tissue_groups, dat$sex, drop = TRUE, sep = "_by_")

message(
  "n = ", nrow(dat),
  ", tissues = ", n_distinct(dat$tissue_groups),
  ", tissue_sex levels = ", nlevels(dat$tissue_sex)
)
message(
  "age_days_scaled range: [",
  round(min(dat$age_days_scaled, na.rm = TRUE), 3), ", ",
  round(max(dat$age_days_scaled, na.rm = TRUE), 3), "]"
)

# ---------------------------------------------------------------------------
# Formulas
# ---------------------------------------------------------------------------
counts_formula <- counts ~
  offset(offset) +
  sex +
  s(
    age_days_scaled,
    by = tissue_sex,
    bs = "tp",
    k = 5,
    m = 1,
    id = "age_by_tissue_sex"
  ) +
  ethnicity_groups +
  assay_groups +
  (1 | dataset_id) +
  (1 + sex || tissue_groups) +
  (0 + ethnicity_groups || tissue_groups)

# Keep a light shape formula for ZINB feasibility (not part of the mean test)
shape_formula <- shape ~ 1 + assay_groups + (1 | tissue_groups)

form <- bf(counts_formula, shape_formula)

# ---------------------------------------------------------------------------
# Priors (production-style student_t)
# ---------------------------------------------------------------------------
i <- mean(log1p(dat$counts / exp(dat$offset)), na.rm = TRUE)
prior <- eval(substitute(
  c(
    prior(student_t(3, i, 1.5), class = Intercept),
    prior(student_t(3, 0, 1), class = Intercept, dpar = shape),
    prior(student_t(3, 0, 5), class = b),
    prior(student_t(3, 0, 2), class = b, dpar = shape)
  ),
  env = list(i = i)
))

out_rds <- file.path(out_dir, paste0("fit_age_smooth_hmc_", gene, ".rds"))

# ---------------------------------------------------------------------------
# Fit — HMC (CPU)
# ---------------------------------------------------------------------------
message("Fitting HMC (CPU) -> ", out_rds)
t0 <- Sys.time()

fit <- brm(
  formula = form,
  data = dat,
  family = zero_inflated_negbinomial(),
  prior = prior,
  backend = "cmdstanr",
  # algorithm = "sampling",
  opencl = NULL,
  chains = 4,
  cores = 8,
  threads = threading(2),
  warmup = 400,
  iter = 900,
  seed = 20260720,
  refresh = 50,
  silent = 0,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  file = out_rds,
  file_refit = "always"
)

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
message("Done in ", round(elapsed, 1), "s. Saved: ", out_rds)
print(fit)
