#!/usr/bin/env Rscript
# ==============================================================================
# SAVI case study: ADRB2 in disease-associated monocytes (PBMC / blood)
#
# Reference script for vignettes/cohort-expression-workflow.Rmd
#
# Pipeline:
#   1. Pseudobulk counts -> harmonise gene ids (ENSG)
#   2. Align to one HCA reference library (monocytic)
#        - core: merge_with_reference_sample + calculate_tmm_offset
#        - wrapper: scale_to_hca_reference
#   3. edgeR QL cohort log(mu) for ADRB2
#        - core: design_from_formula + estimate_logmu_ql
#        - wrapper: estimate_cohort_logmu (filter genes)#   4. Healthy HCA posterior (Normal, blood, 10x Genomics 3)
#   5. Welch test vs HCA
#   6. Plots
# ==============================================================================

suppressPackageStartupMessages({
  library(cli)
  library(Seurat)
  library(dplyr)
  library(brms)
})

pkg_dir <- "/home/a1237163/lab/chen/posteriorHCA"
devtools::load_all(pkg_dir)

cli::cli_h1("posteriorHCA: SAVI ADRB2 Workflow")

# ------------------------------------------------------------------------------
# 1. Load SAVI pseudobulk and harmonise gene ids to ENSG
# ------------------------------------------------------------------------------
cli::cli_h2("1. Preparing SAVI Disease-Associated Monocyte Counts")

savi_path <- Sys.getenv(
  "SAVI_PSEUDOBULK_RDS",
  unset = "/home/a1237163/lab/chen/posteriorHCA_case_studies/SAVI/data/GSE226598_SAVI_pseudobulk_Sample_CellType.rds"
)

if (file.exists(savi_path)) {
  savi <- readRDS(savi_path)
  savi_mono <- subset(savi, subset = CellType == "17. Disease-associated monocytes")
} else {
  data(savi_mono, package = "posteriorHCA", envir = environment())
}

savi_mono <- harmonise_gene_ids(savi_mono, id_type = "auto")
cli::cli_alert_info(
  "Harmonised Seurat object: {nrow(savi_mono)} genes x {ncol(savi_mono)} samples."
)

# ------------------------------------------------------------------------------
# 2. Align to HCA monocytic reference
# ------------------------------------------------------------------------------
cli::cli_h2("2. Aligning to HCA Monocytic Reference Sample")

cell_type <- "monocytic"

cli::cli_h3("2a. Core functions (matrix)")

ref <- load_reference_sample(cell_type)
user_mat <- as.matrix(Seurat::GetAssayData(savi_mono, layer = "counts"))

combined <- merge_with_reference_sample(
  user_mat,
  reference = ref$counts,
  reference_name = ref$sample_id
)
scaling <- calculate_tmm_offset(
  combined,
  reference_name = ref$sample_id,
  method = "TMMwsp"
)

cli::cli_alert_info(
  "Merged matrix: {nrow(combined)} genes x {ncol(combined)} samples."
)
cli::cli_alert_info(
  "Reference `{ref$sample_id}` offset = {scaling$offset[[ref$sample_id]]}."
)

cli::cli_h3("2b. Wrapper scale_to_hca_reference()")

aligned <- scale_to_hca_reference(savi_mono, cell_type)
print(table(aligned$sample_role))

stopifnot(all.equal(
  unname(scaling$offset[colnames(user_mat)]),
  unname(aligned$hca_offset[colnames(user_mat)]),
  tolerance = 1e-10
))
cli::cli_alert_success("Core and wrapper offsets match.")

# ------------------------------------------------------------------------------
# 3. Estimate cohort log(mu) for ADRB2
#
# Absolute log(μ) = QL coefficients with prior.count = 0 on the TMM offset scale.
# Prefer cell-means designs: ~ 0 + Category.
# Core: design_from_formula + estimate_logmu_ql (all genes)
# Wrapper: estimate_cohort_logmu (filter genes)
# ------------------------------------------------------------------------------
cli::cli_h2("3. Estimating Cohort log(mu) for ADRB2")

cli::cli_h3("3a. Core design_from_formula() + estimate_logmu_ql()")

meta_core <- data.frame(
  Category = c(
    as.character(savi_mono$Category[colnames(user_mat)]),
    "reference"
  ),
  row.names = colnames(combined),
  stringsAsFactors = FALSE
)
meta_core$Category <- factor(meta_core$Category)

design <- design_from_formula(~ 0 + Category, meta_core)
est_all <- estimate_logmu_ql(
  counts = combined,
  offset = scaling$offset,
  design = design,
  cell_type = cell_type
)
gene_ensg <- resolve_gene_one("ADRB2", cell_type = cell_type)
est_core <- est_all[est_all$gene == gene_ensg, , drop = FALSE]
est_core$gene_symbol <- "ADRB2"
print(est_core)

cli::cli_h3("3b. Wrapper estimate_cohort_logmu()")

est <- estimate_cohort_logmu(
  aligned,
  formula = ~ 0 + Category,
  genes = "ADRB2"
)
print(est)

stopifnot(all.equal(
  est_core$log_mu[match(est$group, est_core$group)],
  est$log_mu,
  tolerance = 1e-8
))
cli::cli_alert_success("Core and wrapper log(mu) estimates match.")

# ------------------------------------------------------------------------------
# 3c. Demo: one cohort at a time with intercept-only design (~ 1)
#
# Alternative to ~ 0 + Category on the full object. Subset to a single
# Category, fit NB/QL with formula = ~ 1, and treat the intercept as that
# cohort's absolute log(mu). Useful when you want per-cohort fits without a
# multi-level design matrix.
# ------------------------------------------------------------------------------
cli::cli_h3("3c. Demo: single-cohort ~ 1 (intercept = cohort log(mu))")

est_by_level <- purrr::map_dfr(
  levels(aligned$Category),
  .f = function(x) {
    estimate_cohort_logmu(
      subset(aligned, Category == x),
      formula = ~ 1,
      genes = "ADRB2"
    ) %>%
      dplyr::mutate(group = x)
  }
)
print(est_by_level)

# ------------------------------------------------------------------------------
# 4. Healthy HCA baseline draws
# ------------------------------------------------------------------------------
cli::cli_h2("4. Generating Healthy Baseline Draws (Normal, 10x Genomics 3)")

fit <- load_expr_fit(cell_type = cell_type, gene = "ADRB2")

grid <- build_newdata_grid(
  fit,
  disease_groups = "Normal",
  tissue_groups = "blood",
  assay_groups = "10x Genomics 3"
)

hca_res <- expr_draws(
  fit,
  newdata = grid,
  quantity = "linpred",
  collapse = "mean"
)

cli::cli_alert_info(
  "Healthy HCA baseline for {hca_res$gene_symbol}: mean log(mu) = {round(mean(hca_res$draws), 3)}, SD = {round(sd(hca_res$draws), 3)}."
)

hca_pred <- expr_predict(
  fit = fit,
  disease_groups = "Normal",
  tissue_groups = "blood",
  assay_groups = "10x Genomics 3",
  quantity = "linpred",
  collapse = "mean"
)

# ------------------------------------------------------------------------------
# 5. Welch test vs healthy baseline
# ------------------------------------------------------------------------------
cli::cli_h2("5. Testing Cohorts Against Healthy Baseline (QL)")

test_results <- welch_t_test_cohort_hca(
  cohort_est = est,
  hca_draws = hca_pred
)
print(test_results)

cohort <- cohort_estimate_at(est, group = "SAVI")
baseline <- summarize_posterior_draws(hca_res, value = cohort$mu)
welch_test_means(
  cohort$mu, cohort$se,
  baseline$mean, baseline$sd,
  n1 = cohort$n, n2 = baseline$n
)

# ------------------------------------------------------------------------------
# 6. Plots
# ------------------------------------------------------------------------------
cli::cli_h2("6. Plotting Healthy Baseline and Cohort Comparisons")

print(plot_hca_draws(
  draws = hca_res,
  subtitle = "Normal, 10x Genomics 3 healthy baseline"
))

print(plot_cohort_vs_hca(
  hca_draws = hca_res,
  test_results = test_results,
  subtitle = "QL cohort estimates",
  annotate = c("group", "p_value", "direction")
))

print(plot_hca_draws(
  draws = hca_pred,
  title = "ADRB2 predicted count posterior (healthy HCA)"
))
