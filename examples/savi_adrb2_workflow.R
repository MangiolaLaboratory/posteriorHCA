#!/usr/bin/env Rscript
# ==============================================================================
# SAVI case study: ADRB2 in disease-associated monocytes (PBMC / blood)
#
# Reference script for vignettes/cohort-expression-workflow.Rmd
#
# Pipeline:
#   1. Pseudobulk counts -> harmonise gene ids (ENSG)
#   2. TMM-align to one HCA reference library (monocytic)
#   3. edgeR QL cohort log(mu) for ADRB2
#   4. Healthy HCA posterior (Normal, blood, 10x Genomics 3)
#   5. Welch test vs HCA (wrapper + manual helpers)
#   6. Plots
# ==============================================================================

suppressPackageStartupMessages({
  library(cli)
  library(edgeR)
  library(Seurat)
  library(tidyseurat)
  library(brms)
})

pkg_dir <- "/home/a1237163/lab/chen/posteriorHCA"
devtools::load_all(pkg_dir)

cli::cli_h1("posteriorHCA: SAVI ADRB2 Workflow")

# ------------------------------------------------------------------------------
# 1. Load SAVI pseudobulk Seurat object and harmonise gene ids to ENSG
#
# SAVI (GSE226598): PBMC scRNA-seq pseudobulked by sample x cell type.
# We subset disease-associated monocytes (17 samples: CTRL / SAVI / SAVI_treated).
# harmonise_gene_ids() maps symbols -> ENSG so cohort counts match HCA models.
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
# 2. Retrieve atlas reference sample and align
#
# scale_to_hca_reference() merges one HCA library, runs TMMwsp, and stores
# hca_offset / sample_role on each sample. Reference sits at offset 0.
# ------------------------------------------------------------------------------
cli::cli_h2("2. Aligning to HCA Monocytic Reference Sample")

cell_type <- "monocytic"
aligned <- scale_to_hca_reference(savi_mono, cell_type)

cli::cli_alert_info("Reference sample ID: {aligned$hca_reference_name[1]}")
cli::cli_alert_info("Cell type: {aligned$hca_cell_type[1]}")
cli::cli_alert_info("Reference sample offset: {aligned$hca_offset[aligned$hca_reference_name[1]]}")
cli::cli_alert_info("Sample role table:")
print(table(aligned$sample_role))

# ------------------------------------------------------------------------------
# 3. Estimate cohort log(mu) for ADRB2
#
# Primary fit: one edgeR QL model with ~ 0 + Category (group means + reference).
# Secondary loop: re-estimate each Category level alone (~ 1) for comparison.
# Use formula = ~ 0 + Category for downstream testing (includes reference row).
# ------------------------------------------------------------------------------
cli::cli_h2("3. Estimating Cohort log(mu) for ADRB2")

est <- estimate_cohort_logmu(
  aligned,
  formula = ~ 0 + Category,
  genes = "ADRB2"
)
print(est)

est <-
  map_dfr(
    aligned$Category %>% levels(),
    .f = function(x) {
      estimate_cohort_logmu(aligned %>% filter(Category == x),
                            formula = ~ 1,
                            genes = "ADRB2") %>%
        mutate(group = x)
    }
  )
print(est)

# ------------------------------------------------------------------------------
# 4. Generate healthy baseline draws from pre-trained HCA model
#
# load_expr_fit() caches the brms model from Nectar (monocytic / ADRB2).
# Covariates fixed to healthy blood 10x Genomics 3 (matches SAVI tissue).
# hca_res: posterior linpred draws (log mu scale) via expr_draws().
# hca_pred: same fit via expr_predict() convenience wrapper.
# ------------------------------------------------------------------------------
cli::cli_h2("4. Generating Healthy Baseline Draws (Normal, 10x Genomics 3)")

fit <- load_expr_fit(cell_type = cell_type, gene = "ADRB2")

grid <- build_newdata_grid(
  fit,
  disease_groups = "Normal",
  tissue_groups = 'blood',
  assay_groups = "10x Genomics 3"
)

hca_res <- expr_draws(
  fit,
  newdata = grid,
  quantity = "linpred",
  collapse = "mean"
)

cli::cli_alert_info(
  "Healthy HCA baseline for {hca_res$gene_symbol} ({hca_res$gene_ensg}): mean log(mu) = {round(mean(hca_res$draws), 3)}, SD = {round(sd(hca_res$draws), 3)}."
)

hca_pred <- expr_predict(
  fit            = fit,
  disease_groups = "Normal",
  tissue_groups  = 'blood',
  assay_groups   = "10x Genomics 3",
  quantity       = "linpred",
  collapse       = "mean"
)

# ------------------------------------------------------------------------------
# 5. Test cohorts against healthy baseline (QL point estimates)
#
# Wrapper: welch_t_test_cohort_hca() loops cohorts, skips reference by default.
# Manual: extract mu/se with cohort_estimate_at(), summarise HCA draws, then call
# welch_test_means() for the generic heteroscedastic comparison.
# ------------------------------------------------------------------------------
cli::cli_h2("5. Testing Cohorts Against Healthy Baseline (QL)")

test_results <- welch_t_test_cohort_hca(
  cohort_est = est,
  hca_draws = hca_pred
)
print(test_results)

# Manual workflow for a single cohort (SAVI vs HCA)
cohort <- cohort_estimate_at(est, group = "SAVI")
baseline <- summarize_posterior_draws(hca_res, value = cohort$mu)
test <- welch_test_means(
  cohort$mu, cohort$se,
  baseline$mean, baseline$sd,
  n1 = cohort$n, n2 = baseline$n
)

# ------------------------------------------------------------------------------
# 6. Visualise healthy baseline and cohort comparisons
# ------------------------------------------------------------------------------
cli::cli_h2("6. Plotting Healthy Baseline and Cohort Comparisons")

hca_plot <- plot_hca_draws(
  draws = hca_res,
  subtitle = "Normal, 10x Genomics 3 healthy baseline"
)
print(hca_plot)

cohort_plot <- plot_cohort_vs_hca(
  hca_draws = hca_res,
  test_results = test_results,
  subtitle = "QL cohort estimates",
  annotate = c("group", "p_value", "direction")
)
print(cohort_plot)

hca_count_plot <- plot_hca_draws(
  draws = hca_pred,
  title = "ADRB2 predicted count posterior (healthy HCA)"
)
print(hca_count_plot)
