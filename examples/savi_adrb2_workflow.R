#!/usr/bin/env Rscript
# ==============================================================================
# SAVI Case Study: ADRB2 Cohort Estimation & Healthy HCA Baseline Testing
# ==============================================================================

suppressPackageStartupMessages({
  library(cli)
  library(edgeR)
  library(Seurat)
  library(brms)
})

pkg_dir <- "/home/a1237163/lab/chen/posteriorHCA"
devtools::load_all(pkg_dir)

cli::cli_h1("posteriorHCA: SAVI ADRB2 Workflow")

# ------------------------------------------------------------------------------
# 1. Load SAVI pseudobulk Seurat object and harmonise gene ids to ENSG
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

savi_mono <- harmonise_gene_ids(savi_mono, id_type = "symbol")
cli::cli_alert_info(
  "Harmonised Seurat object: {nrow(savi_mono)} genes x {ncol(savi_mono)} samples."
)

# ------------------------------------------------------------------------------
# 2. Retrieve Atlas Reference Sample and Align
# ------------------------------------------------------------------------------
cli::cli_h2("2. Aligning to HCA Monocytic Reference Sample")

cell_type <- "monocytic"
ref <- get_reference_sample_ready(cell_type)
aligned <- scale_to_hca_reference(savi_mono, ref)

cli::cli_alert_info("Reference sample ID: {aligned$hca_reference_name[1]}")
cli::cli_alert_info("Cell type: {aligned$hca_cell_type[1]}")
cli::cli_alert_info("Reference sample offset: {aligned$hca_offset[aligned$hca_reference_name[1]]}")
cli::cli_alert_info("Sample role table:")
print(table(aligned$sample_role))

# ------------------------------------------------------------------------------
# 3. Estimate Cohort log(mu) for ADRB2
# ------------------------------------------------------------------------------
cli::cli_h2("3. Estimating Cohort log(mu) for ADRB2")

est <- estimate_cohort_logmu(
  aligned,
  group = "Category",
  genes = "ADRB2"
)
print(est)

# ------------------------------------------------------------------------------
# 4. Generate Healthy Baseline Draws from Pre-Trained HCA Model
# ------------------------------------------------------------------------------
cli::cli_h2("4. Generating Healthy Baseline Draws (Normal, 10x Genomics 3)")

fit <- load_expr_fit(cell_type = cell_type, gene = "ADRB2")

grid <- build_newdata_grid(
  fit,
  disease_groups = "Normal",
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
  cell_type      = cell_type,
  gene           = "ADRB2",
  disease_groups = "Normal",
  assay_groups   = "10x Genomics 3",
  quantity       = "predict",
  collapse       = "mean"
)

# ------------------------------------------------------------------------------
# 5. Test cohorts against healthy baseline (QL point estimates)
# ------------------------------------------------------------------------------
cli::cli_h2("5. Testing Cohorts Against Healthy Baseline (QL)")

test_results <- welch_t_test_cohort_hca(
  cohort_est = est,
  hca_draws = hca_res
)
print(test_results)

# ------------------------------------------------------------------------------
# 5b. Optional: bootstrap cohort log(mu) and test
# ------------------------------------------------------------------------------
cli::cli_h2("5b. Bootstrap Cohort log(mu) and Test (optional)")

boot_est <- bootstrap_cohort_logmu_batch(
  aligned,
  group = "Category",
  genes = "ADRB2",
  cohort_est = est,
  hca_draws = hca_res,
  n_boot = 2000L,
  seed = 42L
)
print(boot_est)

boot_test <- welch_t_test_cohort_hca(
  cohort_est = boot_est,
  hca_draws = hca_res
)
print(boot_test)

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

boot_plot <- plot_cohort_vs_hca(
  hca_draws = hca_res,
  test_results = boot_test,
  subtitle = "Bootstrap cohort estimates",
  annotate = c("group", "p_value", "direction", "method")
)
print(boot_plot)

hca_count_plot <- plot_hca_draws(
  draws = hca_pred,
  title = "ADRB2 predicted count posterior (healthy HCA)"
)
print(hca_count_plot)
