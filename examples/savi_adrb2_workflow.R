#!/usr/bin/env Rscript
# ==============================================================================
# SAVI Case Study: ADRB2 Cohort Estimation & Healthy HCA Baseline Testing
# ==============================================================================

suppressPackageStartupMessages({
  library(cli)
  library(edgeR)
  library(SummarizedExperiment)
  library(Seurat)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(brms)
  library(readr)
})

# Load posteriorHCA functions
pkg_dir <- "/home/a1237163/lab/chen/posteriorHCA"
for (f in list.files(file.path(pkg_dir, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

cli::cli_h1("posteriorHCA: SAVI ADRB2 Workflow")

# ------------------------------------------------------------------------------
# 1. Load SAVI pseudobulk Seurat object and prepare user counts
# ------------------------------------------------------------------------------
cli::cli_h2("1. Preparing SAVI Disease-Associated Monocyte Counts")

savi_path <- "/home/a1237163/lab/chen/posteriorHCA_case_studies/SAVI/data/GSE226598_SAVI_pseudobulk_Sample_CellType.rds"
savi <- readRDS(savi_path)
savi_mono <- subset(savi, subset = CellType == "17. Disease-associated monocytes")
counts_sym <- as.matrix(LayerData(savi_mono, assay = "RNA", layer = "counts"))

# Map gene symbols to Ensembl IDs
ensg <- mapIds(org.Hs.eg.db, rownames(counts_sym), "ENSEMBL", "SYMBOL", multiVals = "first")
keep <- !is.na(ensg) & !duplicated(ensg)
user <- counts_sym[keep, , drop = FALSE]
rownames(user) <- unname(ensg[keep])

cli::cli_alert_info("User count matrix: {nrow(user)} genes x {ncol(user)} samples.")

# ------------------------------------------------------------------------------
# 2. Retrieve Atlas Reference Sample and Align
# ------------------------------------------------------------------------------
cli::cli_h2("2. Aligning to HCA Monocytic Reference Sample")

ref <- get_reference_sample_ready("monocytic")
aligned <- scale_to_hca_reference(user, ref)

cli::cli_alert_info("Reference sample ID: {aligned$reference_name}")
cli::cli_alert_info("Reference sample offset: {aligned$offset[[aligned$reference_name]]}")
cli::cli_alert_info("Sample role table:")
print(table(aligned$sample_role))

# ------------------------------------------------------------------------------
# 3. Estimate Cohort log(mu) with edgeR / QL Wald SE
# ------------------------------------------------------------------------------
cli::cli_h2("3. Estimating Cohort log(mu) for ADRB2 (ENSG00000169252)")

group <- savi_mono$Category[match(colnames(user), colnames(savi_mono))]
adrb2_id <- "ENSG00000169252" # ADRB2

est <- estimate_cohort_logmu(
  aligned,
  group = group,
  genes = adrb2_id
)
print(est)

# ------------------------------------------------------------------------------
# 4. Generate Healthy Baseline Draws from Pre-Trained HCA Model
# ------------------------------------------------------------------------------
cli::cli_h2("4. Generating Healthy Baseline Draws (Normal, 10x Genomics 3)")

# Load pre-trained brms model for monocytic ADRB2 from Nectar
fit <- load_expr_fit(cell_type = "monocytic", gene_ensg = adrb2_id)

# Build covariate grid (Normal disease, 10x Genomics 3 assay, marginalised demographics)
grid <- build_newdata_grid(
  fit,
  disease_groups = "Normal",
  assay_groups = "10x Genomics 3",
  tissue_groups = 'liver'
)

# Posterior draws of latent log(mu) (linpred without transformation)
hca_res <- expr_draws(
  fit,
  newdata = grid,
  quantity = "linpred",
  collapse = "mean"
)

cli::cli_alert_info("Healthy HCA baseline: Mean log(mu) = {round(mean(hca_res$draws), 3)}, SD = {round(sd(hca_res$draws), 3)}.")

# ------------------------------------------------------------------------------
# 5. Conduct Welch t-Test for Each Cohort vs. Healthy Baseline
# ------------------------------------------------------------------------------
cli::cli_h2("5. Testing Cohorts Against Healthy Baseline")

test_results_list <- list()
for (grp in setdiff(unique(est$group), "reference")) {
  sub_est <- est[est$group == grp, , drop = FALSE]
  res <- test_cohort_vs_hca(
    cohort_est = sub_est,
    hca_draws = hca_res$draws,
    gene = "ADRB2",
    cohort_name = grp
  )
  test_results_list[[grp]] <- res
}

test_results <- do.call(rbind, test_results_list)
print(test_results)
