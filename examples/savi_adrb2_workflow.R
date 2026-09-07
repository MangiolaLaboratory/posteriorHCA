#!/usr/bin/env Rscript
# SAVI case study: ADRB2 (ENSG00000169252) in disease-associated monocytes
#
# user data -> TMM offsets -> model.matrix -> estimate_logmu_ql
#   -> load_expression_fit -> expression_draws
#   -> summarize_posterior_draws -> welch_test_means

suppressPackageStartupMessages({
  library(cli)
  library(Seurat)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(brms)
  library(purrr)
})

pkg_dir <- if (file.exists("DESCRIPTION")) {
  normalizePath(".")
} else if (file.exists("../DESCRIPTION")) {
  normalizePath("..")
} else {
  stop("Run this script from the posteriorHCA package root (or its parent).")
}
devtools::load_all(pkg_dir)

cli::cli_h1("posteriorHCA: SAVI ADRB2 Workflow")

cell_type <- "monocytic"
gene_ensg <- "ENSG00000169252"

# ------------------------------------------------------------------------------
# 1. Load counts and map symbols to ENSG (outside posteriorHCA)
# ------------------------------------------------------------------------------
cli::cli_h2("1. Preparing counts")

savi_path <- Sys.getenv("SAVI_PSEUDOBULK_RDS", unset = "")
if (nzchar(savi_path) && file.exists(savi_path)) {
  savi <- readRDS(savi_path)
  savi_mono <- subset(savi, subset = CellType == "17. Disease-associated monocytes")
} else {
  data(savi_mono, package = "posteriorHCA", envir = environment())
}

user_counts <- as.matrix(Seurat::GetAssayData(savi_mono, layer = "counts"))
mapped <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = rownames(user_counts),
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first"
)
keep <- !is.na(mapped) & !duplicated(mapped)
user_counts <- user_counts[keep, , drop = FALSE]
rownames(user_counts) <- unname(mapped[keep])

sample_metadata <- data.frame(
  Category = factor(savi_mono$Category[colnames(user_counts)]),
  row.names = colnames(user_counts),
  stringsAsFactors = FALSE
)

cli::cli_alert_info(
  "User matrix: {nrow(user_counts)} genes x {ncol(user_counts)} samples."
)

# ------------------------------------------------------------------------------
# 2. Scale to HCA reference (offsets only; estimation uses user libraries)
# ------------------------------------------------------------------------------
cli::cli_h2("2. Aligning to HCA monocytic reference")

reference <- load_reference_sample(cell_type)
combined_counts <- merge_with_reference_sample(
  user_counts,
  reference = reference$counts,
  reference_name = reference$sample_id
)
scaling <- calculate_tmm_offset(
  combined_counts,
  reference_name = reference$sample_id
)
user_offset <- scaling$offset[colnames(user_counts)]

cli::cli_alert_info(
  "Reference `{reference$sample_id}` offset = {scaling$offset[[reference$sample_id]]}."
)

# ------------------------------------------------------------------------------
# 3. Estimate cohort log(mu)
# ------------------------------------------------------------------------------
cli::cli_h2("3. Estimating cohort log(mu)")

design_matrix <- model.matrix(~ 0 + Category, data = sample_metadata)
colnames(design_matrix) <- sub("^Category", "", colnames(design_matrix))

expression_estimates <- estimate_logmu_ql(
  user_counts,
  user_offset,
  design_matrix
)
expression_estimates <- expression_estimates[
  expression_estimates$gene == gene_ensg,
  ,
  drop = FALSE
]
print(expression_estimates)

cli::cli_h3("Intercept-only design (~ 1)")

expression_estimates_by_level <- map_dfr(
  levels(sample_metadata$Category),
  function(category) {
    sample_ids <- rownames(sample_metadata)[sample_metadata$Category == category]
    design_one <- model.matrix(
      ~ 1,
      data = sample_metadata[sample_ids, , drop = FALSE]
    )
    out <- estimate_logmu_ql(
      user_counts[, sample_ids, drop = FALSE],
      user_offset[sample_ids],
      design_one
    )
    out <- out[out$gene == gene_ensg, , drop = FALSE]
    out$group <- category
    out
  }
)
print(expression_estimates_by_level)

# ------------------------------------------------------------------------------
# 4. Healthy HCA baseline draws
# ------------------------------------------------------------------------------
cli::cli_h2("4. Healthy HCA baseline draws")

expression_fit <- load_expression_fit(
  cell_type = cell_type,
  gene_ensg = gene_ensg
)
newdata <- build_newdata_grid(
  expression_fit,
  disease_groups = "Normal",
  tissue_groups = "blood",
  assay_groups = "10x Genomics 3"
)
posterior_draws <- expression_draws(
  expression_fit,
  newdata = newdata,
  quantity = "linpred",
  collapse = "mean"
)

cli::cli_alert_info(
  "HCA mean log(mu) = {round(mean(posterior_draws$draws), 3)}, SD = {round(sd(posterior_draws$draws), 3)}."
)

# ------------------------------------------------------------------------------
# 5. Welch test
# ------------------------------------------------------------------------------
cli::cli_h2("5. Welch test vs healthy baseline")

cohort_estimate <- expression_estimates[
  expression_estimates$group == "SAVI",
  ,
  drop = FALSE
]
posterior_summary <- summarize_posterior_draws(
  posterior_draws,
  value = cohort_estimate$log_mu
)
print(welch_test_means(
  cohort_estimate$log_mu,
  cohort_estimate$se,
  posterior_summary$mean,
  posterior_summary$sd,
  n1 = cohort_estimate$n,
  n2 = posterior_summary$n
))

test_results <- map_dfr(
  expression_estimates$group,
  function(group) {
    cohort_estimate <- expression_estimates[
      expression_estimates$group == group,
      ,
      drop = FALSE
    ]
    posterior_summary <- summarize_posterior_draws(
      posterior_draws,
      value = cohort_estimate$log_mu
    )
    test <- welch_test_means(
      cohort_estimate$log_mu,
      cohort_estimate$se,
      posterior_summary$mean,
      posterior_summary$sd,
      n1 = cohort_estimate$n,
      n2 = posterior_summary$n
    )
    data.frame(
      group = group,
      log_mu = test$mu1,
      se = test$se1,
      p_value = test$p_value,
      stringsAsFactors = FALSE
    )
  }
)
print(test_results)

# ------------------------------------------------------------------------------
# 6. Plots
# ------------------------------------------------------------------------------
cli::cli_h2("6. Plots")

print(plot_hca_draws(
  draws = posterior_draws,
  subtitle = "Normal, 10x Genomics 3 healthy baseline"
))
print(plot_cohort_vs_hca(
  hca_draws = posterior_draws,
  cohort_est = test_results,
  subtitle = "QL cohort estimates",
  annotate = c("group", "p_value")
))
