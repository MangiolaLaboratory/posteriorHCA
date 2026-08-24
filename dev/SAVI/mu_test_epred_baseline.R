# Regenerate the fitted NB *mean parameter* mu for one gene, marginalised over
# the universal healthy-atlas covariate grid.
#
# IMPORTANT: the brms family is zero-inflated negative binomial, so
#   posterior_epred = E[Y] = (1 - zi) * mu     (the expected COUNT)
#   posterior_linpred(transform = FALSE) = log(mu) = b0 + b1*age_decade + ...
# We want mu (the NB location), i.e. exp(linear predictor) = the user's
# "mu = b0 + b1*age_decade + ...". So we use posterior_linpred, NOT epred.
#
# Cached to a writable path for mu_test.R to consume (log_mu draws).

suppressMessages({
  library(brms)
  library(tibble)
  library(tidyr)
  library(dplyr)
})

if (!requireNamespace("poco", quietly = TRUE)) {
  pkgload::load_all("/home/a1237163/lab/chen/poco", quiet = TRUE, export_all = FALSE)
}

data_dir   <- "/home/a1237163/lab/chen/posteriorHCA/dev/SAVI/data"
gene_id    <- "ENSG00000196428" # TSC22D2
out_rds    <- "/home/a1237163/mu_test_epred_baseline.rds"

universal_baseline_new_level_tag <- "__universal_baseline__"

newdata_universal_baseline <- function() {
  tibble::tibble(
    age_decade = NA_character_, sex = NA_character_,
    ethnicity_groups = NA_character_, assay_groups_altered = NA_character_,
    disease_groups_altered = "Normal", tissue_groups = "blood", offset = 0
  )
}

fit_levels <- function(fit, var) {
  d <- fit$data
  if (is.null(d) || !var %in% names(d)) return(character(0))
  v <- d[[var]]
  lv <- if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
  lv[!is.na(lv) & nzchar(lv)]
}

build_newdata_grid <- function(fit, newdata,
                               marginalize_over = c("sex", "age_decade", "ethnicity_groups", "assay_groups_altered"),
                               required = c("tissue_groups", "offset", "disease_groups_altered"),
                               new_level_cols = "dataset_id_altered",
                               new_level_tag = universal_baseline_new_level_tag) {
  newdata <- tibble::as_tibble(newdata)
  is_blank <- function(v) length(v) == 0L || is.na(v[[1L]]) ||
    (is.character(v[[1L]]) && !nzchar(v[[1L]]))
  fit_vars <- colnames(fit$data)
  to_expand <- character(0)
  for (col in marginalize_over) {
    if (!col %in% fit_vars) next
    if (!col %in% names(newdata) || is_blank(newdata[[col]])) to_expand <- c(to_expand, col)
  }
  if (length(to_expand)) {
    levels_list <- lapply(to_expand, function(col) fit_levels(fit, col))
    names(levels_list) <- to_expand
    keep <- vapply(levels_list, length, integer(1L)) > 0L
    levels_list <- levels_list[keep]
    to_expand <- names(levels_list)
  }
  if (length(to_expand)) {
    grid <- do.call(tidyr::expand_grid, levels_list)
    rest <- newdata[, setdiff(names(newdata), to_expand), drop = FALSE]
    out <- dplyr::bind_cols(rest[rep(1L, nrow(grid)), , drop = FALSE], grid)
  } else {
    out <- newdata
  }
  for (col in new_level_cols) if (col %in% fit_vars) out[[col]] <- new_level_tag
  tibble::as_tibble(out)
}

# --- reconstruct fit + epred over the marginalised grid ----------------------
sting_estimates_chunk <- readRDS(file.path(data_dir, "SAVI_estimates_chunk.rds"))
row <- sting_estimates_chunk[sting_estimates_chunk$.feature == gene_id, ]
stopifnot(nrow(row) == 1L)

fit <- poco::reconstruct_brmsfit(row$brms_fit_compressed[[1]])

nd_grid <- build_newdata_grid(fit, newdata_universal_baseline())
cat("grid rows:", nrow(nd_grid), "\n")

# Linear predictor of the mu parameter (log scale) = log(mu). offset = 0 in the
# baseline grid, so this is mu at the reference scale.
log_mu_mat <- brms::posterior_linpred(
  fit, newdata = nd_grid, transform = FALSE, re_formula = NULL,
  allow_new_levels = TRUE, sample_new_levels = "gaussian"
)
log_mu <- as.numeric(log_mu_mat) # pool draws x grid (full predictive spread)
mu     <- exp(log_mu)

# Lever (a): per posterior draw, average log(mu) over the marginalised covariate
# grid -> draws of the healthy GRAND MEAN. Its SD is the uncertainty in the
# *mean* (credible interval), NOT the between-condition spread (prediction
# interval). log_mu_mat is [draws x grid], so rowMeans collapses the grid.
grand_mean_draws <- rowMeans(log_mu_mat)

cat("log(mu) pooled: n =", length(log_mu),
    " mean =", mean(log_mu), " sd =", sd(log_mu), "\n")
cat("grand-mean draws: n =", length(grand_mean_draws),
    " mean =", mean(grand_mean_draws), " sd =", sd(grand_mean_draws), "\n")
cat("mu (natural scale): median =", median(mu), " mean =", mean(mu), "\n")
print(quantile(log_mu, c(.5, .9, .99, 1)))

saveRDS(
  list(gene_id = gene_id,
       log_mu           = log_mu,
       mu               = mu,
       mean_log_mu      = mean(log_mu),
       sd_log_mu        = sd(log_mu),
       median_mu        = median(mu),
       n_draws          = length(log_mu),
       grand_mean_draws = grand_mean_draws,
       mean_grand_mean  = mean(grand_mean_draws),
       sd_grand_mean    = sd(grand_mean_draws),
       n_grand          = length(grand_mean_draws)),
  out_rds
)
cat("Saved ->", out_rds, "\n")
