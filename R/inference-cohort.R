#' Dirichlet weights for Bayesian bootstrap
#'
#' Draws Dirichlet(1, ..., 1) weights scaled so they sum to `n`.
#'
#' @param n Integer number of observations.
#' @return Numeric vector of length `n` summing to `n`.
#' @keywords internal
#' @noRd
draw_dirichlet_weights <- function(n) {
  if (n <= 0L) {
    return(numeric(0))
  }
  if (n == 1L) {
    return(1)
  }
  weights <- stats::rexp(n)
  n * weights / sum(weights)
}

#' Bayesian bootstrap of cohort log(mu)
#'
#' Computes posterior draws of cohort latent log(mu) by repeatedly fitting
#' [edgeR::mglmOneGroup()] with Dirichlet weights.
#'
#' `counts` may be a matrix, `SummarizedExperiment`, `SingleCellExperiment`,
#' `Seurat` object, or an aligned list from [scale_to_hca_reference()].
#'
#' @param counts Count container or matrix.
#' @param offset Numeric vector of sample offsets, or `NULL` if provided in
#'   an aligned `counts` list.
#' @param dispersion Numeric scalar or named vector of negative-binomial
#'   dispersion values. If `NULL`, dispersion is estimated across `counts`.
#' @param group Optional group vector or column name to subset samples. If
#'   `NULL`, all non-reference samples (or all samples) are used.
#' @param gene Character scalar; gene id (must be in `rownames(counts)`).
#' @param n_boot Integer; number of bootstrap iterations (default 2000L).
#' @param seed Optional RNG seed for reproducibility.
#' @param assay Assay name for SE / Seurat input.
#' @return A numeric vector of length `n_boot` containing posterior log(mu) draws.
#' @export
#' @importFrom cli cli_abort
bootstrap_cohort_logmu <- function(
  counts,
  offset = NULL,
  dispersion = NULL,
  group = NULL,
  gene,
  n_boot = 2000L,
  seed = NULL,
  assay = NULL
) {
  if (missing(gene) || is.null(gene) || length(gene) != 1L) {
    cli_abort("`gene` must be a single gene identifier.")
  }

  sample_role <- NULL
  source_obj <- counts
  if (is_aligned_result(counts)) {
    if (is.null(offset)) {
      offset <- counts$offset
    }
    sample_role <- counts$sample_role
    counts <- counts$counts
  } else {
    counts <- extract_count_matrix(source_obj, assay = assay, arg_name = "counts")
  }

  if (!gene %in% rownames(counts)) {
    cli_abort("Gene `{gene}` not found in `counts`.")
  }

  if (is.null(offset)) {
    cli_abort("`offset` is missing.")
  }
  if (length(offset) != ncol(counts)) {
    cli_abort("`offset` length ({length(offset)}) must match `ncol(counts)` ({ncol(counts)}).")
  }
  if (!is.null(names(offset)) && all(names(offset) %in% colnames(counts))) {
    offset <- offset[colnames(counts)]
  }

  cols_idx <- seq_len(ncol(counts))
  if (!is.null(group)) {
    resolved_group <- resolve_group(group, source_obj, ncol(counts), sample_role)
    cols_idx <- which(resolved_group == group[[1]] | resolved_group %in% group)
  } else if (!is.null(sample_role)) {
    cols_idx <- which(sample_role == "user")
  }

  if (length(cols_idx) == 0L) {
    cli_abort("No samples selected for cohort bootstrap.")
  }

  y <- counts[gene, cols_idx, drop = FALSE]
  o <- as.numeric(offset)[cols_idx]
  n_samples <- length(cols_idx)

  if (is.null(dispersion)) {
    samples_df <- data.frame(
      sample_id = colnames(counts),
      group = factor(rep("cohort", ncol(counts))),
      stringsAsFactors = FALSE
    )
    dge <- edgeR::DGEList(counts = counts, samples = samples_df)
    dge$offset <- matrix(as.numeric(offset), nrow = nrow(dge), ncol = ncol(dge), byrow = TRUE)
    design <- stats::model.matrix(~ 1, data = samples_df)
    dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
    disp_vec <- dispersion_from_dge(dge)
    disp <- disp_vec[[gene]]
  } else if (length(dispersion) == 1L && is.numeric(dispersion)) {
    disp <- as.numeric(dispersion)
  } else if (!is.null(names(dispersion)) && gene %in% names(dispersion)) {
    disp <- as.numeric(dispersion[[gene]])
  } else {
    disp <- as.numeric(dispersion[[1]])
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  draws <- vapply(seq_len(n_boot), function(i) {
    w <- draw_dirichlet_weights(n_samples)
    edgeR::mglmOneGroup(
      y,
      offset = o,
      dispersion = disp,
      weights = w
    )[[1]]
  }, numeric(1))

  draws
}

#' Test cohort log(mu) against HCA reference posterior draws
#'
#' Evaluates whether a cohort log(mu) differs significantly from the HCA
#' posterior baseline distribution using a heteroscedastic Welch t-test and
#' empirical posterior rank.
#'
#' @param cohort_est A data frame (e.g. from [estimate_cohort_logmu()]), a
#'   numeric vector of bootstrap draws (from [bootstrap_cohort_logmu()]), or a
#'   list with `log_mu` and `se`.
#' @param hca_draws A numeric vector of posterior draws (e.g. from
#'   [expr_draws()] or `$draws` from [expr_predict()]), or a list containing
#'   `$draws`.
#' @param gene Optional gene name/symbol.
#' @param cohort_name Label for the cohort (default `"cohort"`).
#' @return A data frame with:
#'   \item{gene}{Gene identifier}
#'   \item{cohort}{Cohort label}
#'   \item{cohort_log_mu}{Estimated cohort latent log(mu)}
#'   \item{cohort_se}{Cohort standard error}
#'   \item{hca_mean}{Mean of HCA posterior log(mu) draws}
#'   \item{hca_sd}{Standard deviation of HCA posterior log(mu) draws}
#'   \item{delta_log_mu}{Difference: `cohort_log_mu - hca_mean`}
#'   \item{se_diff}{Combined standard error: `sqrt(cohort_se^2 + hca_sd^2)`}
#'   \item{t_stat}{Welch-style t-statistic}
#'   \item{df}{Welch-Satterthwaite degrees of freedom}
#'   \item{p_value}{Two-sided p-value}
#'   \item{empirical_rank}{Proportion of HCA draws below `cohort_log_mu`}
#'   \item{direction}{Direction category: `"above_hca"`, `"below_hca"`, or `"consistent_with_hca"`}
#' @export
#' @importFrom cli cli_abort
test_cohort_vs_hca <- function(
  cohort_est,
  hca_draws,
  gene = NULL,
  cohort_name = "cohort"
) {
  if (is.list(hca_draws) && !is.data.frame(hca_draws) && "draws" %in% names(hca_draws)) {
    hca_draws <- hca_draws$draws
  }
  if (!is.numeric(hca_draws) || length(hca_draws) < 2L) {
    cli_abort("`hca_draws` must be a numeric vector with at least 2 draws.")
  }

  hca_mean <- mean(hca_draws)
  hca_sd <- stats::sd(hca_draws)
  n_hca <- length(hca_draws)

  if (is.data.frame(cohort_est)) {
    if (nrow(cohort_est) == 0L) {
      cli_abort("`cohort_est` contains 0 rows.")
    }
    if (!all(c("log_mu", "se") %in% names(cohort_est))) {
      cli_abort("`cohort_est` data frame must contain `log_mu` and `se` columns.")
    }
    if (is.null(gene) && "gene" %in% names(cohort_est)) {
      gene <- cohort_est$gene[[1]]
    }
    if ("group" %in% names(cohort_est) && (missing(cohort_name) || identical(cohort_name, "cohort"))) {
      cohort_name <- as.character(cohort_est$group[[1]])
    }
    c_log_mu <- cohort_est$log_mu[[1]]
    c_se <- cohort_est$se[[1]]
    c_n <- if ("n" %in% names(cohort_est)) cohort_est$n[[1]] else 10L
  } else if (is.numeric(cohort_est) && length(cohort_est) > 1L) {
    c_log_mu <- stats::median(cohort_est)
    c_se <- stats::sd(cohort_est)
    c_n <- length(cohort_est)
  } else if (is.list(cohort_est) && all(c("log_mu", "se") %in% names(cohort_est))) {
    c_log_mu <- as.numeric(cohort_est$log_mu[[1]])
    c_se <- as.numeric(cohort_est$se[[1]])
    c_n <- if ("n" %in% names(cohort_est)) as.numeric(cohort_est$n[[1]]) else 10L
  } else {
    cli_abort("Unrecognised format for `cohort_est`.")
  }

  delta <- c_log_mu - hca_mean
  se_diff <- sqrt(c_se^2 + hca_sd^2)
  t_stat <- delta / se_diff

  var1_mean <- c_se^2
  var2_mean <- hca_sd^2

  denom <- (var1_mean^2 / max(1L, c_n - 1L)) + (var2_mean^2 / max(1L, n_hca - 1L))
  df <- if (denom > 0) {
    (var1_mean + var2_mean)^2 / denom
  } else {
    Inf
  }

  p_val <- if (is.finite(df) && df > 0) {
    2 * stats::pt(abs(t_stat), df = df, lower.tail = FALSE)
  } else {
    2 * stats::pnorm(-abs(t_stat))
  }

  empirical_rank <- mean(hca_draws <= c_log_mu)

  direction <- if (p_val < 0.05 && delta > 0) {
    "above_hca"
  } else if (p_val < 0.05 && delta < 0) {
    "below_hca"
  } else {
    "consistent_with_hca"
  }

  data.frame(
    gene = if (!is.null(gene)) as.character(gene) else NA_character_,
    cohort = as.character(cohort_name),
    cohort_log_mu = c_log_mu,
    cohort_se = c_se,
    hca_mean = hca_mean,
    hca_sd = hca_sd,
    delta_log_mu = delta,
    se_diff = se_diff,
    t_stat = t_stat,
    df = df,
    p_value = p_val,
    empirical_rank = empirical_rank,
    direction = direction,
    stringsAsFactors = FALSE
  )
}
