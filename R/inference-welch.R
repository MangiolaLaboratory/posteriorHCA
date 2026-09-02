# Generic Welch-style comparison of two mean estimates
#
# Layering:
#   welch_test_means()           — core test (mu/se/n for both sides)
#   summarize_posterior_draws()  — mean/sd/n (+ rank) from HCA draws
#   cohort_estimate_at()         — mu/se/n from estimate_cohort_logmu output
#   compare_logmu_to_baseline()  — one cohort vs baseline (draws or summary)
#   welch_t_test_cohort_hca()    — batch wrapper (in inference-cohort.R)

#' Welch-Satterthwaite degrees of freedom for two mean estimates
#' @keywords internal
#' @noRd
welch_satterthwaite_df <- function(se1, se2, n1 = NULL, n2 = NULL) {
  v1 <- se1^2
  v2 <- se2^2
  if (is.null(n1) || is.null(n2) || n1 < 2L || n2 < 2L) {
    return(Inf)
  }
  denom <- (v1^2 / (n1 - 1L)) + (v2^2 / (n2 - 1L))
  if (!is.finite(denom) || denom <= 0) {
    return(Inf)
  }
  (v1 + v2)^2 / denom
}

#' Heteroscedastic Welch test for the difference between two mean estimates
#'
#' Compares two estimated means on the same scale (typically latent log(μ))
#' using a Welch-style statistic. When `n1` and `n2` are supplied and both
#' are at least 2, p-values use the Welch–Satterthwaite t distribution;
#' otherwise a normal approximation is used.
#'
#' @param mu1,mu2 Point estimates of the two means.
#' @param se1,se2 Standard errors (or posterior SDs) associated with `mu1`
#'   and `mu2`.
#' @param n1,n2 Optional sample or draw counts used for Welch–Satterthwaite
#'   degrees of freedom.
#' @param alternative `"two.sided"` (default), `"greater"` (mu1 > mu2), or
#'   `"less"` (mu1 < mu2).
#' @return A list with `mu1`, `se1`, `n1`, `mu2`, `se2`, `n2`, `delta`
#'   (`mu1 - mu2`), `se_diff`, `t_stat`, `df`, and `p_value`.
#' @export
#' @importFrom cli cli_abort
welch_test_means <- function(
  mu1,
  se1,
  mu2,
  se2,
  n1 = NULL,
  n2 = NULL,
  alternative = c("two.sided", "greater", "less")
) {
  alternative <- match.arg(alternative)
  if (length(mu1) != 1L || length(mu2) != 1L || length(se1) != 1L || length(se2) != 1L) {
    cli_abort("`mu1`, `se1`, `mu2`, and `se2` must be length-1 numeric values.")
  }
  if (!is.numeric(mu1) || !is.numeric(mu2) || !is.finite(mu1) || !is.finite(mu2)) {
    cli_abort("`mu1` and `mu2` must be finite numeric values.")
  }
  if (!is.numeric(se1) || !is.numeric(se2) || !is.finite(se1) || !is.finite(se2) || se1 < 0 || se2 < 0) {
    cli_abort("`se1` and `se2` must be finite non-negative numeric values.")
  }

  delta <- mu1 - mu2
  se_diff <- sqrt(se1^2 + se2^2)
  t_stat <- if (se_diff > 0) delta / se_diff else if (delta == 0) 0 else Inf * sign(delta)
  df <- welch_satterthwaite_df(se1, se2, n1, n2)

  p_value <- if (is.finite(df) && df > 0 && is.finite(t_stat)) {
    switch(
      alternative,
      two.sided = 2 * stats::pt(abs(t_stat), df = df, lower.tail = FALSE),
      greater = stats::pt(t_stat, df = df, lower.tail = FALSE),
      less = stats::pt(t_stat, df = df, lower.tail = TRUE)
    )
  } else if (is.finite(t_stat)) {
    switch(
      alternative,
      two.sided = 2 * stats::pnorm(-abs(t_stat)),
      greater = stats::pnorm(t_stat, lower.tail = FALSE),
      less = stats::pnorm(t_stat, lower.tail = TRUE)
    )
  } else {
    NA_real_
  }

  list(
    mu1 = unname(as.numeric(mu1)),
    se1 = unname(as.numeric(se1)),
    n1 = n1,
    mu2 = unname(as.numeric(mu2)),
    se2 = unname(as.numeric(se2)),
    n2 = n2,
    delta = unname(as.numeric(delta)),
    se_diff = unname(as.numeric(se_diff)),
    t_stat = unname(as.numeric(t_stat)),
    df = unname(as.numeric(df)),
    p_value = unname(as.numeric(p_value))
  )
}

#' Extract posterior draws from common expression-model outputs
#' @keywords internal
#' @noRd
extract_posterior_draws <- function(draws) {
  if (is.numeric(draws)) {
    return(as.numeric(draws))
  }
  if (is.list(draws) && !is.data.frame(draws)) {
    if ("draws" %in% names(draws)) {
      return(as.numeric(draws$draws))
    }
    if ("pred" %in% names(draws) && is.data.frame(draws$pred) && "value" %in% names(draws$pred)) {
      return(as.numeric(draws$pred$value))
    }
  }
  cli_abort(c(
    "`draws` must be a numeric vector or a list from [expr_draws()] / [expr_predict()].",
    "i" = "Expected `$draws` or `$pred$value`."
  ))
}

#' Summarise posterior draws for comparison with a point estimate
#'
#' @param draws A numeric vector or a list from [expr_draws()] /
#'   [expr_predict()].
#' @param value Optional scalar; when supplied, `empirical_rank` is the
#'   proportion of draws less than or equal to `value`.
#' @return A list with `mean`, `sd`, `n`, and optionally `empirical_rank`.
#' @export
#' @importFrom cli cli_abort
summarize_posterior_draws <- function(draws, value = NULL) {
  x <- extract_posterior_draws(draws)
  if (length(x) < 2L) {
    cli_abort("`draws` must contain at least 2 numeric values.")
  }
  out <- list(
    mean = unname(mean(x)),
    sd = unname(stats::sd(x)),
    n = length(x)
  )
  if (!is.null(value)) {
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value)) {
      cli_abort("`value` must be a single finite numeric value.")
    }
    out$empirical_rank <- mean(x <= value)
  }
  out
}

#' Empirical posterior rank of a value against draws
#'
#' @param draws Numeric vector or expression-model draw object.
#' @param value Scalar value to rank (e.g. cohort log(μ)).
#' @return Proportion of draws less than or equal to `value`.
#' @export
empirical_posterior_rank <- function(draws, value) {
  summarize_posterior_draws(draws, value = value)$empirical_rank
}

#' Extract one cohort estimate row for manual testing
#'
#' @param cohort_est Output from [estimate_cohort_logmu()] or
#'   [bootstrap_cohort_logmu_batch()].
#' @param group Cohort label to extract. Required when `cohort_est` has
#'   more than one row.
#' @param row Optional row index when `group` is not supplied.
#' @return A list with `group`, `mu` (`log_mu`), `se`, `n`, `method`, and
#'   metadata fields from [cohort_est_metadata()].
#' @export
#' @importFrom cli cli_abort
cohort_estimate_at <- function(cohort_est, group = NULL, row = NULL) {
  if (!is.data.frame(cohort_est) || nrow(cohort_est) == 0L) {
    cli_abort("`cohort_est` must be a non-empty data frame.")
  }
  if (!all(c("group", "log_mu", "se") %in% names(cohort_est))) {
    cli_abort("`cohort_est` must contain `group`, `log_mu`, and `se` columns.")
  }

  if (!is.null(group)) {
    sub <- cohort_est[cohort_est$group == group, , drop = FALSE]
    if (nrow(sub) == 0L) {
      cli_abort("No rows in `cohort_est` match `group = \"{group}\"`.")
    }
    if (nrow(sub) > 1L) {
      cli_abort(
        "`cohort_est` has {nrow(sub)} rows for `group = \"{group}\"`; subset to one gene first."
      )
    }
    cohort_row <- sub
  } else if (!is.null(row)) {
    row <- as.integer(row[[1L]])
    if (length(row) != 1L || is.na(row) || row < 1L || row > nrow(cohort_est)) {
      cli_abort("`row` must be a valid row index into `cohort_est`.")
    }
    cohort_row <- cohort_est[row, , drop = FALSE]
  } else if (nrow(cohort_est) == 1L) {
    cohort_row <- cohort_est
  } else {
    cli_abort("Supply `group` or `row` when `cohort_est` has more than one row.")
  }

  meta <- cohort_est_metadata(cohort_row)
  list(
    group = as.character(cohort_row$group[[1]]),
    mu = unname(as.numeric(cohort_row$log_mu[[1]])),
    se = unname(as.numeric(cohort_row$se[[1]])),
    n = if ("n" %in% names(cohort_row)) cohort_row$n[[1]] else NULL,
    method = if ("method" %in% names(cohort_row)) {
      as.character(cohort_row$method[[1]])
    } else {
      "ql"
    },
    gene = meta$gene_ensg,
    gene_symbol = meta$gene_symbol,
    cell_type = meta$cell_type
  )
}

#' Classify a mean difference relative to a reference
#' @keywords internal
#' @noRd
classify_vs_reference <- function(delta, p_value, alpha = 0.05) {
  if (!is.finite(p_value) || p_value >= alpha) {
    return("consistent_with_hca")
  }
  if (delta > 0) {
    "above_hca"
  } else if (delta < 0) {
    "below_hca"
  } else {
    "consistent_with_hca"
  }
}

#' Compare one cohort log(μ) estimate to a baseline
#'
#' Low-level comparison helper combining [welch_test_means()],
#' [summarize_posterior_draws()], and empirical rank. Supply either
#' `baseline_draws` or explicit `baseline_mu`, `baseline_se`, and optionally
#' `baseline_n`.
#'
#' @param mu,se Cohort point estimate and standard error (typically `log_mu`
#'   and `se` from [estimate_cohort_logmu()]).
#' @param n Optional cohort sample size for Welch degrees of freedom.
#' @param baseline_draws Numeric vector or [expr_draws()] output for the
#'   reference baseline.
#' @param baseline_mu,baseline_se,baseline_n Baseline summary when draws are
#'   not supplied.
#' @param cohort Optional cohort label stored in the output.
#' @param method Optional method label (`"ql"`, `"bootstrap"`, etc.).
#' @param gene,gene_symbol,cell_type Optional annotation columns.
#' @param alpha Significance level for the `direction` column.
#' @inheritParams welch_test_means
#' @return A one-row data frame compatible with [welch_t_test_cohort_hca()]
#'   and [plot_cohort_vs_hca()].
#' @export
#' @importFrom cli cli_abort
compare_logmu_to_baseline <- function(
  mu,
  se,
  n = NULL,
  baseline_draws = NULL,
  baseline_mu = NULL,
  baseline_se = NULL,
  baseline_n = NULL,
  cohort = "cohort",
  method = "ql",
  gene = NA_character_,
  gene_symbol = NA_character_,
  cell_type = NA_character_,
  alpha = 0.05,
  alternative = c("two.sided", "greater", "less")
) {
  if (!is.null(baseline_draws)) {
    baseline <- summarize_posterior_draws(baseline_draws, value = mu)
    baseline_mu <- baseline$mean
    baseline_se <- baseline$sd
    baseline_n <- baseline$n
    empirical_rank <- baseline$empirical_rank
  } else {
    if (is.null(baseline_mu) || is.null(baseline_se)) {
      cli_abort(c(
        "Supply `baseline_draws`, or both `baseline_mu` and `baseline_se`.",
        "i" = "Use [summarize_posterior_draws()] on HCA draws for a manual workflow."
      ))
    }
    empirical_rank <- NA_real_
  }

  test <- welch_test_means(
    mu1 = mu,
    se1 = se,
    mu2 = baseline_mu,
    se2 = baseline_se,
    n1 = n,
    n2 = baseline_n,
    alternative = alternative
  )

  data.frame(
    gene = as.character(gene),
    gene_symbol = as.character(gene_symbol),
    cell_type = as.character(cell_type),
    cohort = as.character(cohort),
    method = as.character(method),
    cohort_log_mu = test$mu1,
    cohort_se = test$se1,
    hca_mean = test$mu2,
    hca_sd = test$se2,
    delta_log_mu = test$delta,
    se_diff = test$se_diff,
    t_stat = test$t_stat,
    df = test$df,
    p_value = test$p_value,
    empirical_rank = empirical_rank,
    direction = classify_vs_reference(test$delta, test$p_value, alpha = alpha),
    stringsAsFactors = FALSE
  )
}