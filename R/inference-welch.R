# Welch-style comparison of two mean estimates
#
#   Core:
#     welch_test_means()           — core test (mu/se/n for both sides)
#     summarize_posterior_draws()  — log_mu/se/n from HCA posterior draws
#   Wrapper (compose cores only):
#     compare_cohort_to_hca()      — summarise draws + welch_test_means per row

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
#' @seealso [compare_cohort_to_hca()], [summarize_posterior_draws()]
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

#' Summarise posterior draws for comparison with a point estimate
#'
#' @param draws A numeric vector or a list with a `$draws` element (for
#'   example from [expression_draws()]).
#' @param value Optional scalar; when supplied, `empirical_rank` is the
#'   proportion of draws less than or equal to `value`.
#' @return A list with `log_mu`, `se`, `n`, and optionally `empirical_rank`.
#' @seealso [compare_cohort_to_hca()], [welch_test_means()]
#' @export
#' @importFrom cli cli_abort
summarize_posterior_draws <- function(draws, value = NULL) {
  x <- if (is.numeric(draws)) {
    as.numeric(draws)
  } else if (is.list(draws) && "draws" %in% names(draws)) {
    as.numeric(draws$draws)
  } else {
    cli_abort("`draws` must be a numeric vector or a list with a `$draws` element.")
  }
  if (length(x) < 2L) {
    cli_abort("`draws` must contain at least 2 numeric values.")
  }
  out <- list(
    log_mu = unname(mean(x)),
    se = unname(stats::sd(x)),
    n = length(x)
  )
  if (!is.null(value)) {
    out$empirical_rank <- mean(x <= value)
  }
  out
}

#' Compare cohort log(μ) estimates to an HCA posterior baseline
#'
#' Convenience wrapper around [summarize_posterior_draws()] and
#' [welch_test_means()]. Summarises `hca_draws` once, then runs the Welch
#' comparison for each row of `cohort_estimates`. Does not implement Welch
#' mathematics itself.
#'
#' @param cohort_estimates Data frame from [estimate_logmu_ql()] /
#'   [estimate_cohort_logmu()] with columns `log_mu`, `se`, and preferably
#'   `group`, `gene`, and `n`.
#' @param hca_draws Numeric draws or a list from [expression_draws()] /
#'   [expression_baseline_draws()].
#' @param alternative Passed to [welch_test_means()].
#' @return A data frame with one row per cohort estimate, including
#'   `delta`, `se_diff`, `t_stat`, `df`, and `p_value`.
#' @seealso [summarize_posterior_draws()], [welch_test_means()],
#'   [estimate_cohort_logmu()], [expression_baseline_draws()]
#' @export
#' @importFrom cli cli_abort
compare_cohort_to_hca <- function(
  cohort_estimates,
  hca_draws,
  alternative = c("two.sided", "greater", "less")
) {
  alternative <- match.arg(alternative)
  cohort_estimates <- as.data.frame(cohort_estimates)
  required <- c("log_mu", "se")
  missing <- setdiff(required, names(cohort_estimates))
  if (length(missing)) {
    cli_abort("`cohort_estimates` must contain column{?s}: {missing}.")
  }
  if (!nrow(cohort_estimates)) {
    cli_abort("`cohort_estimates` has zero rows.")
  }

  hca_summary <- summarize_posterior_draws(hca_draws)

  rows <- lapply(seq_len(nrow(cohort_estimates)), function(i) {
    cohort_row <- cohort_estimates[i, , drop = FALSE]
    n1 <- if ("n" %in% names(cohort_row)) cohort_row$n[[1]] else NULL
    test <- welch_test_means(
      mu1 = cohort_row$log_mu[[1]],
      se1 = cohort_row$se[[1]],
      mu2 = hca_summary$log_mu,
      se2 = hca_summary$se,
      n1 = n1,
      n2 = hca_summary$n,
      alternative = alternative
    )
    data.frame(
      gene = if ("gene" %in% names(cohort_row)) {
        as.character(cohort_row$gene[[1]])
      } else {
        NA_character_
      },
      group = if ("group" %in% names(cohort_row)) {
        as.character(cohort_row$group[[1]])
      } else {
        as.character(i)
      },
      log_mu = test$mu1,
      se = test$se1,
      n = if (is.null(test$n1)) NA_integer_ else as.integer(test$n1),
      hca_log_mu = test$mu2,
      hca_se = test$se2,
      hca_n = if (is.null(test$n2)) NA_integer_ else as.integer(test$n2),
      delta = test$delta,
      se_diff = test$se_diff,
      t_stat = test$t_stat,
      df = test$df,
      p_value = test$p_value,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}
