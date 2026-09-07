# Welch-style comparison of two mean estimates
#
#   welch_test_means()           — core test (mu/se/n for both sides)
#   summarize_posterior_draws()  — mean/sd/n from HCA posterior draws

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

#' Summarise posterior draws for comparison with a point estimate
#'
#' @param draws A numeric vector or a list with a `$draws` element (for
#'   example from [expression_draws()]).
#' @param value Optional scalar; when supplied, `empirical_rank` is the
#'   proportion of draws less than or equal to `value`.
#' @return A list with `mean`, `sd`, `n`, and optionally `empirical_rank`.
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
    mean = unname(mean(x)),
    sd = unname(stats::sd(x)),
    n = length(x)
  )
  if (!is.null(value)) {
    out$empirical_rank <- mean(x <= value)
  }
  out
}
