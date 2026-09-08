# Internal bootstrap / mglm user-estimation utilities
#
#   Not part of the supported public cohort-expression API.
#   Retained for validation / benchmarking. Preferred public user helper:
#   estimate_ql(formula, contrast) → estimate + SE.

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

#' Estimate NB dispersions (internal bootstrap path)
#'
#' Internal helper for the experimental bootstrap/mglm estimator. Fits an
#' intercept-only edgeR QL model and returns gene-wise dispersions. Not part
#' of the supported user-facing cohort-expression API; prefer
#' [estimate_ql()] for public workflows.
#'
#' @param counts Gene-by-sample numeric count matrix.
#' @param offset Numeric vector (length `ncol(counts)`) or matrix. Typically
#'   `log(effective_size)` from [calculate_tmm_scaling()] for user samples.
#' @param robust Passed to [edgeR::estimateDisp()] / [edgeR::glmQLFit()].
#' @return Named numeric vector of dispersions (names = gene ids).
#' @keywords internal
#' @noRd
estimate_dispersion_nb <- function(counts, offset, robust = TRUE) {
  counts <- as.matrix(counts)
  design <- matrix(1, nrow = ncol(counts), ncol = 1L)
  colnames(design) <- "(Intercept)"
  rownames(design) <- colnames(counts)
  fit_nb_ql(
    counts = counts,
    offset = offset,
    design = design,
    robust = robust
  )$dispersion
}

#' Bootstrap log(μ) via weighted mglmOneGroup (internal)
#'
#' Internal bootstrap estimator retained for validation/benchmarking; not part
#' of the supported user-facing cohort-expression API. Repeatedly draws
#' Dirichlet(1,…,1) weights and fits [edgeR::mglmOneGroup()] for a single gene.
#' Prefer [estimate_ql()] for public workflows.
#'
#' @param y Numeric vector or one-row matrix of counts for one gene.
#' @param offset Numeric vector of sample offsets (same length as `y`).
#' @param dispersion Positive scalar NB dispersion.
#' @param n_boot Integer number of bootstrap iterations (default `2000L`).
#' @param seed Optional RNG seed.
#' @return Numeric vector of length `n_boot`.
#' @keywords internal
#' @noRd
#' @importFrom cli cli_abort
bootstrap_logmu_mglm <- function(
  y,
  offset,
  dispersion,
  n_boot = 2000L,
  seed = NULL
) {
  n_boot <- as.integer(n_boot[[1L]])
  if (length(n_boot) != 1L || is.na(n_boot) || n_boot < 1L) {
    cli_abort("`n_boot` must be a positive integer.")
  }
  if (is.null(dispersion) || length(dispersion) != 1L ||
      !is.finite(dispersion) || dispersion <= 0) {
    cli_abort("`dispersion` must be a finite positive scalar.")
  }
  if (is.null(dim(y))) {
    y <- matrix(as.numeric(y), nrow = 1L)
  } else {
    y <- as.matrix(y)
    if (nrow(y) != 1L) {
      cli_abort("`y` must be a numeric vector or a one-row matrix.")
    }
  }
  storage.mode(y) <- "double"
  offset <- as.numeric(offset)
  if (length(offset) != ncol(y)) {
    cli_abort(
      "`offset` length ({length(offset)}) must match the number of samples ({ncol(y)})."
    )
  }

  n_samples <- ncol(y)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  vapply(seq_len(n_boot), function(i) {
    w <- draw_dirichlet_weights(n_samples)
    edgeR::mglmOneGroup(
      y,
      offset = offset,
      dispersion = as.numeric(dispersion),
      weights = w
    )[[1]]
  }, numeric(1))
}
