#' Personalized single-sample predictive test against healthy counterfactual
#'
#' Evaluates whether one or more observed counts are surprising under a
#' posterior predictive distribution (`posterior_predict`) from a matched
#' healthy atlas model.
#'
#' @param observed_count Numeric vector of observed counts (one per sample).
#' @param offset Numeric vector of library-size offsets (default `0`).
#' @param predictive_draws A numeric vector of draws (for 1 sample), a draws-by-sample
#'   matrix, or a list containing `$draws` from [expr_draws()] with `quantity = "predict"`.
#' @param linpred_draws Optional numeric vector or matrix of latent log(mu) draws from
#'   `posterior_linpred(transform = FALSE)`.
#' @param epred_draws Optional numeric vector or matrix of expected count draws from
#'   `posterior_epred()`.
#' @param nb_shape Optional scalar or vector of negative-binomial shape parameters
#'   (for latent delta-method SE calculation).
#' @param sample_id Optional character vector of sample identifiers.
#' @param gene Optional gene name or symbol.
#' @return A data frame with one row per sample:
#'   \item{gene}{Gene identifier}
#'   \item{sample_id}{Sample identifier}
#'   \item{observed_count}{Observed count}
#'   \item{offset}{Library size offset}
#'   \item{n_draws}{Number of posterior predictive draws}
#'   \item{healthy_predictive_median}{Median of predicted counts}
#'   \item{healthy_predictive_lo}{2.5\% quantile of predicted counts}
#'   \item{healthy_predictive_hi}{97.5\% quantile of predicted counts}
#'   \item{predictive_p_lower}{Lower predictive tail probability}
#'   \item{predictive_p_upper}{Upper predictive tail probability}
#'   \item{predictive_percentile}{Mid-rank predictive percentile}
#'   \item{predictive_tail_p}{Two-sided predictive tail probability}
#'   \item{predictive_z}{Standard normal deviate score}
#'   \item{absolute_predictive_z}{Absolute Z-score}
#'   \item{direction}{Direction vs healthy baseline}
#'   \item{healthy_epred_median}{Median expected count (if `epred_draws` supplied)}
#'   \item{log2fc_vs_healthy_epred}{Log2 fold change vs expected count (if `epred_draws` supplied)}
#'   \item{patient_log_mu}{Estimated patient log(mu) (if `linpred_draws` supplied)}
#'   \item{healthy_log_mu_mean}{Mean of healthy linpred draws (if `linpred_draws` supplied)}
#'   \item{healthy_log_mu_sd}{SD of healthy linpred draws (if `linpred_draws` supplied)}
#'   \item{welch_style_z}{Heteroscedastic log-mean Z-score (if `linpred_draws` supplied)}
#'   \item{welch_style_p}{Two-sided p-value from log-mean Z-score (if `linpred_draws` supplied)}
#' @export
#' @importFrom cli cli_abort
test_sample_predictive <- function(
  observed_count,
  offset = 0,
  predictive_draws,
  linpred_draws = NULL,
  epred_draws = NULL,
  nb_shape = NULL,
  sample_id = NULL,
  gene = NULL
) {
  if (missing(observed_count) || is.null(observed_count)) {
    cli_abort("`observed_count` is missing.")
  }
  if (is.null(sample_id)) {
    sample_id <- if (!is.null(names(observed_count))) {
      names(observed_count)
    } else if (length(observed_count) == 1L) {
      "sample_1"
    } else {
      paste0("sample_", seq_along(observed_count))
    }
  }

  observed_count <- as.numeric(observed_count)
  n_samples <- length(observed_count)

  # Normalise predictive draws to matrix
  if (is.list(predictive_draws) && !is.data.frame(predictive_draws) && "draws" %in% names(predictive_draws)) {
    predictive_draws <- predictive_draws$draws
  }
  if (is.null(dim(predictive_draws))) {
    predictive_mat <- matrix(as.numeric(predictive_draws), ncol = 1L)
  } else {
    predictive_mat <- as.matrix(predictive_draws)
  }

  if (ncol(predictive_mat) == 1L && n_samples > 1L) {
    predictive_mat <- matrix(
      rep(predictive_mat[, 1L], n_samples),
      ncol = n_samples
    )
  }
  if (ncol(predictive_mat) != n_samples) {
    cli_abort(
      "Number of columns in `predictive_draws` ({ncol(predictive_mat)}) must match length of `observed_count` ({n_samples})."
    )
  }

  # Normalise epred draws if given
  epred_mat <- NULL
  if (!is.null(epred_draws)) {
    if (is.list(epred_draws) && !is.data.frame(epred_draws) && "draws" %in% names(epred_draws)) {
      epred_draws <- epred_draws$draws
    }
    if (is.null(dim(epred_draws))) {
      epred_mat <- matrix(as.numeric(epred_draws), ncol = 1L)
    } else {
      epred_mat <- as.matrix(epred_draws)
    }
    if (ncol(epred_mat) == 1L && n_samples > 1L) {
      epred_mat <- matrix(rep(epred_mat[, 1L], n_samples), ncol = n_samples)
    }
  }

  # Normalise linpred draws if given
  linpred_mat <- NULL
  if (!is.null(linpred_draws)) {
    if (is.list(linpred_draws) && !is.data.frame(linpred_draws) && "draws" %in% names(linpred_draws)) {
      linpred_draws <- linpred_draws$draws
    }
    if (is.null(dim(linpred_draws))) {
      linpred_mat <- matrix(as.numeric(linpred_draws), ncol = 1L)
    } else {
      linpred_mat <- as.matrix(linpred_draws)
    }
    if (ncol(linpred_mat) == 1L && n_samples > 1L) {
      linpred_mat <- matrix(rep(linpred_mat[, 1L], n_samples), ncol = n_samples)
    }
  }

  if (!is.null(nb_shape) && length(nb_shape) == 1L && n_samples > 1L) {
    nb_shape <- rep(nb_shape, n_samples)
  }

  n_draws <- nrow(predictive_mat)

  rows <- lapply(seq_len(n_samples), function(i) {
    y <- observed_count[[i]]
    off <- offset[[i]]
    sid <- sample_id[[i]]
    d <- predictive_mat[, i]

    p_lo <- (1 + sum(d <= y)) / (n_draws + 1)
    p_hi <- (1 + sum(d >= y)) / (n_draws + 1)
    p_percentile <- (0.5 + sum(d < y) + 0.5 * sum(d == y)) / (n_draws + 1)
    p_tail <- min(1, 2 * min(p_lo, p_hi))

    p_clipped <- min(1 - 0.5 / (n_draws + 1), max(0.5 / (n_draws + 1), p_percentile))
    z_score <- stats::qnorm(p_clipped)

    dir_cat <- if (p_lo <= 0.025) {
      "below_personalised_healthy"
    } else if (p_hi <= 0.025) {
      "above_personalised_healthy"
    } else {
      "consistent_with_personalised_healthy"
    }

    res <- data.frame(
      gene = if (!is.null(gene)) as.character(gene) else NA_character_,
      sample_id = as.character(sid),
      observed_count = y,
      offset = off,
      n_draws = n_draws,
      healthy_predictive_median = stats::median(d),
      healthy_predictive_lo = stats::quantile(d, 0.025),
      healthy_predictive_hi = stats::quantile(d, 0.975),
      predictive_p_lower = p_lo,
      predictive_p_upper = p_hi,
      predictive_percentile = p_percentile,
      predictive_tail_p = p_tail,
      predictive_z = z_score,
      absolute_predictive_z = abs(z_score),
      direction = dir_cat,
      stringsAsFactors = FALSE
    )

    if (!is.null(epred_mat)) {
      ep <- epred_mat[, i]
      ep_med <- stats::median(ep)
      res$healthy_epred_median <- ep_med
      res$healthy_epred_lo <- stats::quantile(ep, 0.025)
      res$healthy_epred_hi <- stats::quantile(ep, 0.975)
      res$log2fc_vs_healthy_epred <- log2((y + 0.5) / (ep_med + 0.5))
    }

    if (!is.null(linpred_mat)) {
      lp <- linpred_mat[, i]
      lp_mean <- mean(lp)
      lp_sd <- stats::sd(lp)
      res$healthy_log_mu_mean <- lp_mean
      res$healthy_log_mu_sd <- lp_sd
      res$healthy_log_mu_lo <- stats::quantile(lp, 0.025)
      res$healthy_log_mu_hi <- stats::quantile(lp, 0.975)

      if (y > 0) {
        pat_log_mu <- log(y) - off
        shape_val <- if (!is.null(nb_shape)) nb_shape[[i]] else Inf
        pat_se <- sqrt(1 / y + (if (is.finite(shape_val) && shape_val > 0) 1 / shape_val else 0))
        delta_mu <- pat_log_mu - lp_mean
        w_test <- welch_test_means(
          mu1 = pat_log_mu,
          se1 = pat_se,
          mu2 = lp_mean,
          se2 = lp_sd,
          n2 = length(lp)
        )

        res$patient_log_mu <- pat_log_mu
        res$patient_log_mu_se <- pat_se
        res$delta_log_mu_vs_healthy <- delta_mu
        res$delta_log_mu_se <- w_test$se_diff
        res$welch_style_z <- w_test$t_stat
        res$welch_style_p <- w_test$p_value
      } else {
        res$patient_log_mu <- NA_real_
        res$patient_log_mu_se <- NA_real_
        res$delta_log_mu_vs_healthy <- NA_real_
        res$delta_log_mu_se <- NA_real_
        res$welch_style_z <- NA_real_
        res$welch_style_p <- NA_real_
      }
    }

    res
  })

  do.call(rbind, rows)
}
