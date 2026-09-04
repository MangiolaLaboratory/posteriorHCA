# Composition testing against healthy HCA sccomp draws
#
# Layering (mirrors gene-expression inference):
#   empirical_proportion_confidence() — core EC for one observed vs draws
#   summarize_composition_draws()     — mean / CI by cell type
#   normalize_proportions()           — tidy observed proportion table
#   composition_test()                — batch EC table (one or more samples)
#   composition_posterior_test()      — convenience wrapper (method.R)

#' Empirical confidence of an observed proportion against posterior draws
#'
#' Two-sided tail probability, scaled to \eqn{[0, 1]}:
#' `2 * min(P(obs > draw), P(obs <= draw))`.
#'
#' @param observed Scalar observed proportion.
#' @param draws Numeric vector of posterior predictive proportions.
#' @return Scalar empirical confidence.
#' @export
empirical_proportion_confidence <- function(observed, draws) {
  observed <- as.numeric(observed[[1]])
  draws <- as.numeric(draws)
  2 * min(mean(observed > draws), 1 - mean(observed > draws))
}

#' Summarise composition draws by cell type
#'
#' @param draws Data frame from [composition_draws()]`$draws`, or that list
#'   itself. Must contain `cell_type` and `proportion`.
#' @return Data frame with `cell_type`, `mean`, `lower`, `upper`.
#' @export
summarize_composition_draws <- function(draws) {
  draws <- composition_draws_df(draws)
  cell_types <- unique(draws$cell_type)
  do.call(rbind, lapply(cell_types, function(ct) {
    x <- draws$proportion[draws$cell_type == ct]
    data.frame(
      cell_type = ct,
      mean = mean(x),
      lower = unname(stats::quantile(x, 0.025)),
      upper = unname(stats::quantile(x, 0.975)),
      stringsAsFactors = FALSE
    )
  }))
}

#' Normalise an observed proportions table
#'
#' Accepts either two columns (`cell_type`, `proportion`) or three
#' (`sample_id`, `cell_type`, `proportion`).
#'
#' @param proportions Data frame of observed proportions.
#' @return Data frame with columns `sample_id`, `cell_type`, `proportion`.
#' @export
normalize_proportions <- function(proportions) {
  proportions <- as.data.frame(proportions)
  if (ncol(proportions) == 2L) {
    names(proportions) <- c("cell_type", "proportion")
    proportions$sample_id <- "sample_1"
    proportions <- proportions[, c("sample_id", "cell_type", "proportion")]
  } else {
    names(proportions)[seq_len(3L)] <- c("sample_id", "cell_type", "proportion")
    proportions <- proportions[, c("sample_id", "cell_type", "proportion")]
  }
  proportions$sample_id <- as.character(proportions$sample_id)
  proportions$cell_type <- as.character(proportions$cell_type)
  proportions$proportion <- as.numeric(proportions$proportion)
  proportions
}

#' Extract draws data frame from composition_draws() output
#' @keywords internal
#' @noRd
composition_draws_df <- function(draws) {
  if (is.list(draws) && !is.data.frame(draws) && "draws" %in% names(draws)) {
    draws <- draws$draws
  }
  draws <- as.data.frame(draws)
  if ("L3" %in% names(draws) && !"cell_type" %in% names(draws)) {
    draws$cell_type <- as.character(draws$L3)
  }
  draws$cell_type <- as.character(draws$cell_type)
  draws
}

#' Test observed proportions against healthy composition draws
#'
#' Computes empirical confidence (EC) for each sample × cell type against
#' posterior predictive draws from [composition_draws()].
#'
#' @param proportions Observed proportions (see [normalize_proportions()]).
#' @param draws Output of [composition_draws()], or its `$draws` data frame.
#' @param samples Optional sample ids to keep. Default is all samples in
#'   `proportions`.
#' @return Data frame with `sample_id`, `cell_type`, `proportion_observed`,
#'   `empirical_confidence`, `mean`, `lower`, `upper`.
#' @export
composition_test <- function(proportions, draws, samples = NULL) {
  proportions <- normalize_proportions(proportions)
  if (!is.null(samples)) {
    proportions <- proportions[proportions$sample_id %in% samples, , drop = FALSE]
  }

  draws <- composition_draws_df(draws)
  baseline <- summarize_composition_draws(draws)

  cell_types <- unique(proportions$cell_type)
  draws <- draws[draws$cell_type %in% cell_types, , drop = FALSE]

  rows <- lapply(seq_len(nrow(proportions)), function(i) {
    obs <- proportions[i, , drop = FALSE]
    cell_draws <- draws$proportion[draws$cell_type == obs$cell_type]
    base <- baseline[baseline$cell_type == obs$cell_type, , drop = FALSE]
    data.frame(
      sample_id = obs$sample_id,
      cell_type = obs$cell_type,
      proportion_observed = obs$proportion,
      empirical_confidence = empirical_proportion_confidence(
        obs$proportion,
        cell_draws
      ),
      mean = base$mean[[1]],
      lower = base$lower[[1]],
      upper = base$upper[[1]],
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

#' Test one sample's proportions against healthy composition draws
#'
#' Thin wrapper around [composition_test()] for a single `sample_id`.
#'
#' @inheritParams composition_test
#' @param sample_id Sample to keep. If `NULL` and `proportions` has one sample,
#'   that sample is used.
#' @return Same columns as [composition_test()].
#' @export
composition_test_sample <- function(proportions, draws, sample_id = NULL) {
  proportions <- normalize_proportions(proportions)
  if (is.null(sample_id)) {
    sample_id <- unique(proportions$sample_id)
    if (length(sample_id) != 1L) {
      cli::cli_abort("Provide `sample_id` when `proportions` contains multiple samples.")
    }
  }
  composition_test(proportions, draws, samples = sample_id)
}

#' Test multiple samples against healthy composition draws
#'
#' Thin wrapper around [composition_test()] kept for symmetry with
#' [welch_t_test_cohort_hca()].
#'
#' @inheritParams composition_test
#' @return Same columns as [composition_test()].
#' @export
composition_test_cohort <- function(proportions, draws) {
  composition_test(proportions, draws)
}
