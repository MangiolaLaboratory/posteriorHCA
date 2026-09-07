#' composition_posterior_test: Perform posterior predictive analysis
#'
#' Convenience wrapper that loads the healthy sccomp model (unless `fit` is
#' supplied), draws posterior predictive proportions, tests observed
#' proportions with [composition_test()], and plots with
#' [plot_composition_vs_hca()].
#'
#' Prefer calling [load_sccomp_fit()], [composition_draws()],
#' [composition_test()], and [plot_composition_vs_hca()] directly for
#' modular workflows.
#'
#' @param proportions A data frame of observed cell type proportions.
#' @param sex,age_decade,ethnicity_groups,assay_groups,tissue_groups Metadata
#'   for the healthy query profile. Use `NULL`/`NA` to marginalise.
#' @param disease_groups Deprecated and ignored.
#' @param load_model_to_global_env Deprecated and ignored.
#' @param fit Optional `posteriorHCA_sccomp_fit` from [load_sccomp_fit()].
#' @return A list with `result_table` and `plot`.
#' @export
#' @import ggplot2
composition_posterior_test <- function(
  proportions = NULL,
  sex = NULL,
  age_decade = NULL,
  ethnicity_groups = NULL,
  assay_groups = NULL,
  tissue_groups = NULL,
  disease_groups = NULL,
  load_model_to_global_env = NULL,
  fit = NULL
) {
  if (!is.null(disease_groups)) {
    cli::cli_warn("`disease_groups` is deprecated and ignored (healthy-only model).")
  }
  if (!is.null(load_model_to_global_env)) {
    .Deprecated(msg = "`load_model_to_global_env` is deprecated and ignored.")
  }

  if (is.null(fit)) {
    fit <- load_sccomp_fit()
  }

  sample_ids <- if (is.null(proportions)) {
    "query_sample"
  } else {
    unique(normalize_proportions(proportions)$sample_id)
  }

  newdata <- do.call(
    rbind,
    lapply(sample_ids, function(sid) {
      build_sccomp_newdata(
        fit,
        sample_id = sid,
        age_decade = age_decade,
        sex = sex,
        ethnicity_groups = ethnicity_groups,
        assay_groups = assay_groups,
        tissue_groups = tissue_groups
      )
    })
  )

  draws <- composition_draws(fit, newdata = newdata)

  if (is.null(proportions)) {
    result_table <- summarize_composition_draws(draws)
    plot <- plot_composition_draws(draws)
  } else {
    result_table <- composition_test(proportions, draws)
    plot <- plot_composition_vs_hca(draws, test_results = result_table)
  }

  list(result_table = result_table, plot = plot)
}
