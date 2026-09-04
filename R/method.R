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

#' Query a stored gene-level brms model
#'
#' Builds a covariate grid with [build_newdata_grid()], then draws with
#' [expr_draws()]. Pass a pre-loaded fit from [load_expr_fit()], or supply
#' `cell_type` and `gene` to load from Nectar.
#'
#' Metadata rules (same as [build_newdata_grid()]):
#' * `NA` — marginalise over every level the model knows
#' * one value — fix that level
#' * several values — marginalise over only those levels
#'
#' @param fit Optional `posteriorHCA_expr_fit` object from [load_expr_fit()].
#' @param cell_type Cell type when `fit` is not supplied.
#' @param gene Gene symbol or Ensembl id (e.g. `"ADRB2"` or `"ENSG00000169252"`)
#'   when `fit` is not supplied.
#' @param age_decade,sex,disease_groups,ethnicity_groups,assay_groups,tissue_groups
#'   Metadata choices. See the rules above.
#' @param version `"latest"` or a pinned container/version (e.g. `"V1"`).
#' @param cache_directory,use_cache Passed to [load_expr_fit()] when `fit` is
#'   not supplied.
#' @param quantity `"linpred"` for log(μ), `"predict"` for posterior
#'   predicted counts, `"epred"` for expected counts.
#' @param collapse `"mean"` averages grid profiles within each draw.
#'   `"pool"` stacks every draw x profile. `"sample"` picks one profile
#'   per draw.
#' @param ndraws Number of posterior draws, or `NULL` for all.
#' @param seed Optional RNG seed.
#' @return A list with `pred` (data frame of draws), `summary`, `plot`
#'   (from [plot_hca_draws()]), plus `grid`, `quantity`, and `collapse`.
#'   For `quantity = "predict"` or `"epred"`, the density plot uses a
#'   log1p-scaled x-axis.
#' @export
#'
#' @import ggplot2
expr_predict <- function(
  fit = NULL,
  cell_type = NULL,
  gene = NULL,
  age_decade = NA,
  sex = NA,
  disease_groups = NA,
  ethnicity_groups = NA,
  assay_groups = NA,
  tissue_groups = NA,
  version = "latest",
  cache_directory = get_default_cache_dir(),
  use_cache = TRUE,
  quantity = c("linpred", "predict", "epred"),
  collapse = c("mean", "pool", "sample"),
  ndraws = NULL,
  seed = NULL
) {
  quantity <- match.arg(quantity)
  collapse <- match.arg(collapse)

  fit <- resolve_expr_fit(
    fit = fit,
    cell_type = cell_type,
    gene = gene,
    version = version,
    cache_directory = cache_directory,
    use_cache = use_cache
  )

  grid <- build_newdata_grid(
    fit,
    age_decade = age_decade,
    sex = sex,
    disease_groups = disease_groups,
    ethnicity_groups = ethnicity_groups,
    assay_groups = assay_groups,
    tissue_groups = tissue_groups
  )

  out <- expr_draws(
    fit,
    newdata = grid,
    quantity = quantity,
    collapse = collapse,
    ndraws = ndraws,
    seed = seed
  )

  pred_df <- data.frame(value = out$draws)

  peak_location <- tryCatch(
    {
      dens <- stats::density(out$draws)
      dens$x[which.max(dens$y)]
    },
    error = function(e) NA_real_
  )

  density_plot <- plot_hca_draws(out)

  list(
    pred = pred_df,
    summary = list(
      mean = mean(out$draws),
      median = stats::median(out$draws),
      peak_location = peak_location
    ),
    plot = density_plot,
    draws = out$draws,
    grid = out$grid,
    quantity = out$quantity,
    collapse = out$collapse,
    cell_type = out$cell_type,
    gene_ensg = out$gene_ensg,
    gene_symbol = out$gene_symbol
  )
}