#' composition_posterior_test: Perform posterior predictive analysis
#'
#' This function performs posterior predictive sampling and returns
#' a table and a density plot.
#'
#' @param proportions A data frame of observed cell type proportions.
#' @param sex A string indicating sex ("male", "female", or "unknown").
#' @param age_decade A string for age decade category.
#' @param ethnicity_groups A string indicating ethnicity group.
#' @param assay_groups A string indicating assay group.
#' @param tissue_groups A string indicating tissue group.
#' @param disease_groups Deprecated and ignored. The default healthy sccomp
#'   model is trained on healthy samples only.
#' @param load_model_to_global_env Deprecated and ignored.
#' @param fit Optional `posteriorHCA_sccomp_fit` object from [load_sccomp_fit()].
#'   When `NULL`, the default healthy model is loaded.
#' @return A list with a result table and a density plot.
#' @export
#'
#' @import dplyr tidyr purrr ggplot2 sccomp dittoSeq ggrepel
#' @importFrom scales trans_new
composition_posterior_test <-
  function(
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
      cli::cli_warn(c(
        "`disease_groups` is deprecated and ignored for composition testing.",
        "i" = "The default sccomp model is trained on healthy samples only."
      ))
    }
    if (!is.null(load_model_to_global_env)) {
      .Deprecated(msg = "`load_model_to_global_env` is deprecated and ignored.")
    }

    if (is.null(fit)) {
      fit <- load_sccomp_fit()
    }

    fit_obj <- fit$fit
    count_data <- sccomp_count_data(fit_obj)
    valid_cell_types <- unique(as.character(count_data$L3))

    if (!is.null(proportions)) {
      if (!is.data.frame(proportions) || !ncol(proportions) %in% c(2, 3)) {
        stop("Error: 'proportions' must be a data frame with either two (single sample) or three (multiple samples) columns.")
      }
      if (ncol(proportions) == 2) {
        colnames(proportions) <- c("cell_type", "proportion")
        proportions$sample_id <- "sample_1"
        proportions <- proportions[, c("sample_id", "cell_type", "proportion")]
      } else {
        colnames(proportions) <- c("sample_id", "cell_type", "proportion")
      }

      if (!all(proportions$cell_type %in% valid_cell_types)) {
        stop("Error: The 'cell_type' column must contain valid cell types.")
      }
      if (proportions %>% group_by(sample_id) %>% reframe(d = duplicated(cell_type)) %>% pull(d) %>% any()) {
        stop("Error: The 'cell_type' column must not contain duplicate values.")
      }
      if (!is.numeric(proportions$proportion) || any(proportions$proportion < 0) || any(proportions$proportion >= 1)) {
        stop("Error: The 'proportion' column must be numeric and within the range [0,1).")
      }
    }

    sample_ids <- if (!is.null(proportions)) {
      unique(as.character(proportions$sample_id))
    } else {
      "query_sample"
    }

    newdata <- do.call(
      rbind,
      lapply(sample_ids, function(sid) {
        build_sccomp_newdata(
          fit_obj,
          sample_id = sid,
          age_decade = age_decade,
          sex = sex,
          ethnicity_groups = ethnicity_groups,
          assay_groups = assay_groups,
          tissue_groups = tissue_groups
        )
      })
    )

    predict_res <- composition_draws(
      fit,
      newdata = newdata,
      summary_instead_of_draws = FALSE
    )$draws

    dist_by_cell_type <- predict_res %>%
      group_by(cell_type) %>%
      reframe(
        mean = mean(proportion),
        lower = quantile(proportion, probs = 0.025),
        upper = quantile(proportion, probs = 0.975)
      )

    arcsine_sqrt_trans <- scales::trans_new(
      "arcsine_sqrt",
      transform = function(x) asin(sqrt(x)),
      inverse = function(x) (sin(x))^2
    )

    if (!is.null(proportions)) {
      predict_res <- predict_res %>%
        filter(cell_type %in% unique(proportions$cell_type))
    }

    dist_plot <- ggplot(predict_res, aes(x = proportion)) +
      geom_density(alpha = 0.3) +
      facet_wrap(~cell_type) +
      scale_x_continuous(
        trans = arcsine_sqrt_trans,
        name = "Proportion (Arcsine-Sqrt Scaled)"
      ) +
      labs(
        x = "Proportion",
        y = "Density"
      ) +
      theme_minimal() +
      theme(strip.text = element_text(size = 12, face = "bold"))

    if (!is.null(proportions)) {
      dist_by_cell_type <- proportions %>%
        right_join(
          x = predict_res,
          by = "cell_type",
          suffix = c("_sampled", "_observed"),
          relationship = "many-to-many"
        ) %>%
        group_by(sample_id_observed, cell_type) %>%
        reframe(
          proportion_observed = unique(proportion_observed),
          Empirical_Confidence = 2 * pmin(
            mean(proportion_observed > proportion_sampled),
            1 - mean(proportion_observed > proportion_sampled)
          )
        ) %>%
        left_join(
          y = dist_by_cell_type,
          by = "cell_type",
          relationship = "many-to-many"
        )

      dist_plot <- dist_plot +
        geom_vline(
          data = dist_by_cell_type,
          aes(xintercept = proportion_observed, color = sample_id_observed),
          linetype = "dashed",
          linewidth = 0.5,
          show.legend = TRUE
        ) +
        scale_color_manual(values = dittoSeq::dittoColors()) +
        labs(color = "EC: Empirical Confidence, Sample ID:") +
        theme(
          legend.position = "bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10)
        ) +
        ggrepel::geom_text_repel(
          data = dist_by_cell_type,
          aes(
            x = proportion_observed,
            y = Inf,
            label = paste0("EC:", signif(Empirical_Confidence, 2)),
            color = sample_id_observed
          ),
          size = 4,
          direction = "y",
          segment.color = NA,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
    }

    print(dist_by_cell_type)
    print(dist_plot)

    list(
      result_table = dist_by_cell_type,
      plot = dist_plot
    )
  }

#' Query a stored gene-level brms model
#'
#' Loads the fit, builds a covariate grid with [build_newdata_grid()],
#' then draws with [expr_draws()].
#'
#' Metadata rules (same as [build_newdata_grid()]):
#' * `NA` — marginalise over every level the model knows
#' * one value — fix that level
#' * several values — marginalise over only those levels
#'
#' @param cell_type A string indicating the cell type.
#' @param gene A gene symbol or Ensembl id (e.g. `"ADRB2"` or `"ENSG00000169252"`).
#' @param age_decade,sex,disease_groups,ethnicity_groups,assay_groups,tissue_groups
#'   Metadata choices. See the rules above.
#' @param version `"latest"` or a pinned container/version (e.g. `"V1"`).
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
  cell_type,
  gene,
  age_decade = NA,
  sex = NA,
  disease_groups = NA,
  ethnicity_groups = NA,
  assay_groups = NA,
  tissue_groups = NA,
  version = "latest",
  quantity = c("linpred", "predict", "epred"),
  collapse = c("mean", "pool", "sample"),
  ndraws = NULL,
  seed = NULL
) {
  quantity <- match.arg(quantity)
  collapse <- match.arg(collapse)

  fit <- load_expr_fit(
    cell_type = cell_type,
    gene = gene,
    version = version
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