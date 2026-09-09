# Plotting utilities for posteriorHCA expression outputs

#' Extract posterior draws and metadata from common result objects
#' @param x Output from [expression_draws()], or a numeric vector.
#' @param quantity Override quantity when `x` is a bare numeric vector.
#' @keywords internal
#' @noRd
normalize_draws_input <- function(x, quantity = NULL) {
  if (is.numeric(x)) {
    if (is.null(quantity) || length(quantity) != 1L) {
      cli::cli_abort(
        "`quantity` must be supplied when `draws` is a numeric vector."
      )
    }
    return(list(
      draws = as.numeric(x),
      quantity = match.arg(quantity, c("linpred", "predict", "epred")),
      gene_ensg = NA_character_,
      cell_type = NA_character_
    ))
  }

  if (!is.list(x) || is.data.frame(x)) {
    cli::cli_abort(
      "`draws` must be a numeric vector or a list from [expression_draws()]."
    )
  }

  if (!"draws" %in% names(x)) {
    cli::cli_abort(
      "Could not find posterior draws in `draws` (expected `$draws`)."
    )
  }
  draws_vec <- x$draws

  if (!is.numeric(draws_vec) || length(draws_vec) < 2L) {
    cli::cli_abort("`draws` must contain at least 2 numeric posterior draws.")
  }

  qty <- if (!is.null(quantity)) {
    match.arg(quantity, c("linpred", "predict", "epred"))
  } else if (!is.null(x$quantity)) {
    as.character(x$quantity[[1]])
  } else {
    "linpred"
  }

  meta_gene <- if (!is.null(x$gene_ensg)) as.character(x$gene_ensg[[1]]) else NA_character_
  meta_ct <- if (!is.null(x$cell_type)) as.character(x$cell_type[[1]]) else NA_character_
  list(
    draws = as.numeric(draws_vec),
    quantity = qty,
    gene_ensg = meta_gene,
    cell_type = meta_ct
  )
}

#' X-axis label for an expression-model quantity
#' @keywords internal
#' @noRd
expression_quantity_xlab <- function(quantity) {
  switch(
    quantity,
    linpred = "log(mu)",
    predict = "log1p(predicted count)",
    epred = "log1p(expected count)",
    "value"
  )
}

#' Apply quantity-appropriate x-axis scaling to a ggplot object
#' @keywords internal
#' @noRd
apply_quantity_x_scale <- function(plot, quantity, xlab) {
  if (quantity %in% c("predict", "epred")) {
    plot + ggplot2::scale_x_continuous(
      trans = scales::log1p_trans(),
      name = xlab
    )
  } else {
    plot + ggplot2::labs(x = xlab)
  }
}

#' Build a default title for HCA draw density plots
#' @keywords internal
#' @noRd
default_hca_draws_title <- function(meta, quantity, comparison = FALSE) {
  gene_lab <- if (!is.na(meta$gene_ensg) && nzchar(meta$gene_ensg)) {
    meta$gene_ensg
  } else {
    "Gene"
  }

  ct_lab <- if (!is.na(meta$cell_type) && nzchar(meta$cell_type)) {
    paste0(" (", meta$cell_type, ")")
  } else {
    ""
  }

  baseline_lab <- if (comparison) {
    "Cohort vs healthy HCA log(mu)"
  } else {
    switch(
      quantity,
      linpred = "Healthy HCA log(mu) posterior",
      predict = "Healthy HCA predicted count posterior",
      epred = "Healthy HCA expected count posterior",
      "Healthy HCA posterior"
    )
  }

  paste0(baseline_lab, ": ", gene_lab, ct_lab)
}

#' Build the baseline HCA density ggplot
#' @keywords internal
#' @noRd
build_hca_density_plot <- function(
  draws_vec,
  quantity,
  fill = "#4C78A8"
) {
  xlab <- expression_quantity_xlab(quantity)
  plot_df <- data.frame(value = draws_vec)

  p <- ggplot(plot_df, aes(x = .data$value)) +
    geom_density(fill = fill, colour = NA, alpha = 0.45, linewidth = 1) +
    geom_density(colour = fill, fill = NA, linewidth = 1) +
    theme_minimal() +
    labs(y = "Density")

  apply_quantity_x_scale(p, quantity = quantity, xlab = xlab)
}

#' Assign staggered y positions for cohort markers
#' @keywords internal
#' @noRd
stagger_cohort_y_positions <- function(
  cohort_df,
  draws_vec,
  y_top_frac = 0.92,
  y_bottom_frac = 0.48
) {
  dens <- tryCatch(stats::density(draws_vec), error = function(e) NULL)
  y_max <- if (!is.null(dens)) max(dens$y, na.rm = TRUE) else 0.1
  n <- nrow(cohort_df)
  if (n == 0L) {
    return(cohort_df)
  }

  cohort_df$y_max_dens <- y_max
  cohort_df$y <- if (n == 1L) {
    y_max * ((y_top_frac + y_bottom_frac) / 2)
  } else {
    seq(
      y_max * y_top_frac,
      y_max * y_bottom_frac,
      length.out = n
    )
  }
  cohort_df$label_y <- cohort_df$y + y_max * 0.1
  cohort_df
}

#' Add cohort overlays to an HCA density plot
#' @keywords internal
#' @noRd
add_cohort_overlay <- function(
  plot,
  cohort_df,
  draws_vec,
  show_se = TRUE,
  stagger_heights = TRUE,
  colour_by = c("cohort", "direction"),
  cohort_palette = NULL,
  direction_palette = c(
    above_hca = "#D62728",
    below_hca = "#1F77B4",
    consistent_with_hca = "#7F7F7F"
  ),
  annotate = c("group", "p_value")
) {
  colour_by <- match.arg(colour_by)
  annotate <- intersect(
    annotate,
    c("group", "p_value", "direction", "empirical_rank", "method")
  )

  if (isTRUE(stagger_heights)) {
    cohort_df <- stagger_cohort_y_positions(cohort_df, draws_vec)
  } else {
    dens <- tryCatch(stats::density(draws_vec), error = function(e) NULL)
    y_mark <- if (!is.null(dens)) max(dens$y, na.rm = TRUE) * 0.88 else 0.1
    cohort_df$y <- y_mark
    cohort_df$y_max_dens <- y_mark / 0.88
    cohort_df$label_y <- y_mark + cohort_df$y_max_dens * 0.1
  }

  colour_col <- if (identical(colour_by, "direction") && any(!is.na(cohort_df$direction))) {
    "direction"
  } else {
    "group"
  }

  if (!is.null(cohort_palette) && identical(colour_col, "group")) {
    plot <- plot +
      geom_point(
        data = cohort_df,
        aes(x = .data$log_mu, y = .data$y, colour = .data[[colour_col]]),
        size = 2.8,
        inherit.aes = FALSE
      ) +
      scale_colour_manual(values = cohort_palette, name = "Cohort")
  } else if (identical(colour_col, "direction")) {
    plot <- plot +
      geom_point(
        data = cohort_df,
        aes(x = .data$log_mu, y = .data$y, colour = .data$direction),
        size = 2.8,
        inherit.aes = FALSE
      ) +
      scale_colour_manual(values = direction_palette, name = "Direction")
  } else {
    plot <- plot +
      geom_point(
        data = cohort_df,
        aes(x = .data$log_mu, y = .data$y, colour = .data$group),
        size = 2.8,
        inherit.aes = FALSE
      ) +
      scale_colour_brewer(palette = "Dark2", name = "Cohort")
  }

  if (isTRUE(show_se) && any(is.finite(cohort_df$se))) {
    se_df <- cohort_df[is.finite(cohort_df$se) & cohort_df$se > 0, , drop = FALSE]
    if (nrow(se_df) > 0L) {
      plot <- plot +
        geom_linerange(
          data = se_df,
          aes(
            xmin = .data$log_mu - .data$se,
            xmax = .data$log_mu + .data$se,
            y = .data$y,
            colour = .data[[colour_col]]
          ),
          linewidth = 1,
          inherit.aes = FALSE
        )
    }
  }

  label_parts <- lapply(seq_len(nrow(cohort_df)), function(i) {
    bits <- character(0)
    if ("group" %in% annotate) {
      bits <- c(bits, cohort_df$group[[i]])
    }
    if ("p_value" %in% annotate && is.finite(cohort_df$p_value[[i]])) {
      bits <- c(bits, paste0("p=", signif(cohort_df$p_value[[i]], 2)))
    }
    if ("direction" %in% annotate && !is.na(cohort_df$direction[[i]])) {
      bits <- c(bits, cohort_df$direction[[i]])
    }
    if ("empirical_rank" %in% annotate && is.finite(cohort_df$empirical_rank[[i]])) {
      bits <- c(bits, paste0("rank=", signif(cohort_df$empirical_rank[[i]], 2)))
    }
    if ("method" %in% annotate && !is.na(cohort_df$method[[i]])) {
      bits <- c(bits, cohort_df$method[[i]])
    }
    paste(bits, collapse = "\n")
  })
  cohort_df$label <- unlist(label_parts, use.names = FALSE)

  y_ceiling <- max(c(cohort_df$label_y, cohort_df$y_max_dens), na.rm = TRUE)

  plot +
    geom_text(
      data = cohort_df,
      aes(
        x = .data$log_mu,
        y = .data$label_y,
        label = .data$label,
        colour = .data[[colour_col]]
      ),
      size = 2.9,
      lineheight = 0.9,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    coord_cartesian(clip = "off", ylim = c(0, y_ceiling * 1.12)) +
    theme(
      plot.margin = margin(12, 12, 12, 12)
    )
}

#' Density plot of healthy HCA posterior draws
#'
#' Visualises posterior draws from [expression_draws()]. Optional
#' `query_mu` / `query_SE` / `query_label` overlay cohort point estimates
#' on the same density.
#'
#' @param draws Posterior draws: numeric vector, or list from
#'   [expression_draws()].
#' @param quantity Quantity for bare numeric `draws`. Otherwise inferred from
#'   `draws$quantity`.
#' @param baseline_label Caption label for the HCA density curve.
#' @param fill Fill colour for the HCA density.
#' @param title Plot title. Default is built from gene / cell-type metadata.
#' @param subtitle Optional subtitle.
#' @param query_mu Optional numeric vector of cohort log(mu) estimates to
#'   overlay. Requires `quantity = "linpred"`.
#' @param query_SE Optional numeric vector of SEs for `query_mu` (same length,
#'   or length 1 recycled). Drawn as horizontal error bars when finite.
#' @param query_label Optional character labels for `query_mu` (same length,
#'   or length 1 recycled). Defaults to `"query"`, `"query2"`, ...
#' @return A `ggplot` object.
#' @export
#' @import ggplot2
#' @importFrom cli cli_abort
#' @importFrom scales log1p_trans
plot_hca_draws <- function(
  draws,
  quantity = NULL,
  baseline_label = "Healthy HCA",
  fill = "#4C78A8",
  title = NULL,
  subtitle = NULL,
  query_mu = NULL,
  query_SE = NULL,
  query_label = NULL
) {
  norm <- normalize_draws_input(draws, quantity = quantity)

  if (!is.null(query_mu)) {
    if (!identical(norm$quantity, "linpred")) {
      cli::cli_abort(
        "`query_mu` overlays require `quantity = \"linpred\"` draws."
      )
    }
  }

  p <- build_hca_density_plot(
    draws_vec = norm$draws,
    quantity = norm$quantity,
    fill = fill
  )

  if (!is.null(query_mu)) {
    query_mu <- as.numeric(query_mu)
    n <- length(query_mu)
    if (n < 1L || any(!is.finite(query_mu))) {
      cli::cli_abort("`query_mu` must be a non-empty numeric vector of finite values.")
    }

    if (is.null(query_label)) {
      query_label <- if (n == 1L) "query" else paste0("query", seq_len(n))
    } else {
      query_label <- as.character(query_label)
      if (length(query_label) == 1L && n > 1L) {
        query_label <- rep(query_label, n)
      }
      if (length(query_label) != n) {
        cli::cli_abort("`query_label` must have length 1 or match `query_mu`.")
      }
    }

    if (is.null(query_SE)) {
      se <- rep(NA_real_, n)
    } else {
      se <- as.numeric(query_SE)
      if (length(se) == 1L && n > 1L) {
        se <- rep(se, n)
      }
      if (length(se) != n) {
        cli::cli_abort("`query_SE` must have length 1 or match `query_mu`.")
      }
    }

    cohort_df <- data.frame(
      group = query_label,
      log_mu = query_mu,
      se = se,
      method = NA_character_,
      direction = NA_character_,
      p_value = NA_real_,
      empirical_rank = NA_real_,
      stringsAsFactors = FALSE
    )
    p <- add_cohort_overlay(
      plot = p,
      cohort_df = cohort_df,
      draws_vec = norm$draws,
      show_se = TRUE,
      stagger_heights = TRUE,
      colour_by = "cohort",
      annotate = "group"
    )
  }

  if (is.null(title)) {
    title <- default_hca_draws_title(
      norm,
      norm$quantity,
      comparison = !is.null(query_mu)
    )
  }

  p + labs(title = title, subtitle = subtitle, caption = baseline_label)
}

#' Arcsine square-root transform for composition proportions
#' @keywords internal
#' @noRd
arcsine_sqrt_trans <- function() {
  scales::trans_new(
    "arcsine_sqrt",
    transform = function(x) asin(sqrt(pmax(pmin(x, 1), 0))),
    inverse = function(x) (sin(x))^2
  )
}

#' Density plot of healthy composition posterior draws
#'
#' @param draws Output of [composition_draws()], or its `$draws` data frame.
#' @param cell_types Optional cell types to keep.
#' @param fill Fill colour for the density curves.
#' @param title,subtitle Optional plot labels.
#' @return A `ggplot` object.
#' @export
#' @import ggplot2
plot_composition_draws <- function(
  draws,
  cell_types = NULL,
  fill = "#4C78A8",
  title = NULL,
  subtitle = NULL
) {
  draws <- composition_draws_df(draws)
  if (!is.null(cell_types)) {
    draws <- draws[draws$cell_type %in% cell_types, , drop = FALSE]
  }

  ggplot(draws, aes(x = .data$proportion)) +
    geom_density(alpha = 0.35, fill = fill, colour = NA) +
    facet_wrap(~cell_type) +
    scale_x_continuous(trans = arcsine_sqrt_trans()) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Proportion (arcsine-sqrt)",
      y = "Density",
      caption = "Healthy HCA composition"
    ) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold"))
}

#' Plot observed proportions against healthy composition draws
#'
#' Accepts [composition_draws()] output plus either [composition_test()]
#' results or raw observed proportions.
#'
#' @param draws Output of [composition_draws()], or its `$draws` data frame.
#' @param test_results Optional data frame from [composition_test()].
#' @param proportions Optional observed proportions (used when `test_results`
#'   is not supplied).
#' @param annotate If `TRUE`, label observed values with empirical confidence.
#' @param title,subtitle Optional plot labels.
#' @return A `ggplot` object.
#' @export
#' @import ggplot2
plot_composition_vs_hca <- function(
  draws,
  test_results = NULL,
  proportions = NULL,
  annotate = TRUE,
  title = NULL,
  subtitle = NULL
) {
  if (is.null(test_results) && is.null(proportions)) {
    cli::cli_abort("Provide `test_results` and/or `proportions`.")
  }

  if (is.null(test_results)) {
    test_results <- composition_test(proportions, draws)
  }

  cell_types <- unique(test_results$cell_type)
  p <- plot_composition_draws(
    draws,
    cell_types = cell_types,
    title = title,
    subtitle = subtitle
  )

  p <- p +
    geom_vline(
      data = test_results,
      aes(
        xintercept = .data$proportion_observed,
        colour = .data$sample_id
      ),
      linetype = "dashed",
      linewidth = 0.6
    ) +
    labs(colour = "Sample") +
    theme(legend.position = "bottom")

  if (isTRUE(annotate) && "empirical_confidence" %in% names(test_results)) {
    p <- p +
      ggrepel::geom_text_repel(
        data = test_results,
        aes(
          x = .data$proportion_observed,
          y = Inf,
          label = paste0("EC:", signif(.data$empirical_confidence, 2)),
          colour = .data$sample_id
        ),
        size = 3.5,
        direction = "y",
        segment.color = NA,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }

  p
}
