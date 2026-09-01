# Plotting utilities for posteriorHCA expression outputs

#' Extract posterior draws and metadata from common result objects
#' @param x Output from [expr_draws()], [expr_predict()], or a numeric vector.
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
      gene_symbol = NA_character_,
      cell_type = NA_character_
    ))
  }

  if (!is.list(x) || is.data.frame(x)) {
    cli::cli_abort(
      "`draws` must be a numeric vector or a list from [expr_draws()] / [expr_predict()]."
    )
  }

  draws_vec <- if ("draws" %in% names(x)) {
    x$draws
  } else if ("pred" %in% names(x) && is.data.frame(x$pred) && "value" %in% names(x$pred)) {
    x$pred$value
  } else {
    cli::cli_abort(
      "Could not find posterior draws in `draws` (expected `$draws` or `$pred$value`)."
    )
  }

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

  meta <- expr_metadata(x)
  list(
    draws = as.numeric(draws_vec),
    quantity = qty,
    gene_ensg = meta$gene_ensg,
    gene_symbol = meta$gene_symbol,
    cell_type = meta$cell_type
  )
}

#' Normalise cohort estimates or test results for plotting
#' @keywords internal
#' @noRd
normalize_cohort_plot_df <- function(
  cohort_results,
  exclude_groups = character(0)
) {
  if (is.null(cohort_results)) {
    cli::cli_abort("`cohort_results` is missing.")
  }
  if (!is.data.frame(cohort_results) || nrow(cohort_results) == 0L) {
    cli::cli_abort("`cohort_results` must be a non-empty data frame.")
  }

  if (all(c("cohort", "cohort_log_mu") %in% names(cohort_results))) {
    out <- data.frame(
      group = as.character(cohort_results$cohort),
      log_mu = as.numeric(cohort_results$cohort_log_mu),
      se = if ("cohort_se" %in% names(cohort_results)) {
        as.numeric(cohort_results$cohort_se)
      } else {
        NA_real_
      },
      method = if ("method" %in% names(cohort_results)) {
        as.character(cohort_results$method)
      } else {
        "welch"
      },
      direction = if ("direction" %in% names(cohort_results)) {
        as.character(cohort_results$direction)
      } else {
        NA_character_
      },
      p_value = if ("p_value" %in% names(cohort_results)) {
        as.numeric(cohort_results$p_value)
      } else {
        NA_real_
      },
      empirical_rank = if ("empirical_rank" %in% names(cohort_results)) {
        as.numeric(cohort_results$empirical_rank)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  } else if (all(c("group", "log_mu") %in% names(cohort_results))) {
    out <- data.frame(
      group = as.character(cohort_results$group),
      log_mu = as.numeric(cohort_results$log_mu),
      se = if ("se" %in% names(cohort_results)) {
        as.numeric(cohort_results$se)
      } else {
        NA_real_
      },
      method = if ("method" %in% names(cohort_results)) {
        as.character(cohort_results$method)
      } else {
        "ql"
      },
      direction = if ("direction" %in% names(cohort_results)) {
        as.character(cohort_results$direction)
      } else {
        NA_character_
      },
      p_value = if ("p_value" %in% names(cohort_results)) {
        as.numeric(cohort_results$p_value)
      } else {
        NA_real_
      },
      empirical_rank = if ("empirical_rank" %in% names(cohort_results)) {
        as.numeric(cohort_results$empirical_rank)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  } else {
    cli::cli_abort(
      "`cohort_results` must be from [welch_t_test_cohort_hca()], [estimate_cohort_logmu()], or [bootstrap_cohort_logmu_batch()]."
    )
  }

  if (!is.null(exclude_groups) && length(exclude_groups)) {
    out <- out[!out$group %in% exclude_groups, , drop = FALSE]
  }

  if (nrow(out) == 0L) {
    cli::cli_abort("No cohort rows left to plot after applying `exclude_groups`.")
  }

  out
}

#' @rdname normalize_cohort_plot_df
#' @keywords internal
#' @noRd
normalize_test_results_plot_df <- normalize_cohort_plot_df

#' X-axis label for an expression-model quantity
#' @keywords internal
#' @noRd
expr_quantity_xlab <- function(quantity) {
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
  gene_lab <- if (!is.na(meta$gene_symbol) && nzchar(meta$gene_symbol)) {
    meta$gene_symbol
  } else if (!is.na(meta$gene_ensg) && nzchar(meta$gene_ensg)) {
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
  xlab <- expr_quantity_xlab(quantity)
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
#' Visualises posterior draws from [expr_draws()] or [expr_predict()]. For
#' cohort comparisons against a healthy baseline, use [plot_cohort_vs_hca()].
#'
#' @param draws Posterior draws: numeric vector, or list from [expr_draws()] /
#'   [expr_predict()].
#' @param quantity Quantity for bare numeric `draws`. Otherwise inferred from
#'   `draws$quantity`.
#' @param baseline_label Caption label for the HCA density curve.
#' @param fill Fill colour for the HCA density.
#' @param title Plot title. Default is built from gene / cell-type metadata.
#' @param subtitle Optional subtitle.
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
  subtitle = NULL
) {
  norm <- normalize_draws_input(draws, quantity = quantity)
  p <- build_hca_density_plot(
    draws_vec = norm$draws,
    quantity = norm$quantity,
    fill = fill
  )

  if (is.null(title)) {
    title <- default_hca_draws_title(norm, norm$quantity, comparison = FALSE)
  }

  p + labs(title = title, subtitle = subtitle, caption = baseline_label)
}

#' Plot cohort log(mu) estimates against healthy HCA posterior draws
#'
#' Visualises cohort log(mu) estimates against healthy HCA posterior draws.
#' Accepts outputs from [welch_t_test_cohort_hca()], [estimate_cohort_logmu()],
#' or [bootstrap_cohort_logmu_batch()].
#'
#' @param hca_draws Posterior draws from [expr_draws()] or [expr_predict()].
#'   Must use `quantity = "linpred"` because cohort tests are on log(mu).
#' @param test_results Optional data frame from [welch_t_test_cohort_hca()].
#' @param cohort_est Optional data frame from [estimate_cohort_logmu()] or
#'   [bootstrap_cohort_logmu_batch()]. Used when `test_results` is not supplied.
#' @param exclude_groups Character vector of cohort labels to omit.
#' @param show_se If `TRUE`, draw horizontal error bars for cohort `log_mu +/- se`
#'   (or `cohort_log_mu +/- cohort_se` for Welch output).
#' @param stagger_heights If `TRUE`, place each cohort marker at a different
#'   height above the density to reduce overlap.
#' @param colour_by Colour cohort markers by `"cohort"` or `"direction"`.
#' @param cohort_palette Optional colour vector when `colour_by = "cohort"`.
#' @param annotate Character vector of label fields: `"group"`, `"p_value"`,
#'   `"direction"`, `"empirical_rank"`, `"method"`.
#' @param fill Fill colour for the HCA density curve.
#' @param title Plot title. Default is built from gene / cell-type metadata.
#' @param subtitle Optional subtitle.
#' @return A `ggplot` object.
#' @export
#' @import ggplot2
#' @importFrom cli cli_abort
plot_cohort_vs_hca <- function(
  hca_draws,
  test_results = NULL,
  cohort_est = NULL,
  exclude_groups = character(0),
  show_se = TRUE,
  stagger_heights = TRUE,
  colour_by = c("cohort", "direction"),
  cohort_palette = NULL,
  annotate = c("group", "p_value"),
  fill = "#4C78A8",
  title = NULL,
  subtitle = NULL
) {
  if (is.null(test_results) && is.null(cohort_est)) {
    cli::cli_abort("Provide `test_results` and/or `cohort_est`.")
  }

  norm <- normalize_draws_input(hca_draws)
  if (!identical(norm$quantity, "linpred")) {
    cli::cli_abort(
      "`plot_cohort_vs_hca()` requires `quantity = \"linpred\"` draws from [expr_draws()]."
    )
  }

  cohort_input <- if (!is.null(test_results)) test_results else cohort_est
  cohort_df <- normalize_cohort_plot_df(
    cohort_input,
    exclude_groups = exclude_groups
  )

  p <- build_hca_density_plot(
    draws_vec = norm$draws,
    quantity = norm$quantity,
    fill = fill
  )
  p <- add_cohort_overlay(
    plot = p,
    cohort_df = cohort_df,
    draws_vec = norm$draws,
    show_se = show_se,
    stagger_heights = stagger_heights,
    colour_by = colour_by,
    cohort_palette = cohort_palette,
    annotate = annotate
  )

  if (is.null(title)) {
    title <- default_hca_draws_title(norm, norm$quantity, comparison = TRUE)
  }

  p + labs(
    title = title,
    subtitle = subtitle,
    caption = "Healthy HCA posterior with cohort estimates"
  )
}
