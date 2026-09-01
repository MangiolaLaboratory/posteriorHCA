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

#' Resolve a single finite gene dispersion value
#' @keywords internal
#' @noRd
scalar_gene_dispersion <- function(dispersion, gene_ensg) {
  if (is.null(dispersion)) {
    return(NA_real_)
  }
  if (length(dispersion) == 1L && is.numeric(dispersion)) {
    return(as.numeric(dispersion))
  }
  if (!is.null(names(dispersion)) && gene_ensg %in% names(dispersion)) {
    return(as.numeric(dispersion[[gene_ensg]]))
  }
  as.numeric(dispersion[[1]])
}

#' Estimate negative-binomial dispersion for one gene
#' @keywords internal
#' @noRd
estimate_gene_dispersion <- function(source_obj, gene_ensg, offset = NULL, assay = NULL) {
  if (is_aligned_result(source_obj)) {
    aligned <- aligned_fields(source_obj)
    counts <- aligned$counts
    if (is.null(offset)) {
      offset <- aligned$offset
    }
  } else {
    counts <- extract_count_matrix(source_obj, assay = assay, arg_name = "counts")
    if (is.null(offset)) {
      fields <- aligned_fields(source_obj)
      if (!is.null(fields) && !is.null(fields$offset)) {
        offset <- fields$offset
      }
    }
    if (is.null(offset) && inherits(source_obj, "SummarizedExperiment")) {
      cd <- as.data.frame(SummarizedExperiment::colData(source_obj))
      if ("hca_offset" %in% names(cd)) {
        offset <- cd$hca_offset
        names(offset) <- rownames(cd)
      }
    }
    if (is.null(offset) && inherits(source_obj, "Seurat")) {
      meta <- tryCatch(source_obj[[]], error = function(e) NULL)
      if (!is.null(meta) && "hca_offset" %in% names(meta)) {
        offset <- meta$hca_offset
        names(offset) <- rownames(meta)
      }
    }
  }

  if (is.null(offset)) {
    cli_abort("`offset` is missing.")
  }
  if (!is.null(names(offset)) && all(names(offset) %in% colnames(counts))) {
    offset <- offset[colnames(counts)]
  }

  samples_df <- data.frame(
    sample_id = colnames(counts),
    group = factor(rep("cohort", ncol(counts))),
    stringsAsFactors = FALSE
  )
  dge <- edgeR::DGEList(counts = counts, samples = samples_df)
  dge$offset <- matrix(as.numeric(offset), nrow = nrow(dge), ncol = ncol(dge), byrow = TRUE)
  design <- stats::model.matrix(~ 1, data = samples_df)
  dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
  disp_vec <- dispersion_from_dge(dge)
  as.numeric(disp_vec[[gene_ensg]])
}

#' Bayesian bootstrap of cohort log(mu)
#'
#' Computes posterior draws of cohort latent log(mu) by repeatedly fitting
#' [edgeR::mglmOneGroup()] with Dirichlet weights.
#'
#' `counts` may be a matrix, `SummarizedExperiment`, `SingleCellExperiment`,
#' `Seurat` object, or an aligned list from [scale_to_hca_reference()].
#'
#' @param counts Count container or matrix.
#' @param offset Numeric vector of sample offsets, or `NULL` if provided in
#'   an aligned `counts` list.
#' @param dispersion Numeric scalar or named vector of negative-binomial
#'   dispersion values. If `NULL`, dispersion is estimated across `counts`.
#' @param group Optional group vector or column name to subset samples. If
#'   `NULL`, all non-reference samples (or all samples) are used.
#' @param gene Character scalar; gene id resolved to ENSG (must be in `rownames(counts)`).
#' @param n_boot Integer; number of bootstrap iterations (default 2000L).
#' @param seed Optional RNG seed for reproducibility.
#' @param assay Assay name for SE / Seurat input.
#' @return A numeric vector of length `n_boot` containing posterior log(mu) draws.
#' @export
#' @importFrom cli cli_abort
bootstrap_cohort_logmu <- function(
  counts,
  offset = NULL,
  dispersion = NULL,
  group = NULL,
  gene,
  n_boot = 2000L,
  seed = NULL,
  assay = NULL
) {
  if (missing(gene) || is.null(gene) || length(gene) != 1L) {
    cli_abort("`gene` must be a single gene identifier.")
  }
  n_boot <- as.integer(n_boot[[1L]])
  if (length(n_boot) != 1L || is.na(n_boot) || n_boot < 1L) {
    cli_abort("`n_boot` must be a positive integer.")
  }

  sample_role <- NULL
  source_obj <- counts
  if (is_aligned_result(counts)) {
    aligned <- aligned_fields(counts)
    if (is.null(offset)) {
      offset <- aligned$offset
    }
    sample_role <- aligned$sample_role
    counts <- aligned$counts
  } else {
    counts <- extract_count_matrix(source_obj, assay = assay, arg_name = "counts")
    fields <- aligned_fields(source_obj)
    if (!is.null(fields) && !is.null(fields$sample_role)) {
      sample_role <- fields$sample_role
    } else if (inherits(source_obj, "Seurat")) {
      meta <- tryCatch(source_obj[[]], error = function(e) NULL)
      if (!is.null(meta) && "sample_role" %in% names(meta)) {
        sample_role <- stats::setNames(meta$sample_role, rownames(meta))
      }
    } else if (inherits(source_obj, "SummarizedExperiment")) {
      cd <- as.data.frame(SummarizedExperiment::colData(source_obj))
      if ("sample_role" %in% names(cd)) {
        sample_role <- stats::setNames(cd$sample_role, rownames(cd))
      }
    }
  }

  if (as.character(gene) %in% rownames(counts)) {
    gene_ensg <- as.character(gene)
  } else {
    gene_ensg <- resolve_gene_one(gene, strict = TRUE)
  }

  if (!gene_ensg %in% rownames(counts)) {
    cli_abort("Gene `{gene_ensg}` not found in `counts`.")
  }

  if (is.null(offset)) {
    fields <- aligned_fields(source_obj)
    if (!is.null(fields) && !is.null(fields$offset)) {
      offset <- fields$offset
    }
  }
  if (is.null(offset) && inherits(source_obj, "SummarizedExperiment")) {
    cd <- as.data.frame(SummarizedExperiment::colData(source_obj))
    if ("hca_offset" %in% names(cd)) {
      offset <- cd$hca_offset
      names(offset) <- rownames(cd)
    }
  }
  if (is.null(offset) && inherits(source_obj, "Seurat")) {
    meta <- tryCatch(source_obj[[]], error = function(e) NULL)
    if (!is.null(meta) && "hca_offset" %in% names(meta)) {
      offset <- meta$hca_offset
      names(offset) <- rownames(meta)
    }
  }
  if (is.null(offset)) {
    cli_abort("`offset` is missing.")
  }
  if (length(offset) != ncol(counts)) {
    cli_abort("`offset` length ({length(offset)}) must match `ncol(counts)` ({ncol(counts)}).")
  }
  if (!is.null(names(offset)) && all(names(offset) %in% colnames(counts))) {
    offset <- offset[colnames(counts)]
  }

  cols_idx <- seq_len(ncol(counts))
  if (!is.null(group)) {
    resolved_group <- resolve_group(group, source_obj, ncol(counts), sample_role)
    cols_idx <- which(resolved_group == group[[1]] | resolved_group %in% group)
  } else if (!is.null(sample_role)) {
    cols_idx <- which(sample_role == "user")
  }

  if (length(cols_idx) == 0L) {
    cli_abort("No samples selected for cohort bootstrap.")
  }

  y <- counts[gene_ensg, cols_idx, drop = FALSE]
  y <- as.matrix(y)
  storage.mode(y) <- "double"
  if (anyNA(y)) {
    cli_abort("Gene `{gene_ensg}` has NA counts in the selected samples.")
  }
  o <- as.numeric(offset)[cols_idx]
  n_samples <- length(cols_idx)

  disp <- scalar_gene_dispersion(dispersion, gene_ensg)
  if (!is.finite(disp) || disp <= 0) {
    disp <- estimate_gene_dispersion(source_obj, gene_ensg, offset = offset, assay = assay)
  }
  if (!is.finite(disp) || disp <= 0) {
    cli_abort(
      "Could not resolve a finite positive dispersion for gene `{gene_ensg}`."
    )
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  draws <- vapply(seq_len(n_boot), function(i) {
    w <- draw_dirichlet_weights(n_samples)
    edgeR::mglmOneGroup(
      y,
      offset = o,
      dispersion = disp,
      weights = w
    )[[1]]
  }, numeric(1))

  draws
}

#' Summarise bootstrap log(mu) draws for one cohort
#' @keywords internal
#' @noRd
summarize_bootstrap_draws <- function(draws, hca_draws = NULL) {
  log_mu <- stats::median(draws)
  out <- list(
    log_mu = unname(log_mu),
    mu = unname(exp(log_mu)),
    se = unname(stats::sd(draws)),
    boot_q025 = unname(stats::quantile(draws, 0.025)),
    boot_q975 = unname(stats::quantile(draws, 0.975))
  )
  if (!is.null(hca_draws)) {
    if (is.list(hca_draws) && "draws" %in% names(hca_draws)) {
      hca_draws <- hca_draws$draws
    }
    out$empirical_rank <- mean(hca_draws <= log_mu)
  }
  out
}

#' Bayesian bootstrap of cohort log(mu) for multiple groups
#'
#' Mirrors the ergonomics of [estimate_cohort_logmu()] and
#' [welch_t_test_cohort_hca()]: one function call bootstraps every cohort
#' (except `exclude_groups`) and returns a data frame with `group`, `log_mu`,
#' and `se` columns compatible with downstream testing and plotting.
#'
#' @inheritParams estimate_cohort_logmu
#' @inheritParams bootstrap_cohort_logmu
#' @param cohorts Optional cohort labels to bootstrap. Defaults to all groups
#'   except `exclude_groups`.
#' @param exclude_groups Groups omitted when `cohorts` is `NULL` (default
#'   `"reference"`).
#' @param cohort_est Optional output from [estimate_cohort_logmu()]. When
#'   supplied, `genes` and `dispersion` can be taken from this table.
#' @param hca_draws Optional HCA posterior draws used to compute
#'   `empirical_rank` for each cohort.
#' @param keep_draws If `TRUE`, include a list-column `draws` with the raw
#'   bootstrap vectors.
#' @return A data frame with one row per cohort and columns compatible with
#'   [estimate_cohort_logmu()], plus `method = "bootstrap"`, `boot_q025`,
#'   `boot_q975`, and optionally `empirical_rank` and `draws`.
#' @export
#' @importFrom cli cli_abort
bootstrap_cohort_logmu_batch <- function(
  counts,
  group,
  genes = NULL,
  cohorts = NULL,
  exclude_groups = "reference",
  dispersion = NULL,
  cohort_est = NULL,
  hca_draws = NULL,
  n_boot = 2000L,
  seed = NULL,
  assay = NULL,
  keep_draws = FALSE
) {
  if (!is.null(cohort_est)) {
    if (!is.data.frame(cohort_est) || nrow(cohort_est) == 0L) {
      cli::cli_abort("`cohort_est` must be a non-empty data frame.")
    }
    if (is.null(genes) && "gene" %in% names(cohort_est)) {
      genes <- unique(as.character(cohort_est$gene))
    }
    if (is.null(cohorts) && "group" %in% names(cohort_est)) {
      cohorts <- setdiff(unique(as.character(cohort_est$group)), exclude_groups)
    }
  }

  if (is.null(genes) || length(genes) != 1L) {
    cli::cli_abort("`genes` must be a single gene identifier.")
  }
  n_boot <- as.integer(n_boot[[1L]])
  if (length(n_boot) != 1L || is.na(n_boot) || n_boot < 1L) {
    cli::cli_abort("`n_boot` must be a positive integer.")
  }

  mat_check <- extract_count_matrix(counts, assay = assay, arg_name = "counts")
  if (as.character(genes) %in% rownames(mat_check)) {
    gene_ensg <- as.character(genes)
    gene_symbol <- if (is_ensembl_gene_id(gene_ensg)) NA_character_ else as.character(genes)
  } else {
    gene_ensg <- resolve_gene_one(genes, strict = TRUE)
    gene_symbol <- if (is_ensembl_gene_id(genes)) NA_character_ else as.character(genes)
  }

  if (is.null(dispersion) && !is.null(cohort_est) && "dispersion" %in% names(cohort_est)) {
    disp_vals <- cohort_est$dispersion[cohort_est$gene == gene_ensg]
    disp_vals <- disp_vals[is.finite(disp_vals) & disp_vals > 0]
    if (length(disp_vals)) {
      dispersion <- as.numeric(disp_vals[[1]])
    }
  }

  n_lib <- ncol(extract_count_matrix(counts, assay = assay, arg_name = "counts"))
  sample_role <- NULL
  if (is_aligned_result(counts)) {
    fields <- aligned_fields(counts)
    sample_role <- fields$sample_role
  }
  resolved_group <- resolve_group(
    group,
    counts,
    n_lib,
    sample_role = sample_role,
    assay = assay
  )
  if (is.null(cohorts)) {
    cohorts <- setdiff(unique(as.character(resolved_group)), exclude_groups)
  } else {
    cohorts <- as.character(cohorts)
  }
  cohorts <- cohorts[!is.na(cohorts) & nzchar(cohorts)]
  if (length(cohorts) == 0L) {
    cli::cli_abort("No cohorts selected for bootstrap.")
  }

  if (!is.null(cohort_est) && "gene" %in% names(cohort_est)) {
    gene_ids <- unique(as.character(cohort_est$gene))
    if (length(gene_ids) > 1L) {
      cli::cli_abort("`cohort_est` contains multiple genes; subset to one gene.")
    }
  }

  if (is.null(dispersion) || !is.finite(scalar_gene_dispersion(dispersion, gene_ensg)) ||
      scalar_gene_dispersion(dispersion, gene_ensg) <= 0) {
    dispersion <- estimate_gene_dispersion(counts, gene_ensg, assay = assay)
  }

  out <- lapply(seq_along(cohorts), function(cohort_i) {
    cohort_label <- cohorts[[cohort_i]]
    sub_counts <- subset_aligned_by_group(
      counts,
      resolved_group,
      cohort_label,
      assay = assay
    )
    sub_meta <- extract_sample_metadata(sub_counts)
    n_user <- sum(sub_meta$sample_role == "user", na.rm = TRUE)
    if (n_user < 1L) {
      cli::cli_abort("Cohort `{cohort_label}` has no user samples for bootstrap.")
    }
    cohort_seed <- if (is.null(seed)) NULL else seed + cohort_i - 1L
    boot_draws <- bootstrap_cohort_logmu(
      sub_counts,
      gene = gene_ensg,
      dispersion = dispersion,
      n_boot = n_boot,
      seed = cohort_seed,
      assay = assay
    )
    summary <- summarize_bootstrap_draws(boot_draws, hca_draws = hca_draws)

    meta_row <- NULL
    if (!is.null(cohort_est)) {
      meta_row <- cohort_est[cohort_est$group == cohort_label, , drop = FALSE]
    }

    row <- data.frame(
      gene = gene_ensg,
      gene_symbol = if (!is.null(meta_row) && "gene_symbol" %in% names(meta_row) &&
        !is.na(meta_row$gene_symbol[[1]]) && nzchar(meta_row$gene_symbol[[1]])) {
        as.character(meta_row$gene_symbol[[1]])
      } else {
        gene_symbol
      },
      cell_type = if (!is.null(meta_row) && "cell_type" %in% names(meta_row)) {
        as.character(meta_row$cell_type[[1]])
      } else {
        NA_character_
      },
      group = cohort_label,
      n = n_user,
      log_mu = summary$log_mu,
      mu = summary$mu,
      se = summary$se,
      df = NA_real_,
      dispersion = if (!is.null(meta_row) && "dispersion" %in% names(meta_row)) {
        as.numeric(meta_row$dispersion[[1]])
      } else if (!is.null(dispersion)) {
        if (!is.null(names(dispersion)) && gene_ensg %in% names(dispersion)) {
          as.numeric(dispersion[[gene_ensg]])
        } else {
          as.numeric(dispersion[[1]])
        }
      } else {
        NA_real_
      },
      method = "bootstrap",
      boot_q025 = summary$boot_q025,
      boot_q975 = summary$boot_q975,
      empirical_rank = if (!is.null(summary$empirical_rank)) summary$empirical_rank else NA_real_,
      stringsAsFactors = FALSE
    )
    if (isTRUE(keep_draws)) {
      row$draws <- list(boot_draws)
    }
    row
  })

  do.call(rbind, out)
}

#' Welch t-test for one cohort row against HCA posterior draws
#' @keywords internal
#' @noRd
welch_t_test_one_cohort <- function(
  cohort_row,
  hca_draws,
  hca_meta,
  alpha = 0.05
) {
  if (!all(c("log_mu", "se") %in% names(cohort_row))) {
    cli_abort("`cohort_est` data frame must contain `log_mu` and `se` columns.")
  }

  cohort_meta <- cohort_est_metadata(cohort_row)
  validate_expr_metadata_match(cohort_meta, hca_meta, context = "`cohort_est` and `hca_draws`")

  hca_mean <- mean(hca_draws)
  hca_sd <- stats::sd(hca_draws)
  n_hca <- length(hca_draws)

  c_log_mu <- cohort_row$log_mu[[1]]
  c_se <- cohort_row$se[[1]]
  c_n <- if ("n" %in% names(cohort_row)) cohort_row$n[[1]] else 10L
  cohort_name <- if ("group" %in% names(cohort_row)) {
    as.character(cohort_row$group[[1]])
  } else {
    "cohort"
  }

  delta <- c_log_mu - hca_mean
  se_diff <- sqrt(c_se^2 + hca_sd^2)
  t_stat <- delta / se_diff

  var1_mean <- c_se^2
  var2_mean <- hca_sd^2

  denom <- (var1_mean^2 / max(1L, c_n - 1L)) + (var2_mean^2 / max(1L, n_hca - 1L))
  df <- if (denom > 0) {
    (var1_mean + var2_mean)^2 / denom
  } else {
    Inf
  }

  p_val <- if (is.finite(df) && df > 0) {
    2 * stats::pt(abs(t_stat), df = df, lower.tail = FALSE)
  } else {
    2 * stats::pnorm(-abs(t_stat))
  }

  empirical_rank <- mean(hca_draws <= c_log_mu)

  direction <- if (p_val < alpha && delta > 0) {
    "above_hca"
  } else if (p_val < alpha && delta < 0) {
    "below_hca"
  } else {
    "consistent_with_hca"
  }

  gene_ensg <- cohort_meta$gene_ensg
  if (is.na(gene_ensg) || !nzchar(gene_ensg)) {
    gene_ensg <- hca_meta$gene_ensg
  }
  gene_symbol <- cohort_meta$gene_symbol
  if (is.na(gene_symbol) || !nzchar(gene_symbol)) {
    gene_symbol <- hca_meta$gene_symbol
  }
  cell_type <- cohort_meta$cell_type
  if (is.na(cell_type) || !nzchar(cell_type)) {
    cell_type <- hca_meta$cell_type
  }

  data.frame(
    gene = if (!is.na(gene_ensg)) as.character(gene_ensg) else NA_character_,
    gene_symbol = if (!is.na(gene_symbol)) as.character(gene_symbol) else NA_character_,
    cell_type = if (!is.na(cell_type)) as.character(cell_type) else NA_character_,
    cohort = cohort_name,
    method = if ("method" %in% names(cohort_row)) {
      as.character(cohort_row$method[[1]])
    } else {
      "ql"
    },
    cohort_log_mu = c_log_mu,
    cohort_se = c_se,
    hca_mean = hca_mean,
    hca_sd = hca_sd,
    delta_log_mu = delta,
    se_diff = se_diff,
    t_stat = t_stat,
    df = df,
    p_value = p_val,
    empirical_rank = empirical_rank,
    direction = direction,
    stringsAsFactors = FALSE
  )
}

#' Welch t-test of cohort log(mu) against HCA posterior draws
#'
#' Compares cohort estimates from [estimate_cohort_logmu()] to healthy HCA
#' posterior draws using a heteroscedastic Welch t-test and empirical posterior
#' rank. By default, tests every cohort in `cohort_est` except `exclude_groups`.
#'
#' @param cohort_est A data frame from [estimate_cohort_logmu()] or
#'   [bootstrap_cohort_logmu_batch()] with columns `group`, `log_mu`, and `se`.
#' @param hca_draws A numeric vector of posterior draws, or the list returned
#'   by [expr_draws()] / [expr_predict()].
#' @param cohorts Optional character vector of cohort labels to test. Defaults
#'   to all groups in `cohort_est` except those in `exclude_groups`.
#' @param exclude_groups Groups omitted when `cohorts` is `NULL` (default
#'   `"reference"`).
#' @param alpha Significance level used for the `direction` column.
#' @return A data frame with one row per tested cohort:
#'   \item{gene}{Canonical ENSG id}
#'   \item{gene_symbol}{Gene symbol when available}
#'   \item{cell_type}{Cell type label when available}
#'   \item{cohort}{Cohort label}
#'   \item{cohort_log_mu}{Estimated cohort latent log(mu)}
#'   \item{cohort_se}{Cohort standard error}
#'   \item{hca_mean}{Mean of HCA posterior log(mu) draws}
#'   \item{hca_sd}{Standard deviation of HCA posterior log(mu) draws}
#'   \item{delta_log_mu}{Difference: `cohort_log_mu - hca_mean`}
#'   \item{se_diff}{Combined standard error: `sqrt(cohort_se^2 + hca_sd^2)`}
#'   \item{t_stat}{Welch-style t-statistic}
#'   \item{df}{Welch-Satterthwaite degrees of freedom}
#'   \item{p_value}{Two-sided p-value}
#'   \item{empirical_rank}{Proportion of HCA draws below `cohort_log_mu`}
#'   \item{direction}{`"above_hca"`, `"below_hca"`, or `"consistent_with_hca"`}
#' @export
#' @importFrom cli cli_abort
welch_t_test_cohort_hca <- function(
  cohort_est,
  hca_draws,
  cohorts = NULL,
  exclude_groups = "reference",
  alpha = 0.05
) {
  if (!is.data.frame(cohort_est) || nrow(cohort_est) == 0L) {
    cli_abort("`cohort_est` must be a non-empty data frame from [estimate_cohort_logmu()] or [bootstrap_cohort_logmu_batch()].")
  }
  if (!all(c("group", "log_mu", "se") %in% names(cohort_est))) {
    cli_abort("`cohort_est` must contain `group`, `log_mu`, and `se` columns.")
  }

  hca_meta <- expr_metadata(hca_draws)
  if (is.list(hca_draws) && !is.data.frame(hca_draws) && "draws" %in% names(hca_draws)) {
    hca_draws <- hca_draws$draws
  }
  if (!is.numeric(hca_draws) || length(hca_draws) < 2L) {
    cli_abort("`hca_draws` must be a numeric vector with at least 2 draws.")
  }

  if (is.null(cohorts)) {
    cohorts <- setdiff(unique(as.character(cohort_est$group)), exclude_groups)
  } else {
    cohorts <- as.character(cohorts)
  }
  cohorts <- cohorts[!is.na(cohorts) & nzchar(cohorts)]
  if (length(cohorts) == 0L) {
    cli_abort("No cohorts selected for testing.")
  }

  sub <- cohort_est[cohort_est$group %in% cohorts, , drop = FALSE]
  if (nrow(sub) == 0L) {
    cli_abort(
      "No rows in `cohort_est` match `cohorts`: {paste(cohorts, collapse = ', ')}."
    )
  }

  if ("gene" %in% names(sub)) {
    genes <- unique(as.character(sub$gene))
    if (length(genes) > 1L) {
      cli_abort(
        "`cohort_est` contains multiple genes; subset to one gene or provide one `hca_draws` object per gene."
      )
    }
  }

  validate_expr_metadata_match(cohort_est_metadata(sub[1, , drop = FALSE]), hca_meta)

  out <- lapply(seq_len(nrow(sub)), function(i) {
    welch_t_test_one_cohort(
      cohort_row = sub[i, , drop = FALSE],
      hca_draws = hca_draws,
      hca_meta = hca_meta,
      alpha = alpha
    )
  })

  do.call(rbind, out)
}
