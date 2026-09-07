# Expression-model queries
#
# Metadata choices for each covariate:
#   NA            -> expand over every level the model knows
#   "A"           -> fix that level
#   c("A", "B")   -> expand only over those levels
#
# dataset_id is a random intercept: NA means a new unseen study.

.expression_meta_map <- list(
  age_decade = c("age_decade"),
  sex = c("sex"),
  disease_groups = c(
    "disease_groups_altered",
    "disease_groups___altered",
    "disease_groups"
  ),
  ethnicity_groups = c("ethnicity_groups", "ethnicity_groups_imputed"),
  assay_groups = c(
    "assay_groups_altered",
    "assay_groups___altered",
    "assay_groups"
  ),
  tissue_groups = c("tissue_groups")
)

.dataset_col_candidates <- c("dataset_id_altered", "dataset_id")

#' First matching column name in fit$data
#' @keywords internal
#' @noRd
find_fit_column <- function(fit, candidates) {
  found <- candidates[candidates %in% names(fit$data)]
  if (length(found) == 0L) NA_character_ else found[[1]]
}

#' Allowed levels for one covariate in an expression model
#' @keywords internal
#' @noRd
expression_model_levels <- function(fit, variable) {
  fit <- as_brms_fit(fit)
  column <- if (variable %in% names(.expression_meta_map)) {
    find_fit_column(fit, .expression_meta_map[[variable]])
  } else if (variable %in% names(fit$data)) {
    variable
  } else {
    NA_character_
  }
  if (is.na(column)) {
    cli::cli_abort("This model has no covariate matching `{variable}`.")
  }

  x <- fit$data[[column]]
  from_data <- if (is.factor(x)) {
    levels(x)
  } else {
    sort(unique(as.character(x)))
  }
  from_data <- from_data[!is.na(from_data) & nzchar(from_data)]
  if (is.factor(x) && length(from_data) > 0L) {
    return(from_data)
  }

  coef_names <- tryCatch(rownames(brms::fixef(fit)), error = function(e) character(0))
  if (length(coef_names) == 0L && !is.null(fit$parnames)) {
    coef_names <- as.character(fit$parnames)
  }
  from_coef <- character(0)
  if (length(coef_names) && nzchar(column)) {
    keep <- startsWith(coef_names, column) &
      coef_names != column &
      !grepl(":", coef_names, fixed = TRUE)
    from_coef <- sub(column, "", coef_names[keep], fixed = TRUE)
  }
  unique(c(from_data, from_coef))
}

#' Build a covariate grid for an expression-model query
#'
#' Each metadata argument is one of:
#' * `NA` — marginalise over every level in the model
#' * a single value — fix that level
#' * several values — marginalise over only those levels
#'
#' `dataset_id` is a random intercept. `NA` stamps a new unseen-study label.
#'
#' @param fit A `brmsfit` or [load_expression_fit()] object.
#' @param age_decade,sex,disease_groups,ethnicity_groups,assay_groups,tissue_groups
#'   See the rules above.
#' @param dataset_id Existing dataset id, or `NA` for a new study.
#' @param offset Numeric offset. Atlas queries usually use `0`.
#' @param new_study_id Label used when `dataset_id` is `NA`.
#' @return A data frame, one row per covariate profile.
#' @export
#' @importFrom cli cli_abort cli_alert_info
build_newdata_grid <- function(
  fit,
  age_decade = NA,
  sex = NA,
  disease_groups = NA,
  ethnicity_groups = NA,
  assay_groups = NA,
  tissue_groups = NA,
  dataset_id = NA,
  offset = 0,
  new_study_id = "__new_study__"
) {
  fit <- as_brms_fit(fit)
  if (is.null(fit$data)) {
    cli_abort("`fit` must have a `$data` element.")
  }

  user_choices <- list(
    age_decade = age_decade,
    sex = sex,
    disease_groups = disease_groups,
    ethnicity_groups = ethnicity_groups,
    assay_groups = assay_groups,
    tissue_groups = tissue_groups
  )

  expand_list <- list()
  available_by_col <- list()
  fixed_vars <- character(0)
  marginalised_vars <- character(0)

  for (var in names(user_choices)) {
    column <- find_fit_column(fit, .expression_meta_map[[var]])
    choice <- user_choices[[var]]
    missing_choice <- is.null(choice) || (length(choice) == 1L && is.na(choice))

    if (is.na(column)) {
      if (!missing_choice) {
        cli_abort("This model has no `{var}` column, but `{var}` was specified.")
      }
      next
    }

    available <- expression_model_levels(fit, var)
    if (missing_choice) {
      chosen <- available
    } else {
      chosen <- as.character(choice)
      chosen <- chosen[!is.na(chosen) & nzchar(chosen)]
      if (length(chosen) == 0L) {
        chosen <- available
      }
      unknown <- setdiff(chosen, available)
      if (length(unknown) > 0L) {
        cli_abort(c(
          "Unknown {var} value{?s}: {unknown}.",
          "i" = "This model includes: {available}."
        ))
      }
    }

    expand_list[[column]] <- chosen
    available_by_col[[column]] <- available
    if (missing_choice || length(chosen) > 1L) {
      marginalised_vars <- c(marginalised_vars, var)
    } else {
      fixed_vars <- c(fixed_vars, var)
    }
  }

  if (length(expand_list) == 0L) {
    cli_abort("None of the expected covariates are present in `fit$data`.")
  }

  grid <- tidyr::expand_grid(!!!expand_list)
  grid$offset <- offset

  dataset_col <- find_fit_column(fit, .dataset_col_candidates)
  if (!is.na(dataset_col)) {
    if (is.null(dataset_id) || (length(dataset_id) == 1L && is.na(dataset_id))) {
      grid[[dataset_col]] <- new_study_id
    } else {
      grid[[dataset_col]] <- as.character(dataset_id[[1]])
    }
  }

  if ("counts" %in% names(fit$data) && !"counts" %in% names(grid)) {
    grid$counts <- 1
  }

  for (column in names(available_by_col)) {
    available <- available_by_col[[column]]
    if (is.factor(fit$data[[column]])) {
      grid[[column]] <- factor(as.character(grid[[column]]), levels = available)
    } else {
      grid[[column]] <- as.character(grid[[column]])
    }
  }

  grid <- as.data.frame(grid)
  attr(grid, "fixed") <- fixed_vars
  attr(grid, "marginalised") <- marginalised_vars

  if (nrow(grid) == 1L) {
    cli_alert_info("Covariate grid has 1 profile (all metadata fixed).")
  } else {
    cli_alert_info(
      "Covariate grid has {nrow(grid)} profile{?s} (marginalising: {paste(marginalised_vars, collapse = ', ')})."
    )
  }
  grid
}

#' Posterior draws of gene expression from a fitted brms model
#'
#' @param fit A `brmsfit` or [load_expression_fit()] object.
#' @param newdata Covariate grid from [build_newdata_grid()].
#' @param quantity `"linpred"` (log μ), `"predict"`, or `"epred"`.
#' @param collapse How to combine several grid rows: `"mean"`, `"pool"`,
#'   or `"sample"`.
#' @param ndraws Number of posterior draws, or `NULL` for all.
#' @param transform Passed to `posterior_linpred` only.
#' @param re_formula,allow_new_levels,sample_new_levels Passed to brms.
#' @param seed Optional RNG seed.
#' @return A list with `draws`, `grid`, `quantity`, `collapse`, `n_grid`,
#'   `cell_type`, and `gene_ensg`.
#' @export
#' @importFrom cli cli_abort
expression_draws <- function(
  fit,
  newdata,
  quantity = c("linpred", "predict", "epred"),
  collapse = c("mean", "pool", "sample"),
  ndraws = NULL,
  transform = FALSE,
  re_formula = NULL,
  allow_new_levels = TRUE,
  sample_new_levels = "gaussian",
  seed = NULL
) {
  quantity <- match.arg(quantity)
  collapse <- match.arg(collapse)
  brms_fit <- as_brms_fit(fit)
  newdata <- as.data.frame(newdata)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  pred_args <- list(
    object = brms_fit,
    newdata = newdata,
    summary = FALSE,
    re_formula = re_formula,
    allow_new_levels = allow_new_levels,
    sample_new_levels = sample_new_levels,
    ndraws = ndraws
  )

  draw_matrix <- switch(
    quantity,
    linpred = do.call(
      brms::posterior_linpred,
      c(pred_args, list(transform = transform))
    ),
    predict = do.call(brms::posterior_predict, pred_args),
    epred = do.call(brms::posterior_epred, pred_args)
  )

  if (is.null(dim(draw_matrix)) || ncol(draw_matrix) <= 1L) {
    draws <- as.numeric(draw_matrix)
  } else if (collapse == "mean") {
    draws <- as.numeric(rowMeans(draw_matrix))
  } else if (collapse == "pool") {
    draws <- as.numeric(draw_matrix)
  } else {
    idx <- sample.int(ncol(draw_matrix), nrow(draw_matrix), replace = TRUE)
    draws <- as.numeric(draw_matrix[cbind(seq_len(nrow(draw_matrix)), idx)])
  }

  list(
    draws = draws,
    grid = newdata,
    quantity = quantity,
    collapse = collapse,
    n_grid = nrow(newdata),
    cell_type = if (is_expression_fit(fit)) fit$cell_type else NA_character_,
    gene_ensg = if (is_expression_fit(fit)) fit$gene_ensg else NA_character_
  )
}

#' Load a stored gene-level brms fit
#'
#' @param cell_type Cell type name.
#' @param gene_ensg Ensembl gene id (for example `"ENSG00000169252"`).
#' @inheritParams get_brms_ready
#' @return A `posteriorHCA_expr_fit` list with `fit`, `cell_type`, and
#'   `gene_ensg`.
#' @export
#' @importFrom qs2 qs_read
#' @importFrom cli cli_abort
load_expression_fit <- function(
  cell_type,
  gene_ensg,
  version = "latest",
  cache_directory = get_default_cache_dir(),
  use_cache = TRUE
) {
  gene_ensg <- strip_ensembl_version(as.character(gene_ensg[[1]]))
  if (!is_ensembl_gene_id(gene_ensg)) {
    cli_abort("`gene_ensg` must be an Ensembl gene id (ENSG...).")
  }

  res <- get_brms_ready(
    cell_type = cell_type,
    gene_ensg = gene_ensg,
    version = version,
    cache_directory = cache_directory,
    use_cache = use_cache
  )
  if (!identical(res$status, "success")) {
    extra <- if (!is.null(res$error)) res$error else res$status
    cli_abort("Failed to retrieve brms fit: {extra}")
  }

  obj <- qs2::qs_read(res$path)
  if (!"brms_fit" %in% names(obj)) {
    cli_abort("Cached file does not contain a `brms_fit` column: {res$path}")
  }

  new_expression_fit(
    fit = obj$brms_fit[[1]],
    cell_type = res$cell_type,
    gene_ensg = gene_ensg
  )
}
