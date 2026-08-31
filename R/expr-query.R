# Expression-model queries
#
# User metadata choices for each covariate:
#   NA            -> expand over every level the model knows
#   "A"           -> fix that level (error if it is not in the model)
#   c("A", "B")   -> expand only over those levels (each must be in the model)
#
# dataset_id is different: it is a random intercept for the training studies,
# so NA means "new unseen study", not "average every training dataset".

.expr_meta_map <- list(
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

#' Find the first matching column name in a fitted model's data
#' @keywords internal
#' @noRd
find_fit_column <- function(fit, candidates) {
  nms <- names(fit$data)
  found <- candidates[candidates %in% nms]
  if (length(found) == 0L) {
    return(NA_character_)
  }
  found[[1]]
}

#' Levels stored in a data column (keeps factor order when present)
#' @keywords internal
#' @noRd
levels_from_data_column <- function(x) {
  if (is.factor(x)) {
    lv <- levels(x)
  } else {
    lv <- sort(unique(as.character(x)))
  }
  lv[!is.na(lv) & nzchar(lv)]
}

#' Dummy-variable levels encoded in coefficient names
#'
#' brms names look like `assay_groups_altered10x Genomics 5` (no separator).
#' The reference level does *not* appear here; it only lives in the data.
#' @keywords internal
#' @noRd
levels_from_coef_names <- function(coef_names, column) {
  if (length(coef_names) == 0L || !nzchar(column)) {
    return(character(0))
  }
  keep <- startsWith(coef_names, column) & coef_names != column & !grepl(":", coef_names, fixed = TRUE)
  sub(column, "", coef_names[keep], fixed = TRUE)
}

#' Coefficient names from a brms fit, if available
#' @keywords internal
#' @noRd
coef_names_from_fit <- function(fit) {
  nms <- tryCatch(rownames(brms::fixef(fit)), error = function(e) character(0))
  if (length(nms) == 0L && !is.null(fit$parnames)) {
    nms <- as.character(fit$parnames)
  }
  nms
}

is_missing_choice <- function(x) {
  is.null(x) || (length(x) == 1L && is.na(x))
}

#' Allowed levels for one covariate in an expression model
#'
#' Combines training-data levels (includes the reference category) with
#' dummy names from the fixed-effect coefficients.
#'
#' @param fit A `brmsfit` (or a list with `$data` for testing).
#' @param variable User-facing name (`assay_groups`, `sex`, ...) or the
#'   actual column name in `fit$data`.
#' @return Character vector of allowed levels.
#' @export
expr_model_levels <- function(fit, variable) {
  if (is.null(fit$data)) {
    cli::cli_abort("`fit` must have a `$data` element.")
  }

  if (variable %in% names(.expr_meta_map)) {
    column <- find_fit_column(fit, .expr_meta_map[[variable]])
  } else if (variable %in% names(fit$data)) {
    column <- variable
  } else {
    column <- NA_character_
  }

  if (is.na(column)) {
    cli::cli_abort(
      "This model has no covariate matching `{variable}`."
    )
  }

  from_data <- levels_from_data_column(fit$data[[column]])
  if (is.factor(fit$data[[column]]) && length(from_data) > 0L) {
    return(from_data)
  }
  from_coef <- levels_from_coef_names(coef_names_from_fit(fit), column)
  unique(c(from_data, from_coef))
}

#' Resolve one user choice to a vector of levels
#' @keywords internal
#' @noRd
resolve_choice <- function(user_choice, available, arg_name) {
  if (is_missing_choice(user_choice)) {
    return(available)
  }

  chosen <- as.character(user_choice)
  chosen <- chosen[!is.na(chosen) & nzchar(chosen)]
  if (length(chosen) == 0L) {
    return(available)
  }

  unknown <- setdiff(chosen, available)
  if (length(unknown) > 0L) {
    cli::cli_abort(c(
      "Unknown {arg_name} value{?s}: {unknown}.",
      "i" = "This model includes: {available}."
    ))
  }
  chosen
}

#' Build a covariate grid for an expression-model query
#'
#' Each metadata argument is one of:
#' * `NA` — marginalise over every level in the model
#' * a single value — fix that level (must appear in the model)
#' * several values — marginalise over only those levels
#'
#' `dataset_id` is a random intercept. `NA` (the default) stamps a new
#' unseen-study label rather than expanding over training datasets.
#'
#' @param fit A `brmsfit` with `$data`.
#' @param age_decade,sex,disease_groups,ethnicity_groups,assay_groups,tissue_groups
#'   See the rules above.
#' @param dataset_id Existing dataset id, or `NA` for a new study.
#' @param offset Numeric offset. Atlas queries usually use `0`.
#' @param new_study_id Label used when `dataset_id` is `NA`.
#' @return A data frame, one row per covariate profile. Attributes
#'   `fixed` and `marginalised` name the user-facing variables.
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
  if (is.null(fit$data)) {
    cli_abort("`fit` must have a `$data` element.")
  }
  if (length(offset) != 1L || !is.numeric(offset) || is.na(offset)) {
    cli_abort("`offset` must be a single numeric value.")
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
  column_by_var <- list()
  available_by_col <- list()
  fixed_vars <- character(0)
  marginalised_vars <- character(0)

  for (var in names(user_choices)) {
    column <- find_fit_column(fit, .expr_meta_map[[var]])
    choice <- user_choices[[var]]

    if (is.na(column)) {
      if (!is_missing_choice(choice)) {
        cli_abort(
          "This model has no `{var}` column, but `{var}` was specified."
        )
      }
      next
    }

    available <- expr_model_levels(fit, var)
    if (length(available) == 0L) {
      cli_abort("Model column `{column}` has no usable levels.")
    }

    chosen <- resolve_choice(choice, available, var)
    expand_list[[column]] <- chosen
    column_by_var[[var]] <- column
    available_by_col[[column]] <- available

    if (is_missing_choice(choice) || length(chosen) > 1L) {
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
    if (is_missing_choice(dataset_id)) {
      grid[[dataset_col]] <- new_study_id
    } else {
      grid[[dataset_col]] <- as.character(dataset_id[[1]])
    }
  }

  # Dummy response keeps some brms predict paths happy.
  if ("counts" %in% names(fit$data) && !"counts" %in% names(grid)) {
    grid$counts <- 1
  }

  for (column in names(available_by_col)) {
    available <- available_by_col[[column]]
    if (is.factor(fit$data[[column]])) {
      grid[[column]] <- factor(
        as.character(grid[[column]]),
        levels = available
      )
    } else {
      grid[[column]] <- as.character(grid[[column]])
    }
  }

  grid <- as.data.frame(grid)
  attr(grid, "fixed") <- fixed_vars
  attr(grid, "marginalised") <- marginalised_vars

  n_profile <- nrow(grid)
  if (n_profile == 1L) {
    cli_alert_info("Covariate grid has 1 profile (all metadata fixed).")
  } else {
    marg <- paste(marginalised_vars, collapse = ", ")
    cli_alert_info(
      "Covariate grid has {n_profile} profile{?s} (marginalising: {marg})."
    )
  }

  grid
}

#' Collapse a draws-by-profile matrix over the covariate grid
#'
#' @param y_mat Matrix: rows = posterior draws, columns = grid rows.
#' @param method `"mean"` averages profiles within each draw (grand mean).
#'   `"pool"` stacks every draw x profile (mixture / density). `"sample"`
#'   picks one profile at random for each draw.
#' @return A numeric vector.
#' @export
marginalize_draw_matrix <- function(y_mat, method = c("mean", "pool", "sample")) {
  method <- match.arg(method)

  if (is.null(dim(y_mat))) {
    return(as.numeric(y_mat))
  }
  if (ncol(y_mat) <= 1L) {
    return(as.numeric(y_mat))
  }

  switch(
    method,
    mean = as.numeric(rowMeans(y_mat)),
    pool = as.numeric(y_mat),
    sample = {
      n <- nrow(y_mat)
      idx <- sample.int(ncol(y_mat), n, replace = TRUE)
      as.numeric(y_mat[cbind(seq_len(n), idx)])
    }
  )
}

#' Posterior draws of gene expression from a fitted brms model
#'
#' @param fit A `brmsfit`.
#' @param newdata Covariate grid from [build_newdata_grid()].
#' @param quantity `"linpred"` is log(μ) (`posterior_linpred`,
#'   `transform = FALSE`). `"predict"` is posterior predicted counts.
#'   `"epred"` is the expected count (`posterior_epred`).
#' @param collapse How to combine several grid rows. Use `"mean"` for
#'   location tests; `"pool"` for densities of predicted counts.
#' @param ndraws Number of posterior draws, or `NULL` for all.
#' @param transform Passed to `posterior_linpred` only. Default `FALSE`
#'   keeps log(μ). Set `TRUE` for μ on the count-mean scale.
#' @param re_formula,allow_new_levels,sample_new_levels Passed to brms.
#' @param seed Optional RNG seed.
#' @return A list with `draws`, `grid`, `quantity`, `collapse`, and
#'   `n_grid`.
#' @export
#' @importFrom cli cli_abort
expr_draws <- function(
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

  if (is.null(newdata) || nrow(newdata) < 1L) {
    cli_abort("`newdata` must have at least one row.")
  }
  newdata <- as.data.frame(newdata)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  pred_args <- list(
    object = fit,
    newdata = newdata,
    summary = FALSE,
    re_formula = re_formula,
    allow_new_levels = allow_new_levels,
    sample_new_levels = sample_new_levels,
    ndraws = ndraws
  )

  y_mat <- switch(
    quantity,
    linpred = do.call(
      brms::posterior_linpred,
      c(pred_args, list(transform = transform))
    ),
    predict = do.call(brms::posterior_predict, pred_args),
    epred = do.call(brms::posterior_epred, pred_args)
  )

  draws <- marginalize_draw_matrix(y_mat, method = collapse)

  list(
    draws = draws,
    grid = newdata,
    quantity = quantity,
    collapse = collapse,
    n_grid = nrow(newdata)
  )
}

#' Load a stored gene-level brms fit
#'
#' @inheritParams get_brms_ready
#' @return A `brmsfit`.
#' @export
#' @importFrom qs2 qs_read
#' @importFrom cli cli_abort
load_expr_fit <- function(
  cell_type,
  gene_ensg,
  version = "latest",
  cache_directory = get_default_cache_dir(),
  use_cache = TRUE
) {
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

  obj$brms_fit[[1]]
}
