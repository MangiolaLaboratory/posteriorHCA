# Expression-model queries
#
#   Core:
#     build_newdata_grid() / expression_draws() / load_expression_fit()
#   Wrapper (compose cores only):
#     expression_baseline_draws()
#
# Metadata choices for each covariate:
#   NA / NULL     -> expand over every level the model knows
#   "A"           -> fix that level
#   c("A", "B")   -> expand only over those levels
#   str_detect("10x"), between(6, 8), ...
#                 -> filter levels. The covariate is inserted as the first
#                    argument, then evaluated like dplyr::filter().
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

#' Decide whether a captured metadata expression is "use every level"
#' @keywords internal
#' @noRd
metadata_quo_is_blank <- function(quo) {
  if (rlang::quo_is_missing(quo) || rlang::quo_is_null(quo)) {
    return(TRUE)
  }
  expr <- rlang::quo_get_expr(quo)
  if (is.symbol(expr)) {
    return(as.character(expr) %in% c(
      "NA", "NA_character_", "NA_integer_", "NA_real_", "NA_complex_", "NULL"
    ))
  }
  is.atomic(expr) && length(expr) == 1L && is.na(expr)
}

#' A value the caller already computed, as opposed to a filter expression
#' @keywords internal
#' @noRd
is_concrete_metadata <- function(x) {
  if (is.null(x) || is.factor(x)) {
    return(TRUE)
  }
  is.atomic(x) && !is.language(x)
}

#' Insert the covariate as the first argument of a filter call
#' @keywords internal
#' @noRd
metadata_predicate_expr <- function(expr, arg_name) {
  if (!is.call(expr)) {
    cli::cli_abort(
      "`{arg_name}` must be a value or a filter expression such as `str_detect(...)`."
    )
  }
  if (arg_name %in% all.names(expr)) {
    return(expr)
  }
  as.call(c(list(expr[[1]]), as.symbol(arg_name), as.list(expr)[-1]))
}

#' Keep model levels selected by an exact value or a filter expression
#'
#' Exact values (`"female"`, `c("4", "7")`) are matched literally. A call such
#' as `str_detect("10x")` or `between(6, 8)` is applied like
#' `dplyr::filter()` on the covariate's levels, with the covariate inserted
#' as the first argument. `str_detect()` is case-sensitive, so `"10x"` matches
#' `"10x Genomics 3"` and `"10X"` does not; use `"(?i)10X"` for either case.
#' @keywords internal
#' @noRd
resolve_metadata_levels <- function(quo, levels, arg_name) {
  blank <- list(chosen = levels, missing = TRUE)
  if (metadata_quo_is_blank(quo)) {
    return(blank)
  }

  evaluated <- tryCatch(
    rlang::eval_tidy(quo),
    error = function(e) e
  )
  eval_failed <- inherits(evaluated, "error")
  if (!eval_failed && is_concrete_metadata(evaluated)) {
    if (is.null(evaluated) || (length(evaluated) == 1L && is.na(evaluated))) {
      return(blank)
    }
    chosen <- as.character(evaluated)
    chosen <- chosen[!is.na(chosen) & nzchar(chosen)]
    if (length(chosen) == 0L) {
      return(list(chosen = levels, missing = FALSE))
    }
    unknown <- setdiff(chosen, levels)
    if (length(unknown) > 0L) {
      cli::cli_abort(c(
        "Unknown {arg_name} value{?s}: {unknown}.",
        "i" = "This model includes: {levels}."
      ))
    }
    return(list(chosen = unique(chosen), missing = FALSE))
  }

  expr <- rlang::quo_get_expr(quo)
  if (!eval_failed || !is.call(expr)) {
    cli::cli_abort(
      "`{arg_name}` must be a value or a filter expression such as `str_detect(...)`."
    )
  }

  predicate <- metadata_predicate_expr(expr, arg_name)
  filter_levels <- as.character(levels)
  numeric_levels <- suppressWarnings(as.numeric(filter_levels))
  column <- if (length(filter_levels) > 0L && !anyNA(numeric_levels)) {
    numeric_levels
  } else {
    filter_levels
  }
  data <- stats::setNames(data.frame(column, stringsAsFactors = FALSE), arg_name)
  helper_env <- rlang::new_environment(
    list(str_detect = stringr::str_detect, between = dplyr::between),
    parent = rlang::quo_get_env(quo)
  )
  keep <- tryCatch(
    rlang::eval_tidy(
      rlang::new_quosure(predicate, env = helper_env),
      data = data
    ),
    error = function(e) {
      label <- rlang::expr_deparse(expr)
      cli::cli_abort(c(
        "Could not apply the `{arg_name}` filter `{label}`.",
        "x" = conditionMessage(e)
      ))
    }
  )
  if (!is.logical(keep) || length(keep) != length(levels)) {
    cli::cli_abort(
      "The `{arg_name}` filter must return one TRUE or FALSE per level."
    )
  }
  chosen <- levels[which(keep & !is.na(keep))]
  if (length(chosen) == 0L) {
    cli::cli_abort(c(
      "No `{arg_name}` levels matched that expression.",
      "i" = "This model includes: {levels}."
    ))
  }
  list(chosen = unique(as.character(chosen)), missing = FALSE)
}

#' Build a covariate grid for an expression-model query
#'
#' Each metadata argument is one of:
#' * `NA` — marginalise over every level in the model
#' * a single value — fix that level
#' * several values — marginalise over only those levels
#' * a filter expression such as `str_detect("10x")` or `between(6, 8)` —
#'   keep the levels for which the expression is true. The covariate is
#'   inserted as the first argument and evaluated like `dplyr::filter()`.
#'   `str_detect()` is case-sensitive: `str_detect("10x")` selects
#'   `"10x Genomics 3"` and `"10x Genomics 5"`; `str_detect("10X")` selects
#'   neither. `str_detect("(?i)10X")` ignores case.
#'
#' `dataset_id` is a random intercept. `NA` stamps a new unseen-study label.
#'
#' @param fit A `brmsfit` or [load_expression_fit()] object.
#' @param age_decade,sex,disease_groups,ethnicity_groups,assay_groups,tissue_groups
#'   Level, vector of levels, or a filter expression. See the rules above.
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
  .build_newdata_grid(
    fit = fit,
    metadata = list(
      age_decade = rlang::enquo(age_decade),
      sex = rlang::enquo(sex),
      disease_groups = rlang::enquo(disease_groups),
      ethnicity_groups = rlang::enquo(ethnicity_groups),
      assay_groups = rlang::enquo(assay_groups),
      tissue_groups = rlang::enquo(tissue_groups)
    ),
    dataset_id = dataset_id,
    offset = offset,
    new_study_id = new_study_id
  )
}

#' @keywords internal
#' @noRd
.build_newdata_grid <- function(fit, metadata, dataset_id, offset, new_study_id) {
  fit <- as_brms_fit(fit)
  if (is.null(fit$data)) {
    cli_abort("`fit` must have a `$data` element.")
  }

  expand_list <- list()
  available_by_col <- list()
  fixed_vars <- character(0)
  marginalised_vars <- character(0)

  for (var in names(metadata)) {
    column <- find_fit_column(fit, .expression_meta_map[[var]])
    quo <- metadata[[var]]

    if (is.na(column)) {
      if (!metadata_quo_is_blank(quo)) {
        cli_abort("This model has no `{var}` column, but `{var}` was specified.")
      }
      next
    }

    available <- expression_model_levels(fit, var)
    resolved <- resolve_metadata_levels(quo, available, var)
    chosen <- resolved$chosen
    missing_choice <- resolved$missing

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
#' @param marginalise How to reduce several `newdata` rows (covariate
#'   profiles) to one draw vector. In the Bayesian sense, leaving covariates
#'   free in [build_newdata_grid()] expands a grid over those levels; this
#'   argument then **marginalises** over that grid so the returned draws
#'   represent a baseline that is not conditioned on a single profile.
#'   Options:
#'   \describe{
#'     \item{`"mean"`}{For each posterior draw, average the predicted
#'       quantity across grid rows. Interprets every profile as equally
#'       weighted; typical default for a single comparable baseline.}
#'     \item{`"pool"`}{Concatenate (stack) predictions from all grid rows.
#'       Keeps between-profile spread in the returned vector.}
#'     \item{`"sample"`}{For each posterior draw, pick one grid row at
#'       random. Approximates a discrete mixture over profiles.}
#'   }
#'   When `newdata` has a single row, all three options are equivalent.
#' @param ndraws Number of posterior draws, or `NULL` for all.
#' @param transform Passed to `posterior_linpred` only.
#' @param re_formula,allow_new_levels,sample_new_levels Passed to brms.
#' @param seed Optional RNG seed.
#' @return A list with `draws`, `grid`, `quantity`, `marginalise`, `n_grid`,
#'   `cell_type`, and `gene_ensg`.
#' @seealso [expression_baseline_draws()], [build_newdata_grid()],
#'   [load_expression_fit()]
#' @export
#' @importFrom cli cli_abort
expression_draws <- function(
  fit,
  newdata,
  quantity = c("linpred", "predict", "epred"),
  marginalise = c("mean", "pool", "sample"),
  ndraws = NULL,
  transform = FALSE,
  re_formula = NULL,
  allow_new_levels = TRUE,
  sample_new_levels = "gaussian",
  seed = NULL
) {
  quantity <- match.arg(quantity)
  marginalise <- match.arg(marginalise)
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
  } else if (marginalise == "mean") {
    draws <- as.numeric(rowMeans(draw_matrix))
  } else if (marginalise == "pool") {
    draws <- as.numeric(draw_matrix)
  } else {
    idx <- sample.int(ncol(draw_matrix), nrow(draw_matrix), replace = TRUE)
    draws <- as.numeric(draw_matrix[cbind(seq_len(nrow(draw_matrix)), idx)])
  }

  list(
    draws = draws,
    grid = newdata,
    quantity = quantity,
    marginalise = marginalise,
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

#' Posterior draws for an HCA expression baseline query
#'
#' Convenience wrapper around [load_expression_fit()],
#' [build_newdata_grid()], and [expression_draws()]. Loads the atlas model
#' when `fit` is not supplied, builds the covariate grid, and returns the
#' usual [expression_draws()] list. Does not call brms posterior helpers
#' directly and does not summarise or test draws.
#'
#' @param cell_type Cell type name (required when `fit` is `NULL`).
#' @param gene_ensg Ensembl gene id (required when `fit` is `NULL`).
#' @param fit Optional pre-loaded fit from [load_expression_fit()].
#' @param age_decade,sex,disease_groups,ethnicity_groups,assay_groups,tissue_groups
#'   Passed to [build_newdata_grid()], including filter expressions such as
#'   `str_detect("10x")` and `between(6, 8)`.
#' @param dataset_id,offset,new_study_id Passed to [build_newdata_grid()].
#' @param quantity,marginalise,ndraws,transform,re_formula,allow_new_levels,sample_new_levels,seed
#'   Passed to [expression_draws()].
#' @inheritParams load_expression_fit
#' @return A list from [expression_draws()].
#' @seealso [load_expression_fit()], [build_newdata_grid()],
#'   [expression_draws()], [compare_cohort_to_hca()]
#' @export
#' @importFrom cli cli_abort
expression_baseline_draws <- function(
  cell_type = NULL,
  gene_ensg = NULL,
  fit = NULL,
  age_decade = NA,
  sex = NA,
  disease_groups = NA,
  ethnicity_groups = NA,
  assay_groups = NA,
  tissue_groups = NA,
  dataset_id = NA,
  offset = 0,
  new_study_id = "__new_study__",
  quantity = c("linpred", "predict", "epred"),
  marginalise = c("mean", "pool", "sample"),
  ndraws = NULL,
  transform = FALSE,
  re_formula = NULL,
  allow_new_levels = TRUE,
  sample_new_levels = "gaussian",
  seed = NULL,
  version = "latest",
  cache_directory = get_default_cache_dir(),
  use_cache = TRUE
) {
  quantity <- match.arg(quantity)
  marginalise <- match.arg(marginalise)

  if (is.null(fit)) {
    if (is.null(cell_type) || is.null(gene_ensg)) {
      cli_abort("Provide `fit`, or both `cell_type` and `gene_ensg`.")
    }
    fit <- load_expression_fit(
      cell_type = cell_type,
      gene_ensg = gene_ensg,
      version = version,
      cache_directory = cache_directory,
      use_cache = use_cache
    )
  }

  newdata <- .build_newdata_grid(
    fit = fit,
    metadata = list(
      age_decade = rlang::enquo(age_decade),
      sex = rlang::enquo(sex),
      disease_groups = rlang::enquo(disease_groups),
      ethnicity_groups = rlang::enquo(ethnicity_groups),
      assay_groups = rlang::enquo(assay_groups),
      tissue_groups = rlang::enquo(tissue_groups)
    ),
    dataset_id = dataset_id,
    offset = offset,
    new_study_id = new_study_id
  )

  expression_draws(
    fit,
    newdata = newdata,
    quantity = quantity,
    marginalise = marginalise,
    ndraws = ndraws,
    transform = transform,
    re_formula = re_formula,
    allow_new_levels = allow_new_levels,
    sample_new_levels = sample_new_levels,
    seed = seed
  )
}
