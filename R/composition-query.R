.sccomp_default_release <- "cellNexus_1_0_12"
.sccomp_default_model <- "estimates_age_decade___L3___disease_FALSE___immune_only_TRUE.rds"

.sccomp_user_meta_map <- c(
  age_decade = "age_decade",
  sex = "sex",
  ethnicity_groups = "ethnicity_groups_imputed",
  assay_groups = "assay_groups___altered",
  tissue_groups = "tissue_groups"
)

#' Download a pre-trained sccomp composition model from Nectar
#'
#' Default object path:
#' `sccomp_est/cellNexus_1_0_12/estimates_age_decade___L3___disease_FALSE___immune_only_TRUE.rds`
#'
#' @param model Character. RDS file name inside the release prefix.
#' @param release Character. Nectar prefix under the `sccomp_est` container.
#' @param cache_directory Character. Local cache directory.
#' @param use_cache Logical. Reuse a previously downloaded file.
#' @return Same list as [get_file_ready()], plus `release` and `model`.
#' @export
#' @importFrom cli cli_abort
get_sccomp_ready <- function(
  model = .sccomp_default_model,
  release = .sccomp_default_release,
  cache_directory = get_default_cache_dir(),
  use_cache = TRUE
) {
  res <- get_file_ready(
    cache_directory = cache_directory,
    use_cache = use_cache,
    container = "sccomp_est",
    prefix = release,
    filename = model
  )
  if (!identical(res$status, "success")) {
    extra <- if (!is.null(res$error)) res$error else res$status
    cli_abort("Failed to retrieve sccomp model `{model}`: {extra}")
  }
  c(res, list(release = release, model = model))
}

#' Extract the fitted sccomp object from a loader result
#' @keywords internal
#' @noRd
sccomp_fit_object <- function(x) {
  if (inherits(x, "posteriorHCA_sccomp_fit")) {
    return(x$fit)
  }
  x
}

#' Training metadata columns available in a sccomp fit
#' @keywords internal
#' @noRd
sccomp_count_data <- function(fit) {
  fit <- sccomp_fit_object(fit)
  cd <- attr(fit, "count_data")
  if (is.null(cd)) {
    cli_abort("sccomp fit has no `count_data` attribute.")
  }
  as.data.frame(cd)
}

#' Build sccomp query metadata
#'
#' Each metadata argument is `NA`/`NULL` (leave for sccomp to marginalise),
#' an exact level, or a dplyr filter expression such as `str_detect("10x")`.
#' Several matched levels expand to several rows.
#'
#' @keywords internal
#' @noRd
build_sccomp_newdata <- function(
  fit,
  sample_id = "query_sample",
  age_decade = NA,
  sex = NA,
  ethnicity_groups = NA,
  assay_groups = NA,
  tissue_groups = NA
) {
  fit <- sccomp_fit_object(fit)
  cd <- sccomp_count_data(fit)

  metadata <- list(
    age_decade = list(quo = rlang::enquo(age_decade), column = "age_decade"),
    sex = list(quo = rlang::enquo(sex), column = "sex"),
    ethnicity_groups = list(
      quo = rlang::enquo(ethnicity_groups),
      column = "ethnicity_groups_imputed"
    ),
    assay_groups = list(
      quo = rlang::enquo(assay_groups),
      column = "assay_groups___altered"
    ),
    tissue_groups = list(quo = rlang::enquo(tissue_groups), column = "tissue_groups")
  )

  expand_list <- list()
  for (var in names(metadata)) {
    column <- metadata[[var]]$column
    quo <- metadata[[var]]$quo
    if (!column %in% names(cd)) {
      cli_abort("This sccomp model has no `{column}` column.")
    }

    levels <- unique(as.character(cd[[column]]))
    levels <- levels[!is.na(levels) & nzchar(levels)]
    resolved <- pick_metadata_levels(quo, levels, var)
    # blank → NA for sccomp marginalisation (not expand-all)
    expand_list[[column]] <- if (resolved$NA_expand) NA_character_ else resolved$chosen
  }

  grid <- as.data.frame(tidyr::expand_grid(!!!expand_list))
  n <- nrow(grid)
  sid <- as.character(sample_id[[1]])
  grid$sample_id <- if (n <= 1L) sid else paste0(sid, "__", seq_len(n))
  grid$dataset_id___altered <- NA_character_
  grid[, c(
    "sample_id",
    "age_decade",
    "sex",
    "ethnicity_groups_imputed",
    "assay_groups___altered",
    "dataset_id___altered",
    "tissue_groups"
  )]
}

#' Load a pre-trained healthy-reference sccomp composition model
#'
#' Downloads and caches the default healthy-only CellNexus sccomp model from
#' Nectar, then returns a lightweight wrapper around the fitted object.
#'
#' @inheritParams get_sccomp_ready
#' @return A `posteriorHCA_sccomp_fit` object with elements `fit`, `release`,
#'   `model`, and `path`.
#' @export
#' @importFrom cli cli_abort
load_sccomp_fit <- function(
  model = .sccomp_default_model,
  release = .sccomp_default_release,
  cache_directory = get_default_cache_dir(),
  use_cache = TRUE
) {
  res <- get_sccomp_ready(
    model = model,
    release = release,
    cache_directory = cache_directory,
    use_cache = use_cache
  )
  fit <- readRDS(res$path)
  structure(
    list(
      fit = fit,
      release = res$release,
      model = res$model,
      path = res$path
    ),
    class = c("posteriorHCA_sccomp_fit", "list")
  )
}

#' Draw posterior predictive cell-type proportions from a sccomp fit
#'
#' Uses the sccomp model's full composition formula. Unknown covariates can be
#' left as `NA`; sccomp marginalises over them automatically.
#'
#' @param fit A `posteriorHCA_sccomp_fit` object from [load_sccomp_fit()], or a
#'   fitted sccomp object.
#' @param newdata Optional sample metadata with columns `sample_id`,
#'   `age_decade`, `sex`, `ethnicity_groups_imputed`, `assay_groups___altered`,
#'   `dataset_id___altered`, and `tissue_groups`. When omitted, a one-row table
#'   is built from the scalar metadata arguments below.
#' @param sample_id Sample identifier used when `newdata` is `NULL`.
#' @param age_decade,sex,ethnicity_groups,assay_groups,tissue_groups Metadata
#'   for the query profile: `NA` to marginalise, an exact level, or a dplyr
#'   filter expression such as `str_detect("10x")`.
#' @param summary_instead_of_draws Passed to [sccomp::sccomp_predict()].
#' @param ... Additional arguments passed to [sccomp::sccomp_predict()].
#' @return A list with `draws` (tibble of posterior predictive proportions),
#'   `newdata`, `release`, and `model`.
#' @export
#' @importFrom cli cli_abort
composition_draws <- function(
  fit,
  newdata = NULL,
  sample_id = "query_sample",
  age_decade = NA,
  sex = NA,
  ethnicity_groups = NA,
  assay_groups = NA,
  tissue_groups = NA,
  summary_instead_of_draws = FALSE,
  ...
) {
  if (inherits(fit, "posteriorHCA_sccomp_fit")) {
    release <- fit$release
    model <- fit$model
    fit_obj <- fit$fit
  } else {
    release <- .sccomp_default_release
    model <- .sccomp_default_model
    fit_obj <- sccomp_fit_object(fit)
  }

  if (is.null(newdata)) {
    newdata <- rlang::inject(build_sccomp_newdata(
      fit_obj,
      sample_id = sample_id,
      age_decade = !!rlang::enquo(age_decade),
      sex = !!rlang::enquo(sex),
      ethnicity_groups = !!rlang::enquo(ethnicity_groups),
      assay_groups = !!rlang::enquo(assay_groups),
      tissue_groups = !!rlang::enquo(tissue_groups)
    ))
  } else {
    newdata <- as.data.frame(newdata)
    if (!"sample_id" %in% names(newdata)) {
      cli_abort("`newdata` must contain a `sample_id` column.")
    }
    if (!"dataset_id___altered" %in% names(newdata)) {
      newdata$dataset_id___altered <- NA
    }
  }

  draws <- sccomp::sccomp_predict(
    fit_obj,
    new_data = newdata,
    summary_instead_of_draws = summary_instead_of_draws,
    ...
  )

  if ("L3" %in% names(draws) && !"cell_type" %in% names(draws)) {
    draws$cell_type <- draws$L3
  }

  list(
    draws = draws,
    newdata = newdata,
    release = release,
    model = model
  )
}
