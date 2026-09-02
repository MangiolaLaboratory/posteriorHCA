#' Extract sub-formula from full formula
#'
#' Deprecated. sccomp now marginalises over `NA` metadata directly, so
#' sub-formulas are no longer needed for composition queries.
#'
#' @param full_formula A full formula object.
#' @param factor_names A vector of factor names to retain.
#' @return A sub-formula as a string.
#' @export
#' @import purrr
get_sub_formula <-
  function(full_formula, factor_names){
    .Deprecated(msg = "get_sub_formula() is deprecated because sccomp now marginalises NA metadata directly.")

    terms <-
      full_formula %>%
      terms %>%
      attr('term.labels')

    return(

      terms[
        terms %>% map_lgl(
          .f = function(terms, factor_names){
            ((terms %>%
               strsplit("[+:\\|0-9]") %>%
               unlist %>% trimws %>% unique %>%
               setdiff('')) %in% factor_names) %>% all
          },
          factor_names = factor_names
        )
      ] %>% map_chr(.f = function(f){
        if (grepl(pattern = '\\|', f)){
          return(
            paste('(', f, ')')
          )
        }else{
          return(f)
        }
      }) %>%
        append(
        paste(
          '~',
          full_formula %>% terms %>% attr('intercept')
        ),
        .
      ) %>% paste(collapse = ' + ')

    )

  }

#' Returns the default cache directory with a version number
#' @export
#' @return A length one character vector.
#' @importFrom tools R_user_dir
#' @importFrom utils packageName
get_default_cache_dir <- function() {
  pkg <- tryCatch(utils::packageName(), error = function(e) "posteriorHCA")
  if (is.null(pkg) || !nzchar(pkg)) {
    pkg <- "posteriorHCA"
  }
  tools::R_user_dir(pkg, "cache") |>
    normalizePath() |>
    suppressWarnings()
}

.nectar_auth_token <- "AUTH_b0a86a29c8b74630aac35f471cfe1396"
.nectar_base_url <- "https://object-store.rc.nectar.org.au/v1"

#' Build a Nectar object URL
#'
#' @param container Character. Swift container (`V1`, `meta`, `sccomp_est`, ...).
#' @param ... Path segments joined with `/` (prefix, filename, ...). Empty
#'   segments are dropped.
#' @return Character scalar URL.
#' @keywords internal
#' @noRd
nectar_object_url <- function(container, ...) {
  segments <- c(.nectar_base_url, .nectar_auth_token, container, ...)
  segments <- segments[!is.na(segments) & nzchar(as.character(segments))]
  paste(segments, collapse = "/")
}

#' Normalise a version like `1`, `V1`, or `v1.1` to a container name.
#' @keywords internal
#' @noRd
nectar_container_from_version <- function(version) {
  version <- as.character(version)
  if (identical(version, "latest")) {
    return("latest")
  }
  paste0("V", sub("^[Vv]", "", version))
}

#' Download the Nectar model catalog (`latest.csv` or `all.csv`)
#'
#' `latest.csv` has one row per cell type pointing at the highest version
#' container. `all.csv` lists every cell type × version still published.
#'
#' @param which `"latest"` or `"all"`.
#' @param cache_directory Character. Local cache directory.
#' @param use_cache Logical. Reuse a previously downloaded catalog file.
#' @return A tibble with at least `container`, `ct_name`, and `ver_no`.
#' @export
#' @importFrom readr read_csv cols col_character col_integer
#' @importFrom cli cli_abort
get_nectar_catalog <- function(which = c("latest", "all"),
                               cache_directory = get_default_cache_dir(),
                               use_cache = TRUE) {
  which <- match.arg(which)
  filename <- paste0(which, ".csv")
  res <- get_file_ready(
    cache_directory = cache_directory,
    use_cache = use_cache,
    container = "meta",
    prefix = NULL,
    filename = filename
  )
  if (!identical(res$status, "success")) {
    extra <- if (!is.null(res$error)) res$error else res$status
    cli_abort("Failed to download Nectar catalog `{filename}`: {extra}")
  }
  readr::read_csv(
    res$path,
    show_col_types = FALSE,
    col_types = readr::cols(
      container = readr::col_character(),
      ct_name = readr::col_character(),
      ver_no = readr::col_character(),
      n_genes = readr::col_integer(),
      schema = readr::col_character(),
      updated_at = readr::col_character()
    )
  )
}

#' Resolve which Nectar container holds a cell type
#'
#' @param cell_type Character. Cell type name; passed through `make.names()`.
#' @param version `"latest"` (default), a container (`"V1"`), or a number (`"1"`).
#' @param catalog Optional catalog tibble. If `NULL`, downloaded from Nectar.
#' @param cache_directory,use_cache Passed to [get_nectar_catalog()] when
#'   `catalog` is `NULL`.
#' @return A named list: `cell_type`, `container`, `ver_no`, `n_genes`, `schema`.
#' @export
#' @importFrom cli cli_abort
lookup_cell_type_storage <- function(cell_type,
                                     version = "latest",
                                     catalog = NULL,
                                     cache_directory = get_default_cache_dir(),
                                     use_cache = TRUE) {
  ct <- make.names(cell_type)
  version <- as.character(version)

  if (is.null(catalog)) {
    which <- if (identical(version, "latest")) "latest" else "all"
    catalog <- get_nectar_catalog(
      which = which,
      cache_directory = cache_directory,
      use_cache = use_cache
    )
  }

  if (!all(c("ct_name", "container") %in% names(catalog))) {
    cli_abort("Catalog must have columns `ct_name` and `container`.")
  }

  catalog$ct_name <- as.character(catalog$ct_name)
  catalog$container <- as.character(catalog$container)
  if ("ver_no" %in% names(catalog)) {
    catalog$ver_no <- as.character(catalog$ver_no)
  }

  if (identical(version, "latest")) {
    row <- catalog[catalog$ct_name == ct, , drop = FALSE]
    if (nrow(row) > 1L && "ver_no" %in% names(row)) {
      row <- row[order(numeric_version(row$ver_no), decreasing = TRUE)[[1]], , drop = FALSE]
    }
  } else {
    container <- nectar_container_from_version(version)
    row <- catalog[catalog$ct_name == ct & catalog$container == container, , drop = FALSE]
  }

  if (nrow(row) == 0L) {
    available <- unique(catalog$ct_name)
    cli_abort(
      c(
        "Cell type `{cell_type}` (`{ct}`) is not in the Nectar catalog for version `{version}`.",
        "i" = "Available cell types: {paste(available, collapse = ', ')}"
      )
    )
  }

  list(
    cell_type = ct,
    container = row$container[[1]],
    ver_no = if ("ver_no" %in% names(row)) row$ver_no[[1]] else sub("^V", "", row$container[[1]]),
    n_genes = if ("n_genes" %in% names(row)) as.integer(row$n_genes[[1]]) else NA_integer_,
    schema = if ("schema" %in% names(row)) as.character(row$schema[[1]]) else NA_character_
  )
}

#' Download a gene-level brms fit for a cell type
#'
#' Looks up the version container from `meta/latest.csv` (or a pinned version)
#' and downloads `V{n}/{cell_type}/{gene_ensg}`.
#'
#' @param cell_type Character. Cell type name.
#' @param gene_ensg Character. Ensembl gene ID (object name, no extension).
#' @param version `"latest"` or a container/version pin (e.g. `"V1"`).
#' @param cache_directory Character. Local cache directory.
#' @param use_cache Logical. Reuse a previously downloaded file.
#' @return Same list as [get_file_ready()], plus `cell_type`, `container`, `ver_no`.
#' @export
get_brms_ready <- function(cell_type,
                           gene_ensg,
                           version = "latest",
                           cache_directory = get_default_cache_dir(),
                           use_cache = TRUE) {
  loc <- lookup_cell_type_storage(
    cell_type,
    version = version,
    cache_directory = cache_directory,
    use_cache = use_cache
  )
  res <- get_file_ready(
    cache_directory = cache_directory,
    use_cache = use_cache,
    container = loc$container,
    prefix = loc$cell_type,
    filename = gene_ensg
  )
  res$cell_type <- loc$cell_type
  res$container <- loc$container
  res$ver_no <- loc$ver_no
  res
}

#' Download the reference sample RDS for a cell type
#'
#' Object path: `{container}/reference_samples/{cell_type}.rds`.
#'
#' @inheritParams get_brms_ready
#' @return Same list as [get_file_ready()], plus `cell_type`, `container`, `ver_no`.
#' @export
get_reference_sample_ready <- function(cell_type,
                                       version = "latest",
                                       cache_directory = get_default_cache_dir(),
                                       use_cache = TRUE) {
  loc <- lookup_cell_type_storage(
    cell_type,
    version = version,
    cache_directory = cache_directory,
    use_cache = use_cache
  )
  res <- get_file_ready(
    cache_directory = cache_directory,
    use_cache = use_cache,
    container = loc$container,
    prefix = "reference_samples",
    filename = paste0(loc$cell_type, ".rds")
  )
  res$cell_type <- loc$cell_type
  res$container <- loc$container
  res$ver_no <- loc$ver_no
  res
}

#' @title Get the file path for the pre-trained model
#' @description
#' Deprecated in favour of [get_sccomp_ready()] or [load_sccomp_fit()].
#'
#' @param cache_directory Character. Directory to store the cached file. Defaults to the standard cache directory.
#' @param use_cache Logical. If TRUE, uses the cached version if available.
#' @return Character. The local file path of the model.
#' @export
#' @importFrom cli cli_alert_info cli_abort
get_model_ready <- function(cache_directory = get_default_cache_dir(), use_cache = TRUE) {
  .Deprecated(msg = "get_model_ready() is deprecated; use get_sccomp_ready() or load_sccomp_fit() instead.")
  res <- get_sccomp_ready(
    cache_directory = cache_directory,
    use_cache = use_cache
  )
  cli_alert_info("Using cached model file: {res$path}")
  res$path
}

#' @title Get file from Nectar object storage
#' @description Downloads and caches a file from Nectar object storage using direct URL downloads.
#' The file is synchronized to the cache directory maintaining the same structure as in the object store.
#' @param cache_directory Character. Directory to store the cached file. Defaults to the standard cache directory.
#' @param use_cache Logical. If TRUE, uses the cached version if available.
#' @param container Character. Name of the container in Nectar object storage.
#' @param prefix Character. Prefix (folder) in the container, e.g. a cell type
#'   (`cd8.naive`) or `reference_samples`. `NULL` or `""` for objects at the
#'   container root (e.g. `meta/latest.csv`).
#' @param filename Character. Name of the file to download.
#' @return List. Always returns a list with 'status' (success/not_found/error/empty), 'path' (local file path if successful, NULL otherwise), 'url' (requested URL), and optionally 'error' (error message for error cases).
#' @export
#' @importFrom cli cli_alert_info cli_alert_warning cli_abort
#' @importFrom httr GET write_disk
get_file_ready <- function(cache_directory = get_default_cache_dir(),
                          use_cache = TRUE,
                          container,
                          prefix = NULL,
                          filename) {

  if (missing(container) || is.null(container) || !nzchar(container)) {
    cli_abort("Container name must be specified")
  }
  if (missing(filename) || is.null(filename) || !nzchar(filename)) {
    cli_abort("Filename must be specified")
  }

  if (is.null(prefix) || !nzchar(prefix)) {
    prefix <- NULL
  }

  file_url <- nectar_object_url(container, prefix, filename)
  local_dir <- if (is.null(prefix)) {
    file.path(cache_directory, container)
  } else {
    file.path(cache_directory, container, prefix)
  }
  local_path <- file.path(local_dir, filename)

  if (file.exists(local_path) && isTRUE(use_cache) && isTRUE(file.size(local_path) > 0)) {
    cli_alert_info("Using cached file: {local_path}")
    return(list(status = "success", path = local_path, url = file_url))
  }

  dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)

  cli_alert_info("Downloading file from Nectar object storage...")
  cli_alert_info("URL: {file_url}")
  cli_alert_info("Local path: {local_path}")

  cleanup_local <- function() {
    if (file.exists(local_path)) {
      file.remove(local_path)
    }
  }

  result <- tryCatch(
    {
      response <- httr::GET(file_url, httr::write_disk(local_path, overwrite = TRUE))

      if (response$status_code == 404) {
        cleanup_local()
        cli_alert_info("File not found: {file_url}")
        return(list(status = "not_found", path = NULL, url = file_url))
      }

      if (response$status_code >= 400) {
        cleanup_local()
        cli_alert_info("HTTP error {response$status_code}: {file_url}")
        return(list(
          status = "error",
          path = NULL,
          url = file_url,
          error = paste("HTTP", response$status_code)
        ))
      }

      if (!file.exists(local_path) || file.size(local_path) == 0) {
        cleanup_local()
        cli_alert_info("Downloaded file is empty or does not exist: {local_path}")
        return(list(status = "empty", path = NULL, url = file_url))
      }

      cli_alert_info("Successfully downloaded file: {local_path}")
      return(list(status = "success", path = local_path, url = file_url))
    },
    error = function(e) {
      cleanup_local()
      cli_alert_info("Download failed: {e$message}")
      return(list(status = "error", path = NULL, url = file_url, error = e$message))
    }
  )

  return(result)
}
