#' Extract sub-formula from full formula
#'
#' This function extracts a sub-formula that contains only the specified factor names.
#'
#' @param full_formula A full formula object.
#' @param factor_names A vector of factor names to retain.
#' @return A sub-formula as a string.
#' @export
#' @import purrr
get_sub_formula <-
  function(full_formula, factor_names){

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
  packageName() |>
    R_user_dir(
      "cache"
    ) |>
    normalizePath() |>
    suppressWarnings()
}

#' @title Get the file path for the pre-trained model
#' @description Checks if the pre-trained model file is cached locally. If not, downloads it from the cloud.
#' @param cache_directory Character. Directory to store the cached file. Defaults to the standard cache directory.
#' @param use_cache Logical. If TRUE, uses the cached version if available.
#' @return Character. The local file path of the model.
#' @export
#' @importFrom tools file_path_sans_ext
#' @importFrom httr GET write_disk stop_for_status
#' @importFrom cli cli_alert_info cli_abort
get_model_ready <- function(cache_directory = get_default_cache_dir(), use_cache = TRUE) {
  model_url <- "https://object-store.rc.nectar.org.au/v1/AUTH_06d6e008e3e642da99d806ba3ea629c5/testing/estimates_age_bins___L2.rds"
  model_filename <- basename(model_url)
  model_path <- file.path(cache_directory, model_filename)

  if (file.exists(model_path) && use_cache) {
    cli_alert_info("Using cached model file: {model_path}")
    return(model_path)
  }

  cli_alert_info("Downloading pre-trained model file...")
  dir.create(cache_directory, recursive = TRUE, showWarnings = FALSE)

  tryCatch(
    {
      GET(model_url, write_disk(model_path, overwrite = TRUE)) |> stop_for_status()
    },
    error = function(e) {
      file.remove(model_path)
      cli_abort("Failed to download model file: {e}")
    }
  )

  return(model_path)
}
