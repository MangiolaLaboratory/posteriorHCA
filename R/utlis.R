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
  model_url <- "https://object-store.rc.nectar.org.au/v1/AUTH_b0a86a29c8b74630aac35f471cfe1396/sccomp_est/estimates_age_decade___L3___disease_TRUE___immune_only_TRUE.rds"
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

#' @title Get file from Nectar object storage
#' @description Downloads and caches a file from Nectar object storage using direct URL downloads.
#' The file is synchronized to the cache directory maintaining the same structure as in the object store.
#' @param cache_directory Character. Directory to store the cached file. Defaults to the standard cache directory.
#' @param use_cache Logical. If TRUE, uses the cached version if available.
#' @param container Character. Name of the container in Nectar object storage.
#' @param prefix Character. Prefix (folder name) in the container.
#' @param filename Character. Name of the file to download.
#' @return List. Always returns a list with 'status' (success/not_found/error/empty), 'path' (local file path if successful, NULL otherwise), 'url' (requested URL), and optionally 'error' (error message for error cases).
#' @export
#' @importFrom cli cli_alert_info cli_alert_warning cli_abort
#' @importFrom tools file_path_sans_ext
#' @importFrom httr GET write_disk stop_for_status
get_file_ready <- function(cache_directory = get_default_cache_dir(), 
                          use_cache = TRUE,
                          container,
                          prefix,
                          filename) {
  
  # Validate inputs
  if (missing(container) || is.null(container) || container == "") {
    cli_abort("Container name must be specified")
  }
  if (missing(prefix) || is.null(prefix) || prefix == "") {
    cli_abort("Prefix must be specified")
  }
  if (missing(filename) || is.null(filename) || filename == "") {
    cli_abort("Filename must be specified")
  }
  
  # Construct the object path in the storage
  object_path <- paste(prefix, filename, sep = "/")
  
  # Construct the direct download URL (auth token is fixed)
  auth_token <- "AUTH_b0a86a29c8b74630aac35f471cfe1396"
  base_url <- "https://object-store.rc.nectar.org.au/v1"
  file_url <- paste(base_url, auth_token, container, object_path, sep = "/")
  
  # Construct the local file path maintaining the same structure
  local_dir <- file.path(cache_directory, container, prefix)
  local_path <- file.path(local_dir, filename)
  
  # Check if file exists and use_cache is TRUE
  if (file.exists(local_path) && use_cache) {
    cli_alert_info("Using cached file: {local_path}")
    return(list(status = "success", path = local_path, url = file_url))
  }
  
  # Create local directory structure
  dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
  
  cli_alert_info("Downloading file from Nectar object storage...")
  cli_alert_info("URL: {file_url}")
  cli_alert_info("Local path: {local_path}")
  
  # Download using direct URL (no Swift CLI needed)
  result <- tryCatch(
    {
      response <- GET(file_url, write_disk(local_path, overwrite = TRUE))
      
      # Check HTTP status code
      if (response$status_code == 404) {
        cli_alert_info("File not found: {file_url}")
        return(list(status = "not_found", path = NULL, url = file_url))
      }
      
      # Check for other HTTP errors
      if (response$status_code >= 400) {
        cli_alert_info("HTTP error {response$status_code}: {file_url}")
        return(list(status = "error", path = NULL, url = file_url, error = paste("HTTP", response$status_code)))
      }
      
      # Check if download was successful
      if (!file.exists(local_path) || file.size(local_path) == 0) {
        cli_alert_info("Downloaded file is empty or does not exist: {local_path}")
        return(list(status = "empty", path = local_path, url = file_url))
      }
      
      cli_alert_info("Successfully downloaded file: {local_path}")
      return(list(status = "success", path = local_path, url = file_url))
    },
    error = function(e) {
      # Clean up failed download
      if (file.exists(local_path)) {
        file.remove(local_path)
      }
      cli_alert_info("Download failed: {e$message}")
      return(list(status = "error", path = NULL, url = file_url, error = e$message))
    }
  )
  
  # Return the result (always a list with status and path)
  return(result)
}
