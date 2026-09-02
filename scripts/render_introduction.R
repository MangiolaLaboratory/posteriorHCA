#' Render Introduction vignette and sync README.md
#'
#' `rmarkdown::render()` does not execute the custom `knit:` hook in
#' `vignettes/Introduction.Rmd`, so use this helper when you want both the
#' vignette HTML and an updated package `README.md`.
#'
#' @param input Path to `Introduction.Rmd`. Defaults to the package vignette.
#' @param proj_root Package root directory. Auto-detected when `NULL`.
#' @return Invisibly returns the path to `README.md`.
#' @export
render_introduction <- function(
  input = NULL,
  proj_root = NULL
) {
  if (is.null(proj_root)) {
    proj_root <- find_posteriorhca_root()
  }
  proj_root <- normalizePath(proj_root, mustWork = TRUE)

  if (is.null(input)) {
    input <- file.path(proj_root, "vignettes", "Introduction.Rmd")
  }
  input <- normalizePath(input, mustWork = TRUE)

  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(proj_root, quiet = TRUE)
  } else {
    requireNamespace("posteriorHCA", quietly = FALSE)
  }

  rmarkdown::render(
    input,
    output_format = "rmarkdown::html_vignette",
    knit_root_dir = proj_root,
    params = list(demo_metadata = TRUE, for_github_md = FALSE),
    quiet = FALSE
  )

  rmarkdown::render(
    input,
    output_file = "README.md",
    output_format = "github_document",
    output_dir = proj_root,
    knit_root_dir = proj_root,
    params = list(demo_metadata = FALSE, for_github_md = TRUE),
    quiet = FALSE
  )

  invisible(file.path(proj_root, "README.md"))
}

find_posteriorhca_root <- function() {
  if (requireNamespace("rprojroot", quietly = TRUE)) {
    return(normalizePath(rprojroot::find_package_root_file(), mustWork = TRUE))
  }
  dir <- normalizePath(getwd(), mustWork = FALSE)
  for (candidate in c(dir, dirname(dir))) {
    if (file.exists(file.path(candidate, "DESCRIPTION"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  cli::cli_abort("Could not find package root (no DESCRIPTION file).")
}

if (identical(sys.nframe(), 0L)) {
  render_introduction()
}
