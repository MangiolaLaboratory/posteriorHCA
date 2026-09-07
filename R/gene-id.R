# Gene identifier utilities (canonical ENSG for atlas models)

#' @keywords internal
#' @noRd
is_ensembl_gene_id <- function(x) {
  grepl("^ENSG[0-9]+$", x, ignore.case = TRUE)
}

#' @keywords internal
#' @noRd
strip_ensembl_version <- function(x) {
  sub("\\..*$", "", as.character(x))
}

#' Download the ENSG gene universe for a cell type
#'
#' @inheritParams get_brms_ready
#' @return Character vector of Ensembl gene ids.
#' @export
#' @importFrom cli cli_abort
get_gene_universe <- function(
  cell_type,
  version = "latest",
  cache_directory = get_default_cache_dir(),
  use_cache = TRUE
) {
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
    filename = "genes.csv"
  )
  if (!identical(res$status, "success")) {
    extra <- if (!is.null(res$error)) res$error else res$status
    cli_abort("Failed to retrieve gene universe: {extra}")
  }

  tbl <- readr::read_csv(res$path, show_col_types = FALSE)
  col <- if (".feature" %in% names(tbl)) {
    ".feature"
  } else if ("gene_ensg" %in% names(tbl)) {
    "gene_ensg"
  } else {
    names(tbl)[[1]]
  }
  sort(unique(strip_ensembl_version(as.character(tbl[[col]]))))
}

#' @keywords internal
#' @noRd
is_expression_fit <- function(x) {
  inherits(x, "posteriorHCA_expr_fit")
}

#' @keywords internal
#' @noRd
as_brms_fit <- function(x) {
  if (inherits(x, "brmsfit")) {
    return(x)
  }
  if (is_expression_fit(x)) {
    return(x$fit)
  }
  x
}

#' @keywords internal
#' @noRd
new_expression_fit <- function(fit, cell_type, gene_ensg) {
  structure(
    list(
      fit = fit,
      cell_type = as.character(cell_type),
      gene_ensg = as.character(gene_ensg)
    ),
    class = c("posteriorHCA_expr_fit", "list")
  )
}
