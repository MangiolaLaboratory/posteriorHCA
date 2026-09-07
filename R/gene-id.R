# Gene identifier utilities (canonical ENSG for atlas models)

#' Test whether a string looks like an Ensembl gene id
#' @param x Character vector.
#' @return Logical vector.
#' @keywords internal
#' @noRd
is_ensembl_gene_id <- function(x) {
  grepl("^ENSG[0-9]+$", x, ignore.case = TRUE)
}

#' Strip Ensembl version suffixes (e.g. `ENSG000001.1` -> `ENSG000001`)
#' @param x Character vector.
#' @return Character vector.
#' @keywords internal
#' @noRd
strip_ensembl_version <- function(x) {
  sub("\\..*$", "", as.character(x))
}

#' Download the ENSG gene universe for a cell type
#'
#' Reads `genes.csv` from Nectar object storage (`{container}/{cell_type}/genes.csv`).
#' The file lists one Ensembl gene id per row in a `.feature` column.
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
  inherits(x, "posteriorHCA_expr_fit") || (
    is.list(x) &&
      !is.data.frame(x) &&
      all(c("fit", "cell_type", "gene_ensg") %in% names(x))
  )
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
  if (is.list(x) && "fit" %in% names(x) && inherits(x$fit, "brmsfit")) {
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

#' Extract expression-model metadata from a fit or draws object
#' @keywords internal
#' @noRd
expression_metadata <- function(x) {
  if (is.null(x)) {
    return(list(
      cell_type = NA_character_,
      gene_ensg = NA_character_
    ))
  }

  if (is_expression_fit(x)) {
    return(list(
      cell_type = x$cell_type,
      gene_ensg = x$gene_ensg
    ))
  }

  if (is.list(x) && !is.data.frame(x)) {
    out <- list(
      cell_type = NA_character_,
      gene_ensg = NA_character_
    )
    if (!is.null(x$cell_type)) out$cell_type <- as.character(x$cell_type[[1]])
    if (!is.null(x$gene_ensg)) {
      out$gene_ensg <- as.character(x$gene_ensg[[1]])
    } else if (!is.null(x$gene)) {
      out$gene_ensg <- as.character(x$gene[[1]])
    }
    return(out)
  }

  list(
    cell_type = NA_character_,
    gene_ensg = NA_character_
  )
}
