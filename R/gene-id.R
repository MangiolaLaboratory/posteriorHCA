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

#' Resolve gene symbols or Ensembl ids to canonical ENSG ids
#'
#' Symbols are mapped with `org.Hs.eg.db` when available. When `cell_type` or
#' `universe` is supplied, resolved ids must appear in that universe.
#'
#' @param genes Character vector of symbols and/or Ensembl ids.
#' @param cell_type Optional cell type for universe validation.
#' @param version Nectar version pin passed to [get_gene_universe()].
#' @param universe Optional preloaded ENSG universe. Overrides `cell_type`.
#' @param orgdb Optional `OrgDb` object for symbol mapping. Defaults to
#'   `org.Hs.eg.db` when installed.
#' @param id_type `"auto"` guesses from the input pattern; `"ensembl"` treats
#'   inputs as Ensembl ids; `"symbol"` maps symbols.
#' @param strict If `TRUE`, error when a gene cannot be resolved or is absent
#'   from the universe.
#' @param cache_directory,use_cache Passed to [get_gene_universe()] when needed.
#' @return A named character vector: names are the input `genes`, values are
#'   canonical ENSG ids (or `NA` when unresolved and `strict = FALSE`).
#' @export
#' @importFrom cli cli_abort
resolve_gene <- function(
  genes,
  cell_type = NULL,
  version = "latest",
  universe = NULL,
  orgdb = NULL,
  id_type = c("auto", "ensembl", "symbol"),
  strict = TRUE,
  cache_directory = get_default_cache_dir(),
  use_cache = TRUE
) {
  id_type <- match.arg(id_type)
  genes <- as.character(genes)
  if (!length(genes)) {
    return(setNames(character(0), character(0)))
  }

  input_names <- if (is.null(names(genes)) || !nzchar(names(genes)[[1]])) {
    genes
  } else {
    names(genes)
  }
  genes <- strip_ensembl_version(genes)

  if (id_type == "auto") {
    id_type <- if (mean(is_ensembl_gene_id(genes), na.rm = TRUE) >= 0.5) {
      "ensembl"
    } else {
      "symbol"
    }
  }

  resolved <- rep(NA_character_, length(genes))
  names(resolved) <- input_names

  if (id_type == "ensembl") {
    is_ensg <- is_ensembl_gene_id(genes)
    resolved[is_ensg] <- toupper(genes[is_ensg])
    if (strict && any(!is_ensg)) {
      cli_abort("Not valid Ensembl gene ids: {genes[!is_ensg]}.")
    }
  } else {
    if (is.null(orgdb)) {
      if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
          !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        cli_abort(c(
          "Gene symbol mapping requires AnnotationDbi and org.Hs.eg.db.",
          "i" = "Install them or pass Ensembl ids instead."
        ))
      }
      orgdb <- get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db"))
    }
    mapped <- suppressMessages(AnnotationDbi::mapIds(
      orgdb,
      keys = genes,
      column = "ENSEMBL",
      keytype = "SYMBOL",
      multiVals = "first"
    ))
    resolved <- strip_ensembl_version(unname(mapped))
    names(resolved) <- input_names
    if (strict && any(is.na(resolved))) {
      cli_abort("Could not map gene symbol{?s} to ENSG: {genes[is.na(resolved)]}.")
    }
  }

  if (!is.null(universe) || !is.null(cell_type)) {
    if (is.null(universe)) {
      universe <- get_gene_universe(
        cell_type,
        version = version,
        cache_directory = cache_directory,
        use_cache = use_cache
      )
    }
    missing <- setdiff(resolved[!is.na(resolved)], universe)
    if (length(missing) > 0L) {
      if (strict) {
        label <- if (!is.null(cell_type)) cell_type else "supplied universe"
        cli_abort("Gene{?s} not in `{label}` universe: {missing}.")
      }
      resolved[resolved %in% missing] <- NA_character_
    }
  }

  resolved
}

#' Resolve one gene to a single ENSG id
#' @inheritParams resolve_gene
#' @return Character scalar ENSG id.
#' @export
resolve_gene_one <- function(
  gene,
  cell_type = NULL,
  version = "latest",
  universe = NULL,
  orgdb = NULL,
  id_type = c("auto", "ensembl", "symbol"),
  strict = TRUE,
  cache_directory = get_default_cache_dir(),
  use_cache = TRUE
) {
  if (length(gene) != 1L || is.na(gene) || !nzchar(gene)) {
    cli_abort("`gene` must be a single non-empty gene identifier.")
  }
  out <- resolve_gene(
    gene,
    cell_type = cell_type,
    version = version,
    universe = universe,
    orgdb = orgdb,
    id_type = id_type,
    strict = strict,
    cache_directory = cache_directory,
    use_cache = use_cache
  )
  unname(out[[1]])
}

#' Harmonise gene identifiers to canonical ENSG ids
#'
#' Maps gene symbols to Ensembl ids and drops unmapped or duplicated genes.
#' The output container matches the input class: `Seurat` in gives `Seurat`
#' out, `SummarizedExperiment` in gives `SummarizedExperiment` out, and a
#' matrix or aligned list is returned in the same form.
#'
#' @param x A matrix, `Matrix`, `SummarizedExperiment`, `Seurat` object, or
#'   aligned list from [scale_to_hca_reference()].
#' @param assay Assay name for SE / Seurat input. Default `"counts"` when present.
#' @inheritParams resolve_gene
#' @return An object of the same class as `x` with harmonised ENSG feature names.
#' @export
harmonise_gene_ids <- function(
  x,
  cell_type = NULL,
  version = "latest",
  universe = NULL,
  orgdb = NULL,
  id_type = c("auto", "ensembl", "symbol"),
  strict = FALSE,
  assay = NULL,
  cache_directory = get_default_cache_dir(),
  use_cache = TRUE
) {
  mat <- extract_count_matrix(x, assay = assay, arg_name = "x")

  if (is.null(rownames(mat))) {
    cli_abort("`x` must have gene rownames.")
  }

  resolved <- resolve_gene(
    rownames(mat),
    cell_type = cell_type,
    version = version,
    universe = universe,
    orgdb = orgdb,
    id_type = id_type,
    strict = strict,
    cache_directory = cache_directory,
    use_cache = use_cache
  )
  keep <- !is.na(resolved) & !duplicated(resolved)
  if (!any(keep)) {
    cli_abort("No genes could be mapped to unique ENSG ids.")
  }

  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- resolved[keep]
  storage.mode(mat) <- "double"

  rebuild_count_container(x, mat, assay = assay)
}

#' @keywords internal
#' @noRd
is_expr_fit <- function(x) {
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
  if (is_expr_fit(x)) {
    return(x$fit)
  }
  if (is.list(x) && "fit" %in% names(x) && inherits(x$fit, "brmsfit")) {
    return(x$fit)
  }
  x
}

#' @keywords internal
#' @noRd
new_expr_fit <- function(fit, cell_type, gene_ensg, gene_symbol = NA_character_) {
  structure(
    list(
      fit = fit,
      cell_type = as.character(cell_type),
      gene_ensg = as.character(gene_ensg),
      gene_symbol = as.character(gene_symbol)
    ),
    class = c("posteriorHCA_expr_fit", "list")
  )
}

#' Extract expression-model metadata from a fit or draws object
#' @keywords internal
#' @noRd
expr_metadata <- function(x) {
  if (is.null(x)) {
    return(list(
      cell_type = NA_character_,
      gene_ensg = NA_character_,
      gene_symbol = NA_character_
    ))
  }

  if (is_expr_fit(x)) {
    return(list(
      cell_type = x$cell_type,
      gene_ensg = x$gene_ensg,
      gene_symbol = if (!is.null(x$gene_symbol)) x$gene_symbol else NA_character_
    ))
  }

  if (is.list(x) && !is.data.frame(x)) {
    out <- list(
      cell_type = NA_character_,
      gene_ensg = NA_character_,
      gene_symbol = NA_character_
    )
    if (!is.null(x$cell_type)) out$cell_type <- as.character(x$cell_type[[1]])
    if (!is.null(x$gene_ensg)) {
      out$gene_ensg <- as.character(x$gene_ensg[[1]])
    } else if (!is.null(x$gene)) {
      out$gene_ensg <- as.character(x$gene[[1]])
    }
    if (!is.null(x$gene_symbol)) out$gene_symbol <- as.character(x$gene_symbol[[1]])
    return(out)
  }

  list(
    cell_type = NA_character_,
    gene_ensg = NA_character_,
    gene_symbol = NA_character_
  )
}

#' Ensure two metadata objects refer to the same gene and cell type
#' @keywords internal
#' @noRd
validate_expr_metadata_match <- function(left, right, context = "inputs") {
  for (field in c("gene_ensg", "cell_type")) {
    lval <- left[[field]]
    rval <- right[[field]]
    if (!is.na(lval) && nzchar(lval) && !is.na(rval) && nzchar(rval) && !identical(lval, rval)) {
      cli_abort(c(
        "Mismatched {field} between {context}.",
        "i" = "Left: `{lval}`; right: `{rval}`."
      ))
    }
  }
  invisible(TRUE)
}

#' Extract cohort-estimate metadata
#' @keywords internal
#' @noRd
cohort_est_metadata <- function(cohort_est) {
  if (is.data.frame(cohort_est)) {
    gene_ensg <- if ("gene_ensg" %in% names(cohort_est)) {
      as.character(cohort_est$gene_ensg[[1]])
    } else if ("gene" %in% names(cohort_est)) {
      as.character(cohort_est$gene[[1]])
    } else {
      NA_character_
    }
    return(list(
      cell_type = if ("cell_type" %in% names(cohort_est)) {
        as.character(cohort_est$cell_type[[1]])
      } else {
        NA_character_
      },
      gene_ensg = gene_ensg,
      gene_symbol = if ("gene_symbol" %in% names(cohort_est)) {
        as.character(cohort_est$gene_symbol[[1]])
      } else {
        NA_character_
      }
    ))
  }

  if (is.list(cohort_est) && all(c("log_mu", "se") %in% names(cohort_est))) {
    return(expr_metadata(cohort_est))
  }

  list(
    cell_type = NA_character_,
    gene_ensg = NA_character_,
    gene_symbol = NA_character_
  )
}
