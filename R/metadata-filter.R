#' Resolve one metadata argument to model levels
#'
#' Filter calls such as `str_detect("10x")` get the covariate inserted as the
#' first argument when omitted, then run through `dplyr::filter()`. Exact
#' values (`"A"`, `c("A", "B")`) are matched literally. Blank/`NA` returns
#' `NA_expand = TRUE` so callers can expand (expression) or pass `NA`
#' (composition).
#' @keywords internal
#' @noRd
pick_metadata_levels <- function(quo, levels, name) {
  expr <- rlang::quo_get_expr(quo)
  if (
    rlang::quo_is_missing(quo) ||
      rlang::quo_is_null(quo) ||
      (is.atomic(expr) && length(expr) == 1L && is.na(expr))
  ) {
    return(list(chosen = character(), NA_expand = TRUE))
  }

  # dplyr filter path: patch column name, then filter()
  if (is.call(expr) && !identical(expr[[1]], quote(c))) {
    if (!name %in% all.names(expr)) {
      quo <- rlang::new_quosure(
        as.call(c(list(expr[[1]]), as.symbol(name), as.list(expr)[-1])),
        rlang::quo_get_env(quo)
      )
    }
    chr <- as.character(levels)
    num <- suppressWarnings(as.numeric(chr))
    data <- stats::setNames(
      data.frame(if (length(chr) && !anyNA(num)) num else chr),
      name
    )
    return(list(
      chosen = as.character(dplyr::filter(data, !!quo)[[name]]),
      NA_expand = FALSE
    ))
  }

  # Exact values
  value <- rlang::eval_tidy(quo)
  if (length(value) == 1L && is.na(value)) {
    return(list(chosen = character(), NA_expand = TRUE))
  }
  chosen <- unique(as.character(value))
  chosen <- chosen[!is.na(chosen) & nzchar(chosen)]
  if (!length(chosen)) {
    return(list(chosen = character(), NA_expand = TRUE))
  }
  unknown <- setdiff(chosen, levels)
  if (length(levels) && length(unknown)) {
    cli::cli_abort(c(
      "Unknown {name} value{?s}: {unknown}.",
      "i" = "This model includes: {levels}."
    ))
  }
  list(chosen = chosen, NA_expand = FALSE)
}
