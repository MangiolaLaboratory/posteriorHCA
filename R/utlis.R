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
