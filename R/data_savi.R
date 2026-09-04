#' SAVI disease-associated monocyte pseudobulk (example data)
#'
#' A small Seurat object with pseudobulk RNA counts for disease-associated
#' monocytes from the SAVI cohort ([GSE226598](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226598)).
#'
#' @format A Seurat object with 17 samples (pseudobulk libraries) and sample
#'   metadata columns including `Category`, `Treatment`, `Sample`, and
#'   `CellType`.
#' @source Derived from `GSE226598_SAVI_pseudobulk_Sample_CellType.rds`; see
#'   `vignette("cohort-expression-core", package = "posteriorHCA")` or
#'   `vignette("cohort-expression-wrappers", package = "posteriorHCA")` for
#'   preparation code.
#' @seealso [scale_to_hca_reference()], [estimate_cohort_logmu()],
#'   [welch_t_test_cohort_hca()]
#' @usage data(savi_mono)
#' @examples
#' data(savi_mono)
#' ncol(savi_mono)
"savi_mono"
