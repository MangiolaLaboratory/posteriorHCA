Cohort expression workflow (core)
================
Chen Zhan
2026-09-09

SAVI case study: *ADRB2* (`ENSG00000169252`) in disease-associated
monocytes.

posteriorHCA compares a **user-supplied mean estimate and SE** to an HCA
posterior mean and SE via `welch_test_means()`. The optional edgeR
helper `estimate_ql()` fits a covariate-aware NB QL model and extracts a
user-defined fitted expression quantity (`estimate` + `se`). Any other
valid statistical model that produces a comparable estimate and SE can
be used instead.

This vignette writes the steps explicitly. The concise one-hot wrapper
path is
`vignette("cohort-expression-wrappers", package = "posteriorHCA")`.

Workflow:

1.  Jointly TMM-scale user libraries with one HCA reference
2.  Fit the user cohort with edgeR QL (`formula`), then extract a target
    expression estimand (`contrast`) and SE
3.  If that estimand is absolute, place it on the HCA scale with
    `log(E_HCA)`
4.  Load the HCA expression model and summarise the posterior
5.  Run `welch_test_means()`

## Setup

``` r
library(posteriorHCA)
library(Seurat)
library(AnnotationDbi)
library(org.Hs.eg.db)

cell_type <- "monocytic"
gene_ensg <- "ENSG00000169252"
```

## Prepare counts with Ensembl gene ids

``` r
data(savi_mono)

user_counts <- as.matrix(Seurat::GetAssayData(savi_mono, layer = "counts"))
mapped <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = rownames(user_counts),
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first"
)
keep <- !is.na(mapped) & !duplicated(mapped)
user_counts <- user_counts[keep, , drop = FALSE]
rownames(user_counts) <- unname(mapped[keep])

sample_metadata <- data.frame(
  Category = factor(savi_mono$Category[colnames(user_counts)]),
  Experiment = factor(savi_mono$Experiment[colnames(user_counts)]),
  row.names = colnames(user_counts),
  stringsAsFactors = FALSE
)

dim(user_counts)
#> [1] 17116    17
table(sample_metadata$Category, sample_metadata$Experiment)
#>               
#>                EXP1 EXP2 EXP3
#>   CTRL            3    2    2
#>   SAVI            2    2    1
#>   SAVI_treated    3    2    0
```

## Scale to the HCA reference

Joint TMM / TMMwsp puts user and HCA effective library sizes on one
scale. edgeR is then fitted on **user libraries only** with offset
`log(effective library size)`.

``` r
reference <- load_reference_sample(cell_type)

combined_counts <- merge_with_reference_sample(
  user_counts,
  reference = reference$counts,
  reference_name = reference$sample_id
)
scaling <- calculate_tmm_scaling(
  combined_counts,
  reference_name = reference$sample_id
)

user_offset <- scaling$log_effective_library_size[colnames(user_counts)]
```

## Optional edgeR helper: estimate after accounting for covariates

`estimate_ql()` is an **estimation** helper, not differential-expression
testing. The `formula` controls which effects enter the NB QL model (for
example Category and Experiment). The `contrast` selects which fitted
expression quantity to extract:

```
estimate = c' beta
SE       = sqrt(c' Var(beta) c)
```

posteriorHCA does not interpret the contrast and does not invent a
unique “batch-adjusted mean”; the user owns that estimand. Return value:

```
gene, contrast, estimate, se, df
```

`estimate` is natural-log edgeR scale (`log(mu) = X beta + offset`), not
DE-table log2 `logFC`, and not automatically named `log_mu`.

`estimate_ql()` follows the current edgeR v4 QL workflow, using the NB
dispersion estimated internally by `glmQLFit()` together with
gene-specific moderated quasi-dispersions. Earlier edgeR QL workflows
commonly used trended NB dispersions estimated with `estimateDisp()`; we
currently treat that as a sensitivity-analysis alternative. \### Simple
group mean

Under `~ 0 + Category`, `CategorySAVI` is an absolute normalized group
mean *under that parameterisation*.

``` r
savi_est <- estimate_ql(
  counts = user_counts,
  offset = user_offset,
  metadata = sample_metadata,
  formula = ~ 0 + Category,
  contrast = "CategorySAVI"
)
savi_est <- savi_est[savi_est$gene == gene_ensg, , drop = FALSE]
savi_est
#>                 gene     contrast  estimate        se       df
#> 3946 ENSG00000169252 CategorySAVI -8.208409 0.3762088 19.06436
```

Because this contrast is an absolute mean, add the fixed HCA constant
for comparison on the matched reference scale. SE is unchanged.

``` r
user_log_mu <- savi_est$estimate + scaling$hca_log_effective_library_size
user_se <- savi_est$se
n_savi <- sum(sample_metadata$Category == "SAVI")
c(log_mu = user_log_mu, se = user_se, n = n_savi)
#>    log_mu        se         n 
#> 6.9260840 0.3762088 5.0000000
```

### Accounting for Experiment

Fit Category while including Experiment in the model. Different
contrasts select different fitted conditions of the **same** model —
neither is the unique “correct” adjusted mean.

``` r
exp2 <- paste0("Experiment", levels(sample_metadata$Experiment)[[2]])

# Absolute estimands: SAVI at reference Experiment vs SAVI at EXP2
batch_fit <- estimate_ql(
  counts = user_counts,
  offset = user_offset,
  metadata = sample_metadata,
  formula = ~ 0 + Category + Experiment,
  contrast = c(
    "CategorySAVI",
    paste("CategorySAVI +", exp2)
  )
)
batch_fit <- batch_fit[batch_fit$gene == gene_ensg, , drop = FALSE]
batch_fit
#>                  gene                      contrast  estimate        se
#> 3946  ENSG00000169252                  CategorySAVI -9.962063 0.4974378
#> 21062 ENSG00000169252 CategorySAVI + ExperimentEXP2 -7.902576 0.3239370
#>             df
#> 3946  12.27496
#> 21062 12.27496
```

These are absolute fitted quantities, so add the fixed HCA constant for
comparison (SE unchanged), the same as in the simple group-mean case.

``` r
batch_fit$log_mu <- batch_fit$estimate + scaling$hca_log_effective_library_size
batch_fit[, c("gene", "contrast", "estimate", "log_mu", "se")]
#>                  gene                      contrast  estimate   log_mu
#> 3946  ENSG00000169252                  CategorySAVI -9.962063 5.172430
#> 21062 ENSG00000169252 CategorySAVI + ExperimentEXP2 -7.902576 7.231917
#>              se
#> 3946  0.4974378
#> 21062 0.3239370
```

A relative contrast does **not** get the HCA shift: if both group means
were shifted by the same `log(E_HCA)`, that constant would cancel in the
difference. Use `estimate` / `se` as returned.

``` r
hca_const <- scaling$hca_log_effective_library_size

rel_fit <- estimate_ql(
  counts = user_counts,
  offset = user_offset,
  metadata = sample_metadata,
  formula = ~ 0 + Category + Experiment,
  contrast = c(
    "CategorySAVI",
    "CategoryCTRL",
    "CategorySAVI - CategoryCTRL"
  )
)
rel_fit <- rel_fit[rel_fit$gene == gene_ensg, , drop = FALSE]
rel_fit
#>                  gene                    contrast   estimate        se       df
#> 3946  ENSG00000169252                CategorySAVI  -9.962063 0.4974378 12.27496
#> 21062 ENSG00000169252                CategoryCTRL -13.301385 0.6515197 12.27496
#> 38178 ENSG00000169252 CategorySAVI - CategoryCTRL   3.339321 0.5183784 12.27496

est <- setNames(rel_fit$estimate, rel_fit$contrast)
# Absolute means may be HCA-shifted; their difference is unchanged.
c(
  relative_contrast = unname(est[["CategorySAVI - CategoryCTRL"]]),
  log_mu_CategorySAVI = est[["CategorySAVI"]] + hca_const,
  log_mu_CategoryCTRL = est[["CategoryCTRL"]] + hca_const,
  difference_of_hca_shifted_absolutes =
    (est[["CategorySAVI"]] + hca_const) - (est[["CategoryCTRL"]] + hca_const)
)
#>                   relative_contrast                 log_mu_CategorySAVI 
#>                            3.339321                            5.172430 
#>                 log_mu_CategoryCTRL difference_of_hca_shifted_absolutes 
#>                            1.833109                            3.339321
```

## HCA posterior

``` r
expression_fit <- load_expression_fit(
  cell_type = cell_type,
  gene_ensg = gene_ensg
)

newdata <- build_newdata_grid(
  expression_fit,
  disease_groups = "Normal",
  tissue_groups = "blood",
  assay_groups = "10x Genomics 3"
)

posterior_draws <- expression_draws(
  expression_fit,
  newdata = newdata,
  quantity = "linpred",
  collapse = "mean"
)
```

## Comparison

`welch_test_means()` only needs two means and two SEs on the same scale.
Here we compare the simple SAVI group mean (HCA-shifted) to the HCA
baseline. Estimates from outside posteriorHCA can be passed the same
way.

``` r
posterior_summary <- summarize_posterior_draws(
  posterior_draws,
  value = user_log_mu
)

welch_test_means(
  mu1 = user_log_mu,
  se1 = user_se,
  mu2 = posterior_summary$log_mu,
  se2 = posterior_summary$se,
  n1 = n_savi,
  n2 = posterior_summary$n
)
#> $mu1
#> [1] 6.926084
#> 
#> $se1
#> [1] 0.3762088
#> 
#> $n1
#> [1] 5
#> 
#> $mu2
#> [1] 4.44304
#> 
#> $se2
#> [1] 0.865169
#> 
#> $n2
#> [1] 400
#> 
#> $delta
#> [1] 2.483044
#> 
#> $se_diff
#> [1] 0.9434249
#> 
#> $t_stat
#> [1] 2.631947
#> 
#> $df
#> [1] 123.5459
#> 
#> $p_value
#> [1] 0.009572266
```

``` r
test <- welch_test_means(
  mu1 = user_log_mu,
  se1 = user_se,
  mu2 = posterior_summary$log_mu,
  se2 = posterior_summary$se,
  n1 = n_savi,
  n2 = posterior_summary$n
)
test_results <- data.frame(
  group = "SAVI",
  n = n_savi,
  log_mu = test$mu1,
  se = test$se1,
  hca_log_mu = test$mu2,
  hca_se = test$se2,
  p_value = test$p_value,
  stringsAsFactors = FALSE
)
test_results
#>   group n   log_mu        se hca_log_mu   hca_se     p_value
#> 1  SAVI 5 6.926084 0.3762088    4.44304 0.865169 0.009572266
```

## Plots

``` r
plot_hca_draws(
  draws = posterior_draws,
  subtitle = "Normal, 10x Genomics 3 healthy baseline"
)
```

![](cohort-expression-core_files/figure-gfm/plot-hca-1.png)<!-- -->

``` r
plot_cohort_vs_hca(
  hca_draws = posterior_draws,
  cohort_est = test_results,
  subtitle = "SAVI Category mean via edgeR formula + contrast",
  annotate = c("group", "p_value")
)
```

![](cohort-expression-core_files/figure-gfm/plot-cohort-1.png)<!-- -->

## Session info

``` r
sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: aarch64-apple-darwin23
#> Running under: macOS Tahoe 26.6.2
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: Australia/Adelaide
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] org.Hs.eg.db_3.23.1  AnnotationDbi_1.74.0 IRanges_2.46.0      
#>  [4] S4Vectors_0.50.2     Biobase_2.72.0       BiocGenerics_0.58.1 
#>  [7] generics_0.1.4       Seurat_5.5.1         SeuratObject_5.4.0  
#> [10] sp_2.2-3             posteriorHCA_0.2.0   testthat_3.3.2      
#> 
#> loaded via a namespace (and not attached):
#>   [1] RcppAnnoy_0.0.23            splines_4.6.1              
#>   [3] later_1.4.8                 tibble_3.3.1               
#>   [5] polyclip_1.10-7             brms_2.23.0                
#>   [7] fastDummies_1.7.6           lifecycle_1.0.5            
#>   [9] StanHeaders_2.39.1          edgeR_4.10.4               
#>  [11] rprojroot_2.1.1             vroom_1.7.1                
#>  [13] processx_3.9.0              globals_0.19.1             
#>  [15] sccomp_2.4.0                lattice_0.23-1             
#>  [17] MASS_7.3-66                 backports_1.5.1            
#>  [19] magrittr_2.0.5              limma_3.68.5               
#>  [21] plotly_4.12.1               rmarkdown_2.32             
#>  [23] yaml_2.3.12                 httpuv_1.6.17              
#>  [25] otel_0.2.0                  sctransform_0.4.3          
#>  [27] spam_2.11-4                 sessioninfo_1.2.4          
#>  [29] pkgbuild_1.4.8              spatstat.sparse_3.2-0      
#>  [31] reticulate_1.47.0           cowplot_1.2.0              
#>  [33] pbapply_1.7-5               DBI_1.3.0                  
#>  [35] RColorBrewer_1.1-3          abind_1.4-8                
#>  [37] pkgload_1.5.3               GenomicRanges_1.64.0       
#>  [39] Rtsne_0.17                  purrr_1.2.2                
#>  [41] tensorA_0.36.2.1            inline_0.3.21              
#>  [43] ggrepel_0.9.8               irlba_2.3.7                
#>  [45] listenv_1.0.0               spatstat.utils_3.2-4       
#>  [47] goftest_1.2-3               RSpectra_0.16-2            
#>  [49] spatstat.random_3.5-1       bridgesampling_1.2-1       
#>  [51] fitdistrplus_1.2-6          parallelly_1.48.0          
#>  [53] DelayedArray_0.38.2         codetools_0.2-20           
#>  [55] tidyselect_1.2.1            bayesplot_1.16.0           
#>  [57] farver_2.1.2                matrixStats_1.5.0          
#>  [59] spatstat.explore_3.8-2      Seqinfo_1.2.0              
#>  [61] jsonlite_2.0.0              ellipsis_0.3.3             
#>  [63] progressr_1.0.0             ggridges_0.5.7             
#>  [65] survival_3.8-11             tools_4.6.1                
#>  [67] ica_1.0-3                   Rcpp_1.1.2                 
#>  [69] glue_1.8.1                  SparseArray_1.12.2         
#>  [71] gridExtra_2.3.1             qs2_0.3.1                  
#>  [73] xfun_0.60                   cmdstanr_0.9.0             
#>  [75] MatrixGenerics_1.24.0       distributional_0.8.1       
#>  [77] usethis_3.2.1               dplyr_1.2.1                
#>  [79] withr_3.0.3                 loo_2.10.1                 
#>  [81] instantiate_0.2.3           fastmap_1.2.0              
#>  [83] callr_3.8.0                 digest_0.6.39              
#>  [85] R6_2.6.1                    mime_0.13                  
#>  [87] scattermore_1.2             tensor_1.5.1               
#>  [89] spatstat.data_3.1-9         RSQLite_3.53.3             
#>  [91] tidyr_1.3.2                 data.table_1.18.6.1        
#>  [93] S4Arrays_1.12.0             httr_1.4.9                 
#>  [95] htmlwidgets_1.6.4           uwot_0.2.5                 
#>  [97] pkgconfig_2.0.3             gtable_0.3.6               
#>  [99] blob_1.3.0                  lmtest_0.9-40              
#> [101] S7_0.2.2                    SingleCellExperiment_1.34.0
#> [103] XVector_0.52.0              brio_1.1.5                 
#> [105] htmltools_0.5.9             dotCall64_1.2              
#> [107] scales_1.4.0                png_0.1-9                  
#> [109] posterior_1.7.0             spatstat.univar_3.2-0      
#> [111] knitr_1.52                  rstudioapi_0.19.0          
#> [113] tzdb_0.5.0                  reshape2_1.4.5             
#> [115] coda_0.19-4.1               checkmate_2.3.4            
#> [117] nlme_3.1-171                cachem_1.1.0               
#> [119] zoo_1.9-0                   stringr_1.6.0              
#> [121] KernSmooth_2.23-27          parallel_4.6.1             
#> [123] miniUI_0.1.2                desc_1.4.3                 
#> [125] pillar_1.11.1               grid_4.6.1                 
#> [127] vctrs_0.7.3                 RANN_2.6.3                 
#> [129] promises_1.5.0              stringfish_0.19.2          
#> [131] xtable_1.8-8                cluster_2.1.8.3            
#> [133] evaluate_1.0.5              readr_2.2.0                
#> [135] locfit_1.5-9.12             mvtnorm_1.4-2              
#> [137] cli_3.6.6                   compiler_4.6.1             
#> [139] rlang_1.3.0                 crayon_1.5.3               
#> [141] rstantools_2.7.1            future.apply_1.20.2        
#> [143] labeling_0.4.3              forcats_1.0.1              
#> [145] plyr_1.8.9                  fs_2.1.0                   
#> [147] rstan_2.32.7                stringi_1.8.9              
#> [149] QuickJSR_1.11.0             viridisLite_0.4.3          
#> [151] deldir_2.0-4                Biostrings_2.80.2          
#> [153] devtools_2.5.2              spatstat.geom_3.8-2        
#> [155] Brobdingnag_1.2-9           Matrix_1.7-6               
#> [157] RcppHNSW_0.7.0              hms_1.1.4                  
#> [159] patchwork_1.3.2             bit64_4.8.6                
#> [161] future_1.75.0               ggplot2_4.0.3              
#> [163] statmod_1.5.2               KEGGREST_1.52.2            
#> [165] shiny_1.14.0                SummarizedExperiment_1.42.0
#> [167] ROCR_1.0-12                 igraph_2.3.3               
#> [169] memoise_2.0.1               RcppParallel_6.2.1         
#> [171] bit_4.6.0
```
