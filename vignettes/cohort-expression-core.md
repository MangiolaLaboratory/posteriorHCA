Cohort expression workflow (core)
================
Chen Zhan
2026-09-08

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

``` text
estimate = c' beta
SE       = sqrt(c' Var(beta) c)
```

posteriorHCA does not interpret the contrast and does not invent a
unique “batch-adjusted mean”; the user owns that estimand. Return value:

``` text
gene, contrast, estimate, se, df, dispersion
```

`estimate` is natural-log edgeR scale (`log(mu) = X beta + offset`), not
DE-table log2 `logFC`, and not automatically named `log_mu`.

### Simple group mean

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
#>                 gene     contrast  estimate        se      df dispersion
#> 3946 ENSG00000169252 CategorySAVI -8.208456 0.3763713 19.0649  0.2340799
```

Because this contrast is an absolute mean, add the fixed HCA constant
for comparison on the matched reference scale. SE is unchanged.

``` r
user_log_mu <- savi_est$estimate + scaling$hca_log_effective_library_size
user_se <- savi_est$se
n_savi <- sum(sample_metadata$Category == "SAVI")
c(log_mu = user_log_mu, se = user_se, n = n_savi)
#>    log_mu        se         n 
#> 6.9260372 0.3763713 5.0000000
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
#>                  gene                      contrast  estimate        se     df
#> 3946  ENSG00000169252                  CategorySAVI -9.962119 0.4985965 12.275
#> 21062 ENSG00000169252 CategorySAVI + ExperimentEXP2 -7.902632 0.3247185 12.275
#>       dispersion
#> 3946   0.1740796
#> 21062  0.1740796
```

These are absolute fitted quantities, so add the fixed HCA constant for
comparison (SE unchanged), the same as in the simple group-mean case.

``` r
batch_fit$log_mu <- batch_fit$estimate + scaling$hca_log_effective_library_size
batch_fit[, c("gene", "contrast", "estimate", "log_mu", "se")]
#>                  gene                      contrast  estimate   log_mu
#> 3946  ENSG00000169252                  CategorySAVI -9.962119 5.172374
#> 21062 ENSG00000169252 CategorySAVI + ExperimentEXP2 -7.902632 7.231861
#>              se
#> 3946  0.4985965
#> 21062 0.3247185
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
#>                  gene                    contrast   estimate        se     df
#> 3946  ENSG00000169252                CategorySAVI  -9.962119 0.4985965 12.275
#> 21062 ENSG00000169252                CategoryCTRL -13.301634 0.6530466 12.275
#> 38178 ENSG00000169252 CategorySAVI - CategoryCTRL   3.339515 0.5196051 12.275
#>       dispersion
#> 3946   0.1740796
#> 21062  0.1740796
#> 38178  0.1740796

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
#>                            3.339515                            5.172374 
#>                 log_mu_CategoryCTRL difference_of_hca_shifted_absolutes 
#>                            1.832860                            3.339515
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
#> [1] 6.926037
#> 
#> $se1
#> [1] 0.3763713
#> 
#> $n1
#> [1] 5
#> 
#> $mu2
#> [1] 4.435144
#> 
#> $se2
#> [1] 0.8693214
#> 
#> $n2
#> [1] 400
#> 
#> $delta
#> [1] 2.490894
#> 
#> $se_diff
#> [1] 0.9472988
#> 
#> $t_stat
#> [1] 2.62947
#> 
#> $df
#> [1] 124.8902
#> 
#> $p_value
#> [1] 0.009626487
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
#>   group n   log_mu        se hca_log_mu    hca_se     p_value
#> 1  SAVI 5 6.926037 0.3763713   4.435144 0.8693214 0.009626487
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
#> R version 4.6.1 (2026-06-24 ucrt)
#> Platform: x86_64-w64-mingw32/x64
#> Running under: Windows 11 x64 (build 26200)
#> 
#> Matrix products: default
#>   LAPACK version 3.12.1
#> 
#> locale:
#> [1] LC_COLLATE=English_United Kingdom.utf8 
#> [2] LC_CTYPE=English_United Kingdom.utf8   
#> [3] LC_MONETARY=English_United Kingdom.utf8
#> [4] LC_NUMERIC=C                           
#> [5] LC_TIME=English_United Kingdom.utf8    
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
#>  [4] S4Vectors_0.50.1     Biobase_2.72.0       BiocGenerics_0.58.1 
#>  [7] generics_0.1.4       Seurat_5.5.1         SeuratObject_5.4.0  
#> [10] sp_2.2-3             posteriorHCA_0.2.0   testthat_3.3.2      
#> 
#> loaded via a namespace (and not attached):
#>   [1] RcppAnnoy_0.0.23            splines_4.6.1              
#>   [3] later_1.4.8                 tibble_3.3.1               
#>   [5] polyclip_1.10-7             brms_2.23.0                
#>   [7] fastDummies_1.7.6           lifecycle_1.0.5            
#>   [9] StanHeaders_2.39.1          edgeR_4.10.1               
#>  [11] rprojroot_2.1.1             vroom_1.7.1                
#>  [13] processx_3.9.0              sccomp_2.4.0               
#>  [15] globals_0.19.1              lattice_0.22-9             
#>  [17] MASS_7.3-65                 backports_1.5.1            
#>  [19] magrittr_2.0.5              limma_3.68.4               
#>  [21] plotly_4.12.1               rmarkdown_2.32             
#>  [23] yaml_2.3.12                 httpuv_1.6.17              
#>  [25] otel_0.2.0                  sctransform_0.4.3          
#>  [27] spam_2.11-4                 sessioninfo_1.2.4          
#>  [29] pkgbuild_1.4.8              spatstat.sparse_3.2-0      
#>  [31] reticulate_1.47.0           cowplot_1.2.0              
#>  [33] pbapply_1.7-5               DBI_1.3.0                  
#>  [35] RColorBrewer_1.1-3          multcomp_1.4-32            
#>  [37] abind_1.4-8                 pkgload_1.5.3              
#>  [39] GenomicRanges_1.64.0        Rtsne_0.17                 
#>  [41] purrr_1.2.2                 TH.data_1.1-5              
#>  [43] tensorA_0.36.2.1            sandwich_3.1-3             
#>  [45] inline_0.3.21               ggrepel_0.9.8              
#>  [47] irlba_2.3.7                 listenv_1.0.0              
#>  [49] spatstat.utils_3.2-4        goftest_1.2-3              
#>  [51] RSpectra_0.16-2             spatstat.random_3.5-1      
#>  [53] bridgesampling_1.2-1        fitdistrplus_1.2-6         
#>  [55] parallelly_1.48.0           DelayedArray_0.38.2        
#>  [57] codetools_0.2-20            tidyselect_1.2.1           
#>  [59] bayesplot_1.16.0            farver_2.1.2               
#>  [61] matrixStats_1.5.0           spatstat.explore_3.8-2     
#>  [63] Seqinfo_1.2.0               jsonlite_2.0.0             
#>  [65] ellipsis_0.3.3              progressr_1.0.0            
#>  [67] ggridges_0.5.7              survival_3.8-6             
#>  [69] emmeans_2.0.4               tools_4.6.1                
#>  [71] ica_1.0-3                   Rcpp_1.1.2                 
#>  [73] glue_1.8.1                  SparseArray_1.12.2         
#>  [75] gridExtra_2.3.1             qs2_0.3.1                  
#>  [77] xfun_0.60                   cmdstanr_0.9.0             
#>  [79] MatrixGenerics_1.24.0       distributional_0.8.1       
#>  [81] usethis_3.2.1               dplyr_1.2.1                
#>  [83] withr_3.0.3                 loo_2.10.1                 
#>  [85] instantiate_0.2.3           fastmap_1.2.0              
#>  [87] callr_3.8.0                 digest_0.6.39              
#>  [89] R6_2.6.1                    mime_0.13                  
#>  [91] estimability_2.0.0          scattermore_1.2            
#>  [93] tensor_1.5.1                spatstat.data_3.1-9        
#>  [95] RSQLite_3.53.3              tidyr_1.3.2                
#>  [97] data.table_1.18.6.1         S4Arrays_1.12.0            
#>  [99] httr_1.4.9                  htmlwidgets_1.6.4          
#> [101] uwot_0.2.5                  pkgconfig_2.0.3            
#> [103] gtable_0.3.6                blob_1.3.0                 
#> [105] lmtest_0.9-40               S7_0.2.2                   
#> [107] SingleCellExperiment_1.34.0 XVector_0.52.0             
#> [109] brio_1.1.5                  htmltools_0.5.9            
#> [111] dotCall64_1.2               scales_1.4.0               
#> [113] png_0.1-9                   posterior_1.7.0            
#> [115] spatstat.univar_3.2-0       knitr_1.51                 
#> [117] rstudioapi_0.19.0           tzdb_0.5.0                 
#> [119] reshape2_1.4.5              coda_0.19-4.1              
#> [121] checkmate_2.3.4             nlme_3.1-169               
#> [123] cachem_1.1.0                zoo_1.9-0                  
#> [125] stringr_1.6.0               KernSmooth_2.23-26         
#> [127] parallel_4.6.1              miniUI_0.1.2               
#> [129] desc_1.4.3                  pillar_1.11.1              
#> [131] grid_4.6.1                  vctrs_0.7.3                
#> [133] RANN_2.6.3                  promises_1.5.0             
#> [135] stringfish_0.19.2           xtable_1.8-8               
#> [137] cluster_2.1.8.2             evaluate_1.0.5             
#> [139] readr_2.2.0                 locfit_1.5-9.12            
#> [141] mvtnorm_1.4-2               cli_3.6.6                  
#> [143] compiler_4.6.1              rlang_1.3.0                
#> [145] crayon_1.5.3                rstantools_2.7.1           
#> [147] future.apply_1.20.2         labeling_0.4.3             
#> [149] forcats_1.0.1               plyr_1.8.9                 
#> [151] fs_2.1.0                    rstan_2.32.7               
#> [153] stringi_1.8.9               QuickJSR_1.11.0            
#> [155] viridisLite_0.4.3           deldir_2.0-4               
#> [157] Biostrings_2.80.1           devtools_2.5.2             
#> [159] spatstat.geom_3.8-2         Brobdingnag_1.2-9          
#> [161] Matrix_1.7-5                RcppHNSW_0.7.0             
#> [163] hms_1.1.4                   patchwork_1.3.2            
#> [165] bit64_4.8.6                 future_1.75.0              
#> [167] ggplot2_4.0.3               statmod_1.5.2              
#> [169] KEGGREST_1.52.2             shiny_1.14.0               
#> [171] SummarizedExperiment_1.42.0 ROCR_1.0-12                
#> [173] igraph_2.3.3                memoise_2.0.1              
#> [175] RcppParallel_6.2.1          bit_4.6.0
```
