Cohort expression workflow
================
Chen Zhan
2026-09-07

SAVI case study: *ADRB2* (`ENSG00000169252`) in disease-associated
monocytes.

Workflow:

1.  Scale user libraries to one HCA reference (TMM offsets)
2.  Build a design with `model.matrix()` (`~ 0 + Category` or
    intercept-only `~ 1`)
3.  Estimate user expression means with `estimate_logmu_ql()`
4.  Load the HCA expression model and draw posteriors
5.  Summarise draws and run a Welch test

## Setup

``` r
library(posteriorHCA)
library(Seurat)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(purrr)

cell_type <- "monocytic"
gene_ensg <- "ENSG00000169252"
```

## Prepare counts with Ensembl gene ids

Atlas models are keyed by ENSG. Map gene symbols once with AnnotationDbi
before the posteriorHCA workflow.

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
  row.names = colnames(user_counts),
  stringsAsFactors = FALSE
)

dim(user_counts)
#> [1] 17116    17
table(sample_metadata$Category)
#> 
#>         CTRL         SAVI SAVI_treated 
#>            7            5            5
```

## Scale to the HCA reference

The reference library is used only to compute TMM offsets. User
estimation below uses user libraries and their offsets only.

``` r
reference <- load_reference_sample(cell_type)

combined_counts <- merge_with_reference_sample(
  user_counts,
  reference = reference$counts,
  reference_name = reference$sample_id
)
scaling <- calculate_tmm_offset(
  combined_counts,
  reference_name = reference$sample_id
)

user_offset <- scaling$offset[colnames(user_counts)]
```

## Estimate user expression means

Cell-means design (`~ 0 + Category`):

``` r
design_matrix <- model.matrix(~ 0 + Category, data = sample_metadata)
colnames(design_matrix) <- sub("^Category", "", colnames(design_matrix))

expression_estimates <- estimate_logmu_ql(
  user_counts,
  user_offset,
  design_matrix
)

expression_estimates <- expression_estimates[
  expression_estimates$gene == gene_ensg,
  ,
  drop = FALSE
]
expression_estimates
#>                  gene        group n   log_mu         mu        se       df
#> 3946  ENSG00000169252         CTRL 7 3.301308   27.14812 0.5652071 19.04548
#> 21062 ENSG00000169252         SAVI 5 6.928445 1020.90506 0.3482651 19.04548
#> 38178 ENSG00000169252 SAVI_treated 5 6.283036  535.41145 0.3697428 19.04548
#>       dispersion
#> 3946   0.2340799
#> 21062  0.2340799
#> 38178  0.2340799
```

Intercept-only design (`~ 1`): subset to one `Category`, fit
`model.matrix(~ 1, ...)`, and treat the intercept as that cohort’s
absolute log(μ).

``` r
expression_estimates_by_level <- map_dfr(
  levels(sample_metadata$Category),
  function(category) {
    sample_ids <- rownames(sample_metadata)[sample_metadata$Category == category]
    design_one <- model.matrix(
      ~ 1,
      data = sample_metadata[sample_ids, , drop = FALSE]
    )
    out <- estimate_logmu_ql(
      user_counts[, sample_ids, drop = FALSE],
      user_offset[sample_ids],
      design_one
    )
    out <- out[out$gene == gene_ensg, , drop = FALSE]
    out$group <- category
    out
  }
)
expression_estimates_by_level
#>              gene        group n   log_mu        mu        se        df
#> 1 ENSG00000169252         CTRL 7 3.328751  27.90347 0.3911294 10.902720
#> 2 ENSG00000169252         SAVI 5 6.899132 991.41344 0.2529245  4.015891
#> 3 ENSG00000169252 SAVI_treated 5 6.288969 538.59771 0.2051681  3.957096
#>   dispersion
#> 1  0.2326938
#> 2  0.3595043
#> 3  0.2110409
```

## Load HCA model and draw the healthy baseline

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

## Summarise posterior and Welch test

``` r
cohort_estimate <- expression_estimates[
  expression_estimates$group == "SAVI",
  ,
  drop = FALSE
]

posterior_summary <- summarize_posterior_draws(
  posterior_draws,
  value = cohort_estimate$log_mu
)

welch_test_means(
  cohort_estimate$log_mu,
  cohort_estimate$se,
  posterior_summary$mean,
  posterior_summary$sd,
  n1 = cohort_estimate$n,
  n2 = posterior_summary$n
)
#> $mu1
#> [1] 6.928445
#> 
#> $se1
#> [1] 0.3482651
#> 
#> $n1
#> [1] 5
#> 
#> $mu2
#> [1] 4.388392
#> 
#> $se2
#> [1] 0.8231683
#> 
#> $n2
#> [1] 400
#> 
#> $delta
#> [1] 2.540053
#> 
#> $se_diff
#> [1] 0.893809
#> 
#> $t_stat
#> [1] 2.84183
#> 
#> $df
#> [1] 132.1808
#> 
#> $p_value
#> [1] 0.005197391
```

Compare every cohort the same way with ordinary R:

``` r
test_results <- map_dfr(
  expression_estimates$group,
  function(group) {
    cohort_estimate <- expression_estimates[
      expression_estimates$group == group,
      ,
      drop = FALSE
    ]
    posterior_summary <- summarize_posterior_draws(
      posterior_draws,
      value = cohort_estimate$log_mu
    )
    test <- welch_test_means(
      cohort_estimate$log_mu,
      cohort_estimate$se,
      posterior_summary$mean,
      posterior_summary$sd,
      n1 = cohort_estimate$n,
      n2 = posterior_summary$n
    )
    data.frame(
      group = group,
      log_mu = test$mu1,
      se = test$se1,
      hca_mean = test$mu2,
      hca_sd = test$se2,
      p_value = test$p_value,
      stringsAsFactors = FALSE
    )
  }
)
test_results
#>          group   log_mu        se hca_mean    hca_sd     p_value
#> 1         CTRL 3.301308 0.5652071 4.388392 0.8231683 0.281063450
#> 2         SAVI 6.928445 0.3482651 4.388392 0.8231683 0.005197391
#> 3 SAVI_treated 6.283036 0.3697428 4.388392 0.8231683 0.037975059
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
  subtitle = "QL cohort estimates",
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
#>  [1] purrr_1.2.2          org.Hs.eg.db_3.23.1  AnnotationDbi_1.74.0
#>  [4] IRanges_2.46.0       S4Vectors_0.50.2     Biobase_2.72.0      
#>  [7] BiocGenerics_0.58.1  generics_0.1.4       Seurat_5.5.1        
#> [10] SeuratObject_5.4.0   sp_2.2-3             posteriorHCA_0.2.0  
#> [13] testthat_3.3.2      
#> 
#> loaded via a namespace (and not attached):
#>   [1] RcppAnnoy_0.0.23            splines_4.6.1              
#>   [3] later_1.4.8                 tibble_3.3.1               
#>   [5] polyclip_1.10-7             brms_2.23.0                
#>   [7] fastDummies_1.7.6           lifecycle_1.0.5            
#>   [9] StanHeaders_2.39.1          edgeR_4.10.4               
#>  [11] rprojroot_2.1.1             vroom_1.7.1                
#>  [13] processx_3.9.0              globals_0.19.1             
#>  [15] sccomp_2.4.0                lattice_0.22-9             
#>  [17] MASS_7.3-65                 backports_1.5.1            
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
#>  [39] Rtsne_0.17                  tensorA_0.36.2.1           
#>  [41] inline_0.3.21               ggrepel_0.9.8              
#>  [43] irlba_2.3.7                 listenv_1.0.0              
#>  [45] spatstat.utils_3.2-4        goftest_1.2-3              
#>  [47] RSpectra_0.16-2             spatstat.random_3.5-1      
#>  [49] bridgesampling_1.2-1        fitdistrplus_1.2-6         
#>  [51] parallelly_1.48.0           DelayedArray_0.38.2        
#>  [53] codetools_0.2-20            tidyselect_1.2.1           
#>  [55] bayesplot_1.16.0            farver_2.1.2               
#>  [57] matrixStats_1.5.0           spatstat.explore_3.8-2     
#>  [59] Seqinfo_1.2.0               jsonlite_2.0.0             
#>  [61] ellipsis_0.3.3              progressr_1.0.0            
#>  [63] ggridges_0.5.7              survival_3.8-6             
#>  [65] tools_4.6.1                 ica_1.0-3                  
#>  [67] Rcpp_1.1.2                  glue_1.8.1                 
#>  [69] SparseArray_1.12.2          gridExtra_2.3.1            
#>  [71] qs2_0.3.1                   xfun_0.60                  
#>  [73] cmdstanr_0.9.0              MatrixGenerics_1.24.0      
#>  [75] distributional_0.8.1        usethis_3.2.1              
#>  [77] dplyr_1.2.1                 withr_3.0.3                
#>  [79] loo_2.10.1                  instantiate_0.2.3          
#>  [81] fastmap_1.2.0               callr_3.8.0                
#>  [83] digest_0.6.39               R6_2.6.1                   
#>  [85] mime_0.13                   scattermore_1.2            
#>  [87] tensor_1.5.1                spatstat.data_3.1-9        
#>  [89] RSQLite_3.53.3              tidyr_1.3.2                
#>  [91] data.table_1.18.6.1         S4Arrays_1.12.0            
#>  [93] httr_1.4.9                  htmlwidgets_1.6.4          
#>  [95] uwot_0.2.5                  pkgconfig_2.0.3            
#>  [97] gtable_0.3.6                blob_1.3.0                 
#>  [99] lmtest_0.9-40               S7_0.2.2                   
#> [101] SingleCellExperiment_1.34.0 XVector_0.52.0             
#> [103] brio_1.1.5                  htmltools_0.5.9            
#> [105] dotCall64_1.2               scales_1.4.0               
#> [107] png_0.1-9                   posterior_1.7.0            
#> [109] spatstat.univar_3.2-0       knitr_1.52                 
#> [111] rstudioapi_0.19.0           tzdb_0.5.0                 
#> [113] reshape2_1.4.5              coda_0.19-4.1              
#> [115] checkmate_2.3.4             nlme_3.1-169               
#> [117] cachem_1.1.0                zoo_1.9-0                  
#> [119] stringr_1.6.0               KernSmooth_2.23-26         
#> [121] parallel_4.6.1              miniUI_0.1.2               
#> [123] desc_1.4.3                  pillar_1.11.1              
#> [125] grid_4.6.1                  vctrs_0.7.3                
#> [127] RANN_2.6.3                  promises_1.5.0             
#> [129] stringfish_0.19.2           xtable_1.8-8               
#> [131] cluster_2.1.8.2             evaluate_1.0.5             
#> [133] readr_2.2.0                 locfit_1.5-9.12            
#> [135] mvtnorm_1.4-2               cli_3.6.6                  
#> [137] compiler_4.6.1              rlang_1.3.0                
#> [139] crayon_1.5.3                rstantools_2.7.1           
#> [141] future.apply_1.20.2         labeling_0.4.3             
#> [143] forcats_1.0.1               plyr_1.8.9                 
#> [145] fs_2.1.0                    rstan_2.32.7               
#> [147] stringi_1.8.9               QuickJSR_1.11.0            
#> [149] viridisLite_0.4.3           deldir_2.0-4               
#> [151] Biostrings_2.80.2           devtools_2.5.2             
#> [153] spatstat.geom_3.8-2         Brobdingnag_1.2-9          
#> [155] Matrix_1.7-5                RcppHNSW_0.7.0             
#> [157] hms_1.1.4                   patchwork_1.3.2            
#> [159] bit64_4.8.6                 future_1.75.0              
#> [161] ggplot2_4.0.3               statmod_1.5.2              
#> [163] KEGGREST_1.52.2             shiny_1.14.0               
#> [165] SummarizedExperiment_1.42.0 ROCR_1.0-12                
#> [167] igraph_2.3.3                memoise_2.0.1              
#> [169] RcppParallel_6.2.1          bit_4.6.0
```
