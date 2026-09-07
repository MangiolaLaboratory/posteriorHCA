Cohort expression workflow (wrappers)
================
Chen Zhan
2026-09-07

SAVI case study: *ADRB2* (`ENSG00000169252`) in disease-associated
monocytes.

This vignette is the **high-level wrapper** path. Each wrapper only
composes core functions; it does not reimplement TMM, edgeR QL, brms
draws, or Welch mathematics. The transparent core path is
`vignette("cohort-expression-core", package = "posteriorHCA")`.

``` text
HIGH-LEVEL WRAPPER                  CORE FUNCTIONS CALLED
──────────────────                  ─────────────────────

scale_to_hca_reference()
    ├── load_reference_sample()
    ├── merge_with_reference_sample()
    └── calculate_tmm_offset()

estimate_cohort_logmu()
    ├── extract counts / metadata
    ├── stats::model.matrix()
    └── estimate_logmu_ql()

expression_baseline_draws()
    ├── load_expression_fit()
    ├── build_newdata_grid()
    └── expression_draws()

compare_cohort_to_hca()
    ├── summarize_posterior_draws()
    └── welch_test_means()
```

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

## Map symbols to Ensembl ids

Atlas models are keyed by ENSG. Map once with AnnotationDbi (outside
posteriorHCA), then rebuild a Seurat object that keeps `Category`.

``` r
data(savi_mono)

counts <- as.matrix(Seurat::GetAssayData(savi_mono, layer = "counts"))
mapped <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = rownames(counts),
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first"
)
keep <- !is.na(mapped) & !duplicated(mapped)
counts <- counts[keep, , drop = FALSE]
rownames(counts) <- unname(mapped[keep])

meta <- savi_mono[[]][colnames(counts), , drop = FALSE]
savi_ensg <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta)
```

## Wrapper workflow

`scale_to_hca_reference()` is a convenience wrapper around
`load_reference_sample()` (when `reference` is a cell-type string),
`merge_with_reference_sample()`, and `calculate_tmm_offset()`.

``` r
scaled <- scale_to_hca_reference(
  savi_ensg,
  reference = cell_type
)
```

`estimate_cohort_logmu()` is a convenience wrapper that constructs the
design with `model.matrix()` and delegates estimation to
`estimate_logmu_ql()`.

``` r
cohort_estimates <- estimate_cohort_logmu(
  scaled,
  formula = ~ 0 + Category,
  gene_ensg = gene_ensg
)
cohort_estimates
#>                  gene                group n   log_mu         mu        se
#> 2617  ENSG00000169252         CategoryCTRL 7 3.304671   27.23958 0.6015499
#> 12234 ENSG00000169252         CategorySAVI 5 6.930949 1023.46522 0.3661802
#> 21851 ENSG00000169252 CategorySAVI_treated 5 6.284042  535.95061 0.3896259
#>            df dispersion
#> 2617  19.0239  0.2235231
#> 12234 19.0239  0.2235231
#> 21851 19.0239  0.2235231
```

Intercept-only design (`~ 1`): subset to one `Category`, then call
`estimate_cohort_logmu()` with `formula = ~ 1`. The intercept is that
cohort’s absolute log(μ). Offsets were already computed by
`scale_to_hca_reference()`, so the reference library is not needed in
the subset.

``` r
cohort_estimates_by_level <- map_dfr(
  unique(as.character(scaled$Category[scaled$sample_role == "user"])),
  function(category) {
    out <- estimate_cohort_logmu(
      subset(scaled, sample_role == "user" & Category == category),
      formula = ~ 1,
      gene_ensg = gene_ensg
    )
    out$group <- category
    out
  }
)
cohort_estimates_by_level
#>              gene        group n   log_mu        mu        se        df
#> 1 ENSG00000169252         CTRL 7 3.326279  27.83459 0.3974670 10.928782
#> 2 ENSG00000169252         SAVI 5 6.899991 992.26549 0.2464727  4.013446
#> 3 ENSG00000169252 SAVI_treated 5 6.291712 540.07736 0.2111886  3.954086
#>   dispersion
#> 1  0.2302444
#> 2  0.3381536
#> 3  0.2077055
```

`expression_baseline_draws()` delegates model loading, query-grid
construction and posterior generation to `load_expression_fit()`,
`build_newdata_grid()` and `expression_draws()`.

``` r
hca_draws <- expression_baseline_draws(
  cell_type = cell_type,
  gene_ensg = gene_ensg,
  disease_groups = "Normal",
  tissue_groups = "blood",
  assay_groups = "10x Genomics 3"
)
```

`compare_cohort_to_hca()` summarizes HCA posterior draws with
`summarize_posterior_draws()` and delegates the actual statistical test
to `welch_test_means()`.

``` r
test_results <- compare_cohort_to_hca(
  cohort_estimates,
  hca_draws
)
test_results
#>              gene                group   log_mu        se n hca_mean    hca_sd
#> 1 ENSG00000169252         CategoryCTRL 3.304671 0.6015499 7 4.444192 0.8688719
#> 2 ENSG00000169252         CategorySAVI 6.930949 0.3661802 5 4.444192 0.8688719
#> 3 ENSG00000169252 CategorySAVI_treated 6.284042 0.3896259 5 4.444192 0.8688719
#>   hca_n     delta   se_diff    t_stat        df     p_value
#> 1   400 -1.139521 1.0567879 -1.078287  53.63923 0.285730957
#> 2   400  2.486758 0.9428819  2.637401 133.43381 0.009346491
#> 3   400  1.839850 0.9522325  1.932144 114.35430 0.055815505
```

## Plots

``` r
plot_hca_draws(
  draws = hca_draws,
  subtitle = "Normal, 10x Genomics 3 healthy baseline"
)
```

![](cohort-expression-wrappers_files/figure-gfm/plot-hca-1.png)<!-- -->

``` r
plot_cohort_vs_hca(
  hca_draws = hca_draws,
  cohort_est = test_results,
  subtitle = "QL cohort estimates (wrapper path)",
  annotate = c("group", "p_value")
)
```

![](cohort-expression-wrappers_files/figure-gfm/plot-cohort-1.png)<!-- -->

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
