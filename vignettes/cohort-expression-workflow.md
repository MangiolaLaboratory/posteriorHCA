Cohort expression workflow
================
Chen Zhan
2026-09-04

This vignette follows `examples/savi_adrb2_workflow.R` for the SAVI case
study (*ADRB2* in disease-associated monocytes, PBMC / blood).

At each modelling step we show **two equivalent paths**:

- **Core / manual** — matrix-level helpers you can wire yourself  
- **Wrapper** — same cores, with Seurat / SE / gene-id / metadata
  handling

Pipeline:

1.  Prepare pseudobulk counts (`savi_mono`) and harmonise gene ids  
2.  Align to the HCA monocytic reference
    - core: `merge_with_reference_sample()` + `calculate_tmm_offset()`  
    - wrapper: `scale_to_hca_reference()`  
3.  Estimate cohort log(μ)
    - core: `design_from_formula()` + `estimate_logmu_ql()` (all
      genes)  
    - wrapper: `estimate_cohort_logmu()` (filter genes)  
    - demo: single-cohort `~ 1` (intercept as cohort log(μ))  
4.  Healthy HCA posterior (Normal, blood, 10x Genomics 3)
    - manual: `load_expr_fit()` → `build_newdata_grid()` →
      `expr_draws()`  
    - wrapper: `expr_predict()`  
5.  Welch test vs HCA
    - wrapper: `welch_t_test_cohort_hca()`  
    - manual: `cohort_estimate_at()` + `summarize_posterior_draws()` +
      `welch_test_means()`  
6.  Plot results

## Setup

``` r
library(posteriorHCA)
library(Seurat)
library(dplyr)
library(ggplot2)

cell_type <- "monocytic"
```

## 1. Prepare `savi_mono`

The SAVI study ([de Cevins *et al.*,
2023](https://doi.org/10.1016/j.xcrm.2023.101333)) is PBMC scRNA-seq
pseudobulked by sample and cell type. We subset **disease-associated
monocytes** (17 samples: CTRL, SAVI, SAVI_treated).

``` r
savi_path <- "GSE226598_SAVI_pseudobulk_Sample_CellType.rds"
savi <- readRDS(savi_path)
savi_mono <- subset(savi, subset = CellType == "17. Disease-associated monocytes")
```

``` r
data(savi_mono)
dim(savi_mono)
#> [1] 24431    17
table(savi_mono$Category)
#> 
#>         CTRL         SAVI SAVI_treated 
#>            7            5            5
```

Harmonise gene symbols to Ensembl IDs so cohort counts and HCA models
share row names:

``` r
savi_mono <- harmonise_gene_ids(savi_mono, id_type = "auto")
```

## 2. Align to the HCA monocytic reference - step 0

Reference sits at offset 0 (same convention as atlas brms training:
`offset = log(1 / multiplier)`).

### 2a. Core functions (matrix)

``` r
ref <- load_reference_sample(cell_type)
user_mat <- as.matrix(Seurat::GetAssayData(savi_mono, layer = "counts"))

combined <- merge_with_reference_sample(
  user_mat,
  reference = ref$counts,
  reference_name = ref$sample_id
)
scaling <- calculate_tmm_offset(
  combined,
  reference_name = ref$sample_id,
  method = "TMMwsp"
)

c(
  n_genes = nrow(combined),
  n_samples = ncol(combined),
  reference_offset = unname(scaling$offset[[ref$sample_id]])
)
#>          n_genes        n_samples reference_offset 
#>             9625               18                0
```

### 2b. Wrapper `scale_to_hca_reference()`

``` r
aligned <- scale_to_hca_reference(savi_mono, cell_type)

table(aligned$sample_role)
#> 
#> reference      user 
#>         1        17
aligned$hca_reference_name[1]
#>                     C1_17. Disease-associated monocytes 
#> "e11a0d767c2a97f658791b17cae25860____SC142___monocytic"
```

``` r
all.equal(
  unname(scaling$offset[colnames(user_mat)]),
  unname(aligned$hca_offset[colnames(user_mat)]),
  tolerance = 1e-10
)
#> [1] TRUE
```

## 3. Estimate cohort log(μ) - step 1

Absolute log(μ) comes from QL coefficients with `prior.count = 0` on the
TMM offset scale. Prefer a cell-means design such as `~ 0 + Category`.

### 3a. Core `design_from_formula()` + `estimate_logmu_ql()`

`estimate_logmu_ql()` fits **all genes** once and returns log(μ) + SE
for every gene × design column. Filter afterwards if you only need a
subset.

``` r
meta_core <- data.frame(
  Category = c(
    as.character(savi_mono$Category[colnames(user_mat)]),
    "reference"
  ),
  row.names = colnames(combined),
  stringsAsFactors = FALSE
)
meta_core$Category <- factor(meta_core$Category)

design <- design_from_formula(~ 0 + Category, meta_core)
est_all <- estimate_logmu_ql(
  counts = combined,
  offset = scaling$offset,
  design = design,
  cell_type = cell_type
)
gene_ensg <- resolve_gene_one("ADRB2", cell_type = cell_type)
est_core <- est_all[est_all$gene == gene_ensg, , drop = FALSE]
est_core$gene_symbol <- "ADRB2"
est_core
#>                  gene gene_symbol cell_type        group n   log_mu         mu
#> 2620  ENSG00000169252       ADRB2 monocytic         CTRL 7 3.313065   27.46918
#> 12245 ENSG00000169252       ADRB2 monocytic    reference 1 4.530958   92.84749
#> 21870 ENSG00000169252       ADRB2 monocytic         SAVI 5 6.937536 1030.22874
#> 31495 ENSG00000169252       ADRB2 monocytic SAVI_treated 5 6.286779  537.41932
#>              se       df dispersion
#> 2620  0.6154281 18.98187  0.2226158
#> 12245 0.8046448 18.98187  0.2226158
#> 21870 0.3751294 18.98187  0.2226158
#> 31495 0.3992543 18.98187  0.2226158
```

### 3b. Wrapper `estimate_cohort_logmu()`

Same cores, plus container I/O and gene-id resolution. `genes` only
filters the returned table (dispersion still uses the full matrix).

``` r
est <- estimate_cohort_logmu(
  aligned,
  formula = ~ 0 + Category,
  genes = "ADRB2"
)
est
#>              gene gene_symbol cell_type        group n   log_mu         mu
#> 1 ENSG00000169252       ADRB2 monocytic         CTRL 7 3.313065   27.46918
#> 2 ENSG00000169252       ADRB2 monocytic         SAVI 5 6.937536 1030.22874
#> 3 ENSG00000169252       ADRB2 monocytic SAVI_treated 5 6.286779  537.41932
#> 4 ENSG00000169252       ADRB2 monocytic    reference 1 4.530958   92.84749
#>          se       df dispersion
#> 1 0.6154281 18.98187  0.2226158
#> 2 0.3751294 18.98187  0.2226158
#> 3 0.3992543 18.98187  0.2226158
#> 4 0.8046448 18.98187  0.2226158
```

``` r
all.equal(
  est_core$log_mu[match(est$group, est_core$group)],
  est$log_mu,
  tolerance = 1e-8
)
#> [1] TRUE
```

### 3c. Demo: single-cohort `~ 1` (intercept = cohort log(μ))

Alternative to `~ 0 + Category` on the full object. Subset to one
`Category`, fit NB/QL with `formula = ~ 1`, and treat the intercept as
that cohort’s absolute log(μ). Useful when you want per-cohort fits
without a multi-level design matrix.

``` r
est_by_level <- purrr::map_dfr(
  levels(aligned$Category),
  .f = function(x) {
    estimate_cohort_logmu(
      subset(aligned, Category == x),
      formula = ~ 1,
      genes = "ADRB2"
    ) %>%
      dplyr::mutate(group = x)
  }
)
est_by_level
#>              gene gene_symbol cell_type        group n   log_mu         mu
#> 1 ENSG00000169252       ADRB2 monocytic         CTRL 7 3.334290   28.05845
#> 2 ENSG00000169252       ADRB2 monocytic         SAVI 5 6.908769 1001.01467
#> 3 ENSG00000169252       ADRB2 monocytic SAVI_treated 5 6.292121  540.29813
#>          se        df dispersion
#> 1 0.3984058 10.899835  0.2300399
#> 2 0.2621220  3.994552  0.3382149
#> 3 0.2118680  3.953785  0.2076096
```

## 4. Healthy HCA baseline draws

Covariates fixed to healthy blood 10x Genomics 3 (matches SAVI tissue).
`quantity = "linpred"` is log(μ) on the natural-log scale.

### 4a. Manual coding

``` r
fit <- load_expr_fit(cell_type = cell_type, gene = "ADRB2") # step 2

grid <- build_newdata_grid(
  fit,
  disease_groups = "Normal",
  tissue_groups = "blood",
  assay_groups = "10x Genomics 3"
)

hca_res <- expr_draws(  # step 3
  fit,
  newdata = grid,
  quantity = "linpred",
  collapse = "mean"
)

round(c(mean = mean(hca_res$draws), sd = sd(hca_res$draws)), 3)
#>  mean    sd 
#> 4.444 0.869
```

### 4b. Wrapper `expr_predict()`

``` r
hca_pred <- expr_predict(
  fit = fit,
  disease_groups = "Normal",
  tissue_groups = "blood",
  assay_groups = "10x Genomics 3",
  quantity = "linpred",
  collapse = "mean"
)
```

## 5. Test cohorts against the healthy baseline

### 5a. Wrapper `welch_t_test_cohort_hca()`

``` r
test_results <- welch_t_test_cohort_hca(
  cohort_est = est,
  hca_draws = hca_pred
)
test_results
#>              gene gene_symbol cell_type       cohort method cohort_log_mu
#> 1 ENSG00000169252       ADRB2 monocytic         CTRL     ql      3.313065
#> 2 ENSG00000169252       ADRB2 monocytic         SAVI     ql      6.937536
#> 3 ENSG00000169252       ADRB2 monocytic SAVI_treated     ql      6.286779
#>   cohort_se hca_mean    hca_sd delta_log_mu  se_diff     t_stat        df
#> 1 0.6154281 4.373241 0.9304786    -1.060177 1.115591 -0.9503276  60.06341
#> 2 0.3751294 4.373241 0.9304786     2.564295 1.003251  2.5559858 148.34013
#> 3 0.3992543 4.373241 0.9304786     1.913538 1.012519  1.8898786 127.68952
#>      p_value empirical_rank           direction
#> 1 0.34575423         0.1150 consistent_with_hca
#> 2 0.01159443         0.9975           above_hca
#> 3 0.06104119         0.9900 consistent_with_hca
```

### 5b. Manual helpers (one cohort)

``` r
cohort <- cohort_estimate_at(est, group = "SAVI")
baseline <- summarize_posterior_draws(hca_res, value = cohort$mu) # step 4

welch_test_means( # step 5
  cohort$mu, cohort$se,
  baseline$mean, baseline$sd,
  n1 = cohort$n, n2 = baseline$n
)
#> $mu1
#> [1] 6.937536
#> 
#> $se1
#> [1] 0.3751294
#> 
#> $n1
#> [1] 5
#> 
#> $mu2
#> [1] 4.444192
#> 
#> $se2
#> [1] 0.8688719
#> 
#> $n2
#> [1] 400
#> 
#> $delta
#> [1] 2.493344
#> 
#> $se_diff
#> [1] 0.9463934
#> 
#> $t_stat
#> [1] 2.634575
#> 
#> $df
#> [1] 125.7561
#> 
#> $p_value
#> [1] 0.009483039
```

## 6. Visualise results

``` r
plot_hca_draws(
  draws = hca_res,
  subtitle = "Normal, 10x Genomics 3 healthy baseline"
)
```

![](cohort-expression-workflow_files/figure-gfm/plot-hca-1.png)<!-- -->

``` r
plot_cohort_vs_hca(
  hca_draws = hca_res,
  test_results = test_results,
  subtitle = "QL cohort estimates",
  annotate = c("group", "p_value", "direction")
)
```

![](cohort-expression-workflow_files/figure-gfm/plot-cohort-1.png)<!-- -->

``` r
plot_hca_draws(
  draws = hca_pred,
  title = "ADRB2 predicted count posterior (healthy HCA)"
)
```

![](cohort-expression-workflow_files/figure-gfm/plot-hca-pred-1.png)<!-- -->

## Session info

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-conda-linux-gnu
#> Running under: Red Hat Enterprise Linux 8.4 (Ootpa)
#> 
#> Matrix products: default
#> BLAS/LAPACK: /home/a1237163/miniconda3/envs/R_env/lib/libopenblasp-r0.3.29.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: Australia/Adelaide
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] ggplot2_4.0.3      dplyr_1.2.1        Seurat_5.5.0       SeuratObject_5.4.0
#> [5] sp_2.2-1           posteriorHCA_0.2.0 testthat_3.3.2    
#> 
#> loaded via a namespace (and not attached):
#>   [1] RcppAnnoy_0.0.23            splines_4.5.3              
#>   [3] later_1.4.8                 tibble_3.3.1               
#>   [5] polyclip_1.10-7             brms_2.23.1                
#>   [7] fastDummies_1.7.6           lifecycle_1.0.5            
#>   [9] StanHeaders_2.39.0.9000     edgeR_4.8.2                
#>  [11] rprojroot_2.1.1             vroom_1.7.1                
#>  [13] globals_0.19.1              processx_3.9.0             
#>  [15] sccomp_2.2.0                lattice_0.22-9             
#>  [17] MASS_7.3-65                 backports_1.5.1            
#>  [19] magrittr_2.0.5              limma_3.66.0               
#>  [21] plotly_4.12.0               rmarkdown_2.31             
#>  [23] yaml_2.3.12                 httpuv_1.6.17              
#>  [25] otel_0.2.0                  sctransform_0.4.3          
#>  [27] spam_2.11-3                 spatstat.sparse_3.2-0      
#>  [29] sessioninfo_1.2.3           pkgbuild_1.4.8             
#>  [31] reticulate_1.46.0           DBI_1.3.0                  
#>  [33] cowplot_1.2.0               pbapply_1.7-4              
#>  [35] RColorBrewer_1.1-3          multcomp_1.4-30            
#>  [37] abind_1.4-8                 pkgload_1.5.2              
#>  [39] Rtsne_0.17                  GenomicRanges_1.62.1       
#>  [41] purrr_1.2.2                 BiocGenerics_0.56.0        
#>  [43] TH.data_1.1-5               tensorA_0.36.2.1           
#>  [45] sandwich_3.1-1              IRanges_2.44.0             
#>  [47] S4Vectors_0.48.1            inline_0.3.21              
#>  [49] ggrepel_0.9.8               irlba_2.3.7                
#>  [51] spatstat.utils_3.2-3        listenv_1.0.0              
#>  [53] goftest_1.2-3               RSpectra_0.16-2            
#>  [55] spatstat.random_3.5-0       bridgesampling_1.2-1       
#>  [57] fitdistrplus_1.2-6          parallelly_1.48.0          
#>  [59] codetools_0.2-20            DelayedArray_0.36.1        
#>  [61] tidyselect_1.2.1            bayesplot_1.15.0.9000      
#>  [63] farver_2.1.2                spatstat.explore_3.8-1     
#>  [65] matrixStats_1.5.0           stats4_4.5.3               
#>  [67] Seqinfo_1.0.0               jsonlite_2.0.0             
#>  [69] ellipsis_0.3.3              progressr_0.19.0           
#>  [71] ggridges_0.5.7              survival_3.8-6             
#>  [73] emmeans_2.0.3               tools_4.5.3                
#>  [75] ica_1.0-3                   Rcpp_1.1.2                 
#>  [77] glue_1.8.1                  gridExtra_2.3.1            
#>  [79] SparseArray_1.10.10         qs2_0.2.1                  
#>  [81] xfun_0.57                   cmdstanr_0.9.0             
#>  [83] MatrixGenerics_1.22.0       distributional_0.8.1       
#>  [85] usethis_3.2.1               withr_3.0.3                
#>  [87] loo_2.10.0.9000             instantiate_0.2.3          
#>  [89] fastmap_1.2.0               callr_3.8.0                
#>  [91] digest_0.6.39               R6_2.6.1                   
#>  [93] mime_0.13                   estimability_1.5.1         
#>  [95] scattermore_1.2             tensor_1.5.1               
#>  [97] RSQLite_3.53.1              spatstat.data_3.1-9        
#>  [99] dichromat_2.0-0.1           tidyr_1.3.2                
#> [101] generics_0.1.4              data.table_1.18.4          
#> [103] httr_1.4.8                  htmlwidgets_1.6.4          
#> [105] S4Arrays_1.10.1             uwot_0.2.4                 
#> [107] pkgconfig_2.0.3             gtable_0.3.6               
#> [109] blob_1.3.0                  lmtest_0.9-40              
#> [111] S7_0.2.2                    SingleCellExperiment_1.32.0
#> [113] XVector_0.50.0              brio_1.1.5                 
#> [115] htmltools_0.5.9             dotCall64_1.2              
#> [117] scales_1.4.0                Biobase_2.70.0             
#> [119] png_0.1-9                   posterior_1.7.1            
#> [121] spatstat.univar_3.2-0       knitr_1.51                 
#> [123] rstudioapi_0.18.0           reshape2_1.4.5             
#> [125] tzdb_0.5.0                  curl_7.1.0                 
#> [127] coda_0.19-4.1               checkmate_2.3.4            
#> [129] nlme_3.1-169                org.Hs.eg.db_3.22.0        
#> [131] cachem_1.1.0                zoo_1.8-15                 
#> [133] stringr_1.6.0               KernSmooth_2.23-26         
#> [135] parallel_4.5.3              miniUI_0.1.2               
#> [137] AnnotationDbi_1.72.0        desc_1.4.3                 
#> [139] pillar_1.11.1               grid_4.5.3                 
#> [141] vctrs_0.7.3                 RANN_2.6.2                 
#> [143] promises_1.5.0              stringfish_0.19.0          
#> [145] xtable_1.8-8                cluster_2.1.8.2            
#> [147] evaluate_1.0.5              readr_2.2.0                
#> [149] mvtnorm_1.4-2               cli_3.6.6                  
#> [151] locfit_1.5-9.12             compiler_4.5.3             
#> [153] rlang_1.3.0                 crayon_1.5.3               
#> [155] rstantools_2.6.0.9000       future.apply_1.20.2        
#> [157] labeling_0.4.3              ps_1.9.3                   
#> [159] plyr_1.8.9                  forcats_1.0.1              
#> [161] fs_2.1.0                    rstan_2.39.0.9000          
#> [163] stringi_1.8.9               QuickJSR_1.10.0            
#> [165] deldir_2.0-4                viridisLite_0.4.3          
#> [167] Biostrings_2.78.0           lazyeval_0.2.3             
#> [169] spatstat.geom_3.8-1         devtools_2.5.2             
#> [171] Brobdingnag_1.2-9           Matrix_1.7-5               
#> [173] RcppHNSW_0.7.0              hms_1.1.4                  
#> [175] patchwork_1.3.2             bit64_4.8.2                
#> [177] future_1.75.0               KEGGREST_1.50.0            
#> [179] statmod_1.5.2               shiny_1.14.0               
#> [181] SummarizedExperiment_1.40.0 ROCR_1.0-12                
#> [183] igraph_2.3.1                memoise_2.0.1              
#> [185] RcppParallel_5.1.11-2       bit_4.6.0
```
