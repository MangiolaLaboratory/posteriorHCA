Cohort expression workflow (native edgeR)
================
Chen Zhan
2026-09-09

SAVI case study: *ADRB2* (`ENSG00000169252`) in disease-associated
monocytes.

Same analysis as `vignette("cohort-expression-core")`, but with the
edgeR / limma steps written out. Only `load_reference_sample()` and
`merge_with_reference_sample()` are used from posteriorHCA before the
HCA model query.

``` text
USER + HCA counts
  → TMMwsp
  → offset = log(E_user)
  → glmQLFit → makeContrasts → glmQLFTest
  → estimate + QLF SE
  → + log(E_HCA)
  → compare to HCA posterior
```

## Setup

``` r
library(posteriorHCA)
library(Seurat)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(edgeR)
library(limma)

cell_type <- "monocytic"
gene <- "ENSG00000169252"
ref_name <- "___hca_ref_sample___"
```

## Counts

``` r
data(savi_mono, overwrite = TRUE)

counts <- as.matrix(Seurat::GetAssayData(savi_mono, layer = "counts"))
rn <- rownames(counts)

if (!all(grepl("^ENSG", rn))) {
  known <- rn %in% AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")
  counts <- counts[known, , drop = FALSE]
  ensembl <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = rownames(counts),
    column = "ENSEMBL",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  ok <- !is.na(ensembl) & !duplicated(ensembl)
  counts <- counts[ok, , drop = FALSE]
  rownames(counts) <- unname(ensembl[ok])
}

savi_mono <- Seurat::CreateSeuratObject(
  counts = counts,
  meta.data = savi_mono@meta.data
)
```

## HCA reference

``` r
ref <- load_reference_sample(cell_type)
counts <- merge_with_reference_sample(
  counts,
  reference = ref$counts,
  reference_name = ref_name
)
```

The reference is only for joint TMM. Drop it before fitting.

## TMMwsp scaling

``` r
nf <- calcNormFactors(
  counts,
  refColumn = match(ref_name, colnames(counts)),
  method = "TMMwsp"
)
eff <- colSums(counts) * nf

log_E <- log(eff[colnames(counts) != ref_name])
log_E_hca <- log(eff[[ref_name]])
```

## USER model

``` r
y <- counts[, colnames(counts) != ref_name, drop = FALSE]

offset <- matrix(log_E, nrow = nrow(y), ncol = ncol(y), byrow = TRUE)
design <- model.matrix(
  ~ 0 + Category,
  data = savi_mono@meta.data[colnames(y), , drop = FALSE]
)

fit <- glmQLFit(y, design = design, offset = offset, prior.count = 0, robust = TRUE)

contrast <- makeContrasts(CategorySAVI, levels = design)
```

`CategorySAVI` is the absolute SAVI group log mean under
`~ 0 + Category`.

## Estimate and QLF-implied SE

`glmQLFTest()` reports `logFC` on the log2 scale. Convert to natural
log. The SE below is implied by the one-dimensional QL F statistic, not
a coefficient-covariance SE.

``` r
qlf <- glmQLFTest(fit, contrast = contrast)

estimate <- qlf$table[gene, "logFC"] * log(2)
se_qlf <- abs(qlf$table[gene, "logFC"]) / sqrt(qlf$table[gene, "F"]) * log(2)

# Absolute estimand → place on HCA scale (SE unchanged).
log_mu <- estimate + log_E_hca
```

## HCA posterior

``` r
fit_hca <- load_expression_fit(cell_type = cell_type, gene_ensg = gene)
newdata <- build_newdata_grid(
  fit_hca,
  disease_groups = "Normal",
  tissue_groups = "blood",
  assay_groups = "10x Genomics 3"
)
draws <- expression_draws(
  fit_hca,
  newdata = newdata,
  quantity = "linpred",
  marginalise = "mean"
)
```

## Plot

``` r
plot_hca_draws(
  draws,
  query_mu = log_mu,
  query_SE = se_qlf,
  query_label = "SAVI"
)
```

![](cohort-expression-edger-native_files/figure-gfm/plot-hca-1.png)<!-- -->

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
#>  [1] edgeR_4.10.1         limma_3.68.4         org.Hs.eg.db_3.23.1 
#>  [4] AnnotationDbi_1.74.0 IRanges_2.46.0       S4Vectors_0.50.1    
#>  [7] Biobase_2.72.0       BiocGenerics_0.58.1  generics_0.1.4      
#> [10] Seurat_5.5.1         SeuratObject_5.4.0   sp_2.2-3            
#> [13] posteriorHCA_0.2.0   testthat_3.3.2      
#> 
#> loaded via a namespace (and not attached):
#>   [1] RcppAnnoy_0.0.23            splines_4.6.1              
#>   [3] later_1.4.8                 tibble_3.3.1               
#>   [5] polyclip_1.10-7             brms_2.23.0                
#>   [7] fastDummies_1.7.6           lifecycle_1.0.5            
#>   [9] StanHeaders_2.39.1          rprojroot_2.1.1            
#>  [11] vroom_1.7.1                 processx_3.9.0             
#>  [13] sccomp_2.4.0                globals_0.19.1             
#>  [15] lattice_0.22-9              MASS_7.3-65                
#>  [17] backports_1.5.1             magrittr_2.0.5             
#>  [19] plotly_4.12.1               rmarkdown_2.32             
#>  [21] yaml_2.3.12                 httpuv_1.6.17              
#>  [23] otel_0.2.0                  sctransform_0.4.3          
#>  [25] spam_2.11-4                 sessioninfo_1.2.4          
#>  [27] pkgbuild_1.4.8              spatstat.sparse_3.2-0      
#>  [29] reticulate_1.47.0           cowplot_1.2.0              
#>  [31] pbapply_1.7-5               DBI_1.3.0                  
#>  [33] RColorBrewer_1.1-3          pkgload_1.5.3              
#>  [35] multcomp_1.4-32             abind_1.4-8                
#>  [37] GenomicRanges_1.64.0        Rtsne_0.17                 
#>  [39] purrr_1.2.2                 TH.data_1.1-5              
#>  [41] tensorA_0.36.2.1            sandwich_3.1-3             
#>  [43] inline_0.3.21               ggrepel_0.9.8              
#>  [45] irlba_2.3.7                 listenv_1.0.0              
#>  [47] spatstat.utils_3.2-4        goftest_1.2-3              
#>  [49] RSpectra_0.16-2             spatstat.random_3.5-1      
#>  [51] bridgesampling_1.2-1        fitdistrplus_1.2-6         
#>  [53] parallelly_1.48.0           DelayedArray_0.38.2        
#>  [55] codetools_0.2-20            tidyselect_1.2.1           
#>  [57] bayesplot_1.16.0            farver_2.1.2               
#>  [59] matrixStats_1.5.0           spatstat.explore_3.8-2     
#>  [61] Seqinfo_1.2.0               jsonlite_2.0.0             
#>  [63] ellipsis_0.3.3              progressr_1.0.0            
#>  [65] ggridges_0.5.7              survival_3.8-6             
#>  [67] emmeans_2.0.4               tools_4.6.1                
#>  [69] ica_1.0-3                   Rcpp_1.1.2                 
#>  [71] glue_1.8.1                  SparseArray_1.12.2         
#>  [73] gridExtra_2.3.1             qs2_0.3.1                  
#>  [75] xfun_0.60                   cmdstanr_0.9.0             
#>  [77] MatrixGenerics_1.24.0       usethis_3.2.1              
#>  [79] distributional_0.8.1        dplyr_1.2.1                
#>  [81] withr_3.0.3                 loo_2.10.1                 
#>  [83] instantiate_0.2.3           fastmap_1.2.0              
#>  [85] callr_3.8.0                 digest_0.6.39              
#>  [87] R6_2.6.1                    mime_0.13                  
#>  [89] estimability_2.0.0          scattermore_1.2            
#>  [91] tensor_1.5.1                spatstat.data_3.1-9        
#>  [93] RSQLite_3.53.3              tidyr_1.3.2                
#>  [95] data.table_1.18.6.1         S4Arrays_1.12.0            
#>  [97] httr_1.4.9                  htmlwidgets_1.6.4          
#>  [99] uwot_0.2.5                  pkgconfig_2.0.3            
#> [101] gtable_0.3.6                blob_1.3.0                 
#> [103] lmtest_0.9-40               S7_0.2.2                   
#> [105] SingleCellExperiment_1.34.0 XVector_0.52.0             
#> [107] brio_1.1.5                  htmltools_0.5.9            
#> [109] dotCall64_1.2               scales_1.4.0               
#> [111] png_0.1-9                   posterior_1.7.0            
#> [113] spatstat.univar_3.2-0       knitr_1.51                 
#> [115] rstudioapi_0.19.0           tzdb_0.5.0                 
#> [117] reshape2_1.4.5              coda_0.19-4.1              
#> [119] checkmate_2.3.4             nlme_3.1-169               
#> [121] cachem_1.1.0                zoo_1.9-0                  
#> [123] stringr_1.6.0               KernSmooth_2.23-26         
#> [125] parallel_4.6.1              miniUI_0.1.2               
#> [127] desc_1.4.3                  pillar_1.11.1              
#> [129] grid_4.6.1                  vctrs_0.7.3                
#> [131] RANN_2.6.3                  promises_1.5.0             
#> [133] stringfish_0.19.2           xtable_1.8-8               
#> [135] cluster_2.1.8.2             evaluate_1.0.5             
#> [137] readr_2.2.0                 mvtnorm_1.4-2              
#> [139] cli_3.6.6                   locfit_1.5-9.12            
#> [141] compiler_4.6.1              rlang_1.3.0                
#> [143] crayon_1.5.3                rstantools_2.7.1           
#> [145] future.apply_1.20.2         labeling_0.4.3             
#> [147] forcats_1.0.1               plyr_1.8.9                 
#> [149] fs_2.1.0                    rstan_2.32.7               
#> [151] stringi_1.8.9               QuickJSR_1.11.0            
#> [153] viridisLite_0.4.3           deldir_2.0-4               
#> [155] Biostrings_2.80.1           devtools_2.5.2             
#> [157] spatstat.geom_3.8-2         Brobdingnag_1.2-9          
#> [159] Matrix_1.7-5                RcppHNSW_0.7.0             
#> [161] hms_1.1.4                   patchwork_1.3.2            
#> [163] bit64_4.8.6                 future_1.75.0              
#> [165] ggplot2_4.0.3               KEGGREST_1.52.2            
#> [167] statmod_1.5.2               shiny_1.14.0               
#> [169] SummarizedExperiment_1.42.0 ROCR_1.0-12                
#> [171] igraph_2.3.3                memoise_2.0.1              
#> [173] RcppParallel_6.2.1          bit_4.6.0
```
