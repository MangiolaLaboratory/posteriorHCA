posteriorHCA
================

<!-- badges: start -->

[![Lifecycle:experimental](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

<!-- badges: end -->

`PosteriorHCA` applies **Bayesian inference** to model the comprehensive
**Human Cell Atlas (HCA)**, allowing researchers to probabilistically
analyze both **cellular composition** and **gene expression** patterns
across human samples.

This package provides two main analytical capabilities:

1.  **Cellular Composition Analysis**: Uses the `sccomp` package to
    model variations in cell type proportions across:
    - **Cell Types**
    - **Sex**
    - **Age**
    - **Ethnicity**
    - **Disease conditions**
    - **Tissue Groups**
2.  **Gene Expression Prediction**: Uses `brms` models to predict gene
    expression levels for specific cell types and genes based on sample
    metadata including:
    - **Cell Types**
    - **Sex**
    - **Age Decade**
    - **Disease Groups**
    - **Ethnicity Groups**
    - **Assay Groups**
    - **Tissue Groups**

The probabilistic landscapes are **pre-trained**, meaning they represent
**posterior distributions** estimated from millions of cells.
`PosteriorHCA` makes these **queryable**, allowing researchers to:

- **Test hypotheses** against precomputed reference baselines.
- **Compare new samples** to expected cellular composition
  distributions.
- **Predict gene expression** levels for specific cell types and genes.
- **Assess how cell proportions and gene expression vary** across human
  conditions.

## **Installation**

To install the latest development version of `PosteriorHCA` from GitHub,
run:

``` r
# Install devtools if not already installed
install.packages("devtools")

# Install PosteriorHCA from GitHub
devtools::install_github("MangiolaLaboratory/posteriorHCA")
```

Once installed, load the package:

``` r
library(posteriorHCA)
```

## Usage Examples

### **1. Cellular Composition Analysis**

#### Load relevant packages and example data

The package includes a small example dataset for demonstration:

``` r
library(tidyverse)
library(ggplot2)

data(example_proportions)
print(example_proportions)
#> # A tibble: 26 × 3
#>    sample_id                        cell_type       proportion
#>    <chr>                            <chr>                <dbl>
#>  1 0000c153da22cf963b807c0563aca6a6 b memory          0       
#>  2 0000c153da22cf963b807c0563aca6a6 b naive           0       
#>  3 0000c153da22cf963b807c0563aca6a6 cd14 mono         0       
#>  4 0000c153da22cf963b807c0563aca6a6 cd16 mono         0.00502 
#>  5 0000c153da22cf963b807c0563aca6a6 cd4 fh em         0       
#>  6 0000c153da22cf963b807c0563aca6a6 cd4 naive         0       
#>  7 0000c153da22cf963b807c0563aca6a6 cd4 tcm           0.000359
#>  8 0000c153da22cf963b807c0563aca6a6 cd4 th1 em        0       
#>  9 0000c153da22cf963b807c0563aca6a6 cd4 th1/th17 em   0.000359
#> 10 0000c153da22cf963b807c0563aca6a6 cd4 th17 em       0       
#> # ℹ 16 more rows
```

#### Run Composition Posterior Testing

The `composition_posterior_test()` function allows testing new cellular
composition data against the probabilistic model:

``` r
result <- composition_posterior_test(
    proportions = example_proportions,
    sex = "male",
    age_decade = '4',
    disease_groups = "Normal",
    ethnicity_groups = "European",
    tissue_groups = "blood", 
    load_model_to_global_env = T
)
#> ℹ Using cached model file: /home/a1237163/.cache/R/posteriorHCA/estimates_age_decade___L3___disease_TRUE___immune_only_TRUE.rds
#> Model file loaded: Global Env.
#> All inputs of metadata are valid.
#> Start to generative sampling form posterior...
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 0.0 seconds.
#> # A tibble: 26 × 7
#>    sample_id_observed cell_type proportion_observed Empirical_Confidence    mean
#>    <chr>              <chr>                   <dbl>                <dbl>   <dbl>
#>  1 0000c153da22cf963… b memory             0                           0 0.104  
#>  2 0000c153da22cf963… b naive              0                           0 0.107  
#>  3 0000c153da22cf963… cd14 mono            0                           0 0.216  
#>  4 0000c153da22cf963… cd16 mono            0.00502                     0 0.0195 
#>  5 0000c153da22cf963… cd4 fh em            0                           0 0.00771
#>  6 0000c153da22cf963… cd4 naive            0                           0 0.112  
#>  7 0000c153da22cf963… cd4 tcm              0.000359                    0 0.0199 
#>  8 0000c153da22cf963… cd4 th1 …            0                           0 0.0113 
#>  9 0000c153da22cf963… cd4 th1/…            0.000359                    0 0.0181 
#> 10 0000c153da22cf963… cd4 th17…            0                           0 0.0241 
#> # ℹ 16 more rows
#> # ℹ 2 more variables: lower <dbl>, upper <dbl>
```

![](/README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

#### About `load_model_to_global_env`

In the `composition_posterior_test` function, we added arg -
`load_model_to_global_env` to determine if the model file should be
loaded into global env or local env (scope within the function running)
in R session. By default, it is set to be TRUE -
`load_model_to_global_env = T`, and this is particularly useful when you
wish to repeatedly run the function. Loading model file into global env
prevents re-loading the model file, which usually will cost you minutes
to run.

#### Inspect Composition Results

**Posterior Distribution Table:**

``` r
print(result$result_table)
#> # A tibble: 26 × 7
#>    sample_id_observed cell_type proportion_observed Empirical_Confidence    mean
#>    <chr>              <chr>                   <dbl>                <dbl>   <dbl>
#>  1 0000c153da22cf963… b memory             0                           0 0.104  
#>  2 0000c153da22cf963… b naive              0                           0 0.107  
#>  3 0000c153da22cf963… cd14 mono            0                           0 0.216  
#>  4 0000c153da22cf963… cd16 mono            0.00502                     0 0.0195 
#>  5 0000c153da22cf963… cd4 fh em            0                           0 0.00771
#>  6 0000c153da22cf963… cd4 naive            0                           0 0.112  
#>  7 0000c153da22cf963… cd4 tcm              0.000359                    0 0.0199 
#>  8 0000c153da22cf963… cd4 th1 …            0                           0 0.0113 
#>  9 0000c153da22cf963… cd4 th1/…            0.000359                    0 0.0181 
#> 10 0000c153da22cf963… cd4 th17…            0                           0 0.0241 
#> # ℹ 16 more rows
#> # ℹ 2 more variables: lower <dbl>, upper <dbl>
```

**Density Plot of Posterior Predictions:**

``` r
print(result$plot)
```

![](/README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### **2. Gene Expression Prediction**

The `expr_predict()` function allows predicting gene expression levels
for specific cell types and genes based on sample metadata:

``` r
# Predict gene expression for a specific cell type and gene
expr_result <- expr_predict(
  cell_type = 'cd4 naive',
  gene = "ENSG00000000419",
  age_decade = '7',
  sex = 'female',
  disease_groups = 'Normal',
  ethnicity_groups = 'European',
  assay_groups = '10x Genomics 3',
  tissue_groups = 'blood'
)
#> ℹ Using cached file: /home/a1237163/.cache/R/posteriorHCA/V1/cd4.naive/ENSG00000000419
```

#### Inspect Gene Expression Results

**Summary Statistics:**

``` r
print(expr_result$summary)
#> $mean
#> [1] 97.065
#> 
#> $median
#> [1] 72
#> 
#> $peak_location
#> [1] 49.04881
```

**Density Plot of Predicted Expression:**

``` r
print(expr_result$plot)
```

![](/README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

**Predicted Values:**

``` r
head(expr_result$pred)
#>   value
#> 1    63
#> 2   476
#> 3   110
#> 4    87
#> 5   104
#> 6   239
```

## **Input Arguments and Valid Values**

Both `composition_posterior_test()` and `expr_predict()` functions
accept various metadata inputs corresponding to observed sample
characteristics. These include **sex, age decade, disease condition,
ethnicity, and tissue type, etc**. For composition analysis, users can
also provide observed cell type proportions. However, if some metadata
or proportions are unknown, users can leave these fields empty, and the
functions will still generate meaningful posterior estimates based on
the comprehensive reference.

If users choose to specify metadata, it must be consistent with the
**Human Cell Atlas (HCA) annotations** to ensure accurate comparisons.
Below is a comprehensive list of valid arguments:

- **Cell Types:** “b memory”, “b naive”, “cd14 mono”, “cd16 mono”, “cd4
  fh em”, “cd4 naive”, “cd4 tcm”, “cd4 th1 em”, “cd4 th1/th17 em”, “cd4
  th17 em”, “cd4 th2 em”, “cd8 naive”, “cd8 tcm”, “cd8 tem”, “cdc”,
  “granulocyte”, “ilc”, “macrophage”, “mait”, “mast”, “nk”, “nkt”,
  “pdc”, “plasma”, “tgd”, “treg”

- **Sex:** “female”, “male”, “unknown”

- **Age Decades:** 1, 2, 3, 4, 5, 6, 7, 8, 9, 10

- **Disease Groups:** “Normal”, “COVID-19 related_blood”, “Metabolic and
  Other Disorders: Renal and Urinary Disorders_renal system”, “Systemic
  Lupus Erythematosus_blood”, “Respiratory Conditions: General
  Respiratory Disorders_blood”, “COVID-19 related_nasal, oral, and
  pharyngeal regions”, “COVID-19 related_respiratory system”,
  “Glioblastoma_cerebral lobes and cortical areas”, “Cancer: Hematologic
  Cancer_bone marrow”, “Lung Adenocarcinoma_respiratory system”,
  “Infectious and Immune-related Diseases: Autoimmune and Immune-Related
  Disorders_small intestine”, “Infectious and Immune-related Diseases:
  Respiratory Infections_blood”, “Metabolic and Other Disorders:
  Gastrointestinal Disorders_large intestine”, “Cancer: Renal
  Cancer_renal system”, “Infectious and Immune-related Diseases:
  Autoimmune and Immune-Related Disorders_respiratory system”,
  “Infectious and Immune-related Diseases: Respiratory Infections_nasal,
  oral, and pharyngeal regions”, “Metabolic and Other Disorders: Renal
  and Urinary Disorders_prostate”, “Metabolic and Other
  Disorders_sensory-related structures”, “Cancer: Gastrointestinal
  Cancer_large intestine”, “Other Diseases_brainstem and cerebellar
  structures”, “Metabolic and Other Disorders: Gastrointestinal
  Disorders_oesophagus”, “Cancer: Lung Cancer_respiratory system”, “Lung
  Adenocarcinoma_cerebral lobes and cortical areas”, “Metabolic and
  Other Disorders: Gastrointestinal Disorders_stomach”, “Lung
  Adenocarcinoma_liver”, “Cancer: Lung Cancer_liver”, “Respiratory
  Conditions: Restrictive Lung Diseases_respiratory system”,
  “Cancer_integumentary system (skin)”, “Infectious and Immune-related
  Diseases: Autoimmune and Immune-Related Disorders_large intestine”,
  “Other Diseases_nasal, oral, and pharyngeal regions”, “Cancer: Lung
  Cancer_lymphatic system”, “Neurodegenerative and Neurological
  Disorders: Neurodegenerative Disorders_cerebral lobes and cortical
  areas”, “Cancer: Breast Cancer_respiratory system”, “Metabolic and
  Other Disorders: Metabolic Disorders_blood”, “Cancer: Gastrointestinal
  Cancer_liver”, “Cancer: Renal Cancer_blood”, “Glioblastoma_brainstem
  and cerebellar structures”, “Lung Adenocarcinoma_lymphatic system”,
  “Cancer: Lung Cancer_integumentary system (skin)”, “Respiratory
  Conditions: Obstructive Lung Diseases_respiratory system”, “Metabolic
  and Other Disorders_female reproductive system”, “Metabolic and Other
  Disorders: Gastrointestinal Disorders_liver”, “Cancer: Lung
  Cancer_endocrine system”, “Lung Adenocarcinoma_endocrine system”,
  “COVID-19 related_trachea”, “Neurodegenerative and Neurological
  Disorders_cerebral lobes and cortical areas”, “Other Diseases_renal
  system”, “Cancer: Renal Cancer_endocrine system”, “Cancer: Renal
  Cancer_lymphatic system”

- **Ethnicity Groups:** “African”, “Other/Unknown”, “European”, “East
  Asian”, “Hispanic/Latin American”, “South Asian”, “Native American &
  Pacific Islander”

- **Tissue Groups:** “respiratory system”, “blood”, “renal system”,
  “small intestine”, “nasal, oral, and pharyngeal regions”, “cerebral
  lobes and cortical areas”, “bone marrow”, “breast”, “cardiovascular
  system”, “trachea”, “endocrine system”, “liver”, “large intestine”,
  “lymphatic system”, “prostate”, “sensory-related structures”,
  “oesophagus”, “spleen”, “brainstem and cerebellar structures”,
  “stomach”, “female reproductive system”, “thymus”, “adipose tissue”,
  “integumentary system (skin)”, “digestive system (general)”,
  “gallbladder”, “pancreas”

# Session Info

``` r
sessionInfo()
#> R version 4.4.3 (2025-02-28)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Etc/UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] rstan_2.36.0.9000       StanHeaders_2.36.0.9000 lubridate_1.9.4        
#>  [4] forcats_1.0.1           stringr_1.6.0           dplyr_1.1.4            
#>  [7] purrr_1.2.0             readr_2.1.6             tidyr_1.3.1            
#> [10] tibble_3.3.0            ggplot2_4.0.1           tidyverse_2.0.0        
#> [13] posteriorHCA_0.1.0     
#> 
#> loaded via a namespace (and not attached):
#>   [1] gridExtra_2.3               inline_0.3.21              
#>   [3] rlang_1.1.6                 magrittr_2.0.4             
#>   [5] matrixStats_1.5.0           ggridges_0.5.7             
#>   [7] compiler_4.4.3              loo_2.8.0                  
#>   [9] callr_3.7.6                 vctrs_0.6.5                
#>  [11] pkgconfig_2.0.3             crayon_1.5.3               
#>  [13] fastmap_1.2.0               backports_1.5.0            
#>  [15] XVector_0.46.0              labeling_0.4.3             
#>  [17] utf8_1.2.6                  cmdstanr_0.9.0             
#>  [19] rmarkdown_2.30              tzdb_0.5.0                 
#>  [21] UCSC.utils_1.2.0            ps_1.9.1                   
#>  [23] bit_4.6.0                   xfun_0.53                  
#>  [25] zlibbioc_1.52.0             cachem_1.1.0               
#>  [27] GenomeInfoDb_1.42.3         jsonlite_2.0.0             
#>  [29] DelayedArray_0.32.0         parallel_4.4.3             
#>  [31] R6_2.6.1                    bslib_0.9.0                
#>  [33] stringi_1.8.7               RColorBrewer_1.1-3         
#>  [35] GenomicRanges_1.58.0        jquerylib_0.1.4            
#>  [37] Rcpp_1.1.0                  SummarizedExperiment_1.36.0
#>  [39] knitr_1.50                  IRanges_2.40.1             
#>  [41] bayesplot_1.14.0            Matrix_1.7-2               
#>  [43] timechange_0.3.0            tidyselect_1.2.1           
#>  [45] rstudioapi_0.17.1           dichromat_2.0-0.1          
#>  [47] abind_1.4-8                 yaml_2.3.10                
#>  [49] stringfish_0.17.0           codetools_0.2-20           
#>  [51] curl_7.0.0                  processx_3.8.6             
#>  [53] pkgbuild_1.4.8              lattice_0.22-6             
#>  [55] Biobase_2.66.0              withr_3.0.2                
#>  [57] bridgesampling_1.1-2        sccomp_2.1.20              
#>  [59] S7_0.2.1                    posterior_1.6.1            
#>  [61] coda_0.19-4.1               evaluate_1.0.5             
#>  [63] RcppParallel_5.1.11-1       pillar_1.11.1              
#>  [65] MatrixGenerics_1.18.1       tensorA_0.36.2.1           
#>  [67] checkmate_2.3.3             stats4_4.4.3               
#>  [69] distributional_0.5.0        generics_0.1.4             
#>  [71] vroom_1.6.6                 rprojroot_2.1.0            
#>  [73] S4Vectors_0.44.0            hms_1.1.4                  
#>  [75] rstantools_2.5.0            scales_1.4.0               
#>  [77] qs2_0.1.5                   glue_1.8.0                 
#>  [79] pheatmap_1.0.13             tools_4.4.3                
#>  [81] data.table_1.17.8           dittoSeq_1.18.0            
#>  [83] mvtnorm_1.3-3               fs_1.6.6                   
#>  [85] cowplot_1.2.0               grid_4.4.3                 
#>  [87] QuickJSR_1.8.1              instantiate_0.2.3          
#>  [89] colorspace_2.1-1            SingleCellExperiment_1.28.1
#>  [91] nlme_3.1-167                GenomeInfoDbData_1.2.13    
#>  [93] patchwork_1.3.2             cli_3.6.5                  
#>  [95] S4Arrays_1.6.0              Brobdingnag_1.2-9          
#>  [97] V8_6.0.5                    gtable_0.3.6               
#>  [99] sass_0.4.10                 digest_0.6.38              
#> [101] BiocGenerics_0.52.0         SparseArray_1.6.2          
#> [103] ggrepel_0.9.6               brms_2.23.0                
#> [105] farver_2.1.2                htmltools_0.5.8.1          
#> [107] lifecycle_1.0.4             httr_1.4.7                 
#> [109] bit64_4.6.0-1
```
