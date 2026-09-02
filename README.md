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

To rebuild this vignette **and** sync `README.md`, run
`source("scripts/render_introduction.R"); render_introduction()` (or
`Rscript scripts/render_introduction_cli.R` from the package root).
Plain `rmarkdown::render("vignettes/Introduction.Rmd")` only builds the
vignette HTML.

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

#### Load the healthy sccomp reference

The default composition model is a healthy-only CellNexus sccomp fit.
Unknown covariates can be left as `NA`; sccomp marginalises over them
automatically.

``` r
sccomp_fit <- load_sccomp_fit()
#> ℹ Using cached file: /home/a1237163/.cache/R/posteriorHCA/sccomp_est/cellNexus_1_0_12/estimates_age_decade___L3___disease_FALSE___immune_only_TRUE.rds
```

#### Draw healthy reference proportions

``` r
comp_draws <- composition_draws(
  sccomp_fit,
  sex = "male",
  age_decade = "4",
  ethnicity_groups = "European",
  assay_groups = "10x Genomics 3",
  tissue_groups = "blood"
)
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 4.318 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
head(comp_draws$draws)
#>      sample_id age_decade  sex ethnicity_groups_imputed assay_groups___altered
#> 1 query_sample          4 male                 European         10x Genomics 3
#> 2 query_sample          4 male                 European         10x Genomics 3
#> 3 query_sample          4 male                 European         10x Genomics 3
#> 4 query_sample          4 male                 European         10x Genomics 3
#> 5 query_sample          4 male                 European         10x Genomics 3
#> 6 query_sample          4 male                 European         10x Genomics 3
#>   dataset_id___altered tissue_groups       L3 proportion .draw cell_type
#> 1                   NA         blood b memory 0.05109034     1  b memory
#> 2                   NA         blood b memory 0.04868110     2  b memory
#> 3                   NA         blood b memory 0.04549951     3  b memory
#> 4                   NA         blood b memory 0.04742480     4  b memory
#> 5                   NA         blood b memory 0.04457948     5  b memory
#> 6                   NA         blood b memory 0.04695580     6  b memory
```

#### Compare observed proportions

`composition_posterior_test()` keeps the legacy table + plot interface,
but now uses \[load_sccomp_fit()\] and \[composition_draws()\]
internally.

``` r
result <- composition_posterior_test(
    proportions = example_proportions,
    sex = "male",
    age_decade = "4",
    ethnicity_groups = "European",
    assay_groups = "10x Genomics 3",
    tissue_groups = "blood",
    fit = sccomp_fit
)
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 3.447 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> # A tibble: 26 × 7
#>    sample_id_observed cell_type proportion_observed Empirical_Confidence    mean
#>    <chr>              <chr>                   <dbl>                <dbl>   <dbl>
#>  1 0000c153da22cf963… b memory             0                           0 4.64e-2
#>  2 0000c153da22cf963… b naive              0                           0 5.10e-2
#>  3 0000c153da22cf963… cd14 mono            0                           0 3.73e-2
#>  4 0000c153da22cf963… cd16 mono            0.00502                     0 1.22e-2
#>  5 0000c153da22cf963… cd4 fh em            0                           0 6.82e-4
#>  6 0000c153da22cf963… cd4 naive            0                           0 2.12e-1
#>  7 0000c153da22cf963… cd4 tcm              0.000359                    0 1.85e-1
#>  8 0000c153da22cf963… cd4 th1 …            0                           0 3.74e-3
#>  9 0000c153da22cf963… cd4 th1/…            0.000359                    0 6.91e-3
#> 10 0000c153da22cf963… cd4 th17…            0                           0 1.73e-3
#> # ℹ 16 more rows
#> # ℹ 2 more variables: lower <dbl>, upper <dbl>
```

![](/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

#### Inspect Composition Results

**Posterior Distribution Table:**

``` r
print(result$result_table)
#> # A tibble: 26 × 7
#>    sample_id_observed cell_type proportion_observed Empirical_Confidence    mean
#>    <chr>              <chr>                   <dbl>                <dbl>   <dbl>
#>  1 0000c153da22cf963… b memory             0                           0 4.64e-2
#>  2 0000c153da22cf963… b naive              0                           0 5.10e-2
#>  3 0000c153da22cf963… cd14 mono            0                           0 3.73e-2
#>  4 0000c153da22cf963… cd16 mono            0.00502                     0 1.22e-2
#>  5 0000c153da22cf963… cd4 fh em            0                           0 6.82e-4
#>  6 0000c153da22cf963… cd4 naive            0                           0 2.12e-1
#>  7 0000c153da22cf963… cd4 tcm              0.000359                    0 1.85e-1
#>  8 0000c153da22cf963… cd4 th1 …            0                           0 3.74e-3
#>  9 0000c153da22cf963… cd4 th1/…            0.000359                    0 6.91e-3
#> 10 0000c153da22cf963… cd4 th17…            0                           0 1.73e-3
#> # ℹ 16 more rows
#> # ℹ 2 more variables: lower <dbl>, upper <dbl>
```

**Density Plot of Posterior Predictions:**

``` r
print(result$plot)
```

![](/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

### **2. Gene Expression Prediction**

The `expr_predict()` function allows predicting gene expression levels
for specific cell types and genes based on sample metadata:

``` r
# Predict gene expression for a specific cell type and gene
expr_result <- expr_predict(
  cell_type = "cd4 naive",
  gene = "ENSG00000000419",
  age_decade = "7",
  sex = "female",
  disease_groups = "Normal",
  ethnicity_groups = "European",
  assay_groups = "10x Genomics 3",
  tissue_groups = "blood"
)
#> ℹ Using cached file: /home/a1237163/.cache/R/posteriorHCA/meta/latest.csv
#> ℹ Using cached file: /home/a1237163/.cache/R/posteriorHCA/V1/cd4.naive/genes.csv
#> ℹ Using cached file: /home/a1237163/.cache/R/posteriorHCA/meta/latest.csv
#> ℹ Using cached file: /home/a1237163/.cache/R/posteriorHCA/V1/cd4.naive/ENSG00000000419
#> ℹ Covariate grid has 1 profile (all metadata fixed).
```

#### Inspect Gene Expression Results

**Summary Statistics:**

``` r
print(expr_result$summary)
#> $mean
#> [1] 4.442823
#> 
#> $median
#> [1] 4.45129
#> 
#> $peak_location
#> [1] 4.459412
```

**Density Plot of Predicted Expression:**

``` r
print(expr_result$plot)
```

![](/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

**Predicted Values:**

``` r
head(expr_result$pred)
#>      value
#> 1 3.670794
#> 2 4.554274
#> 3 4.264375
#> 4 4.193046
#> 5 4.535664
#> 6 3.602183
```

## **Input Arguments and Valid Values**

Both `composition_posterior_test()` / \[composition_draws()\] and
\[expr_predict()\] accept metadata inputs corresponding to observed
sample characteristics. For composition, the default healthy sccomp
model uses **sex, age decade, ethnicity, assay, and tissue**. Unknown
covariates can be set to `NA` and are marginalised automatically by
sccomp. For expression queries, metadata can also be left empty to
marginalise over the pre-trained model.

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

<!-- -   **Disease Groups:** "Normal", "COVID-19 related_blood", "Metabolic and Other Disorders: Renal and Urinary Disorders_renal system", "Systemic Lupus Erythematosus_blood", "Respiratory Conditions: General Respiratory Disorders_blood", "COVID-19 related_nasal, oral, and pharyngeal regions", "COVID-19 related_respiratory system", "Glioblastoma_cerebral lobes and cortical areas", "Cancer: Hematologic Cancer_bone marrow", "Lung Adenocarcinoma_respiratory system", "Infectious and Immune-related Diseases: Autoimmune and Immune-Related Disorders_small intestine", "Infectious and Immune-related Diseases: Respiratory Infections_blood", "Metabolic and Other Disorders: Gastrointestinal Disorders_large intestine", "Cancer: Renal Cancer_renal system", "Infectious and Immune-related Diseases: Autoimmune and Immune-Related Disorders_respiratory system", "Infectious and Immune-related Diseases: Respiratory Infections_nasal, oral, and pharyngeal regions", "Metabolic and Other Disorders: Renal and Urinary Disorders_prostate", "Metabolic and Other Disorders_sensory-related structures", "Cancer: Gastrointestinal Cancer_large intestine", "Other Diseases_brainstem and cerebellar structures", "Metabolic and Other Disorders: Gastrointestinal Disorders_oesophagus", "Cancer: Lung Cancer_respiratory system", "Lung Adenocarcinoma_cerebral lobes and cortical areas", "Metabolic and Other Disorders: Gastrointestinal Disorders_stomach", "Lung Adenocarcinoma_liver", "Cancer: Lung Cancer_liver", "Respiratory Conditions: Restrictive Lung Diseases_respiratory system", "Cancer_integumentary system (skin)", "Infectious and Immune-related Diseases: Autoimmune and Immune-Related Disorders_large intestine", "Other Diseases_nasal, oral, and pharyngeal regions", "Cancer: Lung Cancer_lymphatic system", "Neurodegenerative and Neurological Disorders: Neurodegenerative Disorders_cerebral lobes and cortical areas", "Cancer: Breast Cancer_respiratory system", "Metabolic and Other Disorders: Metabolic Disorders_blood", "Cancer: Gastrointestinal Cancer_liver", "Cancer: Renal Cancer_blood", "Glioblastoma_brainstem and cerebellar structures", "Lung Adenocarcinoma_lymphatic system", "Cancer: Lung Cancer_integumentary system (skin)", "Respiratory Conditions: Obstructive Lung Diseases_respiratory system", "Metabolic and Other Disorders_female reproductive system", "Metabolic and Other Disorders: Gastrointestinal Disorders_liver", "Cancer: Lung Cancer_endocrine system", "Lung Adenocarcinoma_endocrine system", "COVID-19 related_trachea", "Neurodegenerative and Neurological Disorders_cerebral lobes and cortical areas", "Other Diseases_renal system", "Cancer: Renal Cancer_endocrine system", "Cancer: Renal Cancer_lymphatic system" -->

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
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
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
#>  [1] lubridate_1.9.5    forcats_1.0.1      stringr_1.6.0      dplyr_1.2.1       
#>  [5] purrr_1.2.2        readr_2.2.0        tidyr_1.3.2        tibble_3.3.1      
#>  [9] ggplot2_4.0.3      tidyverse_2.0.0    posteriorHCA_0.2.0 testthat_3.3.2    
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3          tensorA_0.36.2.1           
#>   [3] rstudioapi_0.18.0           jsonlite_2.0.0             
#>   [5] magrittr_2.0.5              TH.data_1.1-5              
#>   [7] estimability_1.5.1          farver_2.1.2               
#>   [9] rmarkdown_2.31              fs_2.1.0                   
#>  [11] vctrs_0.7.3                 memoise_2.0.1              
#>  [13] htmltools_0.5.9             S4Arrays_1.10.1            
#>  [15] usethis_3.2.1               curl_7.1.0                 
#>  [17] distributional_0.7.0        SparseArray_1.10.10        
#>  [19] StanHeaders_2.39.0.9000     sass_0.4.10                
#>  [21] bslib_0.10.0                desc_1.4.3                 
#>  [23] sandwich_3.1-1              emmeans_2.0.3              
#>  [25] zoo_1.8-15                  cachem_1.1.0               
#>  [27] lifecycle_1.0.5             pkgconfig_2.0.3            
#>  [29] Matrix_1.7-5                R6_2.6.1                   
#>  [31] fastmap_1.2.0               MatrixGenerics_1.22.0      
#>  [33] digest_0.6.39               colorspace_2.1-2           
#>  [35] patchwork_1.3.2             S4Vectors_0.48.1           
#>  [37] ps_1.9.3                    rprojroot_2.1.1            
#>  [39] brms_2.23.0                 pkgload_1.5.2              
#>  [41] GenomicRanges_1.62.1        labeling_0.4.3             
#>  [43] timechange_0.4.0            httr_1.4.8                 
#>  [45] abind_1.4-8                 compiler_4.5.3             
#>  [47] bit64_4.8.0                 withr_3.0.2                
#>  [49] inline_0.3.21               S7_0.2.2                   
#>  [51] backports_1.5.1             QuickJSR_1.9.2             
#>  [53] pkgbuild_1.4.8              MASS_7.3-65                
#>  [55] DelayedArray_0.36.1         sessioninfo_1.2.3          
#>  [57] loo_2.10.1.9000             tools_4.5.3                
#>  [59] otel_0.2.0                  glue_1.8.1                 
#>  [61] callr_3.7.6                 nlme_3.1-169               
#>  [63] cmdstanr_0.9.0              grid_4.5.3                 
#>  [65] instantiate_0.2.3           checkmate_2.3.4            
#>  [67] generics_0.1.4              gtable_0.3.6               
#>  [69] tzdb_0.5.0                  data.table_1.18.2.1        
#>  [71] hms_1.1.4                   utf8_1.2.6                 
#>  [73] stringfish_0.19.0           XVector_0.50.0             
#>  [75] BiocGenerics_0.56.0         ggrepel_0.9.8              
#>  [77] pillar_1.11.1               vroom_1.7.1                
#>  [79] limma_3.66.0                posterior_1.7.1            
#>  [81] splines_4.5.3               lattice_0.22-9             
#>  [83] bit_4.6.0                   survival_3.8-6             
#>  [85] tidyselect_1.2.1            SingleCellExperiment_1.32.0
#>  [87] locfit_1.5-9.12             knitr_1.51                 
#>  [89] gridExtra_2.3               IRanges_2.44.0             
#>  [91] Seqinfo_1.0.0               edgeR_4.8.2                
#>  [93] SummarizedExperiment_1.40.0 stats4_4.5.3               
#>  [95] xfun_0.57                   dittoSeq_1.22.0            
#>  [97] bridgesampling_1.2-1        Biobase_2.70.0             
#>  [99] statmod_1.5.1               devtools_2.5.1             
#> [101] brio_1.1.5                  matrixStats_1.5.0          
#> [103] rstan_2.39.0.9000           pheatmap_1.0.13            
#> [105] stringi_1.8.7               sccomp_2.2.0               
#> [107] yaml_2.3.12                 evaluate_1.0.5             
#> [109] codetools_0.2-20            cli_3.6.6                  
#> [111] RcppParallel_5.1.11-2       xtable_1.8-8               
#> [113] processx_3.9.0              jquerylib_0.1.4            
#> [115] dichromat_2.0-0.1           Rcpp_1.1.1-1               
#> [117] coda_0.19-4.1               parallel_4.5.3             
#> [119] rstantools_2.7.0            ellipsis_0.3.3             
#> [121] bayesplot_1.16.0.9000       Brobdingnag_1.2-9          
#> [123] mvtnorm_1.3-7               scales_1.4.0               
#> [125] ggridges_0.5.7              crayon_1.5.3               
#> [127] rlang_1.2.0                 cowplot_1.2.0              
#> [129] qs2_0.2.0                   multcomp_1.4-30
```
