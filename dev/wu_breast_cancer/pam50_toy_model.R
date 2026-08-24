library(tibble)
library(tidyr)
library(dplyr)
library(tidySummarizedExperiment)

se_breast <- 
  readRDS("/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/Age_Clock/TAR_missing_tissue/se_list_missing.rds") %>% 
  .[['breast']] 

se_breast %>% rowData() %>% as.data.frame() %>% 
  # filter(.abundant_epithelial) %>% 
  rownames_to_column('.feature') %>% 
  dplyr::select(.feature) %>% 
  mutate(
    symb = AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = .feature,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  ) -> ab

pam50_modules_tbl <- tibble(
  module = c(
    "luminal",
    "her2",
    "basal",
    "proliferation",
    "other_pam50"
  ),
  gene = list(
    c(
      "ESR1", "PGR", "FOXA1", "BCL2", "BAG1",
      "NAT1", "MLPH", "SLC39A6", "GPR160", "MAPT"
    ),
    c(
      "ERBB2", "GRB7", "FGFR4"
    ),
    c(
      "KRT5", "KRT14", "KRT17", "EGFR",
      "FOXC1", "CDH3", "MIA", "PHGDH"
    ),
    c(
      "ANLN", "BIRC5", "CCNB1", "CCNE1", "CDC20", "CDC6",
      "CDCA1", "CENPF", "CEP55", "EXO1", "KIF2C", "KNTC2",
      "MELK", "MKI67", "MYBL2", "ORC6L", "PTTG1", "RRM2",
      "TYMS", "UBE2C", "UBE2T"
    ),
    c(
      "ACTR3B", "BLVRA", "CXXC5", "MDM2", "MMP11",
      "MYC", "SFRP1", "TMEM45B"
    )
  )
) %>%
  unnest(gene)

pam50_modules_tbl <-
  pam50_modules_tbl %>% 
  mutate(
    .feature = AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = gene,
      column = "ENSEMBL",
      keytype = "SYMBOL",
      multiVals = "first"
    )
  )

(pam50_modules_tbl$gene %in% ab$symb ) %>% sum()

ab %>% 
  filter(symb %in% pam50_modules_tbl$gene) %>% 
  pull(.feature) -> biomarkers
  
library(targets)
library(dplyr)
library(purrr)
library(furrr)
library(future)
library(lobstr)
library(posterior)
library(progressr)
library(brms)
library(tidyverse)

handlers(global = TRUE)
handlers("progress")

options(future.globals.maxSize = 5 * 1024^3)  # 5 GiB

# plan(multisession, workers = 32)

n_components = 5
corr_threshold = 0.1

with_progress({
  p <- progressor(along = biomarkers)
  
  pam50_estimates_chunk <-
    purrr::map_dfr(
      biomarkers,
      function(ensg) {
        
        p(sprintf("Processing %s", ensg))
        
        cmdstanr::set_cmdstan_path(
          "/hpcfs/users/a1237163/R/.cmdstan/cmdstan-2.36.0"
        )
        
        estimates_chunk <-
          tibble(
            .feature = ensg,
            symb = ab$symb[ab$.feature == ensg],
            cell_type = "epithelial"
          ) %>% 
          mutate(
            se = ensg %>% map(
              ~{
                se_breast %>%
                  filter(cell_type_unified_ensemble == "epithelial") %>%
                  .[.x, , drop = FALSE]
              }
            )
          ) %>% 
          mutate(
            brms_fit = map(se, ~ {
              
              data =
                .x |>
                as_tibble() |>
                mutate(counts = counts |> as.integer()) |>
                droplevels()
              
              n_NAs = data |> filter(counts |> is.na()) |> nrow()
              
              if (n_NAs > 0) {
                warning(glue(
                  "You have {n_NAs} NAs in counts. They have been filtered out"
                ))
                
                data =
                  data |>
                  filter(!counts |> is.na()) |>
                  droplevels()
              }
              
              colnames(data) =
                colnames(data) |>
                stringr::str_replace_all("_+", "_")
              
              formula_chr_counts <-
                "counts ~ 1 + offset(offset) + age_decade * sex + disease_groups_altered + ethnicity_groups + assay_groups_altered + (1 | dataset_id_altered)"
              
              formula_chr_shape <-
                "shape ~ 1 + disease_groups_altered + assay_groups_altered + ethnicity_groups "
              
              if (data |> distinct(disease_groups_altered) |> nrow() == 1) {
                formula_chr_counts =
                  formula_chr_counts |>
                  str_remove_all(fixed("+ disease_groups_altered"))
                
                formula_chr_shape =
                  formula_chr_shape |>
                  str_remove_all(fixed("+ disease_groups_altered"))
              }
              
              if (data |> distinct(ethnicity_groups) |> nrow() == 1) {
                formula_chr_counts =
                  formula_chr_counts |>
                  str_remove_all(fixed("+ ethnicity_groups"))
                
                formula_chr_shape =
                  formula_chr_shape |>
                  str_remove_all(fixed("+ ethnicity_groups"))
              }
              
              if (data |> distinct(assay_groups_altered) |> nrow() == 1) {
                formula_chr_counts =
                  formula_chr_counts |>
                  str_remove_all(fixed("+ assay_groups_altered"))
                
                formula_chr_shape =
                  formula_chr_shape |>
                  str_remove_all(fixed("+ assay_groups_altered"))
              }
              
              if (data |> distinct(dataset_id_altered) |> nrow() == 1) {
                formula_chr_counts =
                  formula_chr_counts |>
                  str_remove_all(fixed("+ (1 | dataset_id_altered)"))
              }
              
              if (data |> distinct(sex) |> nrow() == 1) {
                formula_chr_counts =
                  formula_chr_counts |>
                  str_remove_all(fixed("* sex"))
              }
              
              formula <- bf(
                as.formula(formula_chr_counts),
                as.formula(formula_chr_shape)
              )
              
              if (stringr::str_trim(formula_chr_shape) != "shape ~ 1") {
                
                prior = c(
                  prior(student_t(3, i, 1.5), class = Intercept),
                  prior(student_t(3, 0, 1), class = Intercept, dpar = shape),
                  prior(student_t(3, 0, 5), class = b),
                  prior(student_t(3, 0, 2), class = b, dpar = shape)
                ) |>
                  substitute(env = list(
                    i = mean(log1p(data$counts / exp(data$offset)))
                  )) |>
                  eval()
                
              } else {
                
                prior = c(
                  prior(student_t(3, i, 1.5), class = Intercept),
                  prior(student_t(3, 0, 1), class = Intercept, dpar = shape),
                  prior(student_t(3, 0, 5), class = b)
                ) |>
                  substitute(env = list(
                    i = mean(log1p(data$counts / exp(data$offset)))
                  )) |>
                  eval()
              }
              
              chains <- 2
              
              bterms <- brmsterms(
                formula = brms:::validate_formula(
                  formula,
                  data = data,
                  family = zero_inflated_negbinomial(),
                  autocor = NULL,
                  sparse = NULL,
                  cov_ranef = NULL
                )
              )
              
              bframe <- brms:::brmsframe(bterms, data)
              
              sdata <- brms:::.standata(
                bframe,
                data = data,
                prior = prior,
                data2 = NULL,
                stanvars = NULL,
                threads = NULL
              )
              
              Kc <- sdata$Kc
              Kc_shape <- sdata$Kc_shape
              
              inits <- lapply(seq_len(chains), function(i) {
                list(
                  b = rnorm(Kc, 0, 5),
                  Intercept = rnorm(
                    1,
                    mean(log1p(data$counts / exp(data$offset))),
                    1.5
                  ),
                  b_shape = rnorm(Kc_shape, 0, 2),
                  Intercept_shape = rnorm(1, 0, 1)
                )
              })
              
              res <- brm(
                formula = formula,
                data = data,
                family = zero_inflated_negbinomial(),
                prior = prior,
                chains = chains,
                cores = 2,
                warmup = 400,
                refresh = 10,
                backend = "cmdstanr",
                init = inits,
                iter = 600
              ) %>% 
                poco::strip_brmsfit_envs()
              
              if (nrow(res$data) != nrow(data)) {
                warning(glue(
                  "The number of rows in the data and the number of rows in the brms object are different. The data has been filtered out."
                ))
              }
              
              return(res)
            })
          )
        
        fit_brms <- estimates_chunk$brms_fit[[1]]
        
        posterior_cor <- 
          poco::posterior_correlation(
            fit_brms,
            type = "correlation",
            method = "pearson"
          )
        
        original_draws <- 
          posterior::as_draws_matrix(fit_brms)
        
        partition <- poco::partition_parameters_clusters(
          posterior_cor,
          threshold = corr_threshold,
          min_size = 2L,
          simple_output = FALSE
        )
        
        layout_tbl <- tibble::tibble(
          n_components = n_components,
          corr_threshold = corr_threshold,
          n_blocks = length(partition$blocks),
          block_size_max = if (length(partition$blocks)) {
            max(lengths(partition$blocks))
          } else {
            0L
          },
          block_size_med = if (length(partition$blocks)) {
            as.integer(stats::median(lengths(partition$blocks)))
          } else {
            0L
          },
          n_remainder = length(partition$remainder)
        )
        
        compress_tic <- Sys.time()
        
        compressed_blockwise <- tryCatch(
          poco::compress_brmsfit(
            brmsfit = fit_brms,
            method = "mclust",
            n_components = n_components,
            partition = partition,
            cluster_BPPARAM = 1L,
            verbose = FALSE
          ),
          error = function(e) {
            warning(
              "compress_brmsfit failed (n_components = ", n_components,
              ", corr_threshold = ", corr_threshold, "): ",
              conditionMessage(e),
              call. = FALSE,
              immediate. = TRUE
            )
            NULL
          }
        )
        
        time_used <- as.numeric(
          difftime(Sys.time(), compress_tic, units = "secs")
        )
        
        fid_blockwise <- tryCatch(
          poco::evaluate_compression(
            compressed_blockwise$compressed,
            original_draws,
            seed = 1L
          ),
          error = function(e) {
            warning(
              "evaluate_compression failed (n_components = ", n_components,
              ", corr_threshold = ", corr_threshold, "): ",
              conditionMessage(e),
              call. = FALSE,
              immediate. = TRUE
            )
            NULL
          }
        )
        
        if (is.null(compressed_blockwise)) {
          
          layout_tbl <- tibble::add_column(
            layout_tbl,
            compression_ok = FALSE,
            time_used = NA_real_,
            compression_rate = NA_real_,
            energy_pct = NA_real_,
            c2st_pct = NA_real_,
            reproduction_pct = NA_real_
          )
          
        } else {
          
          estimates_chunk$brms_fit_compressed <- list(compressed_blockwise)
          
          layout_tbl <- tibble::add_column(
            layout_tbl,
            compression_ok = TRUE,
            time_used = time_used,
            compression_rate =
              as.numeric(lobstr::obj_size(compressed_blockwise)) /
              as.numeric(lobstr::obj_size(fit_brms))
          )
        }
        
        if (is.null(fid_blockwise)) {
          
          layout_tbl <- tibble::add_column(
            layout_tbl,
            energy_pct = NA_real_,
            c2st_pct = NA_real_,
            reproduction_pct = NA_real_
          )
          
        } else {
          
          en <- fid_blockwise$metrics$energy$reproduction_pct
          c2 <- fid_blockwise$metrics$c2st$reproduction_pct
          
          layout_tbl <- tibble::add_column(
            layout_tbl,
            energy_pct = if (is.null(en)) NA_real_ else en,
            c2st_pct = if (is.null(c2)) NA_real_ else c2,
            reproduction_pct = fid_blockwise$reproduction_pct
          )
        }
        
        estimates_chunk %>%  
          tibble::add_column(layout_tbl)
      }
    )
})

pam50_estimates_chunk %>% saveRDS('dev/wu_breast_cancer/data/pam50_estimates_chunk.rds')


pam50_estimates_chunk <- readRDS('dev/wu_breast_cancer/data/pam50_estimates_chunk.rds')

se_cancer_epithelial <- readRDS("/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/wu_breast_cancer/data/Wu_etal_2021_BRCA_se_cancer_epithelial.rds")
se_breast_epithelial_ref_sample <-
  se_breast %>%
  filter(cell_type_unified_ensemble == "epithelial") %>%
  filter(offset == 0)

# Harmonise Wu (gene symbols) to ENSEMBL; rebuild SE from counts + colData.
se_cancer_epithelial_ensembl <- {
  sym <- rownames(se_cancer_epithelial)
  ensg <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = sym,
    column = "ENSEMBL",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  keep <- !is.na(ensg) & !duplicated(ensg)
  ensg <- ensg[keep]

  counts_mat <- SummarizedExperiment::assay(se_cancer_epithelial, "counts")[keep, , drop = FALSE]
  rownames(counts_mat) <- ensg
  colnames(counts_mat) <- colnames(se_cancer_epithelial)

  cd <- as.data.frame(SummarizedExperiment::colData(se_cancer_epithelial))
  cd$sample_role <- "cancer_pseudobulk"
  cd$cell_type_unified_ensemble <- "epithelial"
  cd$offset <- 0

  # rownames are ENSEMBL IDs; omit `.feature` in rowData (tidySummarizedExperiment
  # injects `.feature` from rownames and duplicate columns break construction).
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts_mat),
    colData = S4Vectors::DataFrame(cd, row.names = colnames(counts_mat)),
    rowData = S4Vectors::DataFrame(
      gene_symbol = sym[keep],
      row.names = ensg
    )
  )
}

colData(se_breast_epithelial_ref_sample)$sample_role <- "reference"

merge_summarized_experiments <- function(...) {
  ses <- list(...)
  assay_names <- Reduce(intersect, lapply(ses, SummarizedExperiment::assayNames))
  if (length(assay_names) == 0) {
    stop("No shared assays between SummarizedExperiment objects.")
  }

  genes <- Reduce(intersect, lapply(ses, rownames))
  if (length(genes) == 0) {
    stop("No shared rownames between SummarizedExperiment objects.")
  }

  ses <- lapply(ses, function(se) se[genes, , drop = FALSE])

  cds <- lapply(ses, function(se) as.data.frame(SummarizedExperiment::colData(se)))
  all_nms <- Reduce(union, lapply(cds, names))
  cds <- lapply(cds, function(cd) {
    missing <- setdiff(all_nms, names(cd))
    for (nm in missing) {
      cd[[nm]] <- NA
    }
    cd[, all_nms, drop = FALSE]
  })

  assays <- stats::setNames(
    lapply(assay_names, function(an) {
      do.call(
        cbind,
        lapply(ses, function(se) SummarizedExperiment::assay(se, an))
      )
    }),
    assay_names
  )

  SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    colData = S4Vectors::DataFrame(
      do.call(rbind, cds),
      row.names = unlist(lapply(ses, colnames), use.names = FALSE)
    ),
    rowData = SummarizedExperiment::rowData(ses[[1]])
  )
}

se_epithelial <- merge_summarized_experiments(
  se_cancer_epithelial_ensembl,
  se_breast_epithelial_ref_sample
)


rowData(se_epithelial)$.abundant =  rowData(se_breast_epithelial_ref_sample)[ rowData(se_epithelial) %>% rownames(), ".abundant_epithelial"]

se_epithelial <-
  se_epithelial %>% 
  tidybulk::scale_abundance(method = 'TMMwsp', reference_sample = se_epithelial %>% colData() %>% as.data.frame() %>% filter(sample_role == "reference") %>% rownames()) %>% 
  mutate(offset = log(1 / multiplier)) 

i = 27
j = 1

newdata <-
  tibble(
    age_decade = NA,
    # sex = 'female',
    # disease_groups_altered = 'Normal',
    ethnicity_groups = NA,
    assay_groups_altered = "10x Genomics 5",
    # dataset_id_altered = NA,
    offset = se_epithelial$offset[j]
  )

#' Compare an observed count to the healthy posterior predictive distribution.
#'
#' Uses `posterior_predict()` (full ZINB sampling distribution), not `posterior_epred()`.
#' @param fit A `brmsfit` object (or reconstructed from compressed fit).
#' @param newdata Covariate profile for the tumour sample.
#' @param obs Observed count.
#' @param cred_mass Mass for the HPD interval (default 0.95).
#' @param alpha One-sided tail threshold (default 0.025).
gene_deviation_from_fit <- function(
  fit,
  newdata,
  obs,
  cred_mass = 0.95,
  alpha = 0.025
) {
  y_rep <- brms::posterior_predict(
    fit,
    newdata = newdata,
    summary = FALSE,
    re_formula = NA
  ) |>
    as.numeric()

  hpd <- if (requireNamespace("HDInterval", quietly = TRUE)) {
    HDInterval::hdi(y_rep, credMass = cred_mass)
  } else {
    probs <- c((1 - cred_mass) / 2, 1 - (1 - cred_mass) / 2)
    stats::quantile(y_rep, probs = probs, names = FALSE)
  }

  p_above <- mean(y_rep >= obs)
  p_below <- mean(y_rep <= obs)

  tibble(
    observed_count = obs,
    healthy_median = median(y_rep),
    healthy_mean = mean(y_rep),
    healthy_hpd_lo = hpd[1],
    healthy_hpd_hi = hpd[2],
    obs_percentile = mean(y_rep <= obs),
    p_above_healthy = p_above,
    p_below_healthy = p_below,
    posterior_tail_p = 2 * min(p_above, p_below),
    log2fc_vs_median = log2((obs + 1) / (median(y_rep) + 1)),
    direction = dplyr::case_when(
      p_above < alpha ~ "higher_than_healthy",
      p_below < alpha ~ "lower_than_healthy",
      TRUE ~ "consistent_with_healthy"
    ),
    y_rep = list(y_rep)
  )
}

fit_compressed <-
  pam50_estimates_chunk$brms_fit_compressed[[i]] %>%
  poco::reconstruct_brmsfit()

obs <-
  se_epithelial[pam50_estimates_chunk$.feature[i], j] %>%
  assay("counts") %>%
  as.numeric()

gene_deviation_res <- gene_deviation_from_fit(fit_compressed, newdata, obs)
y_rep <- gene_deviation_res$y_rep[[1]]

gene_deviation <-
  gene_deviation_res %>%
  dplyr::select(-y_rep) %>%
  mutate(
    .feature = pam50_estimates_chunk$.feature[i],
    symb = pam50_estimates_chunk$symb[i],
    sample = colnames(se_epithelial)[j],
    .before = observed_count
  )

gene_deviation

se_epithelial[pam50_estimates_chunk$.feature[i], j]  
