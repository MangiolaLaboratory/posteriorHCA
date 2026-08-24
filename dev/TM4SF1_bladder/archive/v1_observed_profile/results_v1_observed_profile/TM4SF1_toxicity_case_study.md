---
title: "TM4SF1 in histological-variant bladder cancer"
subtitle: "Replicated tumour-expression signal and posteriorHCA overlap screen"
author: "posteriorHCA case study"
date: last-modified
format:
  html:
    toc: true
    toc-depth: 3
    code-fold: true
    code-summary: "Show code"
    df-print: paged
    embed-resources: false
execute:
  echo: true
  warning: false
  message: false
  freeze: auto
params:
  case_dir: "."
  model_root: "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE"
  gene_symbol: "TM4SF1"
  gene_ensembl: "ENSG00000169908"
  model_cell_types: ["epithelial", "endothelial"]
  n_boot: 2000
  threshold_quantile: 0.05
  alert_probability: 0.80
  min_donor_keys: 10
  refresh_hca: false
  seed: 293189
---

## Question and interpretation boundary

The bladder-cancer paper identified TM4SF1 as a candidate antigen in histological-variant (HV) tumours and noted that expression in normal human tissues remains a potential safety concern. This report first reconstructs that tumour-expression signal, then asks a narrower safety-screening question: **which Normal/10x3-supported posteriorHCA tissue, cell-type and demographic profiles have modelled TM4SF1 expression that overlaps the HV tumour-epithelial range?**

An overlap flag is not evidence of toxicity. It is a reproducible prioritisation rule for follow-up protein localisation, antigen-density, accessibility and functional-killing experiments. Paired/adjacent normal bladder samples are retained as calibration samples and are never called fully healthy tissue. Three have at least 20 epithelial cells (UC04, VAR07 and VAR10), but only VAR07 has a counterpart in the retained tumour cohort, so this is not a matched-pair analysis.



## Analysis map

The report keeps the biological unit as the patient throughout. Cell-level source data are used only to reproduce the paper figure; all cohort inference uses patient-level summed raw counts.

1. Re-read the authors' released Figure 5 source tables and confirm the TM4SF1 differential-expression rank.
2. Display both the released per-cell contrast and the reconstructed patient-level raw-count contrast.
3. Align the study pseudobulks to each posteriorHCA model's own offset-zero TMMwsp reference.
4. Estimate patient-cohort latent `log(mu)` with the SAVI study's negative-binomial/Bayesian-bootstrap method.
5. Predict posteriorHCA `log(mu)` for observed Normal, 10x Genomics 3 donor–tissue profiles.
6. Define a model-specific HV tumour-expression anchor and identify supported profiles whose posterior overlaps it.

## 1. Reproduce the paper's TM4SF1 highlight


``` r
paper_workbook <- file.path(
  case_dir, "data", "raw", "paper", "41467_2025_59888_MOESM4_ESM.xlsx"
)

paper_de <- readxl::read_excel(paper_workbook, sheet = "Figure5A") |>
  filter(cluster == "Variant") |>
  mutate(
    gene = as.character(gene),
    hv_enriched = avg_log2FC > 0,
    hv_rank = rank(
      if_else(hv_enriched, p_val, NA_real_),
      ties.method = "min",
      na.last = "keep"
    ),
    minus_log10_fdr = -log10(pmax(p_val_adj, .Machine$double.xmin)),
    is_tm4sf1 = gene == tm4_symbol
  )

tm4_paper_de <- paper_de |>
  filter(is_tm4sf1) |>
  transmute(
    source = "paper Figure5A",
    gene,
    avg_log2FC,
    pct_HV = pct.1,
    pct_UC = pct.2,
    p_value = p_val,
    adjusted_p_value = p_val_adj,
    HV_enriched_rank = as.integer(hv_rank)
  )

replication_de <- readr::read_tsv(
  file.path(
    case_dir, "results", "replication_audit", "TM4SF1_original_DE_result.tsv"
  ),
  show_col_types = FALSE
) |>
  transmute(
    analysis,
    n_HV_patients,
    n_UC_patients,
    cells_per_patient = n_cells_per_patient,
    avg_log2FC,
    p_value = p_val,
    adjusted_p_value = coalesce(FDR_BH, p_val_adj),
    HV_enriched_rank = HV_enriched_rank_by_p
  )

write_tsv(
  replication_de,
  file.path(results_dir, "TM4SF1_paper_DE_replication.tsv")
)
```

The released table places TM4SF1 at rank 1 among HV-enriched genes: log2 fold change 2.36, detection 37.6% in HV versus 2.9% in pure UC, and adjusted (p = 4.45e-40). The closest public-DGE reconstruction also ranks TM4SF1 first, but it contains 9 HV and only 2 UC patients at the paper's 150-cell minimum because the deposited UC01 matrix does not contain the published UC01 barcodes.


``` r
replication_de |>
  mutate(
    avg_log2FC = round(avg_log2FC, 2),
    p_value = format(p_value, scientific = TRUE, digits = 3),
    adjusted_p_value = format(adjusted_p_value, scientific = TRUE, digits = 3)
  ) |>
  kable(
    caption = "TM4SF1 in the released paper result and public-DGE sensitivity analyses."
  )
```



Table: TM4SF1 in the released paper result and public-DGE sensitivity analyses.

|analysis                                                            | n_HV_patients| n_UC_patients| cells_per_patient| avg_log2FC|p_value  |adjusted_p_value | HV_enriched_rank|
|:-------------------------------------------------------------------|-------------:|-------------:|-----------------:|----------:|:--------|:----------------|----------------:|
|paper_released_source_DE_9HV_3UC_150_cells_per_patient_unknown_seed |             9|             3|               150|       2.36|1.68e-44 |4.45e-40         |                1|
|closest_public_DGE_9HV_2UC_150_cells_per_patient_seed293189         |             9|             2|               150|       4.30|2.06e-29 |2.38e-26         |                1|
|secondary_public_DGE_9HV_3UC_90_cells_per_patient_seed293189        |             9|             3|                90|       2.42|1.82e-11 |2.26e-09         |               28|


``` r
volcano_cap <- quantile(paper_de$minus_log10_fdr, 0.995, na.rm = TRUE)
paper_de_display <- paper_de |>
  mutate(plot_y = pmin(minus_log10_fdr, volcano_cap))
paper_de_plot <- paper_de_display |>
  ggplot(aes(avg_log2FC, plot_y)) +
  geom_point(colour = "grey72", alpha = 0.45, size = 1) +
  geom_point(
    data = filter(paper_de_display, is_tm4sf1),
    colour = "#B2182B", size = 3
  ) +
  geom_text(
    data = filter(paper_de_display, is_tm4sf1),
    aes(label = gene),
    colour = "#B2182B", nudge_x = 0.25, hjust = 0, fontface = "bold"
  ) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  coord_cartesian(clip = "off") +
  labs(
    x = "HV minus pure-UC average log2 fold change",
    y = expression(-log[10](adjusted~p)),
    title = "TM4SF1 is the strongest HV-enriched signal in the released table",
    subtitle = "Extreme values are capped at the 99.5th percentile for display"
  ) +
  theme_report() +
  theme(plot.margin = margin(5.5, 35, 5.5, 5.5))

save_plot(paper_de_plot, "01_paper_Figure5A_TM4SF1_volcano.png")
paper_de_plot
```

![Released Figure 5A differential-expression table. TM4SF1 is highlighted without recomputing the authors' private object.](/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/figures/paper-volcano-1.png)

## 2. Visualise TM4SF1 in the paper cohort


``` r
paper_cells <- readxl::read_excel(paper_workbook, sheet = "Figure5B") |>
  transmute(
    barcode = as.character(Barcode),
    paper_id = as.character(Name),
    source_expression = as.numeric(TM4SF1),
    histology_group = case_when(
      str_starts(paper_id, "VAR") ~ "HV",
      str_starts(paper_id, "UC") ~ "UC",
      TRUE ~ "other"
    )
  ) |>
  filter(histology_group %in% c("HV", "UC"))

paper_cell_patient <- paper_cells |>
  group_by(paper_id, histology_group) |>
  summarise(
    n_cells = n(),
    mean_source_expression = mean(source_expression),
    median_source_expression = median(source_expression),
    source_nonzero_fraction = mean(source_expression > 1e-8),
    .groups = "drop"
  )

write_tsv(
  paper_cell_patient,
  file.path(results_dir, "TM4SF1_paper_source_patient_summary.tsv")
)
```

The authors released normalized per-cell TM4SF1 values—not raw UMI counts—for 8553 tumour epithelial cells across 9 HV and 3 UC patients. The cell distributions reproduce the figure, but they are descriptive; patients, not cells, are the independent units. Raw DGE counts are used for detection and cohort modelling below.


``` r
paper_cell_plot <- ggplot(
  paper_cells,
  aes(histology_group, source_expression, fill = histology_group)
) +
  geom_violin(scale = "width", trim = TRUE, colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", alpha = 0.8) +
  geom_point(
    data = paper_cell_patient,
    inherit.aes = FALSE,
    mapping = aes(x = histology_group, y = mean_source_expression),
    position = position_jitter(width = 0.07, height = 0),
    shape = 21, fill = "black", colour = "white", size = 2.3
  ) +
  scale_fill_manual(values = c(HV = "#B2182B", UC = "#2166AC")) +
  labs(
    x = NULL,
    y = "Released TM4SF1 expression value",
    title = "TM4SF1 is enriched in HV cells but heterogeneous across patients"
  ) +
  theme_report() +
  guides(fill = "none")

save_plot(paper_cell_plot, "02_paper_Figure5B_TM4SF1_cells.png", 6.5, 5.5)
paper_cell_plot
```

![Released Figure 5B TM4SF1 values. The violin is cell-level source data; black points are patient means.](/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/figures/paper-cell-violin-1.png)


``` r
tm4_patient <- readr::read_tsv(
  file.path(case_dir, "data", "metadata", "TM4SF1_patient_summary.tsv"),
  show_col_types = FALSE
) |>
  filter(
    sample_type == "tumour",
    broad_cell_type == "epithelial",
    paper_id %in% c(sprintf("VAR%02d", 1:9), sprintf("UC%02d", 1:3))
  ) |>
  mutate(
    cohort_status = if_else(
      reconstructed_n150_threshold_met,
      ">=150 reconstructed cells",
      "<150 reconstructed cells"
    ),
    paper_id = factor(paper_id, levels = paper_id[order(CPM)])
  )

patient_cpm_plot <- ggplot(
  tm4_patient,
  aes(paper_id, CPM, colour = histology_group, shape = cohort_status)
) +
  geom_point(size = 3) +
  scale_y_log10(labels = label_number()) +
  scale_colour_manual(values = c(HV = "#B2182B", UC = "#2166AC")) +
  coord_flip() +
  labs(
    x = NULL,
    y = "TM4SF1 CPM (log10 scale)",
    colour = "Histology",
    shape = "Public-DGE cohort status",
    title = "The HV signal is patient-heterogeneous",
    subtitle = "VAR08 contributes most HV molecules; several HV tumours are low"
  ) +
  theme_report()

save_plot(patient_cpm_plot, "03_public_DGE_TM4SF1_patient_CPM.png", 8, 6)
patient_cpm_plot
```

![Patient-level raw-count pseudobulk contrast reconstructed from the deposited DGE matrices. UC01 is retained but marked because it has only 90 reconstructed tumour epithelial cells.](/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/figures/public-dge-patient-contrast-1.png)

### Which broad population carries TM4SF1?


``` r
population_source <- readr::read_tsv(
  file.path(
    case_dir, "data", "metadata", "TM4SF1_population_source_summary.tsv"
  ),
  show_col_types = FALSE
) |>
  mutate(
    sample_type = recode(
      sample_type,
      paired_normal = "paired/adjacent normal",
      tumour = "tumour"
    )
  )

population_plot <- ggplot(
  population_source,
  aes(broad_cell_type, pooled_CPM, fill = sample_type)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_text(
    aes(label = percent(pooled_detection_rate, accuracy = 1)),
    position = position_dodge(width = 0.8),
    vjust = -0.35, size = 3
  ) +
  scale_y_continuous(labels = label_number(), expand = expansion(mult = c(0, 0.14))) +
  scale_fill_manual(values = c(
    "tumour" = "#B2182B",
    "paired/adjacent normal" = "#67A9CF"
  )) +
  labs(
    x = NULL,
    y = "Pooled TM4SF1 CPM",
    fill = "Study compartment",
    title = "Endothelial TM4SF1 is high and nearly universally detected",
    subtitle = "Labels give the fraction of cells with non-zero counts"
  ) +
  theme_report()

save_plot(population_plot, "04_TM4SF1_broad_population_source.png", 8, 5.5)
population_plot
```

![Pooled TM4SF1 abundance and detection across broad reconstructed populations. Endothelium is a key cell class for off-tumour follow-up.](/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/figures/broad-population-source-1.png)

Because pooled CPM is composition-weighted, the same check is shown at patient level:


``` r
population_patient <- readr::read_tsv(
  file.path(case_dir, "data", "metadata", "TM4SF1_patient_summary.tsv"),
  show_col_types = FALSE
) |>
  filter(broad_cell_type %in% c(
    "epithelial", "endothelial", "stromal", "immune"
  )) |>
  mutate(sample_type = recode(
    sample_type,
    paired_normal = "paired/adjacent normal",
    tumour = "tumour"
  ))

population_patient_plot <- ggplot(
  population_patient,
  aes(broad_cell_type, CPM, colour = sample_type)
) +
  geom_boxplot(
    aes(group = interaction(broad_cell_type, sample_type)),
    outlier.shape = NA, width = 0.58, alpha = 0.18,
    position = position_dodge(width = 0.65)
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.12, dodge.width = 0.65
    ),
    alpha = 0.75, size = 1.8
  ) +
  scale_y_continuous(
    trans = scales::transform_log1p(),
    labels = label_number()
  ) +
  scale_colour_manual(values = c(
    "tumour" = "#B2182B",
    "paired/adjacent normal" = "#2166AC"
  )) +
  labs(
    x = NULL,
    y = "Patient TM4SF1 CPM (log1p scale)",
    colour = "Study compartment",
    title = "The endothelial signal is visible across patient pseudobulks"
  ) +
  theme_report()

save_plot(
  population_patient_plot,
  "04b_TM4SF1_broad_population_patient_CPM.png",
  8, 5.5
)
population_patient_plot
```

![Patient-level TM4SF1 CPM by broad population and study compartment. Points are patient pseudobulks; low-cell profiles remain descriptive.](/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/figures/broad-population-patient-distribution-1.png)

Together, the pooled and patient-level checks motivate the comparator strategy. Epithelium is the closest cell-class-matched comparator for tumour epithelium, while endothelium must be analysed separately because its expression may represent a vascular off-tumour exposure. Pooling these classes would hide the relevant biology.

## 3. Cohort-level `log(mu)` and posteriorHCA baseline

### Reusable helpers

The helper functions below deliberately remain in the report. They implement one transparent chain from raw patient counts to model-aligned `log(mu)`. Input validation and cache checks run silently in hidden chunks.


``` r
aggregate_ensembl_counts <- function(se) {
  counts <- as.matrix(SummarizedExperiment::assay(se, "counts"))
  features <- as.data.frame(SummarizedExperiment::rowData(se))
  ensembl <- as.character(features$ensembl_id)
  keep <- !is.na(ensembl) & nzchar(ensembl)

  out <- rowsum(counts[keep, , drop = FALSE], ensembl[keep], reorder = FALSE)
  storage.mode(out) <- "numeric"
  out
}

read_study_pseudobulk <- function(case_dir) {
  tumour_se <- readRDS(file.path(
    case_dir, "data", "processed", "pseudobulk_tumour_epithelial_counts.rds"
  ))
  normal_se <- readRDS(file.path(
    case_dir, "data", "processed",
    "pseudobulk_paired_normal_epithelial_counts.rds"
  ))

  tumour_counts <- aggregate_ensembl_counts(tumour_se)
  normal_counts <- aggregate_ensembl_counts(normal_se)
  tumour_meta <- as.data.frame(colData(tumour_se)) |>
    rownames_to_column("source_column") |>
    mutate(
      analysis_group = histology_group,
      study_state = "tumour epithelial",
      include_for_logmu = TRUE
    )
  normal_meta <- as.data.frame(colData(normal_se)) |>
    rownames_to_column("source_column") |>
    mutate(
      analysis_group = "paired_normal",
      study_state = "paired-normal epithelial",
      include_for_logmu = n_cells >= 20
    )

  common <- intersect(rownames(tumour_counts), rownames(normal_counts))
  normal_keep <- normal_meta$include_for_logmu
  counts <- cbind(
    tumour_counts[common, , drop = FALSE],
    normal_counts[common, normal_keep, drop = FALSE]
  )
  meta <- bind_rows(tumour_meta, normal_meta[normal_keep, , drop = FALSE]) |>
    mutate(
      analysis_column = c(
        paste0("tumour__", tumour_meta$patient_id),
        paste0("paired_normal_epithelial__", normal_meta$patient_id[normal_keep])
      ),
      strict_n150 = study_state != "tumour epithelial" |
        reconstructed_n150_threshold_met
    )

  colnames(counts) <- meta$analysis_column
  rownames(meta) <- meta$analysis_column
  list(counts = counts, metadata = meta)
}

align_to_hca_reference <- function(study_counts, reference_counts, reference_name) {
  shared <- intersect(rownames(study_counts), names(reference_counts))
  combined <- cbind(
    study_counts[shared, , drop = FALSE],
    setNames(reference_counts[shared], shared)
  )
  colnames(combined)[ncol(combined)] <- reference_name

  norm_factors <- edgeR::calcNormFactors(
    combined,
    refColumn = ncol(combined),
    method = "TMMwsp"
  )
  library_size <- colSums(combined)
  effective_size <- library_size * norm_factors
  multiplier <- effective_size[[reference_name]] / effective_size
  offset <- log(1 / multiplier)

  list(
    counts = combined[, colnames(study_counts), drop = FALSE],
    offset = offset[colnames(study_counts)],
    multiplier = multiplier[colnames(study_counts)],
    shared_features = shared
  )
}

fit_shared_dispersion <- function(counts, metadata, offset, feature_id) {
  design_group <- factor(metadata$analysis_group)
  design <- model.matrix(~ 0 + design_group)
  dge <- edgeR::DGEList(counts = counts)
  dge$offset <- matrix(offset, nrow(counts), ncol(counts), byrow = TRUE)

  keep <- edgeR::filterByExpr(dge, design = design)
  keep[match(feature_id, rownames(dge))] <- TRUE
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
  invisible(edgeR::glmQLFit(dge, design, robust = TRUE))

  dispersion <- if (!is.null(dge$trended.dispersion)) {
    dge$trended.dispersion
  } else if (!is.null(dge$tagwise.dispersion)) {
    dge$tagwise.dispersion
  } else {
    rep_len(dge$common.dispersion, nrow(dge))
  }
  names(dispersion) <- rownames(dge)
  unname(dispersion[[feature_id]])
}

draw_dirichlet_weights <- function(n) {
  weights <- rexp(n)
  n * weights / sum(weights)
}

bootstrap_cohort_logmu <- function(
    counts, offset, dispersion, columns, n_boot, seed
) {
  y <- counts[tm4_id, columns, drop = FALSE]
  o <- offset[columns]
  set.seed(seed)
  vapply(seq_len(n_boot), function(i) {
    edgeR::mglmOneGroup(
      y,
      offset = o,
      dispersion = dispersion,
      weights = draw_dirichlet_weights(length(columns))
    )[[1]]
  }, numeric(1))
}

summarise_draws <- function(draws) {
  tibble(
    median = median(draws),
    lower = quantile(draws, 0.025),
    upper = quantile(draws, 0.975)
  )
}

pair_log2_contrast <- function(x, y, n = 2000L, seed = 1L) {
  set.seed(seed)
  delta <- (
    sample(x, n, replace = TRUE) - sample(y, n, replace = TRUE)
  ) / log(2)
  tibble(
    estimate_log2 = median(delta),
    lower_log2 = quantile(delta, 0.025),
    upper_log2 = quantile(delta, 0.975),
    probability_positive = mean(delta > 0)
  )
}
```


``` r
read_hca_bundle <- function(cell_type, model_root, feature_id) {
  store <- file.path(model_root, paste0("V1_", cell_type), "_targets")
  registry <- read.csv(
    file.path(model_root, "_PAK2_TM4SF1_modelled_cell_types.csv"),
    check.names = FALSE
  )
  target_name <- registry |>
    filter(
      cell_type == paste0("V1_", .env$cell_type),
      ensembl_id == feature_id,
      isTRUE(in_feature_df) | in_feature_df
    ) |>
    pull(estimates_branch) |>
    unique()

  chunk <- targets::tar_read_raw(target_name[[1]], store = store)
  fit_index <- match(feature_id, chunk$.feature)
  fit <- chunk$brms_fit[[fit_index]]
  pseudobulk <- targets::tar_read_raw("pseudobulk_sample", store = store)
  sample_meta <- as.data.frame(colData(pseudobulk)) |>
    rownames_to_column("sample_column")

  reference_index <- which(sample_meta$offset == 0)
  reference_name <- sample_meta$sample_column[[reference_index]]
  reference_counts <- as.numeric(
    assay(pseudobulk, "counts")[, reference_index, drop = TRUE]
  )
  names(reference_counts) <- rownames(pseudobulk)

  profiles <- sample_meta |>
    transmute(
      sample_column,
      dataset_id = as.character(dataset_id),
      donor_id = as.character(donor_id),
      tissue = as.character(tissue),
      tissue_groups = as.character(tissue_groups),
      age_decade = as.character(age_decade),
      sex = as.character(sex),
      ethnicity_groups = as.character(ethnicity_groups),
      disease_group = as.character(disease_groups),
      assay_group = as.character(assay_groups)
    ) |>
    filter(
      disease_group == "Normal",
      assay_group == "10x Genomics 3",
      if_all(
        c(
          dataset_id, donor_id, tissue_groups,
          age_decade, sex, ethnicity_groups
        ),
        ~ !is.na(.x) & nzchar(.x)
      )
    ) |>
    mutate(donor_key = paste(dataset_id, donor_id, sep = "___")) |>
    distinct(
      donor_key, tissue_groups, age_decade, sex, ethnicity_groups,
      .keep_all = TRUE
    )

  newdata <- profiles |>
    transmute(
      counts = 1,
      age_decade,
      sex,
      disease_groups_altered = "Normal",
      ethnicity_groups,
      assay_groups_altered = "10x Genomics 3",
      dataset_id_altered = "__TM4SF1_new_study__",
      tissue_groups,
      offset = 0
    )

  fit_draws <- posterior::as_draws_array(fit)
  draw_diagnostics <- posterior::summarise_draws(fit_draws)
  nuts <- brms::nuts_params(fit)
  finite_max <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) max(x) else NA_real_
  }
  finite_min <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) min(x) else NA_real_
  }
  max_rhat_value <- finite_max(draw_diagnostics$rhat)
  min_ess_bulk_value <- finite_min(draw_diagnostics$ess_bulk)
  min_ess_tail_value <- finite_min(draw_diagnostics$ess_tail)
  divergence_count <- sum(
    nuts$Parameter == "divergent__" & nuts$Value == 1,
    na.rm = TRUE
  )
  diagnostics <- tibble(
    cell_type,
    n_model_rows = nrow(fit$data),
    n_posterior_draws = posterior::ndraws(fit_draws),
    max_rhat = max_rhat_value,
    min_ess_bulk = min_ess_bulk_value,
    min_ess_tail = min_ess_tail_value,
    divergences = divergence_count,
    readiness_policy = "Rhat<=1.01; bulk/tail ESS>=400; zero divergences",
    model_ready = all(is.finite(c(
      max_rhat_value, min_ess_bulk_value, min_ess_tail_value,
      divergence_count
    ))) &
      max_rhat_value <= 1.01 &
      min_ess_bulk_value >= 400 &
      min_ess_tail_value >= 400 &
      divergence_count == 0
  )

  list(
    fit = fit,
    newdata = newdata,
    profiles = profiles,
    reference_counts = reference_counts,
    reference_name = reference_name,
    model_store = store,
    target_name = target_name[[1]],
    diagnostics = diagnostics
  )
}

posterior_hca_logmu <- function(fit, newdata, seed) {
  set.seed(seed)
  brms::posterior_linpred(
    fit,
    newdata = newdata,
    transform = FALSE,
    re_formula = NULL,
    allow_new_levels = TRUE,
    sample_new_levels = "gaussian"
  )
}

summarise_hca_groups <- function(
    logmu, profiles, dimension, group_vars, hv_draws,
    threshold, min_donors, seed
) {
  if (!length(group_vars)) {
    groups <- list(comprehensive = seq_len(nrow(profiles)))
  } else {
    key <- do.call(
      interaction,
      c(
        profiles[, group_vars, drop = FALSE],
        list(drop = TRUE, lex.order = TRUE, sep = " | ")
      )
    )
    groups <- split(seq_len(nrow(profiles)), key)
  }

  map2_dfr(groups, seq_along(groups), function(index, group_number) {
    hca_draws <- rowMeans(logmu[, index, drop = FALSE])
    paired <- pair_log2_contrast(
      hca_draws, hv_draws,
      n = max(length(hv_draws), length(hca_draws)),
      seed = seed + group_number
    )
    first <- profiles[index[[1]], , drop = FALSE]
    values <- setNames(
      map(group_vars, ~ first[[.x]][[1]]),
      group_vars
    )
    value_or_na <- function(name) {
      if (name %in% names(values)) as.character(values[[name]]) else NA_character_
    }
    support <- n_distinct(profiles$donor_key[index])

    tibble(
      dimension,
      tissue_groups = value_or_na("tissue_groups"),
      sex = value_or_na("sex"),
      age_decade = value_or_na("age_decade"),
      ethnicity_groups = value_or_na("ethnicity_groups"),
      n_profile_rows = length(index),
      n_donor_keys = support,
      median_log_mu = median(hca_draws),
      lower_log_mu = quantile(hca_draws, 0.025),
      upper_log_mu = quantile(hca_draws, 0.975),
      threshold_log_mu = threshold,
      probability_above_threshold = mean(hca_draws >= threshold),
      probability_hca_ge_random_hv = mean(ecdf(hv_draws)(hca_draws)),
      estimate_log2_vs_hv = paired$estimate_log2,
      lower_log2_vs_hv = paired$lower_log2,
      upper_log2_vs_hv = paired$upper_log2,
      supported = dimension == "comprehensive" | support >= min_donors
    )
  })
}

compute_tissue_rule_sensitivity <- function(
    logmu, profiles, hv_draws, cell_type,
    threshold_quantiles = c(0.01, 0.05, 0.10),
    minimum_donor_keys = c(5L, 10L, 20L),
    alert_probabilities = c(0.50, 0.80, 0.95)
) {
  tissue_index <- split(
    seq_len(nrow(profiles)),
    factor(profiles$tissue_groups)
  )
  tissue_draws <- map(tissue_index, ~ rowMeans(logmu[, .x, drop = FALSE]))
  tissue_support <- map_int(
    tissue_index,
    ~ n_distinct(profiles$donor_key[.x])
  )

  tidyr::crossing(
    threshold_quantile = threshold_quantiles,
    minimum_donor_keys = minimum_donor_keys,
    alert_probability_cutoff = alert_probabilities
  ) |>
    pmap_dfr(function(
        threshold_quantile, minimum_donor_keys, alert_probability_cutoff
    ) {
      threshold <- unname(quantile(hv_draws, threshold_quantile))
      exceedance <- map_dbl(
        tissue_draws,
        ~ mean(.x >= threshold)
      )
      tibble(
        cell_type,
        threshold_quantile,
        threshold_log_mu = threshold,
        minimum_donor_keys,
        alert_probability_cutoff,
        supported_tissues = sum(tissue_support >= minimum_donor_keys),
        prototype_tissue_candidates = sum(
          tissue_support >= minimum_donor_keys &
            exceedance >= alert_probability_cutoff
        )
      )
    })
}
```


``` r
run_one_cell_type <- function(cell_type, study, params, model_root) {
  bundle <- read_hca_bundle(cell_type, model_root, tm4_id)
  aligned <- align_to_hca_reference(
    study$counts,
    bundle$reference_counts,
    bundle$reference_name
  )

  dispersion <- fit_shared_dispersion(
    aligned$counts,
    study$metadata,
    aligned$offset,
    tm4_id
  )

  groups <- list(
    "HV tumour epithelial" = which(
      study$metadata$study_state == "tumour epithelial" &
        study$metadata$histology_group == "HV"
    ),
    "UC tumour epithelial (3 published IDs)" = which(
      study$metadata$study_state == "tumour epithelial" &
        study$metadata$histology_group == "UC"
    ),
    "UC tumour epithelial (2 public-DGE >=150)" = which(
      study$metadata$study_state == "tumour epithelial" &
        study$metadata$histology_group == "UC" &
        study$metadata$strict_n150
    ),
    "paired-normal epithelial (>=20 cells)" = which(
      study$metadata$study_state == "paired-normal epithelial"
    )
  )

  study_draws <- imap_dfr(groups, function(index, group_name) {
    draws <- bootstrap_cohort_logmu(
      aligned$counts,
      aligned$offset,
      dispersion,
      colnames(aligned$counts)[index],
      params$n_boot,
      params$seed + which(names(groups) == group_name)
    )
    tibble(
      cell_type,
      group = group_name,
      draw = seq_along(draws),
      log_mu = draws,
      n_patients = length(index)
    )
  })

  hv_draws <- study_draws |>
    filter(group == "HV tumour epithelial") |>
    pull(log_mu)
  threshold <- unname(quantile(hv_draws, params$threshold_quantile))

  hca_logmu <- posterior_hca_logmu(
    bundle$fit,
    bundle$newdata,
    seed = params$seed + match(cell_type, params$model_cell_types) * 100L
  )

  profile_summary <- bind_rows(
    summarise_hca_groups(
      hca_logmu, bundle$profiles, "comprehensive", character(), hv_draws,
      threshold, params$min_donor_keys, params$seed + 1000L
    ),
    summarise_hca_groups(
      hca_logmu, bundle$profiles, "tissue", "tissue_groups", hv_draws,
      threshold, params$min_donor_keys, params$seed + 2000L
    ),
    summarise_hca_groups(
      hca_logmu, bundle$profiles, "tissue_sex",
      c("tissue_groups", "sex"), hv_draws,
      threshold, params$min_donor_keys, params$seed + 3000L
    ),
    summarise_hca_groups(
      hca_logmu, bundle$profiles, "tissue_age",
      c("tissue_groups", "age_decade"), hv_draws,
      threshold, params$min_donor_keys, params$seed + 4000L
    ),
    summarise_hca_groups(
      hca_logmu, bundle$profiles, "tissue_ethnicity",
      c("tissue_groups", "ethnicity_groups"), hv_draws,
      threshold, params$min_donor_keys, params$seed + 5000L
    )
  ) |>
    mutate(
      cell_type = cell_type,
      model_ready = bundle$diagnostics$model_ready[[1]],
      interpretable_for_demographic_alert = !(
        dimension == "tissue_ethnicity" &
          ethnicity_groups == "Other/Unknown"
      ),
      prototype_overlap_candidate = supported &
        interpretable_for_demographic_alert &
        probability_above_threshold >= params$alert_probability,
      expression_overlap_alert = prototype_overlap_candidate & model_ready,
      .before = 1
    )

  hca_comprehensive <- rowMeans(hca_logmu)
  renal_index <- which(bundle$profiles$tissue_groups == "renal system")
  hca_renal <- if (length(renal_index)) {
    rowMeans(hca_logmu[, renal_index, drop = FALSE])
  } else {
    numeric()
  }

  rule_sensitivity <- compute_tissue_rule_sensitivity(
    hca_logmu,
    bundle$profiles,
    hv_draws,
    cell_type
  ) |>
    mutate(
      model_ready = bundle$diagnostics$model_ready[[1]],
      analysis_id = analysis_id,
      .after = cell_type
    )

  contrast_draws <- study_draws
  if (cell_type != "epithelial") {
    contrast_draws <- contrast_draws |>
      filter(group != "paired-normal epithelial (>=20 cells)")
  }

  cohort_contrasts <- contrast_draws |>
    group_split(group) |>
    map_dfr(function(group_data) {
      comparison_draws <- hca_comprehensive
      comparator <- paste0(
        "posteriorHCA Normal/10x3-supported comprehensive ", cell_type
      )
      if (
        cell_type == "epithelial" &&
        unique(group_data$group) == "paired-normal epithelial (>=20 cells)" &&
        length(hca_renal)
      ) {
        comparison_draws <- hca_renal
        comparator <- paste(
          "posteriorHCA Normal/10x3-supported renal-system",
          "epithelial proxy"
        )
      }
      pair_log2_contrast(
        group_data$log_mu,
        comparison_draws,
        n = params$n_boot,
        seed = params$seed + 6000L + match(
          unique(group_data$group), names(groups)
        )
      ) |>
        mutate(
          cell_type,
          study_group = unique(group_data$group),
          comparator,
          n_study_patients = unique(group_data$n_patients),
          .before = 1
        )
    })

  threshold_registry <- tibble(
    gene_symbol = tm4_symbol,
    ensembl_id = tm4_id,
    cell_type_model = cell_type,
    model_store = bundle$model_store,
    model_target = bundle$target_name,
    normalization_reference = bundle$reference_name,
    normalization_method = "TMMwsp; reference-anchored offset=log(1/multiplier)",
    study_biological_unit = "patient",
    tumour_expression_anchor = "9-patient HV tumour epithelial cohort location",
    threshold_quantile = params$threshold_quantile,
    threshold_log_mu = threshold,
    n_patient_bootstrap_draws = params$n_boot,
    random_seed = params$seed,
    alert_probability_cutoff = params$alert_probability,
    minimum_donor_keys = params$min_donor_keys,
    model_ready = bundle$diagnostics$model_ready[[1]],
    diagnostic_gate = bundle$diagnostics$readiness_policy[[1]],
    input_and_model_fingerprint = input_fingerprint,
    shared_feature_count = length(aligned$shared_features),
    shared_feature_sha256 = digest::digest(
      sort(aligned$shared_features), algo = "sha256"
    ),
    profile_weighting = paste(
      "equal weight per distinct dataset-donor x tissue-group x",
      "age-decade x sex x ethnicity profile"
    ),
    prediction_policy = paste(
      "Normal; 10x Genomics 3; offset=0; observed donor-tissue profiles;",
      "unseen dataset random effect"
    ),
    analysis_version,
    analysis_id
  )

  support_profiles <- bundle$profiles |>
    mutate(cell_type = cell_type, .before = 1)

  list(
    analysis_version = analysis_version,
    analysis_id = analysis_id,
    cell_type = cell_type,
    study_draws = study_draws,
    profile_summary = profile_summary,
    rule_sensitivity = rule_sensitivity,
    cohort_contrasts = cohort_contrasts,
    threshold_registry = threshold_registry,
    diagnostics = bundle$diagnostics,
    support_profiles = support_profiles,
    reference = tibble(
      cell_type,
      reference_sample = bundle$reference_name,
      shared_features = length(aligned$shared_features),
      shared_feature_sha256 = digest::digest(
        sort(aligned$shared_features), algo = "sha256"
      ),
      tm4sf1_dispersion = dispersion
    )
  )
}
```

### Build or load the model-aligned results



### What the baseline represents

The comprehensive reference is an equal-weight empirical standardisation over distinct Normal-labelled donor–tissue–demographic profiles observed with 10x Genomics 3 **within the selected cell type**. Here, "comprehensive" means the eligible atlas support, not complete tissue or population coverage. Epithelial and endothelial biology remain separate. Predictions fix disease to Normal and offset to zero, and draw a new dataset-level random effect as in the SAVI cohort comparison. Here, 10x Genomics 3 is a fixed posteriorHCA reporting scale, not a claim that the bladder study used that assay.

`posterior_linpred(transform = FALSE)` estimates the conditional latent negative-binomial `log(mu)` on the fitted model scale. It is not the marginal zero-inflated expected count, detection probability, surface-protein abundance, antigen density or toxicity. It is used here only for alignment with the SAVI cohort-comparison estimand.


``` r
support_profiles |>
  group_by(cell_type) |>
  summarise(
    donor_keys = n_distinct(donor_key),
    donor_tissue_profiles = n(),
    tissue_groups = n_distinct(tissue_groups),
    age_decades = n_distinct(age_decade),
    sexes = n_distinct(sex),
    ethnicity_groups = n_distinct(ethnicity_groups),
    .groups = "drop"
  ) |>
  kable(
    caption = paste(
      "Normal-labelled, 10x3-supported profiles used for empirical standardisation.",
      "A donor key is dataset_id × donor_id."
    )
  )
```



Table: Normal-labelled, 10x3-supported profiles used for empirical standardisation. A donor key is dataset_id × donor_id.

|cell_type   | donor_keys| donor_tissue_profiles| tissue_groups| age_decades| sexes| ethnicity_groups|
|:-----------|----------:|---------------------:|-------------:|-----------:|-----:|----------------:|
|endothelial |        903|                   979|            29|          10|     2|                8|
|epithelial  |        536|                   573|            21|           9|     2|                7|

Bladder is not a fitted posteriorHCA tissue level. The few Normal-labelled 10x3 `bladder organ` profiles are collapsed into `renal system`, together with renal tissues. Accordingly, the paired-normal comparison below is a calibration against an imperfect renal-system epithelial proxy; this report does not make bladder-specific demographic claims.


``` r
retained_tumour_ids <- study$metadata |>
  filter(
    study_state == "tumour epithelial",
    strict_n150,
    paper_id %in% c(sprintf("VAR%02d", 1:9), sprintf("UC%02d", 1:3))
  ) |>
  pull(paper_id)

paired_normal_calibration <- study$metadata |>
  filter(study_state == "paired-normal epithelial") |>
  transmute(
    paper_id,
    epithelial_cells = n_cells,
    raw_library_size = library_size,
    paired_with_retained_tumour = paper_id %in% retained_tumour_ids
  )

paired_normal_calibration |>
  kable(
    caption = paste(
      "Paired/adjacent-normal epithelial calibration samples with at least",
      "20 cells. Only VAR07 has a retained tumour-cohort counterpart."
    )
  )
```



Table: Paired/adjacent-normal epithelial calibration samples with at least 20 cells. Only VAR07 has a retained tumour-cohort counterpart.

|                                |paper_id | epithelial_cells| raw_library_size|paired_with_retained_tumour |
|:-------------------------------|:--------|----------------:|----------------:|:---------------------------|
|paired_normal_epithelial__UC04  |UC04     |              190|          3807119|FALSE                       |
|paired_normal_epithelial__VAR07 |VAR07    |               30|            64199|TRUE                        |
|paired_normal_epithelial__VAR10 |VAR10    |               29|           153881|FALSE                       |

### Model-readiness note


``` r
model_diagnostics |>
  mutate(
    max_rhat = round(max_rhat, 3),
    min_ess_bulk = round(min_ess_bulk, 1),
    min_ess_tail = round(min_ess_tail, 1)
  ) |>
  kable(
    caption = "Diagnostics of the current local TM4SF1 posteriorHCA fits."
  )
```



Table: Diagnostics of the current local TM4SF1 posteriorHCA fits.

|cell_type   | n_model_rows| n_posterior_draws| max_rhat| min_ess_bulk| min_ess_tail| divergences|readiness_policy                                 |model_ready |
|:-----------|------------:|-----------------:|--------:|------------:|------------:|-----------:|:------------------------------------------------|:-----------|
|epithelial  |         1545|               400|    1.065|         55.2|         33.2|          34|Rhat<=1.01; bulk/tail ESS>=400; zero divergences |FALSE       |
|endothelial |         2308|               400|    1.102|         17.7|         68.3|           6|Rhat<=1.01; bulk/tail ESS>=400; zero divergences |FALSE       |

The current local fits contain 400 retained posterior draws per model. **Neither current fit passes the registered readiness gate.** `TM4SF1_expression_overlap_alerts.tsv` contains only candidates from models that pass this gate and is therefore empty when none is ready. The figures and tables below show ungated prototype overlap candidates, not translational alerts; failed models must be refitted with longer chains before their alerts can be evaluated.

## 4. Paper cohorts versus posteriorHCA


``` r
cohort_contrasts |>
  mutate(
    across(c(estimate_log2, lower_log2, upper_log2), ~ round(.x, 2)),
    probability_positive = percent(probability_positive, accuracy = 0.1)
  ) |>
  select(
    cell_type, study_group, comparator, n_study_patients,
    estimate_log2, lower_log2, upper_log2, probability_positive
  ) |>
  kable(
    caption = paste(
      "Patient-cohort minus posteriorHCA latent log-mu contrasts on a log2 scale.",
      "Intervals independently pair study bootstrap and HCA posterior draws."
    )
  )
```



Table: Patient-cohort minus posteriorHCA latent log-mu contrasts on a log2 scale. Intervals independently pair study bootstrap and HCA posterior draws.

|cell_type   |study_group                               |comparator                                                       | n_study_patients| estimate_log2| lower_log2| upper_log2|probability_positive |
|:-----------|:-----------------------------------------|:----------------------------------------------------------------|----------------:|-------------:|----------:|----------:|:--------------------|
|epithelial  |HV tumour epithelial                      |posteriorHCA Normal/10x3-supported comprehensive epithelial      |                9|          1.86|      -5.10|       8.19|70.0%                |
|epithelial  |UC tumour epithelial (2 public-DGE >=150) |posteriorHCA Normal/10x3-supported comprehensive epithelial      |                2|         -2.70|      -9.41|       4.58|19.8%                |
|epithelial  |UC tumour epithelial (3 published IDs)    |posteriorHCA Normal/10x3-supported comprehensive epithelial      |                3|         -0.92|      -7.76|       5.87|38.9%                |
|epithelial  |paired-normal epithelial (>=20 cells)     |posteriorHCA Normal/10x3-supported renal-system epithelial proxy |                3|          2.32|      -4.45|       8.91|74.4%                |
|endothelial |HV tumour epithelial                      |posteriorHCA Normal/10x3-supported comprehensive endothelial     |                9|          1.36|      -2.45|       5.44|74.2%                |
|endothelial |UC tumour epithelial (2 public-DGE >=150) |posteriorHCA Normal/10x3-supported comprehensive endothelial     |                2|         -3.12|      -7.13|       0.95|7.3%                 |
|endothelial |UC tumour epithelial (3 published IDs)    |posteriorHCA Normal/10x3-supported comprehensive endothelial     |                3|         -1.30|      -5.27|       2.75|27.6%                |


``` r
cohort_plot_data <- cohort_contrasts |>
  filter(
    study_group != "paired-normal epithelial (>=20 cells)" |
      cell_type == "epithelial"
  ) |>
  mutate(
    study_label = recode(
      study_group,
      "HV tumour epithelial" = "HV tumour (9 patients)",
      "UC tumour epithelial (2 public-DGE >=150)" =
        "UC tumour (2 public-DGE patients with >=150 cells)",
      "UC tumour epithelial (3 published IDs)" =
        "UC tumour (3 published IDs; UC01 has 90 cells)",
      "paired-normal epithelial (>=20 cells)" =
        "paired/adjacent normal (3 patients)"
    ),
    reference_label = case_when(
      str_detect(comparator, "renal-system") ~
        "HCA renal-system epithelial proxy",
      cell_type == "epithelial" ~ "HCA comprehensive epithelial reference",
      TRUE ~ "HCA comprehensive endothelial reference"
    ),
    comparison = str_wrap(paste(study_label, "vs", reference_label), 46),
    comparison = factor(comparison, levels = rev(unique(comparison)))
  )

cohort_plot <- ggplot(
  cohort_plot_data,
  aes(estimate_log2, comparison, colour = cell_type)
) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbar(
    aes(xmin = lower_log2, xmax = upper_log2),
    orientation = "y", width = 0.18
  ) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c(
    epithelial = "#B2182B", endothelial = "#2166AC"
  )) +
  labs(
    x = "Study minus posteriorHCA log2 latent expression",
    y = NULL,
    colour = "posteriorHCA model",
    title = "Study cohorts versus separate epithelial/endothelial references",
    subtitle = "Normal-labelled, 10x3-supported posteriorHCA scale"
  ) +
  theme_report()

save_plot(cohort_plot, "05_TM4SF1_cohort_vs_posteriorHCA_logmu.png", 10, 7)
cohort_plot
```

![Cohort-level TM4SF1 contrasts. Positive values indicate higher expression in the study group than the model-aligned posteriorHCA comparator.](/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/figures/cohort-contrast-forest-1.png)

The three-sample paired/adjacent-normal result is a **diagnostic calibration check**, not a healthy-control effect estimate or a matched-pair analysis. A material shift from the posteriorHCA renal proxy can arise from adjacent-tissue biology, bladder-versus-kidney aggregation, platform composition or study effects. Tumour-versus-HCA differences should therefore be read as modular cross-study contrasts rather than a jointly fitted causal effect.

## 5. A reusable threshold contract

The report defines the screening threshold separately for each posteriorHCA cell-type model:

\[
T_c = Q_q\{\log \mu_{\mathrm{HV},c}^{(b)}\},
\]

where \(q =\) 5%, \(c\) is the posteriorHCA model, and the distribution is the nine-patient HV tumour-epithelial Bayesian bootstrap after alignment to that model's TMMwsp reference. A supported Normal-labelled profile qualifies as a prototype candidate when

\[
P(\log \mu_{\mathrm{Normal},c} \geq T_c \mid \mathrm{data}) \geq p_{\min},
\]

with \(p_{\min} =\) 80%. The selected \(q\) is a lower posterior quantile of uncertainty in the nine-patient cohort-location estimator. It is not a quantile of patient expression, a minimum antigen-positive level, or a biological or clinical toxicity threshold. Both the quantile and probability are report parameters, not hidden constants.


``` r
threshold_registry |>
  transmute(
    cell_type_model,
    normalization_reference,
    threshold_quantile,
    threshold_log_mu = round(threshold_log_mu, 3),
    alert_probability_cutoff,
    minimum_donor_keys,
    model_ready,
    diagnostic_gate,
    analysis_version
  ) |>
  kable(
    caption = "Model-specific TM4SF1 threshold registry."
  )
```



Table: Model-specific TM4SF1 threshold registry.

|cell_type_model |normalization_reference                        | threshold_quantile| threshold_log_mu| alert_probability_cutoff| minimum_donor_keys|model_ready |diagnostic_gate                                  |analysis_version           |
|:---------------|:----------------------------------------------|------------------:|----------------:|------------------------:|------------------:|:-----------|:------------------------------------------------|:--------------------------|
|epithelial      |04e410cbab17c7d05877161100a5d1e1___epithelial  |               0.05|            8.367|                      0.8|                 10|FALSE       |Rhat<=1.01; bulk/tail ESS>=400; zero divergences |tm4sf1-logmu-v5-2026-08-18 |
|endothelial     |e9a38836dc09755c74363eec65317e6f___endothelial |               0.05|            6.071|                      0.8|                 10|FALSE       |Rhat<=1.01; bulk/tail ESS>=400; zero divergences |tm4sf1-logmu-v5-2026-08-18 |

The prototype exceedance-probability and donor-support cutoffs are screening choices, so their effect is reported rather than hidden:


``` r
rule_sensitivity |>
  mutate(
    threshold_quantile = percent(threshold_quantile, accuracy = 1),
    probability_rule = paste0(
      "candidates at P>=",
      percent(alert_probability_cutoff, accuracy = 1)
    )
  ) |>
  select(
    cell_type, model_ready, threshold_quantile, threshold_log_mu,
    minimum_donor_keys, supported_tissues, probability_rule,
    prototype_tissue_candidates
  ) |>
  pivot_wider(
    names_from = probability_rule,
    values_from = prototype_tissue_candidates
  ) |>
  kable(
    caption = paste(
      "Number of prototype tissue candidates under alternative HV-anchor",
      "quantiles, support rules and posterior-probability cutoffs."
    )
  )
```



Table: Number of prototype tissue candidates under alternative HV-anchor quantiles, support rules and posterior-probability cutoffs.

|cell_type   |model_ready |threshold_quantile | threshold_log_mu| minimum_donor_keys| supported_tissues| candidates at P>=50%| candidates at P>=80%| candidates at P>=95%|
|:-----------|:-----------|:------------------|----------------:|------------------:|-----------------:|--------------------:|--------------------:|--------------------:|
|epithelial  |FALSE       |1%                 |         7.957641|                  5|                13|                    2|                    0|                    0|
|epithelial  |FALSE       |1%                 |         7.957641|                 10|                11|                    1|                    0|                    0|
|epithelial  |FALSE       |1%                 |         7.957641|                 20|                 5|                    1|                    0|                    0|
|epithelial  |FALSE       |5%                 |         8.367083|                  5|                13|                    0|                    0|                    0|
|epithelial  |FALSE       |5%                 |         8.367083|                 10|                11|                    0|                    0|                    0|
|epithelial  |FALSE       |5%                 |         8.367083|                 20|                 5|                    0|                    0|                    0|
|epithelial  |FALSE       |10%                |         8.532771|                  5|                13|                    0|                    0|                    0|
|epithelial  |FALSE       |10%                |         8.532771|                 10|                11|                    0|                    0|                    0|
|epithelial  |FALSE       |10%                |         8.532771|                 20|                 5|                    0|                    0|                    0|
|endothelial |FALSE       |1%                 |         5.664328|                  5|                18|                   14|                    5|                    1|
|endothelial |FALSE       |1%                 |         5.664328|                 10|                13|                    9|                    3|                    1|
|endothelial |FALSE       |1%                 |         5.664328|                 20|                11|                    7|                    1|                    0|
|endothelial |FALSE       |5%                 |         6.070556|                  5|                18|                   13|                    4|                    0|
|endothelial |FALSE       |5%                 |         6.070556|                 10|                13|                    8|                    2|                    0|
|endothelial |FALSE       |5%                 |         6.070556|                 20|                11|                    6|                    0|                    0|
|endothelial |FALSE       |10%                |         6.246990|                  5|                18|                   13|                    3|                    0|
|endothelial |FALSE       |10%                |         6.246990|                 10|                13|                    8|                    2|                    0|
|endothelial |FALSE       |10%                |         6.246990|                 20|                11|                    6|                    0|                    0|

The raw `log(mu)` cutoff is reusable only with the recorded model store, feature, reference sample, assay policy, offset-zero convention and new-study random-effect policy. Across a different posteriorHCA release or cell-type model, the study counts must be realigned and \(T_c\) recalibrated. The estimand and recomputation procedure are portable; numerical log2 ratios, exceedance probabilities and cutoffs must be recomputed for every model or baseline release.

## 6. Normal-labelled tissue and observed-demographic profile landscape

Only observed Normal/10x3 donor–tissue profiles enter the marginalisation. These are descriptive empirical-profile averages: tissue × sex, age or ethnicity rows can differ in their observed mix of other covariates and datasets, so they are not controlled demographic effects. A tissue-demographic result must contain at least 10 distinct dataset × donor keys before it can qualify as an ungated prototype candidate; a model-ready alert additionally requires the diagnostic gate. Sparse rows remain in the full output with `supported = FALSE` but are suppressed from candidate plots. `Other/Unknown` ethnicity is retained descriptively but cannot qualify as a candidate.

### Tissue and cell type


``` r
tissue_plot_data <- profile_summary |>
  filter(dimension == "tissue", supported) |>
  group_by(cell_type) |>
  mutate(tissue_rank = rank(-median_log_mu, ties.method = "first")) |>
  ungroup() |>
  mutate(
    tissue_label = paste0(tissue_groups, " (n=", n_donor_keys, ")"),
    tissue_label = reorder(tissue_label, median_log_mu)
  )

tissue_plot <- ggplot(
  tissue_plot_data,
  aes(median_log_mu, tissue_label, colour = prototype_overlap_candidate)
) +
  geom_errorbar(
    aes(xmin = lower_log_mu, xmax = upper_log_mu),
    orientation = "y", width = 0.16, colour = "grey65"
  ) +
  geom_point(size = 2.3) +
  geom_vline(
    aes(xintercept = threshold_log_mu),
    linetype = 2, colour = "#B2182B"
  ) +
  facet_wrap(~ cell_type, scales = "free_y") +
  scale_colour_manual(
    values = c(`TRUE` = "#B2182B", `FALSE` = "#4D4D4D"),
    labels = c(
      `TRUE` = "ungated prototype candidate",
      `FALSE` = "below prototype rule"
    )
  ) +
  labs(
    x = "Normal-labelled posteriorHCA latent log(mu), offset 0",
    y = NULL,
    colour = NULL,
    title = "TM4SF1 varies by tissue and cell class"
  ) +
  theme_report()

save_plot(tissue_plot, "06_posteriorHCA_TM4SF1_tissue_profiles.png", 11, 8)
tissue_plot
```

![Normal-labelled posteriorHCA tissue profiles relative to the model-specific HV tumour-expression anchor. Highlighted points pass the donor-support and posterior-probability rules but remain ungated by model readiness.](/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/figures/tissue-profile-plot-1.png)

### Observed sex profiles within tissue


``` r
sex_plot_data <- profile_summary |>
  filter(dimension == "tissue_sex", supported) |>
  group_by(cell_type, tissue_groups) |>
  filter(n_distinct(sex) == 2) |>
  ungroup() |>
  mutate(
    tissue_groups = str_wrap(tissue_groups, 24),
    tissue_groups = reorder(tissue_groups, median_log_mu)
  )

sex_plot <- ggplot(
  sex_plot_data,
  aes(median_log_mu, tissue_groups, colour = sex)
) +
  geom_errorbar(
    aes(xmin = lower_log_mu, xmax = upper_log_mu),
    orientation = "y", width = 0.14,
    position = position_dodge(width = 0.5),
    alpha = 0.65
  ) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_vline(
    aes(xintercept = threshold_log_mu),
    linetype = 2, colour = "grey45"
  ) +
  facet_wrap(~ cell_type, scales = "free_y") +
  scale_colour_manual(values = c(female = "#CC79A7", male = "#0072B2")) +
  labs(
    x = "Normal-labelled posteriorHCA latent log(mu), offset 0",
    y = NULL,
    colour = "Sex",
    title = "TM4SF1 across observed tissue × sex profile subsets"
  ) +
  theme_report()

save_plot(sex_plot, "07_posteriorHCA_TM4SF1_tissue_by_sex.png", 11, 8)
sex_plot
```

![Descriptive observed tissue × sex profile averages. Sex is female/male because those are the levels available in the fitted models; these are not controlled sex effects.](/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/figures/tissue-sex-plot-1.png)

### Observed age-decade profiles within tissue


``` r
top_tissues <- profile_summary |>
  filter(dimension == "tissue", supported) |>
  group_by(cell_type) |>
  slice_max(median_log_mu, n = 8, with_ties = FALSE) |>
  ungroup() |>
  select(cell_type, tissue_groups)

candidate_tissues <- profile_summary |>
  filter(
    prototype_overlap_candidate,
    !is.na(tissue_groups)
  ) |>
  select(cell_type, tissue_groups)

display_tissues <- bind_rows(top_tissues, candidate_tissues) |>
  distinct()

age_plot_data <- profile_summary |>
  filter(dimension == "tissue_age", supported) |>
  inner_join(display_tissues, by = c("cell_type", "tissue_groups")) |>
  mutate(
    age_decade = factor(
      age_decade,
      levels = as.character(sort(unique(as.integer(age_decade))))
    ),
    tissue_groups = str_wrap(tissue_groups, 24)
  )

age_plot <- ggplot(
  age_plot_data,
  aes(age_decade, tissue_groups, fill = probability_above_threshold)
) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_point(
    data = filter(age_plot_data, prototype_overlap_candidate),
    shape = 21, fill = NA, colour = "black", size = 2.6, stroke = 0.8
  ) +
  facet_wrap(~ cell_type, scales = "free_y") +
  scale_fill_viridis_c(
    limits = c(0, 1), labels = label_percent(), option = "C"
  ) +
  labs(
    x = "Age decade (e.g. 6 = 60-69 years)",
    y = NULL,
    fill = "P(Normal profile >= anchor)",
    title = "TM4SF1 across observed tissue × age-decade profile subsets"
  ) +
  theme_report()

save_plot(age_plot, "08_posteriorHCA_TM4SF1_tissue_by_age.png", 11, 7)
age_plot
```

![Posterior probability that a supported Normal-labelled tissue × age-decade profile exceeds the model-specific HV tumour-expression anchor. The eight highest aggregate tissues plus every qualifying prototype tissue are shown; these descriptive averages are not controlled age effects.](/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/figures/tissue-age-plot-1.png)

### Observed ethnicity profiles within tissue


``` r
ethnicity_plot_data <- profile_summary |>
  filter(dimension == "tissue_ethnicity", supported) |>
  inner_join(display_tissues, by = c("cell_type", "tissue_groups")) |>
  mutate(
    tissue_groups = str_wrap(tissue_groups, 24),
    ethnicity_groups = str_wrap(ethnicity_groups, 20)
  )

ethnicity_plot <- ggplot(
  ethnicity_plot_data,
  aes(ethnicity_groups, tissue_groups, fill = probability_above_threshold)
) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_point(
    data = filter(ethnicity_plot_data, prototype_overlap_candidate),
    shape = 21, fill = NA, colour = "black", size = 2.6, stroke = 0.8
  ) +
  facet_wrap(~ cell_type, scales = "free_y") +
  scale_fill_viridis_c(
    limits = c(0, 1), labels = label_percent(), option = "C"
  ) +
  labs(
    x = "Recorded/harmonised ethnicity group",
    y = NULL,
    fill = "P(Normal profile >= anchor)",
    title = "TM4SF1 across observed tissue × ethnicity profile subsets"
  ) +
  theme_report() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_plot(
  ethnicity_plot,
  "09_posteriorHCA_TM4SF1_tissue_by_ethnicity.png",
  12, 8
)
ethnicity_plot
```

![Descriptive observed tissue × recorded/harmonised ethnicity profile averages for the eight highest aggregate tissues plus every qualifying prototype tissue. Sparse strata cannot qualify as candidates; these are not controlled ethnicity effects.](/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/figures/tissue-ethnicity-plot-1.png)

### Ungated prototype candidates for follow-up


``` r
candidate_table <- profile_summary |>
  filter(
    prototype_overlap_candidate,
    dimension %in% c("tissue", "tissue_sex", "tissue_age", "tissue_ethnicity")
  ) |>
  arrange(cell_type, dimension, desc(probability_above_threshold)) |>
  transmute(
    cell_type,
    dimension,
    tissue = tissue_groups,
    sex,
    age_decade,
    ethnicity = ethnicity_groups,
    donor_keys = n_donor_keys,
    median_log_mu = round(median_log_mu, 2),
    log2_vs_HV = round(estimate_log2_vs_hv, 2),
    P_above_threshold = percent(probability_above_threshold, accuracy = 0.1),
    P_at_least_random_HV = percent(
      probability_hca_ge_random_hv, accuracy = 0.1
    )
  )

head(candidate_table, 40) |>
  kable(
    caption = paste(
      "First 40 supported prototype overlap candidates, ungated by model",
      "readiness.",
      "The full machine-readable table is written to results/."
    )
  )
```



Table: First 40 supported prototype overlap candidates, ungated by model readiness. The full machine-readable table is written to results/.

|cell_type   |dimension  |tissue                      |sex    |age_decade |ethnicity | donor_keys| median_log_mu| log2_vs_HV|P_above_threshold |P_at_least_random_HV |
|:-----------|:----------|:---------------------------|:------|:----------|:---------|----------:|-------------:|----------:|:-----------------|:--------------------|
|endothelial |tissue     |female reproductive system  |NA     |NA         |NA        |         10|          8.19|       2.12|91.5%             |81.0%                |
|endothelial |tissue     |integumentary system (skin) |NA     |NA         |NA        |         13|          7.86|       1.42|86.5%             |72.9%                |
|endothelial |tissue_sex |female reproductive system  |female |NA         |NA        |         10|          8.19|       2.11|91.5%             |81.0%                |
|endothelial |tissue_sex |integumentary system (skin) |male   |NA         |NA        |         11|          7.84|       1.32|86.8%             |73.0%                |

The tissue, tissue × sex, tissue × age and tissue × ethnicity rows are nested views of the same atlas profiles; they should not be counted as independent safety findings.

The motivating respiratory-system example is checked explicitly rather than assumed:


``` r
respiratory_sex <- profile_summary |>
  filter(
    dimension == "tissue_sex",
    tissue_groups == "respiratory system",
    supported
  ) |>
  transmute(
    cell_type,
    sex,
    donor_keys = n_donor_keys,
    median_log_mu = round(median_log_mu, 2),
    threshold_log_mu = round(threshold_log_mu, 2),
    P_above_threshold = percent(probability_above_threshold, accuracy = 0.1),
    prototype_candidate = prototype_overlap_candidate,
    model_ready,
    model_ready_alert = expression_overlap_alert
  )

respiratory_sex |>
  kable(
    caption = "Supported respiratory-system sex profiles under the current rule."
  )
```



Table: Supported respiratory-system sex profiles under the current rule.

|cell_type   |sex    | donor_keys| median_log_mu| threshold_log_mu|P_above_threshold |prototype_candidate |model_ready |model_ready_alert |
|:-----------|:------|----------:|-------------:|----------------:|:-----------------|:-------------------|:-----------|:-----------------|
|epithelial  |female |         47|          7.78|             8.37|41.8%             |FALSE               |FALSE       |FALSE             |
|epithelial  |male   |         51|          7.72|             8.37|40.0%             |FALSE               |FALSE       |FALSE             |
|endothelial |female |         40|          6.46|             6.07|61.2%             |FALSE               |FALSE       |FALSE             |
|endothelial |male   |         45|          6.32|             6.07|57.5%             |FALSE               |FALSE       |FALSE             |

No supported respiratory-system sex profile passes the current 80% prototype-overlap rule. Model readiness is reported separately, so the present analysis supports neither a respiratory-system candidate nor a male-lung alert.

## 7. Conclusions and next decisions

The released DE table ranks TM4SF1 first among HV-enriched genes, and the closest public-DGE reconstruction supports HV enrichment, subject to the documented 9-HV/2-UC limitation. Its tumour expression is heterogeneous, so a target-positive HV label should not be interpreted as uniform antigen abundance.

The off-tumour screening conclusion is cell-type specific. Endothelium is a necessary sentinel because the bladder dataset itself shows high endothelial CPM and near-universal detection. At the registered default settings, the ungated tissue-level candidates are endothelial: female reproductive system, integumentary system (skin). These are follow-up candidates, not toxicity calls or model-ready alerts; readiness is reported separately above.

A supported **respiratory-system × male** profile could be prioritised for off-tumour validation only if it passed the registered probability rule, the corresponding female profile were examined, and a refitted posteriorHCA model confirmed the result. The current prototype does not pass that rule and cannot infer vulnerability. The fitted tissue group combines lung and airway compartments and cannot localise risk to a specific pulmonary structure.

The immediate next analysis should be a longer-chain refit of TM4SF1 epithelial and endothelial models, followed by protein-level localisation and antigen-density calibration in the prototype candidate tissues before any toxicity interpretation. This report intentionally stops before defining a clinical toxicity threshold.

## Reproducible outputs


``` r
output_manifest <- tribble(
  ~file, ~contents,
  "TM4SF1_paper_DE_replication.tsv", "Released and public-DGE TM4SF1 DE results",
  "TM4SF1_paper_source_patient_summary.tsv", "Figure5B patient-level source-expression summary",
  "TM4SF1_study_logmu_draws.rds", "Patient Bayesian-bootstrap cohort log(mu) draws",
  "TM4SF1_cohort_logmu_contrasts.tsv", "Study-versus-HCA log2 contrasts",
  "TM4SF1_threshold_registry.tsv", "Cell-type-specific threshold contract and provenance",
  "TM4SF1_alert_rule_sensitivity.tsv", "Ungated prototype tissue-candidate counts while varying anchor quantile, donor support and exceedance probability",
  "posteriorHCA_TM4SF1_profile_summary.tsv", "All tissue/demographic posterior summaries, including unsupported rows",
  "posteriorHCA_TM4SF1_supported_profiles.tsv", "Profiles meeting the minimum donor-key rule",
  "TM4SF1_prototype_overlap_candidates.tsv", "Supported ungated profiles passing the prototype overlap rule",
  "TM4SF1_expression_overlap_alerts.tsv", "Readiness-gated alerts; empty whenever model_ready is false",
  "TM4SF1_model_diagnostics.tsv", "Current brms fit convergence diagnostics",
  "TM4SF1_reference_alignment.tsv", "TMM reference, shared feature count and TM4SF1 dispersion",
  "TM4SF1_toxicity_case_study.md", "Knitted Markdown preview of the executed QMD",
  "output_manifest.tsv", "Machine-readable inventory of report outputs",
  "sessionInfo.txt", "R session and package versions for this report",
  "figures/", "Report figures as reusable PNG files",
  "cache/", "Per-cell-type posterior/cache objects keyed by analysis version"
) |>
  mutate(path = file.path(results_dir, file), .before = 1)

write_tsv(output_manifest, file.path(results_dir, "output_manifest.tsv"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo.txt"))
```

```
## Warning in system("timedatectl", intern = TRUE): running command 'timedatectl'
## had status 1
```

``` r
output_manifest |>
  kable(caption = "Files produced by this report.")
```



Table: Files produced by this report.

|path                                                                                                                                                                |file                                       |contents                                                                                                          |
|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------|:-----------------------------------------------------------------------------------------------------------------|
|/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/TM4SF1_paper_DE_replication.tsv            |TM4SF1_paper_DE_replication.tsv            |Released and public-DGE TM4SF1 DE results                                                                         |
|/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/TM4SF1_paper_source_patient_summary.tsv    |TM4SF1_paper_source_patient_summary.tsv    |Figure5B patient-level source-expression summary                                                                  |
|/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/TM4SF1_study_logmu_draws.rds               |TM4SF1_study_logmu_draws.rds               |Patient Bayesian-bootstrap cohort log(mu) draws                                                                   |
|/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/TM4SF1_cohort_logmu_contrasts.tsv          |TM4SF1_cohort_logmu_contrasts.tsv          |Study-versus-HCA log2 contrasts                                                                                   |
|/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/TM4SF1_threshold_registry.tsv              |TM4SF1_threshold_registry.tsv              |Cell-type-specific threshold contract and provenance                                                              |
|/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/TM4SF1_alert_rule_sensitivity.tsv          |TM4SF1_alert_rule_sensitivity.tsv          |Ungated prototype tissue-candidate counts while varying anchor quantile, donor support and exceedance probability |
|/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/posteriorHCA_TM4SF1_profile_summary.tsv    |posteriorHCA_TM4SF1_profile_summary.tsv    |All tissue/demographic posterior summaries, including unsupported rows                                            |
|/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/posteriorHCA_TM4SF1_supported_profiles.tsv |posteriorHCA_TM4SF1_supported_profiles.tsv |Profiles meeting the minimum donor-key rule                                                                       |
|/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/TM4SF1_prototype_overlap_candidates.tsv    |TM4SF1_prototype_overlap_candidates.tsv    |Supported ungated profiles passing the prototype overlap rule                                                     |
|/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/TM4SF1_expression_overlap_alerts.tsv       |TM4SF1_expression_overlap_alerts.tsv       |Readiness-gated alerts; empty whenever model_ready is false                                                       |
|/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/TM4SF1_model_diagnostics.tsv               |TM4SF1_model_diagnostics.tsv               |Current brms fit convergence diagnostics                                                                          |
|/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/TM4SF1_reference_alignment.tsv             |TM4SF1_reference_alignment.tsv             |TMM reference, shared feature count and TM4SF1 dispersion                                                         |
|/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/TM4SF1_toxicity_case_study.md              |TM4SF1_toxicity_case_study.md              |Knitted Markdown preview of the executed QMD                                                                      |
|/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/output_manifest.tsv                        |output_manifest.tsv                        |Machine-readable inventory of report outputs                                                                      |
|/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/sessionInfo.txt                            |sessionInfo.txt                            |R session and package versions for this report                                                                    |
|/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/figures/                                   |figures/                                   |Report figures as reusable PNG files                                                                              |
|/scratchdata1/groups/phoenix-hpc-mangiola_laboratory/chen/posteriorHCA/dev/TM4SF1_bladder/results/tm4sf1_toxicity_report/cache/                                     |cache/                                     |Per-cell-type posterior/cache objects keyed by analysis version                                                   |

Primary provenance: [GEO GSE293189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293189), the [published study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12174346/), and the authors' [public analysis repository](https://github.com/angelussong/Histological_Variant_Bladder_Cancer_Analysis). The preprocessing discrepancy and exact 67-library mapping are documented in the case-study `results/replication_audit/` and `data/metadata/` directories.
