# Shared helpers: RE tissue_groups vs fixed interaction pilot (V1_nk, 4 TIMES genes).
#
# This revision refits BOTH model families (random-effect and interaction) on the
# login node with a *revised shape formula* and, for each family, also fits the
# previous *rich* shape so we can test whether the shape formula controls the
# posterior SD of the comprehensive liver baseline.
#
#   Counts formulas (unchanged):
#     re          : random slopes inside tissue_groups
#     interaction : fixed tissue_groups x covariate interactions
#
#   Shape formulas:
#     rich   (previous) : includes ethnicity_groups (+ tissue term)
#     simple (revised)  : drops ethnicity/tissue-specific terms so the shape
#                          borrows power across tissues instead of being limited
#                          by small per-tissue sample sizes (e.g. liver n = 71).

formula_specs <- function() {
  re_counts <- "counts ~ 1 + offset(offset) + age_decade * sex + disease_groups_altered + ethnicity_groups + assay_groups_altered + (1 | dataset_id_altered) + (1 + age_decade * sex + ethnicity_groups | tissue_groups)"
  # Same RE mean structure with uncorrelated tissue slopes (brms `||`).
  re_uncorrelated_counts <- "counts ~ 1 + offset(offset) + age_decade * sex + disease_groups_altered + ethnicity_groups + assay_groups_altered + (1 | dataset_id_altered) + (1 + age_decade * sex + ethnicity_groups || tissue_groups)"
  # Split ethnicity slopes out of the correlated tissue RE block (Pathfinder-friendly).
  re_split_counts <- paste(
    "counts ~ 1 + offset(offset) + age_decade * sex + disease_groups_altered",
    "+ ethnicity_groups + assay_groups_altered + (1 | dataset_id_altered)",
    "+ (1 + age_decade * sex | tissue_groups)",
    "+ (0 + ethnicity_groups | tissue_groups)"
  )
  int_counts <- "counts ~ 1 + offset(offset) + disease_groups_altered + assay_groups_altered + (1 | dataset_id_altered) + tissue_groups + age_decade:tissue_groups + sex:tissue_groups + age_decade:sex:tissue_groups + ethnicity_groups:tissue_groups"
  shape_simple <- "shape ~ 1 + disease_groups_altered + assay_groups_altered + (1 | tissue_groups)"

  list(
    re = list(
      counts = re_counts,
      shape_rich = "shape ~ 1 + disease_groups_altered + assay_groups_altered + ethnicity_groups + (1 | tissue_groups)",
      shape_simple = shape_simple
    ),
    re_uncorrelated = list(
      counts = re_uncorrelated_counts,
      shape_simple = shape_simple
    ),
    re_split = list(
      counts = re_split_counts,
      shape_simple = shape_simple
    ),
    interaction = list(
      counts = int_counts,
      shape_rich = "shape ~ 1 + disease_groups_altered + assay_groups_altered + ethnicity_groups + tissue_groups",
      shape_simple = "shape ~ 1 + disease_groups_altered + assay_groups_altered"
    )
  )
}

# Backward-compatible accessors: report the *revised* (simple) production candidate.
formula_re_spec <- function() {
  s <- formula_specs()$re
  list(counts = s$counts, shape = s$shape_simple)
}

formula_interaction_spec <- function() {
  s <- formula_specs()$interaction
  list(counts = s$counts, shape = s$shape_simple)
}

fit_levels <- function(fit, var) {
  d <- fit$data
  if (is.null(d) || !var %in% names(d)) return(character(0))
  v <- d[[var]]
  lv <- if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
  lv[!is.na(lv) & nzchar(lv)]
}

build_comprehensive_grid <- function(
    fit,
    tissue_group = "liver",
    marginalize_over = c("sex", "age_decade", "ethnicity_groups", "assay_groups_altered"),
    new_level_cols = "dataset_id_altered",
    new_level_tag = "__comprehensive_healthy_baseline__"
) {
  newdata <- tibble::tibble(
    age_decade = NA_character_,
    sex = NA_character_,
    ethnicity_groups = NA_character_,
    assay_groups_altered = NA_character_,
    disease_groups_altered = "Normal",
    tissue_groups = tissue_group,
    offset = 0
  )

  fit_vars <- colnames(fit$data)
  levels_list <- lapply(intersect(marginalize_over, fit_vars), function(col) {
    fit_levels(fit, col)
  })
  names(levels_list) <- intersect(marginalize_over, fit_vars)
  levels_list <- levels_list[vapply(levels_list, length, integer(1L)) > 0L]

  if (length(levels_list)) {
    grid <- do.call(tidyr::expand_grid, levels_list)
    rest <- newdata[, setdiff(names(newdata), names(levels_list)), drop = FALSE]
    out <- dplyr::bind_cols(rest[rep(1L, nrow(grid)), , drop = FALSE], grid)
  } else {
    out <- newdata
  }

  for (col in new_level_cols) {
    if (col %in% fit_vars) out[[col]] <- new_level_tag
  }

  for (col in intersect(names(out), fit_vars)) {
    if (is.factor(fit$data[[col]])) {
      out[[col]] <- factor(out[[col]], levels = levels(fit$data[[col]]))
    }
  }

  tibble::as_tibble(out)
}

linpred_baseline_draws <- function(fit, ndraws = NULL, seed = 20260722L) {
  nd_grid <- build_comprehensive_grid(fit)
  if (!is.null(seed)) set.seed(as.integer(seed))
  log_mu_mat <- brms::posterior_linpred(
    fit,
    newdata = nd_grid,
    transform = FALSE,
    re_formula = NULL,
    allow_new_levels = TRUE,
    sample_new_levels = "gaussian",
    ndraws = ndraws
  )
  tibble::tibble(
    draw = seq_len(nrow(log_mu_mat)),
    baseline_log_mu = rowMeans(log_mu_mat),
    n_grid_rows = nrow(nd_grid)
  )
}

build_bf <- function(counts_formula, shape_formula) {
  brms::bf(stats::as.formula(counts_formula), stats::as.formula(shape_formula))
}

# Generic gene-level fitter for any (counts, shape) pair.
# `cores` = parallel chains; `threads_per_chain` = within-chain reduce_sum threads.
# Request Slurm cpus_per_task >= chains * threads_per_chain.
fit_gene_model <- function(
    data, gene, counts_formula, shape_formula,
    cores = 32L, chains = 2L, warmup = 400L, iter = 600L,
    threads_per_chain = NULL,
    adapt_delta = 0.95, max_treedepth = 12L, seed = 42L, tag = "model"
) {
  data <- droplevels(data)
  colnames(data) <- stringr::str_replace_all(colnames(data), "_+", "_")

  i <- mean(log1p(data$counts / exp(data$offset)))
  prior <- eval(substitute(
    c(
      prior(student_t(3, i, 1.5), class = Intercept),
      prior(student_t(3, 0, 1), class = Intercept, dpar = shape),
      prior(student_t(3, 0, 5), class = b),
      prior(student_t(3, 0, 2), class = b, dpar = shape)
    ),
    env = list(i = i)
  ))

  n_chains <- as.integer(chains)
  n_cores_chains <- min(as.integer(cores), n_chains)
  if (is.null(threads_per_chain)) {
    avail <- suppressWarnings(as.integer(parallelly::availableCores()))
    if (!is.finite(avail) || avail < 1L) avail <- n_cores_chains
    threads_per_chain <- max(1L, as.integer(floor(avail / n_cores_chains)))
  }
  threads_per_chain <- as.integer(threads_per_chain)
  threads_arg <- if (threads_per_chain > 1L) {
    brms::threading(threads = threads_per_chain)
  } else {
    NULL
  }

  message(glue::glue(
    "Fitting {tag} model for {gene} (chain_cores = {n_cores_chains}, chains = {n_chains}, ",
    "threads_per_chain = {threads_per_chain}, warmup = {warmup}, iter = {iter}, ",
    "seed = {seed}, adapt_delta = {adapt_delta}, max_treedepth = {max_treedepth}) ..."
  ))
  brms::brm(
    formula = build_bf(counts_formula, shape_formula),
    data = data,
    family = brms::zero_inflated_negbinomial(),
    prior = prior,
    chains = n_chains,
    cores = n_cores_chains,
    threads = threads_arg,
    warmup = warmup,
    iter = iter,
    refresh = 50,
    backend = "cmdstanr",
    silent = 1,
    seed = as.integer(seed),
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  )
}

# The 2 (counts family) x 2 (shape complexity) fit grid.
# `role` distinguishes the revised production candidate (simple shape) from the
# rich-shape reference used only to isolate the shape effect on baseline SD.
fit_grid_spec <- function() {
  specs <- formula_specs()
  tibble::tribble(
    ~model,        ~shape_variant, ~role,        ~counts,                  ~shape,                        ~file_pattern,
    "re",          "simple",       "revised",    specs$re$counts,          specs$re$shape_simple,          "fit_re_{gene}.rds",
    "re",          "rich",         "reference",  specs$re$counts,          specs$re$shape_rich,            "fit_re_rich_{gene}.rds",
    "interaction", "simple",       "revised",    specs$interaction$counts, specs$interaction$shape_simple, "fit_interaction_{gene}.rds",
    "interaction", "rich",         "reference",  specs$interaction$counts, specs$interaction$shape_rich,   "fit_interaction_rich_{gene}.rds"
  )
}

# Stronger MCMC comparison grid: all three mean formulas x simple shape only.
# Used by the targets pipeline (4 chains / 1500 warmup / 3000 iter).
fit_grid_spec_strong <- function() {
  specs <- formula_specs()
  tibble::tribble(
    ~model,             ~shape_variant, ~counts,                            ~shape,                              ~file_pattern,
    "re",               "simple",       specs$re$counts,                    specs$re$shape_simple,               "fit_re_{gene}.rds",
    "interaction",      "simple",       specs$interaction$counts,           specs$interaction$shape_simple,      "fit_interaction_{gene}.rds",
    "re_uncorrelated",  "simple",       specs$re_uncorrelated$counts,       specs$re_uncorrelated$shape_simple,  "fit_re_uncorrelated_{gene}.rds"
  )
}

# Pathfinder grid: original RE + split RE, simple shape, CPU dense.
# Draw budget matched to strong HMC: (iter - warmup) * chains = 1500 * 4 = 6000.
# Pathfinder total draws ≈ num_paths * single_path_draws (psis_resample = FALSE).
fit_grid_spec_pathfinder <- function() {
  specs <- formula_specs()
  tibble::tribble(
    ~model,              ~shape_variant, ~counts,                     ~shape,                       ~file_pattern,
    "pathfinder_re",     "simple",       specs$re$counts,             specs$re$shape_simple,        "fit_pathfinder_re_{gene}.rds",
    "pathfinder_re_split","simple",      specs$re_split$counts,       specs$re_split$shape_simple,  "fit_pathfinder_re_split_{gene}.rds"
  )
}

# Match strong HMC post-warmup draw count.
pathfinder_draw_settings <- function(
    hmc_chains = 4L,
    hmc_warmup = 1500L,
    hmc_iter = 3000L,
    num_paths = 50L
) {
  target_draws <- as.integer(hmc_chains * (hmc_iter - hmc_warmup))
  num_paths <- as.integer(num_paths)
  single_path_draws <- as.integer(ceiling(target_draws / num_paths))
  list(
    target_draws = target_draws,
    num_paths = num_paths,
    single_path_draws = single_path_draws,
    total_draws = num_paths * single_path_draws
  )
}

# Pathfinder fitter: CPU + dense (not sparse, not OpenCL). Successful pilot pattern.
fit_gene_model_pathfinder <- function(
    data, gene, counts_formula, shape_formula,
    num_paths = 50L,
    single_path_draws = 120L,
    threads = 8L,
    history_size = 100L,
    max_lbfgs_iters = 100L,
    seed = 20260720L,
    tag = "pathfinder"
) {
  data <- droplevels(data)
  colnames(data) <- stringr::str_replace_all(colnames(data), "_+", "_")

  i <- mean(log1p(data$counts / exp(data$offset)))
  prior <- eval(substitute(
    c(
      prior(student_t(3, i, 1.5), class = Intercept),
      prior(student_t(3, 0, 1), class = Intercept, dpar = shape),
      prior(student_t(3, 0, 5), class = b),
      prior(student_t(3, 0, 2), class = b, dpar = shape)
    ),
    env = list(i = i)
  ))

  form <- brms::bf(
    stats::as.formula(counts_formula),
    stats::as.formula(shape_formula),
    sparse = FALSE
  )

  # Force a CPU Stan executable (do not reuse an OpenCL-compiled binary).
  force_cpu_stanvars <- brms::stanvar(
    scode = "// force_cpu_recompile_pathfinder_20260722",
    block = "functions"
  )

  threads <- as.integer(threads)
  threads_arg <- if (threads > 1L) brms::threading(threads = threads) else NULL

  message(glue::glue(
    "Fitting {tag} Pathfinder for {gene}: num_paths={num_paths}, ",
    "single_path_draws={single_path_draws}, total≈{as.integer(num_paths) * as.integer(single_path_draws)}, ",
    "threads={threads}, seed={seed} ..."
  ))

  brms::brm(
    formula = form,
    data = data,
    family = brms::zero_inflated_negbinomial(),
    prior = prior,
    stanvars = force_cpu_stanvars,
    backend = "cmdstanr",
    algorithm = "pathfinder",
    opencl = NULL,
    chains = as.integer(num_paths),
    cores = 1L,
    threads = threads_arg,
    history_size = as.integer(history_size),
    max_lbfgs_iters = as.integer(max_lbfgs_iters),
    single_path_draws = as.integer(single_path_draws),
    psis_resample = FALSE,
    seed = as.integer(seed),
    refresh = 50,
    silent = 1
  )
}

# Default TIMES biomarkers present in the V1_nk production store.
pilot_biomarker_genes <- function() {
  c("SPON2", "ZFP36L2", "ZFP36", "VIM")
}

# Timestamped progress line, flushed immediately so it streams to logs / console.
pilot_log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  message(msg)
  utils::flush.console()
  invisible(msg)
}

fmt_dur <- function(secs) {
  secs <- as.numeric(secs)
  if (secs < 90) return(sprintf("%.0fs", secs))
  if (secs < 5400) return(sprintf("%.1fmin", secs / 60))
  sprintf("%.1fh", secs / 3600)
}

run_re_vs_interaction_pilot <- function(
    case_dir,
    out_dir,
    model_store_root = "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE",
    target_cell_type = "nk",
    baseline_tissue = "liver",
    n_draws = 400L,
    refresh_fits = FALSE,
    fit_cores = 32L,
    genes = NULL,
    models = NULL,
    shape_variants = NULL,
    chains = 2L,
    warmup = 400L,
    iter = 600L,
    adapt_delta = 0.95,
    max_treedepth = 12L,
    seed = 42L
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fit_cores <- as.integer(min(fit_cores, parallel::detectCores()))

  store <- file.path(model_store_root, paste0("V1_", target_cell_type), "_targets")

  model_meta <- targets::tar_meta(starts_with("estimates_chunk_"), store = store) |>
    dplyr::filter(.data$type == "branch") |>
    dplyr::mutate(.feature = stringr::str_extract(.data$warnings, "(?<=Gene:___)ENSG\\d+(?=___)"))

  biomarkers <- readr::read_csv(
    file.path(case_dir, "data", "times_biomarkers.csv"),
    show_col_types = FALSE
  ) |>
    dplyr::mutate(in_store = .data$.feature %in% model_meta$.feature) |>
    dplyr::filter(.data$in_store)

  if (!is.null(genes)) {
    biomarkers <- biomarkers |> dplyr::filter(.data$gene %in% genes)
  }

  read_fit_re <- function(feature_id) {
    target_name <- model_meta |>
      dplyr::filter(.data$.feature == feature_id) |>
      dplyr::slice(1) |>
      dplyr::pull(.data$name)
    targets::tar_read_raw(target_name[[1]], store = store)$brms_fit[[1]]
  }

  se <- targets::tar_read(pseudobulk_sample, store = store)
  cd <- tibble::as_tibble(SummarizedExperiment::colData(se))

  readr::write_csv(
    tibble::tibble(
      tissue_groups = names(sort(table(cd$tissue_groups), decreasing = TRUE)),
      n = as.integer(sort(table(cd$tissue_groups), decreasing = TRUE))
    ),
    file.path(out_dir, "nk_pseudobulk_tissue_sample_sizes.csv")
  )

  example_fit <- read_fit_re(biomarkers$.feature[1])
  liver_rows <- example_fit$data |>
    dplyr::filter(.data$tissue_groups == baseline_tissue, .data$disease_groups_altered == "Normal")

  readr::write_csv(
    tibble::tibble(
      metric = c(
        "nk_pseudobulk_total",
        "nk_liver_all_disease",
        "nk_liver_normal_in_fit",
        "n_tissue_levels",
        "n_age_decade_levels",
        "n_ethnicity_levels"
      ),
      value = c(
        nrow(cd),
        sum(cd$tissue_groups == baseline_tissue),
        nrow(liver_rows),
        dplyr::n_distinct(cd$tissue_groups),
        dplyr::n_distinct(example_fit$data$age_decade),
        dplyr::n_distinct(example_fit$data$ethnicity_groups)
      )
    ),
    file.path(out_dir, "liver_sample_size_summary.csv")
  )

  grid <- fit_grid_spec()
  if (!is.null(models)) {
    grid <- grid |> dplyr::filter(.data$model %in% models)
  }
  if (!is.null(shape_variants)) {
    grid <- grid |> dplyr::filter(.data$shape_variant %in% shape_variants)
  }
  if (!nrow(grid)) {
    stop("No fit-grid cells left after models/shape_variants filters.")
  }

  draws_long <- list()
  summary_rows <- list()

  n_genes <- nrow(biomarkers)
  n_cells <- n_genes * nrow(grid)
  cell_i <- 0L
  t_start <- Sys.time()
  pilot_log(glue::glue(
    "Starting pilot: {n_genes} genes x {nrow(grid)} model cells = {n_cells} baselines ",
    "(fit_cores = {fit_cores}, refresh_fits = {isTRUE(refresh_fits)}, ",
    "chains = {chains}, warmup = {warmup}, iter = {iter}, seed = {seed})."
  ))

  for (i in seq_len(n_genes)) {
    gene <- biomarkers$gene[i]
    feature_id <- biomarkers$.feature[i]

    pilot_log(glue::glue("=== Gene {i}/{n_genes}: {gene} — loading production RE fit for data ..."))
    fit_re_prod <- read_fit_re(feature_id)
    data <- droplevels(fit_re_prod$data)
    n_liver_rows <- sum(data$tissue_groups == baseline_tissue)

    for (r in seq_len(nrow(grid))) {
      spec <- grid[r, ]
      cell_i <- cell_i + 1L
      cell_tag <- glue::glue("{spec$model}/{spec$shape_variant}")
      fit_path <- file.path(out_dir, glue::glue(spec$file_pattern, gene = gene))
      t_cell <- Sys.time()

      # The production RE fit already uses the rich RE shape, so it *is* the
      # RE rich-shape reference; reuse it instead of refitting the expensive
      # random-slopes model.
      is_re_rich_reference <- spec$model == "re" && spec$shape_variant == "rich"

      if (is_re_rich_reference) {
        pilot_log(glue::glue("[cell {cell_i}/{n_cells}] {gene} {cell_tag}: reusing production store fit (reference)."))
        fit <- fit_re_prod
        fit_source <- "production_store"
      } else if (file.exists(fit_path) && !isTRUE(refresh_fits)) {
        pilot_log(glue::glue("[cell {cell_i}/{n_cells}] {gene} {cell_tag}: loading cached fit."))
        fit <- readRDS(fit_path)
        fit_source <- "cache"
      } else {
        pilot_log(glue::glue("[cell {cell_i}/{n_cells}] {gene} {cell_tag}: fitting brms (this can take a while) ..."))
        fit <- fit_gene_model(
          data, gene,
          counts_formula = spec$counts,
          shape_formula = spec$shape,
          cores = fit_cores,
          chains = chains,
          warmup = warmup,
          iter = iter,
          adapt_delta = adapt_delta,
          max_treedepth = max_treedepth,
          seed = seed,
          tag = cell_tag
        )
        saveRDS(fit, fit_path, compress = "xz")
        fit_source <- "fit"
        pilot_log(glue::glue(
          "[cell {cell_i}/{n_cells}] {gene} {cell_tag}: fit done in ",
          "{fmt_dur(difftime(Sys.time(), t_cell, units = 'secs'))}; saved {basename(fit_path)}."
        ))
      }

      pilot_log(glue::glue("[cell {cell_i}/{n_cells}] {gene} {cell_tag}: computing baseline draws ..."))
      d <- linpred_baseline_draws(fit, ndraws = n_draws) |>
        dplyr::mutate(
          gene = gene,
          .feature = feature_id,
          model = spec$model,
          shape_variant = spec$shape_variant,
          role = spec$role
        )
      draws_long[[paste(gene, spec$model, spec$shape_variant, sep = "_")]] <- d

      cell_sd <- stats::sd(d$baseline_log_mu)
      summary_rows[[paste(gene, spec$model, spec$shape_variant, sep = "_")]] <- tibble::tibble(
        gene = gene,
        .feature = feature_id,
        model = spec$model,
        shape_variant = spec$shape_variant,
        role = spec$role,
        fit_source = fit_source,
        n_draws = nrow(d),
        mean_log_mu = mean(d$baseline_log_mu),
        sd_log_mu = cell_sd,
        median_mu = stats::median(exp(d$baseline_log_mu)),
        max_rhat = suppressWarnings(max(posterior::rhat(fit), na.rm = TRUE)),
        n_data_rows = nrow(data),
        n_liver_rows = n_liver_rows
      )

      elapsed <- difftime(Sys.time(), t_start, units = "secs")
      eta <- as.numeric(elapsed) / cell_i * (n_cells - cell_i)
      pilot_log(glue::glue(
        "[cell {cell_i}/{n_cells}] {gene} {cell_tag}: done (baseline SD = {round(cell_sd, 3)}). ",
        "Elapsed {fmt_dur(elapsed)}, ETA {fmt_dur(eta)}."
      ))
    }
  }

  pilot_log(glue::glue("All {n_cells} baselines computed in {fmt_dur(difftime(Sys.time(), t_start, units = 'secs'))}. Writing outputs ..."))

  draws_df <- dplyr::bind_rows(draws_long)
  baseline_summary <- dplyr::bind_rows(summary_rows)

  # If this run was a gene/model subset, merge into existing full artefacts
  # so other genes/cells are not wiped.
  comparison_path <- file.path(out_dir, "baseline_comparison_by_formula.csv")
  draws_path <- file.path(out_dir, "baseline_draws_long.csv")
  if (file.exists(comparison_path) && (!is.null(genes) || !is.null(models) || !is.null(shape_variants))) {
    old_sum <- readr::read_csv(comparison_path, show_col_types = FALSE)
    key <- c("gene", "model", "shape_variant")
    old_sum <- old_sum |>
      dplyr::anti_join(baseline_summary |> dplyr::select(dplyr::all_of(key)), by = key)
    baseline_summary <- dplyr::bind_rows(old_sum, baseline_summary) |>
      dplyr::arrange(.data$gene, .data$model, .data$shape_variant)
  }
  if (file.exists(draws_path) && (!is.null(genes) || !is.null(models) || !is.null(shape_variants))) {
    old_draws <- readr::read_csv(draws_path, show_col_types = FALSE)
    key <- c("gene", "model", "shape_variant")
    old_draws <- old_draws |>
      dplyr::anti_join(draws_df |> dplyr::distinct(dplyr::across(dplyr::all_of(key))), by = key)
    draws_df <- dplyr::bind_rows(old_draws, draws_df)
  }

  # Isolate the shape effect: within each (gene, model), how does the simple
  # shape change the baseline mean and SD relative to the rich shape?
  shape_effect <- baseline_summary |>
    dplyr::select(gene, model, shape_variant, mean_log_mu, sd_log_mu, median_mu) |>
    tidyr::pivot_wider(
      names_from = shape_variant,
      values_from = c(mean_log_mu, sd_log_mu, median_mu)
    ) |>
    dplyr::mutate(
      sd_ratio_simple_over_rich = .data$sd_log_mu_simple / .data$sd_log_mu_rich,
      sd_reduction = .data$sd_log_mu_rich - .data$sd_log_mu_simple,
      mean_shift_simple_minus_rich = .data$mean_log_mu_simple - .data$mean_log_mu_rich
    )

  readr::write_csv(baseline_summary, file.path(out_dir, "baseline_comparison_by_formula.csv"))
  readr::write_csv(shape_effect, file.path(out_dir, "baseline_shape_effect.csv"))
  readr::write_csv(draws_df, file.path(out_dir, "baseline_draws_long.csv"))
  saveRDS(draws_df, file.path(out_dir, "baseline_draws_re_vs_interaction.rds"), compress = "xz")

  specs <- formula_specs()
  readr::write_lines(
    c(
      "# Random-effect counts formula",
      specs$re$counts,
      "# RE shape (rich reference / production)",
      specs$re$shape_rich,
      "# RE shape (revised, simple)",
      specs$re$shape_simple,
      "",
      "# Interaction counts formula",
      specs$interaction$counts,
      "# Interaction shape (rich reference)",
      specs$interaction$shape_rich,
      "# Interaction shape (revised, simple)",
      specs$interaction$shape_simple
    ),
    file.path(out_dir, "formula_spec.txt")
  )

  pilot_log(glue::glue("Done. Artefacts written to {out_dir}"))

  invisible(list(
    baseline_summary = baseline_summary,
    shape_effect = shape_effect,
    draws = draws_df
  ))
}
