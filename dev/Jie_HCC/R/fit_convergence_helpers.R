# Convergence / MCMC diagnostics for the RE vs interaction pilot fits.
#
# Inventory covers:
#   - fit_re_<GENE>.rds              : random-effect, simple shape (pilot)
#   - fit_interaction_<GENE>.rds     : interaction, simple shape (pilot)
#   - fit_interaction_rich_<GENE>.rds: interaction, rich shape (pilot; ZFP36 may be strong refit)
#   - production V1_nk store         : random-effect, rich shape (reference)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Shared comprehensive-liver baseline draws (deterministic across report sections).
# sample_new_levels = "gaussian" is RNG; always set seed before linpred.
BASELINE_LINPRED_SEED <- 20260722L

comprehensive_baseline_log_mu <- function(fit, seed = BASELINE_LINPRED_SEED) {
  nd <- build_comprehensive_grid(fit)
  if (!is.null(seed)) set.seed(as.integer(seed))
  lp <- brms::posterior_linpred(
    fit,
    newdata = nd,
    transform = FALSE,
    re_formula = NULL,
    allow_new_levels = TRUE,
    sample_new_levels = "gaussian",
    ndraws = NULL
  )
  rowMeans(lp)
}

# Max R̂ over all parameters — same quantity as posterior::rhat(fit).
# Do NOT call posterior::rhat() on a draws_array: that returns a scalar, not
# a per-parameter vector (posterior ≥1.6 quirk / method dispatch).
max_rhat_from_draws <- function(draws) {
  sm <- posterior::summarise_draws(
    draws,
    posterior::default_convergence_measures()
  )
  rh <- sm$rhat
  rh <- rh[is.finite(rh)]
  if (!length(rh)) {
    return(list(max_rhat = NA_real_, mean_rhat = NA_real_, n_rhat_gt_1_05 = NA_integer_))
  }
  list(
    max_rhat = max(rh),
    mean_rhat = mean(rh),
    n_rhat_gt_1_05 = as.integer(sum(rh > 1.05))
  )
}

pilot_fit_inventory <- function(
    out_dir,
    case_dir = "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC",
    model_store_root = "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE",
    target_cell_type = "nk",
    genes = c("SPON2", "ZFP36L2", "ZFP36", "VIM")
) {
  biomarkers <- readr::read_csv(
    file.path(case_dir, "data", "times_biomarkers.csv"),
    show_col_types = FALSE
  ) |>
    dplyr::filter(.data$gene %in% genes)

  store <- file.path(model_store_root, paste0("V1_", target_cell_type), "_targets")
  model_meta <- targets::tar_meta(starts_with("estimates_chunk_"), store = store) |>
    dplyr::filter(.data$type == "branch") |>
    dplyr::mutate(.feature = stringr::str_extract(.data$warnings, "(?<=Gene:___)ENSG\\d+(?=___)"))

  local_rows <- tidyr::expand_grid(
    gene = biomarkers$gene,
    model = c("re", "interaction"),
    shape_variant = c("simple", "rich")
  ) |>
    dplyr::left_join(biomarkers |> dplyr::select(gene, .feature), by = "gene") |>
    dplyr::mutate(
      role = dplyr::if_else(shape_variant == "simple", "revised", "reference"),
      file_name = dplyr::case_when(
        model == "re" & shape_variant == "simple" ~ paste0("fit_re_", gene, ".rds"),
        model == "re" & shape_variant == "rich" ~ paste0("fit_re_rich_", gene, ".rds"),
        model == "interaction" & shape_variant == "simple" ~ paste0("fit_interaction_", gene, ".rds"),
        model == "interaction" & shape_variant == "rich" ~ paste0("fit_interaction_rich_", gene, ".rds")
      ),
      path = file.path(out_dir, file_name),
      source = dplyr::case_when(
        model == "re" & shape_variant == "rich" ~ "production_store",
        file.exists(path) ~ "local_rds",
        TRUE ~ "missing"
      ),
      available = source != "missing" | (model == "re" & shape_variant == "rich" &
        .feature %in% model_meta$.feature)
    )

  # Mark production RE-rich as available when the gene is in the store,
  # even if no local fit_re_rich_*.rds exists.
  local_rows <- local_rows |>
    dplyr::mutate(
      available = dplyr::if_else(
        model == "re" & shape_variant == "rich" & .feature %in% model_meta$.feature,
        TRUE,
        available
      ),
      source = dplyr::if_else(
        model == "re" & shape_variant == "rich" & !file.exists(path) &
          .feature %in% model_meta$.feature,
        "production_store",
        source
      )
    )

  list(
    inventory = local_rows,
    biomarkers = biomarkers,
    model_meta = model_meta,
    store = store
  )
}

load_pilot_fit <- function(row, inventory_meta) {
  if (isTRUE(row$source == "local_rds") && file.exists(row$path)) {
    return(readRDS(row$path))
  }
  if (isTRUE(row$source == "production_store")) {
    target_name <- inventory_meta$model_meta |>
      dplyr::filter(.data$.feature == row$.feature) |>
      dplyr::slice(1) |>
      dplyr::pull(.data$name)
    if (!length(target_name) || is.na(target_name[[1]])) {
      stop("No production target for feature ", row$.feature)
    }
    return(targets::tar_read_raw(target_name[[1]], store = inventory_meta$store)$brms_fit[[1]])
  }
  stop("Fit not available: ", row$gene, " / ", row$model, " / ", row$shape_variant)
}

fit_stan_settings <- function(fit) {
  a <- fit$fit@stan_args[[1]]
  alg <- as.character(fit$algorithm %||% a$algorithm %||% a$method %||% NA_character_)
  is_pf <- identical(alg, "pathfinder") || identical(a$method, "pathfinder")

  chains <- as.integer(fit$fit@sim$chains %||% NA_integer_)
  iter <- as.integer(fit$fit@sim$iter %||% NA_integer_)
  warmup <- as.integer(fit$fit@sim$warmup %||% NA_integer_)
  thin <- as.integer(fit$fit@sim$thin %||% 1L)
  if (!length(chains)) chains <- NA_integer_
  if (!length(iter)) iter <- NA_integer_
  if (!length(warmup)) warmup <- NA_integer_
  if (!length(thin)) thin <- 1L

  n_draws <- tryCatch(
    as.integer(posterior::ndraws(posterior::as_draws_df(fit))),
    error = function(e) NA_integer_
  )

  post_warmup <- if (is_pf) {
    n_draws
  } else if (is.finite(iter) && is.finite(warmup) && is.finite(thin)) {
    as.integer((iter - warmup) / thin)
  } else {
    NA_integer_
  }

  tibble::tibble(
    chains = chains,
    iter = iter,
    warmup = warmup,
    thin = thin,
    post_warmup_per_chain = post_warmup,
    n_draws = n_draws,
    adapt_delta = as.numeric(a$control$adapt_delta %||% a$adapt_delta %||% NA_real_),
    max_treedepth = as.integer(a$control$max_treedepth %||% a$max_treedepth %||% NA_integer_),
    seed = as.integer(a$seed %||% NA_integer_),
    algorithm = alg,
    backend = as.character(fit$backend %||% NA_character_),
    is_pathfinder = is_pf
  )
}

fit_sampler_summary <- function(fit) {
  np <- tryCatch(brms::nuts_params(fit), error = function(e) NULL)
  if (is.null(np)) {
    return(tibble::tibble(
      n_divergent = NA_integer_,
      pct_divergent = NA_real_,
      max_treedepth_hit = NA_integer_,
      pct_max_treedepth = NA_real_,
      mean_treedepth = NA_real_
    ))
  }
  div <- np |> dplyr::filter(.data$Parameter == "divergent__")
  td <- np |> dplyr::filter(.data$Parameter == "treedepth__")
  settings <- fit_stan_settings(fit)
  cap <- settings$max_treedepth[[1]]
  tibble::tibble(
    n_divergent = as.integer(sum(div$Value)),
    pct_divergent = 100 * mean(div$Value),
    max_treedepth_hit = as.integer(max(td$Value)),
    pct_max_treedepth = if (is.na(cap)) NA_real_ else 100 * mean(td$Value >= cap),
    mean_treedepth = mean(td$Value)
  )
}

fit_rhat_ess_summary <- function(fit) {
  rh <- posterior::rhat(fit)
  ess_b <- tryCatch(posterior::ess_bulk(fit), error = function(e) rep(NA_real_, length(rh)))
  ess_t <- tryCatch(posterior::ess_tail(fit), error = function(e) rep(NA_real_, length(rh)))
  ok <- is.finite(rh)
  tibble::tibble(
    n_params = sum(ok),
    max_rhat = max(rh[ok], na.rm = TRUE),
    mean_rhat = mean(rh[ok], na.rm = TRUE),
    n_rhat_gt_1_01 = sum(rh[ok] > 1.01),
    n_rhat_gt_1_05 = sum(rh[ok] > 1.05),
    n_rhat_gt_1_1 = sum(rh[ok] > 1.1),
    n_rhat_gt_1_5 = sum(rh[ok] > 1.5),
    min_ess_bulk = suppressWarnings(min(ess_b[ok], na.rm = TRUE)),
    min_ess_tail = suppressWarnings(min(ess_t[ok], na.rm = TRUE)),
    median_ess_bulk = suppressWarnings(stats::median(ess_b[ok], na.rm = TRUE))
  )
}

worst_parameters <- function(fit, n = 8L) {
  sm <- posterior::summarise_draws(posterior::as_draws_df(fit))
  sm |>
    dplyr::filter(is.finite(.data$rhat)) |>
    dplyr::arrange(dplyr::desc(.data$rhat), .data$ess_bulk) |>
    dplyr::slice_head(n = n) |>
    dplyr::select(variable, mean, sd, rhat, ess_bulk, ess_tail)
}

diagnose_one_fit <- function(
    row,
    inventory_meta,
    compute_baseline = TRUE,
    n_baseline_draws = NULL
) {
  fit <- load_pilot_fit(row, inventory_meta)
  settings <- fit_stan_settings(fit)
  sampler <- fit_sampler_summary(fit)
  rhat_ess <- fit_rhat_ess_summary(fit)

  baseline <- tibble::tibble(
    baseline_mean = NA_real_,
    baseline_sd = NA_real_,
    baseline_chain_mean_range = NA_real_,
    baseline_max_chain_frac_lt5 = NA_real_,
    baseline_n_draws = NA_integer_
  )
  if (isTRUE(compute_baseline) && exists("build_comprehensive_grid", mode = "function")) {
    n_chains <- settings$chains[[1]]
    n_iter <- settings$post_warmup_per_chain[[1]]
    # NULL / NA / "all" → every post-warmup draw (matches draw-stability endpoint).
    use_all <- is.null(n_baseline_draws) ||
      (is.character(n_baseline_draws) && identical(n_baseline_draws, "all")) ||
      (is.numeric(n_baseline_draws) && !is.finite(n_baseline_draws))
    if (use_all) {
      y <- comprehensive_baseline_log_mu(fit)
    } else {
      nd <- build_comprehensive_grid(fit)
      set.seed(BASELINE_LINPRED_SEED)
      lp <- brms::posterior_linpred(
        fit,
        newdata = nd,
        transform = FALSE,
        re_formula = NULL,
        allow_new_levels = TRUE,
        sample_new_levels = "gaussian",
        ndraws = as.integer(n_baseline_draws)
      )
      y <- rowMeans(lp)
    }
    n_used <- length(y)
    can_label_chains <- use_all &&
      is.finite(n_chains) && is.finite(n_iter) &&
      !isTRUE(settings$is_pathfinder[[1]]) &&
      n_used == as.integer(n_chains * n_iter)
    if (can_label_chains) {
      chain <- rep(seq_len(n_chains), each = n_iter)
      by_chain <- tapply(y, chain, mean)
      frac_lt5 <- tapply(y, chain, function(z) mean(z < 5))
      baseline <- tibble::tibble(
        baseline_mean = mean(y),
        baseline_sd = stats::sd(y),
        baseline_chain_mean_range = as.numeric(diff(range(by_chain))),
        baseline_max_chain_frac_lt5 = max(as.numeric(frac_lt5)),
        baseline_n_draws = as.integer(n_used)
      )
    } else {
      baseline <- tibble::tibble(
        baseline_mean = mean(y),
        baseline_sd = stats::sd(y),
        baseline_chain_mean_range = NA_real_,
        baseline_max_chain_frac_lt5 = mean(y < 5),
        baseline_n_draws = as.integer(n_used)
      )
    }
  }

  # Convergence gate used by this report.
  # Pathfinder is not multi-chain HMC: skip R̂ / divergence / treedepth gates;
  # require finite baseline summaries and a usable draw count.
  if (isTRUE(settings$is_pathfinder[[1]])) {
    gate_pass <- is.finite(baseline$baseline_mean[[1]]) &&
      is.finite(baseline$baseline_sd[[1]]) &&
      isTRUE(settings$n_draws[[1]] >= 100)
  } else {
    gate_pass <- isTRUE(rhat_ess$max_rhat <= 1.05) &&
      isTRUE(sampler$pct_divergent < 1) &&
      (is.na(sampler$pct_max_treedepth) || sampler$pct_max_treedepth < 10) &&
      (is.na(baseline$baseline_chain_mean_range[[1]]) ||
         baseline$baseline_chain_mean_range[[1]] < 1.0)
  }

  dplyr::bind_cols(
    tibble::tibble(
      gene = row$gene,
      .feature = row$.feature,
      model = row$model,
      shape_variant = row$shape_variant,
      role = row$role,
      source = row$source,
      path = row$path,
      gate_pass = gate_pass
    ),
    settings,
    sampler,
    rhat_ess,
    baseline
  )
}

diagnose_all_fits <- function(
    out_dir,
    case_dir = "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC",
    model_store_root = "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE",
    target_cell_type = "nk",
    genes = c("SPON2", "ZFP36L2", "ZFP36", "VIM"),
    compute_baseline = TRUE,
    n_baseline_draws = NULL,
    cache_csv = NULL
) {
  meta <- pilot_fit_inventory(
    out_dir = out_dir,
    case_dir = case_dir,
    model_store_root = model_store_root,
    target_cell_type = target_cell_type,
    genes = genes
  )
  inv <- meta$inventory |> dplyr::filter(.data$available)
  rows <- vector("list", nrow(inv))
  for (i in seq_len(nrow(inv))) {
    message(sprintf(
      "[%d/%d] diagnosing %s / %s / %s ...",
      i, nrow(inv), inv$gene[i], inv$model[i], inv$shape_variant[i]
    ))
    rows[[i]] <- diagnose_one_fit(
      inv[i, ], meta,
      compute_baseline = compute_baseline,
      n_baseline_draws = n_baseline_draws
    )
  }
  out <- dplyr::bind_rows(rows)
  if (!is.null(cache_csv)) {
    readr::write_csv(out, cache_csv)
  }
  list(diagnostics = out, inventory = meta$inventory, meta = meta)
}

# Draws for MCMC plots: keep a small, interpretable parameter set.
select_plot_variables <- function(fit, n_worst = 4L) {
  v <- posterior::variables(fit)
  core <- intersect(
    c("Intercept", "b_Intercept", "Intercept_shape", "b_shape_Intercept", "zi"),
    v
  )
  worst <- worst_parameters(fit, n = n_worst)$variable
  # Prefer liver-related fixed effects if present among variables
  liver <- grep("tissue_groupsliver", v, value = TRUE)
  liver <- head(liver, 2L)
  unique(c(core, liver, worst))
}

baseline_draws_by_chain <- function(fit, ndraws_cap = 2000L) {
  settings <- fit_stan_settings(fit)
  n_chains <- settings$chains[[1]]
  n_iter <- settings$post_warmup_per_chain[[1]]
  use_n <- min(as.integer(n_chains * n_iter), as.integer(ndraws_cap))
  # If capped, still request all and subsample evenly is hard; prefer all when <= cap.
  use_n <- as.integer(n_chains * n_iter)
  if (use_n > ndraws_cap) {
    # Subsample equally per chain via indices after full linpred is too heavy;
    # for large fits use ndraws_cap without chain labels.
    nd <- build_comprehensive_grid(fit)
    lp <- brms::posterior_linpred(
      fit, newdata = nd, transform = FALSE, re_formula = NULL,
      allow_new_levels = TRUE, sample_new_levels = "gaussian",
      ndraws = ndraws_cap
    )
    return(tibble::tibble(
      draw = seq_len(nrow(lp)),
      chain = NA_integer_,
      baseline_log_mu = rowMeans(lp)
    ))
  }
  nd <- build_comprehensive_grid(fit)
  lp <- brms::posterior_linpred(
    fit, newdata = nd, transform = FALSE, re_formula = NULL,
    allow_new_levels = TRUE, sample_new_levels = "gaussian",
    ndraws = use_n
  )
  tibble::tibble(
    draw = seq_len(use_n),
    chain = rep(seq_len(n_chains), each = n_iter),
    iteration = rep(seq_len(n_iter), times = n_chains),
    baseline_log_mu = rowMeans(lp)
  )
}

convergence_gate_label <- function(diag_row) {
  dplyr::case_when(
    isTRUE(diag_row$gate_pass) ~ "PASS",
    isTRUE(diag_row$max_rhat > 1.5) | isTRUE(diag_row$pct_divergent >= 50) ~ "FAIL",
    TRUE ~ "WARN"
  )
}

#------------------------------------------------------------------------------
# Strong formula comparison: production vs RE / interaction / RE-uncorrelated
# (simple shape; strong MCMC for the three local fits).
#------------------------------------------------------------------------------

strong_formula_inventory <- function(
    strong_out_dir,
    case_dir = "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC",
    model_store_root = "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE",
    target_cell_type = "nk",
    genes = c("SPON2", "ZFP36L2", "ZFP36", "VIM")
) {
  biomarkers <- readr::read_csv(
    file.path(case_dir, "data", "times_biomarkers.csv"),
    show_col_types = FALSE
  ) |>
    dplyr::filter(.data$gene %in% genes)

  store <- file.path(model_store_root, paste0("V1_", target_cell_type), "_targets")
  model_meta <- targets::tar_meta(starts_with("estimates_chunk_"), store = store) |>
    dplyr::filter(.data$type == "branch") |>
    dplyr::mutate(
      .feature = stringr::str_extract(.data$warnings, "(?<=Gene:___)ENSG\\d+(?=___)")
    )

  formula_spec <- tibble::tribble(
    ~model,                  ~shape_variant, ~role,          ~file_pattern,                            ~source_pref,
    "production",            "rich",         "production",   NA_character_,                            "production_store",
    "re",                    "simple",       "strong",       "fit_re_{gene}.rds",                      "local_rds",
    "interaction",           "simple",       "strong",       "fit_interaction_{gene}.rds",             "local_rds",
    "re_uncorrelated",       "simple",       "strong",       "fit_re_uncorrelated_{gene}.rds",         "local_rds",
    "pathfinder_re",         "simple",       "pathfinder",   "fit_pathfinder_re_{gene}.rds",           "local_rds",
    "pathfinder_re_split",   "simple",       "pathfinder",   "fit_pathfinder_re_split_{gene}.rds",     "local_rds"
  )

  inv <- tidyr::crossing(
    biomarkers |> dplyr::select(gene, .feature),
    formula_spec
  ) |>
    dplyr::mutate(
      file_name = dplyr::case_when(
        .data$model == "production" ~ NA_character_,
        .data$model == "re" ~ paste0("fit_re_", .data$gene, ".rds"),
        .data$model == "interaction" ~ paste0("fit_interaction_", .data$gene, ".rds"),
        .data$model == "re_uncorrelated" ~ paste0("fit_re_uncorrelated_", .data$gene, ".rds"),
        .data$model == "pathfinder_re" ~ paste0("fit_pathfinder_re_", .data$gene, ".rds"),
        .data$model == "pathfinder_re_split" ~ paste0("fit_pathfinder_re_split_", .data$gene, ".rds"),
        TRUE ~ NA_character_
      ),
      path = dplyr::if_else(
        is.na(.data$file_name),
        NA_character_,
        file.path(strong_out_dir, .data$file_name)
      ),
      source = dplyr::case_when(
        .data$source_pref == "production_store" & .data$.feature %in% model_meta$.feature ~
          "production_store",
        .data$source_pref == "local_rds" & !is.na(.data$path) & file.exists(.data$path) ~
          "local_rds",
        TRUE ~ "missing"
      ),
      available = .data$source != "missing"
    ) |>
    dplyr::select(-file_pattern, -source_pref)

  list(
    inventory = inv,
    biomarkers = biomarkers,
    model_meta = model_meta,
    store = store
  )
}

load_strong_fit <- function(row, inventory_meta) {
  # Reuse the pilot loader (local_rds / production_store).
  load_pilot_fit(row, inventory_meta)
}

diagnose_strong_formula_fits <- function(
    strong_out_dir,
    case_dir = "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC",
    model_store_root = "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE",
    target_cell_type = "nk",
    genes = c("SPON2", "ZFP36L2", "ZFP36", "VIM"),
    compute_baseline = TRUE,
    n_baseline_draws = NULL,
    cache_csv = NULL
) {
  meta <- strong_formula_inventory(
    strong_out_dir = strong_out_dir,
    case_dir = case_dir,
    model_store_root = model_store_root,
    target_cell_type = target_cell_type,
    genes = genes
  )
  inv <- meta$inventory |> dplyr::filter(.data$available)
  rows <- vector("list", nrow(inv))
  for (i in seq_len(nrow(inv))) {
    message(sprintf(
      "[%d/%d] diagnosing %s / %s ...",
      i, nrow(inv), inv$gene[i], inv$model[i]
    ))
    rows[[i]] <- diagnose_one_fit(
      inv[i, ], meta,
      compute_baseline = compute_baseline,
      n_baseline_draws = n_baseline_draws
    )
  }
  out <- dplyr::bind_rows(rows)
  if (!is.null(cache_csv)) {
    dir.create(dirname(cache_csv), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(out, cache_csv)
  }
  list(diagnostics = out, inventory = meta$inventory, meta = meta)
}

collect_strong_baseline_draws <- function(
    strong_out_dir,
    case_dir = "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC",
    model_store_root = "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE",
    target_cell_type = "nk",
    genes = c("SPON2", "ZFP36L2", "ZFP36", "VIM"),
    n_draws = NULL,
    cache_rds = NULL
) {
  if (!is.null(cache_rds) && file.exists(cache_rds)) {
    return(readRDS(cache_rds))
  }
  meta <- strong_formula_inventory(
    strong_out_dir = strong_out_dir,
    case_dir = case_dir,
    model_store_root = model_store_root,
    target_cell_type = target_cell_type,
    genes = genes
  )
  inv <- meta$inventory |> dplyr::filter(.data$available)
  out <- vector("list", nrow(inv))
  for (i in seq_len(nrow(inv))) {
    row <- inv[i, ]
    message(sprintf(
      "[%d/%d] baseline draws %s / %s ...",
      i, nrow(inv), row$gene, row$model
    ))
    fit <- load_strong_fit(row, meta)
    d <- linpred_baseline_draws(fit, ndraws = n_draws) |>
      dplyr::mutate(
        gene = row$gene,
        .feature = row$.feature,
        model = row$model,
        shape_variant = row$shape_variant,
        source = row$source
      )
    out[[i]] <- d
  }
  draws <- dplyr::bind_rows(out)
  if (!is.null(cache_rds)) {
    dir.create(dirname(cache_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(draws, cache_rds, compress = "xz")
  }
  draws
}

#------------------------------------------------------------------------------
# HMC draw-stability curves (exclude production / Pathfinder)
# Cumulative first-k post-warmup iterations × all chains → R̂, baseline μ, SD
#------------------------------------------------------------------------------

hmc_iter_grid_default <- function(n_iter_per_chain) {
  grid <- c(25L, 50L, 75L, 100L, 150L, 200L, 300L, 400L, 600L, 800L, 1000L, 1200L, 1500L)
  grid <- grid[grid <= as.integer(n_iter_per_chain)]
  if (!length(grid) || max(grid) < as.integer(n_iter_per_chain)) {
    grid <- sort(unique(c(grid, as.integer(n_iter_per_chain))))
  }
  as.integer(grid)
}

# First index where value is within tol of the final value, and stays there.
first_stable_index <- function(x, rel_tol = 0.05, abs_tol = 0.05) {
  x <- as.numeric(x)
  if (!length(x) || all(!is.finite(x))) return(NA_integer_)
  xref <- x[length(x)]
  if (!is.finite(xref)) return(NA_integer_)
  tol <- max(abs_tol, rel_tol * abs(xref))
  ok <- is.finite(x) & abs(x - xref) <= tol
  # require all subsequent points also OK
  stay <- rev(dplyr::cumall(rev(ok)))
  idx <- which(stay)
  if (!length(idx)) return(NA_integer_)
  as.integer(idx[[1]])
}

first_rhat_stable_index <- function(max_rhat, gate = 1.05, near_final_tol = 0.02) {
  x <- as.numeric(max_rhat)
  if (!length(x) || all(!is.finite(x))) return(NA_integer_)
  xref <- x[length(x)]
  near_final <- is.finite(x) & abs(x - xref) <= near_final_tol
  under_gate <- is.finite(x) & x <= gate
  ok <- near_final & (under_gate | !is.finite(gate))
  stay <- rev(dplyr::cumall(rev(ok)))
  idx <- which(stay)
  if (!length(idx)) return(NA_integer_)
  as.integer(idx[[1]])
}

hmc_draw_stability_one_fit <- function(
    fit,
    gene,
    model,
    iter_grid = NULL,
    rhat_variables = NULL
) {
  settings <- fit_stan_settings(fit)
  if (isTRUE(settings$is_pathfinder[[1]])) {
    stop("hmc_draw_stability_one_fit is for HMC fits only")
  }
  n_chains <- as.integer(settings$chains[[1]])
  n_iter <- as.integer(settings$post_warmup_per_chain[[1]])
  if (!is.finite(n_chains) || !is.finite(n_iter) || n_iter < 2L) {
    stop("Fit lacks usable post-warmup iterations for stability curves")
  }
  if (is.null(iter_grid)) iter_grid <- hmc_iter_grid_default(n_iter)

  # --- Baseline trajectory: all draws (shared seeded helper), then cumulative ---
  y <- comprehensive_baseline_log_mu(fit)
  draws_df <- posterior::as_draws_df(fit)
  stopifnot(length(y) == nrow(draws_df))
  draws_df$.baseline_log_mu <- y

  # --- R̂ on ALL parameters via summarise_draws (matches rhat(fit) / main table) ---
  if (is.null(rhat_variables)) {
    draws_arr <- posterior::as_draws_array(fit)
  } else {
    draws_arr <- posterior::as_draws_array(fit, variable = rhat_variables)
  }

  rows <- vector("list", length(iter_grid))
  for (i in seq_along(iter_grid)) {
    k <- as.integer(iter_grid[[i]])
    sub <- draws_df |> dplyr::filter(.data$.iteration <= k)
    n_draws <- nrow(sub)
    bl <- sub$.baseline_log_mu
    by_chain <- tapply(bl, sub$.chain, mean)
    d_sub <- posterior::subset_draws(draws_arr, iteration = seq_len(k))
    rh_sum <- max_rhat_from_draws(d_sub)
    rows[[i]] <- tibble::tibble(
      gene = gene,
      model = model,
      n_iter_per_chain = k,
      n_draws = n_draws,
      max_rhat = rh_sum$max_rhat,
      mean_rhat = rh_sum$mean_rhat,
      n_rhat_gt_1_05 = rh_sum$n_rhat_gt_1_05,
      baseline_mean = mean(bl),
      baseline_sd = stats::sd(bl),
      baseline_chain_mean_range = as.numeric(diff(range(by_chain)))
    )
  }
  dplyr::bind_rows(rows)
}

hmc_draw_stability_thresholds <- function(curve_df) {
  curve_df |>
    dplyr::arrange(.data$gene, .data$model, .data$n_draws) |>
    dplyr::group_by(.data$gene, .data$model) |>
    dplyr::group_modify(function(d, keys) {
      idx_rhat <- first_rhat_stable_index(d$max_rhat)
      idx_mean <- first_stable_index(d$baseline_mean, rel_tol = 0.02, abs_tol = 0.05)
      idx_sd <- first_stable_index(d$baseline_sd, rel_tol = 0.05, abs_tol = 0.05)
      pick_n <- function(idx) {
        if (is.na(idx)) NA_integer_ else as.integer(d$n_draws[[idx]])
      }
      tibble::tibble(
        n_draws_full = max(d$n_draws),
        n_iter_full = max(d$n_iter_per_chain),
        max_rhat_full = dplyr::last(d$max_rhat),
        baseline_mean_full = dplyr::last(d$baseline_mean),
        baseline_sd_full = dplyr::last(d$baseline_sd),
        n_draws_rhat_stable = pick_n(idx_rhat),
        n_draws_mean_stable = pick_n(idx_mean),
        n_draws_sd_stable = pick_n(idx_sd)
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      n_draws_all_stable = pmax(
        .data$n_draws_rhat_stable,
        .data$n_draws_mean_stable,
        .data$n_draws_sd_stable,
        na.rm = FALSE
      )
    )
}

collect_hmc_draw_stability <- function(
    strong_out_dir,
    case_dir = "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC",
    model_store_root = "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE",
    target_cell_type = "nk",
    genes = c("SPON2", "ZFP36L2", "ZFP36", "VIM"),
    models = c("re", "interaction", "re_uncorrelated"),
    iter_grid = NULL,
    cache_rds = NULL,
    force = FALSE
) {
  if (!isTRUE(force) && !is.null(cache_rds) && file.exists(cache_rds)) {
    return(readRDS(cache_rds))
  }
  meta <- strong_formula_inventory(
    strong_out_dir = strong_out_dir,
    case_dir = case_dir,
    model_store_root = model_store_root,
    target_cell_type = target_cell_type,
    genes = genes
  )
  inv <- meta$inventory |>
    dplyr::filter(
      .data$available,
      .data$model %in% models,
      !.data$model %in% c("production"),
      !grepl("^pathfinder", .data$model)
    )
  out <- vector("list", nrow(inv))
  for (i in seq_len(nrow(inv))) {
    row <- inv[i, ]
    message(sprintf(
      "[%d/%d] HMC draw-stability %s / %s ...",
      i, nrow(inv), row$gene, row$model
    ))
    fit <- load_strong_fit(row, meta)
    out[[i]] <- hmc_draw_stability_one_fit(
      fit,
      gene = row$gene,
      model = row$model,
      iter_grid = iter_grid
    )
  }
  curves <- dplyr::bind_rows(out)
  thresholds <- hmc_draw_stability_thresholds(curves)
  result <- list(curves = curves, thresholds = thresholds)
  if (!is.null(cache_rds)) {
    dir.create(dirname(cache_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, cache_rds, compress = "xz")
  }
  result
}
