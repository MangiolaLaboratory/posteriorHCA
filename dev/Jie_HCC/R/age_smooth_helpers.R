# Age-smooth (continuous age, tissue×sex thin-plate) helpers for the k-grid
# targets pipeline and local pilots.

age_smooth_ks_default <- function() {
  c(3L, 4L, 5L)
}

# One row per k; crossed with biomarkers in the tar script.
fit_grid_spec_age_smooth <- function(ks = age_smooth_ks_default()) {
  tibble::tibble(
    k = as.integer(ks),
    model = paste0("age_smooth_k", k),
    file_pattern = "fit_age_smooth_hmc_{gene}_k{k}.rds"
  )
}

# Build brms bf() with smooth basis size k (formula, not string — s() needs k).
build_age_smooth_bf <- function(k) {
  k <- as.integer(k)
  counts_formula <- stats::as.formula(
    substitute(
      counts ~
        offset(offset) +
        sex +
        s(
          age_days_scaled,
          by = tissue_sex,
          bs = "tp",
          k = K,
          m = 1,
          id = "age_by_tissue_sex"
        ) +
        ethnicity_groups +
        assay_groups +
        (1 | dataset_id) +
        (1 + sex || tissue_groups) +
        (0 + ethnicity_groups || tissue_groups),
      list(K = k)
    )
  )
  shape_formula <- shape ~ 1 + assay_groups + (1 | tissue_groups)
  brms::bf(counts_formula, shape_formula)
}

age_smooth_formula_text <- function(k) {
  k <- as.integer(k)
  c(
    paste0(
      "counts ~ offset(offset) + sex + ",
      "s(age_days_scaled, by = tissue_sex, bs = \"tp\", k = ", k,
      ", m = 1, id = \"age_by_tissue_sex\") + ",
      "ethnicity_groups + assay_groups + (1 | dataset_id) + ",
      "(1 + sex || tissue_groups) + (0 + ethnicity_groups || tissue_groups)"
    ),
    "shape ~ 1 + assay_groups + (1 | tissue_groups)"
  )
}

# Attach continuous age + assay/dataset + tissue×sex from V1_nk SE colData.
# Keep healthy samples only (disease_groups_altered == "Normal").
prepare_age_smooth_data <- function(
    fit_data,
    col_data,
    healthy_only = TRUE,
    healthy_level = "Normal"
) {
  dat <- droplevels(fit_data)
  colnames(dat) <- stringr::str_replace_all(colnames(dat), "_+", "_")
  stopifnot(nrow(col_data) == nrow(dat))

  if ("age_days_scaled" %in% names(col_data)) {
    dat$age_days_scaled <- as.numeric(col_data$age_days_scaled)
  } else if ("age_days" %in% names(col_data)) {
    dat$age_days <- as.numeric(as.character(col_data$age_days))
    dat$age_days_scaled <- as.numeric(scale(dat$age_days))
  } else {
    stop("col_data needs age_days_scaled or age_days")
  }

  dat$assay_groups <- factor(col_data$assay_groups)
  dat$dataset_id <- factor(col_data$dataset_id)

  if (isTRUE(healthy_only)) {
    if (!"disease_groups_altered" %in% names(dat)) {
      stop("fit data lacks disease_groups_altered; cannot filter healthy samples")
    }
    n_before <- nrow(dat)
    dat <- dat[as.character(dat$disease_groups_altered) == healthy_level, , drop = FALSE]
    dat <- droplevels(dat)
    message(glue::glue(
      "Healthy filter ({healthy_level}): {nrow(dat)} / {n_before} rows retained"
    ))
    if (nrow(dat) < 10L) {
      stop("Too few healthy rows after filter: ", nrow(dat))
    }
  }

  dat$tissue_sex <- interaction(
    dat$tissue_groups, dat$sex,
    drop = TRUE, sep = "_by_"
  )
  dat
}

# HMC fitter matching HCC_age_smooth_one_gene_pilot.R, with control = NULL.
fit_age_smooth_gene <- function(
    data,
    gene,
    k,
    chains = 4L,
    cores = 8L,
    threads_per_chain = 2L,
    warmup = 400L,
    iter = 900L,
    seed = 20260720L,
    control = NULL
) {
  data <- droplevels(data)
  form <- build_age_smooth_bf(k)

  i <- mean(log1p(data$counts / exp(data$offset)), na.rm = TRUE)
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
  n_cores <- min(as.integer(cores), n_chains)
  threads_arg <- if (as.integer(threads_per_chain) > 1L) {
    brms::threading(threads = as.integer(threads_per_chain))
  } else {
    NULL
  }

  message(glue::glue(
    "Fitting age-smooth k={k} for {gene} ",
    "(chains={n_chains}, cores={n_cores}, threads={threads_per_chain}, ",
    "warmup={warmup}, iter={iter}, control=NULL) ..."
  ))

  brms::brm(
    formula = form,
    data = data,
    family = brms::zero_inflated_negbinomial(),
    prior = prior,
    backend = "cmdstanr",
    opencl = NULL,
    chains = n_chains,
    cores = n_cores,
    threads = threads_arg,
    warmup = as.integer(warmup),
    iter = as.integer(iter),
    seed = as.integer(seed),
    refresh = 50,
    silent = 0,
    control = control
  )
}

#------------------------------------------------------------------------------
# Evaluation helpers (report: Jie_HCC_age_smooth_k_grid_eval.qmd)
#------------------------------------------------------------------------------

age_smooth_inventory <- function(
    out_dir,
    genes = c("SPON2", "ZFP36L2", "ZFP36", "VIM"),
    ks = age_smooth_ks_default()
) {
  grid <- tidyr::expand_grid(gene = genes, k = as.integer(ks)) |>
    dplyr::mutate(
      model = paste0("age_smooth_k", k),
      rds_path = file.path(
        out_dir,
        paste0("fit_age_smooth_hmc_", gene, "_k", k, ".rds")
      ),
      available = file.exists(rds_path)
    )
  grid
}

# Comprehensive healthy liver grid for age-smooth fits.
# Margins: sex × ethnicity × assay; age on a short quantile grid (then averaged
# in linpred); tissue fixed to liver; tissue_sex = liver_by_{sex}.
build_comprehensive_grid_age_smooth <- function(
    fit,
    tissue_group = "liver",
    n_age = 5L,
    marginalize_over = c("sex", "ethnicity_groups", "assay_groups"),
    new_level_cols = "dataset_id",
    new_level_tag = "__comprehensive_healthy_baseline__"
) {
  fit_vars <- colnames(fit$data)
  stopifnot("age_days_scaled" %in% fit_vars)

  age_vals <- as.numeric(stats::quantile(
    fit$data$age_days_scaled,
    probs = seq(0.1, 0.9, length.out = as.integer(n_age)),
    na.rm = TRUE,
    names = FALSE
  ))

  levels_list <- lapply(intersect(marginalize_over, fit_vars), function(col) {
    fit_levels(fit, col)
  })
  names(levels_list) <- intersect(marginalize_over, fit_vars)
  levels_list <- levels_list[vapply(levels_list, length, integer(1L)) > 0L]
  levels_list$age_days_scaled <- age_vals

  grid <- do.call(tidyr::expand_grid, levels_list)
  out <- tibble::as_tibble(grid)
  out$tissue_groups <- tissue_group
  out$disease_groups_altered <- "Normal"
  out$offset <- 0

  if ("tissue_sex" %in% fit_vars) {
    out$tissue_sex <- paste(tissue_group, out$sex, sep = "_by_")
  }

  for (col in new_level_cols) {
    if (col %in% fit_vars) out[[col]] <- new_level_tag
  }

  # Align factor levels to the fit
  for (col in intersect(names(out), fit_vars)) {
    if (is.factor(fit$data[[col]])) {
      out[[col]] <- factor(as.character(out[[col]]), levels = levels(fit$data[[col]]))
    } else if (is.character(fit$data[[col]])) {
      out[[col]] <- as.character(out[[col]])
    }
  }
  if ("age_days_scaled" %in% names(out)) {
    out$age_days_scaled <- as.numeric(out$age_days_scaled)
  }
  out
}

linpred_baseline_draws_age_smooth <- function(
    fit,
    ndraws = NULL,
    seed = 20260722L,
    tissue_group = "liver",
    n_age = 5L
) {
  nd_grid <- build_comprehensive_grid_age_smooth(
    fit,
    tissue_group = tissue_group,
    n_age = n_age
  )
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

# Most frequent observed value, with factor-level order used to break ties.
age_smooth_modal_value <- function(x) {
  x_chr <- as.character(x)
  x_chr <- x_chr[!is.na(x_chr)]
  if (!length(x_chr)) {
    stop("Cannot choose a modal value from an all-missing variable")
  }
  counts <- table(x_chr)
  tied <- names(counts)[counts == max(counts)]
  if (is.factor(x)) {
    tied <- levels(x)[levels(x) %in% tied]
  } else {
    tied <- sort(tied)
  }
  tied[[1]]
}

# Blood / liver age curves via posterior_linpred (re_formula = NA → population
# level). Fixed-effect nuisance factors are held at one explicit, reproducible
# reference profile instead of being inherited from an arbitrary training row.
age_smooth_curve_draws <- function(
    fit,
    tissues = c("blood", "liver"),
    n_age = 40L,
    ndraws = 200L,
    re_formula = NA,
    reference_cols = c("ethnicity_groups", "assay_groups"),
    reference_values = NULL,
    seed = 20260722L
) {
  dat <- fit$data
  sex_lv <- sort(unique(as.character(dat$sex)))

  reference_cols <- intersect(reference_cols, names(dat))
  modal_values <- lapply(dat[reference_cols], age_smooth_modal_value)
  if (is.null(reference_values)) reference_values <- list()
  reference_values <- utils::modifyList(modal_values, reference_values)
  unknown_reference_cols <- setdiff(names(reference_values), names(dat))
  if (length(unknown_reference_cols)) {
    stop(
      "Unknown reference_values column(s): ",
      paste(unknown_reference_cols, collapse = ", ")
    )
  }

  template_row <- function(tissue, sex, age) {
    hit <- dat$tissue_groups == tissue & as.character(dat$sex) == sex
    if (!any(hit)) hit <- dat$tissue_groups == tissue
    if (!any(hit)) stop("No training rows for tissue=", tissue)
    row <- dat[which(hit)[[1]], , drop = FALSE]
    row$age_days_scaled <- age
    row$sex <- sex
    row$tissue_groups <- tissue
    if ("tissue_sex" %in% names(row)) {
      row$tissue_sex <- paste(tissue, sex, sep = "_by_")
    }
    row$offset <- 0
    for (col in names(reference_values)) {
      row[[col]] <- reference_values[[col]]
    }
    # Align factors
    for (col in names(row)) {
      if (is.factor(dat[[col]])) {
        row[[col]] <- factor(as.character(row[[col]]), levels = levels(dat[[col]]))
      }
    }
    row
  }

  rows <- vector("list", length(tissues) * length(sex_lv) * as.integer(n_age))
  idx <- 0L
  meta <- vector("list", length(rows))
  for (tissue in tissues) {
    for (sex in sex_lv) {
      support_rows <- dat$tissue_groups == tissue & as.character(dat$sex) == sex
      if (!any(support_rows)) support_rows <- dat$tissue_groups == tissue
      if (!any(support_rows)) stop("No training rows for tissue=", tissue)
      age_rng <- as.numeric(stats::quantile(
        dat$age_days_scaled[support_rows],
        probs = c(0.05, 0.95),
        na.rm = TRUE,
        names = FALSE
      ))
      age_grid <- seq(
        age_rng[[1]], age_rng[[2]], length.out = as.integer(n_age)
      )
      for (age in age_grid) {
        idx <- idx + 1L
        rows[[idx]] <- template_row(tissue, sex, age)
        meta[[idx]] <- tibble::tibble(
          tissue_groups = tissue,
          sex = sex,
          age_days_scaled = age,
          support_lower = age_rng[[1]],
          support_upper = age_rng[[2]],
          reference_profile = paste(
            paste(names(reference_values), unlist(reference_values), sep = "="),
            collapse = "; "
          )
        )
      }
    }
  }
  nd <- dplyr::bind_rows(rows)
  meta_df <- dplyr::bind_rows(meta)

  if (!is.null(seed)) set.seed(as.integer(seed))
  lp <- brms::posterior_linpred(
    fit,
    newdata = nd,
    transform = FALSE,
    re_formula = re_formula,
    allow_new_levels = TRUE,
    sample_new_levels = "gaussian",
    ndraws = ndraws
  )
  # summarise posterior per newdata row
  meta_df |>
    dplyr::mutate(
      estimate = apply(lp, 2L, mean),
      lower = apply(lp, 2L, stats::quantile, probs = 0.025),
      upper = apply(lp, 2L, stats::quantile, probs = 0.975)
    )
}

# Smooth-term-only curves from brms::conditional_smooths(). This deliberately
# matches age_smooth_curve_draws() on tissues, sexes, age range, resolution,
# posterior draw count and interval probability. conditional_smooths() reports
# the posterior median and pointwise quantiles of the selected smooth term; it
# does not include the intercept, other population-level terms, group-level
# terms or the offset.
age_smooth_conditional_smooths <- function(
    fit,
    tissues = c("blood", "liver"),
    n_age = 40L,
    ndraws = 200L,
    prob = 0.95,
    seed = 20260722L
) {
  dat <- fit$data
  required <- c("age_days_scaled", "tissue_groups", "sex", "tissue_sex")
  missing_required <- setdiff(required, names(dat))
  if (length(missing_required)) {
    stop(
      "Fit data lacks required column(s): ",
      paste(missing_required, collapse = ", ")
    )
  }

  sex_lv <- sort(unique(as.character(dat$sex)))
  tissue_sex_lookup <- tidyr::expand_grid(
    tissue_groups = tissues,
    sex = sex_lv
  ) |>
    dplyr::mutate(tissue_sex = paste(tissue_groups, sex, sep = "_by_")) |>
    dplyr::filter(tissue_sex %in% unique(as.character(dat$tissue_sex)))

  if (!nrow(tissue_sex_lookup)) {
    stop("None of the requested tissue-by-sex combinations occur in the fit")
  }

  out <- vector("list", nrow(tissue_sex_lookup))
  for (i in seq_len(nrow(tissue_sex_lookup))) {
    target <- tissue_sex_lookup[i, ]
    support_rows <- as.character(dat$tissue_sex) == target$tissue_sex
    age_rng <- as.numeric(stats::quantile(
      dat$age_days_scaled[support_rows],
      probs = c(0.05, 0.95),
      na.rm = TRUE,
      names = FALSE
    ))
    age_grid <- seq(
      age_rng[[1]], age_rng[[2]], length.out = as.integer(n_age)
    )

    # Resetting to the same seed makes each condition use the same draw subset,
    # which avoids adding draw-selection noise to tissue/sex comparisons.
    if (!is.null(seed)) set.seed(as.integer(seed))
    smooth_list <- brms::conditional_smooths(
      fit,
      int_conditions = list(
        age_days_scaled = age_grid,
        tissue_sex = target$tissue_sex
      ),
      prob = prob,
      resolution = as.integer(n_age),
      ndraws = ndraws
    )
    matching <- which(vapply(
      smooth_list,
      function(x) {
        all(c(
          "age_days_scaled", "tissue_sex", "estimate__", "lower__", "upper__"
        ) %in% names(x))
      },
      logical(1L)
    ))
    if (length(matching) != 1L) {
      stop(
        "Expected one age-by-tissue_sex conditional smooth; found ",
        length(matching)
      )
    }

    out[[i]] <- tibble::as_tibble(smooth_list[[matching]]) |>
      dplyr::transmute(
        tissue_groups = target$tissue_groups,
        sex = target$sex,
        tissue_sex = as.character(tissue_sex),
        age_days_scaled,
        support_lower = age_rng[[1]],
        support_upper = age_rng[[2]],
        estimate = estimate__,
        lower = lower__,
        upper = upper__
      )
  }
  dplyr::bind_rows(out)
}

# Quick compatibility checks for predict / linpred on a smooth fit.
age_smooth_predict_compat <- function(fit, ndraws = 50L, seed = 20260722L) {
  nd <- fit$data[seq_len(min(5L, nrow(fit$data))), , drop = FALSE]
  if (!is.null(seed)) set.seed(as.integer(seed))

  run_one <- function(label, expr) {
    t0 <- Sys.time()
    out <- tryCatch(expr, error = function(e) e)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(out, "error")) {
      tibble::tibble(
        check = label,
        ok = FALSE,
        detail = conditionMessage(out),
        elapsed_sec = elapsed
      )
    } else {
      tibble::tibble(
        check = label,
        ok = TRUE,
        detail = paste0("dim=", paste(dim(out), collapse = "x")),
        elapsed_sec = elapsed
      )
    }
  }

  dplyr::bind_rows(
    run_one("posterior_linpred (re_formula=NULL)", {
      brms::posterior_linpred(
        fit, newdata = nd, transform = FALSE, re_formula = NULL,
        allow_new_levels = TRUE, sample_new_levels = "gaussian",
        ndraws = ndraws
      )
    }),
    run_one("posterior_linpred (re_formula=NA)", {
      brms::posterior_linpred(
        fit, newdata = nd, transform = FALSE, re_formula = NA,
        ndraws = ndraws
      )
    }),
    run_one("posterior_epred", {
      brms::posterior_epred(
        fit, newdata = nd, re_formula = NA, ndraws = ndraws
      )
    }),
    run_one("posterior_predict", {
      brms::posterior_predict(
        fit, newdata = nd, re_formula = NA, ndraws = ndraws
      )
    }),
    run_one("comprehensive liver linpred", {
      d <- linpred_baseline_draws_age_smooth(fit, ndraws = ndraws, seed = seed)
      matrix(d$baseline_log_mu, ncol = 1L)
    })
  )
}

diagnose_age_smooth_fits <- function(
    out_dir,
    genes = c("SPON2", "ZFP36L2", "ZFP36", "VIM"),
    ks = age_smooth_ks_default(),
    compute_baseline = TRUE,
    n_baseline_draws = NULL,
    cache_csv = NULL
) {
  inv <- age_smooth_inventory(out_dir, genes = genes, ks = ks) |>
    dplyr::filter(.data$available)
  rows <- vector("list", nrow(inv))
  for (i in seq_len(nrow(inv))) {
    row <- inv[i, ]
    message(sprintf(
      "[%d/%d] diagnosing %s k=%d ...",
      i, nrow(inv), row$gene, row$k
    ))
    fit <- readRDS(row$rds_path)
    settings <- fit_stan_settings(fit)
    sampler <- fit_sampler_summary(fit)
    rhat_ess <- fit_rhat_ess_summary(fit)

    baseline <- tibble::tibble(
      baseline_mean = NA_real_,
      baseline_sd = NA_real_,
      baseline_n_draws = NA_integer_
    )
    if (isTRUE(compute_baseline)) {
      d <- linpred_baseline_draws_age_smooth(fit, ndraws = n_baseline_draws)
      baseline <- tibble::tibble(
        baseline_mean = mean(d$baseline_log_mu),
        baseline_sd = stats::sd(d$baseline_log_mu),
        baseline_n_draws = nrow(d)
      )
    }

    gate_pass <- isTRUE(rhat_ess$max_rhat <= 1.05) &&
      isTRUE(sampler$pct_divergent < 1)

    rows[[i]] <- dplyr::bind_cols(
      tibble::tibble(
        gene = row$gene,
        k = as.integer(row$k),
        model = row$model,
        rds_path = row$rds_path,
        gate_pass = gate_pass
      ),
      settings,
      sampler,
      rhat_ess,
      baseline
    )
  }
  out <- dplyr::bind_rows(rows)
  if (!is.null(cache_csv)) {
    dir.create(dirname(cache_csv), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(out, cache_csv)
  }
  out
}
