# Minimal helpers for formula_bench.
# 1) prepare_bench_data()  — once from pseudobulk_sample
# 2) fit_bench_row()       — one param_grid row → timing + nested brms fit

bench_genes <- function() {
  c("SPON2", "ZFP36L2", "ZFP36", "VIM")
}

# Map gene symbol → Ensembl id via times_biomarkers.csv
bench_gene_features <- function(case_dir, genes = bench_genes()) {
  readr::read_csv(
    file.path(case_dir, "data", "times_biomarkers.csv"),
    show_col_types = FALSE
  ) |>
    dplyr::filter(.data$gene %in% genes) |>
    dplyr::select(gene, .feature) |>
    dplyr::arrange(.data$gene)
}

# ---------------------------------------------------------------------------
# 1. Prepare data once (healthy samples, biomarker genes only)
# ---------------------------------------------------------------------------
# Returns a list:
#   meta     — sample tibble (factors, offset, tissue_sex, …)
#   counts   — genes × samples matrix (rownames = gene symbols)
#   features — gene / .feature map
prepare_bench_data <- function(
    se,
    features,
    healthy_level = "Normal"
) {
  stopifnot(inherits(se, "SummarizedExperiment"))
  features <- dplyr::distinct(features, .data$gene, .data$.feature)

  # Collapse ___ → _ so names match formula columns (*_altered)
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  names(cd) <- stringr::str_replace_all(names(cd), "_+", "_")

  needed <- c(
    "offset", "age_decade", "age_days_scaled", "sex", "ethnicity_groups",
    "tissue_groups", "assay_groups", "assay_groups_altered",
    "dataset_id", "dataset_id_altered", "disease_groups_altered"
  )
  missing <- setdiff(needed, names(cd))
  if (length(missing)) {
    stop("pseudobulk colData missing: ", paste(missing, collapse = ", "))
  }

  keep_sample <- as.character(cd$disease_groups_altered) == healthy_level
  cd <- cd[keep_sample, , drop = FALSE]
  message(glue::glue(
    "Healthy filter ({healthy_level}): {sum(keep_sample)} / {length(keep_sample)} samples"
  ))

  feat_ids <- features$.feature
  missing_feat <- setdiff(feat_ids, rownames(se))
  if (length(missing_feat)) {
    stop("Features not in SE: ", paste(missing_feat, collapse = ", "))
  }

  # Assay may be HDF5 / DelayedArray — subset SE then realize as dense matrix
  row_idx <- match(feat_ids, rownames(se))
  col_idx <- which(keep_sample)
  se_sub <- se[row_idx, col_idx, drop = FALSE]
  mat <- as.matrix(SummarizedExperiment::assay(se_sub, "counts"))
  storage.mode(mat) <- "double"
  rownames(mat) <- features$gene[match(rownames(se_sub), features$.feature)]
  colnames(mat) <- colnames(se_sub)

  meta <- tibble::as_tibble(cd[, needed, drop = FALSE])
  meta <- dplyr::mutate(
    meta,
    dplyr::across(
      c(
        "age_decade", "sex", "ethnicity_groups", "tissue_groups",
        "assay_groups", "assay_groups_altered", "dataset_id",
        "dataset_id_altered", "disease_groups_altered"
      ),
      \(x) factor(x)
    ),
    tissue_sex = interaction(
      .data$tissue_groups, .data$sex,
      drop = TRUE, sep = "_by_"
    )
  )

  stopifnot(ncol(mat) == nrow(meta))

  list(
    meta = meta,
    counts = mat,
    features = features,
    n_samples = nrow(meta),
    healthy_level = healthy_level
  )
}

# One gene's modelling frame (counts + meta)
gene_model_data <- function(prepared, gene) {
  if (!gene %in% rownames(prepared$counts)) {
    stop("Gene not in prepared counts: ", gene)
  }
  dplyr::mutate(
    prepared$meta,
    counts = as.numeric(prepared$counts[gene, ])
  )
}

# ---------------------------------------------------------------------------
# 2. Fit one param_grid row
# ---------------------------------------------------------------------------
build_bf_from_row <- function(row) {
  family <- as.character(row$family[[1]])
  if (identical(family, "re_string")) {
    return(brms::bf(
      stats::as.formula(row$counts[[1]]),
      stats::as.formula(row$shape[[1]])
    ))
  }
  if (identical(family, "age_smooth")) {
    k <- as.integer(row$k[[1]])
    # Inline age-smooth bf (same as age_smooth_helpers::build_age_smooth_bf)
    counts_f <- stats::as.formula(substitute(
      counts ~
        offset(offset) + sex +
        s(age_days_scaled, by = tissue_sex, bs = "tp", k = K, m = 1,
          id = "age_by_tissue_sex") +
        ethnicity_groups + assay_groups +
        (1 | dataset_id) +
        (1 + sex || tissue_groups) +
        (0 + ethnicity_groups || tissue_groups),
      list(K = k)
    ))
    return(brms::bf(counts_f, shape ~ 1 + assay_groups + (1 | tissue_groups)))
  }
  stop("Unknown family: ", family)
}

fit_bench_row <- function(row, prepared, mcmc = bench_mcmc(), fits_dir = NULL) {
  stopifnot(nrow(row) == 1L)

  gene <- as.character(row$gene[[1]])
  formula_id <- as.character(row$formula_id[[1]])
  setting_id <- as.character(row$setting_id[[1]])
  chains <- as.integer(row$chains[[1]])
  cpus <- as.integer(row$cpus[[1]])
  threads <- bench_resolve_threads(chains, cpus)
  chain_cores <- min(chains, cpus)

  data <- droplevels(gene_model_data(prepared, gene))
  form <- build_bf_from_row(row)

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

  threads_arg <- if (threads > 1L) brms::threading(threads = threads) else NULL

  seed_i <- as.integer(mcmc$seed) +
    sum(utf8ToInt(paste(gene, formula_id, setting_id, sep = "|"))) %% 100000L

  # Prefer scratch for CmdStan temps when present
  if (dir.exists("/scratchdata1/users/a1237163")) {
    tmpdir <- file.path(
      "/scratchdata1/users/a1237163/tmp",
      Sys.getenv("SLURM_JOB_ID", unset = as.character(Sys.getpid()))
    )
    dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
    Sys.setenv(TMPDIR = tmpdir, TMP = tmpdir, TEMP = tmpdir)
  }

  message(glue::glue(
    "FIT {gene} | {formula_id} | {setting_id} | ",
    "chains={chains} cpus={cpus} threads={threads} n={nrow(data)}"
  ))

  t0 <- Sys.time()
  err <- NA_character_
  fit_obj <- tryCatch(
    brms::brm(
      formula = form,
      data = data,
      family = brms::zero_inflated_negbinomial(),
      prior = prior,
      backend = "cmdstanr",
      opencl = NULL,
      chains = chains,
      cores = chain_cores,
      threads = threads_arg,
      warmup = as.integer(mcmc$warmup),
      iter = as.integer(mcmc$iter),
      seed = seed_i,
      refresh = 50,
      silent = 1,
      control = NULL
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      message("FIT ERROR: ", err)
      NULL
    }
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  rds_path <- NA_character_
  if (!is.null(fit_obj) && !is.null(fits_dir)) {
    dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
    rds_path <- file.path(
      fits_dir,
      paste0("fit_", gene, "__", formula_id, "__", setting_id, ".rds")
    )
    saveRDS(fit_obj, rds_path, compress = "xz")
  }

  max_rhat <- NA_real_
  if (!is.null(fit_obj)) {
    rh <- tryCatch(posterior::rhat(fit_obj), error = function(e) numeric())
    rh <- rh[is.finite(rh)]
    if (length(rh)) max_rhat <- max(rh)
  }

  # One result row: metadata + time + nested fit
  dplyr::mutate(
    row,
    n_data_rows = nrow(data),
    threads_per_chain = threads,
    warmup = as.integer(mcmc$warmup),
    iter = as.integer(mcmc$iter),
    fit_ok = !is.null(fit_obj),
    elapsed_sec = elapsed,
    max_rhat = max_rhat,
    rds_path = rds_path,
    error = err,
    fit = list(fit_obj),
    finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )
}
