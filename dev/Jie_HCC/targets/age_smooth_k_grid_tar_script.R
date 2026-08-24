# Targets pipeline: age-smooth k-grid (four TIMES biomarkers × k ∈ {3,4,5}).
#
# Continuous-age thin-plate smooth by tissue×sex (same formula as
# HCC_age_smooth_one_gene_pilot.R), fit on **healthy** samples only
# (disease_groups_altered == "Normal").
# MCMC: 4 chains, 400 warmup / 900 iter, 8 cores, threading(2), control = NULL.
#
# Usage:
#   module load R/4.5.3-gfbf-2025b
#   Rscript targets/run_age_smooth_k_grid.R
#
# Optional:
#   JUST_SCRIPT=true Rscript targets/run_age_smooth_k_grid.R
#   NAMES=fit_SPON2_k3 Rscript targets/run_age_smooth_k_grid.R

library(targets)

tar_script({
  library(tidyverse)
  library(targets)
  library(tarchetypes)
  library(glue)
  library(crew)
  library(crew.cluster)

  case_dir <- "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC"
  helpers_re <- file.path(case_dir, "R", "re_vs_interaction_helpers.R")
  helpers_age <- file.path(case_dir, "R", "age_smooth_helpers.R")
  model_store_root <- "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE"
  nk_store <- file.path(model_store_root, "V1_nk", "_targets")
  out_dir <- file.path(case_dir, "data", "processed", "age_smooth_k_grid_pilot")

  chains <- 4L
  warmup <- 400L
  iter <- 900L
  fit_cores <- 8L
  threads_per_chain <- 2L
  seed <- 20260720L

  # CmdStan compile temp on scratch (node /tmp often too small).
  slurm_tmpdir_exports <- c(
    "#SBATCH -A saigencir003",
    "export TMPDIR=/scratchdata1/users/a1237163/tmp/${SLURM_JOB_ID:-$$}",
    "export TMP=\"$TMPDIR\"",
    "export TEMP=\"$TMPDIR\"",
    "mkdir -p \"$TMPDIR\""
  )

  tar_option_set(
    memory = "transient",
    garbage_collection = 100,
    storage = "worker",
    retrieval = "worker",
    error = "continue",
    cue = tar_cue(mode = "never"),
    workspace_on_error = TRUE,
    # qs is not installed in this R module; use rds for store artefacts.
    format = "rds",
    controller = crew_controller_group(
      crew_controller_slurm(
        name = "elastic_mini",
        workers = 50,
        tasks_max = 100,
        seconds_idle = 300,
        seconds_interval = 10,
        crashes_max = 5,
        options_cluster = crew_options_slurm(
          script_lines = slurm_tmpdir_exports,
          memory_gigabytes_required = c(20, 40, 80),
          cpus_per_task = 1,
          time_minutes = c(60 * 2, 60 * 4, 60 * 8),
          verbose = TRUE
        )
      ),
      # 8 CPUs = 4 chains × 2 within-chain threads (matches one-gene pilot).
      crew_controller_slurm(
        name = "fit_8cores",
        workers = 6,
        tasks_max = 12,
        seconds_idle = 600,
        seconds_interval = 10,
        crashes_max = 5,
        options_cluster = crew_options_slurm(
          script_lines = slurm_tmpdir_exports,
          memory_gigabytes_required = c(80, 160, 200),
          cpus_per_task = 8,
          time_minutes = c(60 * 12, 60 * 24, 60 * 48),
          verbose = TRUE
        )
      )
    )
  )

  source(helpers_re)
  source(helpers_age)

  biomarkers <- readr::read_csv(
    file.path(case_dir, "data", "times_biomarkers.csv"),
    show_col_types = FALSE
  ) |>
    dplyr::filter(.data$gene %in% pilot_biomarker_genes()) |>
    dplyr::select(gene, .feature)

  if (nrow(biomarkers) != length(pilot_biomarker_genes())) {
    stop(
      "Expected ", length(pilot_biomarker_genes()), " biomarker genes; found ",
      nrow(biomarkers), "."
    )
  }

  fit_grid <- tidyr::crossing(biomarkers, fit_grid_spec_age_smooth()) |>
    dplyr::mutate(branch_name = paste0(gene, "_k", k))

  load_gene_fit_data <- function(feature_id, gene_symbol, store, meta) {
    target_name <- meta |>
      dplyr::filter(.data$.feature == feature_id) |>
      dplyr::slice(1) |>
      dplyr::pull(.data$name)
    if (!length(target_name) || is.na(target_name[[1]])) {
      stop(glue::glue(
        "No estimates_chunk branch for {feature_id} ({gene_symbol}) in {store}"
      ))
    }
    targets::tar_read_raw(target_name[[1]], store = store)$brms_fit[[1]]$data
  }

  mapped_fits <- tar_map(
    values = fit_grid,
    names = "branch_name",
    tar_target(
      fit,
      {
        tmpdir <- file.path(
          "/scratchdata1/users/a1237163/tmp",
          Sys.getenv("SLURM_JOB_ID", unset = as.character(Sys.getpid()))
        )
        dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
        Sys.setenv(TMPDIR = tmpdir, TMP = tmpdir, TEMP = tmpdir)

        raw_data <- load_gene_fit_data(.feature, gene, nk_store, nk_model_meta)
        data <- prepare_age_smooth_data(
          raw_data, nk_col_data,
          healthy_only = TRUE,
          healthy_level = "Normal"
        )

        message(glue::glue(
          "Fitting {gene} age-smooth k={k} (healthy Normal only, n={n}): ",
          "chains={chains}, threads={threads_per_chain}, ",
          "warmup={warmup}, iter={iter}, cores={fit_cores}, control=NULL",
          gene = gene,
          k = k,
          n = nrow(data),
          chains = chains,
          threads_per_chain = threads_per_chain,
          warmup = warmup,
          iter = iter,
          fit_cores = fit_cores
        ))

        t0 <- Sys.time()
        fit_obj <- fit_age_smooth_gene(
          data = data,
          gene = gene,
          k = as.integer(k),
          chains = chains,
          cores = fit_cores,
          threads_per_chain = threads_per_chain,
          warmup = warmup,
          iter = iter,
          seed = seed,
          control = NULL
        )
        elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

        rds_path <- file.path(
          out_dir_ready,
          glue::glue(file_pattern, gene = gene, k = k)
        )
        saveRDS(fit_obj, rds_path, compress = "xz")
        message(paste("Saved", rds_path, "in", round(elapsed, 1), "s"))

        rh <- tryCatch(posterior::rhat(fit_obj), error = function(e) NA_real_)
        rh <- rh[is.finite(rh)]
        np <- tryCatch(brms::nuts_params(fit_obj), error = function(e) NULL)
        pct_div <- if (is.null(np)) {
          NA_real_
        } else {
          100 * mean(np$Value[np$Parameter == "divergent__"])
        }

        tibble::tibble(
          gene = gene,
          .feature = .feature,
          k = as.integer(k),
          model = model,
          rds_path = rds_path,
          n_data_rows = nrow(data),
          n_tissue_sex = nlevels(data$tissue_sex),
          chains = chains,
          warmup = warmup,
          iter = iter,
          elapsed_sec = elapsed,
          max_rhat = if (length(rh)) max(rh) else NA_real_,
          n_rhat_gt_1_05 = if (length(rh)) as.integer(sum(rh > 1.05)) else NA_integer_,
          pct_divergent = pct_div
        )
      },
      packages = c(
        "brms", "posterior", "glue", "dplyr", "tibble", "stringr", "targets"
      ),
      resources = tar_resources(crew = tar_resources_crew("fit_8cores"))
    )
  )

  list(
    tar_target(
      out_dir_ready,
      {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        out_dir
      },
      deployment = "main"
    ),

    tar_target(
      nk_model_meta,
      {
        targets::tar_meta(
          starts_with("estimates_chunk_"),
          store = nk_store
        ) |>
          dplyr::filter(.data$type == "branch") |>
          dplyr::mutate(
            .feature = stringr::str_extract(
              .data$warnings,
              "(?<=Gene:___)ENSG\\d+(?=___)"
            )
          ) |>
          dplyr::filter(!is.na(.data$.feature)) |>
          dplyr::select(name, .feature)
      },
      packages = c("targets", "dplyr", "stringr"),
      resources = tar_resources(crew = tar_resources_crew("elastic_mini"))
    ),

    tar_target(
      nk_col_data,
      {
        se <- targets::tar_read(pseudobulk_sample, store = nk_store)
        cd <- as.data.frame(SummarizedExperiment::colData(se))
        # Keep only columns needed for age-smooth prep (smaller qs artefact).
        keep <- intersect(
          c("age_days_scaled", "age_days", "assay_groups", "dataset_id"),
          names(cd)
        )
        cd[, keep, drop = FALSE]
      },
      packages = c("targets", "SummarizedExperiment"),
      resources = tar_resources(crew = tar_resources_crew("elastic_mini"))
    ),

    mapped_fits,

    tar_combine(
      fit_summary,
      mapped_fits,
      command = dplyr::bind_rows(!!!.x) |>
        dplyr::arrange(gene, k),
      packages = c("dplyr")
    ),

    tar_target(
      fit_summary_csv,
      {
        path <- file.path(out_dir_ready, "age_smooth_k_grid_results.csv")
        readr::write_csv(fit_summary, path)
        path
      },
      packages = c("readr"),
      deployment = "main"
    ),

    tar_target(
      formula_spec_txt,
      {
        path <- file.path(out_dir_ready, "formula_spec_age_smooth_k_grid.txt")
        lines <- c(
          "# Age-smooth k-grid: 4 genes × k ∈ {3,4,5}",
          "# Data: healthy only (disease_groups_altered == \"Normal\")",
          "# MCMC: chains=4, warmup=400, iter=900, cores=8, threads=2, control=NULL",
          ""
        )
        for (kk in age_smooth_ks_default()) {
          lines <- c(
            lines,
            paste0("# k = ", kk),
            age_smooth_formula_text(kk),
            ""
          )
        }
        readr::write_lines(lines, path)
        path
      },
      packages = c("readr"),
      deployment = "main"
    )
  )
})
