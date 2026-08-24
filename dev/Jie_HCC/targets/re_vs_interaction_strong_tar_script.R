# Targets pipeline: RE vs interaction vs RE-uncorrelated (simple shape).
#
# Four TIMES biomarkers x three mean formulas, all with the revised *simple*
# shape. MCMC: 4 chains, 1500 warmup, 3000 total iterations, 4 cores per fit.
#
# Usage (login node that can submit Slurm jobs via crew):
#   module load R/4.5.3-gfbf-2025b
#   Rscript targets/run_re_vs_interaction_strong.R
#
# Or manually:
#   cd /hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/re_vs_interaction_strong
#   Rscript /path/to/re_vs_interaction_strong_tar_script.R   # writes _targets.R
#   R -e 'targets::tar_make()'

library(targets)

tar_script({
  library(tidyverse)
  library(targets)
  library(tarchetypes)
  library(glue)
  library(crew)
  library(crew.cluster)

  #-----------------------#
  # Paths / constants -----
  #-----------------------#
  case_dir <- "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC"
  helpers_path <- file.path(case_dir, "R", "re_vs_interaction_helpers.R")
  model_store_root <- "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE"
  nk_store <- file.path(model_store_root, "V1_nk", "_targets")
  out_dir <- file.path(case_dir, "data", "processed", "re_vs_interaction_strong")

  chains <- 4L
  warmup <- 1500L
  iter <- 3000L
  # Parallel chains; within-chain threads set separately (4 x 8 = 32 CPUs).
  fit_cores <- 4L
  threads_per_chain <- 8L
  adapt_delta <- 0.95
  max_treedepth <- 12L
  seed <- 42L

  # CmdStan writes large temp objects under TMPDIR during compile. Node-local
  # /tmp is often too small for the RE random-slopes models, so point workers
  # at user scratch (see failed fit_ZFP36*_re: "ran out of disk space").
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
    # Keep completed fits; only explicitly invalidated targets re-run.
    cue = tar_cue(mode = "never"),
    workspace_on_error = TRUE,
    format = "qs",
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
      # One fit per Slurm task: 32 CPUs = 4 chains x 8 within-chain threads.
      crew_controller_slurm(
        name = "fit_32cores",
        workers = 6,
        tasks_max = 12,
        seconds_idle = 600,
        seconds_interval = 10,
        crashes_max = 5,
        options_cluster = crew_options_slurm(
          script_lines = slurm_tmpdir_exports,
          memory_gigabytes_required = c(80, 160, 200),
          cpus_per_task = 32,
          time_minutes = c(60 * 12, 60 * 24, 60 * 48),
          verbose = TRUE
        )
      )
    )
  )

  source(helpers_path)

  #-----------------------#
  # Static gene x formula grid (4 genes x 3 models = 12 fits) -----
  #-----------------------#
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

  fit_grid <- tidyr::crossing(biomarkers, fit_grid_spec_strong()) |>
    dplyr::mutate(branch_name = paste(gene, model, sep = "_"))

  # Local helper: pull training data for one gene from the production V1_nk store.
  load_gene_training_data <- function(feature_id, gene_symbol, store, meta) {
    target_name <- meta |>
      dplyr::filter(.data$.feature == feature_id) |>
      dplyr::slice(1) |>
      dplyr::pull(.data$name)
    if (!length(target_name) || is.na(target_name[[1]])) {
      stop(glue::glue(
        "No estimates_chunk branch for {feature_id} ({gene_symbol}) in {store}"
      ))
    }
    fit_prod <- targets::tar_read_raw(target_name[[1]], store = store)$brms_fit[[1]]
    data <- droplevels(fit_prod$data)
    colnames(data) <- stringr::str_replace_all(colnames(data), "_+", "_")
    data
  }

  #-----------------------#
  # Pipeline -----
  #-----------------------#
  mapped_fits <- tar_map(
    values = fit_grid,
    names = "branch_name",
    tar_target(
      fit,
      {
        # Belt-and-suspenders: ensure CmdStan compile temp is on scratch.
        tmpdir <- file.path(
          "/scratchdata1/users/a1237163/tmp",
          Sys.getenv("SLURM_JOB_ID", unset = as.character(Sys.getpid()))
        )
        dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
        Sys.setenv(TMPDIR = tmpdir, TMP = tmpdir, TEMP = tmpdir)

        data <- load_gene_training_data(.feature, gene, nk_store, nk_model_meta)

        message(glue::glue(
          "Fitting {gene} / {model} (simple shape): ",
          "chains={chains}, threads_per_chain={threads_per_chain}, ",
          "warmup={warmup}, iter={iter}, chain_cores={fit_cores}",
          gene = gene,
          model = model,
          chains = chains,
          threads_per_chain = threads_per_chain,
          warmup = warmup,
          iter = iter,
          fit_cores = fit_cores
        ))

        fit_obj <- fit_gene_model(
          data = data,
          gene = gene,
          counts_formula = counts,
          shape_formula = shape,
          cores = fit_cores,
          chains = chains,
          threads_per_chain = threads_per_chain,
          warmup = warmup,
          iter = iter,
          adapt_delta = adapt_delta,
          max_treedepth = max_treedepth,
          seed = seed,
          tag = paste(model, "simple", sep = "/")
        )

        rds_path <- file.path(
          out_dir_ready,
          glue::glue(file_pattern, gene = gene)
        )
        saveRDS(fit_obj, rds_path, compress = "xz")
        message(paste("Saved", rds_path))

        # Return metadata only; the brms fit lives in rds_path (and as the
        # side-effect artefact for the Quarto report). Avoids tar_combine
        # loading 12 full fits into memory.
        tibble::tibble(
          gene = gene,
          .feature = .feature,
          model = model,
          shape_variant = shape_variant,
          counts_formula = counts,
          shape_formula = shape,
          rds_path = rds_path,
          n_data_rows = nrow(data),
          chains = chains,
          warmup = warmup,
          iter = iter,
          max_rhat = suppressWarnings(max(posterior::rhat(fit_obj), na.rm = TRUE))
        )
      },
      packages = c(
        "brms", "posterior", "glue", "dplyr", "tibble", "stringr", "targets"
      ),
      resources = tar_resources(crew = tar_resources_crew("fit_32cores"))
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

    mapped_fits,

    tar_combine(
      fit_summary,
      mapped_fits,
      command = dplyr::bind_rows(!!!.x) |>
        dplyr::arrange(gene, model),
      packages = c("dplyr")
    ),

    tar_target(
      fit_summary_csv,
      {
        path <- file.path(out_dir_ready, "fit_summary_strong.csv")
        readr::write_csv(fit_summary, path)
        path
      },
      packages = c("readr"),
      deployment = "main"
    ),

    tar_target(
      formula_spec_txt,
      {
        specs <- formula_specs()
        path <- file.path(out_dir_ready, "formula_spec_strong.txt")
        readr::write_lines(
          c(
            "# MCMC: chains=4, warmup=1500, iter=3000, cores=4",
            "",
            "# 1. Random-effect counts (correlated tissue slopes) + simple shape",
            specs$re$counts,
            specs$re$shape_simple,
            "",
            "# 2. Interaction counts + simple shape",
            specs$interaction$counts,
            specs$interaction$shape_simple,
            "",
            "# 3. Random-effect counts (uncorrelated tissue slopes, `||`) + simple shape",
            specs$re_uncorrelated$counts,
            specs$re_uncorrelated$shape_simple
          ),
          path
        )
        path
      },
      packages = c("readr"),
      deployment = "main"
    )
  )
})
