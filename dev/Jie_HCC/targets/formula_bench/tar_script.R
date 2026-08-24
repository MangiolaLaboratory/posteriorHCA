# formula_bench — simple readable targets pipeline
#
# Flow:
#   1. prepared_data   ← read data/prepared_healthy.rds (built once offline)
#   2. formulas        ← formulas.R
#   3. settings        ← settings.R  (chains × CPUs)
#   4. gene_features   ← hard-coded gene + Ensembl id
#   5. param_grid_8 / _16  ← crossing(...), split by setting$cpus
#   6. fit_result_8 / _16  ← map on matching Slurm controller (8 or 16 CPUs)
#
# Extend: add a row in formulas.R or settings.R, re-run. Only new grid rows fit.
#
#   module load R/4.5.3-gfbf-2025b
#   Rscript targets/formula_bench/run.R

library(targets)

tar_script({
  library(tidyverse)
  library(targets)
  library(tarchetypes)
  library(crew)
  library(crew.cluster)

  # ---- paths ----
  case_dir <- "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC"
  bench_dir <- file.path(case_dir, "targets", "formula_bench")
  root_dir <- file.path(case_dir, "data", "processed", "formula_bench")
  fits_dir <- file.path(root_dir, "fits")

  source(file.path(bench_dir, "formulas.R"))
  source(file.path(bench_dir, "settings.R"))
  source(file.path(bench_dir, "helpers.R"))

  slurm_tmpdir <- c(
    "#SBATCH -A saigencir003",
    "export TMPDIR=/scratchdata1/users/a1237163/tmp/${SLURM_JOB_ID:-$$}",
    "export TMP=\"$TMPDIR\"; export TEMP=\"$TMPDIR\"; mkdir -p \"$TMPDIR\""
  )

  tar_option_set(
    memory = "transient",
    garbage_collection = 100,
    storage = "worker",
    retrieval = "worker",
    error = "continue",
    cue = tar_cue(mode = "thorough"),
    workspace_on_error = TRUE,
    format = "rds",
    controller = crew_controller_group(
      crew_controller_slurm(
        name = "light",
        workers = 128,
        tasks_max = 40,
        seconds_idle = 300,
        crashes_max = 5,
        options_cluster = crew_options_slurm(
          script_lines = slurm_tmpdir,
          memory_gigabytes_required = c(16, 40, 80, 120),
          cpus_per_task = 1,
          time_minutes = c(60 * 2, 60 * 4, 60 * 8),
          verbose = TRUE
        )
      ),
      # One controller per cpus tier — Slurm ask matches setting$cpus
      crew_controller_slurm(
        name = "fit_8",
        workers = 128,
        tasks_max = 20,
        seconds_idle = 600,
        crashes_max = 5,
        options_cluster = crew_options_slurm(
          script_lines = slurm_tmpdir,
          memory_gigabytes_required = c(32, 80, 160, 200),
          cpus_per_task = 8,
          time_minutes = c(60 * 12, 60 * 24, 60 * 48),
          verbose = TRUE
        )
      ),
      crew_controller_slurm(
        name = "fit_16",
        workers = 128,
        tasks_max = 20,
        seconds_idle = 600,
        crashes_max = 5,
        options_cluster = crew_options_slurm(
          script_lines = slurm_tmpdir,
          memory_gigabytes_required = c(32, 80, 160, 200),
          cpus_per_task = 16,
          time_minutes = c(60 * 12, 60 * 24, 60 * 48),
          verbose = TRUE
        )
      )
    )
  )

  list(
    # ------------------------------------------------------------------
    # 1. Prepared data (built once offline — see prepare_data_once.R)
    # ------------------------------------------------------------------
    tar_target(
      root_ready,
      {
        dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(root_dir, "data"), recursive = TRUE, showWarnings = FALSE)
        root_dir
      },
      deployment = "main"
    ),

    # Hard-coded TIMES biomarkers (symbol + Ensembl id)
    tar_target(
      gene_features,
      tibble::tribble(
        ~gene,      ~.feature,
        "SPON2",    "ENSG00000159674",
        "ZFP36L2",  "ENSG00000152518",
        "ZFP36",    "ENSG00000128016",
        "VIM",      "ENSG00000026025"
      ),
      deployment = "main"
    ),

    # Pre-built healthy count table — do not re-extract from V1_nk
    tar_target(
      prepared_data_path,
      {
        path <- file.path(root_ready, "data", "prepared_healthy.rds")
        if (!file.exists(path)) {
          stop(
            "Missing ", path, ". Run once:\n",
            "  Rscript targets/formula_bench/prepare_data_once.R"
          )
        }
        path
      },
      format = "file",
      deployment = "main"
    ),

    tar_target(
      prepared_data,
      readRDS(prepared_data_path),
      deployment = "main"
    ),

    # ------------------------------------------------------------------
    # 2–4. Formula list, settings list, MCMC knobs
    # ------------------------------------------------------------------
    tar_target(formulas, bench_formulas_enabled(), deployment = "main"),
    tar_target(settings, bench_resource_settings_enabled(), deployment = "main"),
    tar_target(mcmc, bench_mcmc(), deployment = "main"),

    # ------------------------------------------------------------------
    # 5. Param grids split by cpus so Slurm cpus_per_task matches the setting
    # ------------------------------------------------------------------
    # Stable tar_group from (gene, formula_id, setting_id) so adding a formula
    # only creates new branches — finished fits are not invalidated.
    tar_target(
      param_grid_8,
      {
        tidyr::crossing(
          gene_features,
          formulas |> dplyr::rename(formula_label = "label"),
          settings |> dplyr::rename(setting_label = "label")
        ) |>
          dplyr::filter(.data$cpus == 8L) |>
          dplyr::mutate(
            tar_group = vapply(
              paste(.data$gene, .data$formula_id, .data$setting_id, sep = "|"),
              digest::digest2int,
              integer(1)
            ) |>
              abs() |>
              pmax(1L)
          ) |>
          dplyr::group_by(tar_group) |>
          targets::tar_group()
      },
      iteration = "group",
      deployment = "main",
      packages = c("tidyr", "tibble", "dplyr", "targets", "digest")
    ),

    tar_target(
      param_grid_16,
      {
        tidyr::crossing(
          gene_features,
          formulas |> dplyr::rename(formula_label = "label"),
          settings |> dplyr::rename(setting_label = "label")
        ) |>
          dplyr::filter(.data$cpus == 16L) |>
          dplyr::mutate(
            tar_group = vapply(
              paste(.data$gene, .data$formula_id, .data$setting_id, sep = "|"),
              digest::digest2int,
              integer(1)
            ) |>
              abs() |>
              pmax(1L)
          ) |>
          dplyr::group_by(tar_group) |>
          targets::tar_group()
      },
      iteration = "group",
      deployment = "main",
      packages = c("tidyr", "tibble", "dplyr", "targets", "digest")
    ),

    # ------------------------------------------------------------------
    # 6. Map fits on matching controllers (8-CPU jobs → fit_8, etc.)
    # ------------------------------------------------------------------
    tar_target(
      fit_result_8,
      fit_bench_row(param_grid_8, prepared_data, mcmc = mcmc, fits_dir = fits_dir),
      pattern = map(param_grid_8),
      iteration = "list",
      packages = c(
        "brms", "cmdstanr", "posterior", "dplyr", "tibble",
        "glue", "stringr", "parallelly"
      ),
      resources = tar_resources(crew = tar_resources_crew("fit_8"))
    ),

    tar_target(
      fit_result_16,
      fit_bench_row(param_grid_16, prepared_data, mcmc = mcmc, fits_dir = fits_dir),
      pattern = map(param_grid_16),
      iteration = "list",
      packages = c(
        "brms", "cmdstanr", "posterior", "dplyr", "tibble",
        "glue", "stringr", "parallelly"
      ),
      resources = tar_resources(crew = tar_resources_crew("fit_16"))
    ),

    # ------------------------------------------------------------------
    # Summaries
    # ------------------------------------------------------------------
    tar_target(
      timing_summary,
      {
        drop_fit <- function(x) {
          if (is.data.frame(x)) {
            return(dplyr::select(x, -dplyr::any_of(c("fit"))))
          }
          dplyr::bind_rows(
            lapply(x, \(d) dplyr::select(d, -dplyr::any_of(c("fit"))))
          )
        }
        dplyr::bind_rows(drop_fit(fit_result_8), drop_fit(fit_result_16)) |>
          dplyr::arrange(.data$gene, .data$formula_id, .data$setting_id)
      },
      packages = c("dplyr"),
      resources = tar_resources(crew = tar_resources_crew("light"))
    ),

    tar_target(
      timing_summary_csv,
      {
        path <- file.path(root_ready, "timing_summary.csv")
        readr::write_csv(timing_summary, path)
        path
      },
      format = "file",
      packages = c("readr"),
      deployment = "main"
    ),

    tar_target(
      grid_manifest_csv,
      {
        path <- file.path(root_ready, "grid_manifest.csv")
        readr::write_csv(
          tidyr::crossing(
            gene_features,
            formulas |> dplyr::select("formula_id", "family", "k"),
            settings |> dplyr::select("setting_id", "chains", "cpus")
          ) |>
            dplyr::mutate(
              expected_threads = pmax(1L, as.integer(floor(.data$cpus / .data$chains))),
              slurm_controller = paste0("fit_", .data$cpus)
            ),
          path
        )
        path
      },
      format = "file",
      packages = c("readr", "tidyr", "tibble", "dplyr"),
      deployment = "main"
    )
  )
})
