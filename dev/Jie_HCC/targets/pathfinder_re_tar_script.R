# Targets pipeline: Pathfinder (CPU dense) for original RE + split RE.
#
# Four TIMES biomarkers x two mean formulas, simple shape.
# Draw budget matched to strong HMC: 4 chains x 1500 post-warmup = 6000 draws
#   -> Pathfinder num_paths=50, single_path_draws=120 (50*120=6000).
#
# Usage:
#   module load R/4.5.3-gfbf-2025b
#   Rscript targets/run_pathfinder_re.R
#
# Optional:
#   JUST_SCRIPT=true Rscript targets/run_pathfinder_re.R
#   NAMES=fit_pathfinder_re_ZFP36 Rscript targets/run_pathfinder_re.R

library(targets)

tar_script({
  library(tidyverse)
  library(targets)
  library(tarchetypes)
  library(glue)
  library(crew)
  library(crew.cluster)

  case_dir <- "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC"
  helpers_path <- file.path(case_dir, "R", "re_vs_interaction_helpers.R")
  model_store_root <- "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE"
  nk_store <- file.path(model_store_root, "V1_nk", "_targets")
  # Save beside the strong HMC fits so the comparison report can pick them up.
  out_dir <- file.path(case_dir, "data", "processed", "re_vs_interaction_strong")

  source(helpers_path)

  pf <- pathfinder_draw_settings(
    hmc_chains = 4L,
    hmc_warmup = 1500L,
    hmc_iter = 3000L,
    num_paths = 50L
  )
  num_paths <- pf$num_paths
  single_path_draws <- pf$single_path_draws
  threads <- 8L
  seed <- 20260720L

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
      # Pathfinder: cores=1, threads=8 within-path.
      crew_controller_slurm(
        name = "fit_pathfinder_8cores",
        workers = 8,
        tasks_max = 8,
        seconds_idle = 300,
        seconds_interval = 10,
        crashes_max = 5,
        options_cluster = crew_options_slurm(
          script_lines = slurm_tmpdir_exports,
          memory_gigabytes_required = c(40, 80, 160),
          cpus_per_task = 8,
          time_minutes = c(60 * 4, 60 * 8, 60 * 12),
          verbose = TRUE
        )
      )
    )
  )

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

  fit_grid <- tidyr::crossing(biomarkers, fit_grid_spec_pathfinder()) |>
    dplyr::mutate(branch_name = paste(gene, model, sep = "_"))

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

        data <- load_gene_training_data(.feature, gene, nk_store, nk_model_meta)

        message(glue::glue(
          "Pathfinder {gene} / {model}: paths={num_paths}, ",
          "draws/path={single_path_draws}, threads={threads}",
          gene = gene,
          model = model,
          num_paths = num_paths,
          single_path_draws = single_path_draws,
          threads = threads
        ))

        fit_obj <- fit_gene_model_pathfinder(
          data = data,
          gene = gene,
          counts_formula = counts,
          shape_formula = shape,
          num_paths = num_paths,
          single_path_draws = single_path_draws,
          threads = threads,
          seed = seed,
          tag = model
        )

        rds_path <- file.path(
          out_dir_ready,
          glue::glue(file_pattern, gene = gene)
        )
        saveRDS(fit_obj, rds_path, compress = "xz")
        message(paste("Saved", rds_path))

        n_draws <- tryCatch(
          nrow(as.data.frame(fit_obj)),
          error = function(e) NA_integer_
        )

        tibble::tibble(
          gene = gene,
          .feature = .feature,
          model = model,
          shape_variant = shape_variant,
          algorithm = "pathfinder",
          counts_formula = counts,
          shape_formula = shape,
          rds_path = rds_path,
          n_data_rows = nrow(data),
          num_paths = as.integer(num_paths),
          single_path_draws = as.integer(single_path_draws),
          n_draws = as.integer(n_draws),
          max_rhat = suppressWarnings({
            rh <- tryCatch(posterior::rhat(fit_obj), error = function(e) NA_real_)
            if (all(is.na(rh))) NA_real_ else max(rh, na.rm = TRUE)
          })
        )
      },
      packages = c(
        "brms", "posterior", "glue", "dplyr", "tibble", "stringr", "targets"
      ),
      resources = tar_resources(crew = tar_resources_crew("fit_pathfinder_8cores"))
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
        path <- file.path(out_dir_ready, "fit_summary_pathfinder.csv")
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
        pf_local <- pathfinder_draw_settings()
        path <- file.path(out_dir_ready, "formula_spec_pathfinder.txt")
        readr::write_lines(
          c(
            glue::glue(
              "# Pathfinder CPU dense: num_paths={pf_local$num_paths}, ",
              "single_path_draws={pf_local$single_path_draws}, ",
              "total_draws={pf_local$total_draws} ",
              "(matched to HMC 4x1500 post-warmup = {pf_local$target_draws})"
            ),
            "",
            "# 1. Original random-effect counts + simple shape",
            specs$re$counts,
            specs$re$shape_simple,
            "",
            "# 2. Split random-effect counts + simple shape",
            specs$re_split$counts,
            specs$re_split$shape_simple
          ),
          path
        )
        path
      },
      packages = c("readr", "glue"),
      deployment = "main"
    )
  )
})
