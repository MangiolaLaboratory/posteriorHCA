library(targets)
library(tidyverse)
library(targets)
library(tarchetypes)
library(glue)
library(crew)
library(crew.cluster)
case_dir <- "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC"
bench_dir <- file.path(case_dir, "targets", "formula_bench")
root_dir <- file.path(case_dir, "data", "processed", "formula_bench")
data_dir <- file.path(root_dir, "data")
fits_dir <- file.path(root_dir, "fits")
model_store_root <- "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE"
nk_store <- file.path(model_store_root, "V1_nk", "_targets")
source(file.path(bench_dir, "formulas.R"))
source(file.path(bench_dir, "settings.R"))
source(file.path(case_dir, "R", "re_vs_interaction_helpers.R"))
source(file.path(case_dir, "R", "age_smooth_helpers.R"))
source(file.path(bench_dir, "helpers.R"))
formulas_tbl <- bench_formulas_enabled()
settings_tbl <- bench_resource_settings_enabled()
mcmc <- bench_mcmc_defaults()
settings_8 <- dplyr::filter(settings_tbl, .data$cpus == 8L)
settings_16 <- dplyr::filter(settings_tbl, .data$cpus == 16L)
slurm_tmpdir_exports <- c("#SBATCH -A saigencir003", "export TMPDIR=/scratchdata1/users/a1237163/tmp/${SLURM_JOB_ID:-$$}", 
    "export TMP=\"$TMPDIR\"", "export TEMP=\"$TMPDIR\"", "mkdir -p \"$TMPDIR\"")
tar_option_set(memory = "transient", garbage_collection = 100, 
    storage = "worker", retrieval = "worker", error = "continue", 
    cue = tar_cue(mode = "thorough"), workspace_on_error = TRUE, 
    format = "rds", controller = crew_controller_group(crew_controller_slurm(name = "elastic_mini", 
        workers = 128, tasks_max = 80, seconds_idle = 300, seconds_interval = 10, 
        crashes_max = 5, options_cluster = crew_options_slurm(script_lines = slurm_tmpdir_exports, 
            memory_gigabytes_required = c(20, 40, 80), cpus_per_task = 1, 
            time_minutes = c(60 * 2, 60 * 4, 60 * 8), verbose = TRUE)), 
        crew_controller_slurm(name = "fit_8cores", workers = 128, 
            tasks_max = 20, seconds_idle = 600, seconds_interval = 10, 
            crashes_max = 5, options_cluster = crew_options_slurm(script_lines = slurm_tmpdir_exports, 
                memory_gigabytes_required = c(80, 160, 200), 
                cpus_per_task = 8, time_minutes = c(60 * 12, 
                  60 * 24, 60 * 48), verbose = TRUE)), crew_controller_slurm(name = "fit_16cores", 
            workers = 128, tasks_max = 12, seconds_idle = 600, 
            seconds_interval = 10, crashes_max = 5, options_cluster = crew_options_slurm(script_lines = slurm_tmpdir_exports, 
                memory_gigabytes_required = c(80, 160, 200), 
                cpus_per_task = 16, time_minutes = c(60 * 12, 
                  60 * 24, 60 * 48), verbose = TRUE))))
list(tar_target(root_dir_ready, {
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(root_dir, "logs"), recursive = TRUE, 
        showWarnings = FALSE)
    root_dir
}, deployment = "main"), tar_target(formulas_registry, formulas_tbl, 
    deployment = "main"), tar_target(settings_registry, settings_tbl, 
    deployment = "main"), tar_target(mcmc_registry, mcmc, deployment = "main"), 
    tar_target(gene_table, {
        dplyr::arrange(dplyr::select(dplyr::filter(readr::read_csv(file.path(case_dir, 
            "data", "times_biomarkers.csv"), show_col_types = FALSE), 
            .data$gene %in% bench_genes_default()), gene, .feature), 
            .data$gene)
    }, packages = c("readr", "dplyr"), deployment = "main"), 
    tar_target(genes, gene_table$gene, deployment = "main"), 
    tar_target(features, gene_table$.feature, deployment = "main"), 
    tar_target(formula_ids, formulas_registry$formula_id, deployment = "main"), 
    tar_target(setting_ids_8, settings_registry$setting_id[settings_registry$cpus == 
        8L], deployment = "main"), tar_target(setting_ids_16, 
        settings_registry$setting_id[settings_registry$cpus == 
            16L], deployment = "main"), tar_target(nk_model_meta, 
        bench_load_nk_model_meta(nk_store), packages = c("targets", 
            "dplyr", "stringr", "tidyselect"), resources = tar_resources(crew = tar_resources_crew("elastic_mini"))), 
    tar_target(nk_col_data, bench_load_nk_col_data(nk_store), 
        packages = c("targets", "SummarizedExperiment"), resources = tar_resources(crew = tar_resources_crew("elastic_mini"))), 
    tar_target(gene_bundle, {
        root_dir_ready
        path <- bench_extract_gene_data_file(gene = genes, feature_id = features, 
            nk_store = nk_store, nk_meta = nk_model_meta, nk_col_data = nk_col_data, 
            data_dir = data_dir, healthy_level = "Normal")
        list(gene = genes, .feature = features, path = path, 
            data = readRDS(path))
    }, pattern = map(genes, features), iteration = "list", packages = c("targets", 
        "dplyr", "stringr", "glue", "tibble", "SummarizedExperiment"), 
        resources = tar_resources(crew = tar_resources_crew("elastic_mini"))), 
    tar_target(fit_8, {
        bench_fit_one(data = gene_bundle$data, gene = gene_bundle$gene, 
            formula_id = formula_ids, setting_id = setting_ids_8, 
            formulas_tbl = formulas_registry, settings_tbl = settings_registry, 
            mcmc = mcmc_registry, fits_dir = fits_dir)
    }, pattern = cross(gene_bundle, formula_ids, setting_ids_8), 
        iteration = "vector", packages = c("brms", "cmdstanr", 
            "posterior", "dplyr", "tibble", "glue", "stringr", 
            "parallelly"), resources = tar_resources(crew = tar_resources_crew("fit_8cores"))), 
    tar_target(fit_16, {
        bench_fit_one(data = gene_bundle$data, gene = gene_bundle$gene, 
            formula_id = formula_ids, setting_id = setting_ids_16, 
            formulas_tbl = formulas_registry, settings_tbl = settings_registry, 
            mcmc = mcmc_registry, fits_dir = fits_dir)
    }, pattern = cross(gene_bundle, formula_ids, setting_ids_16), 
        iteration = "vector", packages = c("brms", "cmdstanr", 
            "posterior", "dplyr", "tibble", "glue", "stringr", 
            "parallelly"), resources = tar_resources(crew = tar_resources_crew("fit_16cores"))), 
    tar_target(timing_summary, {
        dplyr::arrange(dplyr::bind_rows(fit_8, fit_16), .data$gene, 
            .data$formula_id, .data$setting_id)
    }, packages = c("dplyr")), tar_target(timing_summary_csv, 
        {
            path <- file.path(root_dir_ready, "timing_summary.csv")
            readr::write_csv(timing_summary, path)
            path
        }, packages = c("readr"), deployment = "main", format = "file"), 
    tar_target(formula_spec_txt, {
        path <- file.path(root_dir_ready, "formula_spec_bench.txt")
        lines <- c("# formula_bench registry dump", paste0("# MCMC: warmup=", 
            mcmc_registry$warmup, ", iter=", mcmc_registry$iter), 
            paste0("# Healthy only: disease_groups_altered == Normal"), 
            paste0("# Settings: ", paste(settings_registry$setting_id, 
                collapse = ", ")), "")
        for (i in seq_len(nrow(formulas_registry))) {
            row <- formulas_registry[i, , drop = FALSE]
            lines <- c(lines, paste0("# ", row$formula_id, " — ", 
                row$label), bench_formula_text(row), "")
        }
        readr::write_lines(lines, path)
        path
    }, packages = c("readr"), deployment = "main", format = "file"), 
    tar_target(grid_manifest_csv, {
        path <- file.path(root_dir_ready, "grid_manifest.csv")
        grid <- dplyr::arrange(dplyr::mutate(tidyr::crossing(gene_table, 
            dplyr::select(formulas_registry, "formula_id", "family", 
                "k"), dplyr::select(settings_registry, "setting_id", 
                "chains", "cpus")), branch = paste(.data$gene, 
            .data$formula_id, .data$setting_id, sep = "__"), 
            expected_threads = pmax(1L, as.integer(floor(.data$cpus/.data$chains)))), 
            .data$gene, .data$formula_id, .data$setting_id)
        readr::write_csv(grid, path)
        path
    }, packages = c("dplyr", "tidyr", "readr"), deployment = "main", 
        format = "file"))
