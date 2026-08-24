#!/usr/bin/env Rscript
# One-shot: build prepared_healthy.rds for formula_bench.
#
# V1_nk pseudobulk assay is HDF5-backed on a path not available here, so we
# pull counts from production gene brms fits and attach SE colData (no assay).

suppressPackageStartupMessages({
  library(targets)
  library(SummarizedExperiment)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(glue)
  library(readr)
})

case_dir <- "/home/a1237163/lab/chen/posteriorHCA/dev/Jie_HCC"
root_dir <- file.path(case_dir, "data", "processed", "formula_bench")
out_path <- file.path(root_dir, "data", "prepared_healthy.rds")
nk_store <- "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE/V1_nk/_targets"
healthy_level <- "Normal"

gene_features <- tibble::tribble(
  ~gene,      ~.feature,
  "SPON2",    "ENSG00000159674",
  "ZFP36L2",  "ENSG00000152518",
  "ZFP36",    "ENSG00000128016",
  "VIM",      "ENSG00000026025"
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

message("Loading V1_nk model meta ...")
nk_meta <- tar_meta(starts_with("estimates_chunk_"), store = nk_store) |>
  filter(.data$type == "branch") |>
  mutate(
    .feature = str_extract(.data$warnings, "(?<=Gene:___)ENSG\\d+(?=___)")
  ) |>
  filter(!is.na(.data$.feature)) |>
  select(name, .feature)

# ---- counts from production fits (one column per gene) ----
message("Loading production fit$data for ", nrow(gene_features), " genes ...")
count_list <- list()
meta0 <- NULL
for (i in seq_len(nrow(gene_features))) {
  gene <- gene_features$gene[[i]]
  feat <- gene_features$.feature[[i]]
  tname <- nk_meta |>
    filter(.data$.feature == feat) |>
    slice(1) |>
    pull(.data$name)
  if (!length(tname) || is.na(tname[[1]])) {
    stop("No estimates_chunk for ", gene, " (", feat, ")")
  }
  message("  ", gene, " ← ", tname[[1]])
  dat <- tar_read_raw(tname[[1]], store = nk_store)$brms_fit[[1]]$data
  colnames(dat) <- str_replace_all(colnames(dat), "_+", "_")
  count_list[[gene]] <- as.numeric(dat$counts)
  if (is.null(meta0)) {
    meta0 <- dat[, setdiff(names(dat), "counts"), drop = FALSE]
  } else if (nrow(dat) != nrow(meta0)) {
    stop("Row mismatch for ", gene)
  }
}

counts_mat <- do.call(rbind, count_list)
rownames(counts_mat) <- names(count_list)

# ---- attach continuous age / assay_groups / dataset_id from SE colData ----
message("Loading SE colData (no assay realize) ...")
se <- tar_read(pseudobulk_sample, store = nk_store)
cd <- as.data.frame(colData(se))
names(cd) <- str_replace_all(names(cd), "_+", "_")
stopifnot(nrow(cd) == nrow(meta0))

meta0$age_days_scaled <- as.numeric(cd$age_days_scaled)
meta0$assay_groups <- factor(cd$assay_groups)
meta0$dataset_id <- factor(cd$dataset_id)

# ---- healthy filter ----
keep <- as.character(meta0$disease_groups_altered) == healthy_level
message(glue("Healthy filter ({healthy_level}): {sum(keep)} / {length(keep)} samples"))

meta <- as_tibble(meta0[keep, , drop = FALSE])
counts_mat <- counts_mat[, keep, drop = FALSE]

meta <- meta |>
  mutate(
    across(
      any_of(c(
        "age_decade", "sex", "ethnicity_groups", "tissue_groups",
        "assay_groups", "assay_groups_altered", "dataset_id",
        "dataset_id_altered", "disease_groups_altered"
      )),
      \(x) factor(x)
    ),
    tissue_sex = interaction(tissue_groups, sex, drop = TRUE, sep = "_by_")
  )

# Keep only columns needed by formulas
keep_cols <- c(
  "offset", "age_decade", "age_days_scaled", "sex", "ethnicity_groups",
  "tissue_groups", "assay_groups", "assay_groups_altered",
  "dataset_id", "dataset_id_altered", "disease_groups_altered", "tissue_sex"
)
meta <- meta[, keep_cols, drop = FALSE]

out <- list(
  meta = meta,
  counts = counts_mat,
  features = gene_features,
  n_samples = nrow(meta),
  healthy_level = healthy_level,
  source = "production_fit_data + SE_colData"
)

saveRDS(out, out_path, compress = "xz")
message("Wrote ", out_path)
message("  n_samples = ", out$n_samples)
message("  genes = ", paste(rownames(out$counts), collapse = ", "))
message("  size = ", round(file.info(out_path)$size / 1e6, 2), " MB")
