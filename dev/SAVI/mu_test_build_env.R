# Build the minimal environment that mu_test.R depends on, and cache it.
# Reproduces the `merge-savi-reference` chunk of
# SAVI_case_study_universal_baseline.qmd just far enough to run the per-gene
# mean/SE tests (edgeR cohort mu vs brms healthy baseline draws).
#
# Run with the project's conda R:
#   PATH=/home/a1237163/miniconda3/envs/R_env/bin:$PATH \
#   LD_LIBRARY_PATH=/home/a1237163/miniconda3/envs/R_env/lib:$LD_LIBRARY_PATH \
#   Rscript mu_test_build_env.R

suppressMessages({
  library(SummarizedExperiment)
  library(S4Vectors)
  library(SeuratObject)
  library(targets)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
})

data_dir <- "/home/a1237163/lab/chen/posteriorHCA/dev/SAVI/data"
h5ad_path <- "/home/a1237163/lab/Mangiola_ImmuneAtlas/zenodo/pseudobulk_se.h5ad"
model_path <- file.path(
  "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/NEW_CELL_TYPE",
  "V1_monocytic"
)
# data_dir is read-only on this node; cache to a writable location.
cache_rds <- "/home/a1237163/mu_test_se_monocytic.rds"

# --- SAVI pseudobulk (gene-symbol rows) -> Ensembl SummarizedExperiment --------
SAVI_pseudobulk_sobj <-
  readRDS(file.path(data_dir, "GSE226598_SAVI_pseudobulk_Sample_CellType.rds")) |>
  subset(subset = CellType == "17. Disease-associated monocytes")

counts_sym <- as.matrix(
  SeuratObject::LayerData(SAVI_pseudobulk_sobj, assay = "RNA", layer = "counts")
)

sym <- rownames(counts_sym)
ensg <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = sym, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first"
)
keep <- !is.na(ensg) & !duplicated(ensg)
ensg <- ensg[keep]

counts_ens <- counts_sym[keep, , drop = FALSE]
rownames(counts_ens) <- ensg

cd_savi <- as.data.frame(SAVI_pseudobulk_sobj@meta.data)
rownames(cd_savi) <- colnames(counts_ens)
cd_savi$sample_role <- "savi_pseudobulk"
cd_savi$offset <- 0L
if (!"cell_type_unified_ensemble" %in% names(cd_savi)) {
  cd_savi$cell_type_unified_ensemble <- cd_savi$CellType
}

se_savi_ensembl <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = counts_ens),
  colData = S4Vectors::DataFrame(cd_savi, row.names = colnames(counts_ens)),
  rowData = S4Vectors::DataFrame(gene_symbol = sym[keep], row.names = ensg)
)

# --- Healthy atlas reference column (offset == 0) ------------------------------
pseudobulk_sample <- targets::tar_read(
  pseudobulk_sample,
  store = file.path(model_path, "_targets")
)

# Point the HDF5 backend at a path that exists on this node.
pseudobulk_sample@assays@data@listData[["counts"]]@seed@seed@seed@filepath <- h5ad_path

ref_sample <- pseudobulk_sample[, SummarizedExperiment::colData(pseudobulk_sample)$offset == 0]
SummarizedExperiment::colData(ref_sample)$sample_role <- "reference"

# --- rbind-like merge (intersect genes + shared assays) ------------------------
merge_summarized_experiments <- function(...) {
  ses <- list(...)
  assay_names <- Reduce(intersect, lapply(ses, SummarizedExperiment::assayNames))
  genes <- Reduce(intersect, lapply(ses, rownames))
  ses <- lapply(ses, function(se) se[genes, , drop = FALSE])

  cds <- lapply(ses, function(se) as.data.frame(SummarizedExperiment::colData(se)))
  all_nms <- Reduce(union, lapply(cds, names))
  cds <- lapply(cds, function(cd) {
    for (nm in setdiff(all_nms, names(cd))) cd[[nm]] <- NA
    cd[, all_nms, drop = FALSE]
  })

  assays <- stats::setNames(
    lapply(assay_names, function(an) {
      do.call(base::cbind, lapply(ses, function(se) {
        base::as.matrix(SummarizedExperiment::assay(se, an))
      }))
    }),
    assay_names
  )

  SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    colData = S4Vectors::DataFrame(
      do.call(rbind, cds),
      row.names = unlist(lapply(ses, colnames), use.names = FALSE)
    ),
    rowData = SummarizedExperiment::rowData(ses[[1]])
  )
}

se_monocytic <- merge_summarized_experiments(se_savi_ensembl, ref_sample)

ref_sample_name <-
  se_monocytic |>
  SummarizedExperiment::colData() |>
  as.data.frame() |>
  dplyr::filter(sample_role == "reference") |>
  rownames()

cat("se_monocytic:", nrow(se_monocytic), "genes x", ncol(se_monocytic), "samples\n")
cat("Category:\n"); print(table(SummarizedExperiment::colData(se_monocytic)$Category, useNA = "always"))
cat("ref_sample_name:", ref_sample_name, "\n")

saveRDS(
  list(se_monocytic = se_monocytic, ref_sample_name = ref_sample_name),
  cache_rds, compress = "xz"
)
cat("Saved ->", cache_rds, "\n")
