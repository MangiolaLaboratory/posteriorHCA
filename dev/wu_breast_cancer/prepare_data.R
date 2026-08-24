library(Matrix)
library(tidySingleCellExperiment)
library(Seurat)

# Run from the posteriorHCA project root, or set `project_root` explicitly.
project_root <- getwd()
data_dir <- file.path(project_root, "dev", "wu_breast_cancer", "data")
sc_dir <- file.path(data_dir, "Wu_etal_2021_BRCA_scRNASeq")
tar_file <- file.path(data_dir, "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz")

if (!dir.exists(sc_dir)) {
  url <- paste0(
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE176078&format=file&",
    "file=GSE176078%5FWu%5Fetal%5F2021%5FBRCA%5FscRNASeq%2Etar%2Egz"
  )

  if (!file.exists(tar_file)) {
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    download.file(
      url = url,
      destfile = tar_file,
      mode = "wb",
      method = "auto"
    )
  }

  untar(tarfile = tar_file, exdir = data_dir)
}

counts <- Matrix::readMM(file.path(sc_dir, "count_matrix_sparse.mtx"))
genes <- readLines(file.path(sc_dir, "count_matrix_genes.tsv"))
barcodes <- readLines(file.path(sc_dir, "count_matrix_barcodes.tsv"))

rownames(counts) <- genes
colnames(counts) <- barcodes
# counts <- as(counts, "CsparseMatrix")
# storage.mode(counts@x) <- "integer"

metadata <- read.csv(
  file.path(sc_dir, "metadata.csv"),
  row.names = 1,
  check.names = FALSE
)
metadata <- metadata[colnames(counts), , drop = FALSE]

sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = metadata
)

sce_cancer_epithelial <-
  sce %>% filter(celltype_major == "Cancer Epithelial")


library(BiocParallel)

sce_cancer_epithelial_agg <-
  sce_cancer_epithelial |>
  scuttle::aggregateAcrossCells(
    ids = colData(sce_cancer_epithelial)[, c("orig.ident", "celltype_major")],
    statistics = "sum",
    BPPARAM = MulticoreParam(20)
  )

pb_counts <- SummarizedExperiment::assay(sce_cancer_epithelial_agg, "counts")
agg_colData <- as.data.frame(SummarizedExperiment::colData(sce_cancer_epithelial_agg))

subtype_by_patient <- colData(sce_cancer_epithelial) |>
  as.data.frame() |>
  dplyr::distinct(orig.ident, subtype)

sample_ids <- agg_colData$orig.ident
colnames(pb_counts) <- sample_ids
pb_colData <- DataFrame(
  orig.ident = sample_ids,
  celltype_major = "Cancer Epithelial",
  subtype = subtype_by_patient$subtype[
    match(sample_ids, subtype_by_patient$orig.ident)
  ],
  row.names = sample_ids
)

se_cancer_epithelial <- SummarizedExperiment(
  assays = list(counts = pb_counts),
  colData = pb_colData
)

source("/hpcfs/groups/phoenix-hpc-mangiola_laboratory/chen/HPC_posterior/Age_Clock/report/lung_edgeR_figs/keep_identify_abundant_alt.R")



se_cancer_epithelial_ab <-
  se_cancer_epithelial %>% 
  identify_abundant(
    formula_design = ~ subtype
  ) %>% 
  scale_abundance(method = 'TMMwsp') %>%
  reduce_dimensions(method = "PCA",.abundance = 'scaled_counts', .dims = 10 )

se_cancer_epithelial_ab %>%
  colData() %>%
  as.data.frame() %>% 
  ggplot(aes(x = PC1, y = PC2, colour = subtype)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal()


sce_out <- file.path(data_dir, "Wu_etal_2021_BRCA_se_cancer_epithelial.rds")
saveRDS(se_cancer_epithelial, sce_out)

sobj <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

sobj_cancer_epithelial <-
  sobj %>% filter(celltype_major == "Cancer Epithelial")

sobj_pb_cancer_epithelial <-
  AggregateExpression(
    object = sobj_cancer_epithelial,
    assays = "RNA",
    slot = "counts",
    group.by = c("orig.ident", "celltype_major"),
    return.seurat = TRUE
  )

# AggregateExpression only retains group.by columns in meta.data.
subtype_by_patient <- sobj_cancer_epithelial@meta.data |>
  dplyr::distinct(orig.ident, subtype)
sobj_pb_cancer_epithelial$subtype <- subtype_by_patient$subtype[
  match(sobj_pb_cancer_epithelial$orig.ident, subtype_by_patient$orig.ident)
]

sobj_pb_cancer_epithelial <-
  sobj_pb_cancer_epithelial %>%
  NormalizeData(assay = "RNA") %>%
  FindVariableFeatures(assay = "RNA") %>%
  ScaleData(assay = "RNA") %>%
  RunPCA(assay = "RNA", npcs = 10)

DimPlot(
  sobj_pb_cancer_epithelial, reduction = 'pca', group.by = 'subtype'
)
  

sobj_out <- file.path(data_dir, "Wu_etal_2021_BRCA_scRNASeq_seurat.rds")
saveRDS(sobj, sobj_out)

sobj
