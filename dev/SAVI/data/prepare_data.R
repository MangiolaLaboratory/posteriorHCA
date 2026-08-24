library(Seurat)

# Run from the posteriorHCA project root, or set `project_root` explicitly.
project_root <- getwd()

data_dir <- file.path(project_root, "dev", "SAVI", "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

geo_files <- tibble::tribble(
  ~accession,   ~file_name,                              ~url,
  "GSE226598",  "GSE226598_SeuratObject_SAVI.rds.gz",    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE226598&format=file&file=GSE226598%5FSeuratObject%5FSAVI%2Erds%2Egz",
  "GSE226572",  "GSE226572_SeuratObject_IFN.rds.gz",     "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE226572&format=file&file=GSE226572%5FSeuratObject%5FIFN%2Erds%2Egz"
) |>
  dplyr::mutate(
    destfile = file.path(data_dir, file_name)
  )

download_geo_file_wget <- function(url, destfile, overwrite = FALSE) {
  if (file.exists(destfile) && !overwrite) {
    message("Already exists: ", destfile)
    return(invisible(destfile))
  }
  
  dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
  
  cmd <- sprintf(
    "wget -c --tries=10 --timeout=60 --read-timeout=120 -O %s %s",
    shQuote(destfile),
    shQuote(url)
  )
  
  message("Running: ", cmd)
  status <- system(cmd)
  
  if (status != 0) {
    stop("wget failed for: ", basename(destfile))
  }
  
  invisible(destfile)
}

purrr::pwalk(
  list(
    url = geo_files$url,
    destfile = geo_files$destfile
  ),
  download_geo_file_wget
)

# GEO ships these as gzip(gzip(rds)); a single gzfile() leaves an inner gzip stream.
read_geo_rds_gz <- function(path) {
  con <- gzcon(gzfile(path))
  on.exit(close(con), add = TRUE)
  readRDS(con)
}

savi_rds <- geo_files$destfile[geo_files$accession == "GSE226598"]
seurat_savi <- read_geo_rds_gz(savi_rds)

seurat_savi_pb <- AggregateExpression(
  object = seurat_savi,
  assays = "RNA",
  slot = "counts",
  group.by = c("Sample", "CellType"),
  return.seurat = TRUE
)

# colnames are Sample_CellType; Sample ids with "_" become "-" in colnames.
pb_meta <- seurat_savi@meta.data |>
  dplyr::distinct(Sample, CellType, Treatment, Category, Experiment) |>
  dplyr::mutate(
    pb_id = paste(gsub("_", "-", Sample), CellType, sep = "_")
  )

pb_ids <- colnames(seurat_savi_pb)
idx <- match(pb_ids, pb_meta$pb_id)
if (any(is.na(idx))) {
  stop(
    "Could not match pseudobulk ids to metadata: ",
    paste(head(pb_ids[is.na(idx)], 5), collapse = ", ")
  )
}

seurat_savi_pb$Sample <- pb_meta$Sample[idx]
seurat_savi_pb$CellType <- pb_meta$CellType[idx]
seurat_savi_pb$Treatment <- pb_meta$Treatment[idx]
seurat_savi_pb$Category <- pb_meta$Category[idx]
seurat_savi_pb$Experiment <- pb_meta$Experiment[idx]

savi_pb_out <- file.path(data_dir, "GSE226598_SAVI_pseudobulk_Sample_CellType.rds")
saveRDS(seurat_savi_pb, savi_pb_out)


