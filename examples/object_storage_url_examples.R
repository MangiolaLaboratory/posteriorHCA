# Object Storage URL Construction Examples
# ======================================
# 
# This file demonstrates how to construct direct download URLs for Nectar object storage
# instead of using Swift CLI commands.

library(posteriorHCA)

# Example 1: Basic URL Construction Pattern
# =========================================
#
# URL Structure: https://object-store.rc.nectar.org.au/v1/{AUTH_TOKEN}/{CONTAINER}/{OBJECT_PATH}
#
# Your example URL:
# https://object-store.rc.nectar.org.au/v1/AUTH_b0a86a29c8b74630aac35f471cfe1396/V1/cd4.naive/ENSG00000000419
#
# Breakdown:
# - Base URL: https://object-store.rc.nectar.org.au/v1
# - Auth Token: AUTH_b0a86a29c8b74630aac35f471cfe1396
# - Container: V1
# - Object Path: cd4.naive/ENSG00000000419

# Example 2: Using the URL Construction Helper Function
# =====================================================

# Construct URL for your example
auth_token <- "AUTH_b0a86a29c8b74630aac35f471cfe1396"
container <- "V1"
object_path <- "cd4.naive/ENSG00000000419"

url <- construct_object_storage_url(auth_token, container, object_path)
print(paste("Constructed URL:", url))

# Example 3: Different URL Patterns
# =================================

# Pattern 1: Different container and prefix
url1 <- construct_object_storage_url(
  auth_token = "AUTH_b0a86a29c8b74630aac35f471cfe1396",
  container = "V1", 
  object_path = "cd8.memory/ENSG00000000419"
)

# Pattern 2: Different auth token (like in get_model_ready)
url2 <- construct_object_storage_url(
  auth_token = "AUTH_06d6e008e3e642da99d806ba3ea629c5",
  container = "testing",
  object_path = "estimates_age_bins___L2.rds"
)

print(paste("CD8 Memory URL:", url1))
print(paste("Model URL:", url2))

# Example 4: Using the Updated Download Function
# =============================================

# Download a file using direct URL (no Swift CLI needed)
tryCatch({
  local_path <- get_file_ready(
    container = "V1",
    prefix = "cd4.naive",
    filename = "ENSG00000000419"
  )
  print(paste("Downloaded file to:", local_path))
}, error = function(e) {
  print(paste("Download failed:", e$message))
})

# Example 5: Manual URL Construction for Different Scenarios
# ==========================================================

# Scenario 1: Different cell types
cell_types <- c("cd4.naive", "cd4.memory", "cd8.naive", "cd8.memory", "b.cells")
gene_id <- "ENSG00000000419"

for (cell_type in cell_types) {
  url <- construct_object_storage_url(
    auth_token = "AUTH_b0a86a29c8b74630aac35f471cfe1396",
    container = "V1",
    object_path = paste(cell_type, gene_id, sep = "/")
  )
  print(paste(cell_type, "URL:", url))
}

# Scenario 2: Different genes for the same cell type
genes <- c("ENSG00000000419", "ENSG00000000457", "ENSG00000000460")
cell_type <- "cd4.naive"

for (gene in genes) {
  url <- construct_object_storage_url(
    auth_token = "AUTH_b0a86a29c8b74630aac35f471cfe1396",
    container = "V1",
    object_path = paste(cell_type, gene, sep = "/")
  )
  print(paste("Gene", gene, "URL:", url))
}

# Example 6: Batch Download Function
# ==================================

download_multiple_files <- function(container, cell_types, genes, 
                                   cache_directory = get_default_cache_dir()) {
  results <- list()
  
  for (cell_type in cell_types) {
    for (gene in genes) {
      tryCatch({
        local_path <- get_file_ready(
          container = container,
          prefix = cell_type,
          filename = gene,
          cache_directory = cache_directory
        )
        results[[paste(cell_type, gene, sep = "_")]] <- local_path
        print(paste("Successfully downloaded:", cell_type, gene))
      }, error = function(e) {
        print(paste("Failed to download:", cell_type, gene, "-", e$message))
        results[[paste(cell_type, gene, sep = "_")]] <- NULL
      })
    }
  }
  
  return(results)
}

# Example usage of batch download
# cell_types <- c("cd4.naive", "cd8.naive")
# genes <- c("ENSG00000000419", "ENSG00000000457")
# results <- download_multiple_files("V1", cell_types, genes)

print("Object storage URL construction examples completed!")
