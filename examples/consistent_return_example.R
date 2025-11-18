# Consistent Return Format Example
# ================================
# 
# This example shows the consistent list return format from get_file_ready

library(posteriorHCA)

# Example 1: Successful download
# ==============================
cat("Example 1: Successful download\n")
cat("==============================\n")

result <- get_file_ready(
  container = "V1",
  prefix = "cd4.naive",
  filename = "ENSG00000000419"
)

# Always returns a list with status and path
cat("Status:", result$status, "\n")
cat("Path:", result$path, "\n")
cat("URL:", result$url, "\n")

# Example 2: File not found
# =========================
cat("\nExample 2: File not found\n")
cat("=========================\n")

result <- get_file_ready(
  container = "V1",
  prefix = "cd4.naive",
  filename = "NONEXISTENT_GENE"
)

# Always returns a list with status and path
cat("Status:", result$status, "\n")
cat("Path:", result$path, "\n")
cat("URL:", result$url, "\n")

# Example 3: Simple usage pattern
# ==============================
cat("\nExample 3: Simple usage pattern\n")
cat("===============================\n")

download_file <- function(container, prefix, filename) {
  result <- get_file_ready(container = container, prefix = prefix, filename = filename)
  
  if (result$status == "success") {
    cat("✓ File downloaded successfully to:", result$path, "\n")
    return(result$path)
  } else if (result$status == "not_found") {
    cat("✗ File not found:", result$url, "\n")
    return(NULL)
  } else {
    cat("✗ Error:", result$error, "\n")
    return(NULL)
  }
}

# Test the simple usage
file_path <- download_file("V1", "cd4.naive", "ENSG00000000419")
if (!is.null(file_path)) {
  cat("File is ready at:", file_path, "\n")
}

# Example 4: Batch processing with consistent format
# ==================================================
cat("\nExample 4: Batch processing\n")
cat("==========================\n")

filenames <- c("ENSG00000000419", "NONEXISTENT_GENE", "ENSG00000000457")
successful_downloads <- list()
failed_downloads <- list()

for (filename in filenames) {
  result <- get_file_ready(container = "V1", prefix = "cd4.naive", filename = filename)
  
  if (result$status == "success") {
    successful_downloads[[filename]] <- result$path
    cat("✓", filename, "->", result$path, "\n")
  } else {
    failed_downloads[[filename]] <- result$status
    cat("✗", filename, "->", result$status, "\n")
  }
}

cat("\nSummary:\n")
cat("Successful downloads:", length(successful_downloads), "\n")
cat("Failed downloads:", length(failed_downloads), "\n")

cat("\nConsistent return format examples completed!\n")
