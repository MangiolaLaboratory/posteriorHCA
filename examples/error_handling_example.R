# Error Handling Example for get_file_ready
# ==========================================
# 
# This example shows how to handle different return types from get_file_ready

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

cat("Status:", result$status, "\n")
if (result$status == "success") {
  cat("File path:", result$path, "\n")
} else {
  cat("Error:", result$error, "\n")
}

# Example 2: File not found (404 error)
# ====================================
cat("\nExample 2: File not found\n")
cat("=========================\n")

result <- get_file_ready(
  container = "V1",
  prefix = "cd4.naive",
  filename = "NONEXISTENT_GENE"
)

cat("Status:", result$status, "\n")
if (result$status == "not_found") {
  cat("File not found at URL:", result$url, "\n")
} else if (result$status == "success") {
  cat("File path:", result$path, "\n")
} else {
  cat("Error:", result$error, "\n")
}

# Example 3: Robust error handling function
# ========================================
cat("\nExample 3: Robust error handling\n")
cat("=================================\n")

safe_download <- function(container, prefix, filename) {
  result <- get_file_ready(
    container = container,
    prefix = prefix,
    filename = filename
  )
  
  # Always returns a list with status and path
  return(result)
}

# Test the robust function
result <- safe_download("V1", "cd4.naive", "ENSG00000000419")
cat("Robust result status:", result$status, "\n")

# Example 4: Batch download with error handling
# ============================================
cat("\nExample 4: Batch download with error handling\n")
cat("==============================================\n")

batch_download <- function(container, prefix, filenames) {
  results <- list()
  
  for (filename in filenames) {
    cat("Downloading:", filename, "\n")
    result <- get_file_ready(container = container, prefix = prefix, filename = filename)
    
    # Always returns a list with status and path
    results[[filename]] <- result
  }
  
  return(results)
}

# Test batch download
filenames <- c("ENSG00000000419", "NONEXISTENT_GENE", "ENSG00000000457")
batch_results <- batch_download("V1", "cd4.naive", filenames)

# Summary of batch results
cat("\nBatch download summary:\n")
for (filename in names(batch_results)) {
  result <- batch_results[[filename]]
  cat(filename, ":", result$status, "\n")
  if (result$status == "success") {
    cat("  Path:", result$path, "\n")
  } else if (result$status == "not_found") {
    cat("  File not found\n")
  } else {
    cat("  Error:", result$error, "\n")
  }
}

cat("\nError handling examples completed!\n")
