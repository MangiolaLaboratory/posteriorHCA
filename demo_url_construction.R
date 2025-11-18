# Quick Demo: Object Storage URL Construction
# ==========================================
# 
# This script demonstrates the URL construction pattern for Nectar object storage

# Load the package (assuming it's installed)
# library(posteriorHCA)

# Your example URL breakdown:
# https://object-store.rc.nectar.org.au/v1/AUTH_b0a86a29c8b74630aac35f471cfe1396/V1/cd4.naive/ENSG00000000419

# Components:
base_url <- "https://object-store.rc.nectar.org.au/v1"
auth_token <- "AUTH_b0a86a29c8b74630aac35f471cfe1396"
container <- "V1"
prefix <- "cd4.naive"
gene_id <- "ENSG00000000419"

# Construct the URL manually:
object_path <- paste(prefix, gene_id, sep = "/")
full_url <- paste(base_url, auth_token, container, object_path, sep = "/")

cat("URL Construction Demo:\n")
cat("=====================\n")
cat("Base URL:", base_url, "\n")
cat("Auth Token:", auth_token, "\n")
cat("Container:", container, "\n")
cat("Prefix:", prefix, "\n")
cat("Gene ID:", gene_id, "\n")
cat("Object Path:", object_path, "\n")
cat("Full URL:", full_url, "\n\n")

# Compare with your original example:
original_url <- "https://object-store.rc.nectar.org.au/v1/AUTH_b0a86a29c8b74630aac35f471cfe1396/V1/cd4.naive/ENSG00000000419"
cat("Original URL:", original_url, "\n")
cat("Constructed URL:", full_url, "\n")
cat("URLs match:", original_url == full_url, "\n\n")

# Show how to use the simplified function (when package is loaded):
# local_path <- get_file_ready(container = "V1", prefix = "cd4.naive", filename = "ENSG00000000419")
# cat("Downloaded file to:", local_path, "\n")

cat("Pattern Summary:\n")
cat("===============\n")
cat("URL Pattern: {base_url}/{auth_token}/{container}/{prefix}/{gene_id}\n")
cat("Example: https://object-store.rc.nectar.org.au/v1/AUTH_b0a86a29c8b74630aac35f471cfe1396/V1/cd4.naive/ENSG00000000419\n")
