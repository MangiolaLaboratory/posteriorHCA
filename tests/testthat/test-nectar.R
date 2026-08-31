library(testthat)

if (requireNamespace("dittoSeq", quietly = TRUE) &&
    requireNamespace("posteriorHCA", quietly = TRUE)) {
  library(posteriorHCA)
} else {
  suppressPackageStartupMessages({
    library(cli)
    library(httr)
    library(readr)
  })
  utils_file <- file.path(testthat::test_path("..", ".."), "R", "utlis.R")
  if (!file.exists(utils_file)) {
    utils_file <- file.path("R", "utlis.R")
  }
  sys.source(utils_file, envir = environment())
}

fake_catalog <- function() {
  data.frame(
    container = c("V1", "V1", "V2"),
    ct_name = c("cd8.naive", "nk", "cd8.naive"),
    ver_no = c("1", "1", "2"),
    n_genes = c(7920L, 100L, 8000L),
    schema = "expr_v1",
    updated_at = "2026-08-26",
    stringsAsFactors = FALSE
  )
}

nectar_object_url_fn <- function(...) {
  if (exists("nectar_object_url", mode = "function", inherits = TRUE)) {
    nectar_object_url(...)
  } else {
    getFromNamespace("nectar_object_url", "posteriorHCA")(...)
  }
}

test_that("nectar_object_url matches the published layout", {
  expect_equal(
    nectar_object_url_fn("meta", "latest.csv"),
    "https://object-store.rc.nectar.org.au/v1/AUTH_b0a86a29c8b74630aac35f471cfe1396/meta/latest.csv"
  )
  expect_equal(
    nectar_object_url_fn("V1", "cd8.naive", "ENSG00000000419"),
    "https://object-store.rc.nectar.org.au/v1/AUTH_b0a86a29c8b74630aac35f471cfe1396/V1/cd8.naive/ENSG00000000419"
  )
  expect_equal(
    nectar_object_url_fn("V1", "reference_samples", "cd8.naive.rds"),
    "https://object-store.rc.nectar.org.au/v1/AUTH_b0a86a29c8b74630aac35f471cfe1396/V1/reference_samples/cd8.naive.rds"
  )
})

test_that("lookup_cell_type_storage uses make.names and latest version", {
  cat <- fake_catalog()
  loc <- lookup_cell_type_storage("cd8 naive", version = "latest", catalog = cat)
  expect_equal(loc$cell_type, "cd8.naive")
  expect_equal(loc$container, "V2")
  expect_equal(loc$ver_no, "2")

  loc_v1 <- lookup_cell_type_storage("cd8 naive", version = "V1", catalog = cat)
  expect_equal(loc_v1$container, "V1")
  expect_equal(loc_v1$ver_no, "1")

  loc_num <- lookup_cell_type_storage("cd8 naive", version = "1", catalog = cat)
  expect_equal(loc_num$container, "V1")
})

test_that("lookup_cell_type_storage errors for unknown cell types", {
  expect_error(
    lookup_cell_type_storage("not a type", catalog = fake_catalog()),
    "not in the Nectar catalog"
  )
})

skip_if_nectar_unavailable <- function() {
  skip_on_cran()
  url <- nectar_object_url_fn("meta", "latest.csv")
  ok <- tryCatch({
    resp <- httr::GET(url, httr::timeout(20))
    httr::status_code(resp) == 200L
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) {
    skip("Nectar object storage is not reachable")
  }
}

test_that("cd8 naive is in latest catalog and downloads from V1", {
  skip_if_nectar_unavailable()
  skip_if_not_installed("qs2")
  cache <- tempfile("posteriorHCA-cd8-")
  dir.create(cache)

  catalog <- get_nectar_catalog(
    "latest",
    cache_directory = cache,
    use_cache = FALSE
  )
  expect_true("cd8.naive" %in% catalog$ct_name)

  loc <- lookup_cell_type_storage(
    "cd8 naive",
    catalog = catalog
  )
  expect_equal(loc$cell_type, "cd8.naive")
  expect_equal(loc$container, "V1")
  expect_equal(loc$ver_no, "1")

  gene <- "ENSG00000000419"
  fit_res <- get_brms_ready(
    "cd8 naive",
    gene,
    cache_directory = cache,
    use_cache = FALSE
  )
  expect_equal(fit_res$status, "success")
  expect_true(file.exists(fit_res$path))
  expect_gt(file.size(fit_res$path), 0)
  expect_equal(
    fit_res$url,
    nectar_object_url_fn("V1", "cd8.naive", gene)
  )

  obj <- qs2::qs_read(fit_res$path)
  expect_true("brms_fit" %in% names(obj))
  expect_s3_class(obj$brms_fit[[1]], "brmsfit")

  ref_res <- get_reference_sample_ready(
    "cd8 naive",
    cache_directory = cache,
    use_cache = FALSE
  )
  expect_equal(ref_res$status, "success")
  expect_true(file.exists(ref_res$path))
  expect_equal(
    ref_res$url,
    nectar_object_url_fn("V1", "reference_samples", "cd8.naive.rds")
  )

  skip_if_not_installed("SummarizedExperiment")
  ref <- readRDS(ref_res$path)
  expect_s4_class(ref, "SummarizedExperiment")
  expect_equal(ncol(ref), 1L)
  expect_true("counts" %in% SummarizedExperiment::assayNames(ref))
  expect_equal(nrow(ref), 7920L)
})
