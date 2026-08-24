# Chain / CPU settings for formula_bench.
# Add a row (enabled = TRUE) to extend; old param_grid rows stay cached.
#
# threads_per_chain = max(1, floor(min(availableCores(), cpus) / chains))

bench_resource_settings <- function() {
  tibble::tribble(
    ~setting_id, ~chains, ~cpus, ~label, ~enabled,
    "c2_cpu8",   2L,      8L,    "2 chains × 8 CPUs",  TRUE,
    "c4_cpu8",   4L,      8L,    "4 chains × 8 CPUs",  TRUE,
    "c4_cpu16",  4L,      16L,   "4 chains × 16 CPUs", TRUE
  )
}

bench_resource_settings_enabled <- function() {
  dplyr::filter(bench_resource_settings(), .data$enabled) |>
    dplyr::select(-"enabled")
}

bench_mcmc <- function() {
  list(
    warmup = as.integer(Sys.getenv("BENCH_WARMUP", unset = "400")),
    iter   = as.integer(Sys.getenv("BENCH_ITER", unset = "900")),
    seed   = as.integer(Sys.getenv("BENCH_SEED", unset = "20260729"))
  )
}

bench_resolve_threads <- function(chains, cpus) {
  chains <- as.integer(chains)
  cpus <- as.integer(cpus)
  avail <- suppressWarnings(as.integer(parallelly::availableCores()))
  if (!is.finite(avail) || avail < 1L) avail <- cpus
  avail <- min(avail, cpus)
  max(1L, as.integer(floor(avail / chains)))
}
