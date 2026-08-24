# Formula list for formula_bench.
# Add a row (enabled = TRUE) to extend; old param_grid rows stay cached.

bench_formulas <- function() {
  tibble::tribble(
    ~formula_id,  ~family,      ~k,   ~counts, ~shape, ~label, ~enabled,

    "RE-simple",  "re_string",  NA_integer_,
    paste(
      "counts ~ 1 + offset(offset) + age_decade * sex + ethnicity_groups",
      "+ assay_groups_altered + (1 | dataset_id_altered)",
      "+ (1 + age_decade * sex + ethnicity_groups | tissue_groups)"
    ),
    "shape ~ 1 + assay_groups_altered + (1 | tissue_groups)",
    "Correlated tissue RE + simple shape (healthy-adapted)",
    TRUE,

    "RE-split",   "re_string",  NA_integer_,
    paste(
      "counts ~ 1 + offset(offset) + age_decade * sex + ethnicity_groups",
      "+ assay_groups_altered + (1 | dataset_id_altered)",
      "+ (1 + age_decade * sex | tissue_groups)",
      "+ (0 + ethnicity_groups | tissue_groups)"
    ),
    "shape ~ 1 + assay_groups_altered + (1 | tissue_groups)",
    "Split tissue RE blocks + simple shape (healthy-adapted)",
    TRUE,

    "AGE-k3", "age_smooth", 3L, NA_character_, NA_character_,
    "Continuous-age tp smooth by tissue×sex, k=3", TRUE,

    "AGE-k4", "age_smooth", 4L, NA_character_, NA_character_,
    "Continuous-age tp smooth by tissue×sex, k=4", TRUE,

    "AGE-k5", "age_smooth", 5L, NA_character_, NA_character_,
    "Continuous-age tp smooth by tissue×sex, k=5", TRUE
  )
}

bench_formulas_enabled <- function() {
  dplyr::filter(bench_formulas(), .data$enabled) |>
    dplyr::select(-"enabled")
}
