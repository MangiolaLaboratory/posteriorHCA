Cohort expression workflow
================
Chen Zhan
2026-09-04

SAVI case study (*ADRB2* in disease-associated monocytes, PBMC / blood).
The former combined vignette is split into two reports that share the
same marked steps:

| Step | Core | Wrappers |
|---:|----|----|
| 0 | `merge_with_reference_sample()` + `calculate_tmm_offset()` | `scale_to_hca_reference()` |
| 1 | `design_from_formula()` + `estimate_logmu_ql()` | `estimate_cohort_logmu()` |
| 2 | `load_expr_fit()` | `load_expr_fit()` |
| 3 | `build_newdata_grid()` + `expr_draws()` | `expr_predict()` |
| 4 | `summarize_posterior_draws()` | (inside Welch wrapper) |
| 5 | `welch_test_means()` | `welch_t_test_cohort_hca()` |

- Core (matrix API):
  `vignette("cohort-expression-core", package = "posteriorHCA")`
- Wrappers (Seurat / formula):
  `vignette("cohort-expression-wrappers", package = "posteriorHCA")`

Reference script: `examples/savi_adrb2_workflow.R`.
