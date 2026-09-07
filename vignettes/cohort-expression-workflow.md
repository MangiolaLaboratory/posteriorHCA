# Cohort expression workflow (overview)

The full worked example is `vignette("cohort-expression-core", package = "posteriorHCA")`
and `examples/savi_adrb2_workflow.R`.

Intended sequence:

1. Scale user libraries to one HCA reference
2. Build a design with `model.matrix()`
3. Estimate user means with `estimate_logmu_ql()`
4. Load the HCA model with `load_expression_fit()` and draw with `expression_draws()`
5. Summarise with `summarize_posterior_draws()` and test with `welch_test_means()`
