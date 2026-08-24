# Frozen v1 observed-profile snapshot

**Frozen on:** 2026-08-18

This directory is a reproducibility snapshot of the TM4SF1 toxicity case study
**before** the v2 design-standardised baseline redesign.

## Frozen files

| File | Description |
|---|---|
| `TM4SF1_toxicity_case_study_v1_observed_profile.qmd` | Frozen Quarto source (v1 estimand) |
| `TM4SF1_toxicity_case_study_v1_observed_profile.html` | Rendered HTML at freeze time (if present) |
| `results_v1_observed_profile/` | Complete v1 machine-readable outputs |

## Do not modify

The working original `TM4SF1_toxicity_case_study.qmd` at the case-study root was
left unchanged. All v2 development occurs in
`TM4SF1_toxicity_case_study_v2_design_standardised.qmd` with outputs under
`results/tm4sf1_toxicity_report_v2_design_standardised/`.

## v1 baseline estimand

Equal weight per observed Normal/10x3 donor–tissue–demographic profile row.
Comprehensive baseline = mean over all profile-row predictions.
