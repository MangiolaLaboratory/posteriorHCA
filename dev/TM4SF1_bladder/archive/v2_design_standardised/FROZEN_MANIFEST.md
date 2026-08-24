# Frozen v2 design-standardised snapshot

**Frozen on:** 2026-08-18

This directory preserves the v2 design-standardised baseline report before
the v2.1 corrected admissibility revision.

## Known issues (corrected in v2.1)

- Biological sex admissibility was derived from atlas-observed sex per cell type,
  not from anatomical admissibility — caused many both-sex tissues to be labelled
  single-sex.
- Target prediction grid constructed only from observed joint profiles, preventing
  genuine 0.5/0.5 sex weighting for unobserved sex in a tissue.
- No weighted support diagnostics (C_obs, C_k) reported.
- breast classified as hard biological female-only constraint.
- Donor minimum applied inconsistently: described as exclusion but code included
  cells with n > 0.

## Do not modify

v2.1 development occurs in `TM4SF1_toxicity_case_study_v2.1_corrected.qmd`
with outputs under `results/tm4sf1_toxicity_report_v2.1_corrected/`.
