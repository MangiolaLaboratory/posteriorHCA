# Frozen v2.1 corrected admissibility snapshot

**Frozen on:** 2026-08-18

Preserves v2.1 before the v2.2 presentation/support-visibility correction.

## Known issues (corrected in v2.2)

- Demographic plots (tissue × sex, tissue × age, tissue × ethnicity) filtered
  on `supported` (n_unique_donor_keys >= 10), hiding valid model predictions
  with low/zero direct atlas support.
- Many ordinary both-sex tissues appeared male-only or female-only in figures
  despite both sexes being in the target prediction space.
- No visual encoding of empirical support level in descriptive plots.

## Do not modify

v2.2 development in `TM4SF1_toxicity_case_study_v2.2_visibility.qmd`.
