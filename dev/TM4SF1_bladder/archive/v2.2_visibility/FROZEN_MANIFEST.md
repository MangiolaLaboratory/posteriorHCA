# Frozen v2.2 support-visibility snapshot

**Frozen on:** 2026-08-18

## Known issue (corrected in v2.2.1)

- Tissue × sex plotting audit joined `tissue_sex_rules_all` by `tissue_groups`
  only, creating duplicate rows when epithelial and endothelial have different
  `atlas_observed_sex` for the same tissue.

## Do not modify

v2.2.1 development in `TM4SF1_toxicity_case_study_v2.2.1_audit_cleanup.qmd`.
