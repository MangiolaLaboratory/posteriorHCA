# Frozen Fig. 3 figure-development snapshot (slab-envelope version)

**Frozen on:** 2026-08-19

This is the rendered HTML and source of `TM4SF1_Fig3_figure_development`
immediately before the draw-based revision
(`TM4SF1_Fig3_figure_development_v2_hca_draws.qmd`).

## Known presentation issue (corrected in v2)

- Row 2 Normal “density” slabs were quantile-matched visual envelopes through
  the cached posterior median and 95% interval, not kernel densities of the
  tissue-level posteriorHCA draws used to compute overlap probability.
- HV tumour curves used actual Bayesian-bootstrap study draws.
- That mixed visual grammar can make Normal–HV overlap look stronger than the
  reported \(P(\mu_{\mathrm{Normal}} \ge T_{\mathrm{HV}})\).

## Do not modify

Leave this archive unchanged. Development continues in
`TM4SF1_Fig3_figure_development_v2_hca_draws.qmd`.
The live copies
`TM4SF1_Fig3_figure_development.qmd` and
`TM4SF1_Fig3_figure_development.html` are also left in place and must not be
overwritten by the v2 render.
