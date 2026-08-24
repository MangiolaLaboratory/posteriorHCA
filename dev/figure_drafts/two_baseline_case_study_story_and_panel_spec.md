# posteriorHCA case-study figures: source-audited story and analysis specification

## Editorial status

This specification replaces the earlier manual-matching and calibration
blueprints. It is based on the source publications, their supplements, the
local reports and the processed data in this repository. Planned panels contain
no invented results: an empty axis or an `ANALYSIS NEEDED` box means that the
specified analysis must be completed before final artwork.

Both main-figure canvases are fixed at the Nature double-column production
size, **183 mm wide × 170 mm high**, with embedded editable fonts at 5–7 pt.
This layout constraint does not remove the need to report full inferential
detail in the legend, Methods and Extended Data.

The two reusable posteriorHCA query modes remain:

1. **Comprehensive cohort query.** Metadata are averaged over one frozen HCA
   standardisation table that is not reweighted to the target cohort. For
   SAVI, the primary cohort comparator is an unseen-study healthy latent mean:
   it includes one shared new-dataset draw but no individual count noise. For
   HCC cross-platform structural comparisons, the dataset effect is fixed at
   its population mean.
2. **Metadata-conditioned query.** The common operation is to fix available
   target metadata and integrate fields that remain unknown. In SAVI this is a
   predictive query: all libraries share one new-study effect per joint draw
   and observation noise is included to generate a new healthy count. In HCC
   it is a structural query: age and sex are overwritten on the same
   healthy-liver standardisation rows, the dataset effect remains zero, and the
   output is an HCA target-gene posterior percentile—not a GeoMx count.

The comprehensive unseen-study mean, population structural mean and
metadata-conditioned predictive count are saved and labelled separately. They
must not be overlaid as if they were the same probability distribution.

Smooth age is one possible covariate representation inside either query. It is
not a third baseline and is not the story of either case study. Until the
production smooth-age models pass their convergence gate, age-decade is the
provisional production specification and the smooth-age comparison belongs in
Extended Data.

---

## Figure 2 — SAVI: external population-reference audit of a published signature

### Source facts that constrain the claim

- The source study contains five SAVI patients and seven in-study healthy
  controls. Controls were selected for absence of a `STING1` mutation and for
  matching the patients' age range; the paper did not manually select an HCA
  comparator.
- There are five untreated patient libraries and five treated libraries from
  four patients. P1 has two treated time points; P5 has no treated library.
- The 21 genes form a **STING-activation signature**. The paper did not claim that
  every constituent is a standalone biomarker satisfying
  `untreated > treated > control`.
- The source explicitly proposed CCL3, CCL4 and IL-6 as individual blood
  biomarker candidates. TSC22D2, IL1B and ADRB2 are signature constituents.
- The signature was derived from disease-associated monocyte cluster 17, whereas
  the current posteriorHCA fit is the broader `V1_monocytic` model. Current
  gene-density, percentile and class results are therefore pilot results until the
  target and reference populations are aligned.

### Defensible main claim

> posteriorHCA adds an external population-scale healthy reference to a
> small, internally controlled rare-disease study. After matching the target
> and atlas cell populations and treating biological subjects as the
> inferential units, the comprehensive query tests whether the published
> STING-activation signature and its covered constituents remain displaced beyond
> broader healthy variation. The metadata-conditioned query reveals
> target-domain shift and subject-sensitive evidence that cohort averages can
> conceal.

This case does **not** claim that the source paper lacked controls, that it used
a manual HCA pipeline, or that failure of one constituent falsifies the
published signature.

### Figure grammar

- comprehensive HCA reference: charcoal grey;
- metadata-conditioned HCA reference: teal;
- study controls: green/open circle;
- untreated SAVI: vermilion/filled circle;
- treated SAVI: purple/square;
- unresolved or pilot-only evidence: neutral grey;
- a visible result glyph has one panel-specific meaning stated beside that
  panel; cells are never plotted as biological replicates;
- line: 95% interval, unless explicitly labelled a sensitivity range;
- no cell-level point may be presented as an independent replicate.

### a — Observed cohort design and HCA training support

**Purpose.** Replace the unsupported manual-matching multiverse with factual
evidence that defines what posteriorHCA can and cannot condition on.

**Left plot: source cohort.**

- x-axis: age at collection in years;
- y-axis: source sequencing batch (`STING1`, `STING2`, `STING3`), with the
  original 10x 3' v2/v3 distinction stated;
- point: one library;
- shape: control, untreated or treated;
- colour and shape: analysis group/state; filled versus open symbol: female
  versus male, avoiding reuse of the cohort colours for sex;
- connector: repeated libraries from the same patient;
- explicit family key or brackets: P2 is P1's uncle, C8 the paternal aunt and
  C9 the mother;
- C1 and C2: a separate non-numeric `Adult; exact age unavailable` category,
  not an imputed point or a finite invented interval;
- P5: use the 10.5-year age reported in source Fig. 1a, not the local
  eight-year imputation.

**Right plot: frozen `V1_monocytic` healthy-blood support.**

- same x-axis: donor age;
- stacked histogram or rug: one unique training donor, split by sex;
- side text: 569 healthy blood donors, 973 samples, 13 datasets and four assay
  groups in the currently inspected model frame;
- assay support: 46 healthy-blood samples in the `10x Genomics 3` harmonised
  assay group;
- overlay: SAVI collection ages or intervals;
- support flag per target query: observed, marginalised or outside support.

These counts were read from the exact current `V1_monocytic` training frame and
must be regenerated automatically from the frozen production model before
submission. Audit the HCA release for overlap with GSE226598 before calling the
reference independent.

**No-refit analysis.**

1. Filter the exact model frame to healthy blood monocytic samples.
2. Deduplicate by donor ID.
3. Count donors, samples, datasets, assay groups and sex/age support.
4. Treat coarse ages as distributions. Do not replace `Adult` with 35.
5. Marginalise unknown ethnicity and unsupported fixed-effect combinations
   according to the locked model query.

### b — Patient-level recovery of published composition directions

**Plot.** Seven-row forest on one common effect axis.

- y-axis, in source order:
  naive CD4 (increased), effector/memory CD4 (decreased), naive CD8
  (increased), effector/memory CD8 (decreased), NK (decreased), MAIT
  (decreased), and gamma-delta T (decreased);
- open circle: untreated SAVI minus study-control composition effect;
- grey diamond: untreated SAVI minus comprehensive-HCA composition effect;
- effect scale: primary percentage-point difference in posterior expected
  proportion; additive log-ratio is an Extended Data sensitivity;
- open circle and line: untreated-minus-study-control posterior median and 95%
  interval;
- grey diamond and line: untreated-minus-comprehensive-HCA posterior median
  and 95% interval;
- far-left arrow: source-paper direction, shown as context rather than a new
  result.

Treatment observations may be shown as linked context, but treatment
normalisation is not the endpoint because the paper reported no significant
pre/post composition difference.

**Analysis.**

1. Recompute per-subject cell counts from the public Seurat object.
2. Freeze the mapping from the paper's labels to posteriorHCA composition
   labels, including an `other/unmapped` denominator category.
3. Fit `sccomp` 2.2.0 with
   `noise_model = "multi_beta_binomial"`,
   `formula_composition = ~ group + sequencing_batch`,
   `formula_variability = ~1`, HMC inference, raw counts and eight mutually
   exclusive categories: the seven displayed compartments plus `other`.
4. Regenerate the frozen HCA model's full 26-L3 joint draws; aggregate the
   seven compartments and define `other` as all remaining L3 mass within the
   same `.draw`. Every eight-part vector must sum to one. Do not renormalize
   separately fitted compartments or reuse the cached table trimmed to 17 L3
   types.
5. Give P1 one patient-level weight despite two treated libraries.
6. Standardise HCA using the frozen comprehensive HCA donor weights, not the
   SAVI age/sex distribution and not one weight per target library.
7. Sensitivities: inclusion/exclusion of `other`, joint C8/C9 removal,
   related-family removal and additive log-ratio parameterisation.

The current report provides only pooled pilot means. They remain absent from
the blueprint so that placeholders cannot be mistaken for results.

### c — Signature-first external-reference evidence

**Primary top strip: the posteriorHCA-covered signature effect.**

- four columns aligned above the constituent matrix;
- x-axis in every column: mean of the 16 gene-wise log2 effects;
- diamond: covered-subset median;
- line: 95% subject-coupled product-posterior stability envelope; this is
  explicitly not an ordinary credible interval because the genes were fitted
  separately;
- zero: geometric-mean fold change of one;
- printed coverage: 21 source genes, 19 currently mapped to Ensembl and 16
  with current posteriorHCA models;
- label `16-gene covered subset` until the proxy check below passes.

Do not invent a standardized score or overlay the full source score with an
HCA density. Recover the source's exact normalization and scoring procedure,
then calculate the 21- and 16-gene versions on the same SAVI scale. Show their
paired scatter, Spearman correlation, group directions and leave-one-gene-out
range in Extended Data. If the scoring procedure is not recoverable, or the
16-gene subset changes a main source direction, do not call it a proxy for the
21-gene signature.

For each of the four columns, the covered-subset effect is the equal-weight
mean of the 16 gene-wise log2 effects. It is therefore the log2 of their
geometric-mean fold change. The stability envelope uses separate Dirichlet
weights for the five untreated patients and seven controls, one shared
four-patient vector for the paired treatment contrast, one shared HCA-donor
vector, and a conservative common standardized new-study effect across genes.
The exact algorithm and factorized sensitivities are defined in the detailed
analysis plan.

**Constituent matrix below.**

- y-axis: 16 covered constituents;
- aligned columns:
  1. untreated SAVI versus study controls;
  2. untreated SAVI versus comprehensive HCA;
  3. untreated versus treated, using a paired-patient estimand;
  4. study controls versus HCA, as target-domain shift;
- shared x-axis: log2 fold change or another locked common effect scale;
- point: posterior median;
- line: 95% interval;
- zero line: no effect;
- a practical-equivalence band is permitted only after a biological margin is
  justified;
- no significance stars or compressed binary class; retain all four estimates.

**Analysis.**

1. Resolve cell-population alignment before fitting:
   - reconstruct the source signature evidence in the population used by the
     paper (monocytes plus DCs) and keep it separate from the HCA-matched
     external audit;
   - primary reference-aligned target: all monocyte-lineage clusters 5, 12 and
     17 if this matches `V1_monocytic`;
   - CD14-like clusters 5 and 17 as a declared sensitivity;
   - cluster 17 alone only as a source-state analysis, not a matched-HCA
     analysis;
   - monocytes plus DCs only if the estimand is explicitly changed.
2. Fit a donor-aware target count model including group and sequencing batch.
3. Combine P1's two treated libraries hierarchically or average them on the
   estimand scale so P1 has one paired treatment weight.
4. Combine target-model and frozen HCA draws; do not treat HCA draws as
   independent observations in a Welch test.
5. Recover and execute the exact source scoring procedure before any
   21-versus-16 proxy claim.
6. If `untreated > treated > control` is used, label it a new stringent
   standalone-constituent criterion, not the source signature definition.
7. Use the two prespecified covered-subset effects as one joint corroboration
   gate; constituent intervals are descriptive. Any formal constituent screen
   requires a separately frozen multiplicity procedure.

**Separate source-nominated biomarker tier.** If the manuscript retains a
claim about falsifying or challenging *published biomarkers*, audit CCL3, CCL4
and IL6 as a separate prespecified three-gene analysis because these—not
TSC22D2, IL1B or ADRB2—were the individual blood biomarker candidates named by
the source. First verify model availability and cell-population compatibility.
If those three cannot be analysed, the manuscript and Fig. 2 must use
`signature constituent audit`, not `published biomarker falsification`.

### d — Comprehensive and metadata-conditioned density archetypes

**Three non-duplicating columns selected in order after panels c and f are
complete:**

1. `external + internal direction`: largest lower 95% bound of
   untreated-minus-HCA among genes classified direction-positive for both
   untreated-minus-HCA and untreated-minus-study-control;
2. `study-control/HCA shift`: largest absolute control-minus-HCA median among
   the remaining genes whose 95% interval is wholly above or below zero;
3. `largest subject influence`: among remaining genes, first select those whose
   direction category changes under deletion and maximize \(R_g\); otherwise
   maximize \(R_g\).

Remove each selected gene before applying the next rule and break ties
alphabetically. A rule with no eligible gene leaves a labelled blank slot.

**Upper row: comprehensive cohort query.**

- x-axis: `offset-zero log2 marginal expected count`;
- y-axis: posterior density;
- grey filled curve: comprehensive unseen-study healthy latent mean, including
  one shared new-dataset draw but excluding individual count noise;
- green, vermilion and purple curves: study-control, untreated and descriptive
  treated standardised latent-mean posteriors; matched vertical ticks are their
  medians;
- this is a latent-mean distribution for an unseen study, not a distribution
  of observed healthy counts and not merely uncertainty in the structural HCA
  population mean.

**Lower row: metadata-conditioned query.**

- x-axis:
  \([\log(Y+0.5)-\mathrm{offset}]/\log 2\), labelled
  `offset-centred log2 count`;
- y-axis: posterior predictive density;
- five vertically separated teal mini-ridges: one metadata-conditioned
  new-healthy-sample prediction for each untreated SAVI patient;
- five matched vertical target marks: the corresponding observed untreated
  libraries;
- direct labels or a compact key: patient ID, age specification, sex and
  healthy-prediction percentile;
- no single patient is selected or highlighted, preventing an exemplar from
  being chosen by expression outcome.

The current raster crops are not reused in the production-size blueprint:
they distort aspect ratio, omit final axes and come from the cell-misaligned
pilot. Panel d therefore reserves native-vector density axes at final size.
No provisional gene name is printed. Final examples are selected by the rules
above, not visual extremeness.

### e — Metadata-conditioned healthy-prediction percentiles across all libraries

**Plot.**

- x-axis: 17 libraries grouped as seven controls, five untreated and five
  treated;
- y-axis: 16 covered constituents in panel-c order;
- tile: metadata-conditioned healthy-prediction percentile,
  blue–white–vermilion around 0.5;
- right strip: for each gene separately, a horizontal line from the smallest
  to the largest percentile among the seven controls and a diamond at their
  median.

The main panel does not overlay nominal or multiplicity-adjusted significance
dots. A percentile of 0.95 means that the observed count is higher than
approximately 95% of simulated healthy counts for that same gene and
library-metadata query. It is not a 95% probability of disease. These
percentiles remain descriptive until transport is adequate and calibration
has been demonstrated in the core method validation.

The control strip is a **target-domain diagnostic**, not an empirical
calibration curve. The legacy percentile export is not reusable because it imputed
C1/C2 as 35 years, used P5 as eight rather than 10.5 years, and used the
cell-misaligned cluster-17 target. The blueprint therefore leaves the percentile
matrix blank until the corrected analysis is run.

**Analysis.**

1. Correct all source ages; integrate C1/C2 over a declared adult-age
   distribution.
2. Use the reference-aligned pseudobulk target.
3. Draw one shared unseen target-dataset effect per posterior draw when making
   joint cohort claims or simulations.
4. Preserve the repeated-patient structure.
5. Generate three batches of 4,000 predictions per gene × library
   (\(M=12{,}000\)); require batch percentile range no greater than 0.02 and
   double all batches if needed.
6. If target controls are shifted for a gene, suppress disease-specific
   personalised calls or label them untransportable.

### f — Biological-subject and family influence

**Main-figure plot.**

- first y-axis row: the 16-gene covered subset; remaining rows: at most five
  selected constituents;
- left facet x-axis: absolute log2 untreated-SAVI-versus-comprehensive-HCA
  effect;
- right facet x-axis: absolute log2 untreated-SAVI-versus-study-control
  effect;
- vertical zero line: no biological difference;
- black diamond and line: full-data median and interval; for the 16-gene
  aggregate row the line is the 95% subject-coupled stability envelope, while
  constituent rows use their gene-specific 95% credible intervals;
- open points: subject-deleted posterior medians;
- horizontal range: minimum to maximum deleted-data median;
- purple: deletion changes `direction positive`, `direction negative` or
  `direction unresolved`;
- right label: influential subject or family deletion;
- all 16 covered constituents appear in Extended Data.

There is no arbitrary 0.8 threshold.

**Analysis.**

- Run all 12 individual deletions in both facets. A control can change the
  fitted untreated mean through the shared intercept, batch adjustment,
  dispersion and partial pooling even though it is not an arithmetic term in
  untreated-minus-HCA.
- Deleting a patient removes every associated library; deleting P1 removes all
  three P1 libraries.
- Evaluate every deletion fit on the full-data standardisation/batch grid.
- For gene \(g\) and estimand \(q\), define
  \(R_{gq}=\max_d\operatorname{median}(\Delta_{gq,-d})-
  \min_d\operatorname{median}(\Delta_{gq,-d})\), and
  \(R_g=\max(R_{g,U-HCA},R_{g,U-C})\). Selection and tie rules are frozen in
  the detailed plan.
- Add joint removal of C8/C9 and of the connected P1/P2/C8/C9 family in
  Extended Data; promote a family deletion to the main panel only if it alone
  changes a main conclusion.
- Keep every posteriorHCA model fixed; these are target-only fits.

---

## Figure 3 — HCC GeoMx: source-feature-excluded retrieval of HCA reference structure

### Source and capsule facts that constrain the claim

- The paper reports GeoMx WTA profiling of 18,677 genes in 64 source-indexed
  NK-enriched segments from eight HCC tissues.
- The CodeOcean matrix contains 93 geometric-segment columns:
  47 high-coverage columns with 18,393 nonmissing genes and a 46-column
  1,796-gene branch.
- The five TIMES genes occur only in the 18,393-gene branch.
- The union of the source Code Ocean objects `NKatBorder_indx`,
  `NKatStroma_indx` and `NKatTumor_indx` selects 64 segments:
  30 high-coverage and 34 in the 1,796-gene branch. `NKat` is a source object
  name, not a posteriorHCA score.
- Seven A38 segments are called invasive front by the author index but carry
  `stroma-*` labels. The current local extraction follows the labels and the
  source DE follows the author index. This must be shown and resolved.
- There is no healthy GeoMx bridge. Absolute HCA-to-GeoMx count compatibility
  and calibrated GeoMx predictive percentiles are not identifiable.

### Defensible main claim

> A recurrence-label-blind, source-feature-excluded broad-gene test asks
> whether posteriorHCA retrieves liver-NK structure from source-labelled
> NK-enriched GeoMx segments without using spatial coordinates. Segment-label
> reallocations, alternative reference cells and spatial-processing
> sensitivities determine whether retrieval survives the assay shift.
> Comprehensive and age/sex-conditioned structural HCA target-gene posterior
> percentiles provide descriptive reference context, not calibrated GeoMx
> prediction.

The case does not claim general performance for every spatial programme,
GeoMx count generation, healthy compatibility, calibrated tail probabilities,
population-level recurrence prediction or spatial-coordinate modelling.

### a — Auditable data flow and compartment-discrepancy gate

**Plot.** Branched inclusion flow plus the coverage-by-selection cross-tab.

1. Start: 93 CodeOcean geometric-segment columns.
2. Two crossing classifications, shown as a cross-tab rather than a
   sequential filter:
   - author selection: 64 NK-selected segments and 29 context controls;
   - coverage: 47 high-coverage columns with 18,393 genes and 46 columns in
     the 1,796-gene branch;
   - cross-tab: high coverage contains 30 NK-selected plus 17 controls; the
     1,796-gene branch contains 34 NK-selected plus 12 controls;
   - none of the five TIMES genes occurs in the 1,796-gene branch.
3. IF/TC target branch:
   - label-derived assignment: 20 high-coverage segments and five
     patient × compartment profiles;
   - author-index assignment: 27 high-coverage segments and six profiles,
     adding A38 IF;
   - discrepancy: seven A38 segments.
4. Context-specificity branch:
   - all 47 high-coverage segments from four patients.

The matrix prints exact segment counts. Patient is the biological unit; segment
is a spatial sampling unit. Adjacent stroma is a spatial context, not healthy
liver.

**Required resolution.**

- Use author-index assignment for source-reconstruction panel d.
- Use label-derived assignment as a prespecified sensitivity, or obtain the
  raw OMIX005738-03 labels/author clarification.
- Do not call the current 20-to-5 extraction source-faithful while the conflict
  remains.

### b — Recurrence-label-blind, source-feature-excluded reference retrieval

**Fixed, feasible analysis.**

- all 47 high-coverage GeoMx segments from four patients define gene
  completeness and the candidate background universe;
- 27 non-overlap segments enter the four primary contrasts, eight overlap
  segments enter only the overlap sensitivity, and 12 other-context segments
  do not enter the primary contrasts;
- candidate HCA references already confirmed to contain healthy liver support:
  `V1_nk`, `V1_cd8.tem`, `V1_mait`, `V1_tgd`, `V1_monocytic`,
  `V1_macrophage` and `V1_epithelial`;
- exact current intersection: 3,293 genes nonmissing in all 47 segments and
  present in all seven model stores;
- reconstruct the 241 source spatial genes as the union of the `gene` columns
  in `fd_df` (113), `fd_uf` (3), `fu_df` (19) and `fu_uf` (106) from
  `differential_gene_expression.RData`; require 241 unique symbols and checksum
  the manifest;
- 47 of those 241 intersect the 3,293-gene common universe; removing them
  leaves exactly 3,246 recurrence-label-blind, source-feature-excluded genes
  and removes all five TIMES genes;
- add a declared sensitivity that removes the transcript counterparts of the
  segment-selection proteins if they are not already in the excluded set; report
  the resulting universe size rather than continuing to call it 3,246.
- record the identifier-mapping package/version and the one-to-many mapping
  rule because the exact intersection depends on them.
- after the production formula is fixed, require converged usable fits in
  every HCA store needed by a panel and re-freeze the background. If the final
  size is not 3,246, update every denominator and label rather than retaining
  3,246 as a historical name.

**Primary contrast calculation.**

For each of the 3,246 background genes and each matched block, compute

`median log2(Q3) in CD57-rich segments − median log2(Q3) in comparator
segments`.

The exact non-overlap primary blocks are:

- A24 TC: two CD57-rich versus three cold-tumour segments;
- A29 TC: four versus four cold-tumour segments;
- A6 IF: three CD57-rich versus two CD3-rich segments;
- A38 AS: four CD57-rich versus five CD3-rich segments.

Segments labelled both CD57-rich and CD3-rich are excluded in the primary
calculation and included with the CD57-rich set in a named sensitivity.

For each candidate HCA reference and posterior draw, calculate

`healthy-liver candidate-cell latent expression − equal-weight mean latent
expression across all seven candidate cells`,

using an identical metadata-standardisation table and dataset effects fixed at
their population mean. Compare the GeoMx and HCA 3,246-gene contrast vectors
with Spearman correlation. This is a contrast-to-contrast test; direct
Spearman similarity of absolute within-profile cross-gene expression
orderings is platform-sensitive and belongs only in Extended Data.

The common HCA grid uses one canonical
`source-study accession::source donor ID` key, retains only complete joint
metadata tuples supported by all seven stores, restricts age to the
within-tuple common support, gives each retained donor equal weight and uses
the identical rows for every cell type. Fewer than 20 retained donors, loss of
either sex, no common assay level or exclusion of more than 20% of pooled donor
mass stops the analysis. A separate go/no-go audit must confirm identical
pseudobulk response definitions, TMMwsp reference/feature checksums, offsets,
marginal expected-expression scale, identifiers and liver/Normal ontology
across stores.

**Plot.**

- left y-axis: seven candidate healthy-liver HCA references;
- left x-axis: Spearman correlation, including zero and negative values;
- four shaped points per row: posterior medians for A24 TC, A29 TC, A6 IF and
  A38 AS;
- line through each point: 2.5th–97.5th model-plus-gene-resampling stability
  range, not a patient-level confidence interval;
- highlighted row: liver NK, fixed before viewing the result;
- right y-axis: the four matched blocks;
- right x-axis: NK retrieval margin
  `rho(liver NK) − max rho(other six references)`;
- vertical zero: NK ties the strongest alternative;
- negative-margin label: identity of the winning alternative.

No overall diamond is shown because the four blocks mix TC, IF and AS and use
two different comparator definitions; they do not share a pooled estimand.
Calculate the NK margin inside each paired draw/resample with the same
resampled genes for all seven references.

Run 2,000 posterior-only, 2,000 expression-stratified gene-set-only and 2,000
combined replicates, with separately reported ranges. Repeat with new seeds and
require range endpoints within 0.02. Add 50 MSigDB Hallmark
leave-one-set-out sensitivities and enumerate every within-block rich/context
label reallocation preserving the observed group sizes (10, 70, 10 and 126
assignments). Strong seven-reference retrieval requires, for all four blocks,
the lower 2.5% limits of both liver-NK similarity and NK margin above zero, an
observed-label reallocation percentile of at least 0.90, positive sensitivity
margins, and NK retrieval against every eligible liver store—not only the
seven displayed competitors.

No recurrence label or TIMES biomarker enters reference selection.
`A38 AS` is the literal-label context block used only for this specificity
test; the profile-level primary analyses in c–f use author-index `A38 IF`.

### c — Segment aggregation, mapping and compartment robustness

**Plot.**

- primary grid: 24 endpoints, four supported TIMES genes by six author-index
  patient × compartment profiles;
- cell x-axis: target-gene percentile from zero to one;
- black diamond: primary analysis;
- horizontal range: all declared perturbations;
- small symbols: individual perturbations;

The label-derived five-profile, 20-endpoint result is one named sensitivity,
not the primary grid. Because A38 IF disappears under the label-derived
assignment, its four cells carry a categorical `profile absent` flag for that
sensitivity; they must not show a continuous range pretending that the missing
profile produced a numeric estimate.

**Primary target processing.**

- author Q3-normalised continuous values;
- for each target, insert it beside the 3,246 recurrence-label-blind,
  source-feature-excluded background genes and use the tie-adjusted percentile
  with denominator 3,247;
- calculate the percentile separately per segment, then take the median per
  patient × compartment;
- include all 27 author-index IF/TC segments, including CD57/CD3-overlap
  segments. Panel b uses a different primary overlap rule because it requires
  disjoint rich and comparator groups.

A value of 0.90 means that the target is above approximately 90% of the
background genes in that GeoMx segment. It does not mean that the patient is
above 90% of healthy people. The target is inserted as a query and is never
called one of the 3,246 background genes.

**Exact perturbations.**

- mean versus median aggregation;
- every leave-one-segment-out result;
- enumerate every two-segment subset where possible;
- include/exclude CD57/CD3-overlap segments;
- locked-gene-universe variants;
- author-index versus label-derived compartment assignment;
- calculate a target-gene percentile per segment and then aggregate versus
  aggregate Q3 across segments and then calculate one percentile.

Rounding processed values and applying TMMwsp is a legacy sensitivity only, not
the primary cross-platform analysis. There is no arbitrary pass threshold.

### d — Source-pattern reconstruction

**Plot.**

- y-axis: the five TIMES genes;
- x-axis:
  `log2(non-REC group mean / REC group mean)` after first taking patient-level
  segment means on the continuous-Q3 scale;
- circle: IF contrast;
- square: TC contrast;
- zero line: no recurrence-group difference;
- annotation: `technical extraction check: 10/10 source directions recovered`.

Use the author-index compartment assignment. The displayed pilot uses the
audited patient-mean continuous-Q3 reconstruction: IF compares two non-REC
with two REC patients; TC compares one non-REC with one REC patient. These
effective biological counts must be printed. No inferential P-value or AUC is
defensible at this size.

This panel demonstrates that the public-data extraction retains the source
pattern needed before the posteriorHCA transfer analyses in panels b, c, e and
f. It is not itself posteriorHCA validation.

### e — Comprehensive HCA target-gene percentile context

**Plot.**

- y-axis: SPON2, VIM, ZFP36 and ZFP36L2;
- x-axis:
  `HCA target-gene posterior percentile against 3,246 background genes`;
- grey density: posterior uncertainty in the comprehensive healthy-liver-NK
  target-gene percentile;
- six coloured points per gene row: the author-index patient × compartment
  GeoMx target-gene percentiles;
- horizontal bar: panel-c processing and segment-assignment sensitivity range;
- point colour: patient; point shape: IF or TC.

The HCA density is uncertainty in a cross-gene reference percentile, not a healthy
GeoMx predictive distribution. GeoMx bars are perturbation ranges, not
credible intervals.

**Analysis.**

1. Generate HCA cross-gene vectors jointly across the locked background plus
   one inserted query target. For each assembly replicate, sample one retained
   draw independently from every gene fit; never align equal MCMC row numbers
   across separately fitted genes.
2. Use identical gene identifiers and the same frozen background universe in HCA and
   GeoMx.
3. Exclude HLA-DRB1 because no current `V1_nk` model branch is available.
4. Propagate HCA posterior variation and all panel-c perturbations.
5. Use denominator 3,247 and never describe a TIMES target as a member of the
   3,246-gene background.
6. Generate three 2,000-vector batches (6,000 total) per target with frozen
   seeds. Require batch medians to span no more than 0.01 and 95% endpoints no
   more than 0.02; otherwise double all batches.

### f — Age/sex-conditioned structural HCA reference sensitivity

**Primary plot.**

- x-axis:
  `HCA target-gene posterior percentile against 3,246 background genes`,
  shared by the HCA references and the fixed GeoMx observation;
- y-axis: the four supported TIMES genes;
- grey point/interval: comprehensive reference;
- teal point/interval: age/sex-conditioned structural HCA reference;
- arrow: shift caused only by reference conditioning;
- black vertical mark: fixed GeoMx target-gene percentile, so the descriptive GeoMx–HCA gap is
  visible before and after conditioning;
- facets: A24 (41-year-old female), A38 (51-year-old female), and A6/A29
  (57-year-old male); A6 and A29 share one HCA reference but retain separate
  fixed GeoMx marks.

Start from the exact comprehensive standardisation rows and weights and
overwrite age/sex only; retain each row's ethnicity and HCA assay. Set the
dataset effect to zero. Reuse the same gene-specific sampled draw indices for
comprehensive and conditioned vectors so each arrow is a paired age/sex shift.
With the provisional age-decade model, the actual queries are decade 5 for
41F and decade 6 for both 51F and 57M; exact ages are used only if the smooth
model passes production QC.

A24 shows separate IF and TC GeoMx marks; A38 shows author-index IF with a
dagger and states `literal-label IF profile absent`; A6/A29 retains separate A6
IF, A29 IF and A29 TC marks. Every mark carries the panel-c processing range.

**Density inset.**

- x-axis: one prespecified target-gene percentile;
- y-axis: HCA reference density;
- grey curve: comprehensive reference;
- light teal: A24-conditioned reference;
- medium teal: A38-conditioned reference;
- dark teal: A6/A29-conditioned reference;
- black marks: fixed GeoMx SPON2 percentiles;
- SPON2 is fixed before HCA results because
  `dev/Jie_HCC/data/times_biomarkers.csv` labels it the primary/most predictive
  source feature.

This demonstrates sensitivity of the HCA reference to available metadata. It
does not demonstrate calibrated personalised prediction for GeoMx.

**Gate before final rendering.**

- freeze the age formula independently of HCC;
- require converged production fits across the locked background universe;
- the current 12 four-gene smooth-age pilots fail their stated divergence gate,
  so smooth-age values cannot populate this panel yet;
- check age support and extrapolation;
- carry all panel-c mapping sensitivities into the displayed GeoMx
  target-gene percentiles.

---

## Ordered analysis programme

### Gate 1 — Provenance and target definitions

1. Repair the official SAVI Table S6 file and document the final 21-gene list.
2. Audit overlap between posteriorHCA training data and both target studies.
3. Freeze SAVI target-cell mappings and HCC segment/compartment assignments.
4. Resolve the seven-segment A38 discrepancy or retain both assignments as
   labelled analyses.

### Gate 2 — Freeze posteriorHCA independently

1. Fix the production formula, priors, support rules and standardisation
   weights outside these case-study outcomes.
2. Resolve smooth-age convergence before it replaces age-decade.
3. Define exactly when dataset and donor effects are marginalised,
   conditioned or sampled as new levels.
4. Use a shared new-study draw for joint target predictions.

### Gate 3 — SAVI analyses

1. Correct ages and cell-population alignment.
2. Fit patient-aware composition and expression models.
3. Test the signature first; audit constituents second.
4. If retaining a published-biomarker claim, separately fit CCL3, CCL4 and
   IL6; otherwise remove that claim.
5. Recompute comprehensive unseen-study-mean and metadata-conditioned
   predictive densities as native vector plots.
6. Run 12 subject deletions and family sensitivities.
7. Move runtime/resource benchmarking to Extended Data.

### Gate 4 — HCC analyses

1. Lock the 3,246-gene recurrence-label-blind, source-feature-excluded common
   universe.
2. Run the contrast-based correct-reference specificity analysis and
   expression-stratified gene resampling.
3. Run the full segment aggregation and compartment-assignment perturbation
   grid.
4. Reconstruct the author-index source contrast.
5. Build comprehensive and metadata-conditioned target-gene percentile
   posteriors.
6. Keep all outcome labels descriptive at the current patient count.

## Claims prohibited in the current manuscript/figures

- an original SAVI manual-HCA comparator or matching multiverse;
- resource superiority over an unreported source pipeline;
- falsification of the published 21-gene signature by one constituent;
- a claim that every signature constituent should obey
  `untreated > treated > control`;
- equivalence without a justified margin;
- Welch tests that use HCA posterior draws as independent sample size;
- HCC healthy compatibility, GeoMx calibration, GeoMx predictive P-values or
  GeoMx count generation;
- HCC recurrence performance from the current reanalysis;
- any schematic numerical result rendered as observed evidence.
