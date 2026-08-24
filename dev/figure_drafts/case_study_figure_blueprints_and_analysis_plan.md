# posteriorHCA case-study figure blueprints and analysis plan

## Status and authority

This is the single canonical design document for Figures 2 and 3. It replaces
all earlier slide designs, review notes, story drafts and panel specifications
in `dev/figure_drafts`.

The accompanying PDFs are analysis blueprints, not completed result figures.
A dashed frame or orange `NEEDED` label means that the named analysis must be
completed before submission. Empty axes contain no simulated or invented
effect estimates. A populated result is shown only where its source table has
already been audited.

Both PDFs are rendered at 183 mm × 170 mm, with lower-case panel labels,
embedded fonts and vector-native graphics. They contain no figure-level title;
the manuscript legend supplies the title and full caption.

---

## 1. Case-study claims

### Figure 2: SAVI

The SAVI study already contains seven study controls. posteriorHCA is used as
an additional population-scale reference, not as a replacement for those
controls and not as a reconstruction of an imagined manual HCA-matching
pipeline.

The primary question is:

> Does the published SAVI transcriptional or composition difference remain
> unusual relative to healthy variation learned from many HCA cohorts?

The comprehensive baseline assesses the cohort-level effect. The
metadata-conditioned baseline asks whether every observed library is unusual
relative to a hypothetical healthy library with the same available metadata.
Together they can distinguish a coherent external-reference signal from a
target-study shift, an unresolved constituent, or a result dominated by one
person.

The main claim is about the 16 posteriorHCA-covered genes from the published
21-gene STING-activation signature. A constituent that is not supported by the
external audit is described as `not corroborated`; it is not automatically a
falsified published biomarker. CCL3, CCL4 and IL6 require a separate analysis
if the manuscript wishes to discuss the source-nominated individual
biomarkers.

### Figure 3: HCC GeoMx

The HCC case is deliberately narrow. It uses only:

- invasive front (IF) and tumour centre (TC);
- recurrent (REC) and non-recurrent (non-REC) clinical groups;
- 20 NK-enriched GeoMx ROI measurements;
- five patient × compartment profiles from three patients; and
- four posteriorHCA-covered TIMES genes: SPON2, VIM, ZFP36 and ZFP36L2.

The primary question is:

> Can comprehensive and metadata-conditioned healthy liver NK references
> provide useful external context for four spatial-expression features, while
> preserving the source-reported non-REC-versus-REC ordering in IF and TC?

This is evidence that posteriorHCA can be queried using ROI-aggregated spatial
expression profiles without using spatial coordinates. It is not evidence of
general performance across all spatial technologies or tissue programmes.

There is no healthy GeoMx bridge cohort. Therefore Figure 3 does not claim:

- HCA and GeoMx counts are calibrated to the same measurement scale;
- nominal posteriorHCA predictive intervals have GeoMx coverage guarantees;
- a predictive percentile is a disease probability;
- population-level recurrence prediction from three patients; or
- that ROI measurements are independent biological replicates.

---

## 2. Shared statistical terms

### 2.1 Biological subject

One human donor or patient. The subject is the inferential unit. Cells, SAVI
libraries collected repeatedly from one patient, and multiple GeoMx ROI
measurements do not create additional independent people.

### 2.2 GeoMx ROI measurement

One spatial expression measurement represented by one public GeoMx expression
column. The source files and current reports generally call these observations
ROIs. The main figure therefore uses `ROI measurement`, not `segment`.

### 2.3 Patient × compartment profile

All eligible GeoMx ROI measurements from one patient and one compartment are
aggregated into a single profile. Five such profiles are available:

| Profile | Outcome | ROI measurements |
|---|---|---:|
| A24 IF | non-REC | 2 |
| A6 IF | non-REC | 7 |
| A29 IF | REC | 4 |
| A24 TC | non-REC | 3 |
| A29 TC | REC | 4 |

A24 and A29 each contribute both IF and TC profiles. Those profiles remain
linked observations from one patient.

### 2.4 Comprehensive HCA baseline

The comprehensive query averages supported healthy-reference metadata over
one frozen HCA standardisation distribution.

- In SAVI, the cohort comparator is the healthy mean for an unseen study. One
  shared new-study effect is sampled in each posterior draw.
- In HCC, the comparator is the comprehensive healthy liver NK reference for
  each of the four genes. All cohort and patient quantities must use the same
  offset-centred latent-expression definition.

### 2.5 Metadata-conditioned HCA baseline

Available target metadata are fixed and unknown fields are integrated over
the frozen healthy-reference distribution.

- In SAVI, observation noise is included when simulating a new healthy count.
- In HCC, age and sex are conditioned for A24 (41F) and A6/A29 (57M). The
  current predictive-count output is descriptive across platforms and must not
  be labelled calibrated without a healthy GeoMx bridge.

Smooth age is a possible covariate representation inside either baseline. It
is not a third baseline. Age-decade remains provisional until every required
smooth-age fit passes the production convergence gate.

### 2.6 Density curve

A density curve is a smoothed display of many posterior or posterior-predictive
draws. Its area is scaled to one. The x-axis defines the quantity whose draws
are displayed; curve height is relative density, not sample count.

The legend must distinguish:

- a posterior density of a latent expected abundance;
- a posterior-predictive density of a new count; and
- the observed target value shown as a vertical mark.

### 2.7 SAVI healthy-prediction percentile

For observed count \(Y_{ig}^{obs}\) and \(M\) metadata-matched healthy
predictive counts,

\[
r_{ig}^{pred}=
\frac{0.5+\sum_m I(Y_{ig}^{rep}<Y_{ig}^{obs})
+0.5\sum_m I(Y_{ig}^{rep}=Y_{ig}^{obs})}{M+1}.
\]

A value of 0.90 means the observed count is greater than approximately 90% of
same-gene healthy simulations for that metadata query. It is not a 90%
probability of disease.

### 2.8 Comprehensive HCC-versus-HCA effect

For gene \(g\), compartment \(c\) and recurrence group \(r\), define

\[
\Delta_{gcr}
=\log\mu_{gcr}^{HCC}-\log\mu_g^{HCA}.
\]

Zero means the HCC cohort mean equals the comprehensive healthy liver NK
reference on the declared offset-centred latent-expression scale. A negative
value means lower expression in the HCC cohort; a positive value means higher.
The final interval must combine target-cohort and HCA-reference uncertainty in
the same draw table. The existing Welch tests that treat posterior draws as
independent healthy observations are not used.

### 2.9 Intervals and ranges

- `95% credible interval`: 2.5th–97.5th posterior percentiles for one defined
  Bayesian estimand.
- `95% subject-coupled stability envelope`: the SAVI 16-gene aggregate range
  that combines model draws and shared subject/donor reweighting. It is not an
  ordinary joint posterior interval because genes were fitted separately.
- `deletion range`: minimum-to-maximum posterior median over declared
  subject-deleted fits.
- `processing range`: minimum-to-maximum value over declared ROI aggregation
  or preprocessing alternatives.

### 2.10 Production computation gate

Every model contributing a final point, interval or density must satisfy:

- split-chain \(\widehat R\leq1.01\);
- zero divergent transitions;
- zero maximum-treedepth hits;
- bulk and tail effective sample size at least 400; and
- Monte Carlo standard error no greater than 5% of posterior standard
  deviation for every displayed estimand.

A failed fit is labelled `not estimable`; the model is not changed for only
the gene whose result is inconvenient.

---

## 3. Figure 2 — SAVI panel specification

### Figure-level logic

Panel a establishes the study observations and HCA support. Panel b checks
composition directions. Panel c is the signature-first external audit. Panel
d visualises comprehensive versus metadata-conditioned densities. Panel e
shows every library. Panel f tests subject and family influence.

### Panel 2a — Study design and reference support

**Purpose.** Show the actual cohort design, repeated observations, family
relationships and whether required metadata are represented in the frozen HCA
model.

**Left plot.** The x-axis is age at collection. The y-axis is SAVI sequencing
batch and chemistry. One symbol is one library; colour/shape separates study
control, untreated SAVI and JAK-inhibitor-treated SAVI. Lines join libraries
from one patient. C1/C2 are shown in a separate `Adult; exact age unavailable`
box. P5 uses age 10.5. Family annotations identify P1/P2/C8/C9.

**Right plot.** The x-axis is unique healthy HCA donor age. Stacked bars show
donor counts by sex, and orange ticks show exact SAVI ages. Printed totals are
generated from the frozen model frame. The final panel must also state whether
target-study overlap is absent, observed or unresolved.

**Analysis required.** Freeze a 17-library metadata table; verify 12 unique
people; audit HCA donor, sample, dataset, assay, sex and age support; and audit
whether SAVI or derivatives entered posteriorHCA training.

### Panel 2b — Composition directions

**Purpose.** Ask whether seven published composition directions are recovered
using patients as replicates and whether they also hold against HCA.

**Plot.** The y-axis lists naive CD4 T, effector/memory CD4 T, naive CD8 T,
effector/memory CD8 T, NK, MAIT and gamma-delta T. The x-axis is the difference
in expected proportion in percentage points. An open circle and line show
untreated SAVI minus study controls; a grey diamond and line show untreated
SAVI minus the comprehensive unseen-study HCA mean. The arrow beside each row
is the source-reported direction.

**Analysis required.** Construct subject-level counts for eight mutually
exclusive categories: the seven displayed categories plus `other`. Fit one
joint `sccomp` 2.2.0 multi-beta-binomial model with group and sequencing batch.
Regenerate the frozen HCA 26-L3 joint composition draws and aggregate within
the same draw so all displayed categories plus `other` sum to one. Do not
renormalise separately fitted compartments.

**Status.** The inferential axis is blank. Existing pooled means are not valid
patient-level intervals.

### Panel 2c — Covered signature and constituent audit

**Purpose.** Test the 16 posteriorHCA-covered signature genes as one
prespecified aggregate, then show how each constituent behaves.

**Columns.** The four aligned effects are:

1. untreated minus study control;
2. untreated minus comprehensive unseen-study HCA;
3. paired untreated minus treated in P1, P2, P4 and P6; and
4. study control minus comprehensive unseen-study HCA.

The x-axis is log2 fold change in offset-centred marginal expected abundance.
Zero means no difference. The top strip shows the equal-weight mean of the 16
gene-wise log2 effects. Its line is the 95% subject-coupled stability envelope.
Constituent rows show gene-specific posterior medians and 95% credible
intervals.

P1 has two treated observations, which are averaged before P1 contributes one
paired weight. P5 has no treated observation and does not enter the paired
contrast.

**Target model.** Fit one negative-binomial mixed model per gene with group,
sequencing batch and patient random intercept. Use one frozen prior and sampler
configuration for every gene. Evaluate group means over the same
subject-equal batch-standardisation grid. Build one joint draw table containing
`U_std`, `C_std`, `T_std`, `H_newstudy`, `U_paired` and `T_paired`.

**Cell alignment.** Recover the source scoring population and separately
define the HCA-aligned monocyte target. Cluster 17 alone cannot be labelled a
match to broad `V1_monocytic` without annotation evidence. Validate the
21-versus-16 source score before calling the covered subset a proxy for the
full signature.

**Status.** All result positions are blank until cell-aligned target models and
joint target/HCA draws are available.

### Panel 2d — Comprehensive and metadata-conditioned densities

**Purpose.** Make the two query modes visually explicit for three genes
selected by fixed, non-duplicating rules after panels c and f are complete.

The columns are selected in order:

1. strongest gene with positive untreated-minus-HCA and untreated-minus-control
   evidence;
2. largest remaining study-control-versus-HCA shift; and
3. largest remaining subject-deletion influence.

If a category has no eligible gene, its slot remains labelled blank.

**Upper row.** The x-axis is offset-zero log2 marginal expected count; the
y-axis is posterior density. Grey shows the comprehensive unseen-study healthy
mean. Green, vermilion and purple show study controls, untreated patients and
the descriptive treated-group mean on the same standardised scale.

**Lower row.** The x-axis is offset-centred log2 count; the y-axis is
predictive density. Five vertically separated teal ridges show healthy
predictions for P1, P2, P4, P5 and P6. Each ridge has the matched observed
untreated count mark. No single patient is selected as an exemplar.

### Panel 2e — Every-library healthy-prediction percentile

**Plot.** Rows are 16 genes and columns are all 17 libraries: seven controls,
five untreated and five treated. A tile is the SAVI healthy-prediction
percentile defined above. The right strip shows, separately for each gene, the
minimum, median and maximum across the seven study controls.

**Analysis required.** Correct P5 to 10.5 years. Integrate C1/C2 over a declared
adult-age distribution rather than imputing 35. Generate three batches of
4,000 joint predictive draws per gene × library. The batch-specific percentile
range must be no greater than 0.02; otherwise double every batch.

The control strip diagnoses target-study shift. Seven controls cannot estimate
a calibration curve.

### Panel 2f — Subject and family influence

**Plot.** The left facet is untreated/HCA log2 effect; the right facet is
untreated/study-control log2 effect. The first row is the 16-gene aggregate;
up to five constituent rows are selected by a frozen influence rule. A black
diamond and interval show the full-data result, open points show subject-
deleted medians, a thin line shows the deletion range, and purple identifies a
change among positive, negative and unresolved direction categories.

**Analysis required.** Refit the small target model after each of the 12
individual people is removed, for both effects. Add C8/C9 and
P1/P2/C8/C9 family removals in Extended Data. The HCA model remains fixed.
Evaluate every deleted fit on the full-data standardisation grid.

---

## 4. Figure 3 — HCC GeoMx panel specification

### Figure-level logic

Panel a exposes the sparse spatial design. Panel b verifies extraction of the
source non-REC-greater-than-REC direction. Panels c and d use the comprehensive
HCA reference at cohort level. Panel e shows metadata-conditioned densities
for every observed profile. Panel f isolates how much age/sex conditioning
moves the reference.

Only the IF/TC, REC/non-REC and four-gene quantities defined above enter this
case-study figure; unused source compartments and exploratory context branches
are excluded.

### Panel 3a — Actual observations and aggregation

**Purpose.** Show exactly how the spatial measurements used downstream become
biological profiles.

**Flow.** `20 GeoMx ROI measurements` → `aggregate within patient × IF/TC` →
`five spatial profiles from three patients`.

The right-hand matrix shows the number of ROI measurements:

- A24 41F non-REC: IF 2, TC 3;
- A6 57M non-REC: IF 7, no TC profile;
- A29 57M REC: IF 4, TC 4.

Blue denotes non-REC and vermilion denotes REC. The number inside each symbol
is an ROI count, not a patient count. The footer names SPON2, VIM, ZFP36 and
ZFP36L2 and states that patient is the biological unit.

**Required audit.** Freeze the exact ROI-to-profile manifest and verify the
five rows against
`Jie_HCC_IF_TC_patient_pseudobulk_metadata.csv`. Record aggregation, offset and
normalisation rules. HLA-DRB1 remains a source TIMES gene but is excluded from
the four-gene posteriorHCA figure because the required `V1_nk` model branch is
unavailable.

### Panel 3b — Source-direction extraction check

**Purpose.** Verify that the extracted four-gene GeoMx data preserve the
source-reported direction before posteriorHCA is applied.

**Plot.** Rows are the four genes. The x-axis is
`log2(non-REC / REC)` after patient-level aggregation. A blue circle is IF and
an orange square is TC. Zero means no group difference; a positive value is in
the source direction. The annotation reports the number of positive directions
out of eight.

**Sample sizes.** IF has two non-REC patients and one REC patient. TC has one
non-REC patient and one REC patient. No ROI-level P value is displayed.

**Current status.** The eight effects are populated from the audited current
IF/TC result table. This panel is a technical extraction check, not
posteriorHCA validation.

### Panel 3c — Comprehensive reference densities

**Purpose.** Show where each of the four HCC cohorts lies relative to the
comprehensive healthy liver NK reference.

**Plot.** Rows are the four genes. The x-axis is the shared offset-centred log
expected-abundance scale; the y-axis within each row is posterior density. The
grey density is the comprehensive HCA reference. Four directly labelled marks
or thin cohort densities show IF non-REC, IF REC, TC non-REC and TC REC.

All displayed reference and target quantities must use the same offset,
logarithm base, marginalisation policy and latent-mean definition. A posterior
density must be regenerated from saved draws; a Gaussian curve reconstructed
only from a mean and standard deviation is not acceptable final artwork.

**Status.** The blueprint is blank because native-vector density data still
need to be regenerated and scale-audited.

### Panel 3d — Comprehensive HCC-versus-HCA effects

**Purpose.** Provide the formal cohort-level posteriorHCA comparison for all
four genes and all four observed HCC cohorts.

**Plot.** Rows are the four genes. The x-axis is
`HCC cohort minus comprehensive HCA log expected abundance`. Four symbols per
row represent IF non-REC, IF REC, TC non-REC and TC REC. Colour denotes
recurrence and shape denotes IF or TC. Zero means equal to the comprehensive
healthy reference. A point is the joint-posterior median and its line is the
95% credible interval.

The current blueprint shows audited point estimates only. Final lines require
a new joint draw table pairing HCA reference draws with target-cohort
uncertainty. Do not reuse the existing Welch P values, and do not divide by the
posterior standard deviation of a grand mean and call that healthy population
variation.

Because TC is one patient versus one patient and the only REC patient is A29,
the panel is a proof-of-concept external audit, not population recurrence
inference.

### Panel 3e — Metadata-conditioned density grid

**Purpose.** Show all four genes in all five patient × compartment profiles,
without selecting one favourable patient.

**Grid.** Rows are SPON2, VIM, ZFP36 and ZFP36L2. Columns are A24 IF, A6 IF,
A29 IF, A24 TC and A29 TC. Every cell contains:

- a metadata-conditioned healthy predictive density;
- a vertical mark for the observed ROI-aggregated value; and
- optionally the observed percentile as directly labelled descriptive text.

The x-axis is an offset-centred scale declared once for the grid. Colour
identifies REC/non-REC; it does not recolour the healthy reference density.

**Analysis required.** Generate joint predictions using one shared new-dataset
effect per posterior draw; integrate unknown HCA metadata using frozen weights;
retain the same age/sex query for profiles sharing a patient; and save all draw
data before density estimation. Run ROI aggregation and offset sensitivities
without treating individual ROIs as independent patients.

Because no healthy GeoMx bridge exists, the panel is labelled a descriptive
cross-platform query. It contains no calibrated GeoMx P value or coverage
claim.

### Panel 3f — Patient-level predictive deviation and gene influence

**Purpose.** Summarise the personalised baseline across four genes and test
whether one gene alone creates the observed patient ordering.

For profile \(i\), convert each descriptive predictive percentile to its
standard-normal position and calculate

\[
S_i=\frac{1}{4}\sum_g\left|\Phi^{-1}(r_{ig}^{pred})\right|.
\]

**Plot.** Rows are A24 IF, A6 IF, A29 IF, A24 TC and A29 TC. The x-axis is
mean absolute predictive z. A coloured point is the full four-gene score. A
thin line is the minimum-to-maximum score after omitting each one of the four
genes in turn. The line is a leave-one-gene-out sensitivity range, not a
confidence interval. Blue is non-REC and vermilion is REC.

The score is descriptive because the healthy predictive distributions are not
calibrated in GeoMx and there are only three patients with one REC patient. No
AUC or population recurrence claim is shown in the main figure.

The current personalised outputs use age-decade. Smooth-age results replace
them only after the production gate is passed; the existing smooth-age pilot
has no passing configuration and therefore does not populate the figure.

---

## 5. Ordered remaining analysis

### Figure 2

1. Freeze the SAVI library, subject, family, age and cell-annotation manifest.
2. Audit expression and composition HCA support and target-study overlap.
3. Fit the eight-category patient-level composition model.
4. Resolve source and HCA-aligned cell populations.
5. Fit the 16 target expression models with one frozen specification.
6. Assemble the four gene-level effects and subject-coupled aggregate.
7. Generate all-library metadata-conditioned predictions and control
   diagnostics.
8. Run all individual and family deletions.
9. Select panel-2d genes using the frozen rules and generate vector densities.
10. Render from saved tables/draws and independently verify every plotted
    number.

### Figure 3

1. Freeze the 20-ROI-to-five-profile manifest and four covered genes.
2. Verify patient-first IF/TC aggregation, offsets and normalisation.
3. Regenerate the eight panel-3b extraction effects from the current IF/TC
   data.
4. Audit that HCA and HCC quantities used in panel 3c share the declared
   offset-centred latent-expression definition.
5. Generate comprehensive HCA and four-cohort density draws for panel 3c.
6. Build joint target/HCA draw tables and calculate the four cohort effects and
   credible intervals for panel 3d; do not reuse posterior-draw Welch tests.
7. Generate metadata-conditioned predictive draws for all 20 gene × profile
   cells in panel 3e.
8. Calculate the four-gene predictive-deviation score and every
   leave-one-gene-out sensitivity for panel 3f.
9. Run ROI aggregation, offset and age-formula sensitivities; retain
   age-decade until a smooth-age model passes production QC.
10. Render from saved draw tables and verify every plotted value.

---

## 6. Claim-to-evidence limits

| Proposed statement | Evidence required | If unmet |
|---|---|---|
| posteriorHCA adds an external population audit to SAVI | valid HCA support, aligned 16-gene untreated/HCA effect, study-control comparison and subject-deletion stability | say `provides external reference context` |
| the 16-gene covered subset is externally corroborated | both untreated/HCA and untreated/control aggregate stability envelopes above zero with no subject-deletion direction reversal | report the observed unresolved or subject-sensitive pattern |
| a SAVI constituent is not corroborated | its four panel-2c effects and panel-2f deletions fail the prespecified directional pattern | do not call it a falsified published biomarker |
| posteriorHCA is useful for the HCC spatial case | panel 3b preserves the source direction; panels 3c–f provide stable comprehensive and conditioned reference context for the four genes | describe the application as exploratory |
| HCC cohorts differ from comprehensive healthy HCA | panel-3d joint credible intervals using the same latent scale, with the one-REC-patient limitation | report audited point estimates as descriptive only |
| metadata-conditioned querying separates the observed profiles | panel-3e densities and panel-3f scores are stable to gene omission and the final age formula | say the personalised pattern is unresolved or model-sensitive |

No HCC panel alone supports calibrated GeoMx prediction or population-level
recurrence accuracy. The case demonstrates applicability of posteriorHCA
queries to ROI-aggregated spatial expression, within the limitations of three
patients and four model-covered genes.
