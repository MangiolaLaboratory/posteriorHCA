# File guide

The folder is organized by provenance:

- `data/raw/`: downloaded source material; never use these files directly as
  patient replicates without the sample manifest.
- `data/metadata/`: mappings, annotations, QC tables, and descriptions of matrix
  columns.
- `data/processed/`: R objects containing reconstructed single-cell or
  patient-level count data.
- `scripts/`: the reproducible workflow.
- `logs/`: one execution log per workflow stage.

## Which files contain expression data?

### Primary downloaded expression data

`data/raw/dge/*.txt.gz` are the 67 per-library raw digital gene-expression
matrices downloaded from GEO. Each file has gene symbols in rows, cell barcodes
in columns, and raw integer UMI counts as values. A file is a sequencing library,
not necessarily a patient: technical replicates and tumour pieces must be joined
through `GSE293189_sample_manifest.tsv`.

`data/raw/GSE293189_RAW.tar` is the original GEO archive containing the same 67
DGE files. It is retained for provenance; `data/raw/dge/` is the extracted copy
used by the workflow.

### Reconstructed single-cell expression data

| File | R class/content | Dimensions | Values and intended use |
| --- | --- | ---: | --- |
| `GSE293189_raw_counts.rds` | `SingleCellExperiment` | 36,994 genes x 56,192 deposited barcodes | Complete unfiltered raw integer-valued UMI matrix. `rowData` contains gene identifiers and `colData` contains library/patient/QC metadata. This is the archival reconstructed object. |
| `GSE293189_qc_counts.rds` | list: `counts`, `cell_metadata`, `feature_map` | 36,994 x 26,402 cells | Raw counts after gene, UMI, mitochondrial, and doublet filters. It includes the two paraganglioma subjects as well as the bladder cohort. It is deliberately not normalized. |
| `GSE293189_processed_seurat.rds` | Seurat object | 36,994 x 21,782 bladder cells | HV/UC bladder cohort only. The `counts` layer is raw; the `data` layer is LogNormalized. Also contains 2,000 HVGs, PCA, graph clusters, UMAP, module scores, broad annotations, and tumour-epithelial flags. Use this to inspect or reproduce the single-cell analysis. |
| `GSE293189_all_cell_qc_metadata.rds` | data frame | 56,192 rows | QC measurements and pass/fail calls for every deposited barcode. This file contains no expression matrix. |
| `GSE293189_infercnv_input.rds` | list | genes x selected cells | Raw counts and reference/tumour labels prepared for InferCNV. It is input only: no InferCNV HMM result is claimed. |

### Patient-level pseudobulk expression data

These are the preferred inputs for patient-level differential-expression or
TM4SF1 analyses. Each is a `SummarizedExperiment`; `assay(x, "counts")` is a
sparse gene-by-biological-sample matrix produced by summing raw UMIs. No
normalized expression was averaged.

| File | Dimensions | One column represents | Important scope |
| --- | ---: | --- | --- |
| `pseudobulk_tumour_epithelial_counts.rds` | 36,994 genes x 12 patients | patient x tumour x epithelial | The published 9 HV + 3 pure-UC patient IDs. UC01 is retained but flagged as failing the reconstructed 150-cell rule. |
| `pseudobulk_paired_normal_epithelial_counts.rds` | 36,994 x 7 patients | patient x paired-normal x epithelial | All patients with at least one paired-normal epithelial cell. `sufficient_cells_ge20` is true for three patients. These are paired/adjacent normal, not fully healthy samples. |
| `pseudobulk_all_broad_celltypes_counts.rds` | 36,994 x 80 profiles | patient x sample type x broad cell type | Epithelial, endothelial, stromal, and immune profiles from the bladder cohort. Paraganglioma cells are excluded because they were outside the paper's bladder analysis and were not broadly annotated. |

For every pseudobulk object, `colData` records patient ID, original ID, paper ID,
histology, HV/UC group, sample type, broad cell type, number of contributing
cells, and pseudobulk library size. `rowData` maps gene symbols to Ensembl and
Entrez IDs. Although sparse R matrices store their nonzero slot as numeric,
every count was checked to be integer-valued.

## How to open the main data objects

```r
library(SummarizedExperiment)

# Patient-level tumour epithelial raw counts
pb <- readRDS("data/processed/pseudobulk_tumour_epithelial_counts.rds")
counts <- assay(pb, "counts")
sample_information <- as.data.frame(colData(pb))
gene_information <- as.data.frame(rowData(pb))

# Complete reconstructed raw single-cell object
raw <- readRDS("data/processed/GSE293189_raw_counts.rds")
raw_counts <- assay(raw, "counts")

# QC-passing count-only object
qc <- readRDS("data/processed/GSE293189_qc_counts.rds")
qc_counts <- qc$counts
qc_cell_information <- qc$cell_metadata

# Seurat reconstruction: raw versus normalized values
library(SeuratObject)
seu <- readRDS("data/processed/GSE293189_processed_seurat.rds")
raw_counts <- LayerData(seu, assay = "RNA", layer = "counts")
log_normalized <- LayerData(seu, assay = "RNA", layer = "data")
```

## Metadata files

### Sample and patient identity

- `GSE293189_sample_manifest.tsv`: authoritative 67-row mapping from GSM and
  library to patient, original/paper ID, histology, tumour/paired-normal status,
  tumour piece, technical replicate, DGE filename, SRA run, checksum, and local
  path. It preserves both sides of conflicting GEO/filename mappings.
- `GSE293189_patient_manifest.tsv`: one row per biological subject, summarizing
  the 15 bladder patients and two out-of-scope paraganglioma subjects, library
  counts, and published cohort status.
- `metadata_mapping_inconsistencies.tsv`: the eight GSM title, description, and
  filename conflicts that require caution; none were silently resolved.
- `download_manifest.tsv`: download URL, expected/observed size, local path, and
  MD5 checksum for major downloaded artifacts.
- `authors_code_commit.txt`: exact Git commit of the downloaded authors' public
  analysis repository.

### Gene and cell annotation

- `GSE293189_feature_map.tsv`: one row per union gene symbol, with Ensembl ID,
  Entrez ID, and mitochondrial flag. TM4SF1 maps to ENSG00000169908.
- `GSE293189_cell_metadata.tsv.gz`: all 56,192 deposited barcodes, including
  patient/library identity, sequential QC flags, doublet information, and—for
  processed bladder cells—cluster, UMAP, module scores, broad cell class,
  annotation confidence/evidence, and tumour-analysis inclusion status.
- `tumour_epithelial_cell_metadata.tsv.gz`: the 8,021 bladder cells classified
  as tumour epithelial before applying the per-patient 150-cell rule.
- `broad_cluster_annotation.tsv`: one row per Seurat cluster with broad class,
  class-score margin, confidence, marker evidence, and class module scores.

### QC and published-cohort comparison

- `qc_cell_counts_overall.tsv`: overall cell attrition at every sequential QC
  step.
- `qc_cell_counts_by_library.tsv`: the same attrition separately for each of 67
  libraries.
- `qc_cell_counts_by_patient.tsv`: the same attrition after collapsing libraries
  to biological patients.
- `tumour_epithelial_counts_by_patient.tsv`: reconstructed tumour-epithelial
  cells and reconstructed/published inclusion flags for all 15 bladder patients.
- `published_vs_reconstructed_tumour_cells.tsv`: published Figure 1B counts
  versus reconstructed counts for the expected 12 retained patients.
- `validation_checks.tsv`: 31 machine-readable checks of file existence,
  dimensions, count integerness, metadata alignment, library-size sums, and cell
  accounting.
- `validation_passed.txt`: short success stamp from the latest validation run.
- `sessionInfo.txt`: R and package versions used when the report was generated.

### Pseudobulk column descriptions

- `pseudobulk_tumour_epithelial_sample_metadata.tsv`: the 12 columns of the
  tumour-epithelial matrix, including reconstructed and published cell counts.
- `pseudobulk_paired_normal_epithelial_sample_metadata.tsv`: the seven columns
  of the paired-normal epithelial matrix and the cell-sufficiency flag.
- `pseudobulk_sample_metadata.tsv`: the 80 columns of the broad-cell-type matrix.

### TM4SF1 summaries

- `TM4SF1_patient_summary.tsv`: TM4SF1 raw count, total library size, CPM,
  detected-cell count, and detection rate for every patient x sample type x
  broad-cell-type pseudobulk profile.
- `TM4SF1_group_summary.tsv`: descriptive summaries by sample type, broad cell
  class, and histology group across the full annotated bladder cohort.
- `TM4SF1_population_source_summary.tsv`: pooled epithelial, endothelial,
  stromal, and immune TM4SF1 abundance and each population's share of counts,
  separately for tumour and paired normal.
- `TM4SF1_published_tumour_cohort_by_patient.tsv`: the TM4SF1 rows for the exact
  12 published tumour patient IDs.
- `TM4SF1_published_tumour_cohort_group_summary.tsv`: descriptive HV-versus-UC
  totals and pooled CPM for those 12 IDs, with reconstructed threshold counts.

## Other raw source files

- `GSE293189_family.soft.gz` and `GSE293189_family.xml.tgz`: GEO series/GSM
  metadata in SOFT and MINiML formats.
- `GSE293189_filelist.txt`: GEO's official list of supplementary files.
- `PRJNA1243332_SraRunInfo.csv`: complete SRA/FASTQ run metadata. It is metadata,
  not downloaded FASTQ sequence.
- `PRJNA1243332_sra_esearch.json`: SRA project query result used to obtain the
  run table.
- `PMC12174346_fullTextXML.xml`: machine-readable full paper text.
- `PMC12174346_supplementaryFiles.zip`: original paper supplement archive.
- `paper/*.pdf`: paper supplementary information files.
- `paper/*MOESM4_ESM.xlsx`: source-data workbook; Figure5B contains the
  published 8,553 tumour-epithelial cells and per-cell TM4SF1 values, while
  Figure5A contains the released HV-versus-UC differential-expression table.
- `paper/*Fig*_HTML.{jpg,gif}`: downloaded paper figure images.
- `Histological_Variant_Bladder_Cancer_Analysis/`: verbatim clone of the authors'
  public code. Its scripts cover whole-bladder preprocessing (`Bladder_Seqwell`),
  epithelial/tumour analysis (`Epithelial_Seqwell`), InferCNV preparation,
  Cluster 13, correlation, dividing-cell, GO, subtype DEG, basal/luminal, and
  monocle analyses. Several scripts refer to unpublished intermediate RDS files,
  so they are reference logic rather than a directly runnable workflow here.

## Workflow and report files

- `scripts/00_download.R`: downloads and verifies GEO, SRA, paper, and code
  inputs.
- `scripts/01_build_manifest.R`: parses metadata and builds audited library and
  patient mappings.
- `scripts/02_qc_reconstruct.R`: reads DGE matrices, constructs the union raw
  matrix, applies sequential QC, and performs reproducible doublet handling.
- `scripts/03_cluster_annotate.R`: runs Seurat preprocessing/clustering, broad
  multi-marker annotation, tumour selection, and InferCNV-input preparation.
- `scripts/04_pseudobulk_tm4sf1.R`: sums raw counts to patients and produces
  pseudobulk and TM4SF1 summaries.
- `scripts/05_report.R`: regenerates `data/processed/qc_summary.md`.
- `scripts/06_validate_outputs.R`: performs the 31 final integrity checks.
- `scripts/run_all.R`: executes stages 00 through 06 in order.
- `scripts/07_replication_audit.R`: traces UC01 through QC, annotation, GEO/Git
  mappings, source-data barcodes, and public InferCNV feasibility.
- `scripts/08_tumour_replication_de.R`: performs audit-only tumour reclustering,
  Cluster 13 investigation, fixed-seed HV-versus-UC DE, and 100-seed TM4SF1
  downsampling sensitivity analyses.
- `scripts/09_validate_replication.R`: checks the published inclusion/exclusion
  logic, 9+3/8,553 targets, barcode provenance, patient-level pseudobulk, and
  raw-count integrity without hiding known replication failures.
- `scripts/10_replication_report.R`: generates the evidence-backed replication
  discrepancy report and public-method comparison table.
- `scripts/run_replication_audit.R`: executes stages 07 through 10 without
  rewriting the canonical processed data.
- `scripts/utils.R`: common path, download, parsing, sparse-matrix, and mapping
  helpers.
- `data/processed/qc_summary.md`: scientific validation report and main findings.
- `logs/*.log`: concise run logs for the seven workflow stages.

## Publication replication audit

`results/replication_audit/` is a separate, non-destructive evidence layer. Its
main entry point is `replication_discrepancy_report.md`. The folder also contains
UC01 cell-flow and marker tables, source-barcode matching, published-versus-
reconstructed patient counts, InferCNV feasibility, audit-only tumour
reclustering, fixed-seed and 100-seed DE sensitivity results, plots, explicit
replication validation, and checksums of the pre-audit canonical reconstruction.

The audit concludes that the public DGE deposit does not contain the 605 UC01
barcodes used in the paper. Consequently, no corrected canonical pseudobulk was
created and the available reconstruction must not be described as an exact
9 HV + 3 UC / 8,553-cell reproduction.

## Expanded TM4SF1 expression-overlap report

- `TM4SF1_toxicity_case_study.qmd`: the self-contained Quarto analysis. It
  re-reads Figure5A/Figure5B, visualises the public-DGE patient contrast, aligns
  patient pseudobulks to posteriorHCA reference samples with TMMwsp, estimates
  SAVI-style cohort latent `log(mu)`, builds separate Normal/10x3-supported
  epithelial and endothelial references, and applies a support-aware prototype
  overlap rule across tissue, sex, age decade and ethnicity.
- `results/tm4sf1_toxicity_report/TM4SF1_toxicity_case_study.md`: knitted
  Markdown preview of the executed report. The QMD remains the authoritative
  source.
- `results/tm4sf1_toxicity_report/TM4SF1_paper_DE_replication.tsv`: the released
  paper TM4SF1 result and two public-DGE sensitivity reconstructions.
- `results/tm4sf1_toxicity_report/TM4SF1_paper_source_patient_summary.tsv`:
  patient summaries of the released normalized Figure5B values; these are not
  raw UMI counts.
- `results/tm4sf1_toxicity_report/TM4SF1_study_logmu_draws.rds`: patient-level
  Bayesian-bootstrap cohort `log(mu)` draws for both model-reference scales.
- `results/tm4sf1_toxicity_report/TM4SF1_cohort_logmu_contrasts.tsv`: study
  cohort versus posteriorHCA comparisons on the log2 scale.
- `results/tm4sf1_toxicity_report/TM4SF1_threshold_registry.tsv`: the reusable
  threshold contract, including the cell-type model, TMM reference, HV tumour-
  expression anchor, lower cohort-location uncertainty quantile, prototype
  probability, donor support, diagnostic gate, fingerprints and prediction
  policy.
- `results/tm4sf1_toxicity_report/TM4SF1_alert_rule_sensitivity.tsv`: number of
  ungated prototype tissue candidates while varying the HV anchor quantile
  (1%, 5%, 10%), minimum donor support (5, 10, 20) and posterior exceedance
  probability (50%, 80%, 95%).
- `results/tm4sf1_toxicity_report/posteriorHCA_TM4SF1_profile_summary.tsv`: all
  Normal-labelled tissue/demographic posterior summaries. Unsupported rows are
  retained and explicitly marked; subgroup averages are descriptive observed-
  profile summaries rather than controlled demographic effects.
- `results/tm4sf1_toxicity_report/posteriorHCA_TM4SF1_supported_profiles.tsv`:
  the subset meeting the minimum dataset × donor support rule.
- `results/tm4sf1_toxicity_report/TM4SF1_prototype_overlap_candidates.tsv`:
  supported profiles passing the overlap rule before the model-readiness gate.
  These are exploratory follow-up candidates, not toxicity calls.
- `results/tm4sf1_toxicity_report/TM4SF1_expression_overlap_alerts.tsv`:
  readiness-gated alerts. This is empty by design while the current fits have
  `model_ready = FALSE`.
- `results/tm4sf1_toxicity_report/TM4SF1_model_diagnostics.tsv`: convergence and
  sampling diagnostics for the current TM4SF1 brms fits.
- `results/tm4sf1_toxicity_report/TM4SF1_reference_alignment.tsv`: model-specific
  reference samples, shared feature counts and fitted TM4SF1 dispersion.
- `results/tm4sf1_toxicity_report/figures/`: reusable PNG plots from every major
  report section.
- `results/tm4sf1_toxicity_report/cache/`: versioned per-cell-type computation
  caches; set `refresh_hca: true` in the QMD parameters to rebuild them.
- `results/tm4sf1_toxicity_report/output_manifest.tsv`: machine-readable report
  output inventory.
- `results/tm4sf1_toxicity_report/sessionInfo.txt`: R and package versions used
  for the expression-overlap report.

## Restart-only intermediates

- `intermediate/GSE293189_pre_doublet_checkpoint.rds`: uncompressed checkpoint
  of the union raw matrix, cell metadata, and feature map used to restart the
  expensive DGE-ingestion stage. It is redundant with final data and is not an
  analysis deliverable.
- `intermediate/GSE293189_scDblFinder_score_table.rds`: saved artificial-doublet
  score table used to reproduce the final doublet calls.
