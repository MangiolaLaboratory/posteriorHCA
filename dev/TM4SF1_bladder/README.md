# TM4SF1 bladder-cancer case study

This directory contains a restartable reconstruction from the public GEO DGE
matrices through patient-level raw-count pseudobulk. The canonical preprocessing
stops before posteriorHCA integration; the separate report described below adds
an exploratory, model-gated expression-overlap analysis without defining a
biological or clinical toxicity threshold.

Run from this directory with the project R environment:

```bash
/home/a1237163/miniconda3/envs/R_env/bin/Rscript scripts/run_all.R
```

The generated validation report is
[`data/processed/qc_summary.md`](data/processed/qc_summary.md). The exact
67-library mapping is in
[`data/metadata/GSE293189_sample_manifest.tsv`](data/metadata/GSE293189_sample_manifest.tsv),
including eight explicit GEO-title/description/filename conflicts.

See [`FILE_GUIDE.md`](FILE_GUIDE.md) for a description of every data object,
metadata table, source-file group, and workflow script, including R loading
examples.

The publication-discrepancy audit can be rerun after the canonical workflow:

```bash
/home/a1237163/miniconda3/envs/R_env/bin/Rscript scripts/run_replication_audit.R
```

Its main report is
[`results/replication_audit/replication_discrepancy_report.md`](results/replication_audit/replication_discrepancy_report.md).
The audit is non-destructive and does not rewrite the canonical processed
objects or pseudobulk matrices.

All pseudobulk values are sums of raw UMI counts at the patient level. Paired
normal samples mean paired/adjacent normal bladder and are not labelled healthy.
Restart-only QC checkpoints are kept under `data/processed/intermediate/` and
are not analytical deliverables.

## TM4SF1 expression-overlap report

The expanded case-study report is
[`TM4SF1_toxicity_case_study.qmd`](TM4SF1_toxicity_case_study.qmd). It
reproduces the paper's TM4SF1 differential-expression highlight, compares the
patient-level study cohorts with model-aligned posteriorHCA epithelial and
endothelial references using the SAVI latent `log(mu)` method, and writes a
model-specific threshold registry plus descriptive tissue/sex/age/ethnicity
profile summaries.

Render it from this directory with Quarto:

```bash
quarto render TM4SF1_toxicity_case_study.qmd
```

All helper functions are contained in the QMD. Computed tables, figures and
versioned caches are under `results/tm4sf1_toxicity_report/`; a knitted Markdown
preview is saved there as `TM4SF1_toxicity_case_study.md`. The current local
posteriorHCA fits fail the registered readiness gate, so ungated prototype
candidates are exported separately and the model-ready alert file is empty by
design. Neither file represents demonstrated toxicity.
