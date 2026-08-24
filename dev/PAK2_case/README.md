# PAK2 / IPF case-study data layer

Reproducible download and reconstruction of public count data for Watanabe et al., *Eur Respir J* 2025 (DOI [10.1183/13993003.00022-2025](https://doi.org/10.1183/13993003.00022-2025); GEO [GSE285246](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285246); control scRNA from [GSE173896](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173896)).

This stage builds a **lossless, analysis-ready count layer**. It does **not** QC-filter, normalise, integrate, cluster, annotate cells, or run DE.

## How to rerun

On University of Adelaide HPC, from this directory:

```bash
bash -lc 'use_R_env && Rscript scripts/run_all.R'
```

Or step by step: `00_download_geo.R` → `01_build_metadata.R` → `02_build_spatial_seurat.R` → `03_build_scrna_seurat.R` → `04_pseudobulk_spatial.R` → `05_pseudobulk_scrna.R`.

Existing valid files are skipped. The 11.9 Gb `GSE173896_COPD.rds.gz` and the remainder of the GSE173896 cohort are **not** downloaded.

## Data obtained

### GSE285246 spatial (9 Visium sections, 3 IPF donors)

| Donor | Sections | GSM accessions |
| --- | --- | --- |
| JKPF1 | JKPF1_1, JKPF1_2 | GSM8699165, GSM8699166 |
| JKPF2 | JKPF2_1, JKPF2_2, JKPF2_3 | GSM8699167–GSM8699169 |
| JKPF3 | JKPF3_1, JKPF3_2, JKPF3_3, JKPF3_5 | GSM8699170–GSM8699173 |

GEO processed files per section: `matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`, `tissue_positions_list.csv.gz`, `scalefactors_json.json.gz`, hi/low-res tissue images, fiducial and detected-tissue JPEGs.

### GSE285246 IPF scRNA (3 donors)

| Donor | GSM | Cells in GEO MTX |
| --- | --- | --- |
| JKPF1 | GSM8699174 | 5,571 |
| JKPF2 | GSM8699175 | 3,034 |
| JKPF3 | GSM8699176 | 8,163 |

### GSE173896 control scRNA (3 never-smoker donors only)

| Donor | GSM | Cells in GEO MTX |
| --- | --- | --- |
| JK06 | GSM5282546 | 3,608 |
| JK11 | GSM5282547 | 4,721 |
| JK12 | GSM5282548 | 4,775 |

COPD, non-COPD smoker, and later JK19–JK29 samples were not downloaded.

## Files successfully reconstructed

```
9 spatial Seurat objects          data/seurat/spatial/JKPF*.rds
1 spatial list                    data/seurat/spatial/spatial_seurat_list.rds
6 scRNA Seurat objects            data/seurat/scrna/{JKPF1,JKPF2,JKPF3,JK06,JK11,JK12}.rds
1 scRNA list                      data/seurat/scrna/scrna_seurat_list.rds
1 merged count-level scRNA object data/seurat/scrna/scrna_seurat_merged_counts.rds
1 spatial section-level SE        data/pseudobulk/spatial/spatial_by_section_SE.rds
1 scRNA donor-level SE            data/pseudobulk/scrna/scrna_by_donor_unannotated_SE.rds
```

**Not created (annotations unavailable):**

- `data/pseudobulk/spatial/spatial_by_region_SE.rds`
- `data/pseudobulk/scrna/scrna_by_donor_celltype_SE.rds`

The donor-level scRNA `SummarizedExperiment` is a **data-integrity object only**. It is **not** suitable for posteriorHCA cell-type-specific comparison.

Assay `counts` in every pseudobulk object is the raw integer UMI sum. Feature rownames are Ensembl IDs (`ENSG…`). Spatial images are attached on all nine Visium objects.

## Metadata availability

| Item | Located? | Source |
| --- | --- | --- |
| FF/DF spot annotations | **No** | Not in GEO processed files, ERJ supplementary tables, or Figshare Tables 4/6/7 |
| Original scRNA cell-type annotations (Watanabe IPF fibroblast subtypes) | **No** | Not in GSE285246 MTX, ERJ tables, or GSE173896 COPD metadata |
| Age | **Yes** | ERJ Supplementary Table 1 |
| Sex | **Yes** | ERJ Supplementary Table 1 (column is misspelled `Gendar`) |
| Smoking status | **Yes** | ERJ Supplementary Table 1 `Smoking cessation` plus GEO “Never smokers” for JK06/JK11/JK12 |
| Pack-years, BMI | **Yes** | ERJ Supplementary Table 1 |

Donor table (from `data/metadata/donor_metadata.csv`):

| donor_id | condition | disease | source_GSE | sex | age | smoking_status | smoking_cessation | pack_years |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| JKPF1 | IPF | IPF | GSE285246 | Male | 76 | former | 41 years | 5 |
| JKPF2 | IPF | IPF | GSE285246 | Male | 69 | former | 2 years | 47 |
| JKPF3 | IPF | IPF | GSE285246 | Male | 80 | former | 25 years | 35 |
| JK06 | control | healthy | GSE173896 | Male | 58 | never | Never | 0 |
| JK11 | control | healthy | GSE173896 | Female | 56 | never | Never | 0 |
| JK12 | control | healthy | GSE173896 | Male | 33 | never | Never | 0 |

`condition` and `source_GSE` are stored separately. `disease = healthy` for the GSE173896 donors comes from the IPF paper (“three healthy controls”), not from GEO’s “Never smokers” characteristic.

`smoking_status = former` for IPF donors is taken from Table S1 column **Smoking cessation** (years since cessation). The raw string is kept in `smoking_cessation`.

GSE173896 `GSE173896_COPD_meta_new.xlsx` was attached to matching control barcodes as `copd_study_class` / `copd_study_seurat_clusters` / `copd_study_DoubletFinder`. **`Class` is COPD / NC / Never, not a cell type.** Those labels were not used for pseudobulk.

## Missing information

- No barcode-level FF/DF (or other pathological-region) table. The paper reports 97 FF and 412 DF spots identified histologically; those barcode lists were not deposited with GEO or the ERJ supplementary tables.
- No barcode-level Watanabe cell-type labels (alveolar / adventitial / inflammatory / CTHRC1+ / WNT5A+CTHRC1+ myofibroblast, etc.).
- GEO sample records for GSE285246 only list `tissue: IPF lung` (no age/sex there; those come from Table S1).

## Important limitations

1. **Spatial section/spot-level expression** is available (raw counts + images + coordinates). Pathological FF/DF labels are **not**, so region-level FF-vs-DF pseudobulk was not built and must not be inferred from images or expression in this layer.
2. **scRNA is currently donor-level only.** Donor × cell-type pseudobulk requires author (or later) cell annotations. Individual-cell Seurat objects are retained for that later stage.
3. **GSE285246 and GSE173896 are different source studies** (Cell Ranger 6.1.2 vs 3.0.2; IPF vs never-smoker control lungs). Feature IDs happen to be the same 33,538 GRCh38 Ensembl genes, but that does not remove study-level batch or protocol differences.
4. GEO MTX cell/spot counts match ERJ Supplementary Tables S2/S3 exactly. This pipeline did not apply additional QC filters; whatever filtering is already baked into the deposited processed matrices is retained as-is.

## Feature identifiers

Both studies use 10x GRCh38 `features.tsv` with Ensembl ID, gene symbol, and `Gene Expression`.

- Unique Ensembl IDs: **33,538**
- Shared between GSE285246 and GSE173896: **33,538** (no study-private genes)
- Duplicate gene symbols: **24** (kept as distinct Ensembl rows; symbols were not collapsed)

See `data/metadata/feature_map_scrna.csv`. In combined pseudobulk, a gene absent from a source matrix would be stored as `NA` rather than zero; that case did not occur here.

## Validation (structural)

Spatial (all spots kept; images attached):

| section | donor | n_spots | n_features | total UMI |
| --- | --- | --- | --- | --- |
| JKPF1_1 | JKPF1 | 2899 | 33538 | 10,482,320 |
| JKPF1_2 | JKPF1 | 3642 | 33538 | 21,505,200 |
| JKPF2_1 | JKPF2 | 2160 | 33538 | 13,768,584 |
| JKPF2_2 | JKPF2 | 1982 | 33538 | 13,245,679 |
| JKPF2_3 | JKPF2 | 1452 | 33538 | 8,951,519 |
| JKPF3_1 | JKPF3 | 2413 | 33538 | 15,338,392 |
| JKPF3_2 | JKPF3 | 2386 | 33538 | 11,730,543 |
| JKPF3_3 | JKPF3 | 2459 | 33538 | 11,478,140 |
| JKPF3_5 | JKPF3 | 2130 | 33538 | 12,467,614 |

Section-level pseudobulk: `sum(counts)` per section equals that section’s Seurat total UMI. Across sections, 118,967,991 UMIs.

scRNA:

| donor | condition | GEO | n_cells | n_features | total UMI | author cell types |
| --- | --- | --- | --- | --- | --- | --- |
| JKPF1 | IPF | GSE285246 | 5571 | 33538 | 37,526,577 | 0 / 5571 |
| JKPF2 | IPF | GSE285246 | 3034 | 33538 | 28,291,256 | 0 / 3034 |
| JKPF3 | IPF | GSE285246 | 8163 | 33538 | 112,839,081 | 0 / 8163 |
| JK06 | control | GSE173896 | 3608 | 33538 | 32,467,101 | 0 / 3608 |
| JK11 | control | GSE173896 | 4721 | 33538 | 38,090,601 | 0 / 4721 |
| JK12 | control | GSE173896 | 4775 | 33538 | 43,616,850 | 0 / 4775 |

Donor-level scRNA pseudobulk UMI totals match the source objects. No cells were dropped.

## Layout

```
PAK2_case/
├── README.md
├── scripts/          # 00–05 + utils.R + run_all.R
├── data/
│   ├── raw/          # GEO MTX/images, SOFT, ERJ supplements
│   ├── metadata/     # manifests and donor/sample tables
│   ├── seurat/
│   └── pseudobulk/
└── logs/
```

R session details: `logs/sessionInfo.txt`.
