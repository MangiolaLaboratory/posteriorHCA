# GSE293189 replication discrepancy audit

Generated: 2026-08-18 00:17:28 ACST

## Outcome

Across the 12 patient IDs in the published cohort, the public deposit yields **7,862 reconstructed tumour epithelial cells**, including 90 UC01 cells that fail the 150-cell rule. Applying that rule leaves **9 HV + 2 UC patients / 7,772 cells**, not the published 9 HV + 3 UC / 8,553-cell cohort. Across all 15 tumours, including the three expected exclusions, the reconstruction contains 8,021 candidate tumour epithelial cells. The discrepancy is localized primarily to UC01: 90 reconstructed cells versus 605 published cells. No canonical object or pseudobulk matrix was changed, because neither relabelling non-epithelial cells nor lowering the 150-cell threshold would be scientifically justified.

The strongest provenance result is exact: **0 of 605** Figure5B UC01 source barcodes occur in the currently filename-assigned 12049 DGE, and 0 occur in any library assigned to current UC01. Only 3 barcode strings occur anywhere among all 67 deposited matrices, each as an unrelated single-library collision. In contrast, the recoverable source barcodes for UC02 and UC03 map back to their assigned deposited libraries. This is evidence of a missing or different UC01 input matrix, not a routine QC or annotation error.

## 1. Cell flow

### Bladder cohort overall

| stage | cells |
| --- | --- |
| Raw bladder barcodes | 48,842 |
| >=300 detected genes | 23,947 |
| >=500 UMIs | 23,608 |
| <20% mitochondrial reads | 21,970 |
| After doublet filter | 21,782 |
| Processed bladder cells | 21,782 |
| Broad epithelial |  8,299 |
| Candidate tumour epithelial |  8,021 |
| Final tumour epithelial |  8,021 |

The authors report 21,533 post-QC cells globally. The public reconstruction has 21,782 bladder cells after doublet filtering (249 more). Starting from the same 21,970 post-mitochondrial cells, matching 21,533 would require removal of 437 cells rather than 188. Missing DoubletFinder parameters can plausibly explain much of this global difference, but cannot explain why all 605 published UC01 barcodes are absent.

### UC01/current 12049 library

| stage | cells |
| --- | --- |
| Raw bladder barcodes | 1373 |
| >=300 detected genes | 1224 |
| >=500 UMIs | 1198 |
| <20% mitochondrial reads | 1153 |
| After doublet filter | 1142 |
| Processed bladder cells | 1142 |
| Broad epithelial | 90 |
| Candidate tumour epithelial | 90 |
| Final tumour epithelial | 90 |

The 1,142 UC01 QC-passing cells consist of 90 epithelial, 877 immune, 134 stromal, and 41 endothelial cells. Cluster 32 supplies 88 of the 90 epithelial cells and has moderate multi-marker confidence. Its epithelial support includes EPCAM detection in 73.9% of cells and KRT8/KRT18/KRT19 detection in 68.2%/71.6%/65.9%. Major excluded clusters show coherent alternative identities (for example PTPRC in 94.3% of cluster 10, ACTA2 in 100% of cluster 21, and VWF in 96.9% of cluster 14). There is therefore no hidden 515-cell epithelial population that can be recovered by a reasonable broad-marker reinterpretation.

Full evidence: `uc01_cell_flow.tsv`, `uc01_cluster_annotation.tsv`, `uc01_marker_summary.tsv`, `uc01_umap_annotations.pdf`, and `uc01_marker_featureplots.pdf`.

## 2. UC01 identity and metadata conflict

| GSM | library_id | patient_id | original_patient_id | geo_title | geo_original_patient_id | uc01_claim_basis | downloaded |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSM8878368 | PA_U67-12049_Pool_1_2_3_S20_L001 | UC01 | 12049 | Patient 11734 (Plasmacytoid) Tumor | 11734 | DGE_filename_12049 | TRUE |
| GSM8878370 | PA_U67-12050T1_Pool_1_2_3_S72_L003 | VAR10 | 12050 | Patient 12049 (Pure UC) Tumor | 12049 | GEO_title_or_description_12049 | TRUE |
| not_deposited | 12923 | UC01 | 12923 |  |  | authors_pre2026_code_12923 | FALSE |

The original public code mapped UC1 to patient 12923. Commit `d26f0427481c5ea39d1697ec3c83b7921a747504` (16 January 2026) changed only this assignment to 12049, while downstream scripts still refer to `PureUC12923`. Patient 12923 is not represented by a deposited DGE filename. GEO record GSM8878368 carries a 12049 filename but a 11734 title/description; GSM8878370 has a 12049 title/description but carries a 12050T1 DGE filename. These conflicts are preserved, not silently resolved.

## 3. Published versus reconstructed tumour cohort

| patient_id | histology_group | n_final_tumour_epithelial | published_tumour_epithelial_cells | delta_tumour_epithelial | reconstructed_ge150 |
| --- | --- | --- | --- | --- | --- |
| VAR06 | HV | 606 | 645 | -39 | TRUE |
| UC04 | UC | 24 |  |  | FALSE |
| VAR11 | HV | 38 |  |  | FALSE |
| VAR01 | HV | 317 | 358 | -41 | TRUE |
| VAR07 | HV | 699 | 848 | -149 | TRUE |
| VAR09 | HV | 2880 | 2917 | -37 | TRUE |
| VAR05 | HV | 497 | 428 | 69 | TRUE |
| UC03 | UC | 193 | 163 | 30 | TRUE |
| VAR03 | HV | 809 | 812 | -3 | TRUE |
| VAR02 | HV | 478 | 479 | -1 | TRUE |
| VAR08 | HV | 564 | 568 | -4 | TRUE |
| VAR04 | HV | 341 | 342 | -1 | TRUE |
| UC01 | UC | 90 | 605 | -515 | FALSE |
| VAR10 | HV | 97 |  |  | FALSE |
| UC02 | UC | 388 | 388 | 0 | TRUE |

All nine expected HV patients pass the 150-cell threshold, as do UC02 and UC03. UC01 does not. VAR10, VAR11, and UC04 remain correctly excluded. Validation records this as a publication-replication failure rather than changing the threshold or treating libraries as patients.

## 4. Authors' code and method audit

| component | paper_or_public_code | audit_implementation | exact_reproduction | implication |
| --- | --- | --- | --- | --- |
| Raw matrix assembly | Public script starts from private dge_mergedall_032123.rds and omits merge code | Read and union all 67 deposited raw DGE matrices | FALSE | The public DGE deposit is the only reproducible starting point |
| Ambient RNA correction | decontX is run on the private merged object | Not applied to canonical raw counts; omission cannot create missing barcodes | FALSE | May alter expression/annotation, not cell identities or the absent UC01 source barcodes |
| Cell QC | Paper thresholds: >=300 genes, >=500 UMIs, <20% mitochondrial reads | Applied sequentially and recorded per library and patient | TRUE | UC01 loss is not caused by threshold misapplication |
| Doublet removal | Paper names DoubletFinder; public parameters and calls are absent | scDblFinder-style closest reproducible implementation; scores/calls retained | FALSE | Likely contributes to the 249-cell global post-QC difference, but not the 605-barcode mismatch |
| Whole-bladder clustering | Public scripts use 3,000 HVGs, 100 PCs, k=30, final resolution 0.4 | Canonical reconstruction follows paper Methods: 2,000 HVGs, 75-PC graph, resolution 0.5 | FALSE | Paper and public code disagree; both settings are documented |
| Broad annotation | Private label transfer, external marker files, and manual Final_ID cluster overrides | Multi-marker module evidence for epithelial, immune, stromal, endothelial classes | FALSE | Exact authors' labels cannot be recreated from public materials |
| Epithelial selection | cells_epithelial is taken from an unreleased Epithelial object | Marker-supported broad epithelial cells; no access to authors' private membership list | FALSE | Exact 8,553-cell membership cannot be reconstructed from code alone |
| Tumour support by InferCNV | Paper reports InferCNV 1.4.0/HMM; public script only prepares input from private labels | Prepared a reproducible input only; no HMM result or keep/drop list is claimed | FALSE | Cannot adjudicate the missing UC01 sample |
| Tumour clustering | Paper Methods specify 2,000 HVGs, top 75 PCs, resolution 0.5 | Reclustered 7,862 cells using 2,000 HVGs, 75-PC graph, resolution 0.5 | FALSE | Settings were reproduced, but the input cohort and Cluster-13-like state were not |
| UC01 identity | Pre-2026 code uses 12923; commit d26f042 changes it to 12049; other scripts retain PureUC12923 | Audited filename, GEO title/description, Git history, and Figure5B source barcodes | FALSE | The deposited matrix assigned to UC01 is not the matrix behind the published UC01 cells |
| HV-versus-UC DE | Released Figure5A gives the final DE table; exact seed and complete input object are absent | Compared released result; ran fixed-seed and 100-seed public-DGE sensitivity analyses | FALSE | TM4SF1 direction is robust in available data; exact published genome-wide test is not reproducible |

The exact author object cannot be regenerated from the repository because the raw-matrix merge, `dge_mergedall_032123.rds`, reference label object, external marker lists, final epithelial object, manual keep/drop decisions, and InferCNV outputs are not public. DecontX could change gene counts and downstream label scores, but it cannot produce the missing published UC01 barcodes. Exact InferCNV support is likewise not reproducible: the public code prepares an input from private `Final_ID` labels but releases no gene-order file, HMM calls, result object, or epithelial membership decision.

## 5. Tumour reclustering and the reported Cluster 13 state

The 7,862 cells from the expected 12 patient IDs were reclustered with 2,000 HVGs, PCA, the top 75 PCs, k=30, and resolution 0.5. The closest signature-scoring candidate was cluster 3 (801 cells), but VAR07 contributed   86% of it; VAR05 contributed 111 cells and UC01 one cell. This patient-dominated composition does not match a convincing cross-patient recovery of the authors' Cluster 13 state. Candidate markers include KRT6A, KRT5, S100A2, FN1, VIM, and other basal/mesenchymal genes.

See `tumour_recluster_summary.tsv`, `tumour_recluster_composition.tsv`, `tumour_cluster13_candidate_markers.tsv`, `tumour_reclustering_cluster13_candidate.pdf`, and `tumour_signature_expression.pdf`.

## 6. TM4SF1 differential expression and downsampling

| analysis | n_HV_patients | n_UC_patients | n_cells_per_patient | avg_log2FC | p_val | p_val_adj | FDR_BH | pct.1 | pct.2 | HV_enriched_rank_by_p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| paper_released_source_DE_9HV_3UC_150_cells_per_patient_unknown_seed | 9 | 3 | 150 | 2.36002231418336 | 1.67525834324674e-44 | 4.44931863382902e-40 |  | 0.376 | 0.029 | 1 |
| closest_public_DGE_9HV_2UC_150_cells_per_patient_seed293189 | 9 | 2 | 150 | 4.3005528716779 | 2.05664493071056e-29 | 7.60835225667064e-25 | 2.37761008020958e-26 | 0.371851851851852 | 0.0366666666666667 | 1 |
| secondary_public_DGE_9HV_3UC_90_cells_per_patient_seed293189 | 9 | 3 | 90 | 2.41680600033156 | 1.81786148807034e-11 | 6.72499678896742e-07 | 2.26430868315401e-09 | 0.358024691358025 | 0.162962962962963 | 28 |

The released paper source table places TM4SF1 first among HV-enriched genes by p value (log2FC 2.36, p=1.68e-44, source adjusted p=4.45e-40). In the closest feasible public-DGE 150-cell comparison (9 HV + 2 UC), it is also rank 1 (log2FC  4.3, p=2.06e-29, BH FDR=2.38e-26). Across 100 seeds, its median rank is 1 (range 1-2; 79% rank 1 and 100% top 5).

A secondary 9 HV + 3 UC comparison is possible only by lowering every patient to 90 cells. At the fixed seed TM4SF1 ranks 28; over 100 seeds its median rank is 15 (range 5-32; 1% top 5). This is a sensitivity analysis, not a reproduction of the requested 150-cell design. Using the released per-cell TM4SF1 values for the exact 12 published patients, the HV-minus-UC signal remains positive in all 100 resamples, but genome-wide rank sensitivity cannot be calculated because the source workbook releases only per-cell TM4SF1 for this purpose, not the complete underlying count matrix.

| analysis | n_runs | median_rank | rank_min | rank_Q1 | rank_Q3 | rank_max | fraction_top1 | fraction_top5 | fraction_top10 | median_logFC |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| closest_public_DGE_9HV_2UC_150_cells_per_patient | 100 | 1 | 1 | 1 | 1 | 2 | 0.79 | 1 | 1 | 0.792360949730136 |
| secondary_public_DGE_9HV_3UC_90_cells_per_patient | 100 | 15 | 5 | 13 | 18 | 32 | 0 | 0.01 | 0.1 | 0.60592831170955 |

## 7. Validation and disposition

Canonical structural validation: **PASS** (31/31). Replication audit validation: **10/14 pass**. The four failures are the expected UC01/publication-target failures and are retained in `validation_checks_updated.tsv`.

No corrected reconstruction was emitted. The existing patient-level pseudobulk remains valid as a transparent public-DGE reconstruction, with UC01 explicitly flagged below 150 cells; it must not be represented as an exact reproduction of the paper's 9 HV + 3 UC / 8,553-cell cohort. Resolving the discrepancy requires the authors' original UC01/12923 count matrix or the exact private epithelial object and cell-membership list.

## 8. Reproducible outputs

- Audit driver: `scripts/run_replication_audit.R`
- UC01 and provenance audit: `scripts/07_replication_audit.R`
- Reclustering, DE, and sensitivity: `scripts/08_tumour_replication_de.R`
- Validation: `scripts/09_validate_replication.R`
- This report: `scripts/10_replication_report.R`
- Machine-readable method comparison: `results/replication_audit/authors_code_method_audit.tsv`
- Original-object checksums: `results/replication_audit/original_reconstruction/SHA256SUMS.txt`

All audit tables, plots, and the audit-only reclustered object are under `results/replication_audit/`. Canonical data remain under `data/processed/`; no posteriorHCA integration or tumour/toxicity threshold was attempted.
