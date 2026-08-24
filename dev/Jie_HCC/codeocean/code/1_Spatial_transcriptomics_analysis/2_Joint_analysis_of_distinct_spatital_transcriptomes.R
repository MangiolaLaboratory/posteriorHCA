#> sessionInfo()
#R version 4.2.1 (2022-06-23 ucrt)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 22621)

#Matrix products: default

#locale:
#[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8   
#[3] LC_MONETARY=Chinese (Simplified)_China.utf8 LC_NUMERIC=C                               
#[5] LC_TIME=Chinese (Simplified)_China.utf8    

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     

#loaded via a namespace (and not attached):
 #[1] SummarizedExperiment_1.28.0 KEGGREST_1.38.0             locfit_1.5-9.9             
 #[4] tidyselect_1.2.0            lattice_0.20-45             colorspace_2.1-0           
 #[7] vctrs_0.6.3                 generics_0.1.3              stats4_4.2.1               
#[10] utf8_1.2.3                  blob_1.2.4                  XML_3.99-0.16.1            
#[13] rlang_1.1.1                 pillar_1.9.0                glue_1.6.2                 
#[16] DBI_1.2.2                   BiocParallel_1.32.6         BiocGenerics_0.44.0        
#[19] bit64_4.0.5                 RColorBrewer_1.1-3          matrixStats_1.2.0          
#[22] GenomeInfoDbData_1.2.9      lifecycle_1.0.4             zlibbioc_1.44.0            
#[25] MatrixGenerics_1.10.0       Biostrings_2.66.0           munsell_0.5.0              
#[28] gtable_0.3.4                DESeq2_1.38.3               codetools_0.2-18           
#[31] memoise_2.0.1               Biobase_2.58.0              geneplotter_1.76.0         
#[34] IRanges_2.32.0              fastmap_1.1.1               GenomeInfoDb_1.34.9        
#[37] parallel_4.2.1              fansi_1.0.4                 AnnotationDbi_1.60.2       
#[40] Rcpp_1.0.11                 xtable_1.8-4                scales_1.3.0               
#[43] cachem_1.0.8                DelayedArray_0.24.0         S4Vectors_0.36.2           
#[46] annotate_1.76.0             XVector_0.38.0              bit_4.0.5                  
#[49] ggplot2_3.5.0               png_0.1-8                   dplyr_1.1.2                
#[52] GenomicRanges_1.50.2        grid_4.2.1                  cli_3.6.1                  
#[55] tools_4.2.1                 bitops_1.0-7                magrittr_2.0.3             
#[58] RCurl_1.98-1.14             tibble_3.2.1                RSQLite_2.3.5              
#[61] crayon_1.5.2                pkgconfig_2.0.3             Matrix_1.6-5               
#[64] httr_1.4.7                  rstudioapi_0.15.0           R6_2.5.1                   
#[67] compiler_4.2.1  



setwd("/code/1_Spatial_transcriptomics_analysis/")

##read in data
load("/data/Input_1_Spatial_transcriptomics_analysis/transcriptome_differential_analysis_input_3.RData")

####cross-compare all -- appropriate threshold 0.01/4
##use filtering low-mean count
dwgeneidx_1NKatStroma_agst_0NKatBorder <- which(niky_de$l2fc_1NKatStroma_agst_0NKatBorder<0 & niky_de$padj_1NKatStroma_agst_0NKatBorder<0.01/4)     #21
dwgeneidx_1NKatTumor_agst_0NKatBorder <- which(niky_de$l2fc_1NKatTumor_agst_0NKatBorder<0 & niky_de$padj_1NKatTumor_agst_0NKatBorder<0.01/4)    #96   
upgeneidx_1NKatStroma_agst_0NKatBorder <- which(niky_de$l2fc_1NKatStroma_agst_0NKatBorder>0 & niky_de$padj_1NKatStroma_agst_0NKatBorder<0.01/4)     #4
upgeneidx_1NKatTumor_agst_0NKatBorder <- which(niky_de$l2fc_1NKatTumor_agst_0NKatBorder>0 & niky_de$padj_1NKatTumor_agst_0NKatBorder<0.01/4)    #129
nsggeneidx_1NKatStroma_agst_0NKatBorder <- which(niky_de$padj_1NKatStroma_agst_0NKatBorder>=0.01/4)     #6980
nsggeneidx_1NKatTumor_agst_0NKatBorder <- which(niky_de$padj_1NKatTumor_agst_0NKatBorder>=0.01/4)    #15209

##without filtering out low-mean count, lose power, Nsig becomes much smaller. https://support.bioconductor.org/p/80194/
dwgeneidx_1NKatStroma_agst_0NKatBorder <- which(niky_de$l2fc_1NKatStroma_agst_0NKatBorder<0 & niky_de$BH_adjusted_pvalue_1NKatStroma_agst_0NKatBorder<0.01/4)     #8
dwgeneidx_1NKatTumor_agst_0NKatBorder <- which(niky_de$l2fc_1NKatTumor_agst_0NKatBorder<0 & niky_de$BH_adjusted_pvalue_1NKatTumor_agst_0NKatBorder<0.01/4)    #79     
upgeneidx_1NKatStroma_agst_0NKatBorder <- which(niky_de$l2fc_1NKatStroma_agst_0NKatBorder>0 & niky_de$BH_adjusted_pvalue_1NKatStroma_agst_0NKatBorder<0.01/4)     #7
upgeneidx_1NKatTumor_agst_0NKatBorder <- which(niky_de$l2fc_1NKatTumor_agst_0NKatBorder>0 & niky_de$BH_adjusted_pvalue_1NKatTumor_agst_0NKatBorder<0.01/4)    #113
nsggeneidx_1NKatStroma_agst_0NKatBorder <- which(niky_de$BH_adjusted_pvalue_1NKatStroma_agst_0NKatBorder>=0.01/4)     #18377   
nsggeneidx_1NKatTumor_agst_0NKatBorder <- which(niky_de$BH_adjusted_pvalue_1NKatTumor_agst_0NKatBorder>=0.01/4)    #18199 

##our final choice: l2fc & padj
dwgeneidx_1NKatStroma_agst_0NKatBorder <- which(niky_de$l2fc_1NKatStroma_agst_0NKatBorder<0 & niky_de$padj_1NKatStroma_agst_0NKatBorder<0.01/4)     #21
dwgeneidx_1NKatTumor_agst_0NKatBorder <- which(niky_de$l2fc_1NKatTumor_agst_0NKatBorder<0 & niky_de$padj_1NKatTumor_agst_0NKatBorder<0.01/4)    #96   
upgeneidx_1NKatStroma_agst_0NKatBorder <- which(niky_de$l2fc_1NKatStroma_agst_0NKatBorder>0 & niky_de$padj_1NKatStroma_agst_0NKatBorder<0.01/4)     #4
upgeneidx_1NKatTumor_agst_0NKatBorder <- which(niky_de$l2fc_1NKatTumor_agst_0NKatBorder>0 & niky_de$padj_1NKatTumor_agst_0NKatBorder<0.01/4)    #129
nsggeneidx_1NKatStroma_agst_0NKatBorder <- setdiff(which(!is.na(niky_de$pvalue_1NKatStroma_agst_0NKatBorder)),union(dwgeneidx_1NKatStroma_agst_0NKatBorder,upgeneidx_1NKatStroma_agst_0NKatBorder))     #18367
nsggeneidx_1NKatTumor_agst_0NKatBorder <- setdiff(which(!is.na(niky_de$pvalue_1NKatTumor_agst_0NKatBorder)),union(dwgeneidx_1NKatTumor_agst_0NKatBorder,upgeneidx_1NKatTumor_agst_0NKatBorder))     #18166
##96+129+18166=18391=18677-286
#21+4+18367=18392=18677-285

#total trajectories among all counts of AS, IF and TC
idx_concave <- intersect(upgeneidx_1NKatStroma_agst_0NKatBorder,upgeneidx_1NKatTumor_agst_0NKatBorder)   #0

idx_convex <- intersect(dwgeneidx_1NKatStroma_agst_0NKatBorder,dwgeneidx_1NKatTumor_agst_0NKatBorder)    #0
gene_df_convex$gene %in% GeneList$Gene_List   #logical(0)

idx_increase <- intersect(dwgeneidx_1NKatStroma_agst_0NKatBorder,upgeneidx_1NKatTumor_agst_0NKatBorder)   #0
gene_df_increase$gene %in% GeneList$Gene_List   #logical(0)

idx_decrease <- intersect(upgeneidx_1NKatStroma_agst_0NKatBorder,dwgeneidx_1NKatTumor_agst_0NKatBorder)   #0

idx_downflat <- intersect(upgeneidx_1NKatStroma_agst_0NKatBorder,nsggeneidx_1NKatTumor_agst_0NKatBorder)    #4
gene_df_downflat <- niky_de[idx_downflat,c("gene","gene_annotate","gene_annotate_2")]
length(which(gene_df_downflat$gene %in% GeneList$Gene_List))  #2

idx_upflat <- intersect(dwgeneidx_1NKatStroma_agst_0NKatBorder,nsggeneidx_1NKatTumor_agst_0NKatBorder)    #21
gene_df_upflat <- niky_de[idx_upflat,c("gene","gene_annotate","gene_annotate_2")]
length(which(gene_df_upflat$gene %in% GeneList$Gene_List))  #9

idx_flatup <- intersect(nsggeneidx_1NKatStroma_agst_0NKatBorder,upgeneidx_1NKatTumor_agst_0NKatBorder)    #129
length(which(niky_de[idx_flatup,"gene"] %in% GeneList$Gene_List))  #21

idx_flatdown <- intersect(nsggeneidx_1NKatStroma_agst_0NKatBorder,dwgeneidx_1NKatTumor_agst_0NKatBorder)    #96
length(which(niky_de[idx_flatdown,"gene"] %in% GeneList$Gene_List))  #36

idx_flatflat <- intersect(nsggeneidx_1NKatStroma_agst_0NKatBorder,nsggeneidx_1NKatTumor_agst_0NKatBorder)    #18140
length(which(niky_de[idx_flatflat,"gene"] %in% GeneList$Gene_List))  #2082

# 4+21+129+96+18140 = 18390 = 18677-287 (union of (285,286 genes))
# 2+9+21+36+2082 = 2150

display_df <- niky_de[idx_convex,]
display_df <- niky_de[idx_decrease,]
display_df <- niky_de[idx_downflat,]
display_df <- niky_de[idx_upflat,]
display_df <- niky_de[idx_flatup,]
display_df <- niky_de[idx_flatdown,]
display_df <- display_df[order(display_df$l2fc_1NKatStroma_agst_0NKatBorder),c("gene","gene_annotate","gene_annotate_2","in_GeneList",grep("1NKatStroma_agst_0NKatBorder",colnames(display_df),value=T),grep("1NKatTumor_agst_0NKatBorder",colnames(display_df),value=T))]


##use filtering low-mean count
dwgeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder <- which(niky_de$l2fc_1NKatStroma_agst_0Nonrelp_NKatBorder<0 & niky_de$padj_1NKatStroma_agst_0Nonrelp_NKatBorder<0.01/4)     #277
dwgeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder <- which(niky_de$l2fc_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder<0 & niky_de$padj_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder<0.01/4)    #215    
upgeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder <- which(niky_de$l2fc_1NKatStroma_agst_0Nonrelp_NKatBorder>0 & niky_de$padj_1NKatStroma_agst_0Nonrelp_NKatBorder<0.01/4)     #309
upgeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder <- which(niky_de$l2fc_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder>0 & niky_de$padj_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder<0.01/4)    #1175
nsggeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder <- which(niky_de$padj_1NKatStroma_agst_0Nonrelp_NKatBorder>=0.01/4)     #17794
nsggeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder <- which(niky_de$padj_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder>=0.01/4)    #10464  

##without filtering out low-mean count, lose power, Nsig becomes much smaller. https://support.bioconductor.org/p/80194/
dwgeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder <- which(niky_de$l2fc_1NKatStroma_agst_0Nonrelp_NKatBorder<0 & niky_de$BH_adjusted_pvalue_1NKatStroma_agst_0Nonrelp_NKatBorder<0.01/4)     #185
dwgeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder <- which(niky_de$l2fc_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder<0 & niky_de$BH_adjusted_pvalue_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder<0.01/4)    #194
upgeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder <- which(niky_de$l2fc_1NKatStroma_agst_0Nonrelp_NKatBorder>0 & niky_de$BH_adjusted_pvalue_1NKatStroma_agst_0Nonrelp_NKatBorder<0.01/4)     #172
upgeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder <- which(niky_de$l2fc_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder>0 & niky_de$BH_adjusted_pvalue_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder<0.01/4)    #1070
nsggeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder <- which(niky_de$BH_adjusted_pvalue_1NKatStroma_agst_0Nonrelp_NKatBorder>=0.01/4)     #18023
nsggeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder <- which(niky_de$BH_adjusted_pvalue_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder>=0.01/4)    #17123

##our final choice: l2fc & padj
dwgeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder <- which(niky_de$l2fc_1NKatStroma_agst_0Nonrelp_NKatBorder<0 & niky_de$padj_1NKatStroma_agst_0Nonrelp_NKatBorder<0.01/4)     #277
dwgeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder <- which(niky_de$l2fc_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder<0 & niky_de$padj_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder<0.01/4)    #215      
upgeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder <- which(niky_de$l2fc_1NKatStroma_agst_0Nonrelp_NKatBorder>0 & niky_de$padj_1NKatStroma_agst_0Nonrelp_NKatBorder<0.01/4)     #309   
upgeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder <- which(niky_de$l2fc_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder>0 & niky_de$padj_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder<0.01/4)    #1175   
nsggeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder <- setdiff(which(!is.na(niky_de$pvalue_1NKatStroma_agst_0Nonrelp_NKatBorder)),union(dwgeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder,upgeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder))     #17794
nsggeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder <- setdiff(which(!is.na(niky_de$pvalue_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder)),union(dwgeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder,upgeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder))     #16997
##277+309+17794=18380=18677-297
##215+1175+16997=18387=18677-290

#non-recurrent trajectories among AS, IF and TC.
idx_concave_Nonrelp <- intersect(upgeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder,upgeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder)   #1
gene_df_concave_Nonrelp <- niky_de[idx_concave_Nonrelp,c("gene","gene_annotate","gene_annotate_2")]
##gene gene_annotate                                                     gene_annotate_2
##8745 PSAT1          <NA> All Targets, Amino Acid Synthesis & Transport, Glutamine Metabolism
gene_df_concave_Nonrelp$gene %in% GeneList$Gene_List   #FALSE
length(which(niky_de[idx_concave_Nonrelp,"gene"] %in% GeneList$Gene_List))  #0

idx_convex_Nonrelp <- intersect(dwgeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder,dwgeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder)    #5
length(which(niky_de[idx_convex_Nonrelp,"gene"] %in% GeneList$Gene_List))  #3

idx_increase_Nonrelp <- intersect(dwgeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder,upgeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder)   #7
gene_df_increase_Nonrelp <- niky_de[idx_increase_Nonrelp,c("gene","gene_annotate","gene_annotate_2")]
##gene gene_annotate                                                     gene_annotate_2

##306     PYCR2       <NA> All Targets, Oxidative Stress, Amino Acid Synthesis & Transport, Glutamine Metabolism
##1615     H3C2       <NA> All Targets, Epigenetic Modification
##3618    IGF2R       <NA> All Targets, Neutrophil Degranulation, Lysosome, Endocytosis
##5102     CD46       <NA> All Targets, Complement System
##5499      MIF       <NA> All Targets, Amino Acid Synthesis & Transport, Neutrophil Degranulation, Other Interleukin Signaling
##6986    H3-3A       <NA> All Targets, Notch Signaling, Epigenetic Modification, Wnt Signaling, Estrogen Signaling, Senescence, Cell Cycle, Differentiation
##7741 MAPKAPK2       <NA> All Targets, TLR Signaling, Lipid Metabolism, IL-17 Signaling, MAPK Signaling, Senescence, VEGF Signaling
  
gene_df_increase_Nonrelp$gene %in% GeneList$Gene_List   #FALSE FALSE  TRUE  TRUE  TRUE FALSE  TRUE
length(which(niky_de[idx_increase_Nonrelp,"gene"] %in% GeneList$Gene_List))  #4

idx_decrease_Nonrelp <- intersect(upgeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder,dwgeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder)   #0
length(which(niky_de[idx_decrease_Nonrelp,"gene"] %in% GeneList$Gene_List))  #0

idx_downflat_Nonrelp <- intersect(upgeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder,nsggeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder)    #303 
length(which(niky_de[idx_downflat_Nonrelp,"gene"] %in% GeneList$Gene_List))  #143

idx_upflat_Nonrelp <- intersect(dwgeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder,nsggeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder)    #265
length(which(niky_de[idx_upflat_Nonrelp,"gene"] %in% GeneList$Gene_List))  #116

idx_flatup_Nonrelp <- intersect(nsggeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder,upgeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder)    #1163
length(which(niky_de[idx_flatup_Nonrelp,"gene"] %in% GeneList$Gene_List))  #110

idx_flatdown_Nonrelp <- intersect(nsggeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder,dwgeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder)    #209
length(which(niky_de[idx_flatdown_Nonrelp,"gene"] %in% GeneList$Gene_List))  #47

idx_flatflat_Nonrelp <- intersect(nsggeneidx_1NKatStroma_agst_0Nonrelp_NKatBorder,nsggeneidx_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder)    #16421
length(which(niky_de[idx_flatflat_Nonrelp,"gene"] %in% GeneList$Gene_List))  #1721

# 1+5+7+0+303+265+1163+209+16421 = 18374 = 18677-303 (union of (315,295 genes))
# 0+3+4+0+143+116+110+47+1721 = 2144

display_df <- niky_de[idx_concave_Nonrelp,]
display_df <- niky_de[idx_convex_Nonrelp,]
display_df <- niky_de[idx_increase_Nonrelp,]
display_df <- niky_de[idx_decrease_Nonrelp,]
display_df <- niky_de[idx_downflat_Nonrelp,]
display_df <- niky_de[idx_upflat_Nonrelp,]
display_df <- niky_de[idx_flatup_Nonrelp,]
display_df <- niky_de[idx_flatdown_Nonrelp,]
display_df <- display_df[order(display_df$l2fc_1NKatStroma_agst_0Nonrelp_NKatBorder),c("gene","gene_annotate","gene_annotate_2","in_GeneList",grep("1NKatStroma_agst_0Nonrelp_NKatBorder",colnames(display_df),value=T),grep("1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder",colnames(display_df),value=T))]

####cross-compare recurrent patients 
##use filtering low-mean count
dwgeneidx_1NKatStroma_agst_0Relp_NKatBorder <- which(niky_de$l2fc_1NKatStroma_agst_0Relp_NKatBorder<0 & niky_de$padj_1NKatStroma_agst_0Relp_NKatBorder<0.01/4)     #430
dwgeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder <- which(niky_de$l2fc_1Relp_NKatTumor_agst_0Relp_NKatBorder<0 & niky_de$padj_1Relp_NKatTumor_agst_0Relp_NKatBorder<0.01/4)    #0        
upgeneidx_1NKatStroma_agst_0Relp_NKatBorder <- which(niky_de$l2fc_1NKatStroma_agst_0Relp_NKatBorder>0 & niky_de$padj_1NKatStroma_agst_0Relp_NKatBorder<0.01/4)     #377
upgeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder <- which(niky_de$l2fc_1Relp_NKatTumor_agst_0Relp_NKatBorder>0 & niky_de$padj_1Relp_NKatTumor_agst_0Relp_NKatBorder<0.01/4)    #0 
nsggeneidx_1NKatStroma_agst_0Relp_NKatBorder <- which(niky_de$padj_1NKatStroma_agst_0Relp_NKatBorder>=0.01/4)     #11443   
nsggeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder <- which(niky_de$padj_1Relp_NKatTumor_agst_0Relp_NKatBorder>=0.01/4)    #18365

##without filtering out low-mean count, lose power, Nsig becomes much smaller. https://support.bioconductor.org/p/80194/
dwgeneidx_1NKatStroma_agst_0Relp_NKatBorder <- which(niky_de$l2fc_1NKatStroma_agst_0Relp_NKatBorder<0 & niky_de$BH_adjusted_pvalue_1NKatStroma_agst_0Relp_NKatBorder<0.01/4)     #405
dwgeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder <- which(niky_de$l2fc_1Relp_NKatTumor_agst_0Relp_NKatBorder<0 & niky_de$BH_adjusted_pvalue_1Relp_NKatTumor_agst_0Relp_NKatBorder<0.01/4)    #0        
upgeneidx_1NKatStroma_agst_0Relp_NKatBorder <- which(niky_de$l2fc_1NKatStroma_agst_0Relp_NKatBorder>0 & niky_de$BH_adjusted_pvalue_1NKatStroma_agst_0Relp_NKatBorder<0.01/4)     #350
upgeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder <- which(niky_de$l2fc_1Relp_NKatTumor_agst_0Relp_NKatBorder>0 & niky_de$BH_adjusted_pvalue_1Relp_NKatTumor_agst_0Relp_NKatBorder<0.01/4)    #0
nsggeneidx_1NKatStroma_agst_0Relp_NKatBorder <- which(niky_de$BH_adjusted_pvalue_1NKatStroma_agst_0Relp_NKatBorder>=0.01/4)     #17611 
nsggeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder <- which(niky_de$BH_adjusted_pvalue_1Relp_NKatTumor_agst_0Relp_NKatBorder>=0.01/4)    #18365

##our final choice: l2fc & padj
dwgeneidx_1NKatStroma_agst_0Relp_NKatBorder <- which(niky_de$l2fc_1NKatStroma_agst_0Relp_NKatBorder<0 & niky_de$padj_1NKatStroma_agst_0Relp_NKatBorder<0.01/4)     #430
dwgeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder <- which(niky_de$l2fc_1Relp_NKatTumor_agst_0Relp_NKatBorder<0 & niky_de$padj_1Relp_NKatTumor_agst_0Relp_NKatBorder<0.01/4)    #0        
upgeneidx_1NKatStroma_agst_0Relp_NKatBorder <- which(niky_de$l2fc_1NKatStroma_agst_0Relp_NKatBorder>0 & niky_de$padj_1NKatStroma_agst_0Relp_NKatBorder<0.01/4)     #377
upgeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder <- which(niky_de$l2fc_1Relp_NKatTumor_agst_0Relp_NKatBorder>0 & niky_de$padj_1Relp_NKatTumor_agst_0Relp_NKatBorder<0.01/4)    #0 
nsggeneidx_1NKatStroma_agst_0Relp_NKatBorder <- setdiff(which(!is.na(niky_de$pvalue_1NKatStroma_agst_0Relp_NKatBorder)),union(dwgeneidx_1NKatStroma_agst_0Relp_NKatBorder,upgeneidx_1NKatStroma_agst_0Relp_NKatBorder))      #17559  
nsggeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder <- setdiff(which(!is.na(niky_de$pvalue_1Relp_NKatTumor_agst_0Relp_NKatBorder)),union(dwgeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder,upgeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder))    #18365
length(dwgeneidx_1NKatStroma_agst_0Relp_NKatBorder) #430
length(dwgeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder) #0
length(upgeneidx_1NKatStroma_agst_0Relp_NKatBorder) #377
length(upgeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder) #0
length(nsggeneidx_1NKatStroma_agst_0Relp_NKatBorder) #17559
length(nsggeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder) #18365

##430+377+17559=18366=18677-311
##0+0+18365=18365=18677-312

#recurrent trajectories among AS, IF and TC.
idx_concave_Relp <- intersect(upgeneidx_1NKatStroma_agst_0Relp_NKatBorder,upgeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder)   #0

idx_convex_Relp <- intersect(dwgeneidx_1NKatStroma_agst_0Relp_NKatBorder,dwgeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder)    #0

idx_increase_Relp <- intersect(dwgeneidx_1NKatStroma_agst_0Relp_NKatBorder,upgeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder)   #0

idx_decrease_Relp <- intersect(upgeneidx_1NKatStroma_agst_0Relp_NKatBorder,dwgeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder)   #0

idx_downflat_Relp <- intersect(upgeneidx_1NKatStroma_agst_0Relp_NKatBorder,nsggeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder)    #377
length(which(niky_de[idx_downflat_Relp,"gene"] %in% GeneList$Gene_List))  #78

idx_upflat_Relp <- intersect(dwgeneidx_1NKatStroma_agst_0Relp_NKatBorder,nsggeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder)    #430
length(which(niky_de[idx_upflat_Relp,"gene"] %in% GeneList$Gene_List))  #43

idx_flatup_Relp <- intersect(nsggeneidx_1NKatStroma_agst_0Relp_NKatBorder,upgeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder)    #0
length(which(niky_de[idx_flatup_Relp,"gene"] %in% GeneList$Gene_List))  #0

idx_flatdown_Relp <- intersect(nsggeneidx_1NKatStroma_agst_0Relp_NKatBorder,dwgeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder)    #0

idx_flatflat_Relp <- intersect(nsggeneidx_1NKatStroma_agst_0Relp_NKatBorder,nsggeneidx_1Relp_NKatTumor_agst_0Relp_NKatBorder)    #17551
length(which(niky_de[idx_flatflat_Relp,"gene"] %in% GeneList$Gene_List))  #2025

# 0+0+0+0+377+430+0+0+17551 = 18358 = 18677-319 (union of (289,290 genes))
# 78+43+0+2025 = 2146

display_df <- niky_de[idx_downflat_Relp,]
display_df <- niky_de[idx_upflat_Relp,]
display_df <- display_df[order(display_df$l2fc_1NKatStroma_agst_0Relp_NKatBorder),c("gene","gene_annotate","gene_annotate_2","in_GeneList",grep("1NKatStroma_agst_0Relp_NKatBorder",colnames(display_df),value=T),grep("1Relp_NKatTumor_agst_0Relp_NKatBorder",colnames(display_df),value=T))]

##two-way cross-check，need to save the gene names, use BH pval < 0.01/4
gene_inter <- intersect(niky_de[idx_flatflat_Nonrelp,"gene"],niky_de[idx_concave_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_flatflat_Nonrelp,"gene"],niky_de[idx_convex_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_flatflat_Nonrelp,"gene"],niky_de[idx_increase_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_flatflat_Nonrelp,"gene"],niky_de[idx_decrease_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_flatflat_Nonrelp,"gene"],niky_de[idx_downflat_Relp,"gene"])  #245 --> 59  
gene_name <- data.frame(gene=gene_inter,nonrelp="flatflat",relp="downflat",stringsAsFactors=F)
gene_inter <- intersect(niky_de[idx_flatflat_Nonrelp,"gene"],niky_de[idx_upflat_Relp,"gene"])  #313 --> 26 
gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="flatflat",relp="upflat",stringsAsFactors=F))
gene_inter <- intersect(niky_de[idx_flatflat_Nonrelp,"gene"],niky_de[idx_flatup_Relp,"gene"])  #0 --> 0
#gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="flatflat",relp="flatup",stringsAsFactors=F))
#gene_inter <- intersect(niky_de[idx_flatflat_Nonrelp,"gene"],niky_de[idx_flatdown_Relp,"gene"])  
#gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="flatflat",relp="flatdown",stringsAsFactors=F))
gene_inter <- intersect(niky_de[idx_flatflat_Nonrelp,"gene"],niky_de[idx_flatflat_Relp,"gene"])  #15830 --> 1633
#245+313+15830 = 16388
length(intersect(niky_de[idx_flatflat_Nonrelp,"gene"],c(niky_de[idx_flatflat_Relp,"gene"],niky_de[idx_flatdown_Relp,"gene"],niky_de[idx_flatup_Relp,"gene"],niky_de[idx_upflat_Relp,"gene"],niky_de[idx_downflat_Relp,"gene"],niky_de[idx_decrease_Relp,"gene"],niky_de[idx_increase_Relp,"gene"],niky_de[idx_convex_Relp,"gene"],niky_de[idx_concave_Relp,"gene"])))  #16388

#59+26+1+1633 = 1718

gene_inter <- intersect(niky_de[idx_flatdown_Nonrelp,"gene"],niky_de[idx_downflat_Relp,"gene"])  #113 --> 17
gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="flatdown",relp="downflat",stringsAsFactors=F))
gene_inter <- intersect(niky_de[idx_flatdown_Nonrelp,"gene"],niky_de[idx_upflat_Relp,"gene"])  #3 --> 0
gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="flatdown",relp="upflat",stringsAsFactors=F))
gene_inter <- intersect(niky_de[idx_flatdown_Nonrelp,"gene"],niky_de[idx_flatup_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_flatdown_Nonrelp,"gene"],niky_de[idx_flatdown_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_flatdown_Nonrelp,"gene"],niky_de[idx_flatflat_Relp,"gene"])  #91 --> 28
gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="flatdown",relp="flatflat",stringsAsFactors=F))
#113+3+91 = 207
#17+28 = 45

gene_inter <- intersect(niky_de[idx_flatup_Nonrelp,"gene"],niky_de[idx_downflat_Relp,"gene"])  #19 --> 2
gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="flatup",relp="downflat",stringsAsFactors=F))
gene_inter <- intersect(niky_de[idx_flatup_Nonrelp,"gene"],niky_de[idx_upflat_Relp,"gene"])  #106 --> 15
gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="flatup",relp="upflat",stringsAsFactors=F))
gene_inter <- intersect(niky_de[idx_flatup_Nonrelp,"gene"],niky_de[idx_flatup_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_flatup_Nonrelp,"gene"],niky_de[idx_flatdown_Relp,"gene"])  #0 --> 0
#gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="flatup",relp="flatdown",stringsAsFactors=F))
gene_inter <- intersect(niky_de[idx_flatup_Nonrelp,"gene"],niky_de[idx_flatflat_Relp,"gene"])  #1038 --> 93
gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="flatup",relp="flatflat",stringsAsFactors=F))
#19+106+1+1038 = 1164
#2+15+93 = 110

#gene_inter <- intersect(niky_de[idx_upflat_Nonrelp,"gene"],niky_de[idx_downflat_Relp,"gene"])  
#gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="upflat",relp="downflat",stringsAsFactors=F))
gene_inter <- intersect(niky_de[idx_upflat_Nonrelp,"gene"],niky_de[idx_upflat_Relp,"gene"])  #3 --> 1
gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="upflat",relp="upflat",stringsAsFactors=F))
gene_inter <- intersect(niky_de[idx_upflat_Nonrelp,"gene"],niky_de[idx_flatup_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_upflat_Nonrelp,"gene"],niky_de[idx_flatdown_Relp,"gene"])  #0 --> 0
#gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="upflat",relp="flatdown",stringsAsFactors=F))
gene_inter <- intersect(niky_de[idx_upflat_Nonrelp,"gene"],niky_de[idx_flatflat_Relp,"gene"])  #262 --> 115
gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="upflat",relp="flatflat",stringsAsFactors=F))
#3+262 = 265
#1+115 = 116

gene_inter <- intersect(niky_de[idx_downflat_Nonrelp,"gene"],niky_de[idx_downflat_Relp,"gene"])  #0 --> 0
#gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="downflat",relp="downflat",stringsAsFactors=F))
gene_inter <- intersect(niky_de[idx_downflat_Nonrelp,"gene"],niky_de[idx_upflat_Relp,"gene"])
#gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="downflat",relp="upflat",stringsAsFactors=F))
gene_inter <- intersect(niky_de[idx_downflat_Nonrelp,"gene"],niky_de[idx_flatup_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_downflat_Nonrelp,"gene"],niky_de[idx_flatdown_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_downflat_Nonrelp,"gene"],niky_de[idx_flatflat_Relp,"gene"])  #303 --> 143
gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="downflat",relp="flatflat",stringsAsFactors=F))
#0+0+0+303=303
#0+0+0+143=143

gene_inter <- intersect(niky_de[idx_decrease_Nonrelp,"gene"],niky_de[idx_downflat_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_decrease_Nonrelp,"gene"],niky_de[idx_upflat_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_decrease_Nonrelp,"gene"],niky_de[idx_flatup_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_decrease_Nonrelp,"gene"],niky_de[idx_flatdown_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_decrease_Nonrelp,"gene"],niky_de[idx_flatflat_Relp,"gene"])  #0 --> 0
#gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="decrease",relp="flatflat",stringsAsFactors=F))
#0+0+3 = 3
#0+0+1 = 1

gene_inter <- intersect(niky_de[idx_increase_Nonrelp,"gene"],niky_de[idx_downflat_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_increase_Nonrelp,"gene"],niky_de[idx_upflat_Relp,"gene"])  #0 --> 0
#gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="increase",relp="upflat",stringsAsFactors=F))
gene_inter <- intersect(niky_de[idx_increase_Nonrelp,"gene"],niky_de[idx_flatup_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_increase_Nonrelp,"gene"],niky_de[idx_flatdown_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_increase_Nonrelp,"gene"],niky_de[idx_flatflat_Relp,"gene"])  #7 --> 4
gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="increase",relp="flatflat",stringsAsFactors=F))
# 0+0+0+7 = 7
# 0+0+0+4 = 4

gene_inter <- intersect(niky_de[idx_convex_Nonrelp,"gene"],niky_de[idx_downflat_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_convex_Nonrelp,"gene"],niky_de[idx_upflat_Relp,"gene"])  #0 --> 0
#gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="convex",relp="upflat",stringsAsFactors=F))
gene_inter <- intersect(niky_de[idx_convex_Nonrelp,"gene"],niky_de[idx_flatup_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_convex_Nonrelp,"gene"],niky_de[idx_flatdown_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_convex_Nonrelp,"gene"],niky_de[idx_flatflat_Relp,"gene"])  #5 --> 3
gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="convex",relp="flatflat",stringsAsFactors=F))
#0+0+0+0+5 = 5
#0+0+0+0+3 = 3

gene_inter <- intersect(niky_de[idx_concave_Nonrelp,"gene"],niky_de[idx_downflat_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_concave_Nonrelp,"gene"],niky_de[idx_upflat_Relp,"gene"])  #0 --> 0
#gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="concave",relp="upflat",stringsAsFactors=F))
gene_inter <- intersect(niky_de[idx_concave_Nonrelp,"gene"],niky_de[idx_flatup_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_concave_Nonrelp,"gene"],niky_de[idx_flatdown_Relp,"gene"])  #0 --> 0
gene_inter <- intersect(niky_de[idx_concave_Nonrelp,"gene"],niky_de[idx_flatflat_Relp,"gene"])  #1 --> 0
gene_name <- rbind(gene_name,data.frame(gene=gene_inter,nonrelp="concave",relp="flatflat",stringsAsFactors=F))
#0+0+0+0+1 = 1
#0

length(gene_inter)
length(which(gene_inter %in% GeneList$Gene_List))

gene_name$immune_gene <- "F"
gene_name$immune_gene[which(gene_name$gene %in% GeneList$Gene_List)] <- "T"
dim(gene_name) #2509    4

#16388+207+1164+265+303+3+7+5+1=18342
#1718+45+110+116+143+1+4+3+1=2141

#write.csv(gene_name,"gene_name_0.01div4_final.csv")

