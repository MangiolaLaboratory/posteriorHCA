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
#[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
 #[1] pheatmap_1.0.12             DESeq2_1.38.3               SummarizedExperiment_1.28.0
 #[4] Biobase_2.58.0              MatrixGenerics_1.10.0       matrixStats_1.2.0          
 #[7] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9         IRanges_2.32.0             
#[10] S4Vectors_0.36.2            BiocGenerics_0.44.0        

#loaded via a namespace (and not attached):
 #[1] KEGGREST_1.38.0        locfit_1.5-9.9         tidyselect_1.2.0      
 #[4] lattice_0.20-45        colorspace_2.1-0       vctrs_0.6.3           
 #[7] generics_0.1.3         utf8_1.2.3             blob_1.2.4            
#[10] XML_3.99-0.16.1        rlang_1.1.1            pillar_1.9.0          
#[13] glue_1.6.2             DBI_1.2.2              BiocParallel_1.32.6   
#[16] bit64_4.0.5            RColorBrewer_1.1-3     GenomeInfoDbData_1.2.9
#[19] lifecycle_1.0.4        zlibbioc_1.44.0        Biostrings_2.66.0     
#[22] munsell_0.5.0          gtable_0.3.4           codetools_0.2-18      
#[25] memoise_2.0.1          geneplotter_1.76.0     fastmap_1.1.1         
#[28] parallel_4.2.1         fansi_1.0.4            AnnotationDbi_1.60.2  
#[31] Rcpp_1.0.11            xtable_1.8-4           scales_1.3.0          
#[34] cachem_1.0.8           DelayedArray_0.24.0    annotate_1.76.0       
#[37] XVector_0.38.0         bit_4.0.5              ggplot2_3.5.0         
#[40] png_0.1-8              dplyr_1.1.2            grid_4.2.1            
#[43] cli_3.6.1              tools_4.2.1            bitops_1.0-7          
#[46] magrittr_2.0.3         RCurl_1.98-1.14        tibble_3.2.1          
#[49] RSQLite_2.3.5          crayon_1.5.2           pkgconfig_2.0.3       
#[52] Matrix_1.6-5           httr_1.4.7             rstudioapi_0.15.0     
#[55] R6_2.5.1               compiler_4.2.1    


setwd("/code/1_Spatial_transcriptomics_analysis/")

##1. read in data
load("/data/Input_1_Spatial_transcriptomics_analysis/transcriptome_differential_analysis_input_1.RData")

########################Analysis 2. DE############################
####DEseq2 analysis ##to find out the better method of statistic analysis for HCC patients data
library("DESeq2"); library("pheatmap")
packageVersion("DESeq2")  #‘1.38.3’  ##!for newer version

for (i in seq(1,10)){  #pair-to-pair comparisons of nonrelapsed and  relapsed patients' transcriptome data at three locations: border, tumor, stroma
  cat(i,"\n")   #let differern comparison pairs put out in different lines
  
  if (i==1){
    top1_i <-"0NKatBorder"; top2_i <- "1NKatStroma"  #IF as control(0), AS as case(1), compare the differences in gene expression between the two sites
    study1_indx <- NKatBorder_indx
    study1_indx <- paste(rep("...",length(study1_indx)),study1_indx,sep="")  #list out columns of NKatBorder_indx in de_data_all
    study2_indx <- NKatStroma_indx
    study2_indx <- paste(rep("...",length(study2_indx)),study2_indx,sep="")  #list out columns of NKatStroma_indx in de_data_all
    
    niky_de[[paste("ncase_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data framev de_data_all name as ncase_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("nctrl_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as nctrl_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("baseMean_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as baseMean_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("l2fc_",top2_i,"_agst_",top1_i,sep="")]] <- NA   #Insert a new blank column in data frame de_data_all name as l2fc_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("lfcSE_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as pvalue_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("pvalue_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as nctrl_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("padj_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as padj_1NKatStroma_agst_0NKatBorder_indx
  } else if (i==2){
    top1_i <-"0NKatTumor"; top2_i <- "1NKatStroma"  #tumor as control(0), stroma as case(1), compare the differences in gene expression between the two sites
    study1_indx <- NKatTumor_indx
    study1_indx <- paste(rep("...",length(study1_indx)),study1_indx,sep="")  #list out columns of NKatTumor_indx in de_data_all
    study2_indx <- NKatStroma_indx
    study2_indx <- paste(rep("...",length(study2_indx)),study2_indx,sep="")  #list out columns of NKatStroma_indx in de_data_all
    
    niky_de[[paste("ncase_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as ncase_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("nctrl_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as nctrl_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("baseMean_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as baseMean_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("l2fc_",top2_i,"_agst_",top1_i,sep="")]] <- NA   #Insert a new blank column in data frame de_data_all name as l2fc_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("lfcSE_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as pvalue_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("pvalue_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as nctrl_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("padj_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as padj_1NKatStroma_agst_0NKatBorder_indx
  } else if (i==3){
    top1_i <-"0NKatBorder"; top2_i <- "1NKatTumor"  #IF as control(0), TC as case(1), compare the differences in gene expression between the two sites
    study1_indx <- NKatBorder_indx
    study1_indx <- paste(rep("...",length(study1_indx)),study1_indx,sep="")  #list out columns of NKatBorder_indx in de_data_all
    study2_indx <- NKatTumor_indx
    study2_indx <- paste(rep("...",length(study2_indx)),study2_indx,sep="")  #list out columns of NKatTumor_indx in de_data_all
    
    niky_de[[paste("ncase_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as ncase_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("nctrl_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as nctrl_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("baseMean_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as baseMean_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("l2fc_",top2_i,"_agst_",top1_i,sep="")]] <- NA   #Insert a new blank column in data frame de_data_all name as l2fc_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("lfcSE_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as pvalue_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("pvalue_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as nctrl_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("padj_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as padj_1NKatStroma_agst_0NKatBorder_indx
  } else if (i==4){
    top1_i <-"0Relp_NKatBorder"; top2_i <- "1NKatStroma"  #recurrent IF as control(0), AS as case(1), compare the differences in gene expression between the two sites
    study1_indx <- Relp_NKatBorder_indx
    study1_indx <- paste(rep("...",length(study1_indx)),study1_indx,sep="")  #list out columns of NKatBorder_indx in de_data_all
    study2_indx <- NKatStroma_indx
    study2_indx <- paste(rep("...",length(study2_indx)),study2_indx,sep="")  #list out columns of NKatStroma_indx in de_data_all
    
    niky_de[[paste("ncase_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as ncase_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("nctrl_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as nctrl_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("baseMean_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as baseMean_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("l2fc_",top2_i,"_agst_",top1_i,sep="")]] <- NA   #Insert a new blank column in data frame de_data_all name as l2fc_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("lfcSE_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as pvalue_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("pvalue_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as nctrl_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("padj_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as padj_1NKatStroma_agst_0NKatBorder_indx
  } else if (i==5){
    top1_i <-"0Relp_NKatBorder"; top2_i <- "1Relp_NKatTumor"  #recurrent IF as control(0), recurrent TC as case(1), compare the differences in gene expression between the two sites
    study1_indx <- Relp_NKatBorder_indx
    study1_indx <- paste(rep("...",length(study1_indx)),study1_indx,sep="")  #list out columns of NKatBorder_indx in de_data_all
    study2_indx <- Relp_NKatTumor_indx
    study2_indx <- paste(rep("...",length(study2_indx)),study2_indx,sep="")  #list out columns of NKatStroma_indx in de_data_all
    
    niky_de[[paste("ncase_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as ncase_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("nctrl_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as nctrl_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("baseMean_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as baseMean_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("l2fc_",top2_i,"_agst_",top1_i,sep="")]] <- NA   #Insert a new blank column in data frame de_data_all name as l2fc_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("lfcSE_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as pvalue_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("pvalue_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as nctrl_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("padj_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as padj_1NKatStroma_agst_0NKatBorder_indx
  } else if (i==6){
    top1_i <-"0Nonrelp_NKatBorder"; top2_i <- "1NKatStroma"  #non-recurrent IF as control(0), AS as case(1), compare the differences in gene expression between the two sites
    study1_indx <- Nonrelp_NKatBorder_indx
    study1_indx <- paste(rep("...",length(study1_indx)),study1_indx,sep="")  #list out columns of NKatBorder_indx in de_data_all
    study2_indx <- NKatStroma_indx
    study2_indx <- paste(rep("...",length(study2_indx)),study2_indx,sep="")  #list out columns of NKatStroma_indx in de_data_all
    
    niky_de[[paste("ncase_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as ncase_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("nctrl_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as nctrl_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("baseMean_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as baseMean_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("l2fc_",top2_i,"_agst_",top1_i,sep="")]] <- NA   #Insert a new blank column in data frame de_data_all name as l2fc_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("lfcSE_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as pvalue_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("pvalue_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as nctrl_1NKatStroma_agst_0NKatBorder_indx
    niky_de[[paste("padj_",top2_i,"_agst_",top1_i,sep="")]] <- NA  #Insert a new blank column in data frame de_data_all name as padj_1NKatStroma_agst_0NKatBorder_indx
  } else if (i==7){
    top1_i <-"0Nonrelp_NKatBorder"; top2_i <- "1Nonrelp_NKatTumor"  
    study1_indx <- Nonrelp_NKatBorder_indx
    study1_indx <- paste(rep("...",length(study1_indx)),study1_indx,sep="")
    study2_indx <- Nonrelp_NKatTumor_indx
    study2_indx <- paste(rep("...",length(study2_indx)),study2_indx,sep="")
    
    niky_de[[paste("ncase_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("nctrl_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("baseMean_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("l2fc_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("lfcSE_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("pvalue_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("padj_",top2_i,"_agst_",top1_i,sep="")]] <- NA
  } else if (i==8){
    top1_i <-"0Relp_NKatTumor"; top2_i <- "1Nonrelp_NKatTumor"
    study1_indx <- Relp_NKatTumor_indx
    study1_indx <- paste(rep("...",length(study1_indx)),study1_indx,sep="")
    study2_indx <- Nonrelp_NKatTumor_indx
    study2_indx <- paste(rep("...",length(study2_indx)),study2_indx,sep="")
    
    niky_de[[paste("ncase_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("nctrl_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("baseMean_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("l2fc_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("lfcSE_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("pvalue_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("padj_",top2_i,"_agst_",top1_i,sep="")]] <- NA
  } else if (i==9){
    top1_i <-"0Relp_NKatBorder"; top2_i <- "1Nonrelp_NKatBorder"
    study1_indx <- Relp_NKatBorder_indx
    study1_indx <- paste(rep("...",length(study1_indx)),study1_indx,sep="")
    study2_indx <- Nonrelp_NKatBorder_indx
    study2_indx <- paste(rep("...",length(study2_indx)),study2_indx,sep="")
    
    niky_de[[paste("ncase_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("nctrl_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("baseMean_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("l2fc_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("lfcSE_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("pvalue_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("padj_",top2_i,"_agst_",top1_i,sep="")]] <- NA
  } else {
    top1_i <-"0Relp"; top2_i <- "1Nonrelp"
    study1_indx <- Relp_indx
    study1_indx <- paste(rep("...",length(study1_indx)),study1_indx,sep="")
    study2_indx <- Nonrelp_indx
    study2_indx <- paste(rep("...",length(study2_indx)),study2_indx,sep="")
    
    niky_de[[paste("ncase_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("nctrl_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("baseMean_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("l2fc_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("lfcSE_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("pvalue_",top2_i,"_agst_",top1_i,sep="")]] <- NA
    niky_de[[paste("padj_",top2_i,"_agst_",top1_i,sep="")]] <- NA
  }
  
  
  cat("case:",top2_i,length(study2_indx)," VS ","ctrl:",top1_i,length(study1_indx), "\n")  #output every case(1) vs control(0) groups of above for loop in different lines
  topic_2_sort <- as.factor(c(rep(top1_i,length(study1_indx)), rep(top2_i,length(study2_indx))))  ##make sure other--> 0ther, so that it is at frond alphabeta
  de_matrix_2_sort <- data_gene_df_ed[,c(study1_indx,study2_indx)]  #extract the columns of study1_indx and study2_indx in data_gene_df_ed as de_matrix_2_sort
  
  age_2_sort <- as.numeric(t(data_meta_df_ed["age",c(study1_indx,study2_indx)]))  #extract the columns of study1_indx and study2_indx and rows of age in data_meta_df_ed and transfer it in to numeric data and transform it as age_2_sort
  
  gender_num <- t(data_meta_df_ed["gender",c(study1_indx,study2_indx)])  #extract the columns of study1_indx and study2_indx and rows of gender in data_meta_df_ed and transform it as gender_num
  gender_num <- sub("Male","1",gender_num)  #lable the Male in gander_num as 1
  gender_num <- sub("Female","2",gender_num)  #lable the Female in gander_num as 2
  sex_2_sort <- as.factor(gender_num[,1])  #transform the class of the first column in gender_num into factor 
  
  de_matrix_2_sort_int <- matrix(nrow = dim(de_matrix_2_sort)[1], ncol = dim(de_matrix_2_sort)[2])  #creat a matrix that has the same dimension as de_matrix_2_sort
  for (i in 1: dim(de_matrix_2_sort)[2]){
    de_matrix_2_sort_int[,i] <- round(de_matrix_2_sort[,i])   #take the integer value of columns i of de_matrix_2_sort
  }
  
  de_df_2_sort_int_raw <- as.data.frame(de_matrix_2_sort_int)
  row.names(de_df_2_sort_int_raw) <- row.names(data_gene_df_ed)  #make data frames data_gene_df_ed and de_df_2_sort_int_raw have the same row name
  colnames(de_df_2_sort_int_raw) <- c(study1_indx,study2_indx)   #use study1_indx and study2_indx as column names for the data, strat from start from 1 (after minus 1)
  
  ##1788 genes on large collection of sample
  de_df_2_sort_int <- na.omit(de_df_2_sort_int_raw)  #remove na rows or in other words genes (sample number untouched)
  de_df_2_sort_int <- de_df_2_sort_int[rowSums(de_df_2_sort_int)>=5,]  #exrtact column at least 5 rows in de_df_2_sort_int, filter low abundance data
  
  age_2_sort_cut <- cut(age_2_sort,breaks=c(min(age_2_sort)-1,median(age_2_sort),max(age_2_sort)+1))  #divide age_2_sort into three segments with three cricital values
  col_data <- DataFrame(topic_2_sort=topic_2_sort,age_2_sort=age_2_sort_cut,sex_2_sort=sex_2_sort)  #put three different data frames into one
  tab_top <- table(col_data$topic_2_sort)  #output the column topic_2_sort in col_data as table
  print(tab_top)
  colnames(de_df_2_sort_int) <- NULL  #clear the colnames of de_df_2_sort_int
  
  dds <- tryCatch(DESeqDataSetFromMatrix(countData = de_df_2_sort_int, colData = as.data.frame(col_data) ,design = ~ age_2_sort+sex_2_sort+topic_2_sort), error=function(e) "err")  #build the DESeq object for DESeq2
  if (length(dds)==1){  # not using possible_error == "err", because somehow it is very slow
    dds <- tryCatch(DESeqDataSetFromMatrix(countData = de_df_2_sort_int, colData = as.data.frame(col_data) ,design = ~ sex_2_sort+topic_2_sort), error=function(e) "err")
    
    if (length(dds)==1){ 
      dds <- DESeqDataSetFromMatrix(countData = de_df_2_sort_int, colData = as.data.frame(col_data) ,design = ~ topic_2_sort)
    }
  }
  #perform differential expression analysis
  dds <- DESeq(dds,betaPrior=T)  #use wald-test (not LRT) for greater power and smaller p-value; ##!for newer version, need to add betaPrior=T (if F gives MLE estimate of L2FC, T gives MAP estimate of L2FC), in order to be consistent with old runs  --> in order to generate identical pval and l2fc
  #res <- results(dds, alpha = 0.05)
  res <- results(dds, alpha = 0.05/6)  #results of differential expression analysis
  print(mcols(res, use.names=TRUE))
  summary(res)
  # res <- results(dds, alpha = 0.1)
  # summary(res)
  # res <- results(dds, alpha = 0.5)
  # summary(res)
  
  ntop <- 10
  select_id <- intersect(head(order(res$log2FoldChange),ntop), which(res$padj<0.05/6)) #intersect the results of top 10 l2fc and adjust p <0.05/6
  #cat("\tpadj:",res$padj[select_id],"\n")
  #cat("\tpraw:",res$pvalue[select_id],"\n")
  #cat("\tl2fc:",res$log2FoldChange[select_id],"\n")
  identified_de_gene <- row.names(de_df_2_sort_int)[select_id]  #find the rows about results of top 10 l2fc and adjust p <0.05/6
  cat("\tdown_gene:\n\t",paste(identified_de_gene,rep(" --> l2fc=",ntop),res$log2FoldChange[select_id],rep(", p_raw=",ntop),res$pvalue[select_id],sep=" ",collapse="; "),"\n")
  
  select_id <- intersect(tail(order(res$log2FoldChange),ntop), which(res$padj<0.05/6))  #intersect the results of bottom 10 l2fc and adjust p <0.05/6
  identified_de_gene <- row.names(de_df_2_sort_int)[select_id]  #find the rows about results of bottom 10 l2fc and adjust p <0.05/6
  cat("\tup_gene:\n\t",paste(identified_de_gene,rep(" --> l2fc=",ntop),res$log2FoldChange[select_id],rep(", p_raw=",ntop),res$pvalue[select_id],sep=" ",collapse="; "),"\n")
  
  #write
  eid_target <- niky_de$gene
  
  eid_sqc <- row.names(res)
  rw_mat_id <- match(eid_target,eid_sqc)  #match the different results with gene manes
  mat_id <- rw_mat_id[!is.na(rw_mat_id)]
  tag_id <- which(!is.na(rw_mat_id))
  
  niky_de[[paste("ncase_",top2_i,"_agst_",top1_i,sep="")]][tag_id] <- tab_top[2]
  niky_de[[paste("nctrl_",top2_i,"_agst_",top1_i,sep="")]][tag_id] <- tab_top[1]
  niky_de[[paste("baseMean_",top2_i,"_agst_",top1_i,sep="")]][tag_id] <- res$baseMean[mat_id]
  niky_de[[paste("l2fc_",top2_i,"_agst_",top1_i,sep="")]][tag_id] <- res$log2FoldChange[mat_id]
  niky_de[[paste("lfcSE_",top2_i,"_agst_",top1_i,sep="")]][tag_id] <- res$lfcSE[mat_id]
  niky_de[[paste("pvalue_",top2_i,"_agst_",top1_i,sep="")]][tag_id] <- res$pvalue[mat_id]
  niky_de[[paste("padj_",top2_i,"_agst_",top1_i,sep="")]][tag_id] <- res$padj[mat_id]
  
  
  ##10k genes on small collection of sample
  de_df_2_sort_int_2 <- de_df_2_sort_int_raw[setdiff(row.names(de_df_2_sort_int_raw),row.names(de_df_2_sort_int)),]
  keep_col <- !apply(de_df_2_sort_int_2[1,],2,is.na)  #find which columns of first line is not blank 
  de_df_2_sort_int_2 <- de_df_2_sort_int_2[,keep_col]
  de_df_2_sort_int_2 <- na.omit(de_df_2_sort_int_2)
  de_df_2_sort_int_2 <- de_df_2_sort_int_2[rowSums(de_df_2_sort_int_2)>=5,]
  
  age_2_sort_kp_cut <- cut(age_2_sort[keep_col],breaks=c(min(age_2_sort[keep_col])-1,median(age_2_sort[keep_col]),max(age_2_sort[keep_col])+1))
  col_data <- DataFrame(topic_2_sort=topic_2_sort[keep_col],age_2_sort=age_2_sort_kp_cut,sex_2_sort=sex_2_sort[keep_col])
  tab_top <- table(col_data$topic_2_sort)
  print(tab_top)
  colnames(de_df_2_sort_int_2) <- NULL
  
  dds <- tryCatch(DESeqDataSetFromMatrix(countData = de_df_2_sort_int_2, colData = as.data.frame(col_data) ,design = ~ age_2_sort+sex_2_sort+topic_2_sort), error=function(e) "err")
  if (length(dds)==1){  # not using possible_error == "err", because somehow it is very slow.
    dds <- tryCatch(DESeqDataSetFromMatrix(countData = de_df_2_sort_int_2, colData = as.data.frame(col_data) ,design = ~ sex_2_sort+topic_2_sort), error=function(e) "err")
    
    if (length(dds)==1){
      dds <- DESeqDataSetFromMatrix(countData = de_df_2_sort_int_2, colData = as.data.frame(col_data) ,design = ~ topic_2_sort)
    }
  }
  dds <- DESeq(dds,betaPrior=T)  #use wald-test (not LRT) for greater power and smaller p-value; ##!for newer version, need to add betaPrior=T (if F gives MLE estimate of L2FC, T gives MAP estimate of L2FC), in order to be consistent with old runs  --> in order to generate identical pval and l2fc
  #res <- results(dds, alpha = 0.05)
  res <- results(dds, alpha = 0.05/6)
  print(mcols(res, use.names=TRUE))
  summary(res)
  # res <- results(dds, alpha = 0.1)
  # summary(res)
  # res <- results(dds, alpha = 0.5)
  # summary(res)
  
  ntop <- 10
  select_id <- intersect(head(order(res$log2FoldChange),ntop), which(res$padj<0.05/6))
  
  identified_de_gene <- row.names(de_df_2_sort_int_2)[select_id]
  cat("\tdown_gene:\n\t",paste(identified_de_gene,rep(" --> l2fc=",ntop),res$log2FoldChange[select_id],rep(", p_raw=",ntop),res$pvalue[select_id],sep=" ",collapse="; "),"\n")
  
  select_id <- intersect(tail(order(res$log2FoldChange),ntop), which(res$padj<0.05/6))
  identified_de_gene <- row.names(de_df_2_sort_int_2)[select_id]
  cat("\tup_gene:\n\t",paste(identified_de_gene,rep(" --> l2fc=",ntop),res$log2FoldChange[select_id],rep(", p_raw=",ntop),res$pvalue[select_id],sep=" ",collapse="; "),"\n")
  
  #write
  eid_target <- niky_de$gene
  
  eid_sqc <- row.names(res)
  rw_mat_id <- match(eid_target,eid_sqc)
  mat_id <- rw_mat_id[!is.na(rw_mat_id)]
  tag_id <- which(!is.na(rw_mat_id))
  
  niky_de[[paste("ncase_",top2_i,"_agst_",top1_i,sep="")]][tag_id] <- tab_top[2]
  niky_de[[paste("nctrl_",top2_i,"_agst_",top1_i,sep="")]][tag_id] <- tab_top[1]
  niky_de[[paste("baseMean_",top2_i,"_agst_",top1_i,sep="")]][tag_id] <- res$baseMean[mat_id]
  niky_de[[paste("l2fc_",top2_i,"_agst_",top1_i,sep="")]][tag_id] <- res$log2FoldChange[mat_id]
  niky_de[[paste("lfcSE_",top2_i,"_agst_",top1_i,sep="")]][tag_id] <- res$lfcSE[mat_id]
  niky_de[[paste("pvalue_",top2_i,"_agst_",top1_i,sep="")]][tag_id] <- res$pvalue[mat_id]
  niky_de[[paste("padj_",top2_i,"_agst_",top1_i,sep="")]][tag_id] <- res$padj[mat_id]
}

##read in gene annotate data
load("transcriptome_differential_analysis_input_2.RData")

##gene annotate 1
#gene annotate, weather the genes are immune or not
niky_de$gene_annotate <- NA
eid_target <- niky_de$gene
eid_sqc <- gene_anno_db$基因名称
rw_mat_id <- match(eid_target,eid_sqc)  #find out which genes are immune genes
mat_id <- rw_mat_id[!is.na(rw_mat_id)]
tag_id <- which(!is.na(rw_mat_id))

niky_de$gene_annotate[tag_id] <- gene_anno_db$Function[mat_id]
length(which(!is.na(niky_de$gene_annotate)))   #113 immune genes

##gene annotate 2
niky_de$gene_annotate_2 <- NA
eid_target <- niky_de$gene
eid_sqc <- Q3_Normalized$TargetName
rw_mat_id <- match(eid_target,eid_sqc)
mat_id <- rw_mat_id[!is.na(rw_mat_id)]
tag_id <- which(!is.na(rw_mat_id))

niky_de$gene_annotate_2[tag_id] <- Q3_Normalized$TargetGroup[mat_id]
length(which(!is.na(niky_de$gene_annotate_2)))  #1796
length(which(!(is.na(niky_de$gene_annotate_2) & is.na(niky_de$gene_annotate))))  #1892  #total number of annotated genes

##gene annotate 3
niky_de$in_GeneList <- "No"
eid_target <- niky_de$gene
eid_sqc <- GeneList$Gene_List
rw_mat_id <- match(eid_target,eid_sqc)
mat_id <- rw_mat_id[!is.na(rw_mat_id)]
tag_id <- which(!is.na(rw_mat_id))

niky_de$in_GeneList[tag_id] <- "Yes"
length(which(niky_de$in_GeneList=="Yes"))  #2168

##recompute padj by BH to avoid NA values
pval_colname <- colnames(niky_de)[grep("^pval",colnames(niky_de))]
for (colname_i in pval_colname){
  bh_colname_i <- paste("BH_adjusted_",colname_i,sep="")
  niky_de[[bh_colname_i]] <- NA
  idx_na_i <- which(!is.na(niky_de[[colname_i]]))
  niky_de[[bh_colname_i]][idx_na_i] <- p.adjust(niky_de[[colname_i]][idx_na_i], method = "BH", n = length(niky_de[[colname_i]][idx_na_i]))
}

display_df <- niky_de[,c("padj_1NKatStroma_agst_0NKatBorder","BH_adjusted_pvalue_1NKatStroma_agst_0NKatBorder","padj_1NKatTumor_agst_0NKatBorder","BH_adjusted_pvalue_1NKatTumor_agst_0NKatBorder")]
display_df <- niky_de[,c("padj_1NKatStroma_agst_0Nonrelp_NKatBorder","BH_adjusted_pvalue_1NKatStroma_agst_0Nonrelp_NKatBorder","padj_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder","BH_adjusted_pvalue_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder")]
display_df <- niky_de[,c("padj_1NKatStroma_agst_0Relp_NKatBorder","BH_adjusted_pvalue_1NKatStroma_agst_0Relp_NKatBorder","padj_1Relp_NKatTumor_agst_0Relp_NKatBorder","BH_adjusted_pvalue_1Relp_NKatTumor_agst_0Relp_NKatBorder")]
dim(niky_de)  # 18677    84
length(which(is.na(niky_de[,"BH_adjusted_pvalue_1NKatStroma_agst_0NKatBorder"])))
# 285
length(which(is.na(niky_de[,"BH_adjusted_pvalue_1NKatTumor_agst_0NKatBorder"])))
# 286
length(union(which(is.na(niky_de[,"BH_adjusted_pvalue_1NKatStroma_agst_0NKatBorder"])),which(is.na(niky_de[,"BH_adjusted_pvalue_1NKatTumor_agst_0NKatBorder"]))))
# 287

length(which(is.na(niky_de[,"BH_adjusted_pvalue_1NKatStroma_agst_0Nonrelp_NKatBorder"])))
# 297
length(which(is.na(niky_de[,"BH_adjusted_pvalue_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder"])))
# 290
length(union(which(is.na(niky_de[,"BH_adjusted_pvalue_1NKatStroma_agst_0Nonrelp_NKatBorder"])),which(is.na(niky_de[,"BH_adjusted_pvalue_1Nonrelp_NKatTumor_agst_0Nonrelp_NKatBorder"]))))
# 303

length(which(is.na(niky_de[,"BH_adjusted_pvalue_1NKatStroma_agst_0Relp_NKatBorder"])))
# 311
length(which(is.na(niky_de[,"BH_adjusted_pvalue_1Relp_NKatTumor_agst_0Relp_NKatBorder"])))
# 312
length(union(which(is.na(niky_de[,"BH_adjusted_pvalue_1NKatStroma_agst_0Relp_NKatBorder"])),which(is.na(niky_de[,"BH_adjusted_pvalue_1Relp_NKatTumor_agst_0Relp_NKatBorder"]))))
# 319

length(which(is.na(niky_de[,"pvalue_1Relp_NKatTumor_agst_0Relp_NKatBorder"])))
#312

