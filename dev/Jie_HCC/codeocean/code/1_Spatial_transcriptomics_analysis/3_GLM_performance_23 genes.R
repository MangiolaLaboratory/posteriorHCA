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

#other attached packages:
 #[1] pROC_1.16.2         coop_0.6-3          mltools_0.3.5       fmsb_0.7.6         
 #[5] rcompanion_2.4.35   readxl_1.4.3        readr_2.1.5         fitdistrplus_1.1-11
 #[9] survival_3.1-12     MASS_7.3-60.0.1    

#loaded via a namespace (and not attached):
 #[1] nortest_1.0-4       zoo_1.8-12          modeltools_0.2-23   tidyselect_1.2.0   
 #[5] coin_1.4-3          splines_4.2.1       lattice_0.20-45     vctrs_0.6.3        
 #[9] generics_0.1.3      expm_0.999-9        stats4_4.2.1        utf8_1.2.3         
#[13] rlang_1.1.1         e1071_1.7-14        pillar_1.9.0        glue_1.6.2         
#[17] withr_3.0.0         plyr_1.8.9          matrixStats_1.2.0   multcomp_1.4-25    
#[21] lifecycle_1.0.4     rootSolve_1.8.2.4   multcompView_0.1-10 cellranger_1.1.0   
#[25] mvtnorm_1.2-4       codetools_0.2-18    tzdb_0.4.0          lmom_3.0           
#[29] lmtest_0.9-40       parallel_4.2.1      class_7.3-20        fansi_1.0.4        
#[33] TH.data_1.1-2       Rcpp_1.0.11         gld_2.6.6           Exact_3.2          
#[37] hms_1.1.3           dplyr_1.1.2         grid_4.2.1          cli_3.6.1          
#[41] tools_4.2.1         sandwich_3.1-0      magrittr_2.0.3      DescTools_0.99.54  
#[45] proxy_0.4-27        tibble_3.2.1        pkgconfig_2.0.3     libcoin_1.0-10     
#[49] Matrix_1.6-5        data.table_1.15.2   httr_1.4.7          rstudioapi_0.15.0  
#[53] R6_2.5.1            boot_1.3-28         compiler_4.2.1    


setwd("/code/1_Spatial_transcriptomics_analysis/")

library(fitdistrplus);library(readr);library(readxl);library(rcompanion);library(fmsb);library(mltools);library(coop);library(pROC)

##read in data
load("/data/Input_1_Spatial_transcriptomics_analysis/transcriptome_differential_analysis_input_4.RData")

##global parameters
accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}
N_trial <- 100
coef <- qt(.975, df = N_trial-1)
ppv <- function(x){diag(x)[2]/rowSums(x)[2] * 100}
npv <- function(x){diag(x)[1]/rowSums(x)[1] * 100}

######################## ML: genes predict for liver ############################

####ML 1: reurrent vs non-recurrent, predict liver using absolute value, one random sampling####

#training and testing vectors
test_acc_gene <- c()

xpluscova_acc_train_mean <- c(); xpluscova_acc_train_sd <- c() 
xpluscova_acc_test_mean <- c(); xpluscova_acc_test_sd <- c()
xpluscova_nag_r2_train_mean <- c(); xpluscova_nag_r2_train_sd <- c()
xpluscova_nag_pval_train_mean <- c(); xpluscova_nag_pval_train_sd <- c()
xpluscova_kapa_r_train_mean <- c(); xpluscova_kapa_r_train_sd <- c()
xpluscova_kapa_r_test_mean <- c(); xpluscova_kapa_r_test_sd <- c()
xpluscova_kapa_pval_train_mean <- c(); xpluscova_kapa_pval_train_sd <- c()
xpluscova_kapa_pval_test_mean <- c(); xpluscova_kapa_pval_test_sd <- c()
xpluscova_mccr_r_train_mean <- c(); xpluscova_mccr_r_train_sd <- c()
xpluscova_mccr_r_test_mean <- c(); xpluscova_mccr_r_test_sd <- c()

xpluscova_ppv_train_mean <- c();xpluscova_ppv_train_sd <- c()
xpluscova_ppv_test_mean <- c();xpluscova_ppv_test_sd <- c()

xpluscova_npv_train_mean <- c();xpluscova_npv_train_sd <- c()
xpluscova_npv_test_mean <- c();xpluscova_npv_test_sd <- c()

xpluscova_auc_train_mean <- c();xpluscova_auc_train_sd <- c()
xpluscova_auc_test_mean <- c();xpluscova_auc_test_sd <- c()


##use negative binomial distribution to expand gene reading data
for (i in 1:length(gene_show)){
  gene_i <- gene_show[i]
  type_i <- unlist(gene_name[which(gene_name$gene==gene_i),c("nonrelp","relp")])
  
  NKatStroma_i <- na.omit(unlist(data_gene_df_ed[gene_i,paste(rep("...",length(NKatStroma_indx)),NKatStroma_indx,sep="")]))
  
  Relp_NKatBorder_i <- na.omit(unlist(data_gene_df_ed[gene_i,paste(rep("...",length(Relp_NKatBorder_indx)),Relp_NKatBorder_indx,sep="")]))
  Relp_NKatTumor_i <- na.omit(unlist(data_gene_df_ed[gene_i,paste(rep("...",length(Relp_NKatTumor_indx)),Relp_NKatTumor_indx,sep="")]))
  
  Nonrelp_NKatBorder_i <- na.omit(unlist(data_gene_df_ed[gene_i,paste(rep("...",length(Nonrelp_NKatBorder_indx)),Nonrelp_NKatBorder_indx,sep="")]))
  Nonrelp_NKatTumor_i <- na.omit(unlist(data_gene_df_ed[gene_i,paste(rep("...",length(Nonrelp_NKatTumor_indx)),Nonrelp_NKatTumor_indx,sep="")]))
  cat(gene_i,"\tnon-relp vs relp:",type_i,"\tnon-relp border=",Nonrelp_NKatBorder_i,"\tnon-relp tumor=",Nonrelp_NKatTumor_i,"\trelp border=",Relp_NKatBorder_i,"\trelp tumor=",Relp_NKatTumor_i,"\n")
  
  nb_Relp_NKatBorder_i <- fitdist(round(as.numeric(Relp_NKatBorder_i)),"nbinom")
  nb_Relp_NKatTumor_i <- fitdist(round(as.numeric(Relp_NKatTumor_i)), "nbinom")
  nb_Nonrelp_NKatBorder_i <- fitdist(round(as.numeric(Nonrelp_NKatBorder_i)), "nbinom")
  nb_Nonrelp_NKatTumor_i <- fitdist(round(as.numeric(Nonrelp_NKatTumor_i)), "nbinom")
  
  xpluscova_acc_train_vec <- c(); xpluscova_acc_test_vec <- c()
  xpluscova_nag_r2_train_vec <- c(); xpluscova_nag_pval_train_vec <- c()
  xpluscova_kapa_r_train_vec <- c(); xpluscova_kapa_r_test_vec <- c()
  xpluscova_kapa_pval_train_vec <- c(); xpluscova_kapa_pval_test_vec <- c()
  xpluscova_mccr_r_train_vec <- c(); xpluscova_mccr_r_test_vec <- c()
  xpluscova_ppv_train_vec <- c(); xpluscova_ppv_test_vec <- c()
  xpluscova_npv_train_vec <- c(); xpluscova_npv_test_vec <- c()
  xpluscova_auc_train_vec <- c(); xpluscova_auc_test_vec <- c()
  
  
  for (trial_j in 1:N_trial){
    
    set.seed(trial_j)
    nb_bs_Relp_NKatBorder_i <- rnbinom(1000, size=nb_Relp_NKatBorder_i$estimate["size"], mu=nb_Relp_NKatBorder_i$estimate["mu"])
    set.seed(trial_j)
    nb_bs_Relp_NKatTumor_i <- rnbinom(1000, size=nb_Relp_NKatTumor_i$estimate["size"], mu=nb_Relp_NKatTumor_i$estimate["mu"])
    set.seed(trial_j)
    nb_bs_Nonrelp_NKatBorder_i <- rnbinom(1000, size=nb_Nonrelp_NKatBorder_i$estimate["size"], mu=nb_Nonrelp_NKatBorder_i$estimate["mu"])
    set.seed(trial_j)
    nb_bs_Nonrelp_NKatTumor_i <- rnbinom(1000, size=nb_Nonrelp_NKatTumor_i$estimate["size"], mu=nb_Nonrelp_NKatTumor_i$estimate["mu"])
    
    train_df <- data.frame(Relp_y=c(rep("1",750),rep("0",750)),NKatBorder_x=c(nb_bs_Relp_NKatBorder_i[1:750],nb_bs_Nonrelp_NKatBorder_i[1:750]),NKatTumor_x=c(nb_bs_Relp_NKatTumor_i[1:750],nb_bs_Nonrelp_NKatTumor_i[1:750]),stringsAsFactors=T)
    test_df <- data.frame(Relp_y=c(rep("1",250),rep("0",250)),NKatBorder_x=c(nb_bs_Relp_NKatBorder_i[751:1000],nb_bs_Nonrelp_NKatBorder_i[751:1000]),NKatTumor_x=c(nb_bs_Relp_NKatTumor_i[751:1000],nb_bs_Nonrelp_NKatTumor_i[751:1000]),stringsAsFactors=T)
    train_target <- train_df$Relp_y
    test_target <- test_df$Relp_y
    
    ##training: NKatBorder_x + NKatTumor_x
    set.seed(trial_j)
    linearMod <- glm("Relp_y ~ NKatBorder_x + NKatTumor_x ",data = train_df,family=binomial)
    
    pred_prob <- predict(linearMod, train_df, type = "response") #it will look for prs_xpluscova_name features automatically.
    pred_label <- ifelse(pred_prob > 0.5,1,0)
    pred_label <- as.factor(pred_label)
    tb <- table(pred_label, train_target)
    xpluscova_acc_train_vec <- c(xpluscova_acc_train_vec, accuracy(tb))
    
    r2_res <- nagelkerke(linearMod)
    xpluscova_nag_r2_train_vec <- c(xpluscova_nag_r2_train_vec, r2_res$Pseudo.R.squared.for.model.vs.null[3]); xpluscova_nag_pval_train_vec <- c(xpluscova_nag_pval_train_vec, r2_res$Likelihood.ratio.test[4])
    
    cor_test <- Kappa.test(x=train_target, y=pred_label, conf.level=0.95) #calculate Cohen's kappa statistics for agreement and its confidence intervals followed by testing null-hypothesis that the extent of agreement is same as random, kappa statistic equals zero.
    xpluscova_kapa_r_train_vec <- c(xpluscova_kapa_r_train_vec, cor_test$Result$estimate); xpluscova_kapa_pval_train_vec <- c(xpluscova_kapa_pval_train_vec, cor_test$Result$p.value)
    
    xpluscova_mccr_r_train_vec <- c(xpluscova_mccr_r_train_vec, mcc(preds = pred_label, actuals = train_target))
    
    xpluscova_auc_train_vec<- c(xpluscova_auc_train_vec,auc(train_target,pred_prob))
    xpluscova_ppv_train_vec <- c(xpluscova_ppv_train_vec, ppv(tb))
    xpluscova_npv_train_vec <- c(xpluscova_npv_train_vec, npv(tb))
    
    ##testing: NKatBorder_x + NKatTumor_x
    pred_prob <- predict(linearMod, test_df, type = "response") #it will look for prs_xpluscova_name features automatically.
    pred_label <- ifelse(pred_prob > 0.5,1,0)
    pred_label <- as.factor(pred_label)
    tb <- table(pred_label, test_target)
    xpluscova_acc_test_vec <- c(xpluscova_acc_test_vec, accuracy(tb))
    
    cor_test <- Kappa.test(x=test_target, y=pred_label, conf.level=0.95) #calculate Cohen's kappa statistics for agreement and its confidence intervals followed by testing null-hypothesis that the extent of agreement is same as random, kappa statistic equals zero.
    xpluscova_kapa_r_test_vec <- c(xpluscova_kapa_r_test_vec, cor_test$Result$estimate); xpluscova_kapa_pval_test_vec <- c(xpluscova_kapa_pval_test_vec, cor_test$Result$p.value)
    
    xpluscova_mccr_r_test_vec <- c(xpluscova_mccr_r_test_vec, mcc(preds = pred_label, actuals = test_target))
    
    xpluscova_auc_test_vec<- c(xpluscova_auc_test_vec,auc(test_target,pred_prob))
    xpluscova_ppv_test_vec <- c(xpluscova_ppv_test_vec, ppv(tb))
    xpluscova_npv_test_vec <- c(xpluscova_npv_test_vec, npv(tb))
    
  }
  
  
  cat("gene: ",gene_i,"\n")
  cat("\n\tNKatBorder_x + NKatTumor_x: train --> ",mean(xpluscova_acc_train_vec),"+/- ",coef*sd(xpluscova_acc_train_vec)/sqrt(N_trial),"; test --> ",mean(xpluscova_acc_test_vec),"+/- ",coef*sd(xpluscova_acc_test_vec)/sqrt(N_trial),"\n")
  cat("\t\t train nag r2 --> ",mean(xpluscova_nag_r2_train_vec),"+/- ",coef*sd(xpluscova_nag_r2_train_vec)/sqrt(N_trial),"\n")
  cat("\t\t train nag pval --> ",mean(xpluscova_nag_pval_train_vec),"+/- ",coef*sd(xpluscova_nag_pval_train_vec)/sqrt(N_trial),"\n")
  cat("\t\t train kapa r --> ",mean(xpluscova_kapa_r_train_vec),"+/- ",coef*sd(xpluscova_kapa_r_train_vec)/sqrt(N_trial),"; test kapa r --> ",mean(xpluscova_kapa_r_test_vec),"+/- ",coef*sd(xpluscova_kapa_r_test_vec)/sqrt(N_trial),"\n")
  cat("\t\t train kapa pval --> ",mean(xpluscova_kapa_pval_train_vec),"+/- ",coef*sd(xpluscova_kapa_pval_train_vec)/sqrt(N_trial),"; test kapa pval --> ",mean(xpluscova_kapa_pval_test_vec),"+/- ",coef*sd(xpluscova_kapa_pval_test_vec)/sqrt(N_trial),"\n")
  cat("\t\t train mccr r --> ",mean(xpluscova_mccr_r_train_vec),"+/- ",coef*sd(xpluscova_mccr_r_train_vec)/sqrt(N_trial),"; test mccr r --> ",mean(xpluscova_mccr_r_test_vec),"+/- ",coef*sd(xpluscova_mccr_r_test_vec)/sqrt(N_trial),"\n")
  cat("\t\t train ppv --> ",mean(xpluscova_ppv_train_vec),"+/- ",coef*sd(xpluscova_ppv_train_vec)/sqrt(N_trial),"; test ppv --> ",mean(xpluscova_ppv_test_vec),"+/- ",coef*sd(xpluscova_ppv_test_vec)/sqrt(N_trial),"\n")
  cat("\t\t train npv --> ",mean(xpluscova_npv_train_vec),"+/- ",coef*sd(xpluscova_npv_train_vec)/sqrt(N_trial),"; test npv --> ",mean(xpluscova_npv_test_vec),"+/- ",coef*sd(xpluscova_npv_test_vec)/sqrt(N_trial),"\n")
  cat("\t\t train auc --> ",mean(xpluscova_auc_train_vec),"+/- ",coef*sd(xpluscova_auc_train_vec)/sqrt(N_trial),"; test auc --> ",mean(xpluscova_auc_test_vec),"+/- ",coef*sd(xpluscova_auc_test_vec)/sqrt(N_trial),"\n")
  cat("\n##\n")
#RNASE6 	non-relp vs relp: flatdown downflat 	non-relp border= 37.17816 46.17078 25.67298 37.926 27.04288 35.3976 41.03468 27.85385 44.18046 	non-relp tumor= 14.26017 8.63356 8.369875 	relp border= 12.77507 20.5884 13.4848 16.254 11.9196 8.274763 7.585199 30.3408 20.2272 43.8256 10.1136 	relp tumor= 30.3408 33.712 11.3778 17.3376 

  #......

#gene:  S100A4 

	#NKatBorder_x + NKatTumor_x: train -->  97.80667 +/-  0.06709571 ; test -->  97.716 +/-  0.1257689 
		 #train nag r2 -->  0.9545151 +/-  0.001438061 
		 #train nag pval -->  0 +/-  0 
		 #train kapa r -->  0.9561333 +/-  0.001341914 ; test kapa r -->  0.95432 +/-  0.002515379 
		 #train kapa pval -->  0 +/-  0 ; test kapa pval -->  0 +/-  0 
		 #train mccr r -->  0.9562493 +/-  0.001334836 ; test mccr r -->  0.9545279 +/-  0.002485619 
		 #train ppv -->  97.12029 +/-  0.09988053 ; test ppv -->  97.02241 +/-  0.2165756 
		 #train npv -->  98.51623 +/-  0.06402974 ; test npv -->  98.45119 +/-  0.1448398 
		 #train auc -->  0.9963705 +/-  0.0002349166 ; test auc -->  0.9964495 +/-  0.0004281259 

##
#GPR183 	non-relp vs relp: flatdown downflat 	non-relp border= 51.28022 46.17078 29.17384 54.99269 33.63871 42.9828 26.36168 37.80165 39.38981 	non-relp tumor= 20.02493 6.906848 9.677668 	relp border= 22.35638 29.2572 16.856 24.9228 8.668799 22.06603 15.1704 16.856 35.3976 16.856 23.5984 	relp tumor= 28.1736 57.31039 3.7926 21.672 

#gene:  GPR183 

	#NKatBorder_x + NKatTumor_x: train -->  90.602 +/-  0.1389089 ; test -->  90.424 +/-  0.2614948 
		 #train nag r2 -->  0.8081581 +/-  0.00237914 
		 #train nag pval -->  3.664715e-291 +/-  0 
		 #train kapa r -->  0.81204 +/-  0.002778179 ; test kapa r -->  0.80848 +/-  0.005229896 
		 #train kapa pval -->  0 +/-  0 ; test kapa pval -->  0 +/-  0 
		 #train mccr r -->  0.8120773 +/-  0.002776588 ; test mccr r -->  0.8087937 +/-  0.005233904 
		 #train ppv -->  90.72626 +/-  0.1438242 ; test ppv -->  90.56654 +/-  0.3359998 
		 #train npv -->  90.4852 +/-  0.169315 ; test npv -->  90.34423 +/-  0.3530666 
		 #train auc -->  0.9684496 +/-  0.0006258956 ; test auc -->  0.9677255 +/-  0.001411405 

  xpluscova_acc_train_mean <- c(xpluscova_acc_train_mean, mean(xpluscova_acc_train_vec)); xpluscova_acc_train_sd <- c(xpluscova_acc_train_sd, sd(xpluscova_acc_train_vec)) 
  xpluscova_acc_test_mean <- c(xpluscova_acc_test_mean, mean(xpluscova_acc_test_vec)); xpluscova_acc_test_sd <- c(xpluscova_acc_test_sd, sd(xpluscova_acc_test_vec))
  xpluscova_nag_r2_train_mean <- c(xpluscova_nag_r2_train_mean, mean(xpluscova_nag_r2_train_vec)); xpluscova_nag_r2_train_sd <- c(xpluscova_nag_r2_train_sd, sd(xpluscova_nag_r2_train_vec))
  xpluscova_nag_pval_train_mean <- c(xpluscova_nag_pval_train_mean, mean(xpluscova_nag_pval_train_vec)); xpluscova_nag_pval_train_sd <- c(xpluscova_nag_pval_train_sd, sd(xpluscova_nag_pval_train_vec))
  xpluscova_kapa_r_train_mean <- c(xpluscova_kapa_r_train_mean, mean(xpluscova_kapa_r_train_vec)); xpluscova_kapa_r_train_sd <- c(xpluscova_kapa_r_train_sd, sd(xpluscova_kapa_r_train_vec))
  xpluscova_kapa_r_test_mean <- c(xpluscova_kapa_r_test_mean, mean(xpluscova_kapa_r_test_vec)); xpluscova_kapa_r_test_sd <- c(xpluscova_kapa_r_test_sd, sd(xpluscova_kapa_r_test_vec))
  xpluscova_kapa_pval_train_mean <- c(xpluscova_kapa_pval_train_mean, mean(xpluscova_kapa_pval_train_vec)); xpluscova_kapa_pval_train_sd <- c(xpluscova_kapa_pval_train_sd, sd(xpluscova_kapa_pval_train_vec))
  xpluscova_kapa_pval_test_mean <- c(xpluscova_kapa_pval_test_mean, mean(xpluscova_kapa_pval_test_vec)); xpluscova_kapa_pval_test_sd <- c(xpluscova_kapa_pval_test_sd, sd(xpluscova_kapa_pval_test_vec))
  xpluscova_mccr_r_train_mean <- c(xpluscova_mccr_r_train_mean, mean(xpluscova_mccr_r_train_vec)); xpluscova_mccr_r_train_sd <- c(xpluscova_mccr_r_train_sd, sd(xpluscova_mccr_r_train_vec))
  xpluscova_mccr_r_test_mean <- c(xpluscova_mccr_r_test_mean, mean(xpluscova_mccr_r_test_vec)); xpluscova_mccr_r_test_sd <- c(xpluscova_mccr_r_test_sd, sd(xpluscova_mccr_r_test_vec))
 
  xpluscova_ppv_train_mean <- c( xpluscova_ppv_train_mean, mean( xpluscova_ppv_train_vec));  xpluscova_ppv_train_sd <- c( xpluscova_ppv_train_sd, sd( xpluscova_ppv_train_vec))
  xpluscova_ppv_test_mean <- c( xpluscova_ppv_test_mean, mean( xpluscova_ppv_test_vec));  xpluscova_ppv_test_sd <- c( xpluscova_ppv_test_sd, sd( xpluscova_ppv_test_vec))
  xpluscova_npv_train_mean <- c( xpluscova_npv_train_mean, mean( xpluscova_npv_train_vec));  xpluscova_npv_train_sd <- c( xpluscova_npv_train_sd, sd( xpluscova_npv_train_vec))
  xpluscova_npv_test_mean <- c( xpluscova_npv_test_mean, mean( xpluscova_npv_test_vec));  xpluscova_npv_test_sd <- c( xpluscova_npv_test_sd, sd( xpluscova_npv_test_vec))
  xpluscova_auc_train_mean <- c( xpluscova_auc_train_mean, mean( xpluscova_auc_train_vec));  xpluscova_auc_train_sd <- c( xpluscova_auc_train_sd, sd( xpluscova_auc_train_vec))
  xpluscova_auc_test_mean <- c( xpluscova_auc_test_mean, mean( xpluscova_auc_test_vec));  xpluscova_auc_test_sd <- c( xpluscova_auc_test_sd, sd( xpluscova_auc_test_vec)) 
   
}


predLiver_relapse_absoluteX <- data.frame(gene=gene_show,xpluscova_acc_train_mean=xpluscova_acc_train_mean, xpluscova_acc_train_sd=xpluscova_acc_train_sd,  
                                          xpluscova_acc_test_mean=xpluscova_acc_test_mean, xpluscova_acc_test_sd=xpluscova_acc_test_sd, 
                                          xpluscova_nag_r2_train_mean=xpluscova_nag_r2_train_mean, xpluscova_nag_r2_train_sd=xpluscova_nag_r2_train_sd, 
                                          xpluscova_nag_pval_train_mean=xpluscova_nag_pval_train_mean, xpluscova_nag_pval_train_sd=xpluscova_nag_pval_train_sd, 
                                          xpluscova_kapa_r_train_mean=xpluscova_kapa_r_train_mean, xpluscova_kapa_r_train_sd=xpluscova_kapa_r_train_sd, 
                                          xpluscova_kapa_r_test_mean=xpluscova_kapa_r_test_mean, xpluscova_kapa_r_test_sd=xpluscova_kapa_r_test_sd, 
                                          xpluscova_kapa_pval_train_mean=xpluscova_kapa_pval_train_mean, xpluscova_kapa_pval_train_sd=xpluscova_kapa_pval_train_sd, 
                                          xpluscova_kapa_pval_test_mean=xpluscova_kapa_pval_test_mean, xpluscova_kapa_pval_test_sd=xpluscova_kapa_pval_test_sd, 
                                          xpluscova_mccr_r_train_mean=xpluscova_mccr_r_train_mean, xpluscova_mccr_r_train_sd=xpluscova_mccr_r_train_sd, 
                                          xpluscova_mccr_r_test_mean=xpluscova_mccr_r_test_mean, xpluscova_mccr_r_test_sd=xpluscova_mccr_r_test_sd, 
                                          xpluscova_ppv_train_mean=xpluscova_ppv_train_mean,xpluscova_ppv_train_sd=xpluscova_ppv_train_sd,
                                          xpluscova_ppv_test_mean=xpluscova_ppv_test_mean,xpluscova_ppv_test_sd=xpluscova_ppv_test_sd,
                                          xpluscova_npv_train_mean=xpluscova_npv_train_mean,xpluscova_npv_train_sd=xpluscova_npv_train_sd,
                                          xpluscova_npv_test_mean=xpluscova_npv_test_mean,xpluscova_npv_test_sd=xpluscova_npv_test_sd,
                                          xpluscova_auc_train_mean=xpluscova_auc_train_mean,xpluscova_auc_train_sd=xpluscova_auc_train_sd,
                                          xpluscova_auc_test_mean=xpluscova_auc_test_mean,xpluscova_auc_test_sd=xpluscova_auc_test_sd,
                                          stringsAsFactors = F)


#write_csv(predLiver_relapse_absoluteX,"LIVER_TRANSCRIPTOME.csv")