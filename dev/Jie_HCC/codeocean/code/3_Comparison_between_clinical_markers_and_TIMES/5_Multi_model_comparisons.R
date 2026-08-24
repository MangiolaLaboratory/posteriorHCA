# > sessionInfo()
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.5 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ropls_1.30.0    pROC_1.18.5     dplyr_1.1.4     survival_3.1-12
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.12                 XVector_0.38.0              pillar_1.9.0                compiler_4.2.1              plyr_1.8.9                  GenomeInfoDb_1.34.9         zlibbioc_1.44.0            
# [8] MatrixGenerics_1.10.0       bitops_1.0-7                tools_4.2.1                 MultiDataSet_1.26.0         lifecycle_1.0.4             tibble_3.2.1                lattice_0.20-45            
# [15] pkgconfig_2.0.3             rlang_1.1.3                 Matrix_1.6-5                DelayedArray_0.24.0         cli_3.6.2                   rstudioapi_0.16.0           parallel_4.2.1             
# [22] qqman_0.1.9                 GenomeInfoDbData_1.2.9      withr_3.0.0                 generics_0.1.3              vctrs_0.6.5                 S4Vectors_0.36.2            IRanges_2.32.0             
# [29] MultiAssayExperiment_1.24.0 stats4_4.2.1                grid_4.2.1                  tidyselect_1.2.1            calibrate_1.7.7             glue_1.7.0                  Biobase_2.58.0             
# [36] R6_2.5.1                    fansi_1.0.6                 limma_3.54.2                magrittr_2.0.3              MASS_7.3-58                 GenomicRanges_1.50.2        matrixStats_1.2.0          
# [43] splines_4.2.1               BiocGenerics_0.44.0         SummarizedExperiment_1.28.0 utf8_1.2.4                  RCurl_1.98-1.14   

library(survival)
library(dplyr)

setwd("/code/3_Comparison_between_clinical_markers_and_TIMES/")

#######load clinic data####
load("/code/3_Comparison_between_clinical_markers_and_TIMES/results_univariate_multivariate_regression.RData")  #output from 3_Univariate_Multivariate_regression_analysis.R 

#######Stage 1: Pre-process other key clinical factors#######
data_multi  <- all_data %>%
  dplyr::select(DFS, Recurrence, tumor_greatest_dimension, tumor_burden, tumor_area_size, macrovascular_invasion, tumor_differentiation_grade, satellite_lesion, dynamic_enhanced_arterial_phase, Alpha_fetoprotein, lymphadenectasis, AST_ALT_abnormality, MVI_M0_vs_M1M2, BCLC_A_vs_BC, AST_abnormality, capsular_invasion, Ki_67, TNM_2level, neoplastic_thrombus, patient_id, TIMES)

feat_vec <- c("tumor_greatest_dimension", "tumor_burden", "tumor_area_size", "macrovascular_invasion", "tumor_differentiation_grade", "satellite_lesion", "dynamic_enhanced_arterial_phase", "Alpha_fetoprotein", "lymphadenectasis", "AST_ALT_abnormality", "MVI_M0_vs_M1M2", "BCLC_A_vs_BC", "AST_abnormality", "capsular_invasion", "Ki_67", "TNM_2level", "neoplastic_thrombus")
data_multi_fill <- data_multi
for (feat_i in feat_vec) {
  val_i <- data_multi_fill[[feat_i]]
  if (typeof(val_i)=="double"){
    fill_i <- median(val_i, na.rm = T)
  } else {
    temp <- table(val_i)
    fill_i <- names(temp)[which.max(temp)]
  }
  cat(feat_i, fill_i, "\n")
  data_multi_fill[which(is.na(val_i)),feat_i] <- fill_i
}
# tumor_greatest_dimension 5 
# tumor_burden 34.518 
# tumor_area_size 1.9875 
# macrovascular_invasion no 
# tumor_differentiation_grade ≤G2 
# satellite_lesion - 
#   dynamic_enhanced_arterial_phase hyperenhancement 
# Alpha_fetoprotein 62.67 
# lymphadenectasis no 
# AST_ALT_abnormality normal 
# MVI_M0_vs_M1M2 yes 
# BCLC_A_vs_BC high 
# AST_abnormality normal 
# capsular_invasion - 
#   Ki_67 0.3 
# TNM_2level ≥II 
# neoplastic_thrombus - 

table(data_fig3ef$TNM_4level)
# I  II III  IV 
# 56  73  18  17 
summary(data_fig3ef$TIMES)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.02612 0.15755 0.40359 0.40821 0.60635 0.94822      23 
table(data_fig3ef$Recurrence)
# N  Y 
# 98 78 
table(data_fig3ef$TNM_2level)
# <II ≥II 
# 56 108 
table(data_fig3ef$BCLC_A_vs_BC)
# low high 
# 26  104 
table(data_fig3ef$BCLC_3level)
# A  B  C 
# 26 19 85  


#####Stage 2: Compare againist more multi-variate models using ROC Curve Analysis#######
#####Plot receiver operating characteristic (ROC) curves for TIMES and other top performing clinical variables. The bigger the AUC (area_size under the curve), the better the predictive performance. Show TIMES with superior AUC.
library(pROC)
##use data_multi_fill and feat_vec in multi-variate regression above
common_name <- intersect(colnames(data_multi_fill),colnames(data_fig3ef))
#"DFS"          "Recurrence"      "BCLC_A_vs_BC" "TNM_2level"   "patient_id"   "TIMES" 

identical(data_multi_fill$DFS,data_fig3ef$DFS)  #TRUE
identical(data_multi_fill$Recurrence,data_fig3ef$Recurrence)  #TRUE
identical(data_multi_fill$TIMES,data_fig3ef$TIMES)  #TRUE
identical(data_multi_fill$patient_id,data_fig3ef$patient_id)  #TRUE

data_4roc <- data_multi_fill
data_4roc <- cbind(data_4roc, data_fig3ef[,setdiff(colnames(data_fig3ef),common_name)])
data_4roc[["TNM staging"]][which(data_4roc[["TNM staging"]]=="NA")] <- NA

for (feat_i in setdiff(colnames(data_fig3ef),common_name)) {
  val_i <- data_4roc[[feat_i]]
  if (typeof(val_i)=="double"){
    fill_i <- median(val_i, na.rm = T)
  } else {
    temp <- table(val_i)
    fill_i <- names(temp)[which.max(temp)]
  }
  cat(feat_i, fill_i, "\n")
  data_4roc[which(is.na(val_i)),feat_i] <- fill_i
}
# TNM staging II 
# BCLC_3level C 
# TNM_4level II 

multi_model_3size <- glm(Recurrence ~ tumor_greatest_dimension + tumor_area_size + tumor_burden, data = data_4roc,family=binomial)
summary(multi_model_3size)
# Call:
#   glm(formula = Recurrence ~ tumor_greatest_dimension + tumor_area_size + 
#         tumor_burden, family = binomial, data = data_4roc)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.0312  -1.0395  -0.7412   1.2179   2.0704  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)              -1.854823   0.530141  -3.499 0.000467 ***
#   tumor_greatest_dimension  0.443435   0.162797   2.724 0.006453 ** 
#   tumor_area_size          -0.453853   0.337605  -1.344 0.178840    
# tumor_burden              0.005976   0.006121   0.976 0.328849    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 241.71  on 175  degrees of freedom
# Residual deviance: 220.40  on 172  degrees of freedom
# (78 observations deleted due to missingness)
# AIC: 228.4
# 
# Number of Fisher Scoring iterations: 4

multi_model_vascular_invasion <- glm(Recurrence ~ macrovascular_invasion + MVI_M0_vs_M1M2, data = data_4roc,family=binomial)
summary(multi_model_vascular_invasion)
# Call:
#   glm(formula = Recurrence ~ macrovascular_invasion + MVI_M0_vs_M1M2, 
#       family = binomial, data = data_4roc)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.6340  -1.0541  -0.9164   1.3060   1.4632  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)                -0.6506     0.3561  -1.827   0.0677 .
# macrovascular_invasionyes   1.3269     0.5522   2.403   0.0163 *
#   MVI_M0_vs_M1M2yes           0.3533     0.4004   0.882   0.3775  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 241.71  on 175  degrees of freedom
# Residual deviance: 233.37  on 173  degrees of freedom
# (78 observations deleted due to missingness)
# AIC: 239.37
# 
# Number of Fisher Scoring iterations: 4

multi_model_noTIMES <- glm(Recurrence ~ tumor_greatest_dimension + tumor_burden + tumor_area_size + macrovascular_invasion + tumor_differentiation_grade + satellite_lesion + dynamic_enhanced_arterial_phase + 
                             Alpha_fetoprotein + lymphadenectasis + AST_ALT_abnormality + MVI_M0_vs_M1M2 + BCLC_A_vs_BC + AST_abnormality + capsular_invasion + Ki_67 + TNM_2level + neoplastic_thrombus, data = data_4roc,family=binomial)
summary(multi_model_noTIMES)
# Call:
#   glm(formula = Recurrence ~ tumor_greatest_dimension + tumor_burden + 
#         tumor_area_size + macrovascular_invasion + tumor_differentiation_grade + 
#         satellite_lesion + dynamic_enhanced_arterial_phase + Alpha_fetoprotein + 
#         lymphadenectasis + AST_ALT_abnormality + MVI_M0_vs_M1M2 + 
#         BCLC_A_vs_BC + AST_abnormality + capsular_invasion + Ki_67 + 
#         TNM_2level + neoplastic_thrombus, family = binomial, data = data_4roc)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.9233  -0.9021  -0.6860   1.0154   2.3132  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)                                       -2.222e+00  7.170e-01  -3.099  0.00194 **
#   tumor_greatest_dimension                           4.490e-01  1.835e-01   2.447  0.01439 * 
#   tumor_burden                                       3.758e-03  6.285e-03   0.598  0.54989   
# tumor_area_size                                   -4.137e-01  3.494e-01  -1.184  0.23641   
# macrovascular_invasionyes                          8.959e-01  7.425e-01   1.207  0.22758   
# tumor_differentiation_gradeG2-G3                   7.947e-01  5.658e-01   1.405  0.16014   
# tumor_differentiation_grade≥G3                     7.432e-02  5.055e-01   0.147  0.88310   
# satellite_lesion+                                 -4.286e-01  5.340e-01  -0.803  0.42224   
# dynamic_enhanced_arterial_phaseuneven enhancement  6.118e-01  6.293e-01   0.972  0.33096   
# Alpha_fetoprotein                                  7.119e-06  1.106e-05   0.644  0.51982   
# lymphadenectasisyes                                1.502e-01  1.149e+00   0.131  0.89597   
# AST_ALT_abnormalityhigh                            2.895e-01  4.415e-01   0.656  0.51194   
# MVI_M0_vs_M1M2yes                                  4.962e-01  6.843e-01   0.725  0.46838   
# BCLC_A_vs_BChigh                                  -2.916e-01  7.530e-01  -0.387  0.69854   
# AST_abnormalityhigh                               -2.777e-01  4.755e-01  -0.584  0.55924   
# capsular_invasion+                                 4.137e-01  4.443e-01   0.931  0.35178   
# Ki_67                                              7.682e-01  1.313e+00   0.585  0.55841   
# TNM_2level≥II                                     -4.602e-01  6.744e-01  -0.682  0.49498   
# neoplastic_thrombus+                               1.332e-01  6.055e-01   0.220  0.82593   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 241.71  on 175  degrees of freedom
# Residual deviance: 205.63  on 157  degrees of freedom
# (78 observations deleted due to missingness)
# AIC: 243.63
# 
# Number of Fisher Scoring iterations: 8

##plot

ROC_times <- roc(data_4roc$Recurrence, data_4roc$TIMES, na.rm=T, smooth=T,algorithm=1, smooth.method="density")
ROC_times_auc <- auc(ROC_times)
# Area under the curve: 0.8232

ROC_TNM <- roc(data_4roc$Recurrence, as.numeric(data_4roc$TNM_2level)-1, na.rm=T, smooth=T,algorithm=1, smooth.method="density")
ROC_TNM_auc <- auc(ROC_TNM)
# Area under the curve: 0.5148

ROC_BCLC <- roc(data_4roc$Recurrence, as.numeric(data_4roc$BCLC_A_vs_BC)-1, na.rm=T, smooth=T,algorithm=1, smooth.method="density")
ROC_BCLC_auc <- auc(ROC_BCLC)
# Area under the curve: 0.529

ROC_vascular_invasion <- roc(multi_model_vascular_invasion$y, multi_model_vascular_invasion$fitted.values, na.rm=T, smooth=T,algorithm=1, smooth.method="density")
ROC_vascular_invasion_auc <- auc(ROC_vascular_invasion)
#Area under the curve: 0.5854

ROC_3size <- roc(multi_model_3size$y, multi_model_3size$fitted.values, na.rm=T, smooth=T,algorithm=1, smooth.method="density")
ROC_3size_auc <- auc(ROC_3size)
#Area under the curve: 0.6845

ROC_multi <- roc(multi_model_noTIMES$y, multi_model_noTIMES$fitted.values, na.rm=T, smooth=T,algorithm=1, smooth.method="density")
ROC_multi_auc <- auc(ROC_multi)
#Area under the curve: 0.7154

plot(ROC_TNM, col = "#559AC6", main = "ROC of TIMES and other models",lwd=3, xlim=c(1,0), ylim=c(0,1))
lines(ROC_BCLC, col = "#E6AB02",lwd=3)
lines(ROC_vascular_invasion, col = "black",lwd=3)
lines(ROC_3size, col = "#7570B3",lwd=3)
lines(ROC_multi, col = "#66A61E",lwd=3)
lines(ROC_times, col = "#E73334",lwd=3)

legend("bottomright", legend=c(paste("TIMES",round(ROC_times_auc,2),sep=": "),paste("17 clinical factors",round(ROC_multi_auc,2),sep=": "),
                               paste("3 tumor size factors",round(ROC_3size_auc,2),sep=": "),paste("2 vascular invasion factors",round(ROC_vascular_invasion_auc,2),sep=": "),
                               paste("BCLC",round(ROC_BCLC_auc,2),sep=": "),paste("TNM",round(ROC_TNM_auc,2),sep=": ")), 
      col=c("#E73334","#66A61E","#7570B3","black","#E6AB02","#559AC6"), bg=NA, lty=1, lwd=3.5, box.lty=0, cex=1)
###save as /data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/Fig3i_ROCs.pdf: compare ROCs of multiple models

#######Stage 3: Discrimination analysis#######
#BiocManager::install("ropls")
library(ropls)

opls_data <- data_multi_fill[,c("Recurrence","TIMES",feat_vec)]
for (feat_i in feat_vec){
  opls_data[[feat_i]] <- as.numeric(opls_data[[feat_i]])
}
opls_data <- opls_data[!is.na(opls_data$Recurrence),]

oplsda <- opls(opls_data[,c("TIMES",feat_vec)], opls_data$Recurrence, predI = 1, orthoI = NA)
# OPLS-DA
# 176 samples x 18 variables and 1 response
# standard scaling of predictors and response(s)
# 13 (0%) NAs
# R2X(cum) R2Y(cum) Q2(cum) RMSEE pre ort pR2Y  pQ2
# Total     0.32    0.354   0.265 0.403   1   1 0.05 0.05
vip <- getVipVn(oplsda)
#vipVn PLS(-DA): Numerical vector of Variable Importance in Projection; OPLS(-DA): Numerical vector of Variable Importance for Prediction (VIP4,p from Galindo-Prieto et al, 2014)
vip_order <- vip[order(vip, decreasing = TRUE)]
vip_order_df <- data.frame("val"=vip_order)
vip_order_df$nam <- names(vip_order)

x <- barplot(height=vip_order_df$val, col=adjustcolor("olivedrab4", alpha.f=0.7),las=2,xaxt="n",ylim=c(0,3), ylab="Variable importance for differentiating REC and non-REC",cex.axis=1, cex.names=1, cex.lab=1)
text(cex=1, x=x+.2, y=c(vip_order_df$val), labels=gsub("_"," ",vip_order_df$nam), xpd=TRUE, srt=60)
abline(h=1,lwd = 2, lty = 2, col = 'black')

#save as /data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/Fig3h_VIP.pdf : VIP plot for OPLSDA
