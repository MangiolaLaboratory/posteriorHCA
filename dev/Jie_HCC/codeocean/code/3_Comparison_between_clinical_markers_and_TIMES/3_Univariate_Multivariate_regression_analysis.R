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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8      
# [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] metafor_2.4-0   Matrix_1.6-5    dplyr_1.1.4     survival_3.1-12
# 
# loaded via a namespace (and not attached):
#   [1] lattice_0.20-45   fansi_1.0.6       withr_3.0.0       utf8_1.2.4        grid_4.2.1        R6_2.5.1          nlme_3.1-157      lifecycle_1.0.4   magrittr_2.0.3    pillar_1.9.0      rlang_1.1.3      
# [12] cli_3.6.2         rstudioapi_0.16.0 vctrs_0.6.5       generics_0.1.3    splines_4.2.1     tools_4.2.1       glue_1.7.0        parallel_4.2.1    compiler_4.2.1    pkgconfig_2.0.3   tidyselect_1.2.1 
# [23] tibble_3.2.1   

library(survival)
library(dplyr)

setwd("/code/3_Comparison_between_clinical_markers_and_TIMES/")

#######load clinic data####
load("/data/Input_3_Comparison_between_clinical_markers_and_TIMES/clinic_data_input.RData")
##select data that have clinic measurements.
idx_c4 <- which(all_data$`Cohort No.`==4 & all_data$`Slide ID`=="NA")
idx_all <- c(which(all_data$`Cohort No.` %in% c(1,2,3,5)), idx_c4)
all_data <- all_data[idx_all,]

idx_notna <- which(!is.na(all_data$TIMES))
table(all_data$Recurrence[idx_notna])
# N NA  Y 
# 91 68 72 
dim(all_data)
#254 138

all_data$Recurrence <-factor(all_data$Recurrence,levels = c("N","Y"))
colnames(all_data)[which(colnames(all_data)=="Patient ID")] <- "patient_id"
colnames(all_data)[which(colnames(all_data)=="DFS (Months)")] <- "DFS"

#TIMES
summary(all_data$TIMES)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.02612 0.15755 0.40359 0.40821 0.60635 0.94822      23 

data_TIMES  <- all_data %>%
  dplyr::select(DFS, Recurrence, TIMES, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_TIMES, id = patient_id)
res_out <- summary(res.cox)
res_out

cox_res_df <- cbind(cbind(row.names(res_out$cmap)[2],res_out$coefficients), res_out$conf.int)  ##start to store results.
cox_res_df <- as.data.frame(cox_res_df,stringsAsFactors = F); colnames(cox_res_df)[1] <- "X_name"
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_TIMES, 
#         id = patient_id)
# 
# n= 127, number of events= 55 
# (127 observations deleted due to missingness)
# 
# 1:2        coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES  4.4776   88.0275   0.5918 7.566 3.85e-14
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES     88.03    0.01136      27.6     280.8
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.814  (se = 0.034 )
# Likelihood ratio test= 67.63  on 1 df,   p=<2e-16
# Wald test            = 57.25  on 1 df,   p=4e-14
# Score (logrank) test = 71.2  on 1 df,   p=<2e-16

mylogit <- glm(Recurrence ~ TIMES, data = data_TIMES,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_TIMES)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.2341  -0.7240  -0.4267   0.6902   2.2317  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.6909     0.4338  -6.204 5.51e-10 ***
#   TIMES         5.4344     0.8352   6.507 7.68e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 223.75  on 162  degrees of freedom
# Residual deviance: 161.57  on 161  degrees of freedom
# (91 observations deleted due to missingness)
# AIC: 165.57
# 
# Number of Fisher Scoring iterations: 4


#######Stage 1: Univariate regression####
#####TNM staging#####
table(all_data$`TNM staging`)
# I          IA          IB          II        IIIA        IIIB         IVA         IVB   mpT1bNxMx        mpT3 mypT3 Nx Mx          NA         pT1        pT1b  pT1b Nx Mx         pT3      pT3 pN 
# 3           8          30          73           7           4           3           8           1           1           1          90           3           7           3           3           1 
# pT3a         pT4   pT4 Nx Mx      pT4 pN       rmpT4       ypT1b 
# 1           3           1           1           1           1 
all_data$TNM_2level <- NA
all_data$TNM_2level[grep('II',all_data$`TNM staging`)]<-'≥II'
all_data$TNM_2level[grep('IV',all_data$`TNM staging`)]<-'≥II'
all_data$TNM_2level[grep('T[3|4]',all_data$`TNM staging`)]<-'≥II'
all_data$TNM_2level[grep('^I[B|A]',all_data$`TNM staging`)]<-'<II'
all_data$TNM_2level[grep('^I$',all_data$`TNM staging`)]<-'<II'
all_data$TNM_2level[grep('T1',all_data$`TNM staging`)]<-'<II'
all_data$TNM_2level <- factor(all_data$TNM_2level,levels = c("<II","≥II"))
table(all_data$TNM_2level)
# <II ≥II 
# 56 108 

data_TNM  <- all_data %>%
  dplyr::select(DFS, Recurrence, TNM_2level, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TNM_2level, data = data_TNM, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[2],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TNM_2level, data = data_TNM, 
#         id = patient_id)
# 
# n= 104, number of events= 51 
# (147 observations deleted due to missingness)
# 
# 1:2               coef exp(coef) se(coef)     z Pr(>|z|)
# TNM_2level≥II 0.6941    2.0018   0.3225 2.152   0.0314
# 
# 1:2             exp(coef) exp(-coef) lower .95 upper .95
# TNM_2level≥II     2.002     0.4995     1.064     3.766
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.589  (se = 0.035 )
# Likelihood ratio test= 5.02  on 1 df,   p=0.03
# Wald test            = 4.63  on 1 df,   p=0.03
# Score (logrank) test = 4.79  on 1 df,   p=0.03
mylogit <- glm(Recurrence ~ TNM_2level, data = data_TNM,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TNM_2level, family = binomial, data = data_TNM)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -1.159  -1.159  -1.028   1.196   1.335  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)    -0.3629     0.3255  -1.115    0.265
# TNM_2level≥II   0.3203     0.3854   0.831    0.406
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 183.77  on 132  degrees of freedom
# Residual deviance: 183.07  on 131  degrees of freedom
# (121 observations deleted due to missingness)
# AIC: 187.07
# 
# Number of Fisher Scoring iterations: 4


#####BCLC staging#####
#BCLC
all_data$`BCLC staging` <-factor(all_data$`BCLC staging`,levels = c("A","B","C"))
table(all_data$`BCLC staging`)
# A  B  C 
# 26 19 85 

#BCLC 2 level: A vs B+C
all_data$BCLC_A_vs_BC <- NA
all_data$BCLC_A_vs_BC[grep('^[A]',all_data$`BCLC staging`)]<-'low'
all_data$BCLC_A_vs_BC[grep('^[B|C]',all_data$`BCLC staging`)]<-'high'
all_data$BCLC_A_vs_BC <-factor(all_data$BCLC_A_vs_BC,levels = c("low","high"))
table(all_data$BCLC_A_vs_BC)
# low high 
# 26  104

data_BCLC_A_vs_BC  <- all_data %>%
  dplyr::select(DFS, Recurrence, BCLC_A_vs_BC, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ BCLC_A_vs_BC, data = data_BCLC_A_vs_BC, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[2],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ BCLC_A_vs_BC, data = data_BCLC_A_vs_BC, 
#         id = patient_id)
# 
# n= 100, number of events= 48 
# (151 observations deleted due to missingness)
# 
# 1:2                  coef exp(coef) se(coef)     z Pr(>|z|)
# BCLC_A_vs_BChigh 1.1473    3.1497   0.4424 2.593  0.00951
# 
# 1:2                exp(coef) exp(-coef) lower .95 upper .95
# BCLC_A_vs_BChigh      3.15     0.3175     1.323     7.496
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.614  (se = 0.023 )
# Likelihood ratio test= 8.69  on 1 df,   p=0.003
# Wald test            = 6.73  on 1 df,   p=0.01
# Score (logrank) test = 7.44  on 1 df,   p=0.006
mylogit <- glm(Recurrence ~ BCLC_A_vs_BC, data = data_BCLC_A_vs_BC,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ BCLC_A_vs_BC, family = binomial, data = data_BCLC_A_vs_BC)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.1611  -1.1611  -0.9218   1.1938   1.4566  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)       -0.6360     0.4122  -1.543    0.123
# BCLC_A_vs_BChigh   0.5975     0.4565   1.309    0.191
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 179.45  on 129  degrees of freedom
# Residual deviance: 177.68  on 128  degrees of freedom
# (124 observations deleted due to missingness)
# AIC: 181.68
# 
# Number of Fisher Scoring iterations: 4


#MVI (Microvascular invasion) 
colnames(all_data)[which(colnames(all_data)=="Microvascular invasion")] <- "MVI"
all_data$MVI <- factor(all_data$MVI,levels = c("M0","M1","M2"))
table(all_data$MVI)
# M0 M1 M2 
# 35 45 30

#MVI 2 level: M0 vs M1+M2
all_data$MVI_M0_vs_M1M2 <- NA
all_data$MVI_M0_vs_M1M2[grep('^M0',all_data$MVI)]<-'no'
all_data$MVI_M0_vs_M1M2[grep('^M[1|2]',all_data$MVI)]<-'yes'
all_data$MVI_M0_vs_M1M2 <-factor(all_data$MVI_M0_vs_M1M2,levels = c("no","yes"))
table(all_data$MVI_M0_vs_M1M2)
# no yes 
# 35  75 

data_MVI_M0_vs_M1M2  <- all_data %>%
  dplyr::select(DFS, Recurrence, MVI_M0_vs_M1M2, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ MVI_M0_vs_M1M2, data = data_MVI_M0_vs_M1M2, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[2],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ MVI_M0_vs_M1M2, data = data_MVI_M0_vs_M1M2, 
#         id = patient_id)
# 
# n= 86, number of events= 39 
# (168 observations deleted due to missingness)
# 
# 1:2                   coef exp(coef) se(coef)     z Pr(>|z|)
# MVI_M0_vs_M1M2yes 1.1338    3.1074   0.4061 2.792  0.00524
# 
# 1:2                 exp(coef) exp(-coef) lower .95 upper .95
# MVI_M0_vs_M1M2yes     3.107     0.3218     1.402     6.887
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.625  (se = 0.036 )
# Likelihood ratio test= 9.28  on 1 df,   p=0.002
# Wald test            = 7.8  on 1 df,   p=0.005
# Score (logrank) test = 8.55  on 1 df,   p=0.003
mylogit <- glm(Recurrence ~ MVI_M0_vs_M1M2, data = data_MVI_M0_vs_M1M2,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ MVI_M0_vs_M1M2, family = binomial, 
#       data = data_MVI_M0_vs_M1M2)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.1888  -1.1888  -0.9164   1.1661   1.4632  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)        -0.6506     0.3561  -1.827   0.0677 .
# MVI_M0_vs_M1M2yes   0.6773     0.4244   1.596   0.1106  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 151.58  on 109  degrees of freedom
# Residual deviance: 148.96  on 108  degrees of freedom
# (144 observations deleted due to missingness)
# AIC: 152.96
# 
# Number of Fisher Scoring iterations: 4


#macrovascular invasion
colnames(all_data)[which(colnames(all_data)=="Macrovascular invasion")] <- "macrovascular_invasion"
all_data$macrovascular_invasion <- factor(all_data$macrovascular_invasion,levels = c("no","yes"))
table(all_data$macrovascular_invasion)
# no yes 
# 124  19 
data_macrovascular_invasion  <- all_data %>%
  dplyr::select(DFS, Recurrence, macrovascular_invasion, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ macrovascular_invasion, data = data_macrovascular_invasion, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[2],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ macrovascular_invasion, 
#         data = data_macrovascular_invasion, id = patient_id)
# 
# n= 108, number of events= 53 
# (146 observations deleted due to missingness)
# 
# 1:2                           coef exp(coef) se(coef)     z Pr(>|z|)
# macrovascular_invasionyes 1.2872    3.6227   0.3333 3.862 0.000113
# 
# 1:2                         exp(coef) exp(-coef) lower .95 upper .95
# macrovascular_invasionyes     3.623      0.276     1.885     6.963
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.596  (se = 0.031 )
# Likelihood ratio test= 12.15  on 1 df,   p=5e-04
# Wald test            = 14.91  on 1 df,   p=1e-04
# Score (logrank) test = 17.04  on 1 df,   p=4e-05

mylogit <- glm(Recurrence ~ macrovascular_invasion, data = data_macrovascular_invasion,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ macrovascular_invasion, family = binomial, 
#       data = data_macrovascular_invasion)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -1.634  -1.083  -1.083   1.275   1.275  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)                -0.2268     0.1808  -1.255   0.2096  
# macrovascular_invasionyes   1.2564     0.5515   2.278   0.0227 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 198.07  on 142  degrees of freedom
# Residual deviance: 192.22  on 141  degrees of freedom
# (111 observations deleted due to missingness)
# AIC: 196.22
# 
# Number of Fisher Scoring iterations: 4


##multivariate regression analyses for vascular invasion
data_vascular_invasion  <- all_data %>%
  dplyr::select(DFS, Recurrence, MVI_M0_vs_M1M2, macrovascular_invasion, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES + MVI_M0_vs_M1M2 + macrovascular_invasion, data = data_vascular_invasion, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES + MVI_M0_vs_M1M2 + 
#           macrovascular_invasion, data = data_vascular_invasion, id = patient_id)
# 
# n= 86, number of events= 39 
# (168 observations deleted due to missingness)
# 
# 1:2                              coef exp(coef)  se(coef)     z Pr(>|z|)
# TIMES                       5.66239 287.83694   0.91875 6.163 7.13e-10
# MVI_M0_vs_M1M2yes           0.85699   2.35607   0.44032 1.946   0.0516
# macrovascular_invasionyes   0.08257   1.08607   0.39641 0.208   0.8350
# 
# 1:2                         exp(coef) exp(-coef) lower .95 upper .95
# TIMES                       287.837   0.003474   47.5447  1742.574
# MVI_M0_vs_M1M2yes             2.356   0.424436    0.9940     5.585
# macrovascular_invasionyes     1.086   0.920749    0.4994     2.362
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.85  (se = 0.028 )
# Likelihood ratio test= 68.65  on 3 df,   p=8e-15
# Wald test            = 49.85  on 3 df,   p=9e-11
# Score (logrank) test = 75.22  on 3 df,   p=3e-16

mylogit <- glm(Recurrence ~ TIMES + MVI_M0_vs_M1M2 + macrovascular_invasion, data = data_vascular_invasion,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES + MVI_M0_vs_M1M2 + macrovascular_invasion, 
#       family = binomial, data = data_vascular_invasion)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.3920  -0.5558  -0.3358   0.5141   2.4941  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                -3.4168     0.6901  -4.951 7.38e-07 ***
#   TIMES                       6.6624     1.1773   5.659 1.52e-08 ***
#   MVI_M0_vs_M1M2yes           0.3085     0.5664   0.545    0.586    
# macrovascular_invasionyes   0.2908     0.8879   0.328    0.743    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 151.582  on 109  degrees of freedom
# Residual deviance:  93.191  on 106  degrees of freedom
# (144 observations deleted due to missingness)
# AIC: 101.19
# 
# Number of Fisher Scoring iterations: 5


####etiology of HCC (Hep B, C, etc), presence of liver fibrosis, cirrhosis. #####
####(1-1) abnormality of hepatitis B virus surface antigen quantification ####
colnames(all_data)[which(colnames(all_data)=="Abnormality of hepatitis B virus surface antigen quantification")] <- "HPB_surface_antigen"
all_data$HPB_surface_antigen <- factor(all_data$HPB_surface_antigen,levels = c("-","+"))
table(all_data$HPB_surface_antigen)
# -   + 
#   26 109  
# 109/(109+26)
# [1] 0.8074074
data_HPB_surface_antigen  <- all_data %>%
  dplyr::select(DFS, Recurrence, HPB_surface_antigen, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ HPB_surface_antigen, data = data_HPB_surface_antigen, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[2],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ HPB_surface_antigen, 
#         data = data_HPB_surface_antigen, id = patient_id)
# 
# n= 104, number of events= 52 
# (150 observations deleted due to missingness)
# 
# 1:2                       coef exp(coef) se(coef)      z Pr(>|z|)
# HPB_surface_antigen+ -0.4353    0.6471   0.3437 -1.267    0.205
# 
# 1:2                    exp(coef) exp(-coef) lower .95 upper .95
# HPB_surface_antigen+    0.6471      1.545    0.3299     1.269
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.543  (se = 0.031 )
# Likelihood ratio test= 1.48  on 1 df,   p=0.2
# Wald test            = 1.6  on 1 df,   p=0.2
# Score (logrank) test = 1.63  on 1 df,   p=0.2

mylogit <- glm(Recurrence ~ HPB_surface_antigen, data = data_HPB_surface_antigen,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ HPB_surface_antigen, family = binomial, 
#       data = data_HPB_surface_antigen)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -1.244  -1.123  -1.123   1.232   1.232  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)            0.1542     0.3934   0.392    0.695
# HPB_surface_antigen+  -0.2828     0.4377  -0.646    0.518
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 186.96  on 134  degrees of freedom
# Residual deviance: 186.55  on 133  degrees of freedom
# (119 observations deleted due to missingness)
# AIC: 190.55
# 
# Number of Fisher Scoring iterations: 3


####(1-2) abnormality of hepatitis B virus surface antibody quantification ####
colnames(all_data)[which(colnames(all_data)=="Abnormality of hepatitis B virus surface antibody quantification")] <- "HPB_surface_antibody"
all_data$HPB_surface_antibody <- factor(all_data$HPB_surface_antibody,levels = c("-","+"))
table(all_data$HPB_surface_antibody)
# -   + 
#   118  18 

data_HPB_surface_antibody  <- all_data %>%
  dplyr::select(DFS, Recurrence, HPB_surface_antibody, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ HPB_surface_antibody, data = data_HPB_surface_antibody, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[2],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ HPB_surface_antibody, 
#         data = data_HPB_surface_antibody, id = patient_id)
# 
# n= 105, number of events= 52 
# (149 observations deleted due to missingness)
# 
# 1:2                        coef exp(coef) se(coef)     z Pr(>|z|)
# HPB_surface_antibody+ 0.09357   1.09809  0.40823 0.229    0.819
# 
# 1:2                     exp(coef) exp(-coef) lower .95 upper .95
# HPB_surface_antibody+     1.098     0.9107    0.4933     2.444
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.504  (se = 0.027 )
# Likelihood ratio test= 0.05  on 1 df,   p=0.8
# Wald test            = 0.05  on 1 df,   p=0.8
# Score (logrank) test = 0.05  on 1 df,   p=0.8

mylogit <- glm(Recurrence ~ HPB_surface_antibody, data = data_HPB_surface_antibody,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ HPB_surface_antibody, family = binomial, 
#       data = data_HPB_surface_antibody)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -1.177  -1.135  -1.135   1.221   1.221  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)            -0.1018     0.1844  -0.552    0.581
# HPB_surface_antibody+   0.1018     0.5062   0.201    0.841
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 188.27  on 135  degrees of freedom
# Residual deviance: 188.23  on 134  degrees of freedom
# (118 observations deleted due to missingness)
# AIC: 192.23
# 
# Number of Fisher Scoring iterations: 3


####(1-3) abnormality of hepatitis B virus e antigen quantification #### 
colnames(all_data)[which(colnames(all_data)=="Abnormality of hepatitis B virus e antigen quantification")] <- "HPB_e_antigen"
all_data$HPB_e_antigen <- factor(all_data$HPB_e_antigen,levels = c("-","+"))
table(all_data$HPB_e_antigen)
# -   + 
#   114  21 

data_HPB_e_antigen  <- all_data %>%
  dplyr::select(DFS, Recurrence, HPB_e_antigen, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ HPB_e_antigen, data = data_HPB_e_antigen, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[2],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ HPB_e_antigen, data = data_HPB_e_antigen, 
#         id = patient_id)
# 
# n= 104, number of events= 52 
# (150 observations deleted due to missingness)
# 
# 1:2                 coef exp(coef) se(coef)     z Pr(>|z|)
# HPB_e_antigen+ 0.08725   1.09117  0.38528 0.226    0.821
# 
# 1:2              exp(coef) exp(-coef) lower .95 upper .95
# HPB_e_antigen+     1.091     0.9165    0.5128     2.322
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.515  (se = 0.029 )
# Likelihood ratio test= 0.05  on 1 df,   p=0.8
# Wald test            = 0.05  on 1 df,   p=0.8
# Score (logrank) test = 0.05  on 1 df,   p=0.8

mylogit <- glm(Recurrence ~ HPB_e_antigen, data = data_HPB_e_antigen,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ HPB_e_antigen, family = binomial, 
#       data = data_HPB_e_antigen)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -1.148  -1.148  -1.137   1.207   1.218  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)    -0.07020    0.18743  -0.375    0.708
# HPB_e_antigen+ -0.02511    0.47544  -0.053    0.958
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 186.96  on 134  degrees of freedom
# Residual deviance: 186.96  on 133  degrees of freedom
# (119 observations deleted due to missingness)
# AIC: 190.96
# 
# Number of Fisher Scoring iterations: 3


####(1-4) abnormality of hepatitis B virus e antibody quantification  (low, normal as negative, high as positive) ####
colnames(all_data)[which(colnames(all_data)=="Abnormality of hepatitis B virus e antibody quantification")] <- "HPB_e_antibody"
all_data$HPB_e_antibody <- factor(all_data$HPB_e_antibody,levels = c("-","+"))
table(all_data$HPB_e_antibody)
# -   + 
#   119  17   

data_HPB_e_antibody  <- all_data %>%
  dplyr::select(DFS, Recurrence, HPB_e_antibody, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ HPB_e_antibody, data = data_HPB_e_antibody, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[2],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ HPB_e_antibody, data = data_HPB_e_antibody, 
#         id = patient_id)
# 
# n= 105, number of events= 52 
# (149 observations deleted due to missingness)
# 
# 1:2                   coef exp(coef) se(coef)      z Pr(>|z|)
# HPB_e_antibody+ -0.01442   0.98568  0.38565 -0.037     0.97
# 
# 1:2               exp(coef) exp(-coef) lower .95 upper .95
# HPB_e_antibody+    0.9857      1.015    0.4629     2.099
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.523  (se = 0.022 )
# Likelihood ratio test= 0  on 1 df,   p=1
# Wald test            = 0  on 1 df,   p=1
# Score (logrank) test = 0  on 1 df,   p=1

mylogit <- glm(Recurrence ~ HPB_e_antibody, data = data_HPB_e_antibody,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ HPB_e_antibody, family = binomial, 
#       data = data_HPB_e_antibody)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -1.332  -1.114  -1.114   1.242   1.242  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)      -0.1515     0.1839  -0.824    0.410
# HPB_e_antibody+   0.5082     0.5260   0.966    0.334
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 188.27  on 135  degrees of freedom
# Residual deviance: 187.32  on 134  degrees of freedom
# (118 observations deleted due to missingness)
# AIC: 191.32
# 
# Number of Fisher Scoring iterations: 3


####(1-5)  abnormality of hepatitis B virus core antibody quantitation  ####
colnames(all_data)[which(colnames(all_data)=="Abnormality of hepatitis B virus core antibody quantitation")] <- "HPB_core_antibody"
all_data$HPB_core_antibody <- factor(all_data$HPB_core_antibody,levels = c("-","+"))
table(all_data$HPB_core_antibody)
# -  + 
#   84 37  

data_HPB_core_antibody  <- all_data %>%
  dplyr::select(DFS, Recurrence, HPB_core_antibody, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ HPB_core_antibody, data = data_HPB_core_antibody, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[2],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ HPB_core_antibody, data = data_HPB_core_antibody, 
#         id = patient_id)
# 
# n= 90, number of events= 47 
# (164 observations deleted due to missingness)
# 
# 1:2                       coef exp(coef)  se(coef) z Pr(>|z|)
# HPB_core_antibody+ 8.513e-05 1.000e+00 3.214e-01 0        1
# 
# 1:2                  exp(coef) exp(-coef) lower .95 upper .95
# HPB_core_antibody+         1     0.9999    0.5327     1.877
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.491  (se = 0.036 )
# Likelihood ratio test= 0  on 1 df,   p=1
# Wald test            = 0  on 1 df,   p=1
# Score (logrank) test = 0  on 1 df,   p=1

mylogit <- glm(Recurrence ~ HPB_core_antibody, data = data_HPB_core_antibody,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ HPB_core_antibody, family = binomial, 
#       data = data_HPB_core_antibody)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -1.218  -1.218  -1.064   1.137   1.295  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)         0.09531    0.21847   0.436    0.663
# HPB_core_antibody+ -0.36724    0.39729  -0.924    0.355
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 167.73  on 120  degrees of freedom
# Residual deviance: 166.87  on 119  degrees of freedom
# (133 observations deleted due to missingness)
# AIC: 170.87
# 
# Number of Fisher Scoring iterations: 3


####(2) abnormality of hepatitis C antibody-IgG ####
colnames(all_data)[which(colnames(all_data)=="Abnormality of hepatitis C antibody-IgG")] <- "HPC_antibody"
all_data$HPC_antibody <- factor(all_data$HPC_antibody,levels = c("-","+"))
table(all_data$HPC_antibody)
# -   + 
#   127   4  
# > 127/(127+4)
# [1] 0.9694656
data_HPC_antibody  <- all_data %>%
  dplyr::select(DFS, Recurrence, HPC_antibody, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ HPC_antibody, data = data_HPC_antibody, id = patient_id)
# Warning message:
#   In fitter(X, Y, istrat, offset, init, control, weights = weights,  :
#               Loglik converged before variable  1 ; coefficient may be infinite. 
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[2],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ HPC_antibody, data = data_HPC_antibody, 
#         id = patient_id)
# 
# n= 104, number of events= 51 
# (150 observations deleted due to missingness)
# 
# 1:2                   coef  exp(coef)   se(coef)      z Pr(>|z|)
# HPC_antibody+ -1.706e+01  3.915e-08  4.462e+03 -0.004    0.997
# 
# 1:2             exp(coef) exp(-coef) lower .95 upper .95
# HPC_antibody+ 3.915e-08   25541780         0       Inf
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.511  (se = 0.008 )
# Likelihood ratio test= 2.52  on 1 df,   p=0.1
# Wald test            = 0  on 1 df,   p=1
# Score (logrank) test = 1.28  on 1 df,   p=0.3

mylogit <- glm(Recurrence ~ HPC_antibody, data = data_HPC_antibody,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ HPC_antibody, family = binomial, data = data_HPC_antibody)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.1574  -1.1574  -0.7585   1.1975   1.6651  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   -0.04725    0.17752  -0.266    0.790
# HPC_antibody+ -1.05136    1.16827  -0.900    0.368
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 181.41  on 130  degrees of freedom
# Residual deviance: 180.49  on 129  degrees of freedom
# (123 observations deleted due to missingness)
# AIC: 184.49
# 
# Number of Fisher Scoring iterations: 4


####(3) liver cirrhosis####
colnames(all_data)[which(colnames(all_data)=="Liver cirrhosis")] <- "liver_cirrhosis"
all_data$liver_cirrhosis <- factor(all_data$liver_cirrhosis,levels = c("-","+"))
table(all_data$liver_cirrhosis)
# -   + 
#   28 111 
# > 111/(111+28)
# [1] 0.7985612
data_liver_cirrhosis  <- all_data %>%
  dplyr::select(DFS, Recurrence, liver_cirrhosis, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ liver_cirrhosis, data = data_liver_cirrhosis, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[2],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ liver_cirrhosis, data = data_liver_cirrhosis, 
#         id = patient_id)
# 
# n= 107, number of events= 53 
# (147 observations deleted due to missingness)
# 
# 1:2                  coef exp(coef) se(coef)     z Pr(>|z|)
# liver_cirrhosis+ 0.6051    1.8315   0.4066 1.488    0.137
# 
# 1:2                exp(coef) exp(-coef) lower .95 upper .95
# liver_cirrhosis+     1.832      0.546    0.8255     4.064
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.535  (se = 0.03 )
# Likelihood ratio test= 2.56  on 1 df,   p=0.1
# Wald test            = 2.21  on 1 df,   p=0.1
# Score (logrank) test = 2.28  on 1 df,   p=0.1

mylogit <- glm(Recurrence ~ liver_cirrhosis, data = data_liver_cirrhosis,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ liver_cirrhosis, family = binomial, 
#       data = data_liver_cirrhosis)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.2472  -1.2472  -0.8203   1.1092   1.5829  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)       -0.9163     0.4183  -2.190   0.0285 *
#   liver_cirrhosis+   1.0788     0.4596   2.347   0.0189 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 192.63  on 138  degrees of freedom
# Residual deviance: 186.65  on 137  degrees of freedom
# (115 observations deleted due to missingness)
# AIC: 190.65
# 
# Number of Fisher Scoring iterations: 4


#######other factors for univariate regression analysis####

#### group 1: pathology  ####
                                                                         
### tumor differentiation grade: 3 levels ###                                          
table(all_data$`Tumor differentiation grade`)
# G1 G1-G2    G2 G2-G3    G3    G4 
# 2     5    68    21    32     1

all_data$tumor_differentiation_grade <- NA
all_data$tumor_differentiation_grade[grep('G1',all_data$`Tumor differentiation grade`)] <-'≤G2'
all_data$tumor_differentiation_grade[grep('^G2$',all_data$`Tumor differentiation grade`)]<-'≤G2'
all_data$tumor_differentiation_grade[grep('^G3$',all_data$`Tumor differentiation grade`)]<-'≥G3'
all_data$tumor_differentiation_grade[grep('^G4$',all_data$`Tumor differentiation grade`)]<-'≥G3'
all_data$tumor_differentiation_grade[grep('G2-G3',all_data$`Tumor differentiation grade`)]<-'G2-G3'
all_data$tumor_differentiation_grade <- factor(all_data$tumor_differentiation_grade,levels = c("≤G2","G2-G3","≥G3"))
table(all_data$tumor_differentiation_grade)
# ≤G2 G2-G3   ≥G3 
# 75    21    33 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, tumor_differentiation_grade, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ tumor_differentiation_grade, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ tumor_differentiation_grade, 
#         data = data_factor, id = patient_id)
# 
# n= 98, number of events= 48 
# (156 observations deleted due to missingness)
# 
# 1:2                                    coef exp(coef) se(coef)      z Pr(>|z|)
# tumor_differentiation_gradeG2-G3  1.39819   4.04786  0.36399  3.841 0.000122
# tumor_differentiation_grade≥G3   -0.03058   0.96988  0.36450 -0.084 0.933144
# 
# 1:2                                exp(coef) exp(-coef) lower .95 upper .95
# tumor_differentiation_gradeG2-G3    4.0479      0.247    1.9833     8.261
# tumor_differentiation_grade≥G3      0.9699      1.031    0.4747     1.981
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.604  (se = 0.041 )
# Likelihood ratio test= 13.16  on 2 df,   p=0.001
# Wald test            = 16.55  on 2 df,   p=3e-04
# Score (logrank) test = 19.36  on 2 df,   p=6e-05
mylogit <- glm(Recurrence ~ tumor_differentiation_grade, data = data_factor,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ tumor_differentiation_grade, family = binomial, 
#       data = data_factor)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -1.482  -1.055  -1.055   1.203   1.305  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)                       -0.2955     0.2335  -1.266   0.2057  
# tumor_differentiation_gradeG2-G3   0.9886     0.5185   1.907   0.0565 .
# tumor_differentiation_grade≥G3     0.2348     0.4193   0.560   0.5754  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 178.64  on 128  degrees of freedom
# Residual deviance: 174.80  on 126  degrees of freedom
# (125 observations deleted due to missingness)
# AIC: 180.8
# 
# Number of Fisher Scoring iterations: 4


### tumor greatest dimension (cm) ###
colnames(all_data)[which(colnames(all_data)=="Tumor greatest dimension (cm)")] <- "tumor_greatest_dimension"
summary(all_data$tumor_greatest_dimension)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.600   3.725   5.000   6.180   8.000  19.000      96 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, tumor_greatest_dimension, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ tumor_greatest_dimension, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ tumor_greatest_dimension, 
#         data = data_factor, id = patient_id)
# 
# n= 100, number of events= 49 
# (154 observations deleted due to missingness)
# 
# 1:2                           coef exp(coef) se(coef)     z Pr(>|z|)
# tumor_greatest_dimension 0.19282   1.21267  0.03318 5.812 6.17e-09
# 
# 1:2                        exp(coef) exp(-coef) lower .95 upper .95
# tumor_greatest_dimension     1.213     0.8246     1.136     1.294
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.748  (se = 0.036 )
# Likelihood ratio test= 27.74  on 1 df,   p=1e-07
# Wald test            = 33.78  on 1 df,   p=6e-09
# Score (logrank) test = 36.65  on 1 df,   p=1e-09

mylogit <- glm(Recurrence ~ tumor_greatest_dimension, data = data_factor,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ tumor_greatest_dimension, family = binomial, 
#       data = data_factor)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.9229  -0.9689  -0.7503   1.1264   1.7564  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)               -1.4351     0.3848  -3.730 0.000192 ***
#   tumor_greatest_dimension   0.2223     0.0586   3.794 0.000148 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 179.73  on 129  degrees of freedom
# Residual deviance: 161.99  on 128  degrees of freedom
# (124 observations deleted due to missingness)
# AIC: 165.99
# 
# Number of Fisher Scoring iterations: 4


### tumor burden (cm3)###
colnames(all_data)[which(colnames(all_data)=="Tumor burden (cm3)")] <- "tumor_burden"
summary(all_data$tumor_burden)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#   0.0377  10.4077  34.5180 121.1640 131.7960 883.8700      127 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, tumor_burden, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ tumor_burden, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ tumor_burden, data = data_factor, 
#         id = patient_id)
# 
# n= 98, number of events= 49 
# (156 observations deleted due to missingness)
# 
# 1:2                 coef exp(coef)  se(coef)    z Pr(>|z|)
# tumor_burden 0.0027544 1.0027581 0.0005786 4.76 1.94e-06
# 
# 1:2            exp(coef) exp(-coef) lower .95 upper .95
# tumor_burden     1.003     0.9972     1.002     1.004
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.719  (se = 0.04 )
# Likelihood ratio test= 16.91  on 1 df,   p=4e-05
# Wald test            = 22.66  on 1 df,   p=2e-06
# Score (logrank) test = 26.48  on 1 df,   p=3e-07

mylogit <- glm(Recurrence ~ tumor_burden, data = data_factor,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ tumor_burden, family = binomial, data = data_factor)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.1847  -0.9997  -0.9737   1.2727   1.3969  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -0.502759   0.220574  -2.279  0.02265 * 
#   tumor_burden  0.003468   0.001222   2.837  0.00456 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 175.67  on 126  degrees of freedom
# Residual deviance: 165.27  on 125  degrees of freedom
# (127 observations deleted due to missingness)
# AIC: 169.27
# 
# Number of Fisher Scoring iterations: 4


### tumor area_size (cm2) ###
colnames(all_data)[which(colnames(all_data)=="Tumor area size (cm2)")] <- "tumor_area_size"
summary(all_data$tumor_area_size)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0060  0.9008  1.9875  3.7160  4.9500 19.2060     127 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, tumor_area_size, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ tumor_area_size, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ tumor_area_size, data = data_factor, 
#         id = patient_id)
# 
# n= 98, number of events= 49 
# (156 observations deleted due to missingness)
# 
# 1:2                  coef exp(coef) se(coef)     z Pr(>|z|)
# tumor_area_size 0.12903   1.13772  0.02753 4.687 2.78e-06
# 
# 1:2               exp(coef) exp(-coef) lower .95 upper .95
# tumor_area_size     1.138     0.8789     1.078     1.201
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.715  (se = 0.038 )
# Likelihood ratio test= 17.03  on 1 df,   p=4e-05
# Wald test            = 21.96  on 1 df,   p=3e-06
# Score (logrank) test = 24.24  on 1 df,   p=9e-07

mylogit <- glm(Recurrence ~ tumor_area_size, data = data_factor,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ tumor_area_size, family = binomial, 
#       data = data_factor)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.2091  -0.9963  -0.9148   1.2279   1.4738  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)     -0.67770    0.25471  -2.661  0.00780 **
#   tumor_area_size  0.15759    0.05206   3.027  0.00247 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 175.67  on 126  degrees of freedom
# Residual deviance: 164.72  on 125  degrees of freedom
# (127 observations deleted due to missingness)
# AIC: 168.72
# 
# Number of Fisher Scoring iterations: 4


### tumor count ###
colnames(all_data)[which(colnames(all_data)=="Tumor count")] <- "tumor_count"
summary(all_data$tumor_count)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.000   1.000   1.000   1.133   1.000   5.000     126  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, tumor_count, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ tumor_count, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ tumor_count, data = data_factor, 
#         id = patient_id)
# 
# n= 100, number of events= 49 
# (154 observations deleted due to missingness)
# 
# 1:2              coef exp(coef) se(coef)      z Pr(>|z|)
# tumor_count -0.4718    0.6239   0.4749 -0.993    0.321
# 
# 1:2           exp(coef) exp(-coef) lower .95 upper .95
# tumor_count    0.6239      1.603    0.2459     1.583
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.518  (se = 0.019 )
# Likelihood ratio test= 1.45  on 1 df,   p=0.2
# Wald test            = 0.99  on 1 df,   p=0.3
# Score (logrank) test = 1.04  on 1 df,   p=0.3

mylogit <- glm(Recurrence ~ tumor_count, data = data_factor,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ tumor_count, family = binomial, data = data_factor)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.1496  -1.1496  -0.9526   1.2054   1.4203  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.4235     0.5447   0.778    0.437
# tumor_count  -0.4892     0.4654  -1.051    0.293
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 176.95  on 127  degrees of freedom
# Residual deviance: 175.60  on 126  degrees of freedom
# (126 observations deleted due to missingness)
# AIC: 179.6
# 
# Number of Fisher Scoring iterations: 3


### capsular invasion ###                                          
colnames(all_data)[which(colnames(all_data)=="Capsular invasion")] <- "capsular_invasion"
all_data$capsular_invasion <- factor(all_data$capsular_invasion,levels = c("-","+"))
table(all_data$capsular_invasion)
# -  + 
#   79 47 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, capsular_invasion, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ capsular_invasion, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ capsular_invasion, data = data_factor, 
#         id = patient_id)
# 
# n= 99, number of events= 47 
# (155 observations deleted due to missingness)
# 
# 1:2                    coef exp(coef) se(coef)     z Pr(>|z|)
# capsular_invasion+ 0.6588    1.9324   0.2968 2.219   0.0265
# 
# 1:2                  exp(coef) exp(-coef) lower .95 upper .95
# capsular_invasion+     1.932     0.5175      1.08     3.457
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.588  (se = 0.04 )
# Likelihood ratio test= 4.93  on 1 df,   p=0.03
# Wald test            = 4.93  on 1 df,   p=0.03
# Score (logrank) test = 5.1  on 1 df,   p=0.02
mylogit <- glm(Recurrence ~ capsular_invasion, data = data_factor,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ capsular_invasion, family = binomial, 
#       data = data_factor)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.3459  -0.9982  -0.9982   1.0178   1.3678  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)         -0.4372     0.2304  -1.898   0.0578 .
# capsular_invasion+   0.8250     0.3761   2.194   0.0283 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 174.16  on 125  degrees of freedom
# Residual deviance: 169.25  on 124  degrees of freedom
# (128 observations deleted due to missingness)
# AIC: 173.25
# 
# Number of Fisher Scoring iterations: 4


### satellite lesion ###                                          
colnames(all_data)[which(colnames(all_data)=="Satellite lesion")] <- "satellite_lesion"
all_data$satellite_lesion <- factor(all_data$satellite_lesion,levels = c("-","+"))
table(all_data$satellite_lesion)
# -  + 
#   93 33 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, satellite_lesion, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ satellite_lesion, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ satellite_lesion, data = data_factor, 
#         id = patient_id)
# 
# n= 99, number of events= 47 
# (155 observations deleted due to missingness)
# 
# 1:2                   coef exp(coef) se(coef)    z Pr(>|z|)
# satellite_lesion+ 1.1173    3.0566   0.3095 3.61 0.000306
# 
# 1:2                 exp(coef) exp(-coef) lower .95 upper .95
# satellite_lesion+     3.057     0.3272     1.667     5.606
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.644  (se = 0.038 )
# Likelihood ratio test= 11.47  on 1 df,   p=7e-04
# Wald test            = 13.03  on 1 df,   p=3e-04
# Score (logrank) test = 14.39  on 1 df,   p=1e-04
mylogit <- glm(Recurrence ~ satellite_lesion, data = data_factor,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ satellite_lesion, family = binomial, 
#       data = data_factor)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -1.310  -1.060  -1.060   1.299   1.299  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)        -0.2814     0.2094  -1.344    0.179
# satellite_lesion+   0.5868     0.4098   1.432    0.152
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 174.16  on 125  degrees of freedom
# Residual deviance: 172.09  on 124  degrees of freedom
# (128 observations deleted due to missingness)
# AIC: 176.09
# 
# Number of Fisher Scoring iterations: 4


####Multiplicity of HCC tumor nodules also has a high positive predictive accuracy for recurrence after resection and should also be a gold standard for comparison and stratification.####
##multivariate regression analyses for tumor nodules
data_nodule  <- all_data %>%
  dplyr::select(DFS, Recurrence, satellite_lesion, tumor_count, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES + satellite_lesion + tumor_count, data = data_nodule, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES + satellite_lesion + 
#           tumor_count, data = data_nodule, id = patient_id)
# 
# n= 98, number of events= 47 
# (156 observations deleted due to missingness)
# 
# 1:2                     coef exp(coef) se(coef)      z Pr(>|z|)
# TIMES               6.6856  800.8219   0.9442  7.081 1.44e-12
# satellite_lesion+   1.2994    3.6670   0.3366  3.860 0.000114
# tumor_count        -0.8627    0.4220   0.5994 -1.439 0.150114
# 
# 1:2                 exp(coef) exp(-coef) lower .95 upper .95
# TIMES               800.822   0.001249  125.8382  5096.350
# satellite_lesion+     3.667   0.272702    1.8956     7.094
# tumor_count           0.422   2.369490    0.1303     1.366
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.879  (se = 0.026 )
# Likelihood ratio test= 94.34  on 3 df,   p=<2e-16
# Wald test            = 64.48  on 3 df,   p=6e-14
# Score (logrank) test = 93.4  on 3 df,   p=<2e-16

mylogit <- glm(Recurrence ~ TIMES + satellite_lesion + tumor_count, data = data_nodule,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES + satellite_lesion + tumor_count, 
#       family = binomial, data = data_nodule)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.7114  -0.5209  -0.2092   0.5280   2.4964  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        -2.6039     0.9887  -2.634  0.00844 ** 
#   TIMES               7.6471     1.2793   5.978 2.26e-09 ***
#   satellite_lesion+   0.6137     0.6163   0.996  0.31931    
# tumor_count        -1.1433     0.8082  -1.415  0.15721    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 170.116  on 122  degrees of freedom
# Residual deviance:  96.834  on 119  degrees of freedom
# (131 observations deleted due to missingness)
# AIC: 104.83
# 
# Number of Fisher Scoring iterations: 5


### neuroaggression ###                                          
colnames(all_data)[which(colnames(all_data)=="Neuroaggression")] <- "neuroaggression"
all_data$neuroaggression <- factor(all_data$neuroaggression,levels = c("-","+"))
table(all_data$neuroaggression)
# - + 
#   4 1

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, neuroaggression, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ neuroaggression, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
# Warning message:
#   In fitter(X, Y, istrat, offset, init, control, weights = weights,  :
#               Ran out of iterations and did not converge
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ neuroaggression, data = data_factor, 
#         id = patient_id)
# 
# n= 3, number of events= 2 
# (251 observations deleted due to missingness)
# 
# 1:2                     coef exp(coef)  se(coef) z Pr(>|z|)
# neuroaggression+ 2.215e+01 4.171e+09 4.566e+04 0        1
# 
# 1:2                exp(coef) exp(-coef) lower .95 upper .95
# neuroaggression+ 4.171e+09  2.398e-10         0       Inf
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 1  (se = 0 )
# Likelihood ratio test= 2.2  on 1 df,   p=0.1
# Wald test            = 0  on 1 df,   p=1
# Score (logrank) test = 2  on 1 df,   p=0.2


### neoplastic thrombus ###                                          
colnames(all_data)[which(colnames(all_data)=="Neoplastic thrombus")] <- "neoplastic_thrombus"
all_data$neoplastic_thrombus <- factor(all_data$neoplastic_thrombus,levels = c("-","+"))
table(all_data$neoplastic_thrombus)
# -   + 
#   108  19   

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, neoplastic_thrombus, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ neoplastic_thrombus, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ neoplastic_thrombus, 
#         data = data_factor, id = patient_id)
# 
# n= 99, number of events= 48 
# (155 observations deleted due to missingness)
# 
# 1:2                      coef exp(coef) se(coef)     z Pr(>|z|)
# neoplastic_thrombus+ 0.7426    2.1014   0.3598 2.064    0.039
# 
# 1:2                    exp(coef) exp(-coef) lower .95 upper .95
# neoplastic_thrombus+     2.101     0.4759     1.038     4.254
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.56  (se = 0.031 )
# Likelihood ratio test= 3.69  on 1 df,   p=0.05
# Wald test            = 4.26  on 1 df,   p=0.04
# Score (logrank) test = 4.46  on 1 df,   p=0.03
mylogit <- glm(Recurrence ~ neoplastic_thrombus, data = data_factor,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ neoplastic_thrombus, family = binomial, 
#       data = data_factor)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -1.315  -1.100  -1.100   1.257   1.257  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)           -0.1857     0.1933  -0.961    0.337
# neoplastic_thrombus+   0.5042     0.5032   1.002    0.316
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 175.67  on 126  degrees of freedom
# Residual deviance: 174.66  on 125  degrees of freedom
# (127 observations deleted due to missingness)
# AIC: 178.66
# 
# Number of Fisher Scoring iterations: 3


### lymph node metastasis ###                                          
colnames(all_data)[which(colnames(all_data)=="Lymph node metastasis")] <- "lymph_node_metastasis"
all_data$lymph_node_metastasis <- factor(all_data$lymph_node_metastasis,levels = c("-","+"))
table(all_data$lymph_node_metastasis)
# - + 
#   1 2 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, lymph_node_metastasis, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ lymph_node_metastasis, data = data_factor, id = patient_id)
# Error in coxph(Surv(DFS, Recurrence) ~ lymph_node_metastasis, data = data_factor,  : 
#                  No (non-missing) observations
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
res_out


### surgical margin (cancer found?) ###                                          
colnames(all_data)[which(colnames(all_data)=="Surgical margin (cancer found?)")] <- "surgical_margin"
all_data$surgical_margin <- factor(all_data$surgical_margin,levels = c("no","yes"))
table(all_data$surgical_margin)
# no yes 
# 121   0 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, surgical_margin, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ surgical_margin, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ surgical_margin, data = data_factor, 
#         id = patient_id)
# 
# n= 98, number of events= 47 
# (156 observations deleted due to missingness)
# 
# 1:2                  coef exp(coef) se(coef)  z Pr(>|z|)
# surgical_marginyes   NA        NA        0 NA       NA
# 
# 1:2                  exp(coef) exp(-coef) lower .95 upper .95
# surgical_marginyes        NA         NA        NA        NA
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.5  (se = 0 )
# Likelihood ratio test= 0  on 0 df,   p=1
# Wald test            = NA  on 0 df,   p=NA
# Score (logrank) test = 0  on 0 df,   p=1


###group 2: IHC

### AACT ###                                          
all_data$AACT <- factor(all_data$AACT,levels = c("-","+"))
table(all_data$AACT)
# - + 
#   0 1

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, AACT, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ AACT, data = data_factor, id = patient_id)
# Error in dimnames(events) <- list(count = tab1.levels, state = c(ystate,  : 
#                                                                    'dimnames' applied to non-array
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
res_out


### AFP ###                                          
all_data$AFP <- factor(all_data$AFP,levels = c("-","+"))
table(all_data$AFP)
# -  + 
#   59 40 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, AFP, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ AFP, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ AFP, data = data_factor, 
#         id = patient_id)
# 
# n= 76, number of events= 41 
# (178 observations deleted due to missingness)
# 
# 1:2      coef exp(coef) se(coef)     z Pr(>|z|)
# AFP+ 0.3426    1.4087   0.3179 1.078    0.281
# 
# 1:2    exp(coef) exp(-coef) lower .95 upper .95
# AFP+     1.409     0.7099    0.7555     2.626
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.546  (se = 0.043 )
# Likelihood ratio test= 1.16  on 1 df,   p=0.3
# Wald test            = 1.16  on 1 df,   p=0.3
# Score (logrank) test = 1.17  on 1 df,   p=0.3


### Arg1 ###                                          
all_data$Arg1 <- factor(all_data$Arg1,levels = c("-","+"))
table(all_data$Arg1)
# -  + 
#   6 43

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, Arg1, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ Arg1, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ Arg1, data = data_factor, 
#         id = patient_id)
# 
# n= 45, number of events= 21 
# (209 observations deleted due to missingness)
# 
# 1:2       coef exp(coef) se(coef)     z Pr(>|z|)
# Arg1+ 0.7905    2.2045   1.0287 0.768    0.442
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# Arg1+     2.204     0.4536    0.2935     16.56
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.518  (se = 0.033 )
# Likelihood ratio test= 0.75  on 1 df,   p=0.4
# Wald test            = 0.59  on 1 df,   p=0.4
# Score (logrank) test = 0.62  on 1 df,   p=0.4


### CA199 ###                                          
all_data$CA199 <- factor(all_data$CA199,levels = c("-","+"))
table(all_data$CA199)
# - + 
#   5 0  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, CA199, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ CA199, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### CR ###                                          
all_data$CR <- factor(all_data$CR,levels = c("-","+"))
table(all_data$CR)
# - + 
#   1 0 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, CR, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ CR, data = data_factor, id = patient_id)
# Error in dimnames(events) <- list(count = tab1.levels, state = c(ystate,  : 
#                                                                    'dimnames' applied to non-array
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### CK ###                                          
all_data$CK <- factor(all_data$CK,levels = c("-","+"))
table(all_data$CK)
# - + 
#   0 5   

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, CK, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ CK, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### CK-pan ###    
colnames(all_data)[which(colnames(all_data)=="CK-pan")] <- "CK_pan"
all_data$CK_pan <- factor(all_data$CK_pan,levels = c("-","+"))
table(all_data$`CK_pan`)
# - + 
#   0 1  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, CK_pan, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ CK_pan, data = data_factor, id = patient_id)
# Error in dimnames(events) <- list(count = tab1.levels, state = c(ystate,  : 
#                                                                    'dimnames' applied to non-array
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### CK7 ###                                          
all_data$CK7 <- factor(all_data$CK7,levels = c("-","+"))
table(all_data$CK7)
# -  + 
#   73 35   

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, CK7, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ CK7, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ CK7, data = data_factor, 
#         id = patient_id)
# 
# n= 87, number of events= 41 
# (167 observations deleted due to missingness)
# 
# 1:2       coef exp(coef) se(coef)      z Pr(>|z|)
# CK7+ -0.1127    0.8934   0.3460 -0.326    0.745
# 
# 1:2    exp(coef) exp(-coef) lower .95 upper .95
# CK7+    0.8934      1.119    0.4534      1.76
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.516  (se = 0.04 )
# Likelihood ratio test= 0.11  on 1 df,   p=0.7
# Wald test            = 0.11  on 1 df,   p=0.7
# Score (logrank) test = 0.11  on 1 df,   p=0.7


### CK8 ###                                          
all_data$CK8 <- factor(all_data$CK8,levels = c("-","+"))
table(all_data$CK8)
# -  + 
#   2 19

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, CK8, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ CK8, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ CK8, data = data_factor, 
#         id = patient_id)
# 
# n= 19, number of events= 11 
# (235 observations deleted due to missingness)
# 
# 1:2      coef exp(coef) se(coef)     z Pr(>|z|)
# CK8+ 0.4821    1.6195   1.0634 0.453     0.65
# 
# 1:2    exp(coef) exp(-coef) lower .95 upper .95
# CK8+     1.619     0.6175    0.2015     13.02
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.512  (se = 0.052 )
# Likelihood ratio test= 0.23  on 1 df,   p=0.6
# Wald test            = 0.21  on 1 df,   p=0.7
# Score (logrank) test = 0.21  on 1 df,   p=0.6


### p_CEA ### 
colnames(all_data)[which(colnames(all_data)=="p-CEA")] <- "p_CEA"
all_data$p_CEA <- factor(all_data$p_CEA,levels = c("-","+"))
table(all_data$p_CEA)
# -  + 
#   5 37  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, p_CEA, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ p_CEA, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ p_CEA, data = data_factor, 
#         id = patient_id)
# 
# n= 28, number of events= 18 
# (226 observations deleted due to missingness)
# 
# 1:2        coef exp(coef) se(coef)    z Pr(>|z|)
# p_CEA+ 0.3638    1.4388   0.7578 0.48    0.631
# 
# 1:2      exp(coef) exp(-coef) lower .95 upper .95
# p_CEA+     1.439      0.695    0.3258     6.353
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.534  (se = 0.036 )
# Likelihood ratio test= 0.25  on 1 df,   p=0.6
# Wald test            = 0.23  on 1 df,   p=0.6
# Score (logrank) test = 0.23  on 1 df,   p=0.6


### CEA ###                                          
all_data$CEA <- factor(all_data$CEA,levels = c("-","+"))
table(all_data$CEA)
# - + 
#   2 0 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, CEA, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ CEA, data = data_factor, id = patient_id)
# Error in dimnames(events) <- list(count = tab1.levels, state = c(ystate,  : 
#                                                                    'dimnames' applied to non-array
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### CEA_poly ### 
colnames(all_data)[which(colnames(all_data)=="CEA poly")] <- "CEA_poly"
all_data$CEA_poly <- factor(all_data$CEA_poly,levels = c("-","+"))
table(all_data$CEA_poly)
# - + 
#   1 1 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, CEA_poly, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ CEA_poly, data = data_factor, id = patient_id)
# Warning message:
#   In fitter(X, Y, istrat, offset, init, control, weights = weights,  :
#               Ran out of iterations and did not converge
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ CEA_poly, data = data_factor, 
#         id = patient_id)
# 
# n= 2, number of events= 1 
# (252 observations deleted due to missingness)
# 
# 1:2               coef  exp(coef)   se(coef)      z Pr(>|z|)
# CEA_poly+ -2.120e+01  6.190e-10  4.019e+04 -0.001        1
# 
# 1:2         exp(coef) exp(-coef) lower .95 upper .95
# CEA_poly+  6.19e-10  1.615e+09         0       Inf
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 1  (se = 0 )
# Likelihood ratio test= 1.39  on 1 df,   p=0.2
# Wald test            = 0  on 1 df,   p=1
# Score (logrank) test = 1  on 1 df,   p=0.3


### Muc_1 ###   
colnames(all_data)[which(colnames(all_data)=="Muc-1")] <- "Muc_1"
all_data$Muc_1 <- factor(all_data$Muc_1,levels = c("-","+"))
table(all_data$Muc_1)
# -  + 
#   11  3 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, Muc_1, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ Muc_1, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ Muc_1, data = data_factor, 
#         id = patient_id)
# 
# n= 14, number of events= 6 
# (240 observations deleted due to missingness)
# 
# 1:2        coef exp(coef) se(coef)     z Pr(>|z|)
# Muc_1+ 1.0403    2.8299   0.9141 1.138    0.255
# 
# 1:2      exp(coef) exp(-coef) lower .95 upper .95
# Muc_1+      2.83     0.3534    0.4718     16.98
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.603  (se = 0.099 )
# Likelihood ratio test= 1.17  on 1 df,   p=0.3
# Wald test            = 1.3  on 1 df,   p=0.3
# Score (logrank) test = 1.42  on 1 df,   p=0.2


### CK19 ###                                          
all_data$CK19 <- factor(all_data$CK19,levels = c("-","+"))
table(all_data$CK19)
# -  + 
#   63 27 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, CK19, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ CK19, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ CK19, data = data_factor, 
#         id = patient_id)
# 
# n= 70, number of events= 36 
# (184 observations deleted due to missingness)
# 
# 1:2       coef exp(coef) se(coef)    z Pr(>|z|)
# CK19+ 0.6112    1.8426   0.3616 1.69    0.091
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# CK19+     1.843     0.5427     0.907     3.743
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.559  (se = 0.041 )
# Likelihood ratio test= 2.66  on 1 df,   p=0.1
# Wald test            = 2.86  on 1 df,   p=0.09
# Score (logrank) test = 2.94  on 1 df,   p=0.09


### CD31 ###                                          
all_data$CD31 <- factor(all_data$CD31,levels = c("-","+"))
table(all_data$CD31)
# - + 
#   0 2  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, CD31, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ CD31, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### CD10 ###                                          
all_data$CD10 <- factor(all_data$CD10,levels = c("-","+"))
table(all_data$CD10)
# -  + 
#   31 24

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, CD10, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ CD10, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ CD10, data = data_factor, 
#         id = patient_id)
# 
# n= 42, number of events= 24 
# (212 observations deleted due to missingness)
# 
# 1:2        coef exp(coef) se(coef)      z Pr(>|z|)
# CD10+ -0.6829    0.5052   0.4409 -1.549    0.121
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# CD10+    0.5052       1.98    0.2129     1.199
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.611  (se = 0.055 )
# Likelihood ratio test= 2.53  on 1 df,   p=0.1
# Wald test            = 2.4  on 1 df,   p=0.1
# Score (logrank) test = 2.49  on 1 df,   p=0.1


### CK8_CK18 ### 
colnames(all_data)[which(colnames(all_data)=="CK8/CK18")] <- "CK8_CK18"
all_data$CK8_CK18 <- factor(all_data$CK8_CK18,levels = c("-","+"))
table(all_data$CK8_CK18)
# - + 
#   1 1   

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, CK8_CK18, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ CK8_CK18, data = data_factor, id = patient_id)
# Warning message:
#   In fitter(X, Y, istrat, offset, init, control, weights = weights,  :
#               Ran out of iterations and did not converge
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ CK8_CK18, data = data_factor, 
#         id = patient_id)
# 
# n= 2, number of events= 1 
# (252 observations deleted due to missingness)
# 
# 1:2              coef exp(coef)  se(coef)     z Pr(>|z|)
# CK8_CK18+ 2.120e+01 1.615e+09 4.019e+04 0.001        1
# 
# 1:2         exp(coef) exp(-coef) lower .95 upper .95
# CK8_CK18+ 1.615e+09   6.19e-10         0       Inf
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 1  (se = 0 )
# Likelihood ratio test= 1.39  on 1 df,   p=0.2
# Wald test            = 0  on 1 df,   p=1
# Score (logrank) test = 1  on 1 df,   p=0.3


### CK20 ###                                          
all_data$CK20 <- factor(all_data$CK20,levels = c("-","+"))
table(all_data$CK20)
# - + 
#   4 0   

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, CK20, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ CK20, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### GPC3 ###                                          
all_data$GPC3 <- factor(all_data$GPC3,levels = c("-","+"))
table(all_data$GPC3)
# -  + 
#   10 86 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, GPC3, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ GPC3, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ GPC3, data = data_factor, 
#         id = patient_id)
# 
# n= 91, number of events= 45 
# (163 observations deleted due to missingness)
# 
# 1:2        coef exp(coef) se(coef)      z Pr(>|z|)
# GPC3+ -0.2473    0.7809   0.4396 -0.562    0.574
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# GPC3+    0.7809      1.281    0.3299     1.849
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.5  (se = 0.021 )
# Likelihood ratio test= 0.3  on 1 df,   p=0.6
# Wald test            = 0.32  on 1 df,   p=0.6
# Score (logrank) test = 0.32  on 1 df,   p=0.6


### D2_40 ###     
colnames(all_data)[which(colnames(all_data)=="D2-40")] <- "D2_40"
all_data$D2_40 <- factor(all_data$D2_40,levels = c("-","+"))
table(all_data$D2_40)
# - + 
#   3 3  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, D2_40, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ D2_40, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ D2_40, data = data_factor, 
#         id = patient_id)
# 
# n= 6, number of events= 2 
# (248 observations deleted due to missingness)
# 
# 1:2        coef exp(coef) se(coef)     z Pr(>|z|)
# D2_40+ 0.3466    1.4142   1.4355 0.241    0.809
# 
# 1:2      exp(coef) exp(-coef) lower .95 upper .95
# D2_40+     1.414     0.7071   0.08484     23.57
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.643  (se = 0.181 )
# Likelihood ratio test= 0.06  on 1 df,   p=0.8
# Wald test            = 0.06  on 1 df,   p=0.8
# Score (logrank) test = 0.06  on 1 df,   p=0.8


### GPC3 ###                                          
all_data$GPC3 <- factor(all_data$GPC3,levels = c("-","+"))
table(all_data$GPC3)
# -  + 
#   10 86  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, GPC3, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ GPC3, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ GPC3, data = data_factor, 
#         id = patient_id)
# 
# n= 91, number of events= 45 
# (163 observations deleted due to missingness)
# 
# 1:2        coef exp(coef) se(coef)      z Pr(>|z|)
# GPC3+ -0.2473    0.7809   0.4396 -0.562    0.574
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# GPC3+    0.7809      1.281    0.3299     1.849
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.5  (se = 0.021 )
# Likelihood ratio test= 0.3  on 1 df,   p=0.6
# Wald test            = 0.32  on 1 df,   p=0.6
# Score (logrank) test = 0.32  on 1 df,   p=0.6


### Hep1 ###                                          
all_data$Hep1 <- factor(all_data$Hep1,levels = c("-","+"))
table(all_data$Hep1)
# - + 
#   0 1   

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, Hep1, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ Hep1, data = data_factor, id = patient_id)
# Error in dimnames(events) <- list(count = tab1.levels, state = c(ystate,  : 
#                                                                    'dimnames' applied to non-array
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### Hept1 ###                                          
all_data$Hept1 <- factor(all_data$Hept1,levels = c("-","+"))
table(all_data$Hept1)
# - + 
#   0 1  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, Hept1, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ Hept1, data = data_factor, id = patient_id)
# Error in dimnames(events) <- list(count = tab1.levels, state = c(ystate,  : 
#                                                                    'dimnames' applied to non-array
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### HCGβ ###                                          
all_data$HCGβ <- factor(all_data$HCGβ,levels = c("-","+"))
table(all_data$HCGβ)
# - + 
#   1 0  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, HCGβ, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ HCGβ, data = data_factor, id = patient_id)
# Error in dimnames(events) <- list(count = tab1.levels, state = c(ystate,  : 
#                                                                    'dimnames' applied to non-array
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### HSP70 ###                                          
all_data$HSP70 <- factor(all_data$HSP70,levels = c("-","+"))
table(all_data$HSP70)
# -  + 
#   0 9  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, HSP70, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ HSP70, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### HMB45 ###                                          
all_data$HMB45 <- factor(all_data$HMB45,levels = c("-","+"))
table(all_data$HMB45)
# - + 
#   5 0  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, HMB45, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ HMB45, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### Inhibin_a ###    
colnames(all_data)[which(colnames(all_data)=="Inhibin a")] <- "Inhibin_a"
all_data$Inhibin_a <- factor(all_data$Inhibin_a,levels = c("-","+"))
table(all_data$Inhibin_a)
# - + 
#   2 0 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, Inhibin_a, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ Inhibin_a, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### SMA ###                                          
all_data$SMA <- factor(all_data$SMA,levels = c("-","+"))
table(all_data$SMA)
# - + 
#   2 0  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, SMA, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ SMA, data = data_factor, id = patient_id)
# Error in dimnames(events) <- list(count = tab1.levels, state = c(ystate,  : 
#                                                                    'dimnames' applied to non-array
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### CgA ###                                          
all_data$CgA <- factor(all_data$CgA,levels = c("-","+"))
table(all_data$CgA)
# - + 
#   5 1  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, CgA, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ CgA, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ CgA, data = data_factor, 
#         id = patient_id)
# 
# n= 4, number of events= 3 
# (250 observations deleted due to missingness)
# 
# 1:2       coef exp(coef) se(coef)      z Pr(>|z|)
# CgA+ -0.1285    0.8794   1.2535 -0.103    0.918
# 
# 1:2    exp(coef) exp(-coef) lower .95 upper .95
# CgA+    0.8794      1.137   0.07537     10.26
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.583  (se = 0.144 )
# Likelihood ratio test= 0.01  on 1 df,   p=0.9
# Wald test            = 0.01  on 1 df,   p=0.9
# Score (logrank) test = 0.01  on 1 df,   p=0.9


### Syn ###                                          
all_data$Syn <- factor(all_data$Syn,levels = c("-","+"))
table(all_data$Syn)
# - + 
#   8 2  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, Syn, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ Syn, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ Syn, data = data_factor, 
#         id = patient_id)
# 
# n= 8, number of events= 5 
# (246 observations deleted due to missingness)
# 
# 1:2       coef exp(coef) se(coef)      z Pr(>|z|)
# Syn+ -0.8084    0.4456   1.1325 -0.714    0.475
# 
# 1:2    exp(coef) exp(-coef) lower .95 upper .95
# Syn+    0.4456      2.244   0.04841     4.101
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.625  (se = 0.089 )
# Likelihood ratio test= 0.59  on 1 df,   p=0.4
# Wald test            = 0.51  on 1 df,   p=0.5
# Score (logrank) test = 0.54  on 1 df,   p=0.5


### SF_1 ###  
colnames(all_data)[which(colnames(all_data)=="SF-1")] <- "SF_1"
all_data$SF_1 <- factor(all_data$SF_1,levels = c("-","+"))
table(all_data$SF_1)
# - + 
#   0 1

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, SF_1, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ SF_1, data = data_factor, id = patient_id)
# Error in dimnames(events) <- list(count = tab1.levels, state = c(ystate,  : 
#                                                                    'dimnames' applied to non-array
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### PAX_8 ###   
colnames(all_data)[which(colnames(all_data)=="PAX-8")] <- "PAX_8"
all_data$PAX_8 <- factor(all_data$PAX_8,levels = c("-","+"))
table(all_data$PAX_8)
# - + 
#   1 0  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, PAX_8, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ PAX_8, data = data_factor, id = patient_id)
# Error in dimnames(events) <- list(count = tab1.levels, state = c(ystate,  : 
#                                                                    'dimnames' applied to non-array
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### P40 ###                                          
all_data$P40 <- factor(all_data$P40,levels = c("-","+"))
table(all_data$P40)
# - + 
#   1 0 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, P40, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ P40, data = data_factor, id = patient_id)
# Error in dimnames(events) <- list(count = tab1.levels, state = c(ystate,  : 
#                                                                    'dimnames' applied to non-array
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### P53 ###                                          
all_data$P53 <- factor(all_data$P53,levels = c("-","+"))
table(all_data$P53)
# - + 
#   1 1  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, P53, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ P53, data = data_factor, id = patient_id)
# Warning message:
#   In fitter(X, Y, istrat, offset, init, control, weights = weights,  :
#               Ran out of iterations and did not converge
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ P53, data = data_factor, 
#         id = patient_id)
# 
# n= 2, number of events= 1 
# (252 observations deleted due to missingness)
# 
# 1:2          coef  exp(coef)   se(coef)      z Pr(>|z|)
# P53+ -2.120e+01  6.190e-10  4.019e+04 -0.001        1
# 
# 1:2    exp(coef) exp(-coef) lower .95 upper .95
# P53+  6.19e-10  1.615e+09         0       Inf
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 1  (se = 0 )
# Likelihood ratio test= 1.39  on 1 df,   p=0.2
# Wald test            = 0  on 1 df,   p=1
# Score (logrank) test = 1  on 1 df,   p=0.3


### Ki_67 ###   
colnames(all_data)[which(colnames(all_data)=="Ki-67")] <- "Ki_67"
summary(all_data$Ki_67)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0500  0.2000  0.3000  0.3482  0.5000  0.7000     134  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, Ki_67, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ Ki_67, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ Ki_67, data = data_factor, 
#         id = patient_id)
# 
# n= 94, number of events= 47 
# (160 observations deleted due to missingness)
# 
# 1:2       coef exp(coef) se(coef)    z Pr(>|z|)
# Ki_67 1.7200    5.5848   0.7928 2.17     0.03
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# Ki_67     5.585     0.1791     1.181     26.41
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.575  (se = 0.047 )
# Likelihood ratio test= 4.59  on 1 df,   p=0.03
# Wald test            = 4.71  on 1 df,   p=0.03
# Score (logrank) test = 4.8  on 1 df,   p=0.03

mylogit <- glm(Recurrence ~ Ki_67, data = data_factor,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ Ki_67, family = binomial, data = data_factor)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.3861  -1.1178  -0.9864   1.2382   1.4090  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -0.6071     0.3979  -1.526    0.127
# Ki_67         1.5507     1.0119   1.533    0.125
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 166.22  on 119  degrees of freedom
# Residual deviance: 163.82  on 118  degrees of freedom
# (134 observations deleted due to missingness)
# AIC: 167.82
# 
# Number of Fisher Scoring iterations: 4


### CD34 ###                                          
all_data$CD34 <- factor(all_data$CD34,levels = c("-","+"))
table(all_data$CD34)
# -   + 
#   1 113 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, CD34, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ CD34, data = data_factor, id = patient_id)
# Warning message:
#   In fitter(X, Y, istrat, offset, init, control, weights = weights,  :
#               Loglik converged before variable  1 ; coefficient may be infinite. 
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ CD34, data = data_factor, 
#         id = patient_id)
# 
# n= 89, number of events= 45 
# (165 observations deleted due to missingness)
# 
# 1:2          coef exp(coef)  se(coef)     z Pr(>|z|)
# CD34+ 1.605e+01 9.336e+06 3.173e+03 0.005    0.996
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# CD34+   9336048  1.071e-07         0       Inf
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.508  (se = 0.008 )
# Likelihood ratio test= 1.83  on 1 df,   p=0.2
# Wald test            = 0  on 1 df,   p=1
# Score (logrank) test = 0.93  on 1 df,   p=0.3


#### group 3: biochemistry test ####

### Alpha_fetoprotein ###     
colnames(all_data)[which(colnames(all_data)=="Alpha-fetoprotein (ng/ml)")] <- "Alpha_fetoprotein"
summary(all_data$Alpha_fetoprotein)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
#       1.0       6.1      62.7   51958.2     790.7 1936000.0       135 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, Alpha_fetoprotein, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ Alpha_fetoprotein, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ Alpha_fetoprotein, data = data_factor, 
#         id = patient_id)
# 
# n= 92, number of events= 46 
# (162 observations deleted due to missingness)
# 
# 1:2                      coef exp(coef)  se(coef)    z Pr(>|z|)
# Alpha_fetoprotein 1.057e-06 1.000e+00 3.633e-07 2.91  0.00361
# 
# 1:2                 exp(coef) exp(-coef) lower .95 upper .95
# Alpha_fetoprotein         1          1         1         1
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.64  (se = 0.04 )
# Likelihood ratio test= 5.72  on 1 df,   p=0.02
# Wald test            = 8.47  on 1 df,   p=0.004
# Score (logrank) test = 10.35  on 1 df,   p=0.001

mylogit <- glm(Recurrence ~ Alpha_fetoprotein, data = data_factor,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ Alpha_fetoprotein, family = binomial, 
#       data = data_factor)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -1.471  -1.117  -1.117   1.239   1.239  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)       -1.447e-01  1.900e-01  -0.761    0.446
# Alpha_fetoprotein  6.805e-06  9.395e-06   0.724    0.469
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 164.89  on 118  degrees of freedom
# Residual deviance: 158.80  on 117  degrees of freedom
# (135 observations deleted due to missingness)
# AIC: 162.8
# 
# Number of Fisher Scoring iterations: 8


### abnormality_of_alpha_fetoprotein ###    
colnames(all_data)[which(colnames(all_data)=="Abnormality of alpha-fetoprotein")] <- "abnormality_of_alpha_fetoprotein"
all_data$abnormality_of_alpha_fetoprotein <- factor(all_data$abnormality_of_alpha_fetoprotein,levels = c("normal","high"))
table(all_data$abnormality_of_alpha_fetoprotein)
# normal   high 
# 35     84  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, abnormality_of_alpha_fetoprotein, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ abnormality_of_alpha_fetoprotein, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ abnormality_of_alpha_fetoprotein, 
#         data = data_factor, id = patient_id)
# 
# n= 92, number of events= 46 
# (162 observations deleted due to missingness)
# 
# 1:2                                      coef exp(coef) se(coef)     z Pr(>|z|)
# abnormality_of_alpha_fetoproteinhigh 0.3941    1.4831   0.3466 1.137    0.255
# 
# 1:2                                    exp(coef) exp(-coef) lower .95 upper .95
# abnormality_of_alpha_fetoproteinhigh     1.483     0.6743    0.7519     2.925
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.573  (se = 0.03 )
# Likelihood ratio test= 1.38  on 1 df,   p=0.2
# Wald test            = 1.29  on 1 df,   p=0.3
# Score (logrank) test = 1.31  on 1 df,   p=0.3


### carcinoembryonic_antigen ###  
colnames(all_data)[which(colnames(all_data)=="Carcinoembryonic antigen")] <- "carcinoembryonic_antigen"
summary(all_data$carcinoembryonic_antigen)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.620   1.670   2.250   2.956   3.560  12.000     153  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, carcinoembryonic_antigen, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ carcinoembryonic_antigen, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ carcinoembryonic_antigen, 
#         data = data_factor, id = patient_id)
# 
# n= 82, number of events= 44 
# (172 observations deleted due to missingness)
# 
# 1:2                            coef exp(coef) se(coef)      z Pr(>|z|)
# carcinoembryonic_antigen -0.13317   0.87531  0.08574 -1.553     0.12
# 
# 1:2                        exp(coef) exp(-coef) lower .95 upper .95
# carcinoembryonic_antigen    0.8753      1.142    0.7399     1.035
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.562  (se = 0.046 )
# Likelihood ratio test= 2.74  on 1 df,   p=0.1
# Wald test            = 2.41  on 1 df,   p=0.1
# Score (logrank) test = 2.43  on 1 df,   p=0.1


### abnormality_of_carcinoembryonic_antigen ### 
colnames(all_data)[which(colnames(all_data)=="Abnormality of carcinoembryonic antigen")] <- "abnormality_of_carcinoembryonic_antigen"
all_data$abnormality_of_carcinoembryonic_antigen <- factor(all_data$abnormality_of_carcinoembryonic_antigen,levels = c("normal","high"))
table(all_data$abnormality_of_carcinoembryonic_antigen)
# normal   high 
# 98      6  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, abnormality_of_carcinoembryonic_antigen, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ abnormality_of_carcinoembryonic_antigen, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ abnormality_of_carcinoembryonic_antigen, 
#         data = data_factor, id = patient_id)
# 
# n= 85, number of events= 45 
# (169 observations deleted due to missingness)
# 
# 1:2                                              coef exp(coef) se(coef)      z Pr(>|z|)
# abnormality_of_carcinoembryonic_antigenhigh -0.9161    0.4001   1.0129 -0.904    0.366
# 
# 1:2                                           exp(coef) exp(-coef) lower .95 upper .95
# abnormality_of_carcinoembryonic_antigenhigh    0.4001        2.5   0.05495     2.913
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.508  (se = 0.022 )
# Likelihood ratio test= 1.11  on 1 df,   p=0.3
# Wald test            = 0.82  on 1 df,   p=0.4
# Score (logrank) test = 0.88  on 1 df,   p=0.3


### carbohydrate_antigen_CA125 ### 
colnames(all_data)[which(colnames(all_data)=="Carbohydrate antigen CA125")] <- "carbohydrate_antigen_CA125"
summary(all_data$carbohydrate_antigen_CA125)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    1.52    6.60    8.80   30.70   22.30  423.33     223   

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, carbohydrate_antigen_CA125, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ carbohydrate_antigen_CA125, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ carbohydrate_antigen_CA125, 
#         data = data_factor, id = patient_id)
# 
# n= 23, number of events= 11 
# (231 observations deleted due to missingness)
# 
# 1:2                              coef exp(coef) se(coef)     z Pr(>|z|)
# carbohydrate_antigen_CA125 0.001915  1.001917 0.002349 0.815    0.415
# 
# 1:2                          exp(coef) exp(-coef) lower .95 upper .95
# carbohydrate_antigen_CA125     1.002     0.9981    0.9973     1.007
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.589  (se = 0.098 )
# Likelihood ratio test= 0.55  on 1 df,   p=0.5
# Wald test            = 0.67  on 1 df,   p=0.4
# Score (logrank) test = 0.7  on 1 df,   p=0.4


### abnormality_of_carbohydrate_antigen_CA125 ###    
colnames(all_data)[which(colnames(all_data)=="Abnormality of carbohydrate antigen CA125")] <- "abnormality_of_carbohydrate_antigen_CA125"
all_data$abnormality_of_carbohydrate_antigen_CA125 <- factor(all_data$abnormality_of_carbohydrate_antigen_CA125,levels = c("normal","high"))
table(all_data$abnormality_of_carbohydrate_antigen_CA125)
# normal   high 
# 26      5  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, abnormality_of_carbohydrate_antigen_CA125, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ abnormality_of_carbohydrate_antigen_CA125, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ abnormality_of_carbohydrate_antigen_CA125, 
#         data = data_factor, id = patient_id)
# 
# n= 23, number of events= 11 
# (231 observations deleted due to missingness)
# 
# 1:2                                               coef exp(coef) se(coef)     z Pr(>|z|)
# abnormality_of_carbohydrate_antigen_CA125high 1.2199    3.3867   0.6485 1.881     0.06
# 
# 1:2                                             exp(coef) exp(-coef) lower .95 upper .95
# abnormality_of_carbohydrate_antigen_CA125high     3.387     0.2953    0.9501     12.07
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.583  (se = 0.059 )
# Likelihood ratio test= 3.12  on 1 df,   p=0.08
# Wald test            = 3.54  on 1 df,   p=0.06
# Score (logrank) test = 3.99  on 1 df,   p=0.05


### carbohydrate_antigen_CA153 ###   
colnames(all_data)[which(colnames(all_data)=="Carbohydrate antigen CA153")] <- "carbohydrate_antigen_CA153"
summary(all_data$carbohydrate_antigen_CA153)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.160   6.135   8.590  10.560  13.120  29.500     231  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, carbohydrate_antigen_CA153, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ carbohydrate_antigen_CA153, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ carbohydrate_antigen_CA153, 
#         data = data_factor, id = patient_id)
# 
# n= 17, number of events= 9 
# (237 observations deleted due to missingness)
# 
# 1:2                             coef exp(coef) se(coef)     z Pr(>|z|)
# carbohydrate_antigen_CA153 0.01009   1.01014  0.05101 0.198    0.843
# 
# 1:2                          exp(coef) exp(-coef) lower .95 upper .95
# carbohydrate_antigen_CA153      1.01       0.99     0.914     1.116
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.446  (se = 0.111 )
# Likelihood ratio test= 0.04  on 1 df,   p=0.8
# Wald test            = 0.04  on 1 df,   p=0.8
# Score (logrank) test = 0.04  on 1 df,   p=0.8


### abnormality_of_carbohydrate_antigen_CA153 ###   
colnames(all_data)[which(colnames(all_data)=="Abnormality of  carbohydrate antigen CA153")] <- "abnormality_of_carbohydrate_antigen_CA153"
all_data$abnormality_of_carbohydrate_antigen_CA153 <- factor(all_data$abnormality_of_carbohydrate_antigen_CA153,levels = c("normal","high"))
table(all_data$abnormality_of_carbohydrate_antigen_CA153)
# normal   high 
# 23      0 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, abnormality_of_carbohydrate_antigen_CA153, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ abnormality_of_carbohydrate_antigen_CA153, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### carbohydrate_antigen_CA199 ###  
colnames(all_data)[which(colnames(all_data)=="Carbohydrate antigen CA199")] <- "carbohydrate_antigen_CA199"
summary(all_data$carbohydrate_antigen_CA199)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.60    6.63   12.50   20.60   25.49  224.80     149 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, carbohydrate_antigen_CA199, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ carbohydrate_antigen_CA199, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ carbohydrate_antigen_CA199, 
#         data = data_factor, id = patient_id)
# 
# n= 87, number of events= 46 
# (167 observations deleted due to missingness)
# 
# 1:2                              coef exp(coef) se(coef)     z Pr(>|z|)
# carbohydrate_antigen_CA199 0.005630  1.005646 0.004128 1.364    0.173
# 
# 1:2                          exp(coef) exp(-coef) lower .95 upper .95
# carbohydrate_antigen_CA199     1.006     0.9944    0.9975     1.014
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.571  (se = 0.046 )
# Likelihood ratio test= 1.52  on 1 df,   p=0.2
# Wald test            = 1.86  on 1 df,   p=0.2
# Score (logrank) test = 1.89  on 1 df,   p=0.2


### abnormality_of_carbohydrate_antigen_CA199 ### 
colnames(all_data)[which(colnames(all_data)=="Abnormality of  carbohydrate antigen CA199")] <- "abnormality_of_carbohydrate_antigen_CA199"
all_data$abnormality_of_carbohydrate_antigen_CA199 <- factor(all_data$abnormality_of_carbohydrate_antigen_CA199,levels = c("normal","high"))
table(all_data$abnormality_of_carbohydrate_antigen_CA199)
# normal   high 
# 88     18 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, abnormality_of_carbohydrate_antigen_CA199, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ abnormality_of_carbohydrate_antigen_CA199, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ abnormality_of_carbohydrate_antigen_CA199, 
#         data = data_factor, id = patient_id)
# 
# n= 88, number of events= 47 
# (166 observations deleted due to missingness)
# 
# 1:2                                                 coef exp(coef) se(coef)      z Pr(>|z|)
# abnormality_of_carbohydrate_antigen_CA199high -0.05294   0.94843  0.37228 -0.142    0.887
# 
# 1:2                                             exp(coef) exp(-coef) lower .95 upper .95
# abnormality_of_carbohydrate_antigen_CA199high    0.9484      1.054    0.4572     1.967
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.499  (se = 0.032 )
# Likelihood ratio test= 0.02  on 1 df,   p=0.9
# Wald test            = 0.02  on 1 df,   p=0.9
# Score (logrank) test = 0.02  on 1 df,   p=0.9


### tumor_associated_glycoprotein_CA72_4 ###  
colnames(all_data)[which(colnames(all_data)=="Tumor-associated glycoprotein CA72-4")] <- "tumor_associated_glycoprotein_CA72_4"
summary(all_data$tumor_associated_glycoprotein_CA72_4)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.24    0.55    1.00    2.70    3.05   11.22     239 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, tumor_associated_glycoprotein_CA72_4, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ tumor_associated_glycoprotein_CA72_4, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ tumor_associated_glycoprotein_CA72_4, 
#         data = data_factor, id = patient_id)
# 
# n= 10, number of events= 6 
# (244 observations deleted due to missingness)
# 
# 1:2                                      coef exp(coef) se(coef)     z Pr(>|z|)
# tumor_associated_glycoprotein_CA72_4 0.1033    1.1088   0.1858 0.556    0.578
# 
# 1:2                                    exp(coef) exp(-coef) lower .95 upper .95
# tumor_associated_glycoprotein_CA72_4     1.109     0.9019    0.7703     1.596
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.514  (se = 0.105 )
# Likelihood ratio test= 0.28  on 1 df,   p=0.6
# Wald test            = 0.31  on 1 df,   p=0.6
# Score (logrank) test = 0.32  on 1 df,   p=0.6


### abnormality_of_tumor_associated_glycoprotein_CA72_4 ### 
colnames(all_data)[which(colnames(all_data)=="Abnormality of tumor-associated glycoprotein CA72-4")] <- "abnormality_of_tumor_associated_glycoprotein_CA72_4"
all_data$abnormality_of_tumor_associated_glycoprotein_CA72_4 <- factor(all_data$abnormality_of_tumor_associated_glycoprotein_CA72_4,levels = c("normal","high"))
table(all_data$abnormality_of_tumor_associated_glycoprotein_CA72_4)
# normal   high 
# 15      0  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, abnormality_of_tumor_associated_glycoprotein_CA72_4, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ abnormality_of_tumor_associated_glycoprotein_CA72_4, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### carbohydrate_antigen_CA50 ###                                          
colnames(all_data)[which(colnames(all_data)=="Carbohydrate antigen CA50")] <- "carbohydrate_antigen_CA50"
summary(all_data$carbohydrate_antigen_CA50)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.50    3.60    6.61   16.36   16.20  163.72     201 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, carbohydrate_antigen_CA50, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ carbohydrate_antigen_CA50, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ carbohydrate_antigen_CA50, 
#         data = data_factor, id = patient_id)
# 
# n= 40, number of events= 18 
# (214 observations deleted due to missingness)
# 
# 1:2                              coef exp(coef)  se(coef)     z Pr(>|z|)
# carbohydrate_antigen_CA50 0.0006995 1.0006997 0.0083583 0.084    0.933
# 
# 1:2                         exp(coef) exp(-coef) lower .95 upper .95
# carbohydrate_antigen_CA50     1.001     0.9993    0.9844     1.017
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.465  (se = 0.077 )
# Likelihood ratio test= 0.01  on 1 df,   p=0.9
# Wald test            = 0.01  on 1 df,   p=0.9
# Score (logrank) test = 0.01  on 1 df,   p=0.9


### abnormality_of_carbohydrate_antigen_CA50 ###   
colnames(all_data)[which(colnames(all_data)=="Abnormality of carbohydrate antigen CA50")] <- "abnormality_of_carbohydrate_antigen_CA50"
all_data$abnormality_of_carbohydrate_antigen_CA50 <- factor(all_data$abnormality_of_carbohydrate_antigen_CA50,levels = c("normal","high"))
table(all_data$abnormality_of_carbohydrate_antigen_CA50)
# normal   high 
# 43     10 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, abnormality_of_carbohydrate_antigen_CA50, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ abnormality_of_carbohydrate_antigen_CA50, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ abnormality_of_carbohydrate_antigen_CA50, 
#         data = data_factor, id = patient_id)
# 
# n= 40, number of events= 18 
# (214 observations deleted due to missingness)
# 
# 1:2                                               coef exp(coef) se(coef)      z Pr(>|z|)
# abnormality_of_carbohydrate_antigen_CA50high -0.3631    0.6955   0.6389 -0.568     0.57
# 
# 1:2                                            exp(coef) exp(-coef) lower .95 upper .95
# abnormality_of_carbohydrate_antigen_CA50high    0.6955      1.438    0.1988     2.433
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.544  (se = 0.048 )
# Likelihood ratio test= 0.35  on 1 df,   p=0.6
# Wald test            = 0.32  on 1 df,   p=0.6
# Score (logrank) test = 0.33  on 1 df,   p=0.6


### cytokeratin_fragment ###   
colnames(all_data)[which(colnames(all_data)=="Cytokeratin fragment")] <- "cytokeratin_fragment"
summary(all_data$cytokeratin_fragment)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.420   2.090   2.330   2.469   2.750   3.690     243 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, cytokeratin_fragment, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ cytokeratin_fragment, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ cytokeratin_fragment, 
#         data = data_factor, id = patient_id)
# 
# n= 7, number of events= 4 
# (247 observations deleted due to missingness)
# 
# 1:2                     coef exp(coef) se(coef)     z Pr(>|z|)
# cytokeratin_fragment 1.536     4.644    1.104 1.391    0.164
# 
# 1:2                    exp(coef) exp(-coef) lower .95 upper .95
# cytokeratin_fragment     4.644     0.2153    0.5337     40.41
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.667  (se = 0.203 )
# Likelihood ratio test= 2.23  on 1 df,   p=0.1
# Wald test            = 1.94  on 1 df,   p=0.2
# Score (logrank) test = 2.25  on 1 df,   p=0.1


### abnormality_of_cytokeratin_fragment ###   
colnames(all_data)[which(colnames(all_data)=="Abnormality of  cytokeratin fragment")] <- "abnormality_of_cytokeratin_fragment"
all_data$abnormality_of_cytokeratin_fragment <- factor(all_data$abnormality_of_cytokeratin_fragment,levels = c("normal","high"))
table(all_data$abnormality_of_cytokeratin_fragment)
# normal   high 
# 11      1 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, abnormality_of_cytokeratin_fragment, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ abnormality_of_cytokeratin_fragment, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### cytokeratin_19_fragment ###  
colnames(all_data)[which(colnames(all_data)=="Cytokeratin 19 fragment")] <- "cytokeratin_19_fragment"
summary(all_data$cytokeratin_19_fragment)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.240   1.438   1.925   2.312   2.560   5.300     246 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, cytokeratin_19_fragment, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ cytokeratin_19_fragment, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ cytokeratin_19_fragment, 
#         data = data_factor, id = patient_id)
# 
# n= 7, number of events= 3 
# (247 observations deleted due to missingness)
# 
# 1:2                          coef exp(coef) se(coef)      z Pr(>|z|)
# cytokeratin_19_fragment -0.3069    0.7357   0.5399 -0.568     0.57
# 
# 1:2                       exp(coef) exp(-coef) lower .95 upper .95
# cytokeratin_19_fragment    0.7357      1.359    0.2554      2.12
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.571  (se = 0.184 )
# Likelihood ratio test= 0.38  on 1 df,   p=0.5
# Wald test            = 0.32  on 1 df,   p=0.6
# Score (logrank) test = 0.34  on 1 df,   p=0.6


### abnormality_of_cytokeratin_19_fragment ###   
colnames(all_data)[which(colnames(all_data)=="Abnormality of cytokeratin 19 fragment")] <- "abnormality_of_cytokeratin_19_fragment"
all_data$abnormality_of_cytokeratin_19_fragment <- factor(all_data$abnormality_of_cytokeratin_19_fragment,levels = c("normal","high"))
table(all_data$abnormality_of_cytokeratin_19_fragment)
# normal   high 
# 8      0  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, abnormality_of_cytokeratin_19_fragment, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ abnormality_of_cytokeratin_19_fragment, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
#cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out


### hepatitis_B_virus_surface_antigen_quantification ###                                          
colnames(all_data)[which(colnames(all_data)=="Hepatitis B virus surface antigen quantification")] <- "hepatitis_B_virus_surface_antigen_quantification"
summary(all_data$hepatitis_B_virus_surface_antigen_quantification)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.00   11.65   19.30  356.14  116.26 7156.00     118 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, hepatitis_B_virus_surface_antigen_quantification, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ hepatitis_B_virus_surface_antigen_quantification, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ hepatitis_B_virus_surface_antigen_quantification, 
#         data = data_factor, id = patient_id)
# 
# n= 105, number of events= 52 
# (149 observations deleted due to missingness)
# 
# 1:2                                                     coef exp(coef)  se(coef)     z Pr(>|z|)
# hepatitis_B_virus_surface_antigen_quantification 0.0001509 1.0001509 0.0001079 1.398    0.162
# 
# 1:2                                                exp(coef) exp(-coef) lower .95 upper .95
# hepatitis_B_virus_surface_antigen_quantification         1     0.9998    0.9999         1
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.485  (se = 0.044 )
# Likelihood ratio test= 1.53  on 1 df,   p=0.2
# Wald test            = 1.96  on 1 df,   p=0.2
# Score (logrank) test = 2.04  on 1 df,   p=0.2


### hepatitis_B_virus_e_antigen_quantification ###                                          
colnames(all_data)[which(colnames(all_data)=="Hepatitis B virus e antigen quantification")] <- "hepatitis_B_virus_e_antigen_quantification"
summary(all_data$hepatitis_B_virus_e_antigen_quantification)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
#    0.0700    0.3000    0.7150   16.1530    0.9125 1429.8000       118 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, hepatitis_B_virus_e_antigen_quantification, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ hepatitis_B_virus_e_antigen_quantification, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ hepatitis_B_virus_e_antigen_quantification, 
#         data = data_factor, id = patient_id)
# 
# n= 105, number of events= 52 
# (149 observations deleted due to missingness)
# 
# 1:2                                               coef exp(coef)  se(coef)     z Pr(>|z|)
# hepatitis_B_virus_e_antigen_quantification 0.0006993 1.0006996 0.0006822 1.025    0.305
# 
# 1:2                                          exp(coef) exp(-coef) lower .95 upper .95
# hepatitis_B_virus_e_antigen_quantification     1.001     0.9993    0.9994     1.002
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.538  (se = 0.044 )
# Likelihood ratio test= 0.79  on 1 df,   p=0.4
# Wald test            = 1.05  on 1 df,   p=0.3
# Score (logrank) test = 1.14  on 1 df,   p=0.3


### hepatitis_B_virus_surface_antibody_quantification ###                                          
colnames(all_data)[which(colnames(all_data)=="Hepatitis B virus surface antibody quantification")] <- "hepatitis_B_virus_surface_antibody_quantification"
summary(all_data$hepatitis_B_virus_surface_antibody_quantification)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#    0.000    1.915    2.000   20.879    2.000 1000.000      118 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, hepatitis_B_virus_surface_antibody_quantification, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ hepatitis_B_virus_surface_antibody_quantification, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ hepatitis_B_virus_surface_antibody_quantification, 
#         data = data_factor, id = patient_id)
# 
# n= 105, number of events= 52 
# (149 observations deleted due to missingness)
# 
# 1:2                                                       coef  exp(coef)   se(coef)     z Pr(>|z|)
# hepatitis_B_virus_surface_antibody_quantification -0.0007639  0.9992364  0.0028317 -0.27    0.787
# 
# 1:2                                                 exp(coef) exp(-coef) lower .95 upper .95
# hepatitis_B_virus_surface_antibody_quantification    0.9992      1.001    0.9937     1.005
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.484  (se = 0.039 )
# Likelihood ratio test= 0.08  on 1 df,   p=0.8
# Wald test            = 0.07  on 1 df,   p=0.8
# Score (logrank) test = 0.07  on 1 df,   p=0.8


### hepatitis_B_virus_e_antibody_quantification ###                                          
colnames(all_data)[which(colnames(all_data)=="Hepatitis B virus e antibody quantification")] <- "hepatitis_B_virus_e_antibody_quantification"
summary(all_data$hepatitis_B_virus_e_antibody_quantification)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#   0.0000   0.1575   1.0450  12.3624   1.8100 291.5000      144  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, hepatitis_B_virus_e_antibody_quantification, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ hepatitis_B_virus_e_antibody_quantification, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ hepatitis_B_virus_e_antibody_quantification, 
#         data = data_factor, id = patient_id)
# 
# n= 82, number of events= 40 
# (172 observations deleted due to missingness)
# 
# 1:2                                               coef exp(coef) se(coef)     z Pr(>|z|)
# hepatitis_B_virus_e_antibody_quantification 0.002145  1.002147 0.002718 0.789     0.43
# 
# 1:2                                           exp(coef) exp(-coef) lower .95 upper .95
# hepatitis_B_virus_e_antibody_quantification     1.002     0.9979    0.9968     1.008
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.512  (se = 0.051 )
# Likelihood ratio test= 0.53  on 1 df,   p=0.5
# Wald test            = 0.62  on 1 df,   p=0.4
# Score (logrank) test = 0.64  on 1 df,   p=0.4


### hepatitis_B_virus_core_antibody_quantification ###                                          
colnames(all_data)[which(colnames(all_data)=="Hepatitis B virus core antibody quantification")] <- "hepatitis_B_virus_core_antibody_quantification"
all_data$hepatitis_B_virus_core_antibody_quantification <- as.numeric(all_data$hepatitis_B_virus_core_antibody_quantification)
summary(all_data$hepatitis_B_virus_core_antibody_quantification)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.01    0.10    0.10   54.01    5.00  644.90     119 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, hepatitis_B_virus_core_antibody_quantification, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ hepatitis_B_virus_core_antibody_quantification, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ hepatitis_B_virus_core_antibody_quantification, 
#         data = data_factor, id = patient_id)
# 
# n= 104, number of events= 52 
# (150 observations deleted due to missingness)
# 
# 1:2                                                   coef exp(coef)  se(coef)     z Pr(>|z|)
# hepatitis_B_virus_core_antibody_quantification 0.0003199 1.0003200 0.0008369 0.382    0.702
# 
# 1:2                                              exp(coef) exp(-coef) lower .95 upper .95
# hepatitis_B_virus_core_antibody_quantification         1     0.9997    0.9987     1.002
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.522  (se = 0.036 )
# Likelihood ratio test= 0.14  on 1 df,   p=0.7
# Wald test            = 0.15  on 1 df,   p=0.7
# Score (logrank) test = 0.15  on 1 df,   p=0.7


### hepatitis_B_virus_pre_S1_antigen ###                                          
colnames(all_data)[which(colnames(all_data)=="Abnormality of hepatitis B virus pre-S1 antigen")] <- "hepatitis_B_virus_pre_S1_antigen"
all_data$hepatitis_B_virus_pre_S1_antigen <- factor(all_data$hepatitis_B_virus_pre_S1_antigen,levels = c("-","+"))
table(all_data$hepatitis_B_virus_pre_S1_antigen)
# -  + 
#   30 86   

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, hepatitis_B_virus_pre_S1_antigen, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ hepatitis_B_virus_pre_S1_antigen, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ hepatitis_B_virus_pre_S1_antigen, 
#         data = data_factor, id = patient_id)
# 
# n= 92, number of events= 44 
# (162 observations deleted due to missingness)
# 
# 1:2                                    coef exp(coef) se(coef)      z Pr(>|z|)
# hepatitis_B_virus_pre_S1_antigen+ -0.5882    0.5553   0.3192 -1.843   0.0654
# 
# 1:2                                 exp(coef) exp(-coef) lower .95 upper .95
# hepatitis_B_virus_pre_S1_antigen+    0.5553      1.801    0.2971     1.038
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.57  (se = 0.038 )
# Likelihood ratio test= 3.16  on 1 df,   p=0.08
# Wald test            = 3.4  on 1 df,   p=0.07
# Score (logrank) test = 3.49  on 1 df,   p=0.06


### GGT ###                                          
colnames(all_data)[which(colnames(all_data)=="GGT (IU/L)")] <- "GGT"
summary(all_data$GGT)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    1.00   19.05   39.00   79.73   79.00  552.50     183  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, GGT, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ GGT, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ GGT, data = data_factor, 
#         id = patient_id)
# 
# n= 41, number of events= 23 
# (213 observations deleted due to missingness)
# 
# 1:2       coef exp(coef) se(coef)     z Pr(>|z|)
# GGT 0.002356  1.002359 0.002009 1.173    0.241
# 
# 1:2   exp(coef) exp(-coef) lower .95 upper .95
# GGT     1.002     0.9976    0.9984     1.006
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.64  (se = 0.054 )
# Likelihood ratio test= 1.16  on 1 df,   p=0.3
# Wald test            = 1.38  on 1 df,   p=0.2
# Score (logrank) test = 1.42  on 1 df,   p=0.2


### TBIL ###                                          
colnames(all_data)[which(colnames(all_data)=="TBIL (umol/L)")] <- "TBIL"
summary(all_data$TBIL)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.00   10.30   14.50   20.83   20.30  443.40      97 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, TBIL, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TBIL, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TBIL, data = data_factor, 
#         id = patient_id)
# 
# n= 116, number of events= 57 
# (138 observations deleted due to missingness)
# 
# 1:2        coef exp(coef) se(coef)     z Pr(>|z|)
# TBIL 0.002953  1.002958 0.002211 1.336    0.182
# 
# 1:2    exp(coef) exp(-coef) lower .95 upper .95
# TBIL     1.003     0.9971    0.9986     1.007
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.474  (se = 0.04 )
# Likelihood ratio test= 1.26  on 1 df,   p=0.3
# Wald test            = 1.78  on 1 df,   p=0.2
# Score (logrank) test = 1.99  on 1 df,   p=0.2


### IBIL ###                                          
colnames(all_data)[which(colnames(all_data)=="IBIL (umol/L)")] <- "IBIL"
summary(all_data$IBIL)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.00    6.70    8.90   12.51   13.60  215.60     111 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, IBIL, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ IBIL, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ IBIL, data = data_factor, 
#         id = patient_id)
# 
# n= 108, number of events= 53 
# (146 observations deleted due to missingness)
# 
# 1:2        coef exp(coef) se(coef)     z Pr(>|z|)
# IBIL 0.005486  1.005501 0.004979 1.102    0.271
# 
# 1:2    exp(coef) exp(-coef) lower .95 upper .95
# IBIL     1.006     0.9945    0.9957     1.015
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.476  (se = 0.042 )
# Likelihood ratio test= 0.89  on 1 df,   p=0.3
# Wald test            = 1.21  on 1 df,   p=0.3
# Score (logrank) test = 1.32  on 1 df,   p=0.3


### ALT_abnormality ###                                          
colnames(all_data)[which(colnames(all_data)=="ALT abnormality")] <- "ALT_abnormality"
all_data$ALT_abnormality <- factor(all_data$ALT_abnormality,levels = c("normal","high"))
table(all_data$ALT_abnormality)
# normal   high 
# 106     32  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, ALT_abnormality, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ ALT_abnormality, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ ALT_abnormality, data = data_factor, 
#         id = patient_id)
# 
# n= 107, number of events= 52 
# (147 observations deleted due to missingness)
# 
# 1:2                     coef exp(coef) se(coef)     z Pr(>|z|)
# ALT_abnormalityhigh 0.5159    1.6752   0.3024 1.706    0.088
# 
# 1:2                   exp(coef) exp(-coef) lower .95 upper .95
# ALT_abnormalityhigh     1.675     0.5969    0.9261      3.03
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.536  (se = 0.034 )
# Likelihood ratio test= 2.71  on 1 df,   p=0.1
# Wald test            = 2.91  on 1 df,   p=0.09
# Score (logrank) test = 2.98  on 1 df,   p=0.08


### AST_abnormality ###                                          
colnames(all_data)[which(colnames(all_data)=="AST  abnormality")] <- "AST_abnormality"
all_data$AST_abnormality <- factor(all_data$AST_abnormality,levels = c("normal","high"))
table(all_data$AST_abnormality)
# normal   high 
# 91     47  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, AST_abnormality, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ AST_abnormality, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ AST_abnormality, data = data_factor, 
#         id = patient_id)
# 
# n= 107, number of events= 52 
# (147 observations deleted due to missingness)
# 
# 1:2                     coef exp(coef) se(coef)    z Pr(>|z|)
# AST_abnormalityhigh 0.6844    1.9825   0.2875 2.38   0.0173
# 
# 1:2                   exp(coef) exp(-coef) lower .95 upper .95
# AST_abnormalityhigh     1.983     0.5044     1.128     3.483
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.579  (se = 0.037 )
# Likelihood ratio test= 5.48  on 1 df,   p=0.02
# Wald test            = 5.67  on 1 df,   p=0.02
# Score (logrank) test = 5.87  on 1 df,   p=0.02

mylogit <- glm(Recurrence ~ AST_abnormality, data = data_factor,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ AST_abnormality, family = binomial, 
#       data = data_factor)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -1.269  -1.058  -1.058   1.302   1.302  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)          -0.2877     0.2118  -1.358    0.174
# AST_abnormalityhigh   0.5013     0.3619   1.385    0.166
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 190.84  on 137  degrees of freedom
# Residual deviance: 188.91  on 136  degrees of freedom
# (116 observations deleted due to missingness)
# AIC: 192.91
# 
# Number of Fisher Scoring iterations: 4


### AST_ALT_abnormality ###                                          
colnames(all_data)[which(colnames(all_data)=="AST/ALT abnormality")] <- "AST_ALT_abnormality"
all_data$AST_ALT_abnormality <- factor(all_data$AST_ALT_abnormality,levels = c("normal","high"))
table(all_data$AST_ALT_abnormality)
# normal   high 
# 75     45    

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, AST_ALT_abnormality, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ AST_ALT_abnormality, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ AST_ALT_abnormality, 
#         data = data_factor, id = patient_id)
# 
# n= 94, number of events= 44 
# (160 observations deleted due to missingness)
# 
# 1:2                         coef exp(coef) se(coef)     z Pr(>|z|)
# AST_ALT_abnormalityhigh 0.8595    2.3619   0.3074 2.796  0.00517
# 
# 1:2                       exp(coef) exp(-coef) lower .95 upper .95
# AST_ALT_abnormalityhigh     2.362     0.4234     1.293     4.314
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.615  (se = 0.04 )
# Likelihood ratio test= 7.74  on 1 df,   p=0.005
# Wald test            = 7.82  on 1 df,   p=0.005
# Score (logrank) test = 8.3  on 1 df,   p=0.004

mylogit <- glm(Recurrence ~ AST_ALT_abnormality, data = data_factor,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ AST_ALT_abnormality, family = binomial, 
#       data = data_factor)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.3132  -0.9888  -0.9888   1.0474   1.3785  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)              -0.4613     0.2371  -1.946   0.0517 .
# AST_ALT_abnormalityhigh   0.7750     0.3838   2.019   0.0435 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 165.52  on 119  degrees of freedom
# Residual deviance: 161.37  on 118  degrees of freedom
# (134 observations deleted due to missingness)
# AIC: 165.37
# 
# Number of Fisher Scoring iterations: 4


### ALP_abnormality ###                                          
colnames(all_data)[which(colnames(all_data)=="ALP abnormality")] <- "ALP_abnormality"
all_data$ALP_abnormality <- factor(all_data$ALP_abnormality,levels = c("normal","high"))
table(all_data$ALP_abnormality)
# normal   high 
# 90     38 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, ALP_abnormality, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ ALP_abnormality, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ ALP_abnormality, data = data_factor, 
#         id = patient_id)
# 
# n= 98, number of events= 48 
# (156 observations deleted due to missingness)
# 
# 1:2                     coef exp(coef) se(coef)    z Pr(>|z|)
# ALP_abnormalityhigh 0.3656    1.4414   0.3152 1.16    0.246
# 
# 1:2                   exp(coef) exp(-coef) lower .95 upper .95
# ALP_abnormalityhigh     1.441     0.6938    0.7771     2.673
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.539  (se = 0.037 )
# Likelihood ratio test= 1.28  on 1 df,   p=0.3
# Wald test            = 1.35  on 1 df,   p=0.2
# Score (logrank) test = 1.36  on 1 df,   p=0.2


### GGT_abnormality ###                                          
colnames(all_data)[which(colnames(all_data)=="GGT abnormality")] <- "GGT_abnormality"
all_data$GGT_abnormality <- factor(all_data$GGT_abnormality,levels = c("normal","high"))
table(all_data$GGT_abnormality)
# normal   high 
# 35     33 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, GGT_abnormality, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ GGT_abnormality, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ GGT_abnormality, data = data_factor, 
#         id = patient_id)
# 
# n= 38, number of events= 22 
# (216 observations deleted due to missingness)
# 
# 1:2                       coef exp(coef) se(coef)     z Pr(>|z|)
# GGT_abnormalityhigh 0.007055  1.007080 0.442910 0.016    0.987
# 
# 1:2                   exp(coef) exp(-coef) lower .95 upper .95
# GGT_abnormalityhigh     1.007      0.993    0.4227     2.399
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.498  (se = 0.06 )
# Likelihood ratio test= 0  on 1 df,   p=1
# Wald test            = 0  on 1 df,   p=1
# Score (logrank) test = 0  on 1 df,   p=1


### TBIL_abnormality ###                                          
colnames(all_data)[which(colnames(all_data)=="TBIL  abnormality")] <- "TBIL_abnormality"
all_data$TBIL_abnormality <- factor(all_data$TBIL_abnormality,levels = c("normal","high"))
table(all_data$TBIL_abnormality)
# normal   high 
# 104     33  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, TBIL_abnormality, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TBIL_abnormality, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TBIL_abnormality, data = data_factor, 
#         id = patient_id)
# 
# n= 106, number of events= 51 
# (148 observations deleted due to missingness)
# 
# 1:2                      coef exp(coef) se(coef)     z Pr(>|z|)
# TBIL_abnormalityhigh 0.3425    1.4084   0.3168 1.081     0.28
# 
# 1:2                    exp(coef) exp(-coef) lower .95 upper .95
# TBIL_abnormalityhigh     1.408       0.71     0.757      2.62
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.518  (se = 0.031 )
# Likelihood ratio test= 1.11  on 1 df,   p=0.3
# Wald test            = 1.17  on 1 df,   p=0.3
# Score (logrank) test = 1.18  on 1 df,   p=0.3


### DBIL_abnormality ###                                          
colnames(all_data)[which(colnames(all_data)=="DBIL abnormality")] <- "DBIL_abnormality"
all_data$DBIL_abnormality <- factor(all_data$DBIL_abnormality,levels = c("normal","high"))
table(all_data$DBIL_abnormality)
# normal   high 
# 85     52  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, DBIL_abnormality, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ DBIL_abnormality, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ DBIL_abnormality, data = data_factor, 
#         id = patient_id)
# 
# n= 106, number of events= 52 
# (148 observations deleted due to missingness)
# 
# 1:2                       coef exp(coef) se(coef)      z Pr(>|z|)
# DBIL_abnormalityhigh -0.2921    0.7467   0.3077 -0.949    0.342
# 
# 1:2                    exp(coef) exp(-coef) lower .95 upper .95
# DBIL_abnormalityhigh    0.7467      1.339    0.4085     1.365
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.534  (se = 0.035 )
# Likelihood ratio test= 0.94  on 1 df,   p=0.3
# Wald test            = 0.9  on 1 df,   p=0.3
# Score (logrank) test = 0.91  on 1 df,   p=0.3


### IBIL_abnormality ###                                          
colnames(all_data)[which(colnames(all_data)=="IBIL  abnormality")] <- "IBIL_abnormality"
all_data$IBIL_abnormality <- factor(all_data$IBIL_abnormality,levels = c("normal","high"))
table(all_data$IBIL_abnormality)
# normal   high 
# 102     35 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, IBIL_abnormality, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ IBIL_abnormality, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ IBIL_abnormality, data = data_factor, 
#         id = patient_id)
# 
# n= 106, number of events= 52 
# (148 observations deleted due to missingness)
# 
# 1:2                      coef exp(coef) se(coef)     z Pr(>|z|)
# IBIL_abnormalityhigh 0.1571    1.1702   0.3147 0.499    0.617
# 
# 1:2                    exp(coef) exp(-coef) lower .95 upper .95
# IBIL_abnormalityhigh      1.17     0.8546    0.6315     2.168
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.5  (se = 0.03 )
# Likelihood ratio test= 0.24  on 1 df,   p=0.6
# Wald test            = 0.25  on 1 df,   p=0.6
# Score (logrank) test = 0.25  on 1 df,   p=0.6


#### group 4: MR ####

### ADC_signal ###                                          
colnames(all_data)[which(colnames(all_data)=="ADC (Apparent Diffusion Coefficient) signal")] <- "ADC_signal"
all_data$ADC_signal <- factor(all_data$ADC_signal,levels = c("low","high"))
table(all_data$ADC_signal)
# low high 
# 25    3 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, ADC_signal, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ ADC_signal, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ ADC_signal, data = data_factor, 
#         id = patient_id)
# 
# n= 23, number of events= 7 
# (231 observations deleted due to missingness)
# 
# 1:2               coef exp(coef) se(coef)     z Pr(>|z|)
# ADC_signalhigh 0.802     2.230    1.156 0.694    0.488
# 
# 1:2              exp(coef) exp(-coef) lower .95 upper .95
# ADC_signalhigh      2.23     0.4484    0.2314     21.49
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.541  (se = 0.087 )
# Likelihood ratio test= 0.42  on 1 df,   p=0.5
# Wald test            = 0.48  on 1 df,   p=0.5
# Score (logrank) test = 0.51  on 1 df,   p=0.5


### DWI_signal ###                                          
colnames(all_data)[which(colnames(all_data)=="DWI (Diffusion Weighted Imaging) signal")] <- "DWI_signal"
all_data$DWI_signal <- factor(all_data$DWI_signal,levels = c("limited","high"))
table(all_data$DWI_signal)
# limited    high 
# 1      38 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, DWI_signal, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ DWI_signal, data = data_factor, id = patient_id)
# Warning message:
#   In fitter(X, Y, istrat, offset, init, control, weights = weights,  :
#               Loglik converged before variable  1 ; coefficient may be infinite.  
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ DWI_signal, data = data_factor, 
#         id = patient_id)
# 
# n= 31, number of events= 12 
# (223 observations deleted due to missingness)
# 
# 1:2                   coef exp(coef)  se(coef)     z Pr(>|z|)
# DWI_signalhigh 1.819e+01 7.926e+07 9.518e+03 0.002    0.998
# 
# 1:2              exp(coef) exp(-coef) lower .95 upper .95
# DWI_signalhigh  79255486  1.262e-08         0       Inf
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.526  (se = 0.025 )
# Likelihood ratio test= 1.65  on 1 df,   p=0.2
# Wald test            = 0  on 1 df,   p=1
# Score (logrank) test = 0.87  on 1 df,   p=0.4


### tumor_border ###                                          
colnames(all_data)[which(colnames(all_data)=="Tumor border")] <- "tumor_border"
all_data$tumor_border <- factor(all_data$tumor_border,levels = c("blurry","clear"))
table(all_data$tumor_border)
# blurry  clear 
# 10     19 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, tumor_border, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ tumor_border, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ tumor_border, data = data_factor, 
#         id = patient_id)
# 
# n= 24, number of events= 11 
# (230 observations deleted due to missingness)
# 
# 1:2                    coef exp(coef) se(coef)      z Pr(>|z|)
# tumor_borderclear -0.4383    0.6451   0.6322 -0.693    0.488
# 
# 1:2                 exp(coef) exp(-coef) lower .95 upper .95
# tumor_borderclear    0.6451       1.55    0.1869     2.227
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.609  (se = 0.086 )
# Likelihood ratio test= 0.46  on 1 df,   p=0.5
# Wald test            = 0.48  on 1 df,   p=0.5
# Score (logrank) test = 0.49  on 1 df,   p=0.5


### dynamic_enhanced_arterial_phase ###                                          
colnames(all_data)[which(colnames(all_data)=="Dynamic enhanced arterial phase")] <- "dynamic_enhanced_arterial_phase"
all_data$dynamic_enhanced_arterial_phase <- factor(all_data$dynamic_enhanced_arterial_phase,levels = c("hyperenhancement","uneven enhancement"))
table(all_data$dynamic_enhanced_arterial_phase)
# hyperenhancement uneven enhancement 
# 23                 18  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, dynamic_enhanced_arterial_phase, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ dynamic_enhanced_arterial_phase, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ dynamic_enhanced_arterial_phase, 
#         data = data_factor, id = patient_id)
# 
# n= 32, number of events= 13 
# (222 observations deleted due to missingness)
# 
# 1:2                                                   coef exp(coef) se(coef)     z Pr(>|z|)
# dynamic_enhanced_arterial_phaseuneven enhancement 2.0585    7.8345   0.6869 2.997  0.00273
# 
# 1:2                                                 exp(coef) exp(-coef) lower .95 upper .95
# dynamic_enhanced_arterial_phaseuneven enhancement     7.834     0.1276     2.039     30.11
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.759  (se = 0.066 )
# Likelihood ratio test= 10.83  on 1 df,   p=0.001
# Wald test            = 8.98  on 1 df,   p=0.003
# Score (logrank) test = 11.96  on 1 df,   p=5e-04

mylogit <- glm(Recurrence ~ dynamic_enhanced_arterial_phase, data = data_factor,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ dynamic_enhanced_arterial_phase, family = binomial, 
#       data = data_factor)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.4823  -0.7775  -0.7775   0.9005   1.6394  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)                                        -1.0415     0.4749  -2.193   0.0283 *
#   dynamic_enhanced_arterial_phaseuneven enhancement   1.7346     0.6896   2.516   0.0119 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 56.227  on 40  degrees of freedom
# Residual deviance: 49.317  on 39  degrees of freedom
# (213 observations deleted due to missingness)
# AIC: 53.317
# 
# Number of Fisher Scoring iterations: 4


### portal_venous_phase ###                                          
colnames(all_data)[which(colnames(all_data)=="Portal venous phase")] <- "portal_venous_phase"
all_data$portal_venous_phase <- factor(all_data$portal_venous_phase,levels = c("enhancement","attenuation"))
table(all_data$portal_venous_phase)
# enhancement attenuation 
# 4          37

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, portal_venous_phase, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ portal_venous_phase, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ portal_venous_phase, 
#         data = data_factor, id = patient_id)
# 
# n= 33, number of events= 13 
# (221 observations deleted due to missingness)
# 
# 1:2                                 coef exp(coef) se(coef)      z Pr(>|z|)
# portal_venous_phaseattenuation -0.7889    0.4544   0.8034 -0.982    0.326
# 
# 1:2                              exp(coef) exp(-coef) lower .95 upper .95
# portal_venous_phaseattenuation    0.4544      2.201   0.09409     2.194
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.549  (se = 0.061 )
# Likelihood ratio test= 0.83  on 1 df,   p=0.4
# Wald test            = 0.96  on 1 df,   p=0.3
# Score (logrank) test = 1.01  on 1 df,   p=0.3


### delayed_phase ###                                          
colnames(all_data)[which(colnames(all_data)=="Delayed phase")] <- "delayed_phase"
all_data$delayed_phase <- factor(all_data$delayed_phase,levels = c("enhancement","attenuation"))
table(all_data$delayed_phase)
# enhancement attenuation 
# 1          39

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, delayed_phase, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ delayed_phase, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ delayed_phase, data = data_factor, 
#         id = patient_id)
# 
# n= 32, number of events= 13 
# (222 observations deleted due to missingness)
# 
# 1:2                           coef exp(coef) se(coef)      z Pr(>|z|)
# delayed_phaseattenuation -1.3193    0.2673   1.0703 -1.233    0.218
# 
# 1:2                        exp(coef) exp(-coef) lower .95 upper .95
# delayed_phaseattenuation    0.2673      3.741   0.03281     2.178
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.528  (se = 0.029 )
# Likelihood ratio test= 1.1  on 1 df,   p=0.3
# Wald test            = 1.52  on 1 df,   p=0.2
# Score (logrank) test = 1.75  on 1 df,   p=0.2


### lymph_node ###                                          
colnames(all_data)[which(colnames(all_data)=="Lymph node")] <- "lymph_node"
all_data$lymph_node <- factor(all_data$lymph_node,levels = c("unenlarged","enlarged"))
table(all_data$lymph_node)
# unenlarged   enlarged 
# 33         9  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, lymph_node, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ lymph_node, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ lymph_node, data = data_factor, 
#         id = patient_id)
# 
# n= 33, number of events= 14 
# (221 observations deleted due to missingness)
# 
# 1:2                     coef exp(coef) se(coef)      z Pr(>|z|)
# lymph_nodeenlarged -0.3107    0.7329   1.0626 -0.292     0.77
# 
# 1:2                  exp(coef) exp(-coef) lower .95 upper .95
# lymph_nodeenlarged    0.7329      1.364   0.09133     5.882
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.53  (se = 0.035 )
# Likelihood ratio test= 0.09  on 1 df,   p=0.8
# Wald test            = 0.09  on 1 df,   p=0.8
# Score (logrank) test = 0.09  on 1 df,   p=0.8


#### group 5: CT ####

### shadow ###                                          
colnames(all_data)[which(colnames(all_data)=="Shadow")] <- "shadow"
all_data$shadow <- factor(all_data$shadow,levels = c("low-density","high-density"))
table(all_data$shadow)
# low-density high-density 
# 54        11  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, shadow, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ shadow, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ shadow, data = data_factor, 
#         id = patient_id)
# 
# n= 49, number of events= 24 
# (205 observations deleted due to missingness)
# 
# 1:2                    coef exp(coef) se(coef)     z Pr(>|z|)
# shadowhigh-density 0.3086    1.3615   0.5128 0.602    0.547
# 
# 1:2                  exp(coef) exp(-coef) lower .95 upper .95
# shadowhigh-density     1.361     0.7345    0.4983      3.72
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.531  (se = 0.048 )
# Likelihood ratio test= 0.34  on 1 df,   p=0.6
# Wald test            = 0.36  on 1 df,   p=0.5
# Score (logrank) test = 0.36  on 1 df,   p=0.5


### shadow_border ###                                          
colnames(all_data)[which(colnames(all_data)=="Shadow border")] <- "shadow_border"
all_data$shadow_border <- factor(all_data$shadow_border,levels = c("blurry","clear"))
table(all_data$shadow_border)
# blurry  clear 
# 34     15  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, shadow_border, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ shadow_border, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ shadow_border, data = data_factor, 
#         id = patient_id)
# 
# n= 36, number of events= 18 
# (218 observations deleted due to missingness)
# 
# 1:2                     coef exp(coef) se(coef)      z Pr(>|z|)
# shadow_borderclear -0.3781    0.6852   0.5348 -0.707     0.48
# 
# 1:2                  exp(coef) exp(-coef) lower .95 upper .95
# shadow_borderclear    0.6852      1.459    0.2402     1.954
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.549  (se = 0.058 )
# Likelihood ratio test= 0.52  on 1 df,   p=0.5
# Wald test            = 0.5  on 1 df,   p=0.5
# Score (logrank) test = 0.51  on 1 df,   p=0.5


### lymphadenectasis ###                                          
colnames(all_data)[which(colnames(all_data)=="Lymphadenectasis")] <- "lymphadenectasis"
all_data$lymphadenectasis <- factor(all_data$lymphadenectasis,levels = c("no","yes"))
table(all_data$lymphadenectasis)
# no yes 
# 51   8 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, lymphadenectasis, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ lymphadenectasis, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ lymphadenectasis, data = data_factor, 
#         id = patient_id)
# 
# n= 46, number of events= 22 
# (208 observations deleted due to missingness)
# 
# 1:2                     coef exp(coef) se(coef)     z Pr(>|z|)
# lymphadenectasisyes 1.3937    4.0297   0.4869 2.862   0.0042
# 
# 1:2                   exp(coef) exp(-coef) lower .95 upper .95
# lymphadenectasisyes      4.03     0.2482     1.552     10.46
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.611  (se = 0.05 )
# Likelihood ratio test= 6.5  on 1 df,   p=0.01
# Wald test            = 8.19  on 1 df,   p=0.004
# Score (logrank) test = 9.58  on 1 df,   p=0.002

mylogit <- glm(Recurrence ~ lymphadenectasis, data = data_factor,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ lymphadenectasis, family = binomial, 
#       data = data_factor)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -1.665  -1.030  -1.030   1.332   1.332  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)          -0.3567     0.2845  -1.254   0.2100  
# lymphadenectasisyes   1.4553     0.8646   1.683   0.0924 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 81.367  on 58  degrees of freedom
# Residual deviance: 78.102  on 57  degrees of freedom
# (195 observations deleted due to missingness)
# AIC: 82.102
# 
# Number of Fisher Scoring iterations: 4


#### group 6: ultrasound ####

### internal_echo ###                                          
colnames(all_data)[which(colnames(all_data)=="Internal echo")] <- "internal_echo"
all_data$internal_echo <- factor(all_data$internal_echo,levels = c("low-level","high-level"))
table(all_data$internal_echo)
# low-level high-level 
# 29         17 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, internal_echo, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ internal_echo, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ internal_echo, data = data_factor, 
#         id = patient_id)
# 
# n= 36, number of events= 20 
# (218 observations deleted due to missingness)
# 
# 1:2                          coef exp(coef) se(coef)      z Pr(>|z|)
# internal_echohigh-level -0.5855    0.5568   0.5646 -1.037      0.3
# 
# 1:2                       exp(coef) exp(-coef) lower .95 upper .95
# internal_echohigh-level    0.5568      1.796    0.1841     1.684
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.542  (se = 0.062 )
# Likelihood ratio test= 1.19  on 1 df,   p=0.3
# Wald test            = 1.08  on 1 df,   p=0.3
# Score (logrank) test = 1.11  on 1 df,   p=0.3


### internal_echo_border ###                                          
colnames(all_data)[which(colnames(all_data)=="Internal echo border")] <- "internal_echo_border"
all_data$internal_echo_border <- factor(all_data$internal_echo_border,levels = c("blurry","clear"))
table(all_data$internal_echo_border)
# blurry  clear 
# 15     36  

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, internal_echo_border, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ internal_echo_border, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ internal_echo_border, 
#         data = data_factor, id = patient_id)
# 
# n= 41, number of events= 24 
# (213 observations deleted due to missingness)
# 
# 1:2                            coef exp(coef) se(coef)      z Pr(>|z|)
# internal_echo_borderclear -0.6222    0.5368   0.4398 -1.415    0.157
# 
# 1:2                         exp(coef) exp(-coef) lower .95 upper .95
# internal_echo_borderclear    0.5368      1.863    0.2267     1.271
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.566  (se = 0.051 )
# Likelihood ratio test= 1.86  on 1 df,   p=0.2
# Wald test            = 2  on 1 df,   p=0.2
# Score (logrank) test = 2.07  on 1 df,   p=0.2


### blood_flow_characteristics ###                                          
colnames(all_data)[which(colnames(all_data)=="Blood flow characteristics")] <- "blood_flow_characteristics"
all_data$blood_flow_characteristics <- factor(all_data$blood_flow_characteristics,levels = c("no","yes"))
table(all_data$blood_flow_characteristics)
# no yes 
# 16  40 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, blood_flow_characteristics, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ blood_flow_characteristics, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
temp_df <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int); colnames(temp_df)[1] <- "X_name"
temp_df
cox_res_df <- rbind(cox_res_df,temp_df)  ##store results.
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ blood_flow_characteristics, 
#         data = data_factor, id = patient_id)
# 
# n= 44, number of events= 24 
# (210 observations deleted due to missingness)
# 
# 1:2                               coef exp(coef) se(coef)     z Pr(>|z|)
# blood_flow_characteristicsyes 1.0306    2.8027   0.7408 1.391    0.164
# 
# 1:2                             exp(coef) exp(-coef) lower .95 upper .95
# blood_flow_characteristicsyes     2.803     0.3568    0.6562     11.97
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.558  (se = 0.042 )
# Likelihood ratio test= 2.57  on 1 df,   p=0.1
# Wald test            = 1.94  on 1 df,   p=0.2
# Score (logrank) test = 2.11  on 1 df,   p=0.1


##post-process
cox_res_df$`Pr(>|z|)` <- as.numeric(cox_res_df$`Pr(>|z|)`)
cox_res_df$coef <- as.numeric(cox_res_df$coef)
cox_res_df_sig <- cox_res_df[which(cox_res_df$`Pr(>|z|)`<0.05),]
cox_res_df_sig_sort <- cox_res_df_sig[order(cox_res_df_sig$`Pr(>|z|)`),]

#write.table(cox_res_df_sig_sort, file = "/data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/cox_univariate_regression.txt",sep="\t",na="NA",row.names=F, col.names=T)

unlist(cox_res_df_sig_sort$X_name)  #18
# "TIMES"                          "tumor_greatest_dimension"                                      "tumor_burden"                                   "tumor_area_size" 
# "macrovascular_invasionyes"                  "tumor_differentiation_gradeG2-G3"                                 "satellite_lesion+" "dynamic_enhanced_arterial_phaseuneven enhancement" 
# "Alpha_fetoprotein"                               "lymphadenectasisyes"                           "AST_ALT_abnormalityhigh"                                 "MVI_M0_vs_M1M2yes" 
# "BCLC_A_vs_BChigh"                               "AST_abnormalityhigh"                                "capsular_invasion+"                                             "Ki_67" 
# "TNM_2level≥II"                              "neoplastic_thrombus+" 

#######Stage 2: Multi-variate regression analysis####

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

res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES + tumor_greatest_dimension + tumor_burden + tumor_area_size + macrovascular_invasion + tumor_differentiation_grade + satellite_lesion + dynamic_enhanced_arterial_phase + 
                   Alpha_fetoprotein + lymphadenectasis + AST_ALT_abnormality + MVI_M0_vs_M1M2 + BCLC_A_vs_BC + AST_abnormality + capsular_invasion + Ki_67 + TNM_2level + neoplastic_thrombus, data = data_multi_fill, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES + tumor_greatest_dimension + 
#           tumor_burden + tumor_area_size + macrovascular_invasion + 
#           tumor_differentiation_grade + satellite_lesion + dynamic_enhanced_arterial_phase + 
#           Alpha_fetoprotein + lymphadenectasis + AST_ALT_abnormality + 
#           MVI_M0_vs_M1M2 + BCLC_A_vs_BC + AST_abnormality + capsular_invasion + 
#           Ki_67 + TNM_2level + neoplastic_thrombus, data = data_multi_fill, 
#         id = patient_id)
# 
# n= 127, number of events= 55 
# (127 observations deleted due to missingness)
# 
# 1:2                                                       coef  exp(coef)   se(coef)      z Pr(>|z|)
# TIMES                                              4.844e+00  1.269e+02  7.326e-01  6.611 3.82e-11
# tumor_greatest_dimension                           2.574e-01  1.294e+00  1.435e-01  1.794 0.072883
# tumor_burden                                      -1.554e-03  9.984e-01  4.043e-03 -0.384 0.700639
# tumor_area_size                                   -1.128e-01  8.934e-01  2.384e-01 -0.473 0.636250
# macrovascular_invasionyes                         -3.337e-01  7.163e-01  5.859e-01 -0.570 0.568957
# tumor_differentiation_gradeG2-G3                  -1.878e-02  9.814e-01  5.108e-01 -0.037 0.970671
# tumor_differentiation_grade≥G3                     1.660e-01  1.181e+00  4.395e-01  0.378 0.705606
# satellite_lesion+                                  1.439e+00  4.216e+00  4.259e-01  3.378 0.000729
# dynamic_enhanced_arterial_phaseuneven enhancement  8.854e-01  2.424e+00  5.393e-01  1.642 0.100662
# Alpha_fetoprotein                                 -5.747e-07  1.000e+00  4.911e-07 -1.170 0.241956
# lymphadenectasisyes                                6.535e-01  1.922e+00  7.072e-01  0.924 0.355499
# AST_ALT_abnormalityhigh                            7.323e-01  2.080e+00  3.605e-01  2.031 0.042209
# MVI_M0_vs_M1M2yes                                  1.337e-01  1.143e+00  6.274e-01  0.213 0.831287
# BCLC_A_vs_BChigh                                   6.108e-01  1.842e+00  6.541e-01  0.934 0.350390
# AST_abnormalityhigh                                3.137e-01  1.368e+00  3.844e-01  0.816 0.414487
# capsular_invasion+                                 9.270e-01  2.527e+00  3.745e-01  2.475 0.013324
# Ki_67                                              3.472e-01  1.415e+00  1.012e+00  0.343 0.731510
# TNM_2level≥II                                     -5.460e-01  5.793e-01  5.236e-01 -1.043 0.297046
# neoplastic_thrombus+                              -2.493e-02  9.754e-01  4.749e-01 -0.053 0.958128
# 
# 1:2                                                 exp(coef) exp(-coef) lower .95 upper .95
# TIMES                                              126.9217   0.007879   30.1932   533.536
# tumor_greatest_dimension                             1.2935   0.773079    0.9764     1.714
# tumor_burden                                         0.9984   1.001556    0.9906     1.006
# tumor_area_size                                      0.8934   1.119357    0.5599     1.425
# macrovascular_invasionyes                            0.7163   1.396137    0.2272     2.258
# tumor_differentiation_gradeG2-G3                     0.9814   1.018959    0.3606     2.671
# tumor_differentiation_grade≥G3                       1.1806   0.847009    0.4989     2.794
# satellite_lesion+                                    4.2155   0.237218    1.8295     9.713
# dynamic_enhanced_arterial_phaseuneven enhancement    2.4239   0.412561    0.8423     6.976
# Alpha_fetoprotein                                    1.0000   1.000001    1.0000     1.000
# lymphadenectasisyes                                  1.9222   0.520244    0.4806     7.687
# AST_ALT_abnormalityhigh                              2.0799   0.480791    1.0261     4.216
# MVI_M0_vs_M1M2yes                                    1.1430   0.874875    0.3342     3.910
# BCLC_A_vs_BChigh                                     1.8420   0.542892    0.5111     6.639
# AST_abnormalityhigh                                  1.3684   0.730766    0.6442     2.907
# capsular_invasion+                                   2.5268   0.395752    1.2128     5.265
# Ki_67                                                1.4150   0.706689    0.1948    10.280
# TNM_2level≥II                                        0.5793   1.726346    0.2076     1.616
# neoplastic_thrombus+                                 0.9754   1.025248    0.3845     2.474
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.874  (se = 0.03 )
# Likelihood ratio test= 117.3  on 19 df,   p=4e-16
# Wald test            = 84.04  on 19 df,   p=4e-10
# Score (logrank) test = 143.6  on 19 df,   p=<2e-16

mylogit <- glm(Recurrence ~ TIMES + tumor_greatest_dimension + tumor_burden + tumor_area_size + macrovascular_invasion + tumor_differentiation_grade + satellite_lesion + dynamic_enhanced_arterial_phase + 
                 Alpha_fetoprotein + lymphadenectasis + AST_ALT_abnormality + MVI_M0_vs_M1M2 + BCLC_A_vs_BC + AST_abnormality + capsular_invasion + Ki_67 + TNM_2level + neoplastic_thrombus, data = data_multi_fill,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES + tumor_greatest_dimension + 
#         tumor_burden + tumor_area_size + macrovascular_invasion + 
#         tumor_differentiation_grade + satellite_lesion + dynamic_enhanced_arterial_phase + 
#         Alpha_fetoprotein + lymphadenectasis + AST_ALT_abnormality + 
#         MVI_M0_vs_M1M2 + BCLC_A_vs_BC + AST_abnormality + capsular_invasion + 
#         Ki_67 + TNM_2level + neoplastic_thrombus, family = binomial, 
#       data = data_multi_fill)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.1105  -0.6824  -0.3523   0.6332   2.3348  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                       -4.144e+00  9.159e-01  -4.524 6.06e-06 ***
#   TIMES                                              5.209e+00  9.535e-01   5.463 4.69e-08 ***
#   tumor_greatest_dimension                           3.860e-01  2.000e-01   1.930   0.0536 .  
# tumor_burden                                      -2.080e-03  7.476e-03  -0.278   0.7809    
# tumor_area_size                                   -1.404e-01  4.080e-01  -0.344   0.7308    
# macrovascular_invasionyes                          4.725e-01  9.315e-01   0.507   0.6120    
# tumor_differentiation_gradeG2-G3                   3.797e-01  7.126e-01   0.533   0.5942    
# tumor_differentiation_grade≥G3                     5.126e-01  5.794e-01   0.885   0.3762    
# satellite_lesion+                                 -2.900e-01  6.348e-01  -0.457   0.6478    
# dynamic_enhanced_arterial_phaseuneven enhancement  1.883e-01  8.179e-01   0.230   0.8179    
# Alpha_fetoprotein                                  2.428e-06  5.618e-06   0.432   0.6657    
# lymphadenectasisyes                                1.106e+00  1.600e+00   0.691   0.4893    
# AST_ALT_abnormalityhigh                            5.963e-01  5.375e-01   1.109   0.2672    
# MVI_M0_vs_M1M2yes                                  1.650e-01  7.652e-01   0.216   0.8293    
# BCLC_A_vs_BChigh                                  -3.408e-01  8.510e-01  -0.400   0.6888    
# AST_abnormalityhigh                               -3.197e-01  5.646e-01  -0.566   0.5712    
# capsular_invasion+                                 3.302e-01  5.371e-01   0.615   0.5387    
# Ki_67                                              2.516e-01  1.590e+00   0.158   0.8743    
# TNM_2level≥II                                     -4.261e-01  8.253e-01  -0.516   0.6057    
# neoplastic_thrombus+                               4.885e-01  7.441e-01   0.657   0.5115    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 223.75  on 162  degrees of freedom
# Residual deviance: 144.20  on 143  degrees of freedom
# (91 observations deleted due to missingness)
# AIC: 184.2
# 
# Number of Fisher Scoring iterations: 7


#######Stage 3: generate forest plot####
library(metafor)

## mutivariate regression ##
res_out <- summary(res.cox)
res_out_all <- cbind(cbind(row.names(res_out$cmap)[c(-1)],res_out$coefficients), res_out$conf.int)
colnames(res_out_all)[1] <- "X_name"
res_out_all <- as.data.frame(res_out_all, stringsAsFactors = F)
res_out_all$X_name <- as.character(res_out_all$X_name)
res_out_all[,-1] <- lapply(res_out_all[,-1],as.numeric)

res_out_all$coef_lw95 <- res_out_all$coef - 1.96*res_out_all$`se(coef)`
res_out_all$coef_up95 <- res_out_all$coef + 1.96*res_out_all$`se(coef)`
res_out_all$coef_lw95_2 <- log(res_out_all$`lower .95`)
res_out_all$coef_up95_2 <- log(res_out_all$`upper .95`)

summary(abs((res_out_all$coef_lw95-res_out_all$coef_lw95_2)/res_out_all$coef_lw95_2))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 7.743e-06 1.508e-05 2.228e-05 6.354e-05 3.493e-05 5.038e-04 
summary(abs((res_out_all$coef_up95-res_out_all$coef_up95_2)/res_out_all$coef_up95_2))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 4.202e-06 9.797e-06 1.541e-05 1.730e-05 2.087e-05 4.559e-05  

res_out_all <- res_out_all[order(res_out_all$`Pr(>|z|)`),]

p.vals <- res_out_all$`Pr(>|z|)`
p.vals <- ifelse(p.vals < 0.01, 
                 format(p.vals,digits = 3,scientific = TRUE,trim = TRUE),
                 format(round(p.vals, 2), nsmall=2, trim=TRUE))
#pdf('multi_cox_nw.pdf',width = 9,height = 8)
forest(res_out_all$coef, ci.lb = res_out_all$coef_lw95, ci.ub = res_out_all$coef_up95, slab = gsub("_"," ",res_out_all$X_name), 
       ilab= p.vals, ilab.pos=4, ilab.xpos = 2, xlab = "Effect Size", xlim = c(-2,7), main = "Multivariate Analysis", shade=TRUE, top=1, cex=1.2, lwd=2)
#dev.off()
##save as /data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/Fig3g_multivariate_cox.pdf : forest plot for Multivariate analysis


## Univariate regression ##
cox_res_df_sig_sort[,-1] <- lapply(cox_res_df_sig_sort[,-1],as.numeric)

cox_res_df_sig_sort$coef_lw95 <- cox_res_df_sig_sort$coef - 1.96*cox_res_df_sig_sort$`se(coef)`
cox_res_df_sig_sort$coef_up95 <- cox_res_df_sig_sort$coef + 1.96*cox_res_df_sig_sort$`se(coef)`
cox_res_df_sig_sort$coef_lw95_2 <- log(cox_res_df_sig_sort$`lower .95`)
cox_res_df_sig_sort$coef_up95_2 <- log(cox_res_df_sig_sort$`upper .95`)

summary(abs((cox_res_df_sig_sort$coef_lw95-cox_res_df_sig_sort$coef_lw95_2)/cox_res_df_sig_sort$coef_lw95_2))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 6.424e-06 1.899e-05 3.890e-05 7.158e-05 7.850e-05 3.473e-04 
summary(abs((cox_res_df_sig_sort$coef_up95-cox_res_df_sig_sort$coef_up95_2)/cox_res_df_sig_sort$coef_up95_2))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 3.781e-06 6.192e-06 7.431e-06 7.033e-06 8.201e-06 8.951e-06 

p.vals <- cox_res_df_sig_sort$`Pr(>|z|)`
p.vals <- ifelse(p.vals < 0.01, 
                 format(p.vals,digits = 3,scientific = TRUE,trim = TRUE),
                 format(round(p.vals, 3), nsmall=3, trim=TRUE))

forest(cox_res_df_sig_sort$coef, ci.lb = cox_res_df_sig_sort$coef_lw95, ci.ub = cox_res_df_sig_sort$coef_up95, slab = gsub("_"," ",cox_res_df_sig_sort$X_name), 
       ilab= p.vals, ilab.pos=4, ilab.xpos = 2, xlab = "Effect Size", xlim = c(-2,6), main = "Single-variate Analysis", shade=TRUE, top=1, cex=1.2, lwd=2)
##save as /data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/ExFig5d_univariate_cox.pdf : forest plot for univariate analysis

#save.image("/code/3_Comparison_between_clinical_markers_and_TIMES/results_univariate_multivariate_regression.RData")
