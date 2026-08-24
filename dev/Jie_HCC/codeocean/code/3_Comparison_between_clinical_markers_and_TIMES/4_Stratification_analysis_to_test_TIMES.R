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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] dplyr_1.1.4     survival_3.1-12
# 
# loaded via a namespace (and not attached):
#   [1] lattice_0.20-45   fansi_1.0.6       withr_3.0.0       utf8_1.2.4        grid_4.2.1        R6_2.5.1          lifecycle_1.0.4   magrittr_2.0.3    pillar_1.9.0     
# [10] rlang_1.1.3       cli_3.6.2         rstudioapi_0.16.0 Matrix_1.6-5      vctrs_0.6.5       generics_0.1.3    splines_4.2.1     tools_4.2.1       glue_1.7.0       
# [19] parallel_4.2.1    compiler_4.2.1    pkgconfig_2.0.3   tidyselect_1.2.1  tibble_3.2.1  

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


##1. BCLC 2 level: A vs B+C
all_data$`BCLC staging` <-factor(all_data$`BCLC staging`,levels = c("A","B","C"))
table(all_data$`BCLC staging`)
# A  B  C 
# 26 19 85 

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


##stratified test: BCLC=A group
all_data_no <- all_data[which(all_data$BCLC_A_vs_BC=="low"),]
data_BCLC_A_vs_BC_no  <- all_data_no %>%
  dplyr::select(DFS, Recurrence, BCLC_A_vs_BC, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_BCLC_A_vs_BC_no, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_BCLC_A_vs_BC_no, 
#         id = patient_id)
# 
# n= 21, number of events= 6 
# (5 observations deleted due to missingness)
# 
# 1:2         coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES    7.105  1218.180    3.287 2.161   0.0307
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES      1218  0.0008209     1.939    765412
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.758  (se = 0.159 )
# Likelihood ratio test= 5.92  on 1 df,   p=0.01
# Wald test            = 4.67  on 1 df,   p=0.03
# Score (logrank) test = 5.7  on 1 df,   p=0.02

mylogit <- glm(Recurrence ~ TIMES, data = data_BCLC_A_vs_BC_no,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_BCLC_A_vs_BC_no)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.4457  -0.8370  -0.5491   0.8511   2.3593  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   -3.023      1.257  -2.405   0.0162 *
#   TIMES          5.741      2.614   2.196   0.0281 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 33.542  on 25  degrees of freedom
# Residual deviance: 27.212  on 24  degrees of freedom
# AIC: 31.212
# 
# Number of Fisher Scoring iterations: 4

##stratified test: BCLC=B+C group
all_data_yes <- all_data[which(all_data$BCLC_A_vs_BC=="high"),]
data_BCLC_A_vs_BC_yes  <- all_data_yes %>%
  dplyr::select(DFS, Recurrence, BCLC_A_vs_BC, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_BCLC_A_vs_BC_yes, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_BCLC_A_vs_BC_yes, 
#         id = patient_id)
# 
# n= 79, number of events= 42 
# (25 observations deleted due to missingness)
# 
# 1:2         coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES   6.1879  486.8045   0.9569 6.466    1e-10
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES     486.8   0.002054     74.61      3176
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.84  (se = 0.028 )
# Likelihood ratio test= 70.02  on 1 df,   p=<2e-16
# Wald test            = 41.81  on 1 df,   p=1e-10
# Score (logrank) test = 65.52  on 1 df,   p=6e-16

mylogit <- glm(Recurrence ~ TIMES, data = data_BCLC_A_vs_BC_yes,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_BCLC_A_vs_BC_yes)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.4923  -0.4390  -0.2541   0.4798   2.5229  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -3.7529     0.7158  -5.243 1.58e-07 ***
#   TIMES         7.6801     1.3350   5.753 8.77e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 144.136  on 103  degrees of freedom
# Residual deviance:  78.493  on 102  degrees of freedom
# AIC: 82.493
# 
# Number of Fisher Scoring iterations: 5


##2. MVI (Microvascular invasion) 
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


##stratified test: no group
all_data_no <- all_data[which(all_data$MVI_M0_vs_M1M2=="no"),] #35
data_MVI_M0_vs_M1M2_no  <- all_data_no %>%
  dplyr::select(DFS, Recurrence, MVI_M0_vs_M1M2, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_MVI_M0_vs_M1M2_no, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_MVI_M0_vs_M1M2_no, 
#         id = patient_id)
# 
# n= 29, number of events= 8 
# (6 observations deleted due to missingness)
# 
# 1:2         coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES    8.507  4950.524    2.553 3.332 0.000861
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES      4951   0.000202     33.23    737451
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.891  (se = 0.07 )
# Likelihood ratio test= 16.7  on 1 df,   p=4e-05
# Wald test            = 11.1  on 1 df,   p=9e-04
# Score (logrank) test = 16.57  on 1 df,   p=5e-05
mylogit <- glm(Recurrence ~ TIMES, data = data_MVI_M0_vs_M1M2_no,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_MVI_M0_vs_M1M2_no)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.5639  -0.6053  -0.3253   0.6504   2.5551  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)   -3.599      1.149  -3.132  0.00174 **
#   TIMES          7.066      2.335   3.025  0.00248 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 45.004  on 34  degrees of freedom
# Residual deviance: 30.638  on 33  degrees of freedom
# AIC: 34.638
# 
# Number of Fisher Scoring iterations: 5

##stratified test: yes group
all_data_yes <- all_data[which(all_data$MVI_M0_vs_M1M2=="yes"),] #75
data_MVI_M0_vs_M1M2_yes  <- all_data_yes %>%
  dplyr::select(DFS, Recurrence, MVI_M0_vs_M1M2, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_MVI_M0_vs_M1M2_yes, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_MVI_M0_vs_M1M2_yes, 
#         id = patient_id)
# 
# n= 57, number of events= 31 
# (18 observations deleted due to missingness)
# 
# 1:2         coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES   5.0271  152.4868   0.9116 5.515 3.49e-08
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES     152.5   0.006558     25.54     910.3
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.819  (se = 0.036 )
# Likelihood ratio test= 42.69  on 1 df,   p=6e-11
# Wald test            = 30.41  on 1 df,   p=3e-08
# Score (logrank) test = 43.84  on 1 df,   p=4e-11
mylogit <- glm(Recurrence ~ TIMES, data = data_MVI_M0_vs_M1M2_yes,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_MVI_M0_vs_M1M2_yes)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.4011  -0.5351   0.2787   0.5067   2.2742  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -3.0347     0.6895  -4.401 1.08e-05 ***
#   TIMES         6.6056     1.3437   4.916 8.83e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 103.959  on 74  degrees of freedom
# Residual deviance:  62.632  on 73  degrees of freedom
# AIC: 66.632
# 
# Number of Fisher Scoring iterations: 5


##3. macrovascular invasion
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


##stratified test: no group
all_data_no <- all_data[which(all_data$macrovascular_invasion=="no"),] #124
data_macrovascular_invasion_no  <- all_data_no %>%
  dplyr::select(DFS, Recurrence, macrovascular_invasion, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_macrovascular_invasion_no, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_macrovascular_invasion_no, 
#         id = patient_id)
# 
# n= 91, number of events= 40 
# (33 observations deleted due to missingness)
# 
# 1:2        coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES   6.898   990.140    1.024 6.739  1.6e-11
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES     990.1    0.00101     133.2      7362
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.844  (se = 0.032 )
# Likelihood ratio test= 65.18  on 1 df,   p=7e-16
# Wald test            = 45.41  on 1 df,   p=2e-11
# Score (logrank) test = 62.34  on 1 df,   p=3e-15

mylogit <- glm(Recurrence ~ TIMES, data = data_macrovascular_invasion_no,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_macrovascular_invasion_no)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.3003  -0.6834  -0.3616   0.6732   2.3930  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -3.1451     0.5706  -5.512 3.54e-08 ***
#   TIMES         6.4449     1.1244   5.732 9.93e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 170.32  on 123  degrees of freedom
# Residual deviance: 117.40  on 122  degrees of freedom
# AIC: 121.4
# 
# Number of Fisher Scoring iterations: 4

##stratified test: yes group
all_data_yes <- all_data[which(all_data$macrovascular_invasion=="yes"),] #19
data_macrovascular_invasion_yes  <- all_data_yes %>%
  dplyr::select(DFS, Recurrence, macrovascular_invasion, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_macrovascular_invasion_yes, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_macrovascular_invasion_yes, 
#         id = patient_id)
# 
# n= 17, number of events= 13 
# (2 observations deleted due to missingness)
# 
# 1:2       coef exp(coef) se(coef)   z Pr(>|z|)
# TIMES  4.132    62.315    1.968 2.1   0.0358
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES     62.31    0.01605     1.316      2950
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.767  (se = 0.06 )
# Likelihood ratio test= 9.16  on 1 df,   p=0.002
# Wald test            = 4.41  on 1 df,   p=0.04
# Score (logrank) test = 6.9  on 1 df,   p=0.009

mylogit <- glm(Recurrence ~ TIMES, data = data_macrovascular_invasion_yes,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_macrovascular_invasion_yes)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.8270  -0.1260   0.3402   0.3726   2.0267  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   -2.429      1.517  -1.601   0.1095  
# TIMES          6.420      2.654   2.419   0.0156 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 21.901  on 18  degrees of freedom
# Residual deviance: 12.485  on 17  degrees of freedom
# AIC: 16.485
# 
# Number of Fisher Scoring iterations: 5


##4. abnormality of hepatitis B virus surface antigen quantification ####
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


##stratified test: "-" negative group
all_data_negative <- all_data[which(all_data$HPB_surface_antigen=="-"),] #26
data_HPB_surface_antigen_negative  <- all_data_negative %>%
  dplyr::select(DFS, Recurrence, HPB_surface_antigen, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_surface_antigen_negative, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_surface_antigen_negative, 
#         id = patient_id)
# 
# n= 20, number of events= 11 
# (6 observations deleted due to missingness)
# 
# 1:2        coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES   4.646   104.188    1.532 3.032  0.00243
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES     104.2   0.009598     5.171      2099
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.789  (se = 0.054 )
# Likelihood ratio test= 11.29  on 1 df,   p=8e-04
# Wald test            = 9.19  on 1 df,   p=0.002
# Score (logrank) test = 11.2  on 1 df,   p=8e-04

mylogit <- glm(Recurrence ~ TIMES, data = data_HPB_surface_antigen_negative,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_HPB_surface_antigen_negative)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.6782  -0.4119   0.1017   0.3653   1.8599  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   -4.911      2.146  -2.288   0.0221 *
#   TIMES         11.151      4.592   2.428   0.0152 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 35.890  on 25  degrees of freedom
# Residual deviance: 17.853  on 24  degrees of freedom
# AIC: 21.853
# 
# Number of Fisher Scoring iterations: 6


##stratified test: "+" positive group
all_data_positive <- all_data[which(all_data$HPB_surface_antigen=="+"),] #109
data_HPB_surface_antigen_positive  <- all_data_positive %>%
  dplyr::select(DFS, Recurrence, HPB_surface_antigen, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_surface_antigen_positive, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_surface_antigen_positive, 
#         id = patient_id)
# 
# n= 84, number of events= 41 
# (25 observations deleted due to missingness)
# 
# 1:2         coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES    7.543  1887.937    1.171 6.441 1.19e-10
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES      1888  0.0005297     190.1     18746
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.855  (se = 0.034 )
# Likelihood ratio test= 71.63  on 1 df,   p=<2e-16
# Wald test            = 41.48  on 1 df,   p=1e-10
# Score (logrank) test = 63.24  on 1 df,   p=2e-15

mylogit <- glm(Recurrence ~ TIMES, data = data_HPB_surface_antigen_positive,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_HPB_surface_antigen_positive)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.2086  -0.6845  -0.3767   0.6308   2.3527  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -3.0228     0.5877  -5.143 2.70e-07 ***
#   TIMES         6.0540     1.0860   5.575 2.48e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 150.66  on 108  degrees of freedom
# Residual deviance: 102.30  on 107  degrees of freedom
# AIC: 106.3
# 
# Number of Fisher Scoring iterations: 4


####5. abnormality of hepatitis B virus surface antibody quantification ####
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

##stratified test: "-" negative group
all_data_negative <- all_data[which(all_data$HPB_surface_antibody=="-"),] #118
data_HPB_surface_antibody_negative  <- all_data_negative %>%
  dplyr::select(DFS, Recurrence, HPB_surface_antibody, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_surface_antibody_negative, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_surface_antibody_negative, 
#         id = patient_id)
# 
# n= 90, number of events= 45 
# (28 observations deleted due to missingness)
# 
# 1:2        coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES   6.888   980.823    1.024 6.724 1.76e-11
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES     980.8    0.00102     131.7      7304
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.845  (se = 0.031 )
# Likelihood ratio test= 72.61  on 1 df,   p=<2e-16
# Wald test            = 45.22  on 1 df,   p=2e-11
# Score (logrank) test = 65.86  on 1 df,   p=5e-16

mylogit <- glm(Recurrence ~ TIMES, data = data_HPB_surface_antibody_negative,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_HPB_surface_antibody_negative)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.2886  -0.6468  -0.3556   0.5922   2.3842  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -3.1196     0.5797  -5.382 7.38e-08 ***
#   TIMES         6.3835     1.0917   5.848 4.99e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 163.28  on 117  degrees of freedom
# Residual deviance: 107.58  on 116  degrees of freedom
# AIC: 111.58
# 
# Number of Fisher Scoring iterations: 4


##stratified test: "+" positive group
all_data_positive <- all_data[which(all_data$HPB_surface_antibody=="+"),] #18
data_HPB_surface_antibody_positive  <- all_data_positive %>%
  dplyr::select(DFS, Recurrence, HPB_surface_antibody, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_surface_antibody_positive, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_surface_antibody_positive, 
#         id = patient_id)
# 
# n= 15, number of events= 7 
# (3 observations deleted due to missingness)
# 
# 1:2        coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES   6.586   725.208    2.692 2.447   0.0144
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES     725.2   0.001379     3.706    141921
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.841  (se = 0.06 )
# Likelihood ratio test= 9.72  on 1 df,   p=0.002
# Wald test            = 5.99  on 1 df,   p=0.01
# Score (logrank) test = 8.93  on 1 df,   p=0.003

mylogit <- glm(Recurrence ~ TIMES, data = data_HPB_surface_antibody_positive,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_HPB_surface_antibody_positive)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.76741  -0.52804  -0.03958   0.46884   1.80648  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   -3.926      1.841  -2.133   0.0329 *
#   TIMES          8.298      3.696   2.245   0.0248 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 24.953  on 17  degrees of freedom
# Residual deviance: 14.885  on 16  degrees of freedom
# AIC: 18.885
# 
# Number of Fisher Scoring iterations: 5


####6. abnormality of hepatitis B virus e antigen quantification #### 
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

##stratified test: "-" negative group
all_data_negative <- all_data[which(all_data$HPB_e_antigen=="-"),] #114
data_HPB_e_antigen_negative  <- all_data_negative %>%
  dplyr::select(DFS, Recurrence, HPB_e_antigen, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_e_antigen_negative, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_e_antigen_negative, 
#         id = patient_id)
# 
# n= 88, number of events= 44 
# (26 observations deleted due to missingness)
# 
# 1:2         coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES   6.3347  563.8219   0.9667 6.553 5.65e-11
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES     563.8   0.001774     84.77      3750
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.833  (se = 0.032 )
# Likelihood ratio test= 64.11  on 1 df,   p=1e-15
# Wald test            = 42.94  on 1 df,   p=6e-11
# Score (logrank) test = 58.85  on 1 df,   p=2e-14

mylogit <- glm(Recurrence ~ TIMES, data = data_HPB_e_antigen_negative,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_HPB_e_antigen_negative)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.3233  -0.6085  -0.3479   0.6270   2.3898  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -3.1401     0.5988  -5.244 1.57e-07 ***
#   TIMES         6.5035     1.1334   5.738 9.58e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 157.90  on 113  degrees of freedom
# Residual deviance: 103.51  on 112  degrees of freedom
# AIC: 107.51
# 
# Number of Fisher Scoring iterations: 4


##stratified test: "+" positive group
all_data_positive <- all_data[which(all_data$HPB_e_antigen=="+"),] #21
data_HPB_e_antigen_positive  <- all_data_positive %>%
  dplyr::select(DFS, Recurrence, HPB_e_antigen, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_e_antigen_positive, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_e_antigen_positive, 
#         id = patient_id)
# 
# n= 16, number of events= 8 
# (5 observations deleted due to missingness)
# 
# 1:2         coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES    7.772  2372.899    2.830 2.746  0.00603
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES      2373  0.0004214      9.25    608693
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.876  (se = 0.05 )
# Likelihood ratio test= 15.28  on 1 df,   p=9e-05
# Wald test            = 7.54  on 1 df,   p=0.006
# Score (logrank) test = 13.64  on 1 df,   p=2e-04

mylogit <- glm(Recurrence ~ TIMES, data = data_HPB_e_antigen_positive,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_HPB_e_antigen_positive)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.2985  -0.7321  -0.3364   0.5072   2.3450  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   -3.319      1.422  -2.335   0.0196 *
#   TIMES          6.648      2.709   2.454   0.0141 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 29.065  on 20  degrees of freedom
# Residual deviance: 18.900  on 19  degrees of freedom
# AIC: 22.9
# 
# Number of Fisher Scoring iterations: 4


####7. abnormality of hepatitis B virus e antibody quantification  (low, normal as negative, high as positive) ####
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

##stratified test: "-" negative group
all_data_negative <- all_data[which(all_data$HPB_e_antibody=="-"),] #119
data_HPB_e_antibody_negative  <- all_data_negative %>%
  dplyr::select(DFS, Recurrence, HPB_e_antibody, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_e_antibody_negative, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_e_antibody_negative, 
#         id = patient_id)
# 
# n= 91, number of events= 44 
# (28 observations deleted due to missingness)
# 
# 1:2        coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES   6.229   507.373    0.925 6.734 1.65e-11
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES     507.4   0.001971     82.79      3109
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.837  (se = 0.029 )
# Likelihood ratio test= 67.25  on 1 df,   p=2e-16
# Wald test            = 45.35  on 1 df,   p=2e-11
# Score (logrank) test = 63.3  on 1 df,   p=2e-15

mylogit <- glm(Recurrence ~ TIMES, data = data_HPB_e_antibody_negative,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_HPB_e_antibody_negative)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.1875  -0.7137  -0.3650   0.6347   2.3652  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -3.0530     0.5656  -5.397 6.76e-08 ***
#   TIMES         6.0305     1.0422   5.786 7.20e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 164.29  on 118  degrees of freedom
# Residual deviance: 111.86  on 117  degrees of freedom
# AIC: 115.86
# 
# Number of Fisher Scoring iterations: 4


##stratified test: "+" positive group
all_data_positive <- all_data[which(all_data$HPB_e_antibody=="+"),] #17
data_HPB_e_antibody_positive  <- all_data_positive %>%
  dplyr::select(DFS, Recurrence, HPB_e_antibody, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_e_antibody_positive, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_e_antibody_positive, 
#         id = patient_id)
# 
# n= 14, number of events= 8 
# (3 observations deleted due to missingness)
# 
# 1:2          coef exp(coef)  se(coef)     z Pr(>|z|)
# TIMES 1.176e+01 1.283e+05 4.997e+00 2.354   0.0186
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES    128285  7.795e-06      7.16 2.298e+09
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.897  (se = 0.054 )
# Likelihood ratio test= 14.6  on 1 df,   p=1e-04
# Wald test            = 5.54  on 1 df,   p=0.02
# Score (logrank) test = 10.48  on 1 df,   p=0.001

mylogit <- glm(Recurrence ~ TIMES, data = data_HPB_e_antibody_positive,family=binomial)
# Warning messages:
#   1: glm.fit: algorithm did not converge 
# 2: glm.fit: fitted probabilities numerically 0 or 1 occurred 
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_HPB_e_antibody_positive)
# 
# Deviance Residuals: 
#   Min          1Q      Median          3Q         Max  
# -1.136e-04  -2.100e-08   2.100e-08   2.100e-08   1.158e-04  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)    -1711     480468  -0.004    0.997
# TIMES           4200    1179150   0.004    0.997
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 2.3035e+01  on 16  degrees of freedom
# Residual deviance: 2.6307e-08  on 15  degrees of freedom
# AIC: 4
# 
# Number of Fisher Scoring iterations: 25


####8. abnormality of hepatitis B virus core antibody quantitation  ####
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

##stratified test: "-" negative group
all_data_negative <- all_data[which(all_data$HPB_core_antibody=="-"),] #84
data_HPB_core_antibody_negative  <- all_data_negative %>%
  dplyr::select(DFS, Recurrence, HPB_core_antibody, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_core_antibody_negative, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_core_antibody_negative, 
#         id = patient_id)
# 
# n= 65, number of events= 33 
# (19 observations deleted due to missingness)
# 
# 1:2         coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES    7.189  1324.707    1.287 5.585 2.34e-08
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES      1325  0.0007549     106.3     16512
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.869  (se = 0.029 )
# Likelihood ratio test= 54.97  on 1 df,   p=1e-13
# Wald test            = 31.19  on 1 df,   p=2e-08
# Score (logrank) test = 46.97  on 1 df,   p=7e-12

mylogit <- glm(Recurrence ~ TIMES, data = data_HPB_core_antibody_negative,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_HPB_core_antibody_negative)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.2702  -0.6310   0.3758   0.6020   2.2977  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.8864     0.6663  -4.332 1.48e-05 ***
#   TIMES         6.0697     1.2171   4.987 6.13e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 116.258  on 83  degrees of freedom
# Residual deviance:  77.167  on 82  degrees of freedom
# AIC: 81.167
# 
# Number of Fisher Scoring iterations: 4


##stratified test: "+" positive group
all_data_positive <- all_data[which(all_data$HPB_core_antibody=="+"),] #37
data_HPB_core_antibody_positive  <- all_data_positive %>%
  dplyr::select(DFS, Recurrence, HPB_core_antibody, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_core_antibody_positive, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_HPB_core_antibody_positive, 
#         id = patient_id)
# 
# n= 25, number of events= 14 
# (12 observations deleted due to missingness)
# 
# 1:2        coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES   6.358   577.163    1.740 3.654 0.000258
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES     577.2   0.001733     19.07     17470
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.819  (se = 0.056 )
# Likelihood ratio test= 20.41  on 1 df,   p=6e-06
# Wald test            = 13.35  on 1 df,   p=3e-04
# Score (logrank) test = 19.57  on 1 df,   p=1e-05

mylogit <- glm(Recurrence ~ TIMES, data = data_HPB_core_antibody_positive,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_HPB_core_antibody_positive)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.5156  -0.4884  -0.2460   0.7013   1.7140  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)   -4.094      1.314  -3.116  0.00183 **
#   TIMES          8.539      2.667   3.202  0.00136 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 50.615  on 36  degrees of freedom
# Residual deviance: 29.030  on 35  degrees of freedom
# AIC: 33.03
# 
# Number of Fisher Scoring iterations: 5


####9. abnormality of hepatitis C antibody-IgG ####
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

##stratified test: "-" negative group
all_data_negative <- all_data[which(all_data$HPC_antibody=="-"),] #127
data_HPC_antibody_negative  <- all_data_negative %>%
  dplyr::select(DFS, Recurrence, HPC_antibody, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_HPC_antibody_negative, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_HPC_antibody_negative, 
#         id = patient_id)
# 
# n= 102, number of events= 51 
# (25 observations deleted due to missingness)
# 
# 1:2          coef exp(coef)  se(coef)     z Pr(>|z|)
# TIMES    7.1573 1283.4687    0.9694 7.383 1.55e-13
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES      1283  0.0007791       192      8582
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.854  (se = 0.025 )
# Likelihood ratio test= 88.42  on 1 df,   p=<2e-16
# Wald test            = 54.51  on 1 df,   p=2e-13
# Score (logrank) test = 80.05  on 1 df,   p=<2e-16

mylogit <- glm(Recurrence ~ TIMES, data = data_HPC_antibody_negative,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_HPC_antibody_negative)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.4333  -0.5433  -0.2858   0.5529   2.4825  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -3.6218     0.6396  -5.663 1.49e-08 ***
#   TIMES         7.3602     1.1994   6.136 8.44e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 175.99  on 126  degrees of freedom
# Residual deviance: 107.23  on 125  degrees of freedom
# AIC: 111.23
# 
# Number of Fisher Scoring iterations: 5


##stratified test: "+" positive group
all_data_positive <- all_data[which(all_data$HPC_antibody=="+"),] #4
data_HPC_antibody_positive  <- all_data_positive %>%
  dplyr::select(DFS, Recurrence, HPC_antibody, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_HPC_antibody_positive, id = patient_id)
# Error in dimnames(events) <- list(count = tab1.levels, state = c(ystate,  : 
#                                                                    'dimnames' applied to non-array
summary(res.cox)

mylogit <- glm(Recurrence ~ TIMES, data = data_HPC_antibody_positive,family=binomial)
# Warning message:
#   glm.fit: fitted probabilities numerically 0 or 1 occurred 
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_HPC_antibody_positive)
# 
# Deviance Residuals: 
#   1           2           3           4  
# -2.110e-08  -1.549e-05   1.289e-05  -2.110e-08  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   -142.2   263202.0  -0.001        1
# TIMES          220.4   414217.2   0.001        1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 4.4987e+00  on 3  degrees of freedom
# Residual deviance: 4.0608e-10  on 2  degrees of freedom
# AIC: 4
# 
# Number of Fisher Scoring iterations: 24


####10. liver cirrhosis####
colnames(all_data)[which(colnames(all_data)=="Liver cirrhosis")] <- "liver_cirrhosis"
all_data$liver_cirrhosis <- factor(all_data$liver_cirrhosis,levels = c("-","+"))
table(all_data$liver_cirrhosis)
# -   + 
#   28 111 
# > 111/(111+28)
# [1] 0.7985612

###no/mild fibrosis (F0-F1), moderate fibrosis (F2-F3), and cirrhosis (F4) based on the METAVIR scoring system
table(all_data$`METAVIR score`)
# 0   1   2   3   4  
# 16  21  35  37 116  
all_data$fibr_3class <- NA
all_data$fibr_3class[grep('^["0"|"1"]$',all_data$`METAVIR score`)]<-'I'
all_data$fibr_3class[grep('^["2"|"3"]$',all_data$`METAVIR score`)]<-'II'
all_data$fibr_3class[grep('^["4"]$',all_data$`METAVIR score`)]<-'III'
all_data$fibr_3class <- factor(all_data$fibr_3class,levels = c("I","II","III"))
table(all_data$fibr_3class)
# I  II III 
# 37  72 116 
#116/(37+  72+ 116)=0.5155556
#72/(37+  72+ 116)=0.32

##stratified test: "I"  group
all_data_classI <- all_data[which(all_data$fibr_3class=="I"),] #37
data_fibr_3class_classI  <- all_data_classI %>%
  dplyr::select(DFS, Recurrence, fibr_3class, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_fibr_3class_classI, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_fibr_3class_classI, 
#         id = patient_id)
# 
# n= 15, number of events= 4 
# (22 observations deleted due to missingness)
# 
# 1:2        coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES   6.494   660.938    2.984 2.176   0.0296
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES     660.9   0.001513     1.905    229312
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.82  (se = 0.105 )
# Likelihood ratio test= 9.68  on 1 df,   p=0.002
# Wald test            = 4.73  on 1 df,   p=0.03
# Score (logrank) test = 13.42  on 1 df,   p=2e-04

mylogit <- glm(Recurrence ~ TIMES, data = data_fibr_3class_classI,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_fibr_3class_classI)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.0576  -0.5961  -0.3714  -0.3294   2.3555  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)   -3.132      1.104  -2.838  0.00454 **
#   TIMES          5.296      2.393   2.213  0.02688 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 24.085  on 22  degrees of freedom
# Residual deviance: 17.462  on 21  degrees of freedom
# (14 observations deleted due to missingness)
# AIC: 21.462
# 
# Number of Fisher Scoring iterations: 5


##stratified test: "II"  group
all_data_classI <- all_data[which(all_data$fibr_3class=="II"),] #72
data_fibr_3class_classI  <- all_data_classI %>%
  dplyr::select(DFS, Recurrence, fibr_3class, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_fibr_3class_classI, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_fibr_3class_classI, 
#         id = patient_id)
# 
# n= 22, number of events= 3 
# (50 observations deleted due to missingness)
# 
# 1:2          coef exp(coef)  se(coef)     z Pr(>|z|)
# TIMES 1.252e+01 2.732e+05 6.207e+00 2.017   0.0437
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES    273192   3.66e-06     1.423 5.246e+10
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.981  (se = 0.017 )
# Likelihood ratio test= 12.45  on 1 df,   p=4e-04
# Wald test            = 4.07  on 1 df,   p=0.04
# Score (logrank) test = 11.76  on 1 df,   p=6e-04

mylogit <- glm(Recurrence ~ TIMES, data = data_fibr_3class_classI,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_fibr_3class_classI)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.18051  -0.35239  -0.05560  -0.02089   2.03117  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   -9.258      4.671  -1.982   0.0475 *
#   TIMES         14.635      8.124   1.801   0.0716 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 21.2537  on 22  degrees of freedom
# Residual deviance:  9.7387  on 21  degrees of freedom
# (49 observations deleted due to missingness)
# AIC: 13.739
# 
# Number of Fisher Scoring iterations: 7

##stratified test: "III"  group
all_data_classI <- all_data[which(all_data$fibr_3class=="III"),] #116
data_fibr_3class_classI  <- all_data_classI %>%
  dplyr::select(DFS, Recurrence, fibr_3class, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_fibr_3class_classI, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_fibr_3class_classI, 
#         id = patient_id)
# 
# n= 85, number of events= 46 
# (31 observations deleted due to missingness)
# 
# 1:2         coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES   6.3952  598.9440   0.9687 6.602 4.07e-11
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES     598.9    0.00167      89.7      3999
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.844  (se = 0.028 )
# Likelihood ratio test= 68.13  on 1 df,   p=<2e-16
# Wald test            = 43.58  on 1 df,   p=4e-11
# Score (logrank) test = 62.11  on 1 df,   p=3e-15

mylogit <- glm(Recurrence ~ TIMES, data = data_fibr_3class_classI,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_fibr_3class_classI)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.5368  -0.5608   0.3154   0.5905   2.3442  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   -3.053      0.610  -5.004 5.61e-07 ***
#   TIMES          7.022      1.225   5.731 1.00e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 151.582  on 109  degrees of freedom
# Residual deviance:  94.571  on 108  degrees of freedom
# (6 observations deleted due to missingness)
# AIC: 98.571
# 
# Number of Fisher Scoring iterations: 5


####11. viral HCC####
all_data$virus <- NA
temp <- which((all_data$HPB_surface_antigen =="+") | (all_data$HPC_antibody=="+"))
temp2 <- which((all_data$HPB_surface_antigen =="-") & (all_data$HPC_antibody=="-"))
all_data$virus[temp] <- "+"
all_data$virus[temp2] <- "-"
all_data$virus <- factor(all_data$virus,levels = c("-","+"))
table(all_data$virus)
# -   + 
#   23 113 
#113/(23+113) = 0.8308824

data_virus  <- all_data %>%
  dplyr::select(DFS, Recurrence, virus, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ virus, data = data_virus, id = patient_id)
res_out <- summary(res.cox)
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ virus, data = data_virus, 
#         id = patient_id)
# 
# n= 105, number of events= 52 
# (149 observations deleted due to missingness)
# 
# 1:2         coef exp(coef) se(coef)     z Pr(>|z|)
# virus+ -0.6159    0.5401   0.3479 -1.77   0.0767
# 
# 1:2      exp(coef) exp(-coef) lower .95 upper .95
# virus+    0.5401      1.851    0.2731     1.068
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.551  (se = 0.031 )
# Likelihood ratio test= 2.8  on 1 df,   p=0.09
# Wald test            = 3.13  on 1 df,   p=0.08
# Score (logrank) test = 3.23  on 1 df,   p=0.07

mylogit <- glm(Recurrence ~ virus, data = data_virus,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ virus, family = binomial, data = data_virus)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -1.291  -1.110  -1.110   1.246   1.246  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.2624     0.4206   0.624    0.533
# virus+       -0.4220     0.4610  -0.915    0.360
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 188.27  on 135  degrees of freedom
# Residual deviance: 187.43  on 134  degrees of freedom
# (118 observations deleted due to missingness)
# AIC: 191.43
# 
# Number of Fisher Scoring iterations: 3

##stratified test: "-" negative group
all_data_negative <- all_data[which(all_data$virus=="-"),] #23
data_virus_negative  <- all_data_negative %>%
  dplyr::select(DFS, Recurrence, virus, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_virus_negative, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_virus_negative, 
#         id = patient_id)
# 
# n= 19, number of events= 11 
# (4 observations deleted due to missingness)
# 
# 1:2       coef exp(coef) se(coef)    z Pr(>|z|)
# TIMES  4.586    98.059    1.498 3.06  0.00221
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES     98.06     0.0102       5.2      1849
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.809  (se = 0.052 )
# Likelihood ratio test= 11.66  on 1 df,   p=6e-04
# Wald test            = 9.36  on 1 df,   p=0.002
# Score (logrank) test = 11.75  on 1 df,   p=6e-04

mylogit <- glm(Recurrence ~ TIMES, data = data_virus_negative,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_virus_negative)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.55954  -0.35352   0.06685   0.45362   1.84352  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   -5.468      2.669  -2.049   0.0405 *
#   TIMES         13.113      6.043   2.170   0.0300 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 31.492  on 22  degrees of freedom
# Residual deviance: 14.617  on 21  degrees of freedom
# AIC: 18.617
# 
# Number of Fisher Scoring iterations: 6


##stratified test: "+" positive group
all_data_positive <- all_data[which(all_data$virus=="+"),] #113
data_virus_positive  <- all_data_positive %>%
  dplyr::select(DFS, Recurrence, virus, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_virus_positive, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_virus_positive, 
#         id = patient_id)
# 
# n= 86, number of events= 41 
# (27 observations deleted due to missingness)
# 
# 1:2         coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES    7.597  1992.368    1.186 6.407 1.48e-10
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES      1992  0.0005019       195     20355
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.856  (se = 0.034 )
# Likelihood ratio test= 71.14  on 1 df,   p=<2e-16
# Wald test            = 41.05  on 1 df,   p=1e-10
# Score (logrank) test = 62.42  on 1 df,   p=3e-15

mylogit <- glm(Recurrence ~ TIMES, data = data_virus_positive,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_virus_positive)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.2107  -0.6670  -0.3624   0.6313   2.3850  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -3.1095     0.5884  -5.285 1.26e-07 ***
#   TIMES         6.1575     1.0810   5.696 1.23e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 155.93  on 112  degrees of freedom
# Residual deviance: 104.80  on 111  degrees of freedom
# AIC: 108.8
# 
# Number of Fisher Scoring iterations: 4


##multi-variate regression
data_virus_inf  <- all_data %>%
  dplyr::select(DFS, Recurrence, HPB_surface_antigen, HPC_antibody, fibr_3class, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES + HPB_surface_antigen + HPC_antibody + fibr_3class, data = data_virus_inf, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES + HPB_surface_antigen + 
#           HPC_antibody + fibr_3class, data = data_virus_inf, id = patient_id)
# 
# n= 103, number of events= 51 
# (151 observations deleted due to missingness)
# 
# 1:2                          coef  exp(coef)   se(coef)      z Pr(>|z|)
# TIMES                 6.921e+00  1.014e+03  9.373e-01  7.384 1.53e-13
# HPB_surface_antigen+ -6.181e-01  5.390e-01  3.783e-01 -1.634    0.102
# HPC_antibody+        -1.774e+01  1.972e-08  4.104e+03 -0.004    0.997
# fibr_3classII        -5.793e-01  5.603e-01  7.940e-01 -0.730    0.466
# fibr_3classIII        1.712e-01  1.187e+00  5.455e-01  0.314    0.754
# 
# 1:2                    exp(coef) exp(-coef) lower .95 upper .95
# TIMES                1.014e+03  9.865e-04  161.4605  6364.450
# HPB_surface_antigen+ 5.390e-01  1.855e+00    0.2568     1.131
# HPC_antibody+        1.972e-08  5.072e+07    0.0000       Inf
# fibr_3classII        5.603e-01  1.785e+00    0.1182     2.656
# fibr_3classIII       1.187e+00  8.426e-01    0.4074     3.457
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.866  (se = 0.023 )
# Likelihood ratio test= 95.1  on 5 df,   p=<2e-16
# Wald test            = 62.31  on 5 df,   p=4e-12
# Score (logrank) test = 91.23  on 5 df,   p=<2e-16

mylogit <- glm(Recurrence ~ TIMES + HPB_surface_antigen + HPC_antibody + fibr_3class, data = data_virus_inf,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES + HPB_surface_antigen + HPC_antibody + 
#         fibr_3class, family = binomial, data = data_virus_inf)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.6346  -0.4889  -0.1572   0.4467   2.7568  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           -3.7222     1.1909  -3.126  0.00177 ** 
#   TIMES                  8.5226     1.4288   5.965 2.45e-09 ***
#   HPB_surface_antigen+  -0.7351     0.6661  -1.104  0.26972    
# HPC_antibody+         -2.1312     1.6818  -1.267  0.20508    
# fibr_3classII         -1.2285     1.2231  -1.004  0.31517    
# fibr_3classIII         0.3361     1.0062   0.334  0.73837    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 177.321  on 127  degrees of freedom
# Residual deviance:  94.586  on 122  degrees of freedom
# (126 observations deleted due to missingness)
# AIC: 106.59
# 
# Number of Fisher Scoring iterations: 5


###11. tumor count ###
colnames(all_data)[which(colnames(all_data)=="Tumor count")] <- "tumor_count"
summary(all_data$tumor_count)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.000   1.000   1.000   1.133   1.000   5.000     126   

all_data$tumor_count_binary <- NA
all_data$tumor_count_binary[which(all_data$tumor_count<=1)] <-'≤1'
all_data$tumor_count_binary[which(all_data$tumor_count>1)]<-'>1'
all_data$tumor_count_binary <- factor(all_data$tumor_count_binary,levels = c("≤1",">1"))
table(all_data$tumor_count_binary)
# ≤1  >1 
# 115  13 

data_factor  <- all_data %>%
  dplyr::select(DFS, Recurrence, tumor_count_binary, patient_id)
res.cox <- coxph(Surv(DFS, Recurrence) ~ tumor_count_binary, data = data_factor, id = patient_id)
res_out <- summary(res.cox)
res_out
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ tumor_count_binary, data = data_factor, 
#         id = patient_id)
# 
# n= 100, number of events= 49 
# (154 observations deleted due to missingness)
# 
# 1:2                       coef exp(coef) se(coef)     z Pr(>|z|)
# tumor_count_binary>1 -0.4240    0.6544   0.5973 -0.71    0.478
# 
# 1:2                    exp(coef) exp(-coef) lower .95 upper .95
# tumor_count_binary>1    0.6544      1.528     0.203      2.11
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.517  (se = 0.019 )
# Likelihood ratio test= 0.57  on 1 df,   p=0.4
# Wald test            = 0.5  on 1 df,   p=0.5
# Score (logrank) test = 0.51  on 1 df,   p=0.5

##stratified test: "-" negative group
all_data_negative <- all_data[which(all_data$tumor_count_binary=="≤1"),] #115
data_factor_negative  <- all_data_negative %>%
  dplyr::select(DFS, Recurrence, tumor_count_binary, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_factor_negative, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_factor_negative, 
#         id = patient_id)
# 
# n= 92, number of events= 45 
# (23 observations deleted due to missingness)
# 
# 1:2         coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES   6.5249  681.9035   0.9217 7.079 1.45e-12
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES     681.9   0.001466       112      4153
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.839  (se = 0.028 )
# Likelihood ratio test= 76.64  on 1 df,   p=<2e-16
# Wald test            = 50.11  on 1 df,   p=1e-12
# Score (logrank) test = 73.55  on 1 df,   p=<2e-16

mylogit <- glm(Recurrence ~ TIMES, data = data_factor_negative,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_factor_negative)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.4866  -0.5281  -0.2893   0.4909   2.4756  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -3.6155     0.6556  -5.515 3.49e-08 ***
#   TIMES         7.5083     1.2712   5.906 3.50e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 159.206  on 114  degrees of freedom
# Residual deviance:  93.124  on 113  degrees of freedom
# AIC: 97.124
# 
# Number of Fisher Scoring iterations: 5

##stratified test: "+" positive group
all_data_positive <- all_data[which(all_data$tumor_count_binary==">1"),] #13
data_factor_positive  <- all_data_positive %>%
  dplyr::select(DFS, Recurrence, tumor_count_binary, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_factor_positive, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_factor_positive, 
#         id = patient_id)
# 
# n= 8, number of events= 4 
# (5 observations deleted due to missingness)
# 
# 1:2          coef exp(coef)  se(coef)     z Pr(>|z|)
# TIMES 6.644e+01 7.188e+28 8.501e+01 0.782    0.434
# 
# 1:2     exp(coef) exp(-coef) lower .95  upper .95
# TIMES 7.188e+28  1.391e-29  3.12e-44 1.656e+101
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.944  (se = 0.057 )
# Likelihood ratio test= 8.72  on 1 df,   p=0.003
# Wald test            = 0.61  on 1 df,   p=0.4
# Score (logrank) test = 4.24  on 1 df,   p=0.04

mylogit <- glm(Recurrence ~ TIMES, data = data_factor_positive,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_factor_positive)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.8518  -0.3487  -0.1383   0.6947   1.4506  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   -5.929      3.604  -1.645   0.1000 .
# TIMES          9.516      5.411   1.759   0.0786 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 17.3232  on 12  degrees of freedom
# Residual deviance:  9.5934  on 11  degrees of freedom
# AIC: 13.593
# 
# Number of Fisher Scoring iterations: 6


###12. satellite lesion ###                                          
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

##stratified test: "-" negative group
all_data_negative <- all_data[which(all_data$satellite_lesion=="-"),] #93
data_factor_negative  <- all_data_negative %>%
  dplyr::select(DFS, Recurrence, satellite_lesion, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_factor_negative, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_factor_negative, 
#         id = patient_id)
# 
# n= 74, number of events= 30 
# (19 observations deleted due to missingness)
# 
# 1:2         coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES    7.501  1809.379    1.278 5.868  4.4e-09
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES      1809  0.0005527     147.8     22156
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.854  (se = 0.039 )
# Likelihood ratio test= 54.36  on 1 df,   p=2e-13
# Wald test            = 34.44  on 1 df,   p=4e-09
# Score (logrank) test = 52.28  on 1 df,   p=5e-13

mylogit <- glm(Recurrence ~ TIMES, data = data_factor_negative,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_factor_negative)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.2873  -0.5269  -0.2907   0.5428   2.5663  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -3.6419     0.7202  -5.057 4.26e-07 ***
#   TIMES         7.3166     1.3799   5.302 1.14e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 127.102  on 92  degrees of freedom
# Residual deviance:  76.701  on 91  degrees of freedom
# AIC: 80.701
# 
# Number of Fisher Scoring iterations: 5

##stratified test: "+" positive group
all_data_positive <- all_data[which(all_data$satellite_lesion=="+"),] #33
data_factor_positive  <- all_data_positive %>%
  dplyr::select(DFS, Recurrence, satellite_lesion, patient_id, TIMES)
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES, data = data_factor_positive, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES, data = data_factor_positive, 
#         id = patient_id)
# 
# n= 25, number of events= 17 
# (8 observations deleted due to missingness)
# 
# 1:2       coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES  4.568    96.357    1.323 3.453 0.000554
# 
# 1:2     exp(coef) exp(-coef) lower .95 upper .95
# TIMES     96.36    0.01038      7.21      1288
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.817  (se = 0.046 )
# Likelihood ratio test= 19.45  on 1 df,   p=1e-05
# Wald test            = 11.93  on 1 df,   p=6e-04
# Score (logrank) test = 17.45  on 1 df,   p=3e-05

mylogit <- glm(Recurrence ~ TIMES, data = data_factor_positive,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Recurrence ~ TIMES, family = binomial, data = data_factor_positive)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.3490  -0.6059   0.3873   0.4666   2.0913  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -2.5382     0.9647  -2.631  0.00851 **
#   TIMES         5.8974     1.8501   3.188  0.00143 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 44.987  on 32  degrees of freedom
# Residual deviance: 28.863  on 31  degrees of freedom
# AIC: 32.863
# 
# Number of Fisher Scoring iterations: 4

