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


data_fig3ef  <- all_data %>%
  dplyr::select(DFS, Recurrence, TIMES, BCLC_A_vs_BC, TNM_2level, `TNM staging`, `BCLC staging`, patient_id)
colnames(data_fig3ef)[7] <- "BCLC_3level"

table(data_fig3ef$`TNM staging`)
# I          IA          IB          II        IIIA        IIIB         IVA         IVB   mpT1bNxMx        mpT3 mypT3 Nx Mx          NA         pT1        pT1b  pT1b Nx Mx         pT3      pT3 pN 
# 3           8          30          73           7           4           3           8           1           1           1          90           3           7           3           3           1 
# pT3a         pT4   pT4 Nx Mx      pT4 pN       rmpT4       ypT1b 
# 1           3           1           1           1           1 

data_fig3ef$TNM_4level <- NA
data_fig3ef$TNM_4level[grep('^I[B|A]$',data_fig3ef$`TNM staging`)]<-'I'
data_fig3ef$TNM_4level[grep('^I$',data_fig3ef$`TNM staging`)]<-'I'
data_fig3ef$TNM_4level[grep('T1',data_fig3ef$`TNM staging`)]<-'I'
data_fig3ef$TNM_4level[grep('^II$',data_fig3ef$`TNM staging`)]<-'II'
data_fig3ef$TNM_4level[grep('^III[B|A]$',data_fig3ef$`TNM staging`)]<-'III'
data_fig3ef$TNM_4level[grep('T3',data_fig3ef$`TNM staging`)]<-'III'
data_fig3ef$TNM_4level[grep('^IV[B|A]$',data_fig3ef$`TNM staging`)]<-'IV'
data_fig3ef$TNM_4level[grep('T4',data_fig3ef$`TNM staging`)]<-'IV'
data_fig3ef$TNM_4level <- factor(data_fig3ef$TNM_4level,levels = c("I","II","III","IV"))
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

####### To generate survival curves: Fig. 3f and Extended Data Fig. 5c #######
#### Fig. 3f ####
data_fig3ef$TIMES_2level <- NA
data_fig3ef$TIMES_2level[which(data_fig3ef$TIMES < 0.5)] <- "N"
data_fig3ef$TIMES_2level[which(data_fig3ef$TIMES >= 0.5)] <- "Y"
data_fig3ef$TIMES_2level <- factor(data_fig3ef$TIMES_2level,levels = c("N","Y"))
table(data_fig3ef$TIMES_2level)
# N   Y 
# 143  88 

hcc_sur <- Surv(data_fig3ef$DFS, as.numeric(data_fig3ef$Recurrence))

##1. TNM and TIMES stratified
plot(survfit(formula = hcc_sur ~ TNM_2level + TIMES_2level, data=data_fig3ef), lwd=3, conf.int = F, mark.time = T, bty = "n", las = 1, xlim=c(0,40), cex = 2,
     ylab = "Disease free survival", xlab = "Time (months)", main = "", lty =c("solid", "dashed","solid", "dashed"),
     col=c(adjustcolor("#559AC6", alpha.f=1),adjustcolor("#559AC6", alpha.f=1),adjustcolor("#E73334", alpha.f=1),adjustcolor("#E73334", alpha.f=1)))
legend("bottomleft", legend=c("TNM_2level<II, TIMES_2level=N", "TNM_2level<II, TIMES_2level=Y","TNM_2level≥II, TIMES_2level=N","TNM_2level≥II, TIMES_2level=Y"), 
       col=c(adjustcolor("#559AC6", alpha.f=1),adjustcolor("#559AC6", alpha.f=1),adjustcolor("#E73334", alpha.f=1),adjustcolor("#E73334", alpha.f=1)), 
       bg=NA, lwd=4, box.lty=0, cex=1.1,lty=c("solid", "dashed","solid", "dashed"))
##save as /data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/Fig3f_survival_TNM_TIMES.pdf  : survival curves in TNM and TIMES stratified subgroups

##(1.1) in TNM <II subgroup ##
grp_idx <- which((data_fig3ef$TNM_2level=="<II") & !is.na(data_fig3ef$TIMES))  #53
data_fig3ef_grp <- data_fig3ef[grp_idx,]
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES_2level, data = data_fig3ef_grp, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES_2level, data = data_fig3ef_grp, 
#         id = patient_id)
# 
# n= 34, number of events= 14 
# (19 observations deleted due to missingness)
# 
# 1:2                coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES_2levelY  2.9369   18.8581   0.7758 3.786 0.000153
# 
# 1:2             exp(coef) exp(-coef) lower .95 upper .95
# TIMES_2levelY     18.86    0.05303     4.123     86.26
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.829  (se = 0.038 )
# Likelihood ratio test= 22.61  on 1 df,   p=2e-06
# Wald test            = 14.33  on 1 df,   p=2e-04
# Score (logrank) test = 26.1  on 1 df,   p=3e-07 (report in the figure)

##(1.2) in TNM ≥II subgroup ##
grp_idx <- which((data_fig3ef$TNM_2level=="≥II") & !is.na(data_fig3ef$TIMES))  #102
data_fig3ef_grp <- data_fig3ef[grp_idx,]
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES_2level, data = data_fig3ef_grp, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES_2level, data = data_fig3ef_grp, 
#         id = patient_id)
# 
# n= 64, number of events= 34 
# (38 observations deleted due to missingness)
# 
# 1:2                coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES_2levelY  2.4132   11.1693   0.5353 4.508 6.53e-06
# 
# 1:2             exp(coef) exp(-coef) lower .95 upper .95
# TIMES_2levelY     11.17    0.08953     3.912     31.89
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.743  (se = 0.035 )
# Likelihood ratio test= 33.18  on 1 df,   p=8e-09
# Wald test            = 20.32  on 1 df,   p=7e-06
# Score (logrank) test = 31.81  on 1 df,   p=2e-08 (report in the figure)

##2. BCLC and TIMES stratified
plot(survfit(formula = hcc_sur ~ BCLC_A_vs_BC + TIMES_2level, data=data_fig3ef), lwd=3, conf.int = F, mark.time = T, bty = "n", las = 1, xlim=c(0,40), cex = 2,
     ylab = "Disease free survival", xlab = "Time (months)", main = "", lty =c("solid", "dashed","solid", "dashed"),
     col=c(adjustcolor("#559AC6", alpha.f=1),adjustcolor("#559AC6", alpha.f=1),adjustcolor("#E73334", alpha.f=1),adjustcolor("#E73334", alpha.f=1)))
legend("bottomleft", legend=c("BCLC_A, TIMES_2level=N", "BCLC_A, TIMES_2level=Y","BCLC_BC, TIMES_2level=N","BCLC_BC, TIMES_2level=Y"), 
       col=c(adjustcolor("#559AC6", alpha.f=1),adjustcolor("#559AC6", alpha.f=1),adjustcolor("#E73334", alpha.f=1),adjustcolor("#E73334", alpha.f=1)), 
       bg=NA, lwd=4, box.lty=0, cex=1.1,lty=c("solid", "dashed","solid", "dashed"))
#save as /data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/Fig3f_survival_BCLC_TIMES.pdf  : survival curves in BCLC and TIMES stratified subgroups

##(2.1) in BCLC: A subgroup ##
grp_idx <- which((data_fig3ef$BCLC_A_vs_BC=="low") & !is.na(data_fig3ef$TIMES))  #26
data_fig3ef_grp <- data_fig3ef[grp_idx,]
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES_2level, data = data_fig3ef_grp, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES_2level, data = data_fig3ef_grp, 
#         id = patient_id)
# 
# n= 21, number of events= 6 
# (5 observations deleted due to missingness)
# 
# 1:2              coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES_2levelY 1.932     6.900    0.877 2.202   0.0276
# 
# 1:2             exp(coef) exp(-coef) lower .95 upper .95
# TIMES_2levelY       6.9     0.1449     1.237     38.49
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.75  (se = 0.092 )
# Likelihood ratio test= 5.2  on 1 df,   p=0.02
# Wald test            = 4.85  on 1 df,   p=0.03
# Score (logrank) test = 6.43  on 1 df,   p=0.01  (report in the figure)

##(2.2) in BCLC: B+C subgroup ##
grp_idx <- which((data_fig3ef$BCLC_A_vs_BC=="high") & !is.na(data_fig3ef$TIMES))  #104
data_fig3ef_grp <- data_fig3ef[grp_idx,]
res.cox <- coxph(Surv(DFS, Recurrence) ~ TIMES_2level, data = data_fig3ef_grp, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TIMES_2level, data = data_fig3ef_grp, 
#         id = patient_id)
# 
# n= 79, number of events= 42 
# (25 observations deleted due to missingness)
# 
# 1:2                coef exp(coef) se(coef)     z Pr(>|z|)
# TIMES_2levelY  2.8038   16.5079   0.5285 5.305 1.13e-07
# 
# 1:2             exp(coef) exp(-coef) lower .95 upper .95
# TIMES_2levelY     16.51    0.06058     5.859     46.51
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.765  (se = 0.03 )
# Likelihood ratio test= 51.98  on 1 df,   p=6e-13
# Wald test            = 28.15  on 1 df,   p=1e-07
# Score (logrank) test = 50.75  on 1 df,   p=1e-12
 
##combine the two figures above and generate Fig. 3f


#### Extended Data Fig. 5c ####
hcc_sur <- Surv(data_fig3ef$DFS, as.numeric(data_fig3ef$Recurrence))

##(1) TNM stratified ##
plot(survfit(formula = hcc_sur ~ TNM_4level, data=data_fig3ef), lwd=3, conf.int = F, mark.time = T, bty = "n", las = 1, xlim=c(0,40), cex = 2,
     ylab = "Disease free survival", xlab = "Time (months)", main = "TNM", col=c(adjustcolor("#1B9E77", alpha.f=1), adjustcolor("#D95F02", alpha.f=1), adjustcolor("#7570B3", alpha.f=1), adjustcolor("#E6AB02", alpha.f=1)))
legend("bottomleft", legend=c("I", "II", "III", "IV"), col=c(adjustcolor("#1B9E77", alpha.f=1), adjustcolor("#D95F02", alpha.f=1), adjustcolor("#7570B3", alpha.f=1), adjustcolor("#E6AB02", alpha.f=1)), bg=NA, lty=1, lwd=4, box.lty=0, cex=1.3)
#save as /data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/ExFig5c_survival_TNM.pdf : survival curves in TNM stratified subgroups

##cox test
grp_idx <- which(((data_fig3ef$TNM_4level=="I") | (data_fig3ef$TNM_4level=="II")) & !is.na(data_fig3ef$TIMES))  #124
data_fig3ef_grp <- data_fig3ef[grp_idx,]
data_fig3ef_grp$TNM_4level <- factor(data_fig3ef_grp$TNM_4level,levels = c("I","II"))
table(data_fig3ef_grp$TNM_4level)
# I II 
# 53 71
res.cox <- coxph(Surv(DFS, Recurrence) ~ TNM_4level, data = data_fig3ef_grp, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TNM_4level, data = data_fig3ef_grp, 
#         id = patient_id)
# 
# n= 88, number of events= 40 
# (36 observations deleted due to missingness)
# 
# 1:2              coef exp(coef) se(coef)     z Pr(>|z|)
# TNM_4levelII 0.5256    1.6914   0.3437 1.529    0.126
# 
# 1:2            exp(coef) exp(-coef) lower .95 upper .95
# TNM_4levelII     1.691     0.5912    0.8623     3.318
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.575  (se = 0.041 )
# Likelihood ratio test= 2.43  on 1 df,   p=0.1
# Wald test            = 2.34  on 1 df,   p=0.1
# Score (logrank) test = 2.38  on 1 df,   p=0.1

grp_idx <- which(((data_fig3ef$TNM_4level=="I") | (data_fig3ef$TNM_4level=="III")) & !is.na(data_fig3ef$TIMES))  #69
data_fig3ef_grp <- data_fig3ef[grp_idx,]
data_fig3ef_grp$TNM_4level <- factor(data_fig3ef_grp$TNM_4level,levels = c("I","III"))
table(data_fig3ef_grp$TNM_4level)
# I III 
# 53  16 
res.cox <- coxph(Surv(DFS, Recurrence) ~ TNM_4level, data = data_fig3ef_grp, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TNM_4level, data = data_fig3ef_grp, 
#         id = patient_id)
# 
# n= 41, number of events= 20 
# (28 observations deleted due to missingness)
# 
# 1:2               coef exp(coef) se(coef)     z Pr(>|z|)
# TNM_4levelIII 1.4757    4.3742   0.5222 2.826  0.00471
# 
# 1:2             exp(coef) exp(-coef) lower .95 upper .95
# TNM_4levelIII     4.374     0.2286     1.572     12.17
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.619  (se = 0.051 )
# Likelihood ratio test= 6.67  on 1 df,   p=0.01
# Wald test            = 7.99  on 1 df,   p=0.005
# Score (logrank) test = 9.44  on 1 df,   p=0.002

grp_idx <- which(((data_fig3ef$TNM_4level=="I") | (data_fig3ef$TNM_4level=="IV")) & !is.na(data_fig3ef$TIMES))  #68
data_fig3ef_grp <- data_fig3ef[grp_idx,]
data_fig3ef_grp$TNM_4level <- factor(data_fig3ef_grp$TNM_4level,levels = c("I","IV"))
table(data_fig3ef_grp$TNM_4level)
# I IV 
# 53 15 
res.cox <- coxph(Surv(DFS, Recurrence) ~ TNM_4level, data = data_fig3ef_grp, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TNM_4level, data = data_fig3ef_grp, 
#         id = patient_id)
# 
# n= 37, number of events= 16 
# (31 observations deleted due to missingness)
# 
# 1:2              coef exp(coef) se(coef)     z Pr(>|z|)
# TNM_4levelIV 1.5034    4.4967   0.8113 1.853   0.0639
# 
# 1:2            exp(coef) exp(-coef) lower .95 upper .95
# TNM_4levelIV     4.497     0.2224    0.9169     22.05
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.564  (se = 0.049 )
# Likelihood ratio test= 2.57  on 1 df,   p=0.1
# Wald test            = 3.43  on 1 df,   p=0.06
# Score (logrank) test = 4.12  on 1 df,   p=0.04

grp_idx <- which(((data_fig3ef$TNM_4level=="II") | (data_fig3ef$TNM_4level=="III")) & !is.na(data_fig3ef$TIMES))  #87
data_fig3ef_grp <- data_fig3ef[grp_idx,]
data_fig3ef_grp$TNM_4level <- factor(data_fig3ef_grp$TNM_4level,levels = c("II","III"))
table(data_fig3ef_grp$TNM_4level)
# II III 
# 71  16 
res.cox <- coxph(Surv(DFS, Recurrence) ~ TNM_4level, data = data_fig3ef_grp, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TNM_4level, data = data_fig3ef_grp, 
#         id = patient_id)
# 
# n= 61, number of events= 32 
# (26 observations deleted due to missingness)
# 
# 1:2               coef exp(coef) se(coef)    z Pr(>|z|)
# TNM_4levelIII 0.8889    2.4324   0.4605 1.93   0.0536
# 
# 1:2             exp(coef) exp(-coef) lower .95 upper .95
# TNM_4levelIII     2.432     0.4111    0.9864     5.998
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.554  (se = 0.037 )
# Likelihood ratio test= 3.11  on 1 df,   p=0.08
# Wald test            = 3.73  on 1 df,   p=0.05
# Score (logrank) test = 3.97  on 1 df,   p=0.05

grp_idx <- which(((data_fig3ef$TNM_4level=="II") | (data_fig3ef$TNM_4level=="IV")) & !is.na(data_fig3ef$TIMES))  #86
data_fig3ef_grp <- data_fig3ef[grp_idx,]
data_fig3ef_grp$TNM_4level <- factor(data_fig3ef_grp$TNM_4level,levels = c("II","IV"))
table(data_fig3ef_grp$TNM_4level)
# II IV 
# 71 15
res.cox <- coxph(Surv(DFS, Recurrence) ~ TNM_4level, data = data_fig3ef_grp, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TNM_4level, data = data_fig3ef_grp, 
#         id = patient_id)
# 
# n= 57, number of events= 28 
# (29 observations deleted due to missingness)
# 
# 1:2              coef exp(coef) se(coef)     z Pr(>|z|)
# TNM_4levelIV 0.8556    2.3529   0.7392 1.157    0.247
# 
# 1:2            exp(coef) exp(-coef) lower .95 upper .95
# TNM_4levelIV     2.353      0.425    0.5525     10.02
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.532  (se = 0.032 )
# Likelihood ratio test= 1.06  on 1 df,   p=0.3
# Wald test            = 1.34  on 1 df,   p=0.2
# Score (logrank) test = 1.42  on 1 df,   p=0.2

grp_idx <- which(((data_fig3ef$TNM_4level=="III") | (data_fig3ef$TNM_4level=="IV")) & !is.na(data_fig3ef$TIMES))  #31
data_fig3ef_grp <- data_fig3ef[grp_idx,]
data_fig3ef_grp$TNM_4level <- factor(data_fig3ef_grp$TNM_4level,levels = c("III","IV"))
table(data_fig3ef_grp$TNM_4level)
# III  IV 
# 16  15 
res.cox <- coxph(Surv(DFS, Recurrence) ~ TNM_4level, data = data_fig3ef_grp, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ TNM_4level, data = data_fig3ef_grp, 
#         id = patient_id)
# 
# n= 10, number of events= 8 
# (21 observations deleted due to missingness)
# 
# 1:2              coef exp(coef) se(coef)     z Pr(>|z|)
# TNM_4levelIV 0.2700    1.3099   0.8666 0.312    0.755
# 
# 1:2            exp(coef) exp(-coef) lower .95 upper .95
# TNM_4levelIV      1.31     0.7634    0.2397      7.16
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.541  (se = 0.109 )
# Likelihood ratio test= 0.09  on 1 df,   p=0.8
# Wald test            = 0.1  on 1 df,   p=0.8
# Score (logrank) test = 0.1  on 1 df,   p=0.8


##(2) BCLC stratified ##
table(data_fig3ef$BCLC_3level)
# A  B  C 
# 26 19 85
l1_idx <- which((data_fig3ef$BCLC_3level=="A") & !is.na(data_fig3ef$TIMES))  #26
l2_idx <- which((data_fig3ef$BCLC_3level=="B") & !is.na(data_fig3ef$TIMES))  #19
l3_idx <- which((data_fig3ef$BCLC_3level=="C") & !is.na(data_fig3ef$TIMES))  #85
l1_dat <- data_fig3ef$TIMES[l1_idx]
l2_dat <- data_fig3ef$TIMES[l2_idx]
l3_dat <- data_fig3ef$TIMES[l3_idx]

plot(survfit(formula = hcc_sur ~ BCLC_3level, data=data_fig3ef), lwd=3, conf.int = F, mark.time = T, bty = "n", las = 1, xlim=c(0,40), cex = 2,
     ylab = "Disease free survival (%)", xlab = "Time (months)", main = "BCLC", col=c(adjustcolor("#1B9E77", alpha.f=1), adjustcolor("#D95F02", alpha.f=1), adjustcolor("#7570B3", alpha.f=1)))
legend("bottomleft", legend=c("A", "B", "C"), col=c(adjustcolor("#1B9E77", alpha.f=1), adjustcolor("#D95F02", alpha.f=1), adjustcolor("#7570B3", alpha.f=1)), bg=NA, lty=1, lwd=4, box.lty=0, cex=1.3)
#save as /data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/ExFig5c_survival_BCLC.pdf : survival curves in BCLC stratified subgroups

##cox test
grp_idx <- which(((data_fig3ef$BCLC_3level=="A") | (data_fig3ef$BCLC_3level=="B")) & !is.na(data_fig3ef$TIMES))  #45
data_fig3ef_grp <- data_fig3ef[grp_idx,]
data_fig3ef_grp$BCLC_3level <- factor(data_fig3ef_grp$BCLC_3level,levels = c("A","B"))
table(data_fig3ef_grp$BCLC_3level)
# A  B 
# 26 19 
res.cox <- coxph(Surv(DFS, Recurrence) ~ BCLC_3level, data = data_fig3ef_grp, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ BCLC_3level, data = data_fig3ef_grp, 
#         id = patient_id)
# 
# n= 39, number of events= 16 
# (6 observations deleted due to missingness)
# 
# 1:2              coef exp(coef) se(coef)     z Pr(>|z|)
# BCLC_3levelB 1.1136    3.0452   0.5222 2.132    0.033
# 
# 1:2            exp(coef) exp(-coef) lower .95 upper .95
# BCLC_3levelB     3.045     0.3284     1.094     8.475
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.701  (se = 0.052 )
# Likelihood ratio test= 4.75  on 1 df,   p=0.03
# Wald test            = 4.55  on 1 df,   p=0.03
# Score (logrank) test = 5  on 1 df,   p=0.03

grp_idx <- which(((data_fig3ef$BCLC_3level=="A") | (data_fig3ef$BCLC_3level=="C")) & !is.na(data_fig3ef$TIMES))  #111
data_fig3ef_grp <- data_fig3ef[grp_idx,]
data_fig3ef_grp$BCLC_3level <- factor(data_fig3ef_grp$BCLC_3level,levels = c("A","C"))
table(data_fig3ef_grp$BCLC_3level)
# A  C 
# 26 85 
res.cox <- coxph(Surv(DFS, Recurrence) ~ BCLC_3level, data = data_fig3ef_grp, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ BCLC_3level, data = data_fig3ef_grp, 
#         id = patient_id)
# 
# n= 82, number of events= 38 
# (29 observations deleted due to missingness)
# 
# 1:2             coef exp(coef) se(coef)    z Pr(>|z|)
# BCLC_3levelC 1.213     3.363    0.456 2.66  0.00782
# 
# 1:2            exp(coef) exp(-coef) lower .95 upper .95
# BCLC_3levelC     3.363     0.2973     1.376     8.221
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.637  (se = 0.026 )
# Likelihood ratio test= 8.87  on 1 df,   p=0.003
# Wald test            = 7.07  on 1 df,   p=0.008
# Score (logrank) test = 7.85  on 1 df,   p=0.005

grp_idx <- which(((data_fig3ef$BCLC_3level=="B") | (data_fig3ef$BCLC_3level=="C")) & !is.na(data_fig3ef$TIMES))  #104
data_fig3ef_grp <- data_fig3ef[grp_idx,]
data_fig3ef_grp$BCLC_3level <- factor(data_fig3ef_grp$BCLC_3level,levels = c("B","C"))
table(data_fig3ef_grp$BCLC_3level)
# B  C 
# 19 85 
res.cox <- coxph(Surv(DFS, Recurrence) ~ BCLC_3level, data = data_fig3ef_grp, id = patient_id)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Recurrence) ~ BCLC_3level, data = data_fig3ef_grp, 
#         id = patient_id)
# 
# n= 79, number of events= 42 
# (25 observations deleted due to missingness)
# 
# 1:2               coef exp(coef) se(coef)      z Pr(>|z|)
# BCLC_3levelC -0.0135    0.9866   0.3668 -0.037    0.971
# 
# 1:2            exp(coef) exp(-coef) lower .95 upper .95
# BCLC_3levelC    0.9866      1.014    0.4807     2.025
# 
# States: 1= (s0), 2= Y 
# 
# Concordance= 0.51  (se = 0.038 )
# Likelihood ratio test= 0  on 1 df,   p=1
# Wald test            = 0  on 1 df,   p=1
# Score (logrank) test = 0  on 1 df,   p=1

##combine the two figures above and generate Extended Data Fig. 5c

