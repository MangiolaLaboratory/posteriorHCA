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

#######To generate barplot: Fig. 3e and Extended Data Fig. 5b #######
#### Fig. 3e: TIMES in REC and non-REC in three types of stratification ####
##(1) unstratified ##
no_idx <- which((data_fig3ef$Recurrence=="N") & !is.na(data_fig3ef$TIMES))  #91
yes_idx <- which((data_fig3ef$Recurrence=="Y") & !is.na(data_fig3ef$TIMES))  #72
no_dat <- data_fig3ef$TIMES[no_idx]
yes_dat <- data_fig3ef$TIMES[yes_idx]

t_res_rec <- t.test(no_dat, yes_dat, alternative = c("two.sided"), mu=0, paired=F, var.equal=F)
# Welch Two Sample t-test
# 
# data:  no_dat and yes_dat
# t = -9.0618, df = 147.31, p-value = 7.295e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.3959901 -0.2541973
# sample estimates:
#   mean of x mean of y 
# 0.2911220 0.6162157 

sum_dat <- data.frame(name=c("no","yes"),mean=c(mean(no_dat), mean(yes_dat)),sd=c(sd(no_dat), sd(yes_dat)),n=c(length(no_dat), length(yes_dat)),stringsAsFactors = F)
sum_dat$se <- sum_dat$sd / sqrt(sum_dat$n)
barCenters <- barplot(height = sum_dat$mean, width=0.5, space=0.1,
                      las = 2, col=c(adjustcolor("#559AC6", alpha.f=0.5), adjustcolor("#E73334", alpha.f=0.5)),
                      ylim = c(0, 1), xlim=c(0,2),
                      cex.names = 0.75, xaxt = "n",
                      ylab = "TIMES",
                      border = "black", axes = TRUE)
text(x = c(0.65,1.05), y = -0.05, srt = 0, adj = 1, labels = c("non-REC","REC"), xpd = TRUE)
stripchart(x=no_dat,pch=19, cex=1.5, at=c(0.3), add=T,vertical=T,method = "jitter",col=adjustcolor("#559AC6", alpha.f=0.95))
stripchart(x=yes_dat,pch=19, cex=1.5, at=c(0.85), add=T,vertical=T,method = "jitter",col=adjustcolor("#E73334", alpha.f=0.95))
segments(barCenters, sum_dat$mean - sum_dat$se * 1.96, barCenters,sum_dat$mean + sum_dat$se * 1.96, lwd = 1.5)
arrows(barCenters, sum_dat$mean - sum_dat$se * 1.96, barCenters,sum_dat$mean + sum_dat$se * 1.96, lwd = 1.5, angle = 90, code = 3, length = 0.05)
##save as /data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/Fig3e_unstratified.pdf : barplot of times in unstratified group

##(2.1) in TNM <II subgroup ##
grp_idx <- which((data_fig3ef$TNM_2level=="<II") & !is.na(data_fig3ef$TIMES))  #53
data_fig3ef_grp <- data_fig3ef[grp_idx,]

no_idx <- which((data_fig3ef_grp$Recurrence=="N") & !is.na(data_fig3ef_grp$TIMES))  #22
yes_idx <- which((data_fig3ef_grp$Recurrence=="Y") & !is.na(data_fig3ef_grp$TIMES))  #16
no_dat <- data_fig3ef_grp$TIMES[no_idx]
yes_dat <- data_fig3ef_grp$TIMES[yes_idx]

t_res_rec <- t.test(no_dat, yes_dat, alternative = c("two.sided"), mu=0, paired=F, var.equal=F)
# Welch Two Sample t-test
# 
# data:  no_dat and yes_dat
# t = -5.602, df = 25.768, p-value = 7.132e-06
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.4907071 -0.2271798
# sample estimates:
#   mean of x mean of y 
# 0.2538666 0.6128100 

sum_dat <- data.frame(name=c("no","yes"),mean=c(mean(no_dat), mean(yes_dat)),sd=c(sd(no_dat), sd(yes_dat)),n=c(length(no_dat), length(yes_dat)),stringsAsFactors = F)
sum_dat$se <- sum_dat$sd / sqrt(sum_dat$n)
barCenters <- barplot(height = sum_dat$mean, width=0.5, space=0.1,
                      las = 2, col=c(adjustcolor("#559AC6", alpha.f=0.5), adjustcolor("#E73334", alpha.f=0.5)),
                      ylim = c(0, 1), xlim=c(0,2),
                      cex.names = 0.75, xaxt = "n",
                      ylab = "TIMES",
                      border = "black", axes = TRUE)
text(x = c(0.65,1.05), y = -0.05, srt = 0, adj = 1, labels = c("non-REC","REC"), xpd = TRUE)
stripchart(x=no_dat,pch=19, cex=1.5, at=c(0.3), add=T,vertical=T,method = "jitter",col=adjustcolor("#559AC6", alpha.f=0.95))
stripchart(x=yes_dat,pch=19, cex=1.5, at=c(0.85), add=T,vertical=T,method = "jitter",col=adjustcolor("#E73334", alpha.f=0.95))
segments(barCenters, sum_dat$mean - sum_dat$se * 1.96, barCenters,sum_dat$mean + sum_dat$se * 1.96, lwd = 1.5)
arrows(barCenters, sum_dat$mean - sum_dat$se * 1.96, barCenters,sum_dat$mean + sum_dat$se * 1.96, lwd = 1.5, angle = 90, code = 3, length = 0.05)
##save as /data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/Fig3e_TNM_low.pdf : barplot of times in TNM <II subgroup

##(2.2) in TNM ≥II subgroup ##
grp_idx <- which((data_fig3ef$TNM_2level=="≥II") & !is.na(data_fig3ef$TIMES))  #102
data_fig3ef_grp <- data_fig3ef[grp_idx,]

no_idx <- which((data_fig3ef_grp$Recurrence=="N") & !is.na(data_fig3ef_grp$TIMES))  #46
yes_idx <- which((data_fig3ef_grp$Recurrence=="Y") & !is.na(data_fig3ef_grp$TIMES))  #43
no_dat <- data_fig3ef_grp$TIMES[no_idx]
yes_dat <- data_fig3ef_grp$TIMES[yes_idx]

t_res_rec <- t.test(no_dat, yes_dat, alternative = c("two.sided"), mu=0, paired=F, var.equal=F)
# Welch Two Sample t-test
# 
# data:  no_dat and yes_dat
# t = -8.8945, df = 86.98, p-value = 7.4e-14
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.4839378 -0.3071551
# sample estimates:
#   mean of x mean of y 
# 0.2924734 0.6880198 

sum_dat <- data.frame(name=c("no","yes"),mean=c(mean(no_dat), mean(yes_dat)),sd=c(sd(no_dat), sd(yes_dat)),n=c(length(no_dat), length(yes_dat)),stringsAsFactors = F)
sum_dat$se <- sum_dat$sd / sqrt(sum_dat$n)
barCenters <- barplot(height = sum_dat$mean, width=0.5, space=0.1,
                      las = 2, col=c(adjustcolor("#559AC6", alpha.f=0.5), adjustcolor("#E73334", alpha.f=0.5)),
                      ylim = c(0, 1), xlim=c(0,2),
                      cex.names = 0.75, xaxt = "n",
                      ylab = "TIMES",
                      border = "black", axes = TRUE)
text(x = c(0.65,1.05), y = -0.05, srt = 0, adj = 1, labels = c("non-REC","REC"), xpd = TRUE)
stripchart(x=no_dat,pch=19, cex=1.5, at=c(0.3), add=T,vertical=T,method = "jitter",col=adjustcolor("#559AC6", alpha.f=0.95))
stripchart(x=yes_dat,pch=19, cex=1.5, at=c(0.85), add=T,vertical=T,method = "jitter",col=adjustcolor("#E73334", alpha.f=0.95))
segments(barCenters, sum_dat$mean - sum_dat$se * 1.96, barCenters,sum_dat$mean + sum_dat$se * 1.96, lwd = 1.5)
arrows(barCenters, sum_dat$mean - sum_dat$se * 1.96, barCenters,sum_dat$mean + sum_dat$se * 1.96, lwd = 1.5, angle = 90, code = 3, length = 0.05)
##save as /data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/Fig3e_TNM_high.pdf: barplot of times in TNM ≥II subgroup

##(3.1) in BCLC: A subgroup ##
grp_idx <- which((data_fig3ef$BCLC_A_vs_BC=="low") & !is.na(data_fig3ef$TIMES))  #26
data_fig3ef_grp <- data_fig3ef[grp_idx,]

no_idx <- which((data_fig3ef_grp$Recurrence=="N") & !is.na(data_fig3ef_grp$TIMES))  #17
yes_idx <- which((data_fig3ef_grp$Recurrence=="Y") & !is.na(data_fig3ef_grp$TIMES))  #9
no_dat <- data_fig3ef_grp$TIMES[no_idx]
yes_dat <- data_fig3ef_grp$TIMES[yes_idx]

t_res_rec <- t.test(no_dat, yes_dat, alternative = c("two.sided"), mu=0, paired=F, var.equal=F)
# Welch Two Sample t-test
# 
# data:  no_dat and yes_dat
# t = -2.3517, df = 12.324, p-value = 0.03609
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.39690685 -0.01572557
# sample estimates:
#   mean of x mean of y 
# 0.3088296 0.5151458 

sum_dat <- data.frame(name=c("no","yes"),mean=c(mean(no_dat), mean(yes_dat)),sd=c(sd(no_dat), sd(yes_dat)),n=c(length(no_dat), length(yes_dat)),stringsAsFactors = F)
sum_dat$se <- sum_dat$sd / sqrt(sum_dat$n)
barCenters <- barplot(height = sum_dat$mean, width=0.5, space=0.1,
                      las = 2, col=c(adjustcolor("#559AC6", alpha.f=0.5), adjustcolor("#E73334", alpha.f=0.5)),
                      ylim = c(0, 1), xlim=c(0,2),
                      cex.names = 0.75, xaxt = "n",
                      ylab = "TIMES",
                      border = "black", axes = TRUE)
text(x = c(0.65,1.05), y = -0.05, srt = 0, adj = 1, labels = c("non-REC","REC"), xpd = TRUE)
stripchart(x=no_dat,pch=19, cex=1.5, at=c(0.3), add=T,vertical=T,method = "jitter",col=adjustcolor("#559AC6", alpha.f=0.95))
stripchart(x=yes_dat,pch=19, cex=1.5, at=c(0.85), add=T,vertical=T,method = "jitter",col=adjustcolor("#E73334", alpha.f=0.95))
segments(barCenters, sum_dat$mean - sum_dat$se * 1.96, barCenters,sum_dat$mean + sum_dat$se * 1.96, lwd = 1.5)
arrows(barCenters, sum_dat$mean - sum_dat$se * 1.96, barCenters,sum_dat$mean + sum_dat$se * 1.96, lwd = 1.5, angle = 90, code = 3, length = 0.05)
##save as /data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/Fig3e_BCLC_low.pdf: barplot of times in BCLC =A subgroup

##(3.2) in BCLC: B+C subgroup ##
grp_idx <- which((data_fig3ef$BCLC_A_vs_BC=="high") & !is.na(data_fig3ef$TIMES))  #104
data_fig3ef_grp <- data_fig3ef[grp_idx,]

no_idx <- which((data_fig3ef_grp$Recurrence=="N") & !is.na(data_fig3ef_grp$TIMES))  #53
yes_idx <- which((data_fig3ef_grp$Recurrence=="Y") & !is.na(data_fig3ef_grp$TIMES))  #51
no_dat <- data_fig3ef_grp$TIMES[no_idx]
yes_dat <- data_fig3ef_grp$TIMES[yes_idx]

t_res_rec <- t.test(no_dat, yes_dat, alternative = c("two.sided"), mu=0, paired=F, var.equal=F)
# Welch Two Sample t-test
# 
# data:  no_dat and yes_dat
# t = -10.637, df = 101.93, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.4971681 -0.3408944
# sample estimates:
#   mean of x mean of y 
# 0.2699268 0.6889580 

sum_dat <- data.frame(name=c("no","yes"),mean=c(mean(no_dat), mean(yes_dat)),sd=c(sd(no_dat), sd(yes_dat)),n=c(length(no_dat), length(yes_dat)),stringsAsFactors = F)
sum_dat$se <- sum_dat$sd / sqrt(sum_dat$n)
barCenters <- barplot(height = sum_dat$mean, width=0.5, space=0.1,
                      las = 2, col=c(adjustcolor("#559AC6", alpha.f=0.5), adjustcolor("#E73334", alpha.f=0.5)),
                      ylim = c(0, 1), xlim=c(0,2),
                      cex.names = 0.75, xaxt = "n",
                      ylab = "TIMES",
                      border = "black", axes = TRUE)
text(x = c(0.65,1.05), y = -0.05, srt = 0, adj = 1, labels = c("non-REC","REC"), xpd = TRUE)
stripchart(x=no_dat,pch=19, cex=1.5, at=c(0.3), add=T,vertical=T,method = "jitter",col=adjustcolor("#559AC6", alpha.f=0.95))
stripchart(x=yes_dat,pch=19, cex=1.5, at=c(0.85), add=T,vertical=T,method = "jitter",col=adjustcolor("#E73334", alpha.f=0.95))
segments(barCenters, sum_dat$mean - sum_dat$se * 1.96, barCenters,sum_dat$mean + sum_dat$se * 1.96, lwd = 1.5)
arrows(barCenters, sum_dat$mean - sum_dat$se * 1.96, barCenters,sum_dat$mean + sum_dat$se * 1.96, lwd = 1.5, angle = 90, code = 3, length = 0.05)
##save as /data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/Fig3e_BCLC_high.pdf: barplot of times in BCLC =B+C subgroup

##combine the five subfigures above and generate Fig. 3e


#### Extended Data Fig. 5b: TIMES in TNM subgroups ####
##TNM stratified ##
l1_idx <- which((data_fig3ef$TNM_4level=="I") & !is.na(data_fig3ef$TIMES))  #53
l2_idx <- which((data_fig3ef$TNM_4level=="II") & !is.na(data_fig3ef$TIMES))  #71
l3_idx <- which((data_fig3ef$TNM_4level=="III") & !is.na(data_fig3ef$TIMES))  #16
l4_idx <- which((data_fig3ef$TNM_4level=="IV") & !is.na(data_fig3ef$TIMES))  #15
l12_idx <- c(l1_idx, l2_idx)  #124
l34_idx <- c(l3_idx, l4_idx)  #31
l12_dat <- data_fig3ef$TIMES[l12_idx]
l34_dat <- data_fig3ef$TIMES[l34_idx]
t.test(l12_dat, l34_dat, alternative = c("two.sided"), mu=0, paired=F, var.equal=F)  #use
# Welch Two Sample t-test
# 
# data:  l12_dat and l34_dat
# t = -2.107, df = 48.644, p-value = 0.0403
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.207526398 -0.004892161
# sample estimates:
#   mean of x mean of y 
# 0.4488561 0.5550654 

sum_dat <- data.frame(name=c("Stage1-2","Stage3-4"),mean=c(mean(l12_dat), mean(l34_dat)),sd=c(sd(l12_dat), sd(l34_dat)),n=c(length(l12_dat), length(l34_dat)),stringsAsFactors = F)
sum_dat$se <- sum_dat$sd / sqrt(sum_dat$n)
barCenters <- barplot(height = sum_dat$mean, width=0.5, space=0.1,
                      las = 2, col=c(adjustcolor("#559AC6", alpha.f=0.5), adjustcolor("#E73334", alpha.f=0.5)),
                      ylim = c(0, 1), xlim=c(0,2),
                      cex.names = 0.75, xaxt = "n",
                      ylab = "TIMES",
                      border = "black", axes = TRUE)
text(x = c(0.65,1.05), y = -0.05, srt = 0, adj = 1, labels = c("Stage1-2","Stage3-4"), xpd = TRUE)
stripchart(x=l12_dat,pch=19, cex=1.5, at=c(0.3), add=T,vertical=T,method = "jitter",col=adjustcolor("#559AC6", alpha.f=0.95))
stripchart(x=l34_dat,pch=19, cex=1.5, at=c(0.85), add=T,vertical=T,method = "jitter",col=adjustcolor("#E73334", alpha.f=0.95))
segments(barCenters, sum_dat$mean - sum_dat$se * 1.96, barCenters,sum_dat$mean + sum_dat$se * 1.96, lwd = 1.5)
arrows(barCenters, sum_dat$mean - sum_dat$se * 1.96, barCenters,sum_dat$mean + sum_dat$se * 1.96, lwd = 1.5, angle = 90, code = 3, length = 0.05)

##save as /data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/ExFig5b_barplot_times_in_TNM12vs34.pdf: barplot of times in TNM subgroups

