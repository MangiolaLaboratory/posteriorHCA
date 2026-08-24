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
#   [1] pROC_1.18.5     dplyr_1.1.4     survival_3.1-12
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
library(pROC)

setwd("/code/3_Comparison_between_clinical_markers_and_TIMES/")

#######load clinic data####
load("/input/Input_3_Comparison_between_clinical_markers_and_TIMES/clinic_data_input.RData")

####to generate Fig3d. non-REC vs REC ####
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

b4x <- all_data$TIMES[which((all_data$Recurrence=="N") & !is.na(all_data$TIMES))]  #91
x <- all_data$TIMES[which((all_data$Recurrence=="Y") & !is.na(all_data$TIMES))]  #72

wil_t <- wilcox.test(x, b4x,alternative = c("two.sided"),mu = 0,paired=F,exact=F,correct=F,conf.int=T,conf.level=0.95)
# Wilcoxon rank sum test
# 
# data:  x and b4x
# W = 5457, p-value = 3.135e-13
# alternative hypothesis: true location shift is not equal to 0
# 95 percent confidence interval:
#   0.3070558 0.4388399
# sample estimates:
#   difference in location 
# 0.3748984 

xmin <- -0.1; xmax <- 1.1; break_seq <- seq(from=xmin,to=xmax,by=0.1); ymin <- 0; ymax <- 3.5
my.hist <- hist(x, xlim=c(xmin,xmax), ylim=c(ymin,ymax),main="", xlab="", ylab="",breaks=break_seq,freq=F,lty="blank",las=1,col=NULL)
lines(c(my.hist$breaks, max(my.hist$breaks)),c(0,my.hist$density,0), type='S',col="#CC4B40",lwd=4)
my.hist <- hist(b4x, add=T, breaks=break_seq,freq=F,lty="blank",col=NULL)
lines(c(my.hist$breaks, max(my.hist$breaks)),c(0,my.hist$density,0), type='S',col="#2367AC",lwd=4)

text(0.5,3,labels=paste("P = ",formatC(wil_t$p.value,format = "e",digits = 2),sep=""))
text(0.5,3.3,labels=paste("MD = ",round(wil_t$estimate,3),sep=""))
##save as /data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/Fig3d_recurrence_wilcox.pdf

##compute accuracy and specificity
length(which(b4x<0.5))  #78
length(which(b4x>=0.5))  #13
length(which(x<0.5))  #16
length(which(x>=0.5))  #56

(78+56)/(78+13+16+56)
#0.8220859  (Accuracy)

78/(78+13)
#0.8571429 (Specificity)

56/(56+16)
#0.7777778 (Sensitivity)

56/(56+13)
#0.8115942 (Positive Predictive Value)

78/(78+16)
#0.8297872 (Negative Predictive Value)

##to generate Extended Figure 5e: Distributions of TIMES scores for anti-PD1 immunotherapy recipients
load("/input/Input_3_Comparison_between_clinical_markers_and_TIMES/clinic_data_input.RData")

table(all_data$Immunotherapy)
# NA  PD (non-response)    PR (response)
# 270  12                  13 

b4x <- all_data$TIMES[which(all_data$Immunotherapy=="PD")]
x <- all_data$TIMES[which(all_data$Immunotherapy=="PR")]

wil_t <- wilcox.test(x, b4x,alternative = c("two.sided"),mu = 0,paired=F,exact=F,correct=F,conf.int=T,conf.level=0.95)
# Wilcoxon rank sum test
# 
# data:  x and b4x
# W = 136, p-value = 0.001606
# alternative hypothesis: true location shift is not equal to 0
# 95 percent confidence interval:
#   0.1240579 0.4038181
# sample estimates:
#   difference in location 
# 0.2673768 

xmin <- -0.1; xmax <- 1.1; break_seq <- seq(from=xmin,to=xmax,by=0.15); ymin <- 0; ymax <- 4
my.hist <- hist(x, xlim=c(xmin,xmax), ylim=c(ymin,ymax),main="", xlab="", ylab="",breaks=break_seq,freq=F,lty="blank",las=1,col=NULL)
lines(c(my.hist$breaks, max(my.hist$breaks)),c(0,my.hist$density,0), type='S',col="#CC4B40",lwd=4)
my.hist <- hist(b4x, add=T, breaks=break_seq,freq=F,lty="blank",col=NULL)
lines(c(my.hist$breaks, max(my.hist$breaks)),c(0,my.hist$density,0), type='S',col="#2367AC",lwd=4)

text(0.5,3.5,labels=paste("P = ",formatC(wil_t$p.value,format = "e",digits = 2),sep=""))
text(0.5,4,labels=paste("MD = ",round(wil_t$estimate,3),sep=""))
##/data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/ExFig5e_immunotherpy_wilcox.pdf


##to generate Extended Figure 5f: Receiver operating characteristic (ROC) curves for predicting recurrence outcome
load("/input/Input_3_Comparison_between_clinical_markers_and_TIMES/nk_data_input.RData")

summary(dat_stain_multi$TIMES)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.06032 0.15871 0.49840 0.48257 0.80545 0.94822 


#1. "(CD3-CD56+)%"
##univariate analysis
##cd57 at IF
res.cox <- coxph(Surv(DFS, Status) ~ `(CD3-CD56+)%`, data = dat_stain_multi, id = `patient ID`)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Status) ~ `(CD3-CD56+)%`, data = dat_stain_multi, 
#         id = `patient ID`)
# 
# n= 61, number of events= 30 
# 
# 1:2                 coef exp(coef) se(coef)      z Pr(>|z|)
# `(CD3-CD56+)%` -0.5285    0.5895   0.2848 -1.855   0.0635
# 
# 1:2              exp(coef) exp(-coef) lower .95 upper .95
# `(CD3-CD56+)%`    0.5895      1.696    0.3373      1.03
# 
# States: 1= (s0), 2= REC 
# 
# Concordance= 0.599  (se = 0.054 )
# Likelihood ratio test= 4.08  on 1 df,   p=0.04
# Wald test            = 3.44  on 1 df,   p=0.06
# Score (logrank) test = 3.56  on 1 df,   p=0.06

mylogit <- glm(Status ~ `(CD3-CD56+)%`, data = dat_stain_multi,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Status ~ `(CD3-CD56+)%`, family = binomial, data = dat_stain_multi)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.5032  -1.2004  -0.6135   1.1010   1.7669  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)      0.7837     0.4721   1.660   0.0969 .
# `(CD3-CD56+)%`  -0.7520     0.3735  -2.013   0.0441 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 84.548  on 60  degrees of freedom
# Residual deviance: 79.901  on 59  degrees of freedom
# AIC: 83.901
# 
# Number of Fisher Scoring iterations: 4


##cd57 at TC
res.cox <- coxph(Surv(DFS, Status) ~ `(CD3-CD56+)%-tc`, data = dat_stain_multi, id = `patient ID`)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Status) ~ `(CD3-CD56+)%-tc`, data = dat_stain_multi, 
#         id = `patient ID`)
# 
# n= 61, number of events= 30 
# 
# 1:2                    coef exp(coef) se(coef)      z Pr(>|z|)
# `(CD3-CD56+)%-tc` -0.2886    0.7493   0.2340 -1.233    0.217
# 
# 1:2                 exp(coef) exp(-coef) lower .95 upper .95
# `(CD3-CD56+)%-tc`    0.7493      1.335    0.4737     1.185
# 
# States: 1= (s0), 2= REC 
# 
# Concordance= 0.602  (se = 0.057 )
# Likelihood ratio test= 2.85  on 1 df,   p=0.09
# Wald test            = 1.52  on 1 df,   p=0.2
# Score (logrank) test = 1.45  on 1 df,   p=0.2

mylogit <- glm(Status ~ `(CD3-CD56+)%-tc`, data = dat_stain_multi,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Status ~ `(CD3-CD56+)%-tc`, family = binomial, 
#       data = dat_stain_multi)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.28876  -1.22346  -0.07129   1.10285   1.59949  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)         0.2651     0.3221   0.823    0.410
# `(CD3-CD56+)%-tc`  -0.3717     0.2813  -1.322    0.186
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 84.548  on 60  degrees of freedom
# Residual deviance: 81.216  on 59  degrees of freedom
# AIC: 85.216
# 
# Number of Fisher Scoring iterations: 5


##multivariate analysis
res.cox <- coxph(Surv(DFS, Status) ~ TIMES+`(CD3-CD56+)%`+`(CD3-CD56+)%-tc`, data = dat_stain_multi, id = `patient ID`)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Status) ~ TIMES + `(CD3-CD56+)%` + 
#           `(CD3-CD56+)%-tc`, data = dat_stain_multi, id = `patient ID`)
# 
# n= 61, number of events= 30 
# 
# 1:2                       coef  exp(coef)   se(coef)      z Pr(>|z|)
# TIMES                 9.4308 12466.0093     1.9940  4.730 2.25e-06
# `(CD3-CD56+)%`       -0.3699     0.6908     0.3615 -1.023    0.306
# `(CD3-CD56+)%-tc`    -0.2379     0.7883     0.2747 -0.866    0.386
# 
# 1:2                 exp(coef) exp(-coef) lower .95 upper .95
# TIMES             1.247e+04  8.022e-05  250.2793 6.209e+05
# `(CD3-CD56+)%`    6.908e-01  1.448e+00    0.3402 1.403e+00
# `(CD3-CD56+)%-tc` 7.883e-01  1.269e+00    0.4601 1.350e+00
# 
# States: 1= (s0), 2= REC 
# 
# Concordance= 0.879  (se = 0.029 )
# Likelihood ratio test= 80.13  on 3 df,   p=<2e-16
# Wald test            = 23.78  on 3 df,   p=3e-05
# Score (logrank) test = 68.02  on 3 df,   p=1e-14


#2. "(CD3-CD57+)%"
##univariate analysis
##cd57 at IF
res.cox <- coxph(Surv(DFS, Status) ~ `(CD3-CD57+)%`, data = dat_stain_multi, id = `patient ID`)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Status) ~ `(CD3-CD57+)%`, data = dat_stain_multi, 
#         id = `patient ID`)
# 
# n= 61, number of events= 30 
# 
# 1:2                 coef exp(coef) se(coef)      z Pr(>|z|)
# `(CD3-CD57+)%` -1.2300    0.2923   0.4201 -2.928  0.00341
# 
# 1:2              exp(coef) exp(-coef) lower .95 upper .95
# `(CD3-CD57+)%`    0.2923      3.421    0.1283    0.6658
# 
# States: 1= (s0), 2= REC 
# 
# Concordance= 0.652  (se = 0.054 )
# Likelihood ratio test= 12.96  on 1 df,   p=3e-04
# Wald test            = 8.57  on 1 df,   p=0.003
# Score (logrank) test = 8.14  on 1 df,   p=0.004

mylogit <- glm(Status ~ `(CD3-CD57+)%`, data = dat_stain_multi,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Status ~ `(CD3-CD57+)%`, family = binomial, data = dat_stain_multi)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.88410  -1.02113  -0.02935   0.89971   1.52684  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)      2.0177     0.7290   2.768  0.00564 **
#   `(CD3-CD57+)%`  -2.2664     0.7949  -2.851  0.00436 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 84.548  on 60  degrees of freedom
# Residual deviance: 68.708  on 59  degrees of freedom
# AIC: 72.708
# 
# Number of Fisher Scoring iterations: 5


##cd57 at TC
res.cox <- coxph(Surv(DFS, Status) ~ `(CD3-CD57+)%-tc`, data = dat_stain_multi, id = `patient ID`)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Status) ~ `(CD3-CD57+)%-tc`, data = dat_stain_multi, 
#         id = `patient ID`)
# 
# n= 61, number of events= 30 
# 
# 1:2                    coef exp(coef) se(coef)      z Pr(>|z|)
# `(CD3-CD57+)%-tc` -0.1486    0.8619   0.1052 -1.412    0.158
# 
# 1:2                 exp(coef) exp(-coef) lower .95 upper .95
# `(CD3-CD57+)%-tc`    0.8619       1.16    0.7012     1.059
# 
# States: 1= (s0), 2= REC 
# 
# Concordance= 0.534  (se = 0.054 )
# Likelihood ratio test= 2.28  on 1 df,   p=0.1
# Wald test            = 1.99  on 1 df,   p=0.2
# Score (logrank) test = 2.03  on 1 df,   p=0.2

mylogit <- glm(Status ~ `(CD3-CD57+)%-tc`, data = dat_stain_multi,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Status ~ `(CD3-CD57+)%-tc`, family = binomial, 
#       data = dat_stain_multi)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.3058  -1.1988  -0.7242   1.1343   1.4156  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)         0.3084     0.3654   0.844    0.399
# `(CD3-CD57+)%-tc`  -0.1781     0.1381  -1.290    0.197
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 84.548  on 60  degrees of freedom
# Residual deviance: 82.771  on 59  degrees of freedom
# AIC: 86.771
# 
# Number of Fisher Scoring iterations: 4


##multivariate analysis
res.cox <- coxph(Surv(DFS, Status) ~ TIMES+`(CD3-CD57+)%`+`(CD3-CD57+)%-tc`, data = dat_stain_multi, id = `patient ID`)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Status) ~ TIMES + `(CD3-CD57+)%` + 
#           `(CD3-CD57+)%-tc`, data = dat_stain_multi, id = `patient ID`)
# 
# n= 61, number of events= 30 
# 
# 1:2                       coef  exp(coef)   se(coef)      z Pr(>|z|)
# TIMES              9.232e+00  1.022e+04  1.968e+00  4.692  2.7e-06
# `(CD3-CD57+)%`    -9.818e-02  9.065e-01  5.295e-01 -0.185    0.853
# `(CD3-CD57+)%-tc` -1.323e-02  9.869e-01  1.370e-01 -0.097    0.923
# 
# 1:2                 exp(coef) exp(-coef) lower .95 upper .95
# TIMES             1.022e+04  9.784e-05  216.1006 4.834e+05
# `(CD3-CD57+)%`    9.065e-01  1.103e+00    0.3211 2.559e+00
# `(CD3-CD57+)%-tc` 9.869e-01  1.013e+00    0.7544 1.291e+00
# 
# States: 1= (s0), 2= REC 
# 
# Concordance= 0.869  (se = 0.032 )
# Likelihood ratio test= 77.81  on 3 df,   p=<2e-16
# Wald test            = 22.15  on 3 df,   p=6e-05
# Score (logrank) test = 68.02  on 3 df,   p=1e-14


#3. "(CD3-CD16+CD56+)%"
##univariate analysis
##cd57 at IF
res.cox <- coxph(Surv(DFS, Status) ~ `(CD3-CD16+CD56+)%`, data = dat_stain_multi, id = `patient ID`)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Status) ~ `(CD3-CD16+CD56+)%`, data = dat_stain_multi, 
#         id = `patient ID`)
# 
# n= 61, number of events= 30 
# 
# 1:2                        coef exp(coef)  se(coef)      z Pr(>|z|)
# `(CD3-CD16+CD56+)%` -5.358721  0.004707  2.206067 -2.429   0.0151
# 
# 1:2                   exp(coef) exp(-coef) lower .95 upper .95
# `(CD3-CD16+CD56+)%`  0.004707      212.5 6.236e-05    0.3553
# 
# States: 1= (s0), 2= REC 
# 
# Concordance= 0.648  (se = 0.05 )
# Likelihood ratio test= 8.43  on 1 df,   p=0.004
# Wald test            = 5.9  on 1 df,   p=0.02
# Score (logrank) test = 6.27  on 1 df,   p=0.01

mylogit <- glm(Status ~ `(CD3-CD16+CD56+)%`, data = dat_stain_multi,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Status ~ `(CD3-CD16+CD56+)%`, family = binomial, 
#       data = dat_stain_multi)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.5520  -1.0738  -0.2236   0.9946   1.6240  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)           0.8865     0.4155   2.134  0.03286 * 
#   `(CD3-CD16+CD56+)%`  -7.5634     2.8887  -2.618  0.00884 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 84.548  on 60  degrees of freedom
# Residual deviance: 74.491  on 59  degrees of freedom
# AIC: 78.491
# 
# Number of Fisher Scoring iterations: 4


##cd57 at TC
res.cox <- coxph(Surv(DFS, Status) ~ `(CD3-CD16+CD56+)%-tc`, data = dat_stain_multi, id = `patient ID`)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Status) ~ `(CD3-CD16+CD56+)%-tc`, data = dat_stain_multi, 
#         id = `patient ID`)
# 
# n= 61, number of events= 30 
# 
# 1:2                         coef exp(coef) se(coef)      z Pr(>|z|)
# `(CD3-CD16+CD56+)%-tc` -1.1131    0.3286   1.3061 -0.852    0.394
# 
# 1:2                      exp(coef) exp(-coef) lower .95 upper .95
# `(CD3-CD16+CD56+)%-tc`    0.3286      3.044    0.0254      4.25
# 
# States: 1= (s0), 2= REC 
# 
# Concordance= 0.611  (se = 0.059 )
# Likelihood ratio test= 0.97  on 1 df,   p=0.3
# Wald test            = 0.73  on 1 df,   p=0.4
# Score (logrank) test = 0.74  on 1 df,   p=0.4

mylogit <- glm(Status ~ `(CD3-CD16+CD56+)%-tc`, data = dat_stain_multi,family=binomial)
summary(mylogit)
# Call:
#   glm(formula = Status ~ `(CD3-CD16+CD56+)%-tc`, family = binomial, 
#       data = dat_stain_multi)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.2428  -1.2168  -0.4607   1.1259   1.4333  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)              0.1649     0.3087   0.534    0.593
# `(CD3-CD16+CD56+)%-tc`  -1.7097     1.5944  -1.072    0.284
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 84.548  on 60  degrees of freedom
# Residual deviance: 83.051  on 59  degrees of freedom
# AIC: 87.051
# 
# Number of Fisher Scoring iterations: 3


##multivariate analysis
res.cox <- coxph(Surv(DFS, Status) ~ TIMES+`(CD3-CD16+CD56+)%`+`(CD3-CD16+CD56+)%-tc`, data = dat_stain_multi, id = `patient ID`)
summary(res.cox)
# Call:
#   coxph(formula = Surv(DFS, Status) ~ TIMES + `(CD3-CD16+CD56+)%` + 
#           `(CD3-CD16+CD56+)%-tc`, data = dat_stain_multi, id = `patient ID`)
# 
# n= 61, number of events= 30 
# 
# 1:2                            coef  exp(coef)   se(coef)      z Pr(>|z|)
# TIMES                   9.598e+00  1.473e+04  2.110e+00  4.548 5.41e-06
# `(CD3-CD16+CD56+)%`    -3.129e+00  4.378e-02  2.704e+00 -1.157    0.247
# `(CD3-CD16+CD56+)%-tc`  8.538e-02  1.089e+00  1.636e+00  0.052    0.958
# 
# 1:2                      exp(coef) exp(-coef) lower .95 upper .95
# TIMES                  1.473e+04  6.789e-05 2.355e+02 9.212e+05
# `(CD3-CD16+CD56+)%`    4.378e-02  2.284e+01 2.184e-04 8.774e+00
# `(CD3-CD16+CD56+)%-tc` 1.089e+00  9.182e-01 4.408e-02 2.691e+01
# 
# States: 1= (s0), 2= REC 
# 
# Concordance= 0.876  (se = 0.03 )
# Likelihood ratio test= 79.34  on 3 df,   p=<2e-16
# Wald test            = 21.11  on 3 df,   p=1e-04
# Score (logrank) test = 68.02  on 3 df,   p=1e-14


###compute ROC
ROC_times_stain <- roc(dat_stain_multi$Status, dat_stain_multi$TIMES, na.rm=T, smooth=T,algorithm=1, smooth.method="density")  #factor vs numeric
ROC_times_stain_auc <- auc(ROC_times_stain)
# Area under the curve: 0.9961
ROC_matureIF_stain <- roc(dat_stain_multi$Status, dat_stain_multi$`(CD3-CD57+)%`, na.rm=T, smooth=T,algorithm=1, smooth.method="density")  #factor vs numeric
ROC_matureIF_stain_auc <- auc(ROC_matureIF_stain)
# Area under the curve: 0.7449
ROC_cytotIF_stain <- roc(dat_stain_multi$Status, dat_stain_multi$`(CD3-CD16+CD56+)%`, na.rm=T, smooth=T,algorithm=1, smooth.method="density")  #factor vs numeric
ROC_cytotIF_stain_auc <- auc(ROC_cytotIF_stain)
# Area under the curve: 0.6821
ROC_ordinIF_stain <- roc(dat_stain_multi$Status, dat_stain_multi$`(CD3-CD56+)%`, na.rm=T, smooth=T,algorithm=1, smooth.method="density")  #factor vs numeric
ROC_ordinIF_stain_auc <- auc(ROC_ordinIF_stain)
# Area under the curve: 0.6141

ROC_matureTC_stain <- roc(dat_stain_multi$Status, dat_stain_multi$`(CD3-CD57+)%-tc`, na.rm=T, smooth=T,algorithm=1, smooth.method="density")  #factor vs numeric
ROC_matureTC_stain_auc <- auc(ROC_matureTC_stain)
# Area under the curve: 0.5409
ROC_cytotTC_stain <- roc(dat_stain_multi$Status, dat_stain_multi$`(CD3-CD16+CD56+)%-tc`, na.rm=T, smooth=T,algorithm=1, smooth.method="density")  #factor vs numeric
ROC_cytotTC_stain_auc <- auc(ROC_cytotTC_stain)
# Area under the curve: 0.5517
ROC_ordinTC_stain <- roc(dat_stain_multi$Status, dat_stain_multi$`(CD3-CD56+)%-tc`, na.rm=T, smooth=T,algorithm=1, smooth.method="density")  #factor vs numeric
ROC_ordinTC_stain_auc <- auc(ROC_ordinTC_stain)
# Area under the curve: 0.5798

plot(ROC_ordinTC_stain, col = "#559AC6", main = "",lwd=3, xlim=c(1,0), ylim=c(0,1))
lines(ROC_cytotTC_stain, col = "#E6AB02",lwd=3)
lines(ROC_matureTC_stain, col = "blue",lwd=3)
lines(ROC_ordinIF_stain, col = "black",lwd=3)
lines(ROC_cytotIF_stain, col = "#7570B3",lwd=3)
lines(ROC_matureIF_stain, col = "#66A61E",lwd=3)
lines(ROC_times_stain, col = "#E73334",lwd=3)

legend("bottomright", legend=c(paste("TIMES",round(ROC_times_stain_auc,2),sep=": "),paste("_matureIF",round(ROC_matureIF_stain_auc,2),sep=": "),
                               paste("_cytotIF",round(ROC_cytotIF_stain_auc,2),sep=": "),paste("_ordinIF",round(ROC_ordinIF_stain_auc,2),sep=": "),
                               paste("_matureTC",round(ROC_matureTC_stain_auc,2),sep=": "),paste("_cytotTC",round(ROC_cytotTC_stain_auc,2),sep=": "),paste("_ordinTC",round(ROC_ordinTC_stain_auc,2),sep=": ")), 
       col=c("#E73334","#66A61E","#7570B3","black","blue","#E6AB02","#559AC6"), bg=NA, lty=1, lwd=3.5, box.lty=0, cex=1)

#save as /data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/ExFig5f_ROCs_NK.pdf

