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

#######load data####
load("/data/Input_3_Comparison_between_clinical_markers_and_TIMES/stability_analysis_data_input.RData")

#######stability of results####
dim(times_data)
#1193   16

weighted.geomean <- function(x, w, ...) exp(weighted.mean(log(x), w, ...))

idx_id <- which(!is.na(times_data$mean_spon2))
##original
times_org <- c()
mean_org <- times_data[idx_id,c("mean_spon2","mean_zfp36l2","mean_zfp36","mean_vimentin","mean_hladr")]
sd_org <- times_data[idx_id,c("sd_spon2","sd_zfp36l2","sd_zfp36","sd_vimentin","sd_hladr")]
for (i in 1:dim(mean_org)[1]){
  weight_i <- 1/sd_org[i,]^2  
  times_i <- weighted.geomean(mean_org[i,], weight_i)
  times_org <- c(times_org, times_i)
}

idx_id_2 <- c(idx_id[-1]-1,dim(times_data)[1])

##2 samples
n=2
times_2sam <- c()
for (i in 1:length(idx_id)){
  idx_i <- idx_id[i]
  idx_i_2 <- idx_id_2[i]
  set.seed(i)
  select_idx <- sample(seq(idx_i,idx_i_2),n,replace = F)
  score_data <- times_data[select_idx,c("score_spon2","score_zfp36l2","score_zfp36","score_vimentin","score_hladr")]
  sd_i <- apply(score_data, 2, function(x) sqrt(var(x) * (n-1)/n))
  mean_i <- apply(score_data, 2, mean)
  
  weight_i <- 1/sd_i^2  
  times_i <- weighted.geomean(mean_i, weight_i)
  times_2sam <- c(times_2sam,times_i)
}

##3 samples
n=3
times_3sam <- c()
for (i in 1:length(idx_id)){
  idx_i <- idx_id[i]
  idx_i_2 <- idx_id_2[i]
  set.seed(i)
  select_idx <- sample(seq(idx_i,idx_i_2),n,replace = F)
  score_data <- times_data[select_idx,c("score_spon2","score_zfp36l2","score_zfp36","score_vimentin","score_hladr")]
  sd_i <- apply(score_data, 2, function(x) sqrt(var(x) * (n-1)/n))
  mean_i <- apply(score_data, 2, mean)
  
  weight_i <- 1/sd_i^2  
  times_i <- weighted.geomean(mean_i, weight_i)
  times_3sam <- c(times_3sam,times_i)
}


##correlation test
cor.test(times_org,times_2sam,method = "spearman")
cor.test(times_org,times_3sam,method = "spearman")

# Spearman's rank correlation rho
# 
# data:  times_org and times_2sam
# S = 365700, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.7729367 

# Spearman's rank correlation rho
# 
# data:  times_org and times_3sam
# S = 146922, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.9160705 


##plot##
##switch on the following comments if output is needed.
#pdf("/data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/ExFig6bc_sampleNo_correlation_HCC.pdf",width=9,height=4)
#par(mfrow = c(1,2),mar=c(4,4,2,2))

cor_test <- cor.test(times_org,times_2sam,method = "spearman")
df_1 <- data.frame(x=times_org,y=times_2sam,stringsAsFactors = F)
linearMod <- lm(y~x,df_1)
eq = paste0("y = ", round(linearMod$coefficients[2],2), "x + ", round(linearMod$coefficients[1],2))
plot(x=times_org,y=times_2sam, pch=20, col=adjustcolor("dodgerblue4", alpha.f=0.4),xlab="all samples are used for a patient",ylab="2 samples are used",cex.lab=1, cex.axis=1, bty="l",xlim=c(0,1), ylim=c(0,1.1), cex=1)
abline(linearMod, col=adjustcolor("red2", alpha.f=0.8),lty = 'dashed', lwd=3)
mtext(paste("spearman's rho = ",round(cor_test$estimate,4),"\np = ",formatC(cor_test$p.value, format = "e", digits = 2),sep=""), adj=1, cex=1, line=-2)
mtext(eq, adj=1, cex=1, side=1, line=-2)
summary(linearMod)
# Call:
#   lm(formula = y ~ x, data = df_1)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.73019 -0.09474 -0.01475  0.08152  0.71350 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.08497    0.02380   3.571 0.000441 ***
#   x            0.87291    0.04944  17.655  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1891 on 211 degrees of freedom
# (27 observations deleted due to missingness)
# Multiple R-squared:  0.5963,	Adjusted R-squared:  0.5944 
# F-statistic: 311.7 on 1 and 211 DF,  p-value: < 2.2e-16

cor_test <- cor.test(times_org,times_3sam,method = "spearman")
df_1 <- data.frame(x=times_org,y=times_3sam,stringsAsFactors = F)
linearMod <- lm(y~x,df_1)
eq = paste0("y = ", round(linearMod$coefficients[2],2), "x + ", round(linearMod$coefficients[1],2))
plot(x=times_org,y=times_3sam, pch=20, col=adjustcolor("dodgerblue4", alpha.f=0.4),xlab="all samples are used for a patient",ylab="3 samples are used",cex.lab=1, cex.axis=1, bty="l",xlim=c(0,1), ylim=c(0,1.1), cex=1)
abline(linearMod, col=adjustcolor("red2", alpha.f=0.8),lty = 'dashed', lwd=3)
mtext(paste("spearman's rho = ",round(cor_test$estimate,4),"\np = ",formatC(cor_test$p.value, format = "e", digits = 2),sep=""), adj=1, cex=1, line=-2)
mtext(eq, adj=1, cex=1, side=1, line=-2)
summary(linearMod)
# Call:
#   lm(formula = y ~ x, data = df_1)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.55879 -0.03088  0.00397  0.03233  0.47470 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.003019   0.013524   0.223    0.824    
# x           0.991011   0.028220  35.118   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1129 on 217 degrees of freedom
# (21 observations deleted due to missingness)
# Multiple R-squared:  0.8504,	Adjusted R-squared:  0.8497 
# F-statistic:  1233 on 1 and 217 DF,  p-value: < 2.2e-16

#dev.off()

###saved as /data/Output_3_Comparison_between_Clinical_Markers_and_TIMES/ExFig6bc_sampleNo_correlation_HCC.pdf : correlation for stability analysis

