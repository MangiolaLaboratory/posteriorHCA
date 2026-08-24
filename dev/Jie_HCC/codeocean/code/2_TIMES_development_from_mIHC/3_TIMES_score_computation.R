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
 #[1] caret_6.0-94         lattice_0.20-45      ggplot2_3.5.0        rpart.plot_3.1.2    
 #[5] RColorBrewer_1.1-3   lhs_1.1.6            rpart_4.1.16         gbm_2.1.9           
 #[9] randomForest_4.7-1.1 e1071_1.7-14         pROC_1.16.2          coop_0.6-3          
#[13] mltools_0.3.5        fmsb_0.7.6           readxl_1.4.3         readr_2.1.5         
#[17] fitdistrplus_1.1-11  survival_3.1-12      MASS_7.3-60.0.1     

#loaded via a namespace (and not attached):
 #[1] Rcpp_1.0.11          lubridate_1.9.3      listenv_0.9.1        class_7.3-20        
 #[5] digest_0.6.33        ipred_0.9-14         foreach_1.5.2        utf8_1.2.3          
 #[9] parallelly_1.37.1    R6_2.5.1             cellranger_1.1.0     plyr_1.8.9          
#[13] stats4_4.2.1         hardhat_1.3.1        pillar_1.9.0         rlang_1.1.1         
#[17] rstudioapi_0.15.0    data.table_1.15.2    Matrix_1.6-5         splines_4.2.1       
#[21] stringr_1.5.1        gower_1.0.1          munsell_0.5.0        proxy_0.4-27        
#[25] compiler_4.2.1       pkgconfig_2.0.3      globals_0.16.3       nnet_7.3-19         
#[29] tidyselect_1.2.0     tibble_3.2.1         prodlim_2023.08.28   codetools_0.2-18    
#[33] fansi_1.0.4          future_1.33.1        dplyr_1.1.2          tzdb_0.4.0          
#[37] withr_3.0.0          ModelMetrics_1.2.2.2 recipes_1.0.10       grid_4.2.1          
#[41] nlme_3.1-164         gtable_0.3.4         lifecycle_1.0.4      magrittr_2.0.3      
#[45] scales_1.3.0         stringi_1.7.12       future.apply_1.11.1  cli_3.6.1           
#[49] reshape2_1.4.4       timeDate_4032.109    generics_0.1.3       vctrs_0.6.3         
#[53] lava_1.8.0           iterators_1.0.14     tools_4.2.1          glue_1.6.2          
#[57] purrr_1.0.2          hms_1.1.3            parallel_4.2.1       timechange_0.3.0    
#[61] colorspace_2.1-0   



setwd("/code/2_TIMES_development_from_mIHC/")

library(fitdistrplus);library(readr);library(readxl);library(fmsb);library(mltools);library(coop);library(pROC);library(e1071);library(randomForest);library(gbm);library(rpart);library(lhs);library(RColorBrewer);library(rpart.plot);library(caret);library(smallstuff) 

##read in data
load("/data/Input_2_TIMES_development_from_mIHC/TIMES_SCORE_input.RData")

###mark the row where the sample is located
A4_2<-c(1:5);A16<-c(6:10);A33<-c(11:15);A34_1<-c(16:20);A36<-c(21:25);A31<-c(26:30);A3_2<-c(31:35);A39_2<-c(36:40);A35<-c(41:45)
A38<-c(46:50);A39_1<-c(51:55);A4_1<-c(56:60);A5<-c(61:65);A11<-c(66:70);A12<-c(71:75);A14_1<-c(76:79);A23<-c(80:84);A32<-c(85:89);A26<-c(90:94);A27<-c(95:98)
A30<-c(99:102);A34_2<-c(103:106);A40<-c(107:111);A59<-c(112:114);A61<-c(115:118);A62<-c(119:122);A63<-c(123:126);A64<-c(127:130);A7<-c(131:135);A8<-c(136:140);A17<-c(141:145)
A18_2<-c(146:150);A21<-c(151:155);A25<-c(156:160);A37_1<-c(161:165);A6<-c(166:170);A9_2<-c(171:175);A10_2<-c(176:180);A13<-c(181:185);A24<-c(186:190);A1<-c(191:195)
A2<-c(196:200);A15<-c(201:205);A19<-c(206:210);A20<-c(211:215);A22<-c(216:220);A37_2<-c(221:225);A41<-c(226:230);A42<-c(231:235);A43<-c(236:240);A44<-c(241:245);A46<-c(246:250)
A47<-c(251:255);A48<-c(256:260);A49<-c(261:265);A50<-c(266:270);A53<-c(271:275);A54<-c(276:278);A55<-c(279:283);A56<-c(284:288);A57<-c(289:293)

###according to the model training sequence,12 is SPON2,23 is ZFP36L2,34 is ZFP36,45 is VIM, and 56 is HLA-DRB1

weighted.geomean <- function(x, w, ...) exp(weighted.mean(log(x), w, ...))

#A4_2
myvariable_A4_2 <- c(mean(D_DATA[A4_2,12]),mean(D_DATA[A4_2,23]), mean(D_DATA[A4_2,34]),mean(D_DATA[A4_2,45]),
                     mean(D_DATA[A4_2,56]))##Some data to average  ## 0.8681768 0.7968680 0.7901432 0.8949941 0.7925568
myweight_A4_2 <- c(1/pop.sd(D_DATA[A4_2,12])^2,1/pop.sd(D_DATA[A4_2,23])^2,1/pop.sd(D_DATA[A4_2,34])^2,1/pop.sd(D_DATA[A4_2,45])^2, 1/pop.sd(D_DATA[A4_2,56])^2)   ##Weights  ## 93.25589  99.20930 185.11201 866.97516 147.99853
A4_2_score <- weighted.geomean(myvariable_A4_2, myweight_A4_2)  ## 0.8600793

#A16
myvariable_A16 <- c(mean(D_DATA[A16,12]),mean(D_DATA[A16,23]), mean(D_DATA[A16,34]),mean(D_DATA[A16,45]),
                    mean(D_DATA[A16,56]))##Some data to average  ## 0.8268878 0.8636770 0.8829983 0.7089639 0.7732907
myweight_A16 <- c(1/pop.sd(D_DATA[A16,12])^2,1/pop.sd(D_DATA[A16,23])^2,1/pop.sd(D_DATA[A16,34])^2,1/pop.sd(D_DATA[A16,45])^2,
                  1/pop.sd(D_DATA[A16,56])^2)   ##Weights  ## 51.09097 161.07069  50.21516  64.08914 419.46710
A16_score <- weighted.geomean(myvariable_A16, myweight_A16)  ## 0.7967842

#A33
myvariable_A33 <- c(mean(D_DATA[A33,12]),mean(D_DATA[A33,23]), mean(D_DATA[A33,34]),mean(D_DATA[A33,45]),
                    mean(D_DATA[A33,56]))##Some data to average  ## 0.8943496 0.7460463 0.8256020 0.6239864 0.7024668
myweight_A33 <- c(1/pop.sd(D_DATA[A33,12])^2,1/pop.sd(D_DATA[A33,23])^2,1/pop.sd(D_DATA[A33,34])^2,1/pop.sd(D_DATA[A33,45])^2,
                  1/pop.sd(D_DATA[A33,56])^2)   ##Weights  ## 271.09226  28.52018  80.55335  28.79599  96.56742
A33_score <- weighted.geomean(myvariable_A33, myweight_A33)  ## 0.8176928

#A34_1
myvariable_A34_1 <- c(mean(D_DATA[A34_1,12]),mean(D_DATA[A34_1,23]), mean(D_DATA[A34_1,34]),mean(D_DATA[A34_1,45]),
                      mean(D_DATA[A34_1,56]))##Some data to average  ## 0.8222004 0.8101577 0.7207990 0.8029562 0.7782979
myweight_A34_1 <- c(1/pop.sd(D_DATA[A34_1,12])^2,1/pop.sd(D_DATA[A34_1,23])^2,1/pop.sd(D_DATA[A34_1,34])^2,1/pop.sd(D_DATA[A34_1,45])^2,1/pop.sd(D_DATA[A34_1,56])^2)   ##Weights  ## 101.90452 148.11113  60.06387  39.76954 107.89602
A34_1_score <- weighted.geomean(myvariable_A34_1, myweight_A34_1)  ## 0.7923048

#A36
myvariable_A36 <- c(mean(D_DATA[A36,12]),mean(D_DATA[A36,23]), mean(D_DATA[A36,34]),mean(D_DATA[A36,45]),
                    mean(D_DATA[A36,56]))##Some data to average  ## 0.8355849 0.8446092 0.7703158 0.6964975 0.8178388
myweight_A36 <- c(1/pop.sd(D_DATA[A36,12])^2,1/pop.sd(D_DATA[A36,23])^2,1/pop.sd(D_DATA[A36,34])^2,1/pop.sd(D_DATA[A36,45])^2,
                  1/pop.sd(D_DATA[A36,56])^2)   ##Weights ## 126.51815 103.86272  43.87684  31.48162  89.01506
A36_score <- weighted.geomean(myvariable_A36, myweight_A36) ## 0.8144864

#A31
myvariable_A31 <- c(mean(D_DATA[A31,12]),mean(D_DATA[A31,23]), mean(D_DATA[A31,34]),mean(D_DATA[A31,45]),
                    mean(D_DATA[A31,56]))##Some data to average ## 0.8558510 0.8966402 0.7789311 0.7859570 0.7629398
myweight_A31 <- c(1/pop.sd(D_DATA[A31,12])^2,1/pop.sd(D_DATA[A31,23])^2,1/pop.sd(D_DATA[A31,34])^2,1/pop.sd(D_DATA[A31,45])^2,
                  1/pop.sd(D_DATA[A31,56])^2)   ##Weights ## 108.79097 201.23947  55.56764 558.27664 233.39390
A31_score <- weighted.geomean(myvariable_A31, myweight_A31) ## 0.8054459

#A3_2
myvariable_A3_2 <- c(mean(D_DATA[A3_2,12]),mean(D_DATA[A3_2,23]), mean(D_DATA[A3_2,34]),mean(D_DATA[A3_2,45]),
                     mean(D_DATA[A3_2,56]))##Some data to average  ## 0.9175715 0.7898069 0.6883359 0.9498124 0.8243918
myweight_A3_2 <- c(1/pop.sd(D_DATA[A3_2,12])^2,1/pop.sd(D_DATA[A3_2,23])^2,1/pop.sd(D_DATA[A3_2,34])^2,1/pop.sd(D_DATA[A3_2,45])^2,1/pop.sd(D_DATA[A3_2,56])^2)   ##Weights  ##   642.53581    70.72637    68.91075 50337.96025   201.64791
A3_2_score <- weighted.geomean(myvariable_A3_2, myweight_A3_2) ## 0.9482224

#A39_2
myvariable_A39_2 <- c(mean(D_DATA[A39_2,12]),mean(D_DATA[A39_2,23]), mean(D_DATA[A39_2,34]),mean(D_DATA[A39_2,45]),
                      mean(D_DATA[A39_2,56]))##Some data to average  ## 0.8887775 0.8019247 0.8043101 0.8666472 0.8955388
myweight_A39_2 <- c(1/pop.sd(D_DATA[A39_2,12])^2,1/pop.sd(D_DATA[A39_2,23])^2,1/pop.sd(D_DATA[A39_2,34])^2,1/pop.sd(D_DATA[A39_2,45])^2,1/pop.sd(D_DATA[A39_2,56])^2)   ##Weights  ## 108.7575  109.6507  138.8190  151.1437 1924.7339
A39_2_score <- weighted.geomean(myvariable_A39_2, myweight_A39_2)  ## 0.8835465

#A35
myvariable_A35 <- c(mean(D_DATA[A35,12]),mean(D_DATA[A35,23]), mean(D_DATA[A35,34]),mean(D_DATA[A35,45]),
                    mean(D_DATA[A35,56]))##Some data to average  ## 0.9030043 0.8546803 0.7167555 0.7924343 0.8035799
myweight_A35 <- c(1/pop.sd(D_DATA[A35,12])^2,1/pop.sd(D_DATA[A35,23])^2,1/pop.sd(D_DATA[A35,34])^2,1/pop.sd(D_DATA[A35,45])^2,
                  1/pop.sd(D_DATA[A35,56])^2)   ##Weights  ## 1328.87111  213.26679   25.73492  103.65737  402.00851
A35_score <- weighted.geomean(myvariable_A35, myweight_A35) ## 0.8696217

#A38
myvariable_A38 <- c(mean(D_DATA[A38,12]),mean(D_DATA[A38,23]), mean(D_DATA[A38,34]),mean(D_DATA[A38,45]),
                    mean(D_DATA[A38,56]))##Some data to average  ## 0.7952590 0.7568380 0.7958425 0.8111118 0.7262411
myweight_A38 <- c(1/pop.sd(D_DATA[A38,12])^2,1/pop.sd(D_DATA[A38,23])^2,1/pop.sd(D_DATA[A38,34])^2,1/pop.sd(D_DATA[A38,45])^2,
                  1/pop.sd(D_DATA[A38,56])^2)   ##Weights  ## 250.49888  56.26287 210.49630 243.58056 238.15693
A38_score <- weighted.geomean(myvariable_A38, myweight_A38)  ## 0.7799292

#A39_1
myvariable_A39_1 <- c(mean(D_DATA[A39_1,12]),mean(D_DATA[A39_1,23]), mean(D_DATA[A39_1,34]),mean(D_DATA[A39_1,45]),
                      mean(D_DATA[A39_1,56]))##Some data to average  ## 0.8651548 0.7853747 0.7957800 0.8060260 0.8619805
myweight_A39_1 <- c(1/pop.sd(D_DATA[A39_1,12])^2,1/pop.sd(D_DATA[A39_1,23])^2,1/pop.sd(D_DATA[A39_1,34])^2,1/pop.sd(D_DATA[A39_1,45])^2, 1/pop.sd(D_DATA[A39_1,56])^2)   ##Weights  ## 190.05040  54.82881  83.57095  28.49997 199.11035
A39_1_score <- weighted.geomean(myvariable_A39_1, myweight_A39_1)  ## 0.8420692

#A4_1
myvariable_A4_1 <- c(mean(D_DATA[A4_1,12]),mean(D_DATA[A4_1,23]), mean(D_DATA[A4_1,34]),mean(D_DATA[A4_1,45]),
                     mean(D_DATA[A4_1,56]))##Some data to average  ## 0.8634495 0.8741069 0.6649354 0.7969981 0.8955932
myweight_A4_1 <- c(1/pop.sd(D_DATA[A4_1,12])^2,1/pop.sd(D_DATA[A4_1,23])^2,1/pop.sd(D_DATA[A4_1,34])^2,1/pop.sd(D_DATA[A4_1,45])^2,1/pop.sd(D_DATA[A4_1,56])^2)   ##Weights  ## 925.74986 186.91610  30.81964  43.06486 375.03157
A4_1_score <- weighted.geomean(myvariable_A4_1, myweight_A4_1) ## 0.8659415

#A5
myvariable_A5 <- c(mean(D_DATA[A5,12]),mean(D_DATA[A5,23]), mean(D_DATA[A5,34]),mean(D_DATA[A5,45]),
                   mean(D_DATA[A5,56]))##Some data to average  ## 0.8653078 0.7733259 0.7995267 0.8292710 0.8930929
myweight_A5 <- c(1/pop.sd(D_DATA[A5,12])^2,1/pop.sd(D_DATA[A5,23])^2,1/pop.sd(D_DATA[A5,34])^2,1/pop.sd(D_DATA[A5,45])^2,
                 1/pop.sd(D_DATA[A5,56])^2)   ##Weights  ## 198.62913  40.59706  55.39164 187.13742 461.25105
A5_score <- weighted.geomean(myvariable_A5, myweight_A5)  ## 0.8631774

#A11
myvariable_A11 <- c(mean(D_DATA[A11,12]),mean(D_DATA[A11,23]), mean(D_DATA[A11,34]),mean(D_DATA[A11,45]),
                    mean(D_DATA[A11,56]))##Some data to average  ## 0.7604015 0.8007394 0.8404768 0.8565877 0.8688481
myweight_A11 <- c(1/pop.sd(D_DATA[A11,12])^2,1/pop.sd(D_DATA[A11,23])^2,1/pop.sd(D_DATA[A11,34])^2,1/pop.sd(D_DATA[A11,45])^2,
                  1/pop.sd(D_DATA[A11,56])^2)   ##Weights  ## 33.04966  88.22900 618.35908 117.83524 198.41748
A11_score <- weighted.geomean(myvariable_A11, myweight_A11)  ## 0.8414663

#A12
myvariable_A12 <- c(mean(D_DATA[A12,12]),mean(D_DATA[A12,23]), mean(D_DATA[A12,34]),mean(D_DATA[A12,45]),
                    mean(D_DATA[A12,56]))##Some data to average  ## 0.8119736 0.7706018 0.7972130 0.7471124 0.8599676
myweight_A12 <- c(1/pop.sd(D_DATA[A12,12])^2,1/pop.sd(D_DATA[A12,23])^2,1/pop.sd(D_DATA[A12,34])^2,1/pop.sd(D_DATA[A12,45])^2,
                  1/pop.sd(D_DATA[A12,56])^2)   ##Weights  ## 99.42561 200.50200 585.01524  40.50842  98.43362
A12_score <- weighted.geomean(myvariable_A12, myweight_A12)  ## 0.7970934

#A14_1
myvariable_A14_1 <- c(mean(D_DATA[A14_1,12]),mean(D_DATA[A14_1,23]), mean(D_DATA[A14_1,34]),mean(D_DATA[A14_1,45]),
                      mean(D_DATA[A14_1,56]))##Some data to average  ## 0.8840616 0.7213683 0.8302680 0.7403093 0.8733166
myweight_A14_1 <- c(1/pop.sd(D_DATA[A14_1,12])^2,1/pop.sd(D_DATA[A14_1,23])^2,1/pop.sd(D_DATA[A14_1,34])^2,1/pop.sd(D_DATA[A14_1,45])^2,1/pop.sd(D_DATA[A14_1,56])^2)   ##Weights  ## 245.42854  44.10401 523.19577  39.73137 232.62700
A14_1_score <- weighted.geomean(myvariable_A14_1, myweight_A14_1)  ## 0.8429193

#A23
myvariable_A23 <- c(mean(D_DATA[A23,12]),mean(D_DATA[A23,23]), mean(D_DATA[A23,34]),mean(D_DATA[A23,45]),
                    mean(D_DATA[A23,56]))##Some data to average  ## 0.8032287 0.7723234 0.7655707 0.8650906 0.8293853
myweight_A23 <- c(1/pop.sd(D_DATA[A23,12])^2,1/pop.sd(D_DATA[A23,23])^2,1/pop.sd(D_DATA[A23,34])^2,1/pop.sd(D_DATA[A23,45])^2,
                  1/pop.sd(D_DATA[A23,56])^2)   ##Weights  ## 108.78128  86.34288  92.31813 271.65808 186.27502
A23_score <- weighted.geomean(myvariable_A23, myweight_A23)  ## 0.8231978

#A32
myvariable_A32 <- c(mean(D_DATA[A32,12]),mean(D_DATA[A32,23]), mean(D_DATA[A32,34]),mean(D_DATA[A32,45]),
                    mean(D_DATA[A32,56]))##Some data to average  ## 0.9088981 0.7732209 0.7642211 0.9412641 0.8147877
myweight_A32 <- c(1/pop.sd(D_DATA[A32,12])^2,1/pop.sd(D_DATA[A32,23])^2,1/pop.sd(D_DATA[A32,34])^2,1/pop.sd(D_DATA[A32,45])^2,
                  1/pop.sd(D_DATA[A32,56])^2)   ##Weights  ## 268.40985    38.58266   179.21669 10657.43360   407.40904
A32_score <- weighted.geomean(myvariable_A32, myweight_A32)  ## 0.9320921

#A26
myvariable_A26 <- c(mean(D_DATA[A26,12]),mean(D_DATA[A26,23]), mean(D_DATA[A26,34]),mean(D_DATA[A26,45]),
                    mean(D_DATA[A26,56]))##Some data to average  ## 0.6602089 0.7190649 0.7533375 0.6166749 0.7991225
myweight_A26 <- c(1/pop.sd(D_DATA[A26,12])^2,1/pop.sd(D_DATA[A26,23])^2,1/pop.sd(D_DATA[A26,34])^2,1/pop.sd(D_DATA[A26,45])^2,
                  1/pop.sd(D_DATA[A26,56])^2)   ##Weights  ## 78.05789  22.39576 120.63008  98.25148 206.96504
A26_score <- weighted.geomean(myvariable_A26, myweight_A26)  ## 0.7269037

#A27
myvariable_A27 <- c(mean(D_DATA[A27,12]),mean(D_DATA[A27,23]), mean(D_DATA[A27,34]),mean(D_DATA[A27,45]),
                    mean(D_DATA[A27,56]))##Some data to average  ## 0.8225011 0.7476877 0.7917618 0.8163084 0.9128848
myweight_A27 <- c(1/pop.sd(D_DATA[A27,12])^2,1/pop.sd(D_DATA[A27,23])^2,1/pop.sd(D_DATA[A27,34])^2,1/pop.sd(D_DATA[A27,45])^2,
                  1/pop.sd(D_DATA[A27,56])^2)   ##Weights  ## 77.07806  52.76121  55.72327 176.70653 997.88276
A27_score <- weighted.geomean(myvariable_A27, myweight_A27)  ## 0.8823591

#A30
myvariable_A30 <- c(mean(D_DATA[A30,12]),mean(D_DATA[A30,23]), mean(D_DATA[A30,34]),mean(D_DATA[A30,45]),
                    mean(D_DATA[A30,56]))##Some data to average  ## 0.5296273 0.7164445 0.6189021 0.8589963 0.6217063
myweight_A30 <- c(1/pop.sd(D_DATA[A30,12])^2,1/pop.sd(D_DATA[A30,23])^2,1/pop.sd(D_DATA[A30,34])^2,1/pop.sd(D_DATA[A30,45])^2,
                  1/pop.sd(D_DATA[A30,56])^2)   ##Weights  ## 103.3401  118.4975   78.4189 1957.0972  168.2283
A30_score <- weighted.geomean(myvariable_A30, myweight_A30)  ## 0.8069612

#A34_2
myvariable_A34_2 <- c(mean(D_DATA[A34_2,12]),mean(D_DATA[A34_2,23]), mean(D_DATA[A34_2,34]),mean(D_DATA[A34_2,45]),
                      mean(D_DATA[A34_2,56]))##Some data to average  ## 0.6188620 0.7190700 0.7490197 0.7483039 0.8027732
myweight_A34_2 <- c(1/pop.sd(D_DATA[A34_2,12])^2,1/pop.sd(D_DATA[A34_2,23])^2,1/pop.sd(D_DATA[A34_2,34])^2,1/pop.sd(D_DATA[A34_2,45])^2,1/pop.sd(D_DATA[A34_2,56])^2)   ##Weights  ## 46.60305  38.85205 190.90252  70.17503 305.83780
A34_2_score <- weighted.geomean(myvariable_A34_2, myweight_A34_2)  ## 0.761346

#A40
myvariable_A40 <- c(mean(D_DATA[A40,12]),mean(D_DATA[A40,23]), mean(D_DATA[A40,34]),mean(D_DATA[A40,45]),
                    mean(D_DATA[A40,56]))##Some data to average  ## 0.7448856 0.7938306 0.7272111 0.8647772 0.8282013
myweight_A40 <- c(1/pop.sd(D_DATA[A40,12])^2,1/pop.sd(D_DATA[A40,23])^2,1/pop.sd(D_DATA[A40,34])^2,1/pop.sd(D_DATA[A40,45])^2,
                  1/pop.sd(D_DATA[A40,56])^2)   ##Weights  ## 40.91214 294.00783  26.57555 143.31946 522.46950
A40_score <- weighted.geomean(myvariable_A40, myweight_A40)  ## 0.8169419

#A59
myvariable_A59 <- c(mean(D_DATA[A59,12]),mean(D_DATA[A59,23]), mean(D_DATA[A59,34]),mean(D_DATA[A59,45]),
                    mean(D_DATA[A59,56]))##Some data to average  ## 0.8221278 0.6400914 0.7922060 0.8011191 0.8947018
myweight_A59 <- c(1/pop.sd(D_DATA[A59,12])^2,1/pop.sd(D_DATA[A59,23])^2,1/pop.sd(D_DATA[A59,34])^2,1/pop.sd(D_DATA[A59,45])^2,
                  1/pop.sd(D_DATA[A59,56])^2)   ##Weights  ## 168.10953   10.21676 5248.67739 1061.55056  273.18590
A59_score <- weighted.geomean(myvariable_A59, myweight_A59)  ## 0.7979877

#A61
myvariable_A61 <- c(mean(D_DATA[A61,12]),mean(D_DATA[A61,23]), mean(D_DATA[A61,34]),mean(D_DATA[A61,45]),
                    mean(D_DATA[A61,56]))##Some data to average  ## 0.7623067 0.7153039 0.6765456 0.6472364 0.7421734
myweight_A61 <- c(1/pop.sd(D_DATA[A61,12])^2,1/pop.sd(D_DATA[A61,23])^2,1/pop.sd(D_DATA[A61,34])^2,1/pop.sd(D_DATA[A61,45])^2,
                  1/pop.sd(D_DATA[A61,56])^2)   ##Weights  ## 35.53857 438.93050 122.63868  38.17133  78.09678
A61_score <- weighted.geomean(myvariable_A61, myweight_A61)  ## 0.709803

#A62
myvariable_A62 <- c(mean(D_DATA[A62,12]),mean(D_DATA[A62,23]), mean(D_DATA[A62,34]),mean(D_DATA[A62,45]),
                    mean(D_DATA[A62,56]))##Some data to average  ## 0.8092858 0.7847566 0.7330935 0.7209147 0.7535627
myweight_A62 <- c(1/pop.sd(D_DATA[A62,12])^2,1/pop.sd(D_DATA[A62,23])^2,1/pop.sd(D_DATA[A62,34])^2,1/pop.sd(D_DATA[A62,45])^2,
                  1/pop.sd(D_DATA[A62,56])^2)   ##Weights  ## 162.43476 199.54608 293.27065 403.59014  54.60996
A62_score <- weighted.geomean(myvariable_A62, myweight_A62)  ## 0.7493312

#A63
myvariable_A63 <- c(mean(D_DATA[A63,12]),mean(D_DATA[A63,23]), mean(D_DATA[A63,34]),mean(D_DATA[A63,45]),
                    mean(D_DATA[A63,56]))##Some data to average  # 0.6180983 0.6850619 0.7829468 0.7696193 0.8323585
myweight_A63 <- c(1/pop.sd(D_DATA[A63,12])^2,1/pop.sd(D_DATA[A63,23])^2,1/pop.sd(D_DATA[A63,34])^2,1/pop.sd(D_DATA[A63,45])^2,
                  1/pop.sd(D_DATA[A63,56])^2)   ##Weights  ## 95.53736 113.31870 149.55868 118.40639 276.08207
A63_score <- weighted.geomean(myvariable_A63, myweight_A63)  ## 0.7595338

#A64
myvariable_A64 <- c(mean(D_DATA[A64,12]),mean(D_DATA[A64,23]), mean(D_DATA[A64,34]),mean(D_DATA[A64,45]),
                    mean(D_DATA[A64,56]))##Some data to average  ## 0.8156098 0.7133044 0.7919957 0.7481662 0.6269729
myweight_A64 <- c(1/pop.sd(D_DATA[A64,12])^2,1/pop.sd(D_DATA[A64,23])^2,1/pop.sd(D_DATA[A64,34])^2,1/pop.sd(D_DATA[A64,45])^2,
                  1/pop.sd(D_DATA[A64,56])^2)   ##Weights  ## 600.05077   34.54602   77.86330 1344.54086   20.30005
A64_score <- weighted.geomean(myvariable_A64, myweight_A64)  ## 0.7667574

#A7
myvariable_A7 <- c(mean(D_DATA[A7,12]),mean(D_DATA[A7,23]), mean(D_DATA[A7,34]),mean(D_DATA[A7,45]),
                   mean(D_DATA[A7,56]))##Some data to average  ## 0.8138257 0.8061658 0.7569045 0.7325619 0.7633649
myweight_A7 <- c(1/pop.sd(D_DATA[A7,12])^2,1/pop.sd(D_DATA[A7,23])^2,1/pop.sd(D_DATA[A7,34])^2,1/pop.sd(D_DATA[A7,45])^2,
                 1/pop.sd(D_DATA[A7,56])^2)   ##Weights  ## 120.70633  45.73522 125.65919  53.88522  22.19100
A7_score <- weighted.geomean(myvariable_A7, myweight_A7)  ## 0.7778778

#A8
myvariable_A8 <- c(mean(D_DATA[A8,12]),mean(D_DATA[A8,23]), mean(D_DATA[A8,34]),mean(D_DATA[A8,45]),
                   mean(D_DATA[A8,56]))##Some data to average  ## 0.8864592 0.8959826 0.8755712 0.8386150 0.765094
myweight_A8 <- c(1/pop.sd(D_DATA[A8,12])^2,1/pop.sd(D_DATA[A8,23])^2,1/pop.sd(D_DATA[A8,34])^2,1/pop.sd(D_DATA[A8,45])^2,
                 1/pop.sd(D_DATA[A8,56])^2)   ##Weights  ## 110.81440 424.10406 315.95727 217.57407  47.62144
A8_score <- weighted.geomean(myvariable_A8, myweight_A8)  ## 0.8719206

#A17
myvariable_A17 <- c(mean(D_DATA[A17,12]),mean(D_DATA[A17,23]), mean(D_DATA[A17,34]),mean(D_DATA[A17,45]),
                    mean(D_DATA[A17,56]))##Some data to average  ## 0.8676519 0.6582666 0.7769788 0.8329133 0.8943050
myweight_A17 <- c(1/pop.sd(D_DATA[A17,12])^2,1/pop.sd(D_DATA[A17,23])^2,1/pop.sd(D_DATA[A17,34])^2,1/pop.sd(D_DATA[A17,45])^2,
                  1/pop.sd(D_DATA[A17,56])^2)   ##Weights  ## 76.61060  62.03987 569.18087 272.31448 120.37211
A17_score <- weighted.geomean(myvariable_A17, myweight_A17)  ## 0.8013786

#A18_2
myvariable_A18_2 <- c(mean(D_DATA[A18_2,12]),mean(D_DATA[A18_2,23]), mean(D_DATA[A18_2,34]),mean(D_DATA[A18_2,45]),
                      mean(D_DATA[A18_2,56]))##Some data to average  ## 0.2546821 0.2536929 0.2486658 0.2199824 0.1919727
myweight_A18_2 <- c(1/pop.sd(D_DATA[A18_2,12])^2,1/pop.sd(D_DATA[A18_2,23])^2,1/pop.sd(D_DATA[A18_2,34])^2,1/pop.sd(D_DATA[A18_2,45])^2,1/pop.sd(D_DATA[A18_2,56])^2)   ##Weights  ## 316.85990  87.74356 164.60431  41.88833 152.38142
A18_2_score <- weighted.geomean(myvariable_A18_2, myweight_A18_2)  ## 0.2374513

#A21
myvariable_A21 <- c(mean(D_DATA[A21,12]),mean(D_DATA[A21,23]), mean(D_DATA[A21,34]),mean(D_DATA[A21,45]),
                    mean(D_DATA[A21,56]))##Some data to average  ## 0.3447483 0.3005808 0.2188586 0.1481435 0.1567003
myweight_A21 <- c(1/pop.sd(D_DATA[A21,12])^2,1/pop.sd(D_DATA[A21,23])^2,1/pop.sd(D_DATA[A21,34])^2,1/pop.sd(D_DATA[A21,45])^2,
                  1/pop.sd(D_DATA[A21,56])^2)   ##Weights  ## 113.45014  60.54381 104.73596 567.15621 158.93763
A21_score <- weighted.geomean(myvariable_A21, myweight_A21)  ## 0.1787047

#A25
myvariable_A25 <- c(mean(D_DATA[A25,12]),mean(D_DATA[A25,23]), mean(D_DATA[A25,34]),mean(D_DATA[A25,45]),
                    mean(D_DATA[A25,56]))##Some data to average  ## 0.3217887 0.2079028 0.3103731 0.2323987 0.2703097
myweight_A25 <- c(1/pop.sd(D_DATA[A25,12])^2,1/pop.sd(D_DATA[A25,23])^2,1/pop.sd(D_DATA[A25,34])^2,1/pop.sd(D_DATA[A25,45])^2,
                  1/pop.sd(D_DATA[A25,56])^2)   ##Weights  ## 273.60398 185.62943  18.02207  70.46502 312.34233
A25_score <- weighted.geomean(myvariable_A25, myweight_A25)  ## 0.2674373

#A37_1
myvariable_A37_1 <- c(mean(D_DATA[A37_1,12]),mean(D_DATA[A37_1,23]), mean(D_DATA[A37_1,34]),mean(D_DATA[A37_1,45]),
                      mean(D_DATA[A37_1,56]))##Some data to average  ## 0.2399122 0.2395951 0.2586685 0.2234571 0.1901271
myweight_A37_1 <- c(1/pop.sd(D_DATA[A37_1,12])^2,1/pop.sd(D_DATA[A37_1,23])^2,1/pop.sd(D_DATA[A37_1,34])^2,1/pop.sd(D_DATA[A37_1,45])^2,1/pop.sd(D_DATA[A37_1,56])^2)   ##Weights  ## 84.68808 338.70335  48.88268 137.55970  63.94557
A37_1_score <- weighted.geomean(myvariable_A37_1, myweight_A37_1)  ## 0.2324069

#A6
myvariable_A6 <- c(mean(D_DATA[A6,12]),mean(D_DATA[A6,23]), mean(D_DATA[A6,34]),mean(D_DATA[A6,45]),
                   mean(D_DATA[A6,56]))##Some data to average  ## 0.1675562 0.1647110 0.2421950 0.1612614 0.1990028
myweight_A6 <- c(1/pop.sd(D_DATA[A6,12])^2,1/pop.sd(D_DATA[A6,23])^2,1/pop.sd(D_DATA[A6,34])^2,1/pop.sd(D_DATA[A6,45])^2,
                 1/pop.sd(D_DATA[A6,56])^2)   ##Weights  ## 90.92686  73.60324 197.20357  30.11867  41.41548
A6_score <- weighted.geomean(myvariable_A6, myweight_A6)  ## 0.2003125

#A9_2
myvariable_A9_2 <- c(mean(D_DATA[A9_2,12]),mean(D_DATA[A9_2,23]), mean(D_DATA[A9_2,34]),mean(D_DATA[A9_2,45]),
                     mean(D_DATA[A9_2,56]))##Some data to average  ## 0.2254188 0.1755901 0.1842855 0.1673071 0.2356756
myweight_A9_2 <- c(1/pop.sd(D_DATA[A9_2,12])^2,1/pop.sd(D_DATA[A9_2,23])^2,1/pop.sd(D_DATA[A9_2,34])^2,1/pop.sd(D_DATA[A9_2,45])^2,1/pop.sd(D_DATA[A9_2,56])^2)   ##Weights  ## 109.51124 212.18743 186.23785 109.45049  70.44189
A9_2_score <- weighted.geomean(myvariable_A9_2, myweight_A9_2)  ## 0.1893246

#A10_2
myvariable_A10_2 <- c(mean(D_DATA[A10_2,12]),mean(D_DATA[A10_2,23]), mean(D_DATA[A10_2,34]),mean(D_DATA[A10_2,45]),
                      mean(D_DATA[A10_2,56]))##Some data to average  ## 0.2789501 0.1903157 0.2768448 0.1118891 0.2363435
myweight_A10_2 <- c(1/pop.sd(D_DATA[A10_2,12])^2,1/pop.sd(D_DATA[A10_2,23])^2,1/pop.sd(D_DATA[A10_2,34])^2,1/pop.sd(D_DATA[A10_2,45])^2,1/pop.sd(D_DATA[A10_2,56])^2)   ##Weights  ## 79.85793 228.25560  45.72578 350.17885  26.30490
A10_2_score <- weighted.geomean(myvariable_A10_2, myweight_A10_2)  ## 0.1587083

#A13
myvariable_A13 <- c(mean(D_DATA[A13,12]),mean(D_DATA[A13,23]), mean(D_DATA[A13,34]),mean(D_DATA[A13,45]),
                    mean(D_DATA[A13,56]))##Some data to average  ## 0.2373651 0.1685988 0.1992689 0.1546391 0.2358687
myweight_A13 <- c(1/pop.sd(D_DATA[A13,12])^2,1/pop.sd(D_DATA[A13,23])^2,1/pop.sd(D_DATA[A13,34])^2,1/pop.sd(D_DATA[A13,45])^2,
                  1/pop.sd(D_DATA[A13,56])^2)   ##Weights  ## 155.77914 274.79615 145.67433 668.63035  51.87152
A13_score <- weighted.geomean(myvariable_A13, myweight_A13)  ## 0.1735139

#A24
myvariable_A24 <- c(mean(D_DATA[A24,12]),mean(D_DATA[A24,23]), mean(D_DATA[A24,34]),mean(D_DATA[A24,45]),
                    mean(D_DATA[A24,56]))##Some data to average  ## 0.2831552 0.1312572 0.2766321 0.2999274 0.1050931
myweight_A24 <- c(1/pop.sd(D_DATA[A24,12])^2,1/pop.sd(D_DATA[A24,23])^2,1/pop.sd(D_DATA[A24,34])^2,1/pop.sd(D_DATA[A24,45])^2,
                  1/pop.sd(D_DATA[A24,56])^2)   ##Weights  ## 55.12813 298.46813  37.93137  74.18400 156.61114
A24_score <- weighted.geomean(myvariable_A24, myweight_A24)  ## 0.1534338

#A1
myvariable_A1 <- c(mean(D_DATA[A1,12]),mean(D_DATA[A1,23]), mean(D_DATA[A1,34]),mean(D_DATA[A1,45]),
                   mean(D_DATA[A1,56]))##Some data to average  ## 0.2213982 0.2562292 0.2471515 0.2367353 0.1738497
myweight_A1 <- c(1/pop.sd(D_DATA[A1,12])^2,1/pop.sd(D_DATA[A1,23])^2,1/pop.sd(D_DATA[A1,34])^2,1/pop.sd(D_DATA[A1,45])^2,
                 1/pop.sd(D_DATA[A1,56])^2)   ##Weights  ## 124.08019  88.40406  46.24544  33.37113 294.42709
A1_score <- weighted.geomean(myvariable_A1, myweight_A1)  ## 0.2029799

#A2
myvariable_A2 <- c(mean(D_DATA[A2,12]),mean(D_DATA[A2,23]), mean(D_DATA[A2,34]),mean(D_DATA[A2,45]),
                   mean(D_DATA[A2,56]))##Some data to average  ## 0.2585977 0.2710504 0.3066784 0.2117230 0.1981043
myweight_A2 <- c(1/pop.sd(D_DATA[A2,12])^2,1/pop.sd(D_DATA[A2,23])^2,1/pop.sd(D_DATA[A2,34])^2,1/pop.sd(D_DATA[A2,45])^2,
                 1/pop.sd(D_DATA[A2,56])^2)   ##Weights  ## 300.72073  41.36757  19.81514  63.45320  59.86080
A2_score <- weighted.geomean(myvariable_A2, myweight_A2)  ## 0.2464645

#A15
myvariable_A15 <- c(mean(D_DATA[A15,12]),mean(D_DATA[A15,23]), mean(D_DATA[A15,34]),mean(D_DATA[A15,45]),
                    mean(D_DATA[A15,56]))##Some data to average  ## 0.1883103 0.1824289 0.1641991 0.2686014 0.1716222
myweight_A15 <- c(1/pop.sd(D_DATA[A15,12])^2,1/pop.sd(D_DATA[A15,23])^2,1/pop.sd(D_DATA[A15,34])^2,1/pop.sd(D_DATA[A15,45])^2,
                  1/pop.sd(D_DATA[A15,56])^2)   ##Weights  ## 42.79975 127.05211 130.89675  50.06624  59.17674
A15_score <- weighted.geomean(myvariable_A15, myweight_A15)  ## 0.1839178

#A19
myvariable_A19 <- c(mean(D_DATA[A19,12]),mean(D_DATA[A19,23]), mean(D_DATA[A19,34]),mean(D_DATA[A19,45]),
                    mean(D_DATA[A19,56]))##Some data to average  ## 0.2025356 0.1439984 0.2325813 0.2401293 0.2729259
myweight_A19 <- c(1/pop.sd(D_DATA[A19,12])^2,1/pop.sd(D_DATA[A19,23])^2,1/pop.sd(D_DATA[A19,34])^2,1/pop.sd(D_DATA[A19,45])^2,
                  1/pop.sd(D_DATA[A19,56])^2)   ##Weights  ## 106.74028  50.19648  43.87770  61.79958  70.52771
A19_score <- weighted.geomean(myvariable_A19, myweight_A19)  ## 0.2153941

#A20
myvariable_A20 <- c(mean(D_DATA[A20,12]),mean(D_DATA[A20,23]), mean(D_DATA[A20,34]),mean(D_DATA[A20,45]),
                    mean(D_DATA[A20,56]))##Some data to average  ## 0.06875204 0.08557863 0.15317728 0.19970806 0.30290846
myweight_A20 <- c(1/pop.sd(D_DATA[A20,12])^2,1/pop.sd(D_DATA[A20,23])^2,1/pop.sd(D_DATA[A20,34])^2,1/pop.sd(D_DATA[A20,45])^2,
                  1/pop.sd(D_DATA[A20,56])^2)   ##Weights  ## 1122.03142  160.13938  101.14939   79.88345  153.69553
A20_score <- weighted.geomean(myvariable_A20, myweight_A20)  ## 0.08965307

#A22
myvariable_A22 <- c(mean(D_DATA[A22,12]),mean(D_DATA[A22,23]), mean(D_DATA[A22,34]),mean(D_DATA[A22,45]),
                    mean(D_DATA[A22,56]))##Some data to average  ## 0.05931511 0.25418445 0.15254946 0.23939072 0.26796535
myweight_A22 <- c(1/pop.sd(D_DATA[A22,12])^2,1/pop.sd(D_DATA[A22,23])^2,1/pop.sd(D_DATA[A22,34])^2,1/pop.sd(D_DATA[A22,45])^2,
                  1/pop.sd(D_DATA[A22,56])^2)   ##Weights  ## 615.22076  58.17365 918.05530 281.56245  67.30834
A22_score <- weighted.geomean(myvariable_A22, myweight_A22)  ## 0.1249875

#A37_2
myvariable_A37_2 <- c(mean(D_DATA[A37_2,12]),mean(D_DATA[A37_2,23]), mean(D_DATA[A37_2,34]),mean(D_DATA[A37_2,45]),
                      mean(D_DATA[A37_2,56]))##Some data to average  ## 0.1475600 0.2057827 0.1705636 0.2159861 0.2700320
myweight_A37_2 <- c(1/pop.sd(D_DATA[A37_2,12])^2,1/pop.sd(D_DATA[A37_2,23])^2,1/pop.sd(D_DATA[A37_2,34])^2,1/pop.sd(D_DATA[A37_2,45])^2,1/pop.sd(D_DATA[A37_2,56])^2)   ##Weights  ## 104.92950 142.88866 338.32845  86.34762 134.61771
A37_2_score <- weighted.geomean(myvariable_A37_2, myweight_A37_2)  ## 0.1915969

#A41
myvariable_A41 <- c(mean(D_DATA[A41,12]),mean(D_DATA[A41,23]), mean(D_DATA[A41,34]),mean(D_DATA[A41,45]),
                    mean(D_DATA[A41,56]))##Some data to average  ## 0.07144199 0.10490452 0.16624875 0.19595165 0.11146058
myweight_A41 <- c(1/pop.sd(D_DATA[A41,12])^2,1/pop.sd(D_DATA[A41,23])^2,1/pop.sd(D_DATA[A41,34])^2,1/pop.sd(D_DATA[A41,45])^2,
                  1/pop.sd(D_DATA[A41,56])^2)   ##Weights  ## 360.04077 191.88992 224.90034  81.48027 500.34161
A41_score <- weighted.geomean(myvariable_A41, myweight_A41)  ## 0.1085553

#A42
myvariable_A42 <- c(mean(D_DATA[A42,12]),mean(D_DATA[A42,23]), mean(D_DATA[A42,34]),mean(D_DATA[A42,45]),
                    mean(D_DATA[A42,56]))##Some data to average  ## 0.1277249 0.1607335 0.1669829 0.2394749 0.2208910
myweight_A42 <- c(1/pop.sd(D_DATA[A42,12])^2,1/pop.sd(D_DATA[A42,23])^2,1/pop.sd(D_DATA[A42,34])^2,1/pop.sd(D_DATA[A42,45])^2,
                  1/pop.sd(D_DATA[A42,56])^2)   ##Weights  ## 215.60843 199.29011 127.59012 210.78282  47.10151
A42_score <- weighted.geomean(myvariable_A42, myweight_A42)  ## 0.1720207

#A43
myvariable_A43 <- c(mean(D_DATA[A43,12]),mean(D_DATA[A43,23]), mean(D_DATA[A43,34]),mean(D_DATA[A43,45]),
                    mean(D_DATA[A43,56]))##Some data to average  ## 0.09145923 0.21913124 0.24557669 0.20298743 0.09199722
myweight_A43 <- c(1/pop.sd(D_DATA[A43,12])^2,1/pop.sd(D_DATA[A43,23])^2,1/pop.sd(D_DATA[A43,34])^2,1/pop.sd(D_DATA[A43,45])^2,
                  1/pop.sd(D_DATA[A43,56])^2)   ##Weights  ## 511.0288 238.3876 241.1592 203.3707 210.4586
A43_score <- weighted.geomean(myvariable_A43, myweight_A43)  ## 0.1411956

#A44
myvariable_A44 <- c(mean(D_DATA[A44,12]),mean(D_DATA[A44,23]), mean(D_DATA[A44,34]),mean(D_DATA[A44,45]),
                    mean(D_DATA[A44,56]))##Some data to average  ## 0.0508713 0.1143693 0.1789831 0.1727069 0.1476455
myweight_A44 <- c(1/pop.sd(D_DATA[A44,12])^2,1/pop.sd(D_DATA[A44,23])^2,1/pop.sd(D_DATA[A44,34])^2,1/pop.sd(D_DATA[A44,45])^2,
                  1/pop.sd(D_DATA[A44,56])^2)   ##Weights  ## 1249.58310  214.52108  455.58873  165.53939   35.35426
A44_score <- weighted.geomean(myvariable_A44, myweight_A44)  ## 0.08102046F

#A46
myvariable_A46 <- c(mean(D_DATA[A46,12]),mean(D_DATA[A46,23]), mean(D_DATA[A46,34]),mean(D_DATA[A46,45]),
                    mean(D_DATA[A46,56]))##Some data to average  ## 0.0529067 0.2292708 0.1589254 0.2247381 0.1510699
myweight_A46 <- c(1/pop.sd(D_DATA[A46,12])^2,1/pop.sd(D_DATA[A46,23])^2,1/pop.sd(D_DATA[A46,34])^2,1/pop.sd(D_DATA[A46,45])^2,
                  1/pop.sd(D_DATA[A46,56])^2)   ##Weights  ## 1368.9978  344.7830  108.4540   78.1747  377.4486
A46_score <- weighted.geomean(myvariable_A46, myweight_A46)  ## 0.08703855

#A47
myvariable_A47 <- c(mean(D_DATA[A47,12]),mean(D_DATA[A47,23]), mean(D_DATA[A47,34]),mean(D_DATA[A47,45]),
                    mean(D_DATA[A47,56]))##Some data to average  ## 0.05195973 0.22167032 0.21249949 0.18560341 0.09098257
myweight_A47 <- c(1/pop.sd(D_DATA[A47,12])^2,1/pop.sd(D_DATA[A47,23])^2,1/pop.sd(D_DATA[A47,34])^2,1/pop.sd(D_DATA[A47,45])^2,
                  1/pop.sd(D_DATA[A47,56])^2)   ##Weights  ## 2168.17857   79.79030  111.81033   45.72423 1112.88640
A47_score <- weighted.geomean(myvariable_A47, myweight_A47)  ## 0.06816067

#A48
myvariable_A48 <- c(mean(D_DATA[A48,12]),mean(D_DATA[A48,23]), mean(D_DATA[A48,34]),mean(D_DATA[A48,45]),
                    mean(D_DATA[A48,56]))##Some data to average  ## 0.09026491 0.20141991 0.12631981 0.09821492 0.24280074
myweight_A48 <- c(1/pop.sd(D_DATA[A48,12])^2,1/pop.sd(D_DATA[A48,23])^2,1/pop.sd(D_DATA[A48,34])^2,1/pop.sd(D_DATA[A48,45])^2,
                  1/pop.sd(D_DATA[A48,56])^2)   ##Weights  ## 3127.1786   23.0742  334.8031 1918.7712  176.7863
A48_score <- weighted.geomean(myvariable_A48, myweight_A48)  ## 0.09815989

#A49
myvariable_A49 <- c(mean(D_DATA[A49,12]),mean(D_DATA[A49,23]), mean(D_DATA[A49,34]),mean(D_DATA[A49,45]),
                    mean(D_DATA[A49,56]))##Some data to average  ## 0.07143354 0.09756379 0.27427277 0.20343374 0.07695978
myweight_A49 <- c(1/pop.sd(D_DATA[A49,12])^2,1/pop.sd(D_DATA[A49,23])^2,1/pop.sd(D_DATA[A49,34])^2,1/pop.sd(D_DATA[A49,45])^2,
                  1/pop.sd(D_DATA[A49,56])^2)   ##Weights  ## 484.04658  201.74569   41.79254  101.58492 1169.07642
A49_score <- weighted.geomean(myvariable_A49, myweight_A49)  ## 0.08352773

#A50
myvariable_A50 <- c(mean(D_DATA[A50,12]),mean(D_DATA[A50,23]), mean(D_DATA[A50,34]),mean(D_DATA[A50,45]),
                    mean(D_DATA[A50,56]))##Some data to average  ## 0.13482024 0.25826412 0.16053095 0.11925756 0.04715907
myweight_A50 <- c(1/pop.sd(D_DATA[A50,12])^2,1/pop.sd(D_DATA[A50,23])^2,1/pop.sd(D_DATA[A50,34])^2,1/pop.sd(D_DATA[A50,45])^2,
                  1/pop.sd(D_DATA[A50,56])^2)   ##Weights  ## 98.31271   50.70869  303.23636  184.79851 2338.88550
A50_score <- weighted.geomean(myvariable_A50, myweight_A50)  ## 0.06031813

#A53
myvariable_A53 <- c(mean(D_DATA[A53,12]),mean(D_DATA[A53,23]), mean(D_DATA[A53,34]),mean(D_DATA[A53,45]),
                    mean(D_DATA[A53,56]))##Some data to average  ## 0.13970640 0.07507459 0.14719960 0.28730958 0.10267459
myweight_A53 <- c(1/pop.sd(D_DATA[A53,12])^2,1/pop.sd(D_DATA[A53,23])^2,1/pop.sd(D_DATA[A53,34])^2,1/pop.sd(D_DATA[A53,45])^2,
                  1/pop.sd(D_DATA[A53,56])^2)   ##Weights   ## 48.56415 883.55376 681.25325  26.44531 144.69544
A53_score <- weighted.geomean(myvariable_A53, myweight_A53)  ## 0.1033064

#A54
myvariable_A54 <- c(mean(D_DATA[A54,12]),mean(D_DATA[A54,23]), mean(D_DATA[A54,34]),mean(D_DATA[A54,45]),
                    mean(D_DATA[A54,56]))##Some data to average  ## 0.24651061 0.31656640 0.26466509 0.27470232 0.03778484
myweight_A54 <- c(1/pop.sd(D_DATA[A54,12])^2,1/pop.sd(D_DATA[A54,23])^2,1/pop.sd(D_DATA[A54,34])^2,1/pop.sd(D_DATA[A54,45])^2,
                  1/pop.sd(D_DATA[A54,56])^2)   ##Weights  ## 16.00247   26.84269 2246.58227   51.05905 6988.5127
A54_score <- weighted.geomean(myvariable_A54, myweight_A54)  ## 0.06161239

#A55
myvariable_A55 <- c(mean(D_DATA[A55,12]),mean(D_DATA[A55,23]), mean(D_DATA[A55,34]),mean(D_DATA[A55,45]),
                    mean(D_DATA[A55,56]))##Some data to average  ## 0.1295407 0.2233013 0.2510932 0.1280947 0.1079854
myweight_A55 <- c(1/pop.sd(D_DATA[A55,12])^2,1/pop.sd(D_DATA[A55,23])^2,1/pop.sd(D_DATA[A55,34])^2,1/pop.sd(D_DATA[A55,45])^2,
                  1/pop.sd(D_DATA[A55,56])^2)   ##Weights  ## 139.91129 197.52254  35.07134 301.30302 641.02793
A55_score <- weighted.geomean(myvariable_A55, myweight_A55)  ## 0.1306011

#A56
myvariable_A56 <- c(mean(D_DATA[A56,12]),mean(D_DATA[A56,23]), mean(D_DATA[A56,34]),mean(D_DATA[A56,45]),
                    mean(D_DATA[A56,56]))##Some data to average  ## 0.1466967 0.1593175 0.2235012 0.1631370 0.1123265
myweight_A56 <- c(1/pop.sd(D_DATA[A56,12])^2,1/pop.sd(D_DATA[A56,23])^2,1/pop.sd(D_DATA[A56,34])^2,1/pop.sd(D_DATA[A56,45])^2,
                  1/pop.sd(D_DATA[A56,56])^2)   ##Weights  ## 66.12919 220.74583  63.46608  99.11107 217.24686
A56_score <- weighted.geomean(myvariable_A56, myweight_A56)  ## 0.1461418

#A57
myvariable_A57 <- c(mean(D_DATA[A57,12]),mean(D_DATA[A57,23]), mean(D_DATA[A57,34]),mean(D_DATA[A57,45]),
                    mean(D_DATA[A57,56]))##Some data to average  ## 0.2904914 0.3055187 0.2436036 0.1236616 0.3057642
myweight_A57 <- c(1/pop.sd(D_DATA[A57,12])^2,1/pop.sd(D_DATA[A57,23])^2,1/pop.sd(D_DATA[A57,34])^2,1/pop.sd(D_DATA[A57,45])^2,
                  1/pop.sd(D_DATA[A57,56])^2)   ##Weights  ## 42.58238  32.93576 206.96028  99.42595  37.28392
A57_score <- weighted.geomean(myvariable_A57, myweight_A57)  ## 0.2193395


A_score<-data.frame(A4_2_score,A16_score,A33_score,A34_1_score,A36_score,A31_score,A3_2_score,A39_2_score,A35_score,A38_score,
                    A39_1_score,A4_1_score,A5_score,A11_score,A12_score,A14_1_score,A23_score,A32_score,A26_score,A27_score,
                    A30_score,A34_2_score,A40_score,A59_score,A61_score,A62_score,A63_score,A64_score,A7_score,A8_score,
                    A17_score,A18_2_score,A21_score,A25_score,A37_1_score,A6_score,A9_2_score,A10_2_score,A13_score,
                    A24_score,A1_score,A2_score,A15_score,A19_score,A20_score,A22_score,A37_2_score,A41_score,A42_score,
                    A43_score,A44_score,A46_score,A47_score,A48_score,A49_score,A50_score,A53_score,A54_score,A55_score,
                    A56_score,A57_score,stringsAsFactors = F)
#write.csv(A_score,file="LIVER_score_TIMES_seed1_LIVER_ALL.csv")
