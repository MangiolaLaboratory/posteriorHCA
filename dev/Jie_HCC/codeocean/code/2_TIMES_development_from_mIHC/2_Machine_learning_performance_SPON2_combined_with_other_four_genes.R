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
 #[1] rpart.plot_3.1.2     RColorBrewer_1.1-3   randomForest_4.7-1.1 caret_6.0-94         lattice_0.20-45     
 #[6] ggplot2_3.5.0        e1071_1.7-14         rpart_4.1.16         gbm_2.1.9            pROC_1.18.5         
#[11] mltools_0.3.5        lhs_1.1.6            readxl_1.4.3        

#loaded via a namespace (and not attached):
 #[1] Rcpp_1.0.11          lubridate_1.9.3      listenv_0.9.1        class_7.3-20         digest_0.6.33       
 #[6] ipred_0.9-14         foreach_1.5.2        utf8_1.2.3           parallelly_1.37.1    R6_2.5.1            
#[11] cellranger_1.1.0     plyr_1.8.9           hardhat_1.3.1        stats4_4.2.1         pillar_1.9.0        
#[16] rlang_1.1.1          rstudioapi_0.15.0    data.table_1.15.2    Matrix_1.6-5         splines_4.2.1       
#[21] gower_1.0.1          stringr_1.5.1        munsell_0.5.0        proxy_0.4-27         compiler_4.2.1      
#[26] pkgconfig_2.0.3      globals_0.16.3       nnet_7.3-19          tidyselect_1.2.0     tibble_3.2.1        
#[31] prodlim_2023.08.28   codetools_0.2-18     fansi_1.0.4          future_1.33.1        dplyr_1.1.2         
#[36] withr_3.0.0          MASS_7.3-60.0.1      recipes_1.0.10       ModelMetrics_1.2.2.2 grid_4.2.1          
#[41] nlme_3.1-164         gtable_0.3.4         lifecycle_1.0.4      magrittr_2.0.3       scales_1.3.0        
#[46] future.apply_1.11.1  cli_3.6.1            stringi_1.7.12       reshape2_1.4.4       timeDate_4032.109   
#[51] generics_0.1.3       vctrs_0.6.3          lava_1.8.0           iterators_1.0.14     tools_4.2.1         
#[56] glue_1.6.2           purrr_1.0.2          parallel_4.2.1       survival_3.5-8       timechange_0.3.0    
#[61] colorspace_2.1-0 

library(lhs);library(mltools);library(pROC);library(gbm);library(rpart);library(e1071);library(caret);library(randomForest);library(rpart);library(RColorBrewer);library(rpart.plot)

setwd("/code/2_TIMES_development_from_mIHC/")
##read in data
load("/data/Input_2_TIMES_development_from_mIHC/SPON2_combine_with_other_four_genes_mIHC.RData")


##global parameters
accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}
N_trial <- 100
coef <- qt(.975, df = N_trial-1)
ppv <- function(x){diag(x)[2]/rowSums(x)[2] * 100}
npv <- function(x){diag(x)[1]/rowSums(x)[1] * 100}

#training and testing vectors
acc_train_vec <- c();acc_train_sd <- c() 
ppv_train_vec <- c();npv_train_vec <- c()
ppv_train_vec <- c();npv_train_vec <- c()
auc_train_vec <- c();auc_train_vec <- c()

acc_test_vec <- c();acc_test_sd <- c() 
ppv_test_vec <- c();npv_test_vec <- c()
ppv_test_vec <- c();npv_test_vec <- c()
auc_test_vec <- c();auc_test_vec <- c()

acc_train_vec_SVM_R <- c();acc_train_sd_SVM_R <- c() 
ppv_train_vec_SVM_R <- c();npv_train_vec_SVM_R <- c()
ppv_train_vec_SVM_R <- c();npv_train_vec_SVM_R <- c()
auc_train_vec_SVM_R <- c();auc_train_vec_SVM_R <- c()

acc_test_vec_SVM_R <- c();acc_test_sd_SVM_R <- c() 
ppv_test_vec_SVM_R <- c();npv_test_vec_SVM_R <- c()
ppv_test_vec_SVM_R <- c();npv_test_vec_SVM_R <- c()
auc_test_vec_SVM_R <- c();auc_test_vec_SVM_R <- c()

acc_train_vec_SVM_L <- c();acc_train_sd_SVM_L <- c() 
ppv_train_vec_SVM_L <- c();npv_train_vec_SVM_L <- c()
ppv_train_vec_SVM_L <- c();npv_train_vec_SVM_L <- c()
auc_train_vec_SVM_L <- c();auc_train_vec_SVM_L <- c()

acc_test_vec_SVM_L <- c();acc_test_sd_SVM_L <- c() 
ppv_test_vec_SVM_L <- c();npv_test_vec_SVM_L <- c()
ppv_test_vec_SVM_L <- c();npv_test_vec_SVM_L <- c()
auc_test_vec_SVM_L <- c();auc_test_vec_SVM_L <- c()

acc_train_vec_RF <- c();acc_train_sd_RF <- c() 
ppv_train_vec_RF <- c();npv_train_vec_RF <- c()
ppv_train_vec_RF <- c();npv_train_vec_RF <- c()
auc_train_vec_RF <- c();auc_train_vec_RF <- c()

acc_test_vec_RF <- c();acc_test_sd_RF <- c() 
ppv_test_vec_RF <- c();npv_test_vec_RF <- c()
ppv_test_vec_RF <- c();npv_test_vec_RF <- c()
auc_test_vec_RF <- c();auc_test_vec_RF <- c()

acc_train_vec_DT <- c();acc_train_sd_DT <- c() 
ppv_train_vec_DT <- c();npv_train_vec_DT <- c()
ppv_train_vec_DT <- c();npv_train_vec_DT <- c()
auc_train_vec_DT <- c();auc_train_vec_DT <- c()

acc_test_vec_DT <- c();acc_test_sd_DT <- c() 
ppv_test_vec_DT <- c();npv_test_vec_DT <- c()
ppv_test_vec_DT <- c();npv_test_vec_DT <- c()
auc_test_vec_DT <- c();auc_test_vec_DT <- c()

acc_train_vec_bgm <- c();acc_train_sd_bgm <- c() 
ppv_train_vec_bgm <- c();npv_train_vec_bgm <- c()
ppv_train_vec_bgm <- c();npv_train_vec_bgm <- c()
auc_train_vec_bgm <- c();auc_train_vec_bgm <- c()

acc_test_vec_bgm <- c();acc_test_sd_bgm <- c() 
ppv_test_vec_bgm <- c();npv_test_vec_bgm <- c()
ppv_test_vec_bgm <- c();npv_test_vec_bgm <- c()
auc_test_vec_bgm <- c();auc_test_vec_bgm <- c()


for (trial_i in 1:100){
  
  set.seed(trial_i)
  train_ind_r_1 <- sample.int(n=length(liver_r_data_gene_ed_1[,1]), size = floor(0.75*length(liver_r_data_gene_ed_1[,1])), replace = F)
  train_df_r_1 <-liver_r_data_gene_ed_1[train_ind_r_1,]
  test_df_r_1 <-liver_r_data_gene_ed_1[-train_ind_r_1,]
  set.seed(trial_i)
  train_ind_nr_1 <- sample.int(n=length(liver_nr_data_gene_ed_1[,1]), size = floor(0.75*length(liver_nr_data_gene_ed_1[,1])), replace = F)
  train_df_nr_1 <-liver_nr_data_gene_ed_1[train_ind_nr_1,]
  test_df_nr_1 <-liver_nr_data_gene_ed_1[-train_ind_nr_1,]
  
  train_df_1 <- rbind(train_df_r_1,train_df_nr_1)
  train_df_1$Relp_y <- as.factor(train_df_1$Relp_y)
  test_df_1 <- rbind(test_df_r_1,test_df_nr_1)
  test_df_1$Relp_y <- as.factor(test_df_1$Relp_y)
  train_target_1 <- train_df_1$Relp_y
  test_target_1 <- test_df_1$Relp_y

  #combination of SPON2 and VIM
  linearMod <- glm("Relp_y ~ NKatBorder_x_SPON2 + NKatTumor_x_SPON2 + NKatBorder_x_VIM + NKatTumor_x_VIM",data = train_df_1,family=binomial,control=list(maxit=100)) #linear model
  model_SVM_R <- svm(Relp_y ~ NKatBorder_x_SPON2 + NKatTumor_x_SPON2 + NKatBorder_x_VIM + NKatTumor_x_VIM,data=train_df_1, kernel = "radial",probability = TRUE,decision.values = TRUE) #support vector machines - radial
  model_SVM_L <- svm(Relp_y ~ NKatBorder_x_SPON2 + NKatTumor_x_SPON2 + NKatBorder_x_VIM + NKatTumor_x_VIM,data=train_df_1, kernel = "linear",probability = TRUE,decision.values = TRUE) #support vector machines - linear
  model_RF <- randomForest(Relp_y ~ NKatBorder_x_SPON2 + NKatTumor_x_SPON2 + NKatBorder_x_VIM + NKatTumor_x_VIM,data = train_df_1,ntree =40,importance=TRUE,proximity=TRUE,replace=F) #random forest
  model_DT <- rpart("Relp_y ~ NKatBorder_x_SPON2 + NKatTumor_x_SPON2 + NKatBorder_x_VIM + NKatTumor_x_VIM",method="class",data = train_df_1,maxdepth=9) #decision tree
  gbm_tunepara <- data.frame(n.trees=1200,interaction.depth=11,shrinkage=0.01,n.minobsinnode=3)
  set.seed(trial_i)
  fn_x <- "as.factor(Relp_y) ~ NKatBorder_x_SPON2 + NKatTumor_x_SPON2 + NKatBorder_x_VIM + NKatTumor_x_VIM"
  xgbMod_x <- train(as.formula(fn_x),data = train_df_1,method = "gbm",trControl = trainControl(method="none",seeds = trial_i),verbose=0,tuneGrid=gbm_tunepara) #xgboost
  
  #training
  pred_prob <- predict(linearMod,train_df_1,type="response")
  pred_label <- ifelse(pred_prob > 0.5,1,0)
  pred_label <- as.factor(pred_label)
  tb <- table(pred_label, train_target_1) #confusion matrix, if pred_prob>0.5, we regarded it as 1, if pred_prob<0.5, we regarded it as 0.
  
  pred_prob_SVM_R <- predict(model_SVM_R, train_df_1,decision.values = TRUE, probability = TRUE)
  pred_prob_SVM_R_p <- attr(pred_prob_SVM_R, "probabilities")
  pred_prob_SVM_R_d <- attr(pred_prob_SVM_R, "decision.values")
  pred_label_SVM_R <- ifelse(pred_prob_SVM_R_p[,1] > 0.5,1,0)
  pred_label_SVM_R <- as.factor(pred_label_SVM_R)
  tb_SVM_R <- table(pred_label_SVM_R, train_target_1) #confusion matrix, if pred_label_SVM_R>0.5, we regarded it as 1, if pred_label_SVM_R<0.5, we regarded it as 0.
  
  pred_prob_SVM_L <- predict(model_SVM_L, train_df_1, decision.values = TRUE, probability = TRUE)
  pred_prob_SVM_L_p <- attr(pred_prob_SVM_L, "probabilities")
  pred_prob_SVM_L_d <- attr(pred_prob_SVM_L, "decision.values")
  pred_label_SVM_L <- ifelse(pred_prob_SVM_L_p[,1] > 0.5,1,0)
  pred_label_SVM_L <- as.factor(pred_label_SVM_L)
  tb_SVM_L <- table(pred_label_SVM_L, train_target_1) #confusion matrix, if pred_label_SVM_L>0.5, we regarded it as 1, if pred_label_SVM_L<0.5, we regarded it as 0.
  
  pred_prob_RF <- predict(model_RF,train_df_1,type="prob")
  pred_label_RF <- ifelse(pred_prob_RF[,2] > 0.5,1,0)
  pred_label_RF <- as.factor(pred_label_RF)
  tb_RF <- table(pred_label_RF, train_target_1) #confusion matrix, if pred_label_RF>0.5, we regarded it as 1, if pred_label_RF<0.5, we regarded it as 0.
  
  pred_prob_DT <- predict(model_DT,train_df_1,type="prob")
  pred_label_DT <- ifelse(pred_prob_DT[,2] > 0.5,1,0)
  pred_label_DT <- as.factor(pred_label_DT)
  tb_DT <- table(pred_label_DT, train_target_1) #confusion matrix, if pred_label_DT>0.5, we regarded it as 1, if pred_label_DT<0.5, we regarded it as 0.
  
  pred_prob_bgm <- predict(xgbMod_x,train_df_1,type="prob")
  pred_label_bgm <- ifelse(pred_prob_bgm[,2] > 0.5,1,0)
  pred_label_bgm <- as.factor(pred_label_bgm)
  tb_bgm <- table(pred_label_bgm, train_target_1) #confusion matrix, if pred_label_bgm>0.5, we regarded it as 1, if pred_label_bmg<0.5, we regarded it as 0.
  
  
  auc_train_vec<- c(auc_train_vec,auc(train_target_1,pred_prob)) #AUC: area under ROC curve
  ppv_train_vec <- c(ppv_train_vec, ppv(tb)) #PPV: positive predictive value
  npv_train_vec <- c(npv_train_vec, npv(tb)) #NPV: negative predictive value
  acc_train_vec <- c(acc_train_vec, accuracy(tb)) #ACC: accuracy
  
  auc_train_vec_SVM_R<- c(auc_train_vec_SVM_R,auc(train_target_1,pred_prob_SVM_R_p[,1]))
  ppv_train_vec_SVM_R <- c(ppv_train_vec_SVM_R, ppv(tb_SVM_R))
  npv_train_vec_SVM_R <- c(npv_train_vec_SVM_R, npv(tb_SVM_R))
  acc_train_vec_SVM_R <- c(acc_train_vec_SVM_R, accuracy(tb_SVM_R))
  
  auc_train_vec_SVM_L<- c(auc_train_vec_SVM_L,auc(train_target_1,pred_prob_SVM_L_p[,1]))
  ppv_train_vec_SVM_L <- c(ppv_train_vec_SVM_L, ppv(tb_SVM_L))
  npv_train_vec_SVM_L <- c(npv_train_vec_SVM_L, npv(tb_SVM_L))
  acc_train_vec_SVM_L <- c(acc_train_vec_SVM_L, accuracy(tb_SVM_L))
  
  auc_train_vec_RF<- c(auc_train_vec_RF,auc(train_target_1,pred_prob_RF[,2]))
  ppv_train_vec_RF <- c(ppv_train_vec_RF, ppv(tb_RF))
  npv_train_vec_RF <- c(npv_train_vec_RF, npv(tb_RF))
  acc_train_vec_RF <- c(acc_train_vec_RF, accuracy(tb_RF))
  
  auc_train_vec_DT<- c(auc_train_vec_DT,auc(train_target_1,pred_prob_DT[,2]))
  ppv_train_vec_DT <- c(ppv_train_vec_DT, ppv(tb_DT))
  npv_train_vec_DT <- c(npv_train_vec_DT, npv(tb_DT))
  acc_train_vec_DT <- c(acc_train_vec_DT, accuracy(tb_DT))
  
  auc_train_vec_bgm<- c(auc_train_vec_bgm,auc(train_target_1,pred_prob_bgm[,2]))
  ppv_train_vec_bgm <- c(ppv_train_vec_bgm, ppv(tb_bgm))
  npv_train_vec_bgm <- c(npv_train_vec_bgm, npv(tb_bgm))
  acc_train_vec_bgm <- c(acc_train_vec_bgm, accuracy(tb_bgm))
  
  #testing
  pred_prob_2 <- predict(linearMod,test_df_1,type="response")
  pred_label_2 <- ifelse(pred_prob_2 > 0.5,1,0)
  pred_label_2 <- as.factor(pred_label_2)
  tb_2 <- table(pred_label_2, test_target_1)
  
  
  pred_prob_SVM_R_2 <- predict(model_SVM_R, test_df_1, decision.values = TRUE, probability = TRUE)
  pred_prob_SVM_R_p_2 <- attr(pred_prob_SVM_R_2, "probabilities")
  pred_prob_SVM_R_d_2 <- attr(pred_prob_SVM_R_2, "decision.values")
  pred_label_SVM_R_2 <- ifelse(pred_prob_SVM_R_p_2[,1] > 0.5,1,0)
  pred_label_SVM_R_2 <- as.factor(pred_label_SVM_R_2)
  tb_SVM_R_2 <- table(pred_label_SVM_R_2, test_target_1)
  
  
  pred_prob_SVM_L_2 <- predict(model_SVM_L, test_df_1, decision.values = TRUE, probability = TRUE)
  pred_prob_SVM_L_p_2 <- attr(pred_prob_SVM_L_2, "probabilities")
  pred_prob_SVM_L_d_2 <- attr(pred_prob_SVM_L_2, "decision.values")
  pred_label_SVM_L_2 <- ifelse(pred_prob_SVM_L_p_2[,1] > 0.5,1,0)
  pred_label_SVM_L_2 <- as.factor(pred_label_SVM_L_2)
  tb_SVM_L_2 <- table(pred_label_SVM_L_2, test_target_1)
  
  pred_prob_RF_2 <- predict(model_RF,test_df_1,type="prob")
  pred_label_RF_2 <- ifelse(pred_prob_RF_2[,2] > 0.5,1,0)
  pred_label_RF_2 <- as.factor(pred_label_RF_2)
  tb_RF_2 <- table(pred_label_RF_2, test_target_1)
  
  
  pred_prob_DT_2 <- predict(model_DT,test_df_1,type="prob")
  pred_label_DT_2 <- ifelse(pred_prob_DT_2[,2] > 0.5,1,0)
  pred_label_DT_2 <- as.factor(pred_label_DT_2)
  tb_DT_2 <- table(pred_label_DT_2, test_target_1)
  
  
  pred_prob_bgm_2 <- predict(xgbMod_x,test_df_1,type="prob")
  pred_label_bgm_2 <- ifelse(pred_prob_bgm_2[,2] > 0.5,1,0)
  pred_label_bgm_2 <- as.factor(pred_label_bgm_2)
  tb_bgm_2 <- table(pred_label_bgm_2, test_target_1)
  
  
  auc_test_vec<- c(auc_test_vec,auc(test_target_1,pred_prob_2))
  ppv_test_vec <- c(ppv_test_vec, ppv(tb_2))
  npv_test_vec <- c(npv_test_vec, npv(tb_2))
  acc_test_vec <- c(acc_test_vec, accuracy(tb_2))
  
  auc_test_vec_SVM_R<- c(auc_test_vec_SVM_R,auc(test_target_1,pred_prob_SVM_R_p_2[,1]))
  ppv_test_vec_SVM_R <- c(ppv_test_vec_SVM_R, ppv(tb_SVM_R_2))
  npv_test_vec_SVM_R <- c(npv_test_vec_SVM_R, npv(tb_SVM_R_2))
  acc_test_vec_SVM_R <- c(acc_test_vec_SVM_R, accuracy(tb_SVM_R_2))
  
  auc_test_vec_SVM_L<- c(auc_test_vec_SVM_L,auc(test_target_1,pred_prob_SVM_L_p_2[,1]))
  ppv_test_vec_SVM_L <- c(ppv_test_vec_SVM_L, ppv(tb_SVM_L_2))
  npv_test_vec_SVM_L <- c(npv_test_vec_SVM_L, npv(tb_SVM_L_2))
  acc_test_vec_SVM_L <- c(acc_test_vec_SVM_L, accuracy(tb_SVM_L_2))
  
  auc_test_vec_RF<- c(auc_test_vec_RF,auc(test_target_1,pred_prob_RF_2[,2]))
  ppv_test_vec_RF <- c(ppv_test_vec_RF, ppv(tb_RF_2))
  npv_test_vec_RF <- c(npv_test_vec_RF, npv(tb_RF_2))
  acc_test_vec_RF <- c(acc_test_vec_RF, accuracy(tb_RF_2))
  
  
  auc_test_vec_DT<- c(auc_test_vec_DT,auc(test_target_1,pred_prob_DT_2[,2]))
  ppv_test_vec_DT <- c(ppv_test_vec_DT, ppv(tb_DT_2))
  npv_test_vec_DT <- c(npv_test_vec_DT, npv(tb_DT_2))
  acc_test_vec_DT <- c(acc_test_vec_DT, accuracy(tb_DT_2))
  
  
  auc_test_vec_bgm<- c(auc_test_vec_bgm,auc(test_target_1,pred_prob_bgm_2[,2]))
  ppv_test_vec_bgm <- c(ppv_test_vec_bgm, ppv(tb_bgm_2))
  npv_test_vec_bgm <- c(npv_test_vec_bgm, npv(tb_bgm_2))
  acc_test_vec_bgm <- c(acc_test_vec_bgm, accuracy(tb_bgm_2))
}

acc_train_mean <-mean(acc_train_vec,na.rm=T); acc_train_sd <-sd(acc_train_vec,na.rm=T)
acc_test_mean <- mean(acc_test_vec,na.rm=T); acc_test_sd <-sd(acc_test_vec,na.rm=T)
acc_train_mean_SVM_R <-mean(acc_train_vec_SVM_R,na.rm=T); acc_train_sd_SVM_R <-sd(acc_train_vec_SVM_R,na.rm=T)
acc_test_mean_SVM_R <- mean(acc_test_vec_SVM_R,na.rm=T); acc_test_sd_SVM_R <-sd(acc_test_vec_SVM_R,na.rm=T)
acc_train_mean_SVM_L <-mean(acc_train_vec_SVM_L,na.rm=T); acc_train_sd_SVM_L <-sd(acc_train_vec_SVM_L,na.rm=T)
acc_test_mean_SVM_L <- mean(acc_test_vec_SVM_L,na.rm=T); acc_test_sd_SVM_L <-sd(acc_test_vec_SVM_L,na.rm=T)
acc_train_mean_RF <-mean(acc_train_vec_RF,na.rm=T); acc_train_sd_RF <-sd(acc_train_vec_RF,na.rm=T)
acc_test_mean_RF <- mean(acc_test_vec_RF,na.rm=T); acc_test_sd_RF <-sd(acc_test_vec_RF,na.rm=T)
acc_train_mean_DT <-mean(acc_train_vec_DT,na.rm=T); acc_train_sd_DT <-sd(acc_train_vec_DT,na.rm=T)
acc_test_mean_DT <- mean(acc_test_vec_DT,na.rm=T); acc_test_sd_DT <-sd(acc_test_vec_DT,na.rm=T)
acc_train_mean_bgm <-mean(acc_train_vec_bgm,na.rm=T); acc_train_sd_bgm <-sd(acc_train_vec_bgm,na.rm=T)
acc_test_mean_bgm <- mean(acc_test_vec_bgm,na.rm=T); acc_test_sd_bgm <-sd(acc_test_vec_bgm,na.rm=T)

ppv_train_mean <- mean(ppv_train_vec,na.rm=T);ppv_train_sd <- sd(ppv_train_vec,na.rm=T)
ppv_test_mean <- mean(ppv_test_vec,na.rm=T);ppv_test_sd <- sd(ppv_test_vec,na.rm=T)
ppv_train_mean_SVM_R <-mean(ppv_train_vec_SVM_R,na.rm=T); ppv_train_sd_SVM_R <-sd(ppv_train_vec_SVM_R,na.rm=T)
ppv_test_mean_SVM_R <- mean(ppv_test_vec_SVM_R,na.rm=T); ppv_test_sd_SVM_R <-sd(ppv_test_vec_SVM_R,na.rm=T)
ppv_train_mean_SVM_L <-mean(ppv_train_vec_SVM_L,na.rm=T); ppv_train_sd_SVM_L <-sd(ppv_train_vec_SVM_L,na.rm=T)
ppv_test_mean_SVM_L <- mean(ppv_test_vec_SVM_L,na.rm=T); ppv_test_sd_SVM_L <-sd(ppv_test_vec_SVM_L,na.rm=T)
ppv_train_mean_RF <-mean(ppv_train_vec_RF,na.rm=T); ppv_train_sd_RF <-sd(ppv_train_vec_RF,na.rm=T)
ppv_test_mean_RF <- mean(ppv_test_vec_RF,na.rm=T); ppv_test_sd_RF <-sd(ppv_test_vec_RF,na.rm=T)
ppv_train_mean_DT <-mean(ppv_train_vec_DT,na.rm=T); ppv_train_sd_DT <-sd(ppv_train_vec_DT,na.rm=T)
ppv_test_mean_DT <- mean(ppv_test_vec_DT,na.rm=T); ppv_test_sd_DT <-sd(ppv_test_vec_DT,na.rm=T)
ppv_train_mean_bgm <-mean(ppv_train_vec_bgm,na.rm=T); ppv_train_sd_bgm <-sd(ppv_train_vec_bgm,na.rm=T)
ppv_test_mean_bgm <- mean(ppv_test_vec_bgm,na.rm=T); ppv_test_sd_bgm <-sd(ppv_test_vec_bgm,na.rm=T)
npv_train_mean <- mean(npv_train_vec,na.rm=T);npv_train_sd <- sd(npv_train_vec,na.rm=T)
npv_test_mean <- mean(npv_test_vec,na.rm=T);npv_test_sd <-sd(npv_test_vec,na.rm=T)
npv_train_mean_SVM_R <-mean(npv_train_vec_SVM_R,na.rm=T); npv_train_sd_SVM_R <-sd(npv_train_vec_SVM_R,na.rm=T)
npv_test_mean_SVM_R <- mean(npv_test_vec_SVM_R,na.rm=T); npv_test_sd_SVM_R <-sd(npv_test_vec_SVM_R,na.rm=T)
npv_train_mean_SVM_L <-mean(npv_train_vec_SVM_L,na.rm=T); npv_train_sd_SVM_L <-sd(npv_train_vec_SVM_L,na.rm=T)
npv_test_mean_SVM_L <- mean(npv_test_vec_SVM_L,na.rm=T); npv_test_sd_SVM_L <-sd(npv_test_vec_SVM_L,na.rm=T)
npv_train_mean_RF <-mean(npv_train_vec_RF,na.rm=T); npv_train_sd_RF <-sd(npv_train_vec_RF,na.rm=T)
npv_test_mean_RF <- mean(npv_test_vec_RF,na.rm=T); npv_test_sd_RF <-sd(npv_test_vec_RF,na.rm=T)
npv_train_mean_DT <-mean(npv_train_vec_DT,na.rm=T); npv_train_sd_DT <-sd(npv_train_vec_DT,na.rm=T)
npv_test_mean_DT <- mean(npv_test_vec_DT,na.rm=T); npv_test_sd_DT <-sd(npv_test_vec_DT,na.rm=T)
npv_train_mean_bgm <-mean(npv_train_vec_bgm,na.rm=T); npv_train_sd_bgm <-sd(npv_train_vec_bgm,na.rm=T)
npv_test_mean_bgm <- mean(npv_test_vec_bgm,na.rm=T); npv_test_sd_bgm <-sd(npv_test_vec_bgm,na.rm=T)
auc_train_mean <- mean(auc_train_vec,na.rm=T);auc_train_sd <-sd(auc_train_vec,na.rm=T)
auc_test_mean <- mean(auc_test_vec,na.rm=T);auc_test_sd <-sd(auc_test_vec,na.rm=T)
auc_train_mean_SVM_R <-mean(auc_train_vec_SVM_R,na.rm=T); auc_train_sd_SVM_R <-sd(auc_train_vec_SVM_R,na.rm=T)
auc_test_mean_SVM_R <- mean(auc_test_vec_SVM_R,na.rm=T); auc_test_sd_SVM_R <-sd(auc_test_vec_SVM_R,na.rm=T)
auc_train_mean_SVM_L <-mean(auc_train_vec_SVM_L,na.rm=T); auc_train_sd_SVM_L <-sd(auc_train_vec_SVM_L,na.rm=T)
auc_test_mean_SVM_L <- mean(auc_test_vec_SVM_L,na.rm=T); auc_test_sd_SVM_L <-sd(auc_test_vec_SVM_L,na.rm=T)
auc_train_mean_RF <-mean(auc_train_vec_RF,na.rm=T); auc_train_sd_RF <-sd(auc_train_vec_RF,na.rm=T)
auc_test_mean_RF <- mean(auc_test_vec_RF,na.rm=T); auc_test_sd_RF <-sd(auc_test_vec_RF,na.rm=T)
auc_train_mean_DT <-mean(auc_train_vec_DT,na.rm=T); auc_train_sd_DT <-sd(auc_train_vec_DT,na.rm=T)
auc_test_mean_DT <- mean(auc_test_vec_DT,na.rm=T); auc_test_sd_DT <-sd(auc_test_vec_DT,na.rm=T)
auc_train_mean_bgm <-mean(auc_train_vec_bgm,na.rm=T); auc_train_sd_bgm <-sd(auc_train_vec_bgm,na.rm=T)
auc_test_mean_bgm <- mean(auc_test_vec_bgm,na.rm=T); auc_test_sd_bgm <-sd(auc_test_vec_bgm,na.rm=T)

predLiver_relapse_absoluteX <- data.frame(acc_train_mean=acc_train_mean, acc_train_sd=acc_train_sd,
                                          acc_test_mean=acc_test_mean, acc_test_sd=acc_test_sd,
                                          acc_train_mean_SVM_R=acc_train_mean_SVM_R, acc_train_sd_SVM_R=acc_train_sd_SVM_R,
                                          acc_test_mean_SVM_R=acc_test_mean_SVM_R, acc_test_sd_SVM_R=acc_test_sd_SVM_R,
                                          acc_train_mean_SVM_L=acc_train_mean_SVM_L, acc_train_sd_SVM_L=acc_train_sd_SVM_L,
                                          acc_test_mean_SVM_L=acc_test_mean_SVM_L, acc_test_sd_SVM_L=acc_test_sd_SVM_L,
                                          acc_train_mean_RF=acc_train_mean_RF, acc_train_sd_RF=acc_train_sd_RF,
                                          acc_test_mean_RF=acc_test_mean_RF, acc_test_sd_RF=acc_test_sd_RF,
                                          acc_train_mean_DT=acc_train_mean_DT, acc_train_sd_DT=acc_train_sd_DT,
                                          acc_test_mean_DT=acc_test_mean_DT, acc_test_sd_DT=acc_test_sd_DT,
                                          acc_train_mean_bgm=acc_train_mean_bgm, acc_train_sd_bgm=acc_train_sd_bgm,
                                          acc_test_mean_bgm=acc_test_mean_bgm, acc_test_sd_bgm=acc_test_sd_bgm,
                                          ppv_train_mean=ppv_train_mean, ppv_train_sd=ppv_train_sd,
                                          ppv_test_mean=ppv_test_mean, ppv_test_sd=ppv_test_sd,
                                          ppv_train_mean_SVM_R=ppv_train_mean_SVM_R, ppv_train_sd_SVM_R=ppv_train_sd_SVM_R,
                                          ppv_test_mean_SVM_R=ppv_test_mean_SVM_R, ppv_test_sd_SVM_R=ppv_test_sd_SVM_R,
                                          ppv_train_mean_SVM_L=ppv_train_mean_SVM_L, ppv_train_sd_SVM_L=ppv_train_sd_SVM_L,
                                          ppv_test_mean_SVM_L=ppv_test_mean_SVM_L, ppv_test_sd_SVM_L=ppv_test_sd_SVM_L,
                                          ppv_train_mean_RF=ppv_train_mean_RF, ppv_train_sd_RF=ppv_train_sd_RF,
                                          ppv_test_mean_RF=ppv_test_mean_RF, ppv_test_sd_RF=ppv_test_sd_RF,
                                          ppv_train_mean_DT=ppv_train_mean_DT, ppv_train_sd_DT=ppv_train_sd_DT,
                                          ppv_test_mean_DT=ppv_test_mean_DT, ppv_test_sd_DT=ppv_test_sd_DT,
                                          ppv_train_mean_bgm=ppv_train_mean_bgm, ppv_train_sd_bgm=ppv_train_sd_bgm,
                                          ppv_test_mean_bgm=ppv_test_mean_bgm, ppv_test_sd_bgm=ppv_test_sd_bgm,
                                          npv_train_mean=npv_train_mean, npv_train_sd=npv_train_sd,
                                          npv_test_mean=npv_test_mean, npv_test_sd=npv_test_sd,
                                          npv_train_mean_SVM_R=npv_train_mean_SVM_R, npv_train_sd_SVM_R=npv_train_sd_SVM_R,
                                          npv_test_mean_SVM_R=npv_test_mean_SVM_R, npv_test_sd_SVM_R=npv_test_sd_SVM_R,
                                          npv_train_mean_SVM_L=npv_train_mean_SVM_L, npv_train_sd_SVM_L=npv_train_sd_SVM_L,
                                          npv_test_mean_SVM_L=npv_test_mean_SVM_L, npv_test_sd_SVM_L=npv_test_sd_SVM_L,
                                          npv_train_mean_RF=npv_train_mean_RF, npv_train_sd_RF=npv_train_sd_RF,
                                          npv_test_mean_RF=npv_test_mean_RF, npv_test_sd_RF=npv_test_sd_RF,
                                          npv_train_mean_DT=npv_train_mean_DT, npv_train_sd_DT=npv_train_sd_DT,
                                          npv_test_mean_DT=npv_test_mean_DT, npv_test_sd_DT=npv_test_sd_DT,
                                          npv_train_mean_bgm=npv_train_mean_bgm, npv_train_sd_bgm=npv_train_sd_bgm,
                                          npv_test_mean_bgm=npv_test_mean_bgm, npv_test_sd_bgm=npv_test_sd_bgm,
                                          auc_train_mean=auc_train_mean, auc_train_sd=auc_train_sd,
                                          auc_test_mean=auc_test_mean, auc_test_sd=auc_test_sd,
                                          auc_train_mean_SVM_R=auc_train_mean_SVM_R, auc_train_sd_SVM_R=auc_train_sd_SVM_R,
                                          auc_test_mean_SVM_R=auc_test_mean_SVM_R, auc_test_sd_SVM_R=auc_test_sd_SVM_R,
                                          auc_train_mean_SVM_L=auc_train_mean_SVM_L, auc_train_sd_SVM_L=auc_train_sd_SVM_L,
                                          auc_test_mean_SVM_L=auc_test_mean_SVM_L, auc_test_sd_SVM_L=auc_test_sd_SVM_L,
                                          auc_train_mean_RF=auc_train_mean_RF, auc_train_sd_RF=auc_train_sd_RF,
                                          auc_test_mean_RF=auc_test_mean_RF, auc_test_sd_RF=auc_test_sd_RF,
                                          auc_train_mean_DT=auc_train_mean_DT, auc_train_sd_DT=auc_train_sd_DT,
                                          auc_test_mean_DT=auc_test_mean_DT, auc_test_sd_DT=auc_test_sd_DT,
                                          auc_train_mean_bgm=auc_train_mean_bgm, auc_train_sd_bgm=auc_train_sd_bgm,
                                          auc_test_mean_bgm=auc_test_mean_bgm, auc_test_sd_bgm=auc_test_sd_bgm, stringsAsFactors = F)

#write.csv(predLiver_relapse_absoluteX,file="SPON2_VIM_LIVER_NEW_performace_FINAL.csv")


