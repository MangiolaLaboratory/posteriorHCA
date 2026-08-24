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

#loaded via a namespace (and not attached):
#[1] compiler_4.2.1    tools_4.2.1       rstudioapi_0.15.0


setwd("../data/Input_1_Spatial_transcriptomics_analysis/")

###read in data
load("Examined_their_capabilities.RData")

##graph about performance of transcriptome 23 genes
pdf(paste("hist_red_chart",".pdf",sep=""),width=8.25,height=11.75)
par(oma = c(1,1,3,1), mfrow = c(10,1),mar=c(3,3,2,1))

liver_re_acc_Mean_vec <- data2$xpluscova_acc_test_mean #23 genes accuracy
liver_re_acc_SD_vec <- data2$xpluscova_acc_test_sd #23 genes standard deviation of accuracy
liverE_re_acc_vec <- 1.960*liver_re_acc_SD_vec/sqrt(23) #23 genes standard error of accuracy

tabbedMeans1 <- matrix(data =NA,nrow=length(liver_re_acc_Mean_vec),ncol=1)
tabbedMeans1[,1] <- liver_re_acc_Mean_vec

tabbedE1 <- matrix(data=NA,nrow=length(liver_re_acc_Mean_vec),ncol=1)
tabbedE1[,1] <- liverE_re_acc_vec

tabbedcol1 <- matrix(data=NA,nrow=length(liver_re_acc_Mean_vec),ncol=1)
tabbedcol1[,1] <- adjustcolor("brown", alpha.f=0.59)

barCenter1 <- barplot(height = tabbedMeans1,beside = TRUE, las = 2,ylim = c(0.5, 1),cex.names = 1, 
                      main = " Liver  Recurrent & Non-recurrent",ylab = "Perfomance",border = "white", axes = TRUE,legend.text = TRUE,
                      col=tabbedcol1,args.legend = list(x = "x",cex = .7),angle=30,names.arg = top_vec,width=5,space=0.5)

segments(barCenter1, tabbedMeans1 - tabbedE1, barCenter1,
         tabbedMeans1 + tabbedE1, lwd = 1)

arrows(barCenter1, tabbedMeans1 - tabbedE1, barCenter1,
       tabbedMeans1 + tabbedE1, lwd = 1, angle = 90, code = 3, length = 0.03)

dev.off()
