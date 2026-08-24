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


setwd("../data/Input_2_TIMES_development_from_mIHC/")

###read in data
load("prediction_accuracy_input.RData") 

##graph about accuracy performance of 5 gene ROIs
pdf(paste("hist_protein_per_chart",".pdf",sep=""),width=8.25,height=11.75)
par(oma = c(1,1,3,1), mfrow = c(5,3),mar=c(3,3,2,1))
barplot(height=c,names.arg = gene,ylim = c(-0.02,0.08),col = "lightblue")
dev.off()
