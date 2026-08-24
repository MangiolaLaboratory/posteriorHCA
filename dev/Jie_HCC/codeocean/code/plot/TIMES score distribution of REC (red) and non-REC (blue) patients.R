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
#[1] compiler_4.2.1 tools_4.2.1   


setwd("../data/Input_2_TIMES_development_from_mIHC/")

##read in data
load("TIMES_distribution_input.RData")

#graph about performance of wilcoxon test between recurrent and non-recurrent ROIs
condi_shift <- c();pval <- c()

pdf("hist_B256.pdf")  
par(mfrow=c(5,4),mar=c(2,4,1,1)) 

b4x <- liver_r_data_gene_ed #recurrent TIMES
x <- liver_nr_data_gene_ed  #non-recurrent TIMES

wil_t <- wilcox.test(x, b4x,alternative = c("two.sided"),mu = 0,paired=F,exact=F,correct=F,conf.int=T,conf.level=0.95) #wilcoxon test

xmin <- -0.2; xmax <- 1.2; break_seq <- seq(from=xmin,to=xmax,by=0.1); ymin <- 0; ymax <- 6
if (wil_t$estimate<0){ #non-recurrent distribution < recurrent distribution
  my.hist <- hist(x, xlim=c(xmin,xmax), ylim=c(ymin,ymax),main="", xlab="", ylab="PD&PR",breaks=break_seq,freq=F,lty="blank",las=1)
  lines(c(my.hist$breaks, max(my.hist$breaks)),c(0,my.hist$density,0), type='S',col="#2367AC",lwd=1.5)
  
  my.hist <- hist(b4x, add=T, breaks=break_seq,freq=F,lty="blank")
  lines(c(my.hist$breaks, max(my.hist$breaks)),c(0,my.hist$density,0), type='S',col="#CC4B40",lwd=1.5)
  
  text(1,2.2,labels=paste("p = ",formatC(wil_t$p.value,format = "e",digits = 2),sep=""))
  text(1,1.9,labels=paste("shift = ",formatC(wil_t$estimate,format = "e",digits = 2),sep=""))
  text(-0.2,0.3,labels="NR",col="#2367AC",font=2)
  text(1,0.3,labels="R",col="#CC4B40",font=2)
} else { #recurrent distribution < non-recurrent distribution
  my.hist <- hist(b4x, xlim=c(xmin,xmax), ylim=c(ymin,ymax),main="", xlab="", ylab="PD&PR_",breaks=break_seq,freq=F,lty="blank",las=1)
  lines(c(my.hist$breaks, max(my.hist$breaks)),c(0,my.hist$density,0), type='S',col="#2367AC",lwd=1.5)
  
  my.hist <- hist(x, add=T, breaks=break_seq,freq=F,lty="blank")
  lines(c(my.hist$breaks, max(my.hist$breaks)),c(0,my.hist$density,0), type='S',col="#CC4B40",lwd=1.5)
  
  text(1,2.2,labels=paste("p = ",formatC(wil_t$p.value,format = "e",digits = 2),sep=""))
  text(1,1.9,labels=paste("shift = ",formatC(wil_t$estimate,format = "e",digits = 2),sep=""))
  text(-0.2,0.3,labels="R",col="#2367AC",font=2)
  text(1,0.3,labels="NR",col="#CC4B40",font=2)
}

dev.off()


