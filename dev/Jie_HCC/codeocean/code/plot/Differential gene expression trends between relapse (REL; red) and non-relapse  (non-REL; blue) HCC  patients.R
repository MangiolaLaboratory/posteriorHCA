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


#####################liver nike plot###########################
setwd("../data/Input_1_Spatial_transcriptomics_analysis/")

###read in data
load("differential_gene_expression.RData")


pdf(paste("line_chart",".pdf",sep=""),width=8.25,height=11.75)
par(oma = c(1,1,3,1), mfrow = c(2,2),mar=c(3,3,2,1))

##graph about flat-up(non-recurrent)_down-flat(recurrent) trajectory####

gene_show <- fu_df$gene

#recurrent
Relp_NKatStroma <- unlist(data_gene_df_ed[gene_show,paste(rep("...",length(NKatStroma_indx)),NKatStroma_indx,sep="")])
Relp_NKatBorder <- unlist(data_gene_df_ed[gene_show,paste(rep("...",length(Relp_NKatBorder_indx)),Relp_NKatBorder_indx,sep="")])
Relp_NKatTumor <- unlist(data_gene_df_ed[gene_show,paste(rep("...",length(Relp_NKatTumor_indx)),Relp_NKatTumor_indx,sep="")])

#non-recurrent
Nonrelp_NKatStroma<- unlist(data_gene_df_ed[gene_show,paste(rep("...",length(NKatStroma_indx)),NKatStroma_indx,sep="")])
Nonrelp_NKatBorder <- unlist(data_gene_df_ed[gene_show,paste(rep("...",length(Nonrelp_NKatBorder_indx)),Nonrelp_NKatBorder_indx,sep="")])
Nonrelp_NKatTumor <- unlist(data_gene_df_ed[gene_show,paste(rep("...",length(Nonrelp_NKatTumor_indx)),Nonrelp_NKatTumor_indx,sep="")])

x_vec <- c(1,2,3)
y_vec_Relp <- c(mean(Relp_NKatStroma,na.rm=T), mean(Relp_NKatBorder,na.rm=T), mean(Relp_NKatTumor,na.rm=T))
y_up_Relp <- c(mean(Relp_NKatStroma,na.rm=T)+1.96/sqrt(19)*sd(Relp_NKatStroma,na.rm=T), mean(Relp_NKatBorder,na.rm=T)+1.96/sqrt(19)*sd(Relp_NKatBorder,na.rm=T), mean(Relp_NKatTumor,na.rm=T)+1.96/sqrt(19)*sd(Relp_NKatTumor,na.rm=T))
y_lw_Relp <- c(mean(Relp_NKatStroma,na.rm=T)-1.96/sqrt(19)*sd(Relp_NKatStroma,na.rm=T), mean(Relp_NKatBorder,na.rm=T)-1.96/sqrt(19)*sd(Relp_NKatBorder,na.rm=T), mean(Relp_NKatTumor,na.rm=T)-1.96/sqrt(19)*sd(Relp_NKatTumor,na.rm=T))

y_vec_Nonrelp <- c(mean(Nonrelp_NKatStroma,na.rm=T), mean(Nonrelp_NKatBorder,na.rm=T), mean(Nonrelp_NKatTumor,na.rm=T))
y_lw_Nonrelp <- c(mean(Nonrelp_NKatStroma,na.rm=T)-1.96/sqrt(19)*sd(Nonrelp_NKatStroma,na.rm=T), mean(Nonrelp_NKatBorder,na.rm=T)-1.96/sqrt(19)*sd(Nonrelp_NKatBorder,na.rm=T), mean(Nonrelp_NKatTumor,na.rm=T)-1.96/sqrt(19)*sd(Nonrelp_NKatTumor,na.rm=T))
y_up_Nonrelp <- c(mean(Nonrelp_NKatStroma,na.rm=T)+1.96/sqrt(19)*sd(Nonrelp_NKatStroma,na.rm=T), mean(Nonrelp_NKatBorder,na.rm=T)+1.96/sqrt(19)*sd(Nonrelp_NKatBorder,na.rm=T),mean(Nonrelp_NKatTumor,na.rm=T)+1.96/sqrt(19)*sd(Nonrelp_NKatTumor,na.rm=T))

ylim_max <- max(max(y_up_Nonrelp,na.rm=T),max(y_up_Relp,na.rm=T))
if (ylim_max == 'NaN') {ylim_max <- 0.3}
if (ylim_max == Inf) {ylim_max <- 0.3}
plot(x = x_vec, y = y_vec_Relp, type = "l", col="#E73334",xlim = c(0.9,3.1), ylim = c(0,ylim_max*1.03),xaxs="i",yaxs="i",xlab="",ylab="",cex.lab=0.75, cex.axis=0.75, cex.main=0.75)
polygon(c(x_vec,rev(x_vec)),c(y_lw_Relp,rev(y_up_Relp)),col = adjustcolor("#E73334", alpha.f=0.8), border = FALSE)
lines(x = x_vec, y = y_vec_Nonrelp, col="#559AC6",lty = 'dashed',xaxs="i",yaxs="i")
polygon(c(x_vec,rev(x_vec)),c(y_lw_Nonrelp,rev(y_up_Nonrelp)),col = adjustcolor("#559AC6", alpha.f=0.8), border = FALSE)
legend("topright", legend=c(paste("Recurrent",sep=": "), paste("Non-recurrent",sep=": ")), col=c("#E73334","#559AC6"), lty=c(1,2), bty="n",cex=1.5)


##graph about flat-up(non-recurrent)_up-flat(recurrent) trajectory####

gene_show_2 <- fu_uf$gene

#recurrent
Relp_NKatStroma<- unlist(data_gene_df_ed[gene_show_2,paste(rep("...",length(NKatStroma_indx)),NKatStroma_indx,sep="")])
Relp_NKatBorder<- unlist(data_gene_df_ed[gene_show_2,paste(rep("...",length(Relp_NKatBorder_indx)),Relp_NKatBorder_indx,sep="")])
Relp_NKatTumor<- unlist(data_gene_df_ed[gene_show_2,paste(rep("...",length(Relp_NKatTumor_indx)),Relp_NKatTumor_indx,sep="")])

#non-recurrent
Nonrelp_NKatStroma<- unlist(data_gene_df_ed[gene_show_2,paste(rep("...",length(NKatStroma_indx)),NKatStroma_indx,sep="")])
Nonrelp_NKatBorder<- unlist(data_gene_df_ed[gene_show_2,paste(rep("...",length(Nonrelp_NKatBorder_indx)),Nonrelp_NKatBorder_indx,sep="")])
Nonrelp_NKatTumor<- unlist(data_gene_df_ed[gene_show_2,paste(rep("...",length(Nonrelp_NKatTumor_indx)),Nonrelp_NKatTumor_indx,sep="")])

x_vec <- c(1,2,3)
y_vec_Relp <- c(mean(Relp_NKatStroma,na.rm=T), mean(Relp_NKatBorder,na.rm=T), mean(Relp_NKatTumor,na.rm=T))
y_up_Relp <- c(mean(Relp_NKatStroma,na.rm=T)+1.96/sqrt(106)*sd(Relp_NKatStroma,na.rm=T), mean(Relp_NKatBorder,na.rm=T)+1.96/sqrt(106)*sd(Relp_NKatBorder,na.rm=T), mean(Relp_NKatTumor,na.rm=T)+1.96/sqrt(106)*sd(Relp_NKatTumor,na.rm=T))
y_lw_Relp <- c(mean(Relp_NKatStroma,na.rm=T)-1.96/sqrt(106)*sd(Relp_NKatStroma,na.rm=T), mean(Relp_NKatBorder,na.rm=T)-1.96/sqrt(106)*sd(Relp_NKatBorder,na.rm=T), mean(Relp_NKatTumor,na.rm=T)-1.96/sqrt(106)*sd(Relp_NKatTumor,na.rm=T))

y_vec_Nonrelp <- c(mean(Nonrelp_NKatStroma,na.rm=T), mean(Nonrelp_NKatBorder,na.rm=T), mean(Nonrelp_NKatTumor,na.rm=T))
y_lw_Nonrelp <- c(mean(Nonrelp_NKatStroma,na.rm=T)-1.96/sqrt(106)*sd(Nonrelp_NKatStroma,na.rm=T), mean(Nonrelp_NKatBorder,na.rm=T)-1.96/sqrt(106)*sd(Nonrelp_NKatBorder,na.rm=T), mean(Nonrelp_NKatTumor,na.rm=T)-1.96/sqrt(106)*sd(Nonrelp_NKatTumor,na.rm=T))
y_up_Nonrelp <- c(mean(Nonrelp_NKatStroma,na.rm=T)+1.96/sqrt(106)*sd(Nonrelp_NKatStroma,na.rm=T), mean(Nonrelp_NKatBorder,na.rm=T)+1.96/sqrt(106)*sd(Nonrelp_NKatBorder,na.rm=T),mean(Nonrelp_NKatTumor,na.rm=T)+1.96/sqrt(106)*sd(Nonrelp_NKatTumor,na.rm=T))

ylim_max <- max(max(y_up_Nonrelp,na.rm=T),max(y_up_Relp,na.rm=T))
if (ylim_max == 'NaN') {ylim_max <- 0.3}
if (ylim_max == Inf) {ylim_max <- 0.3}
plot(x = x_vec, y = y_vec_Relp, type = "l", col="#E73334",xlim = c(0.9,3.1), ylim = c(0,ylim_max*1.03),xaxs="i",yaxs="i",xlab="",ylab="",cex.lab=0.75, cex.axis=0.75, cex.main=0.75)
polygon(c(x_vec,rev(x_vec)),c(y_lw_Relp,rev(y_up_Relp)),col = adjustcolor("#E73334", alpha.f=0.8), border = FALSE)
lines(x = x_vec, y = y_vec_Nonrelp, col="#559AC6",lty = 'dashed',xaxs="i",yaxs="i")
polygon(c(x_vec,rev(x_vec)),c(y_lw_Nonrelp,rev(y_up_Nonrelp)),col = adjustcolor("#559AC6", alpha.f=0.8), border = FALSE)
legend("topright", legend=c(paste("Recurrent",sep=": "), paste("Non-recurrent",sep=": ")), col=c("#E73334","#559AC6"), lty=c(1,2), bty="n",cex=1.5)


##graph about flat-down(non-recurrent)_down-flat(recurrent) trajectory####

gene_show_3 <- fd_df$gene

#recurrent
Relp_NKatStroma<- unlist(data_gene_df_ed[gene_show_3,paste(rep("...",length(NKatStroma_indx)),NKatStroma_indx,sep="")])
Relp_NKatBorder<- unlist(data_gene_df_ed[gene_show_3,paste(rep("...",length(Relp_NKatBorder_indx)),Relp_NKatBorder_indx,sep="")])
Relp_NKatTumor<- unlist(data_gene_df_ed[gene_show_3,paste(rep("...",length(Relp_NKatTumor_indx)),Relp_NKatTumor_indx,sep="")])

#non-recurrent
Nonrelp_NKatStroma<- unlist(data_gene_df_ed[gene_show_3,paste(rep("...",length(NKatStroma_indx)),NKatStroma_indx,sep="")])
Nonrelp_NKatBorder<- unlist(data_gene_df_ed[gene_show_3,paste(rep("...",length(Nonrelp_NKatBorder_indx)),Nonrelp_NKatBorder_indx,sep="")])
Nonrelp_NKatTumor<- unlist(data_gene_df_ed[gene_show_3,paste(rep("...",length(Nonrelp_NKatTumor_indx)),Nonrelp_NKatTumor_indx,sep="")])

x_vec <- c(1,2,3)

#recurrent
y_vec_Relp <- c(mean(Relp_NKatStroma,na.rm=T), mean(Relp_NKatBorder,na.rm=T), mean(Relp_NKatTumor,na.rm=T))
y_up_Relp <- c(mean(Relp_NKatStroma,na.rm=T)+1.96/sqrt(113)*sd(Relp_NKatStroma,na.rm=T), mean(Relp_NKatBorder,na.rm=T)+1.96/sqrt(113)*sd(Relp_NKatBorder,na.rm=T), mean(Relp_NKatTumor,na.rm=T)+1.96/sqrt(113)*sd(Relp_NKatTumor,na.rm=T))
y_lw_Relp <- c(mean(Relp_NKatStroma,na.rm=T)-1.96/sqrt(113)*sd(Relp_NKatStroma,na.rm=T), mean(Relp_NKatBorder,na.rm=T)-1.96/sqrt(113)*sd(Relp_NKatBorder,na.rm=T), mean(Relp_NKatTumor,na.rm=T)-1.96/sqrt(113)*sd(Relp_NKatTumor,na.rm=T))

#non-recurrent
y_vec_Nonrelp <- c(mean(Nonrelp_NKatStroma,na.rm=T), mean(Nonrelp_NKatBorder,na.rm=T), mean(Nonrelp_NKatTumor,na.rm=T))
y_lw_Nonrelp <- c(mean(Nonrelp_NKatStroma,na.rm=T)-1.96/sqrt(113)*sd(Nonrelp_NKatStroma,na.rm=T), mean(Nonrelp_NKatBorder,na.rm=T)-1.96/sqrt(113)*sd(Nonrelp_NKatBorder,na.rm=T), mean(Nonrelp_NKatTumor,na.rm=T)-1.96/sqrt(113)*sd(Nonrelp_NKatTumor,na.rm=T))
y_up_Nonrelp <- c(mean(Nonrelp_NKatStroma,na.rm=T)+1.96/sqrt(113)*sd(Nonrelp_NKatStroma,na.rm=T), mean(Nonrelp_NKatBorder,na.rm=T)+1.96/sqrt(113)*sd(Nonrelp_NKatBorder,na.rm=T),mean(Nonrelp_NKatTumor,na.rm=T)+1.96/sqrt(113)*sd(Nonrelp_NKatTumor,na.rm=T))

ylim_max <- max(max(y_up_Nonrelp,na.rm=T),max(y_up_Relp,na.rm=T))
if (ylim_max == 'NaN') {ylim_max <- 0.3}
if (ylim_max == Inf) {ylim_max <- 0.3}
plot(x = x_vec, y = y_vec_Relp, type = "l", col="#E73334",xlim = c(0.9,3.1), ylim = c(0,ylim_max*1.03),xaxs="i",yaxs="i",xlab="",ylab="",cex.lab=0.75, cex.axis=0.75, cex.main=0.75)
polygon(c(x_vec,rev(x_vec)),c(y_lw_Relp,rev(y_up_Relp)),col = adjustcolor("#E73334", alpha.f=0.8), border = FALSE)
lines(x = x_vec, y = y_vec_Nonrelp, col="#559AC6",lty = 'dashed',xaxs="i",yaxs="i")
polygon(c(x_vec,rev(x_vec)),c(y_lw_Nonrelp,rev(y_up_Nonrelp)),col = adjustcolor("#559AC6", alpha.f=0.8), border = FALSE)
legend("topright", legend=c(paste("Recurrent",sep=": "), paste("Non-recurrent",sep=": ")), col=c("#E73334","#559AC6"), lty=c(1,2), bty="n",cex=1.5)


##graph about flat-down(non-recurrent)_up-flat(recurrent) trajectory####

gene_show <- fd_uf$gene

## line of gene UGT2B15  ##
type_1 <- unlist(fd_uf[which(fd_uf$gene=="UGT2B15"),c("relp","nonrelp")])

#recurrent
Relp_NKatStroma_1 <- unlist(data_gene_df_ed["UGT2B15",paste(rep("...",length(NKatStroma_indx)),NKatStroma_indx,sep="")])
Relp_NKatBorder_1 <- unlist(data_gene_df_ed["UGT2B15",paste(rep("...",length(Relp_NKatBorder_indx)),Relp_NKatBorder_indx,sep="")])
Relp_NKatTumor_1 <- unlist(data_gene_df_ed["UGT2B15",paste(rep("...",length(Relp_NKatTumor_indx)),Relp_NKatTumor_indx,sep="")])

#non-recurrent
Nonrelp_NKatStroma_1 <- unlist(data_gene_df_ed["UGT2B15",paste(rep("...",length(NKatStroma_indx)),NKatStroma_indx,sep="")])
Nonrelp_NKatBorder_1 <- unlist(data_gene_df_ed["UGT2B15",paste(rep("...",length(Nonrelp_NKatBorder_indx)),Nonrelp_NKatBorder_indx,sep="")])
Nonrelp_NKatTumor_1 <- unlist(data_gene_df_ed["UGT2B15",paste(rep("...",length(Nonrelp_NKatTumor_indx)),Nonrelp_NKatTumor_indx,sep="")])
cat("UGT2B15","\t",type_1,"\n\t",Nonrelp_NKatStroma_1,"\t",Nonrelp_NKatBorder_1,"\t",Nonrelp_NKatTumor_1,"\n\t",Relp_NKatStroma_1,"\t",Relp_NKatBorder_1,"\t",Relp_NKatTumor_1,"\n")

x_vec <- c(1,2,3)

#recurrent
y_vec_Relp_1 <- c(mean(Relp_NKatStroma_1,na.rm=T), mean(Relp_NKatBorder_1,na.rm=T), mean(Relp_NKatTumor_1,na.rm=T))
y_up_Relp_1 <- c(mean(Relp_NKatStroma_1,na.rm=T)+1.96*sd(Relp_NKatStroma_1,na.rm=T), mean(Relp_NKatBorder_1,na.rm=T)+1.96*sd(Relp_NKatBorder_1,na.rm=T), mean(Relp_NKatTumor_1,na.rm=T)+1.96*sd(Relp_NKatTumor_1,na.rm=T))
y_lw_Relp_1 <- c(mean(Relp_NKatStroma_1,na.rm=T)-1.96*sd(Relp_NKatStroma_1,na.rm=T), mean(Relp_NKatBorder_1,na.rm=T)-sd(Relp_NKatBorder_1,na.rm=T), mean(Relp_NKatTumor_1,na.rm=T)-1.96*sd(Relp_NKatTumor_1,na.rm=T))

#non-recurrent
y_vec_Nonrelp_1 <- c(mean(Nonrelp_NKatStroma_1,na.rm=T), mean(Nonrelp_NKatBorder_1,na.rm=T), mean(Nonrelp_NKatTumor_1,na.rm=T))
y_lw_Nonrelp_1 <- c(mean(Nonrelp_NKatStroma_1,na.rm=T)-1.96*sd(Nonrelp_NKatStroma_1,na.rm=T), mean(Nonrelp_NKatBorder_1,na.rm=T)-1.96*sd(Nonrelp_NKatBorder_1,na.rm=T), mean(Nonrelp_NKatTumor_1,na.rm=T)-1.96*sd(Nonrelp_NKatTumor_1,na.rm=T))
y_up_Nonrelp_1 <- c(mean(Nonrelp_NKatStroma_1,na.rm=T)+1.96*sd(Nonrelp_NKatStroma_1,na.rm=T), mean(Nonrelp_NKatBorder_1,na.rm=T)+1.96*sd(Nonrelp_NKatBorder_1,na.rm=T), mean(Nonrelp_NKatTumor_1,na.rm=T)+1.96*sd(Nonrelp_NKatTumor_1,na.rm=T))

ylim_max <- max(max(y_up_Nonrelp_1,na.rm=T),max(y_up_Relp_1,na.rm=T))
if (ylim_max == 'NaN') {ylim_max <- 0.3}
if (ylim_max == Inf) {ylim_max <- 0.3}

## line of gene SPP2  ##
type_2 <- unlist(fd_uf[which(fd_uf$gene=="SPP2"),c("relp","nonrelp")])

#recurrent
Relp_NKatStroma_2 <- unlist(data_gene_df_ed["SPP2",paste(rep("...",length(NKatStroma_indx)),NKatStroma_indx,sep="")])
Relp_NKatBorder_2 <- unlist(data_gene_df_ed["SPP2",paste(rep("...",length(Relp_NKatBorder_indx)),Relp_NKatBorder_indx,sep="")])
Relp_NKatTumor_2 <- unlist(data_gene_df_ed["SPP2",paste(rep("...",length(Relp_NKatTumor_indx)),Relp_NKatTumor_indx,sep="")])

#non-recurrent
Nonrelp_NKatStroma_2 <- unlist(data_gene_df_ed["SPP2",paste(rep("...",length(NKatStroma_indx)),NKatStroma_indx,sep="")])
Nonrelp_NKatBorder_2 <- unlist(data_gene_df_ed["SPP2",paste(rep("...",length(Nonrelp_NKatBorder_indx)),Nonrelp_NKatBorder_indx,sep="")])
Nonrelp_NKatTumor_2 <- unlist(data_gene_df_ed["SPP2",paste(rep("...",length(Nonrelp_NKatTumor_indx)),Nonrelp_NKatTumor_indx,sep="")])
cat("SPP2","\t",type_2,"\n\t",Nonrelp_NKatStroma_2,"\t",Nonrelp_NKatBorder_2,"\t",Nonrelp_NKatTumor_2,"\n\t",Relp_NKatStroma_2,"\t",Relp_NKatBorder_2,"\t",Relp_NKatTumor_2,"\n")

x_vec <- c(1,2,3)

#recurrent
y_vec_Relp_2 <- c(mean(Relp_NKatStroma_2,na.rm=T), mean(Relp_NKatBorder_2,na.rm=T), mean(Relp_NKatTumor_2,na.rm=T))
y_up_Relp_2 <- c(mean(Relp_NKatStroma_2,na.rm=T)+1.96*sd(Relp_NKatStroma_2,na.rm=T), mean(Relp_NKatBorder_2,na.rm=T)+1.96*sd(Relp_NKatBorder_2,na.rm=T), mean(Relp_NKatTumor_2,na.rm=T)+1.96*sd(Relp_NKatTumor_2,na.rm=T))
y_lw_Relp_2 <- c(mean(Relp_NKatStroma_2,na.rm=T)-1.96*sd(Relp_NKatStroma_2,na.rm=T), mean(Relp_NKatBorder_2,na.rm=T)-1.96*sd(Relp_NKatBorder_2,na.rm=T), mean(Relp_NKatTumor_2,na.rm=T)-1.96*sd(Relp_NKatTumor_2,na.rm=T))

#non-recurrent
y_vec_Nonrelp_2 <- c(mean(Nonrelp_NKatStroma_2,na.rm=T), mean(Nonrelp_NKatBorder_2,na.rm=T), mean(Nonrelp_NKatTumor_2,na.rm=T))
y_up_Nonrelp_2 <- c(mean(Nonrelp_NKatStroma_2,na.rm=T)+1.96*sd(Nonrelp_NKatStroma_2,na.rm=T), mean(Nonrelp_NKatBorder_2,na.rm=T)+1.96*sd(Nonrelp_NKatBorder_2,na.rm=T), mean(Nonrelp_NKatTumor_2,na.rm=T)+1.96*sd(Nonrelp_NKatTumor_2,na.rm=T))
y_lw_Nonrelp_2 <- c(mean(Nonrelp_NKatStroma_2,na.rm=T)-1.96*sd(Nonrelp_NKatStroma_2,na.rm=T), mean(Nonrelp_NKatBorder_2,na.rm=T)-1.96*sd(Nonrelp_NKatBorder_2,na.rm=T), mean(Nonrelp_NKatTumor_2,na.rm=T)-1.96*sd(Nonrelp_NKatTumor_2,na.rm=T))


## line of gene SDS  ##
type_3 <- unlist(fd_uf[which(fd_uf$gene=="SDS"),c("relp","nonrelp")])

#recurrent
Relp_NKatStroma_3 <- unlist(data_gene_df_ed["SDS",paste(rep("...",length(NKatStroma_indx)),NKatStroma_indx,sep="")])
Relp_NKatBorder_3 <- unlist(data_gene_df_ed["SDS",paste(rep("...",length(Relp_NKatBorder_indx)),Relp_NKatBorder_indx,sep="")])
Relp_NKatTumor_3 <- unlist(data_gene_df_ed["SDS",paste(rep("...",length(Relp_NKatTumor_indx)),Relp_NKatTumor_indx,sep="")])

#non-recurrent
Nonrelp_NKatStroma_3 <- unlist(data_gene_df_ed["SDS",paste(rep("...",length(NKatStroma_indx)),NKatStroma_indx,sep="")])
Nonrelp_NKatBorder_3 <- unlist(data_gene_df_ed["SDS",paste(rep("...",length(Nonrelp_NKatBorder_indx)),Nonrelp_NKatBorder_indx,sep="")])
Nonrelp_NKatTumor_3 <- unlist(data_gene_df_ed["SDS",paste(rep("...",length(Nonrelp_NKatTumor_indx)),Nonrelp_NKatTumor_indx,sep="")])
cat("SDS","\t",type_3,"\n\t",Nonrelp_NKatStroma_3,"\t",Nonrelp_NKatBorder_3,"\t",Nonrelp_NKatTumor_3,"\n\t",Relp_NKatStroma_3,"\t",Relp_NKatBorder_3,"\t",Relp_NKatTumor_3,"\n")

x_vec <- c(1,2,3)

#recurrent
y_vec_Relp_3 <- c(mean(Relp_NKatStroma_3,na.rm=T), mean(Relp_NKatBorder_3,na.rm=T), mean(Relp_NKatTumor_3,na.rm=T))
y_up_Relp_3 <- c(mean(Relp_NKatStroma_3,na.rm=T)+1.96*sd(Relp_NKatStroma_3,na.rm=T), mean(Relp_NKatBorder_3,na.rm=T)+1.96*sd(Relp_NKatBorder_3,na.rm=T), mean(Relp_NKatTumor_3,na.rm=T)+1.96*sd(Relp_NKatTumor_3,na.rm=T))
y_lw_Relp_3 <- c(mean(Relp_NKatStroma_3,na.rm=T)-1.96*sd(Relp_NKatStroma_3,na.rm=T), mean(Relp_NKatBorder_3,na.rm=T)-1.96*sd(Relp_NKatBorder_3,na.rm=T), mean(Relp_NKatTumor_3,na.rm=T)-1.96*sd(Relp_NKatTumor_3,na.rm=T))

#non-recurrent
y_vec_Nonrelp_3 <- c(mean(Nonrelp_NKatStroma_3,na.rm=T), mean(Nonrelp_NKatBorder_3,na.rm=T), mean(Nonrelp_NKatTumor_3,na.rm=T))
y_up_Nonrelp_3 <- c(mean(Nonrelp_NKatStroma_3,na.rm=T)+1.96*sd(Nonrelp_NKatStroma_3,na.rm=T), mean(Nonrelp_NKatBorder_3,na.rm=T)+1.96*sd(Nonrelp_NKatBorder_3,na.rm=T), mean(Nonrelp_NKatTumor_3,na.rm=T)+1.96*sd(Nonrelp_NKatTumor_3,na.rm=T))
y_lw_Nonrelp_3 <- c(mean(Nonrelp_NKatStroma_3,na.rm=T)-1.96*sd(Nonrelp_NKatStroma_3,na.rm=T), mean(Nonrelp_NKatBorder_3,na.rm=T)-1.96*sd(Nonrelp_NKatBorder_3,na.rm=T), mean(Nonrelp_NKatTumor_3,na.rm=T)-1.96*sd(Nonrelp_NKatTumor_3,na.rm=T))

## add labels  ##
plot(x = x_vec, y = y_vec_Relp_1, type = "l", col="chocolate3",xlim = c(0.9,3.1), ylim = c(0,ylim_max*1.03),xaxs="i",yaxs="i",xlab="",ylab="",cex.lab=1, cex.axis=1, cex.main=1)
lines(x = x_vec, y = y_vec_Relp_2, col="brown3",xaxs="i",yaxs="i")
lines(x = x_vec, y = y_vec_Relp_3, col="mediumpurple2",xaxs="i",yaxs="i")
lines(x = x_vec, y = y_vec_Nonrelp_1, col="cyan2",lty = 'dashed')
lines(x = x_vec, y = y_vec_Nonrelp_2, col="mediumaquamarine",lty = 'dashed')
lines(x = x_vec, y = y_vec_Nonrelp_3, col="lightskyblue2",lty = 'dashed')
legend("topleft", legend=c(paste("UGT2B15_Recurrent",sep=": "), paste("UGT2B15_Non-recurrent",sep=": "), paste("SPP2_Recurrent",sep=": "), paste("SPP2_Non-recurrent",sep=": "), paste("SDS_Recurrent",sep=": "), paste("SDS_Non-recurrent",sep=": ")) ,col=c("chocolate3","cyan2","brown3","mediumaquamarine","mediumpurple2","lightskyblue2"), lty=c(1,2,3), bty="n",cex=1.5)

dev.off()
