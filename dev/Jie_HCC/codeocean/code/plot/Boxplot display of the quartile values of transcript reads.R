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
load("boxplot_input.RData")

##plot of special genens
length(gene_show)/10 #2.3
temp_max <- c(); temp_min <- c()

for (fig_i in 1:3){
  
  pdf(paste("FDR_61gene_nonrelap_VS_relap_figset",fig_i,".pdf",sep=""),width=8.25,height=11.75)
  par(oma = c(1,1,3,1), mfrow = c(5,2),mar=c(3,3,2,1))
  
  if(fig_i==3){
    gene_set <- gene_show[((fig_i-1)*10+1):length(gene_show)]
  } else {
    gene_set <- gene_show[((fig_i-1)*10+1):((fig_i-1)*10+10)]
  }
  
  for (gene_i in gene_set){
    idx_de_i <- which(liver_de_result$gene==gene_i)
    if (length(idx_de_i)==0) cat("error")
    padj_border <- liver_de_result$padj_1Relp_NKatBorder_agst_0Nonrelp_NKatBorder[idx_de_i] #draw non-recurrent lines based on this value
    padj_tumor <- liver_de_result$padj_1Relp_NKatTumor_agst_0Nonrelp_NKatTumor[idx_de_i]
    #draw recurrent lines based on this value
    type_i <- unlist(gene_name[which(gene_name$gene==gene_i),c("nonrelp","relp")])
    
    NKatStroma_i <- unlist(data_gene_df_ed[gene_i,paste(rep("...",length(NKatStroma_indx)),NKatStroma_indx,sep="")]) #AS 
     
    #recurrent values
    Relp_NKatBorder_i <- unlist(data_gene_df_ed[gene_i,paste(rep("...",length(Relp_NKatBorder_indx)),Relp_NKatBorder_indx,sep="")]) #IF
    Relp_NKatTumor_i <- unlist(data_gene_df_ed[gene_i,paste(rep("...",length(Relp_NKatTumor_indx)),Relp_NKatTumor_indx,sep="")]) #TC
    
    #non-recurrent values
    Nonrelp_NKatBorder_i <- unlist(data_gene_df_ed[gene_i,paste(rep("...",length(Nonrelp_NKatBorder_indx)),Nonrelp_NKatBorder_indx,sep="")])
    Nonrelp_NKatTumor_i <- unlist(data_gene_df_ed[gene_i,paste(rep("...",length(Nonrelp_NKatTumor_indx)),Nonrelp_NKatTumor_indx,sep="")])
    
    x_vec <- c(1,2,3)

    #recurrent y axis
    y_vec_Relp <- c(mean(NKatStroma_i,na.rm=T), mean(Relp_NKatBorder_i,na.rm=T), mean(Relp_NKatTumor_i,na.rm=T))
    y_lw_Relp <- c(min(NKatStroma_i,na.rm=T), min(Relp_NKatBorder_i,na.rm=T), min(Relp_NKatTumor_i,na.rm=T))
    y_up_Relp <- c(max(NKatStroma_i,na.rm=T), max(Relp_NKatBorder_i,na.rm=T), max(Relp_NKatTumor_i,na.rm=T))
    
    #non-recurrent y axis
    y_vec_Nonrelp <- c(mean(NKatStroma_i,na.rm=T), mean(Nonrelp_NKatBorder_i,na.rm=T), mean(Nonrelp_NKatTumor_i,na.rm=T))
    y_lw_Nonrelp <- c(min(NKatStroma_i,na.rm=T), min(Nonrelp_NKatBorder_i,na.rm=T), min(Nonrelp_NKatTumor_i,na.rm=T))
    y_up_Nonrelp <- c(max(NKatStroma_i,na.rm=T), max(Nonrelp_NKatBorder_i,na.rm=T), max(Nonrelp_NKatTumor_i,na.rm=T))

    ylim_min <- min(min(y_lw_Nonrelp,na.rm=T),min(y_lw_Relp,na.rm=T))
    ylim_max <- max(max(y_up_Nonrelp,na.rm=T),max(y_up_Relp,na.rm=T)) * 1.1
    
    #draw boxplots
    boxplot(NKatStroma_i,Relp_NKatBorder_i,Relp_NKatTumor_i, Nonrelp_NKatBorder_i, Nonrelp_NKatTumor_i, main=gene_i,at=c(1,2,3,2,3),names=c("Stroma","Border","Tumor","",""),notch=,
            col=c(adjustcolor("gray8", alpha.f=0.1),adjustcolor("#E73334", alpha.f=0.2),adjustcolor("#E73334", alpha.f=0.2),adjustcolor("#559AC6", alpha.f=0.2),adjustcolor("#559AC6", alpha.f=0.2)),
            outpch = NA, boxwex = 0.5,xlim=c(0.7,3.5),ylim=c(ylim_min,ylim_max))
    stripchart(x=NKatStroma_i,pch=19, cex=0.75, at=c(1), add=T,vertical=T,method = "jitter",col=adjustcolor("gray8", alpha.f=0.5))
    stripchart(x=Relp_NKatBorder_i,pch=19, cex=0.75, at=c(2), add=T,vertical=T,method = "jitter",col=adjustcolor("#E73334", alpha.f=0.8))
    stripchart(x=Relp_NKatTumor_i,pch=19, cex=0.75, at=c(3), add=T,vertical=T,method = "jitter",col=adjustcolor("#E73334", alpha.f=0.8))
    stripchart(x=Nonrelp_NKatBorder_i,pch=19, cex=0.75, at=c(2), add=T,vertical=T,method = "jitter",col=adjustcolor("#559AC6", alpha.f=0.8))
    stripchart(x=Nonrelp_NKatTumor_i,pch=19, cex=0.75, at=c(3), add=T,vertical=T,method = "jitter",col=adjustcolor("#559AC6", alpha.f=0.8))
    lines(x = x_vec, y = y_vec_Relp, col=adjustcolor("#E73334", alpha.f=0.8),xaxs="i",yaxs="i",lwd=6)
    lines(x = x_vec, y = y_vec_Nonrelp, col=adjustcolor("#559AC6", alpha.f=0.8),xaxs="i",yaxs="i",lwd=6)
    legend("top", legend=c(paste(names(type_i),type_i,sep="_",collapse = "_vs_"),paste("FDR border = ",padj_border,sep=""),paste("FDR tumor = ",padj_tumor,sep="")), bty="n",cex=0.8)
    
    temp_max <- c(temp_max,ylim_max); temp_min <- c(temp_min,ylim_min)
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend("topright", inset = 0.01, legend=c("Relapsed","Non-relapsed"), col=c(adjustcolor("#E73334", alpha.f=0.7),adjustcolor("#559AC6", alpha.f=0.7)), lty=c(1,1), bty="n",cex=1.2,lwd=5,xpd = TRUE, horiz = TRUE)
  dev.off()
}
