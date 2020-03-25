#' Pipeline Analysis 
library(MAESTRO)
library(Seurat)
library(ggplot2)
library(Matrix)
library(future)
library(future.apply)
library(pbapply)
library(aricode)
options(bitmapType='cairo')
plan("multiprocess", workers = 8)
options(future.globals.maxSize = 10*1024^3)

# load datasets using different peak-RP models
PBMC1 <- readRDS("integration_RPmodel_benchmark/10X_PBMC_gene_score_TSS_10K_integration.rds")
PBMC2 <- readRDS("integration_RPmodel_benchmark/10X_PBMC_gene_score_TSS_1K_integration.rds")
PBMC3 <- readRDS("integration_RPmodel_benchmark/10X_PBMC_gene_score_GB_10K_integration.rds")

BMMC_PBMC1 <- readRDS("integration_RPmodel_benchmark/10X_BMMC_PBMC_gene_score_TSS_10K_integration.rds")
BMMC_PBMC2 <- readRDS("integration_RPmodel_benchmark/10X_BMMC_PBMC_gene_score_TSS_1K_integration.rds")
BMMC_PBMC3 <- readRDS("integration_RPmodel_benchmark/10X_BMMC_PBMC_gene_score_GB_10K_integration.rds")

BMMC1 <- readRDS("integration_RPmodel_benchmark/10X_BMMC_gene_score_TSS_10K_integration.rds")
BMMC2 <- readRDS("integration_RPmodel_benchmark/10X_BMMC_gene_score_TSS_1K_integration.rds")
BMMC3 <- readRDS("integration_RPmodel_benchmark/10X_BMMC_gene_score_GB_10K_integration.rds")

PBMC1_score <- PBMC1@meta.data[which(PBMC1@meta.data$tech=="ATAC"),"prediction.score.max"]
PBMC2_score <- PBMC2@meta.data[which(PBMC2@meta.data$tech=="ATAC"),"prediction.score.max"]
PBMC3_score <- PBMC3@meta.data[which(PBMC3@meta.data$tech=="ATAC"),"prediction.score.max"]

BMMC_PBMC1_score <- BMMC_PBMC1@meta.data[which(BMMC_PBMC1@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC_PBMC2_score <- BMMC_PBMC2@meta.data[which(BMMC_PBMC2@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC_PBMC3_score <- BMMC_PBMC3@meta.data[which(BMMC_PBMC3@meta.data$tech=="ATAC"),"prediction.score.max"]

BMMC1_score <- BMMC1@meta.data[which(BMMC1@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC2_score <- BMMC2@meta.data[which(BMMC2@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC3_score <- BMMC3@meta.data[which(BMMC3@meta.data$tech=="ATAC"),"prediction.score.max"]

# label prediction score distribution
library(RColorBrewer)
ncol <- brewer.pal(8,"Set1")
HP=0.5
pdf("ATAC_RNA_prediction_score_density.pdf",width=12,height=3.5)
par(mfrow=c(1,3), mar=c(5,5,3,3))
plot(density(PBMC1_score),xlim=c(0,1),lwd=2,col=ncol[1],main="PBMC Different Donors",ylim=c(0,4))
lines(density(PBMC2_score),lwd=2,col=ncol[2])
lines(density(PBMC3_score),lwd=2,col=ncol[3])
abline(v=0.5,lty=2)
legend("topleft",c(paste("TSS 10K mean = ", round(mean(PBMC1_score),3)),
                   paste("TSS 1K mean =", round(mean(PBMC2_score),3)),
                   paste("Genebody 10K mean =", round(mean(PBMC3_score),3))),col=ncol[1:3],lwd=2,box.lty=0,cex=0.7)
legend("left",c(paste("TSS 10K HP = ", length(PBMC1_score[PBMC1_score>=HP])),
                   paste("TSS 1K HP =", length(PBMC2_score[PBMC2_score>HP])),
                   paste("Genebody 10K HP =", length(PBMC3_score[PBMC3_score>HP]))),col=ncol[1:3],lwd=2,box.lty=0,cex=0.7)
plot(density(BMMC_PBMC1_score),xlim=c(0,1),lwd=2,col=ncol[1],main="PBMC Same Donor",ylim=c(0,6))
lines(density(BMMC_PBMC2_score),lwd=2,col=ncol[2])
lines(density(BMMC_PBMC3_score),lwd=2,col=ncol[3])
abline(v=0.5,lty=2)
legend("topleft",c(paste("TSS 10K mean = ", round(mean(BMMC_PBMC1_score),3)),
                   paste("TSS 1K mean =", round(mean(BMMC_PBMC2_score),3)),
                   paste("Genebody 10K mean =", round(mean(BMMC_PBMC3_score),3))),col=ncol[1:3],lwd=2,box.lty=0,cex=0.7)
legend("left",c(paste("TSS 10K HP = ", length(BMMC_PBMC1_score[BMMC_PBMC1_score>=HP])),
                   paste("TSS 1K HP =", length(BMMC_PBMC2_score[BMMC_PBMC2_score>HP])),
                   paste("Genebody 10K HP =", length(BMMC_PBMC3_score[BMMC_PBMC3_score>HP]))),col=ncol[1:3],lwd=2,box.lty=0,cex=0.7)
plot(density(BMMC1_score),xlim=c(0,1),lwd=2,col=ncol[1],main="BMMC Same Donor",ylim=c(0,4))
lines(density(BMMC2_score),lwd=2,col=ncol[2])
lines(density(BMMC3_score),lwd=2,col=ncol[3])
abline(v=0.5,lty=2)
legend("topleft",c(paste("TSS 10K mean = ", round(mean(BMMC1_score),3)),
                   paste("TSS 1K mean =", round(mean(BMMC2_score),3)),
                   paste("Genebody 10K mean =", round(mean(BMMC3_score),3))),col=ncol[1:3],lwd=2,box.lty=0,cex=0.7)
legend("left",c(paste("TSS 10K HP = ", length(BMMC1_score[BMMC1_score>=HP])),
                   paste("TSS 1K HP =", length(BMMC2_score[BMMC2_score>HP])),
                   paste("Genebody 10K HP =", length(BMMC3_score[BMMC3_score>HP]))),col=ncol[1:3],lwd=2,box.lty=0,cex=0.7)
dev.off()

# NMI compare and correlation between RNA and ATAC
correlation_compare <- function(combined, name, method = "pearson", score.filter = 0.0, expr.filter = 0.0)
{
  # split the RNA and ATAC into two object
  combined_RNA <- subset(combined, cells = rownames(combined@meta.data[which(combined@meta.data[,'tech']=='RNA'),]))
  combined_ATAC <- subset(combined, cells = rownames(combined@meta.data[which(combined@meta.data[,'tech']=='ATAC'),]))
  ATAC_clusters <- combined_ATAC@meta.data[which(combined_ATAC@meta.data$prediction.score.max>=score.filter),]$assign.ident
  ATAC_labels <- combined_ATAC@meta.data[which(combined_ATAC@meta.data$prediction.score.max>=score.filter),]$ATAC_snn_res.0.6
  combined_ATAC_NMI <- NULL
  combined_ATAC_NMI <- NMI(ATAC_labels,ATAC_clusters)
  
  # calculate the mean expression/activity at cluster level
  combined_mean <- NULL; combined_mean_name <- NULL
  shared_gene <- intersect(rownames(combined$RNA), rownames(combined$ACTIVITY))
  for(celltype in unique(combined@meta.data$assign.ident))
  {  
    RNA_cells <- rownames(combined@meta.data[which(combined@meta.data[,'tech']=='RNA'&combined@meta.data[,'assign.ident']==celltype),])
    ATAC_cells <- rownames(combined@meta.data[which(combined@meta.data[,'tech']=='ATAC'&combined@meta.data[,'assign.ident']==celltype&combined@meta.data$prediction.score.max>=score.filter),])
    if(length(RNA_cells)>=10&length(ATAC_cells)>=10){
    combined_mean <- cbind(combined_mean, apply(log2(combined$RNA@data[shared_gene,RNA_cells]+1),1,mean))
    combined_mean <- cbind(combined_mean, apply(log2(combined$ACTIVITY@data[shared_gene,ATAC_cells]+1),1,mean))
    combined_mean_name <- c(combined_mean_name, paste0(celltype,"_RNA"), paste0(celltype,"_ATAC"))}
  }
  colnames(combined_mean) <- combined_mean_name
 
  cor_out <- NULL;cor_name <- NULL
  pdf(paste0(name,".pdf"),width=14,height=5.5)
  par(mfrow=c(2,6))
  for(i in seq(1,ncol(combined_mean),2))
  { 
    if((!is.na(mean(combined_mean[,i])))&(!is.na(mean(combined_mean[,i+1]))))
      {tmp <- combined_mean[,c(i,i+1)]
       tmp <- tmp[which(tmp[,1]>expr.filter),]
       smoothScatter(tmp[,1], tmp[,2], xlab="log2(RNA count+1)", ylab="log2(ATAC score+1)",pch='.',col="blue",main=gsub("_RNA","",colnames(combined_mean)[i]), colramp = colorRampPalette(c("blue", "orange", "red")))
       legend("topright", paste0("R = ",round(cor(tmp[,1], tmp[,2], method=method),2)),box.lty=0)
       cor_out <- c(cor_out, round(cor(tmp[,1], tmp[,2], method=method),2))
       cor_name <- c(cor_name, gsub("_RNA","",colnames(combined_mean)[i]))}
  }
  dev.off()
  names(cor_out) <- cor_name
  return(list(cor=cor_out, NMI = combined_ATAC_NMI))
}

PBMC1_result <- correlation_compare(PBMC1, "10X_PBMC_gene_score_TSS_10K_cor")
PBMC2_result <- correlation_compare(PBMC2, "10X_PBMC_gene_score_TSS_1K_cor")
PBMC3_result <- correlation_compare(PBMC3, "10X_PBMC_gene_score_GB_10K_cor")

PBMC1_result$cor['CD4 Naive T'] = NA;PBMC1_result$cor['Macrophage'] = NA;
PBMC2_result$cor['CD4 Naive T'] = NA;PBMC2_result$cor['Macrophage'] = NA;PBMC2_result$cor['pDC'] = NA;
PBMC3_result$cor['CD4 Naive T'] = NA

BMMC_PBMC1_result <- correlation_compare(BMMC_PBMC1, "10X_BMMC_PBMC_gene_score_TSS_10K_cor")
BMMC_PBMC2_result <- correlation_compare(BMMC_PBMC2, "10X_BMMC_PBMC_gene_score_TSS_1K_cor")
BMMC_PBMC3_result <- correlation_compare(BMMC_PBMC3, "10X_BMMC_PBMC_gene_score_GB_10K_cor")

BMMC1_result <- correlation_compare(BMMC1, "10X_BMMC_gene_score_TSS_10K_cor")
BMMC2_result <- correlation_compare(BMMC2, "10X_BMMC_gene_score_TSS_1K_cor")
BMMC3_result <- correlation_compare(BMMC3, "10X_BMMC_gene_score_GB_10K_cor")

pdf("ATAC_RNA_pearson_cor_allcells_compare.pdf",width=12,height=3.5)
par(mfrow=c(1,3),mar=c(8,5,3,3))
barplot(rbind(PBMC1_result$cor[sort(names(PBMC1_result$cor))], PBMC2_result$cor[sort(names(PBMC1_result$cor))], PBMC3_result$cor[sort(names(PBMC1_result$cor))]), 
        beside=T, col=ncol[1:3], border=NA, names=sort(names(PBMC1_result$cor)), legend.text=c("TSS 10K","TSS 1K","Genebody 10K"),main="PBMC Different Donors", las=2, args.legend = list(x="topleft",box.lty=0,ncol=3), ylab="Pearson's correlation between RNA and ATAC", ylim=c(0,0.4))
barplot(rbind(BMMC_PBMC1_result$cor[sort(names(BMMC_PBMC1_result$cor))], BMMC_PBMC2_result$cor[sort(names(BMMC_PBMC1_result$cor))], BMMC_PBMC3_result$cor[sort(names(BMMC_PBMC1_result$cor))]), 
        beside=T, col=ncol[1:3], border=NA, names=sort(names(BMMC_PBMC1_result$cor)), legend.text=c("TSS 10K","TSS 1K","Genebody 10K"),main="PBMC Same Donor", las=2, args.legend = list(x="topleft",box.lty=0,ncol=3), ylab="Pearson's correlation between RNA and ATAC", ylim=c(0,0.4))
barplot(rbind(BMMC1_result$cor[sort(names(BMMC1_result$cor))], BMMC2_result$cor[sort(names(BMMC1_result$cor))], BMMC3_result$cor[sort(names(BMMC1_result$cor))]), 
        beside=T, col=ncol[1:3], border=NA, names=sort(names(BMMC1_result$cor)), legend.text=c("TSS 10K","TSS 1K","Genebody 10K"),main="BMMC Same Donor", las=2, args.legend = list(x="topleft",box.lty=0,ncol=3), ylab="Pearson's correlation between RNA and ATAC", ylim=c(0,0.4))
dev.off()
pdf("ATAC_RNA_NMI_allcells_compare.pdf",width=8,height=3.5)
par(mfrow=c(1,3),mar=c(8,5,3,3))
bp <- barplot(c(PBMC1_result$NMI, PBMC2_result$NMI, PBMC3_result$NMI), beside=T, col=ncol[1:3], border=NA, names=c("TSS 10K","TSS 1K","Genebody 10K"), main="PBMC Different Donors", las=2, args.legend = list(x="topleft",box.lty=0), ylab="NMI of scATAC-seq cluster\nto transferred celltypes", ylim=c(0,1))
text(bp, c(PBMC1_result$NMI, PBMC2_result$NMI, PBMC3_result$NMI)+0.05, round(c(PBMC1_result$NMI, PBMC2_result$NMI, PBMC3_result$NMI),3), cex=0.7)
bp <- barplot(c(BMMC_PBMC1_result$NMI, BMMC_PBMC2_result$NMI, BMMC_PBMC3_result$NMI), beside=T, col=ncol[1:3], border=NA, names=c("TSS 10K","TSS 1K","Genebody 10K"), main="PBMC Same Donor", las=2, args.legend = list(x="topleft",box.lty=0), ylab="NMI of scATAC-seq cluster\nto transferred celltypes", ylim=c(0,1))
text(bp, c(BMMC_PBMC1_result$NMI, BMMC_PBMC2_result$NMI, BMMC_PBMC3_result$NMI)+0.05, round(c(BMMC_PBMC1_result$NMI, BMMC_PBMC2_result$NMI, BMMC_PBMC3_result$NMI),3), cex=0.7)
bp <- barplot(c(BMMC1_result$NMI, BMMC2_result$NMI, BMMC3_result$NMI), beside=T, col=ncol[1:3], border=NA, names=c("TSS 10K","TSS 1K","Genebody 10K"), main="BMMC Same Donor", las=2, args.legend = list(x="topleft",box.lty=0), ylab="NMI of scATAC-seq cluster\nto transferred celltypes", ylim=c(0,1))
text(bp, c(BMMC1_result$NMI, BMMC2_result$NMI, BMMC3_result$NMI)+0.05, round(c(BMMC1_result$NMI, BMMC2_result$NMI, BMMC3_result$NMI),3), cex=0.7)
dev.off()

































