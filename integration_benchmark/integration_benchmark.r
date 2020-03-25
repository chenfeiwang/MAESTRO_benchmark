library(MAESTRO)
library(Seurat)
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(aricode)
ncol <- brewer.pal(8,"Set1")

# load integration files
PBMC_MAESTRO <- readRDS('integration_RPmodel_benchmark/10X_PBMC_gene_score_TSS_10K_integration.rds')
PBMC_snapATAC <- readRDS('integration_benchmark/10X_PBMC_snapATAC_SeuratObj.rds')
PBMC_Seurat <- readRDS('integration_benchmark/10X_PBMC_Seurat_SeuratObj.rds')
PBMC_cicero <- readRDS('integration_benchmark/10X_PBMC_cicero_SeuratObj.rds')

BMMC_PBMC_MAESTRO <- readRDS('integration_RPmodel_benchmark/10X_BMMC_PBMC_gene_score_TSS_10K_integration.rds')
BMMC_PBMC_snapATAC <- readRDS('integration_benchmark/10X_BMMC_PBMC_snapATAC_SeuratObj.rds')
BMMC_PBMC_Seurat <- readRDS('integration_benchmark/10X_BMMC_PBMC_Seurat_SeuratObj.rds')
BMMC_PBMC_cicero <- readRDS('integration_benchmark/10X_BMMC_PBMC_cicero_SeuratObj.rds')

BMMC_MAESTRO <- readRDS('integration_RPmodel_benchmark/10X_BMMC_gene_score_TSS_10K_integration.rds')
BMMC_snapATAC <- readRDS('integration_benchmark/10X_BMMC_Merged_snapATAC_SeuratObj.rds')
BMMC_Seurat <- readRDS('integration_benchmark/10X_BMMC_Merged_Seurat_SeuratObj.rds')
BMMC_cicero <- readRDS('integration_benchmark/10X_BMMC_Merged_cicero_SeuratObj.rds')

# label prediction score distribution
PBMC_MAESTRO_score <- PBMC_MAESTRO@meta.data[which(PBMC_MAESTRO@meta.data$tech=="ATAC"),"prediction.score.max"]
PBMC_Seurat_score <- PBMC_Seurat@meta.data[which(PBMC_Seurat@meta.data$tech=="ATAC"),"prediction.score.max"]
PBMC_snapATAC_score <- PBMC_snapATAC@meta.data[which(PBMC_snapATAC@meta.data$tech=="ATAC"),"prediction.score.max"]
PBMC_cicero_score <- PBMC_cicero@meta.data[which(PBMC_cicero@meta.data$tech=="ATAC"),"prediction.score.max"]

BMMC_PBMC_MAESTRO_score <- BMMC_PBMC_MAESTRO@meta.data[which(BMMC_PBMC_MAESTRO@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC_PBMC_Seurat_score <- BMMC_PBMC_Seurat@meta.data[which(BMMC_PBMC_Seurat@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC_PBMC_snapATAC_score <- BMMC_PBMC_snapATAC@meta.data[which(BMMC_PBMC_snapATAC@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC_PBMC_cicero_score <- BMMC_PBMC_cicero@meta.data[which(BMMC_PBMC_cicero@meta.data$tech=="ATAC"),"prediction.score.max"]

BMMC_MAESTRO_score <- BMMC_MAESTRO@meta.data[which(BMMC_MAESTRO@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC_Seurat_score <- BMMC_Seurat@meta.data[which(BMMC_Seurat@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC_snapATAC_score <- BMMC_snapATAC@meta.data[which(BMMC_snapATAC@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC_cicero_score <- BMMC_cicero@meta.data[which(BMMC_cicero@meta.data$tech=="ATAC"),"prediction.score.max"]

HP = 0.5
pdf("ATAC_RNA_prediction_score_density.pdf",width=12,height=3.5)
par(mfrow=c(1,3),mar=c(5,5,3,3))
plot(density(PBMC_MAESTRO_score),xlim=c(0,1),lwd=2,col=ncol[1],main="PBMC Different Donors",ylim=c(0,4))
lines(density(PBMC_Seurat_score),lwd=2,col=ncol[2])
lines(density(PBMC_snapATAC_score),lwd=2,col=ncol[3])
lines(density(PBMC_cicero_score),lwd=2,col=ncol[4])
legend("topleft",c(paste("MAESTRO mean = ", round(mean(PBMC_MAESTRO_score),3)),
                   paste("Seurat mean =", round(mean(PBMC_Seurat_score),3)),
                   paste("snapATAC mean =", round(mean(PBMC_snapATAC_score),3)),
                   paste("cicero mean =",round(mean(PBMC_cicero_score),3))),col=ncol[1:4],lwd=2,box.lty=0,cex=0.8)
legend("left",c(paste("MAESTRO HP = ", length(PBMC_MAESTRO_score[PBMC_MAESTRO_score>=HP])),
                   paste("Seurat HP =", length(PBMC_Seurat_score[PBMC_Seurat_score>HP])),
                   paste("snapATAC HP =", length(PBMC_snapATAC_score[PBMC_snapATAC_score>HP])),
                   paste("cicero HP =",length(PBMC_cicero_score[PBMC_cicero_score>HP]))),col=ncol[1:4],lwd=2,box.lty=0,cex=0.8)
abline(v=HP,lty=2)
plot(density(BMMC_PBMC_MAESTRO_score),xlim=c(0,1),lwd=2,col=ncol[1],main="PBMC Same Donor",ylim=c(0,6))
lines(density(BMMC_PBMC_Seurat_score),lwd=2,col=ncol[2])
lines(density(BMMC_PBMC_snapATAC_score),lwd=2,col=ncol[3])
lines(density(BMMC_PBMC_cicero_score),lwd=2,col=ncol[4])
legend("topleft",c(paste("MAESTRO mean =", round(mean(BMMC_PBMC_MAESTRO_score),3)),
                   paste("Seurat mean =", round(mean(BMMC_PBMC_Seurat_score),3)),
                   paste("snapATAC mean =", round(mean(BMMC_PBMC_snapATAC_score),3)),
                   paste("cicero mean =", round(mean(BMMC_PBMC_cicero_score),3))),col=ncol[1:4],lwd=2,box.lty=0,cex=0.8)
legend("left",c(paste("MAESTRO HP = ", length(BMMC_PBMC_MAESTRO_score[BMMC_PBMC_MAESTRO_score>HP])),
                   paste("Seurat HP =", length(BMMC_PBMC_Seurat_score[BMMC_PBMC_Seurat_score>HP])),
                   paste("snapATAC HP =", length(BMMC_PBMC_snapATAC_score[BMMC_PBMC_snapATAC_score>HP])),
                   paste("cicero HP =",length(BMMC_PBMC_cicero_score[BMMC_PBMC_cicero_score>HP]))),col=ncol[1:4],lwd=2,box.lty=0,cex=0.8)
abline(v=HP,lty=2)
plot(density(BMMC_MAESTRO_score),xlim=c(0,1),lwd=2,col=ncol[1],main="BMMC Same Donor",ylim=c(0,4))
lines(density(BMMC_Seurat_score),lwd=2,col=ncol[2])
lines(density(BMMC_snapATAC_score),lwd=2,col=ncol[3])
lines(density(BMMC_cicero_score),lwd=2,col=ncol[4])
legend("topleft",c(paste("MAESTRO mean =", round(mean(BMMC_MAESTRO_score),3)),
                   paste("Seurat mean =", round(mean(BMMC_Seurat_score),3)),
                   paste("snapATAC mean =", round(mean(BMMC_snapATAC_score),3)),
                   paste("cicero mean =", round(mean(BMMC_cicero_score),3))),col=ncol[1:4],lwd=2,box.lty=0,cex=0.8)
legend("left",c(paste("MAESTRO HP = ", length(BMMC_MAESTRO_score[BMMC_MAESTRO_score>HP])),
                   paste("Seurat HP =", length(BMMC_Seurat_score[BMMC_Seurat_score>HP])),
                   paste("snapATAC HP =", length(BMMC_snapATAC_score[BMMC_snapATAC_score>HP])),
                   paste("cicero HP =",length(BMMC_cicero_score[BMMC_cicero_score>HP]))),col=ncol[1:4],lwd=2,box.lty=0,cex=0.8)
abline(v=HP,lty=2)
dev.off()

# count the number of cells that have successfully recieved labels from scRNA-seq
table(BMMC_MAESTRO@meta.data[which(BMMC_MAESTRO@meta.data$tech=="RNA"),"assign.ident"])
table(BMMC_MAESTRO@meta.data[which(BMMC_MAESTRO@meta.data$tech=="ATAC"),"assign.ident"])
table(BMMC_Seurat@meta.data[which(BMMC_Seurat@meta.data$tech=="ATAC"),"assign.ident"])
table(BMMC_snapATAC@meta.data[which(BMMC_snapATAC@meta.data$tech=="ATAC"),"assign.ident"])
table(BMMC_cicero@meta.data[which(BMMC_cicero@meta.data$tech=="ATAC"),"assign.ident"])

correlation_compare <- function(combined, name, method = "pearson", score.filter = 0.0, expr.filter = 0.5)
{
  # split the RNA and ATAC into two object
  combined_RNA <- subset(combined, cells = rownames(combined@meta.data[which(combined@meta.data[,'tech']=='RNA'),]))
  combined_ATAC <- subset(combined, cells = rownames(combined@meta.data[which(combined@meta.data[,'tech']=='ATAC'),]))
  ATAC_clusters <- combined_ATAC@meta.data[which(combined_ATAC@meta.data$prediction.score.max>=score.filter),]$assign.ident
  ATAC_labels <- combined_ATAC@meta.data[which(combined_ATAC@meta.data$prediction.score.max>=score.filter),]$seurat_clusters
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

# NMI compare, correlation between RNA and ATAC
PBMC_MAESTRO_result <- correlation_compare(PBMC_MAESTRO, "10X_PBMC_MAESTRO")
PBMC_snapATAC_result <- correlation_compare(PBMC_snapATAC,  "10X_PBMC_snapATAC")
PBMC_Seurat_result <- correlation_compare(PBMC_Seurat, "10X_PBMC_Seurat")
PBMC_cicero_result <- correlation_compare(PBMC_cicero, "10X_PBMC_cicero")
PBMC_snapATAC_result$cor['CD4 Naive T'] = NA;PBMC_snapATAC_result$cor['pDC'] = NA;PBMC_Seurat_result$cor['CD4 Naive T'] = NA;
PBMC_cicero_result$cor['CD4 Naive T'] = NA;PBMC_cicero_result$cor['DC'] = NA;PBMC_cicero_result$cor['pDC'] = NA

BMMC_PBMC_MAESTRO_result <- correlation_compare(BMMC_PBMC_MAESTRO, "10X_BMMC_PBMC_MAESTRO")
BMMC_PBMC_snapATAC_result <- correlation_compare(BMMC_PBMC_snapATAC,  "10X_BMMC_PBMC_snapATAC")
BMMC_PBMC_Seurat_result <- correlation_compare(BMMC_PBMC_Seurat, "10X_BMMC_PBMC_Seurat")
BMMC_PBMC_cicero_result <- correlation_compare(BMMC_PBMC_cicero, "10X_BMMC_PBMC_cicero")
BMMC_PBMC_cicero_result$cor['CD16 Mono'] = NA;BMMC_PBMC_cicero_result$cor['DC'] = NA;BMMC_PBMC_cicero_result$cor['pDC'] = NA

BMMC_MAESTRO_result <- correlation_compare(BMMC_MAESTRO, "10X_BMMC_Merged_MAESTRO")
BMMC_snapATAC_result <- correlation_compare(BMMC_snapATAC, "10X_BMMC_Merged_snapATAC")
BMMC_Seurat_result <- correlation_compare(BMMC_Seurat, "10X_BMMC_Merged_Seurat")
BMMC_cicero_result <- correlation_compare(BMMC_cicero, "10X_BMMC_Merged_cicero")

pdf("ATAC_RNA_pearson_cor_high_allcells_compare.pdf",width=12,height=3.5)
par(mfrow=c(1,3),mar=c(8,5,3,3))
barplot(rbind(PBMC_MAESTRO_result$cor[sort(names(PBMC_MAESTRO_result$cor))], PBMC_Seurat_result$cor[sort(names(PBMC_MAESTRO_result$cor))], PBMC_snapATAC_result$cor[sort(names(PBMC_MAESTRO_result$cor))], PBMC_cicero_result$cor[sort(names(PBMC_MAESTRO_result$cor))]), 
        beside=T, col=ncol[1:4], border=NA, names=sort(names(PBMC_MAESTRO_result$cor)), legend.text=c("MAESTRO","Seurat","snapATAC","cicero"),main="PBMC Different Donors", las=2, args.legend = list(x="topleft",box.lty=0), ylab="Spearman's correlation between RNA and ATAC", ylim=c(-0.4,0.4))
barplot(rbind(BMMC_PBMC_MAESTRO_result$cor[sort(names(BMMC_PBMC_MAESTRO_result$cor))], BMMC_PBMC_Seurat_result$cor[sort(names(BMMC_PBMC_MAESTRO_result$cor))], BMMC_PBMC_snapATAC_result$cor[sort(names(BMMC_PBMC_MAESTRO_result$cor))], BMMC_PBMC_cicero_result$cor[sort(names(BMMC_PBMC_MAESTRO_result$cor))]),
        beside=T, col=ncol[1:4], border=NA, names=sort(names(BMMC_PBMC_MAESTRO_result$cor)), legend.text=c("MAESTRO","Seurat","snapATAC","cicero"),main="PBMC Same Donor", las=2, args.legend = list(x="topleft",box.lty=0), ylab="Spearman's correlation between RNA and ATAC", ylim=c(-0.4,0.4))
barplot(rbind(BMMC_MAESTRO_result$cor[sort(names(BMMC_MAESTRO_result$cor))], BMMC_Seurat_result$cor[sort(names(BMMC_MAESTRO_result$cor))], BMMC_snapATAC_result$cor[sort(names(BMMC_MAESTRO_result$cor))], BMMC_cicero_result$cor[sort(names(BMMC_MAESTRO_result$cor))]),
        beside=T, col=ncol[1:4], border=NA, names=sort(names(BMMC_MAESTRO_result$cor)), legend.text=c("MAESTRO","Seurat","snapATAC","cicero"),main="BMMC Same Donor", las=2, args.legend = list(x="topleft",box.lty=0), ylab="Spearman's correlation between RNA and ATAC", ylim=c(-0.4,0.4))
dev.off()
pdf("ATAC_RNA_NMI_allcells_compare.pdf",width=8,height=4.5)
par(mfrow=c(1,3),mar=c(8,5,3,3))
bp <- barplot(c(PBMC_MAESTRO_result$NMI, PBMC_Seurat_result$NMI, PBMC_snapATAC_result$NMI, PBMC_cicero_result$NMI), beside=T, col=ncol[1:4], border=NA, names=c("MAESTRO","Seurat","snapATAC","cicero"), main="PBMC Different Donors", las=2, args.legend = list(x="topleft",box.lty=0), ylab="NMI of scATAC-seq cluster\nto transferred celltypes", ylim=c(0,1))
text(bp, c(PBMC_MAESTRO_result$NMI, PBMC_Seurat_result$NMI, PBMC_snapATAC_result$NMI, PBMC_cicero_result$NMI)+0.05, round(c(PBMC_MAESTRO_result$NMI, PBMC_Seurat_result$NMI, PBMC_snapATAC_result$NMI, PBMC_cicero_result$NMI),3),cex=0.85)
bp <- barplot(c(BMMC_PBMC_MAESTRO_result$NMI, BMMC_PBMC_Seurat_result$NMI, BMMC_PBMC_snapATAC_result$NMI, BMMC_PBMC_cicero_result$NMI), beside=T, col=ncol[1:4], border=NA, names=c("MAESTRO","Seurat","snapATAC","cicero"), main="PBMC Same Donor", las=2, args.legend = list(x="topleft",box.lty=0), ylab="NMI of scATAC-seq cluster\nto transferred celltypes", ylim=c(0,1))
text(bp, c(BMMC_PBMC_MAESTRO_result$NMI, BMMC_PBMC_Seurat_result$NMI, BMMC_PBMC_snapATAC_result$NMI, BMMC_PBMC_cicero_result$NMI)+0.05, round(c(BMMC_PBMC_MAESTRO_result$NMI, BMMC_PBMC_Seurat_result$NMI, BMMC_PBMC_snapATAC_result$NMI, BMMC_PBMC_cicero_result$NMI),3),cex=0.85)
bp <- barplot(c(BMMC_MAESTRO_result$NMI, BMMC_Seurat_result$NMI, BMMC_snapATAC_result$NMI, BMMC_cicero_result$NMI), beside=T, col=ncol[1:4], border=NA, names=c("MAESTRO","Seurat","snapATAC","cicero"), main="BMMC Same Donor", las=2, args.legend = list(x="topleft",box.lty=0), ylab="NMI of scATAC-seq cluster\nto transferred celltypes", ylim=c(0,1))
text(bp, c(BMMC_MAESTRO_result$NMI, BMMC_Seurat_result$NMI, BMMC_snapATAC_result$NMI, BMMC_cicero_result$NMI)+0.05, round(c(BMMC_MAESTRO_result$NMI, BMMC_Seurat_result$NMI, BMMC_snapATAC_result$NMI, BMMC_cicero_result$NMI),3),cex=0.85)
dev.off()




















