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
library(RColorBrewer)
library(aricode)
ncol <- brewer.pal(9,"Set1")

# function declare
correlation_compare <- function(combined, name, method = "spearman", score.filter = 0.0, expr.filter = 0.0)
{
  # split the RNA and ATAC into two object
  combined_RNA <- subset(combined, cells = rownames(combined@meta.data[which(combined@meta.data[,'tech']=='RNA'),]))
  combined_ATAC <- subset(combined, cells = rownames(combined@meta.data[which(combined@meta.data[,'tech']=='ATAC'),]))
  ATAC_clusters <- combined_ATAC@meta.data[which(combined_ATAC@meta.data$prediction.score.max>=score.filter),]$assign.ident
  ATAC_labels <- combined_ATAC@meta.data[which(combined_ATAC@meta.data$prediction.score.max>=score.filter),]$seurat_clusters
  combined_ATAC_NMI <- NMI(ATAC_labels,ATAC_clusters)
  
  # calculate the mean expression/activity at cluster level
  combined_mean <- NULL; combined_mean_name <- NULL
#  shared_gene <- intersect(rownames(combined$RNA), rownames(combined$ACTIVITY))        # all shared genes
  shared_gene <- intersect(combined$RunPCA.RNA$features, rownames(combined$ACTIVITY))   # only focused on HVGs
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
       legend("topright", paste0("R = ",round(cor(tmp[,1], tmp[,2], method=method),2),", p = ", round(cor.test(tmp[,1], tmp[,2], method=method)$p.value,2)),box.lty=0,text.col="white")
       cor_out <- c(cor_out, round(cor(tmp[,1], tmp[,2], method=method),2))
       cor_name <- c(cor_name, gsub("_RNA","",colnames(combined_mean)[i]))}
  }
  dev.off()
  names(cor_out) <- cor_name
  return(list(cor=cor_out, NMI = combined_ATAC_NMI))
}

Boxplot_List <- function(plot_list,MAIN,NAMES,YLAB,YLIM, COL){
  n <- length(plot_list)
  xpos <- 0:(n-1)+1.5
  p_value <- c()
  for (i in seq(n-1)){
    tryCatch({p_value <- c(p_value,wilcox.test(plot_list[[i]],plot_list[[i+1]],alternative="greater")$p.value)},
      error = function(err){p_value <- c(p_value,1)})
  }
  mark <- symnum(p_value, cutpoints=c(0,0.001,0.01,0.05,1), symbols=c("***","**","*","-"))
  bp <- boxplot(plot_list,at=0:(n-1)+1, xlim=c(0.5,n+0.5),outline=F,plot=F)
  boxplot(plot_list,at=0:(n-1)+1, xlim=c(0.5,n+0.5),ylim=YLIM,border=COL,outline=F,names=NAMES,main=MAIN,lwd=2,ylab=YLAB,las=2)
  box(lwd=2)
}

# load datasets using different peak-RP models
PBMC1 <- readRDS("integration_RPmodel_benchmark/PBMC/10X_PBMC_gene_score_TSS_10K_integration.rds")
# PBMC2 <- readRDS("integration_RPmodel_benchmark/PBMC/10X_PBMC_gene_score_TSS_1K_integration.rds")
PBMC3 <- readRDS("integration_RPmodel_benchmark/PBMC/10X_PBMC_gene_score_GB_10K_integration.rds")
PBMC4 <- readRDS("integration_RPmodel_benchmark/PBMC/10X_PBMC_gene_score_RMP_integration.rds")
PBMC5 <- readRDS("integration_RPmodel_benchmark/PBMC/10X_PBMC_gene_score_RMPE_integration.rds")
PBMC6 <- readRDS("integration_RPmodel_benchmark/PBMC/10X_PBMC_gene_score_RMPG_integration.rds")
PBMC7 <- readRDS("integration_RPmodel_benchmark/PBMC/10X_PBMC_gene_score_EN_integration.rds")
PBMC8 <- readRDS("integration_RPmodel_benchmark/PBMC/10X_PBMC_gene_score_GBN_integration.rds")
PBMC9 <- readRDS("integration_RPmodel_benchmark/PBMC/10X_PBMC_gene_score_EN_RMPE_integration.rds")
PBMC10 <- readRDS("integration_RPmodel_benchmark/PBMC/10X_PBMC_gene_score_GBN_RMPG_integration.rds")

BMMC_PBMC1 <- readRDS("integration_RPmodel_benchmark/BMMC_PBMC/10X_BMMC_PBMC_gene_score_TSS_10K_integration.rds")
# BMMC_PBMC2 <- readRDS("integration_RPmodel_benchmark/BMMC_PBMC/10X_BMMC_PBMC_gene_score_TSS_1K_integration.rds")
BMMC_PBMC3 <- readRDS("integration_RPmodel_benchmark/BMMC_PBMC/10X_BMMC_PBMC_gene_score_GB_10K_integration.rds")
BMMC_PBMC4 <- readRDS("integration_RPmodel_benchmark/BMMC_PBMC/10X_BMMC_PBMC_gene_score_RMP_integration.rds")
BMMC_PBMC5 <- readRDS("integration_RPmodel_benchmark/BMMC_PBMC/10X_BMMC_PBMC_gene_score_RMPE_integration.rds")
BMMC_PBMC6 <- readRDS("integration_RPmodel_benchmark/BMMC_PBMC/10X_BMMC_PBMC_gene_score_RMPG_integration.rds")
BMMC_PBMC7 <- readRDS("integration_RPmodel_benchmark/BMMC_PBMC/10X_BMMC_PBMC_gene_score_EN_integration.rds")
BMMC_PBMC8 <- readRDS("integration_RPmodel_benchmark/BMMC_PBMC/10X_BMMC_PBMC_gene_score_GBN_integration.rds")
BMMC_PBMC9 <- readRDS("integration_RPmodel_benchmark/BMMC_PBMC/10X_BMMC_PBMC_gene_score_EN_RMPE_integration.rds")
BMMC_PBMC10 <- readRDS("integration_RPmodel_benchmark/BMMC_PBMC/10X_BMMC_PBMC_gene_score_GBN_RMPG_integration.rds")

BMMC1 <- readRDS("integration_RPmodel_benchmark/BMMC/10X_BMMC_gene_score_TSS_10K_integration.rds")
# BMMC2 <- readRDS("integration_RPmodel_benchmark/BMMC/10X_BMMC_gene_score_TSS_1K_integration.rds")
BMMC3 <- readRDS("integration_RPmodel_benchmark/BMMC/10X_BMMC_gene_score_GB_10K_integration.rds")
BMMC4 <- readRDS("integration_RPmodel_benchmark/BMMC/10X_BMMC_gene_score_RMP_integration.rds")
BMMC5 <- readRDS("integration_RPmodel_benchmark/BMMC/10X_BMMC_gene_score_RMPE_integration.rds")
BMMC6 <- readRDS("integration_RPmodel_benchmark/BMMC/10X_BMMC_gene_score_RMPG_integration.rds")
BMMC7 <- readRDS("integration_RPmodel_benchmark/BMMC/10X_BMMC_gene_score_EN_integration.rds")
BMMC8 <- readRDS("integration_RPmodel_benchmark/BMMC/10X_BMMC_gene_score_GBN_integration.rds")
BMMC9 <- readRDS("integration_RPmodel_benchmark/BMMC/10X_BMMC_gene_score_EN_RMPE_integration.rds")
BMMC10 <- readRDS("integration_RPmodel_benchmark/BMMC/10X_BMMC_gene_score_GBN_RMPG_integration.rds")

# label prediction score distribution
PBMC1_score <- PBMC1@meta.data[which(PBMC1@meta.data$tech=="ATAC"),"prediction.score.max"]
# PBMC2_score <- PBMC2@meta.data[which(PBMC2@meta.data$tech=="ATAC"),"prediction.score.max"]
PBMC3_score <- PBMC3@meta.data[which(PBMC3@meta.data$tech=="ATAC"),"prediction.score.max"]
PBMC4_score <- PBMC4@meta.data[which(PBMC4@meta.data$tech=="ATAC"),"prediction.score.max"]
PBMC5_score <- PBMC5@meta.data[which(PBMC5@meta.data$tech=="ATAC"),"prediction.score.max"]
PBMC6_score <- PBMC6@meta.data[which(PBMC6@meta.data$tech=="ATAC"),"prediction.score.max"]
PBMC7_score <- PBMC7@meta.data[which(PBMC7@meta.data$tech=="ATAC"),"prediction.score.max"]
PBMC8_score <- PBMC8@meta.data[which(PBMC8@meta.data$tech=="ATAC"),"prediction.score.max"]
PBMC9_score <- PBMC9@meta.data[which(PBMC9@meta.data$tech=="ATAC"),"prediction.score.max"]
PBMC10_score <- PBMC10@meta.data[which(PBMC10@meta.data$tech=="ATAC"),"prediction.score.max"]

BMMC_PBMC1_score <- BMMC_PBMC1@meta.data[which(BMMC_PBMC1@meta.data$tech=="ATAC"),"prediction.score.max"]
# BMMC_PBMC2_score <- BMMC_PBMC2@meta.data[which(BMMC_PBMC2@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC_PBMC3_score <- BMMC_PBMC3@meta.data[which(BMMC_PBMC3@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC_PBMC4_score <- BMMC_PBMC4@meta.data[which(BMMC_PBMC4@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC_PBMC5_score <- BMMC_PBMC5@meta.data[which(BMMC_PBMC5@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC_PBMC6_score <- BMMC_PBMC6@meta.data[which(BMMC_PBMC6@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC_PBMC7_score <- BMMC_PBMC7@meta.data[which(BMMC_PBMC7@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC_PBMC8_score <- BMMC_PBMC8@meta.data[which(BMMC_PBMC8@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC_PBMC9_score <- BMMC_PBMC9@meta.data[which(BMMC_PBMC9@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC_PBMC10_score <- BMMC_PBMC10@meta.data[which(BMMC_PBMC10@meta.data$tech=="ATAC"),"prediction.score.max"]

BMMC1_score <- BMMC1@meta.data[which(BMMC1@meta.data$tech=="ATAC"),"prediction.score.max"]
# BMMC2_score <- BMMC2@meta.data[which(BMMC2@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC3_score <- BMMC3@meta.data[which(BMMC3@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC4_score <- BMMC4@meta.data[which(BMMC4@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC5_score <- BMMC5@meta.data[which(BMMC5@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC6_score <- BMMC6@meta.data[which(BMMC6@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC7_score <- BMMC7@meta.data[which(BMMC7@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC8_score <- BMMC8@meta.data[which(BMMC8@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC9_score <- BMMC9@meta.data[which(BMMC9@meta.data$tech=="ATAC"),"prediction.score.max"]
BMMC10_score <- BMMC10@meta.data[which(BMMC10@meta.data$tech=="ATAC"),"prediction.score.max"]

pdf("ATAC_RNA_prediction_score_boxplot.pdf",width=10,height=3.5)
par(mfrow=c(1,3),mar=c(10,5,3,3))
Boxplot_List(list(PBMC1_score, PBMC3_score, PBMC4_score, PBMC5_score, PBMC6_score, PBMC7_score, PBMC8_score,PBMC9_score, PBMC10_score), "PBMC Different Donors",
           c("P","G","P & Remove adjacent P","P & Remove adjacent P+E","P & Remove adjacent P+G","P & Add normalized E","P & Add normalized G","Remove adjacent P+E & \nAdd normalized E","Remove adjacent P+G & \nAdd normalized G"), "Max celltype label prediction score", c(0.3,1), ncol[1:9])
Boxplot_List(list(BMMC_PBMC1_score, BMMC_PBMC3_score, BMMC_PBMC4_score, BMMC_PBMC5_score, BMMC_PBMC6_score, BMMC_PBMC7_score, BMMC_PBMC8_score,BMMC_PBMC9_score, BMMC_PBMC10_score), "PBMC Same Donors",
           c("P","G","P & Remove adjacent P","P & Remove adjacent P+E","P & Remove adjacent P+G","P & Add normalized E","P & Add normalized G","Remove adjacent P+E & \nAdd normalized E","Remove adjacent P+G & \nAdd normalized G"), "Max celltype label prediction score", c(0.3,1), ncol[1:9])
Boxplot_List(list(BMMC1_score, BMMC3_score, BMMC4_score, BMMC5_score, BMMC6_score, BMMC7_score, BMMC8_score, BMMC9_score, BMMC10_score), "BMMC Same Donors",
           c("P","G","P & Remove adjacent P","P & Remove adjacent P+E","P & Remove adjacent P+G","P & Add normalized E","P & Add normalized G","Remove adjacent P+E & \nAdd normalized E","Remove adjacent P+G & \nAdd normalized G"), "Max celltype label prediction score", c(0.3,1), ncol[1:9])
dev.off()

# Spearman correlation between RNA and ATAC
PBMC1_result <- correlation_compare(PBMC1, "10X_PBMC_gene_score_TSS_10K_cor")
# PBMC2_result <- correlation_compare(PBMC2, "10X_PBMC_gene_score_TSS_1K_cor")
PBMC3_result <- correlation_compare(PBMC3, "10X_PBMC_gene_score_GB_10K_cor")
PBMC4_result <- correlation_compare(PBMC4, "10X_PBMC_gene_score_RMP_cor")
PBMC5_result <- correlation_compare(PBMC5, "10X_PBMC_gene_score_RMPE_cor")
PBMC6_result <- correlation_compare(PBMC6, "10X_PBMC_gene_score_RMPG_cor")
PBMC7_result <- correlation_compare(PBMC7, "10X_PBMC_gene_score_EN_cor")
PBMC8_result <- correlation_compare(PBMC8, "10X_PBMC_gene_score_GBN_cor")
PBMC9_result <- correlation_compare(PBMC9, "10X_PBMC_gene_score_EN_RMPE_cor")
PBMC10_result <- correlation_compare(PBMC10, "10X_PBMC_gene_score_GBN_RMPG_cor")

BMMC_PBMC1_result <- correlation_compare(BMMC_PBMC1, "10X_BMMC_PBMC_gene_score_TSS_10K_cor")
# BMMC_PBMC2_result <- correlation_compare(BMMC_PBMC2, "10X_BMMC_PBMC_gene_score_TSS_1K_cor")
BMMC_PBMC3_result <- correlation_compare(BMMC_PBMC3, "10X_BMMC_PBMC_gene_score_GB_10K_cor")
BMMC_PBMC4_result <- correlation_compare(BMMC_PBMC4, "10X_BMMC_PBMC_gene_score_RMP_cor")
BMMC_PBMC5_result <- correlation_compare(BMMC_PBMC5, "10X_BMMC_PBMC_gene_score_RMPE_cor")
BMMC_PBMC6_result <- correlation_compare(BMMC_PBMC6, "10X_BMMC_PBMC_gene_score_RMPG_cor")
BMMC_PBMC7_result <- correlation_compare(BMMC_PBMC7, "10X_BMMC_PBMC_gene_score_EN_cor")
BMMC_PBMC8_result <- correlation_compare(BMMC_PBMC8, "10X_BMMC_PBMC_gene_score_GBN_cor")
BMMC_PBMC9_result <- correlation_compare(BMMC_PBMC9, "10X_BMMC_PBMC_gene_score_EN_RMPE_cor")
BMMC_PBMC10_result <- correlation_compare(BMMC_PBMC10, "10X_BMMC_PBMC_gene_score_GBN_RMPG_cor")

BMMC1_result <- correlation_compare(BMMC1, "10X_BMMC_gene_score_TSS_10K_cor")
# BMMC2_result <- correlation_compare(BMMC2, "10X_BMMC_gene_score_TSS_1K_cor")
BMMC3_result <- correlation_compare(BMMC3, "10X_BMMC_gene_score_GB_10K_cor")
BMMC4_result <- correlation_compare(BMMC4, "10X_BMMC_gene_score_RMP_cor")
BMMC5_result <- correlation_compare(BMMC5, "10X_BMMC_gene_score_RMPE_cor")
BMMC6_result <- correlation_compare(BMMC6, "10X_BMMC_gene_score_RMPG_cor")
BMMC7_result <- correlation_compare(BMMC7, "10X_BMMC_gene_score_EN_cor")
BMMC8_result <- correlation_compare(BMMC8, "10X_BMMC_gene_score_GBN_cor")
BMMC9_result <- correlation_compare(BMMC9, "10X_BMMC_gene_score_EN_RMPE_cor")
BMMC10_result <- correlation_compare(BMMC10, "10X_BMMC_gene_score_GBN_RMPG_cor")

PBMC_cor <- rbind(PBMC1_result$cor[sort(names(PBMC1_result$cor))], PBMC3_result$cor[sort(names(PBMC1_result$cor))], PBMC4_result$cor[sort(names(PBMC1_result$cor))], 
                  PBMC5_result$cor[sort(names(PBMC1_result$cor))], PBMC6_result$cor[sort(names(PBMC1_result$cor))], PBMC7_result$cor[sort(names(PBMC1_result$cor))], PBMC8_result$cor[sort(names(PBMC1_result$cor))], 
                  PBMC9_result$cor[sort(names(PBMC1_result$cor))], PBMC10_result$cor[sort(names(PBMC1_result$cor))])
rownames(PBMC_cor) <- c("P","G","P & Remove adjacent P","P & Remove adjacent P+E","P & Remove adjacent P+G","P & Add normalized E","P & Add normalized G","Remove adjacent P+E & \nAdd normalized E","Remove adjacent P+G & \nAdd normalized G")
BMMC_PBMC_cor <- rbind(BMMC_PBMC1_result$cor[sort(names(BMMC_PBMC1_result$cor))], BMMC_PBMC3_result$cor[sort(names(BMMC_PBMC1_result$cor))], BMMC_PBMC4_result$cor[sort(names(BMMC_PBMC1_result$cor))], 
                       BMMC_PBMC5_result$cor[sort(names(BMMC_PBMC1_result$cor))], BMMC_PBMC6_result$cor[sort(names(BMMC_PBMC1_result$cor))], BMMC_PBMC7_result$cor[sort(names(BMMC_PBMC1_result$cor))], BMMC_PBMC8_result$cor[sort(names(BMMC_PBMC1_result$cor))], 
                       BMMC_PBMC9_result$cor[sort(names(BMMC_PBMC1_result$cor))], BMMC_PBMC10_result$cor[sort(names(BMMC_PBMC1_result$cor))])
rownames(BMMC_PBMC_cor) <- c("P","G","P & Remove adjacent P","P & Remove adjacent P+E","P & Remove adjacent P+G","P & Add normalized E","P & Add normalized G","Remove adjacent P+E & \nAdd normalized E","Remove adjacent P+G & \nAdd normalized G")
BMMC_cor <- rbind(BMMC1_result$cor[sort(names(BMMC1_result$cor))], BMMC3_result$cor[sort(names(BMMC1_result$cor))], BMMC4_result$cor[sort(names(BMMC1_result$cor))], 
                  BMMC5_result$cor[sort(names(BMMC1_result$cor))], BMMC6_result$cor[sort(names(BMMC1_result$cor))], BMMC7_result$cor[sort(names(BMMC1_result$cor))], BMMC8_result$cor[sort(names(BMMC1_result$cor))], 
                  BMMC9_result$cor[sort(names(BMMC1_result$cor))], BMMC10_result$cor[sort(names(BMMC1_result$cor))])
rownames(BMMC_cor) <- c("P","G","P & Remove adjacent P","P & Remove adjacent P+E","P & Remove adjacent P+G","P & Add normalized E","P & Add normalized G","Remove adjacent P+E & \nAdd normalized E","Remove adjacent P+G & \nAdd normalized G")

pdf("ATAC_RNA_spearman_cor_var_compare_boxplot.pdf",width=10,height=4)
par(mfrow=c(1,3),mar=c(12,5,3,3))
Boxplot_List(list(PBMC_cor[1,], PBMC_cor[2,], PBMC_cor[3,], PBMC_cor[4,], PBMC_cor[5,], PBMC_cor[6,], PBMC_cor[7,],PBMC_cor[8,], PBMC_cor[9,]), "PBMC Different Donors",
           c("P","G","P & Remove adjacent P","P & Remove adjacent P+E","P & Remove adjacent P+G","P & Add normalized E","P & Add normalized G","Remove adjacent P+E & \nAdd normalized E","Remove adjacent P+G & \nAdd normalized G"), "Spearman's correlation between RNA and ATAC", c(0.2,0.65), ncol[1:9])
Boxplot_List(list(BMMC_PBMC_cor[1,], BMMC_PBMC_cor[2,], BMMC_PBMC_cor[3,], BMMC_PBMC_cor[4,], BMMC_PBMC_cor[5,], BMMC_PBMC_cor[6,], BMMC_PBMC_cor[7,],BMMC_PBMC_cor[8,], BMMC_PBMC_cor[9,]), "PBMC Same Donors",
           c("P","G","P & Remove adjacent P","P & Remove adjacent P+E","P & Remove adjacent P+G","P & Add normalized E","P & Add normalized G","Remove adjacent P+E & \nAdd normalized E","Remove adjacent P+G & \nAdd normalized G"), "Spearman's correlation between RNA and ATAC", c(0.2,0.65), ncol[1:9])
Boxplot_List(list(BMMC_cor[1,], BMMC_cor[2,], BMMC_cor[3,], BMMC_cor[4,], BMMC_cor[5,], BMMC_cor[6,], BMMC_cor[7,], BMMC_cor[8,], BMMC_cor[9,]), "BMMC Same Donors",
           c("P","G","P & Remove adjacent P","P & Remove adjacent P+E","P & Remove adjacent P+G","P & Add normalized E","P & Add normalized G","Remove adjacent P+E & \nAdd normalized E","Remove adjacent P+G & \nAdd normalized G"), "Spearman's correlation between RNA and ATAC", c(0.2,0.65), ncol[1:9])
dev.off()

































