library(MAESTRO)
library(ggplot2)
library(aricode)
library(cluster)
library(scales)
library(lawstat)
library(RColorBrewer)
ncol <- brewer.pal(8,"Set1")

# cellline data analysis
{
  cellline_scABC <- readRDS('clustering_benchmark_methods_public_data/GSE65360_cellline_scABC_cluster.rds')
  cellline_scABC_umap <- readRDS('clustering_benchmark_methods_public_data/GSE65360_cellline_scABC_umap.rds')
  cellline_cisTopic <- readRDS('clustering_benchmark_methods_public_data/GSE65360_cellline_cisTopic_cluster.rds')
  cellline_cisTopic_umap <- readRDS('clustering_benchmark_methods_public_data/GSE65360_cellline_cisTopic_umap.rds')
  cellline_snapATAC <- readRDS('clustering_benchmark_methods_public_data/GSE65360_cellline_snapATAC_cluster.rds')
  cellline_snapATAC_umap <- readRDS('clustering_benchmark_methods_public_data/GSE65360_cellline_snapATAC_umap.rds')
  cellline_MAESTRO <- readRDS('clustering_benchmark_methods_public_data/GSE65360_cellline_MAESTRO_SeuratObj.rds')

  get_orign_index <- function(orig, umap) {x<-NULL;y<-NULL;
                                 for(celltype in unique(orig)) {x <- c(x,mean(umap[orig==celltype,1])); y <- c(y,mean(umap[orig==celltype,2]))}
                                 return(cbind(x,y))}
 
  cellline_scABC_orig <-  gsub(".well","",sapply(strsplit(names(cellline_scABC),"\\."),function(x) paste0(x[2],'.',x[3])))
  cellline_cisTopic_orig <- gsub(".well","",sapply(strsplit(names(cellline_cisTopic),"\\."),function(x) paste0(x[2],'.',x[3])))
  cellline_snapATAC_orig <- gsub(".WELL","",sapply(strsplit(names(cellline_snapATAC),"\\-"),function(x) paste0(x[2],'.',x[3])))
  cellline_MAESTRO_orig <- gsub(".well","",sapply(strsplit(rownames(cellline_MAESTRO@meta.data),"\\."),function(x) paste0(x[2],'.',x[3])))
  cellline_scABC_score <- NMI(cellline_scABC_orig, as.numeric(cellline_scABC))
  cellline_cisTopic_score <- NMI(cellline_cisTopic_orig, as.numeric(cellline_cisTopic))   
  cellline_snapATAC_score <- NMI(cellline_snapATAC_orig, as.numeric(cellline_snapATAC))
  cellline_MAESTRO_score <- NMI(cellline_MAESTRO_orig, cellline_MAESTRO@meta.data[,"seurat_clusters"])
  clusters <- rev(hue_pal()(16))
  scABC_col <- plyr::mapvalues(x = cellline_scABC_orig, from = sort(unique(cellline_scABC_orig)), to = clusters)
  cisTopic_col <- plyr::mapvalues(x = cellline_cisTopic_orig, from = sort(unique(cellline_cisTopic_orig)), to = clusters)
  snapATAC_col <- plyr::mapvalues(x = cellline_snapATAC_orig, from = sort(unique(cellline_snapATAC_orig)), to = clusters)
  MAESTRO_col <- plyr::mapvalues(x = cellline_MAESTRO_orig, from = sort(unique(cellline_MAESTRO_orig)), to = clusters)

  pdf("GSE65360_cellline_NMI_comparison.pdf",width=3.5,height=5)
  par(mar=c(8,5,3,3))
  s <-  c(cellline_scABC_score,cellline_cisTopic_score, cellline_snapATAC_score, cellline_MAESTRO_score)
  p <- barplot(s,col=ncol[2], ylab="NMI",las=2, names=c("scABC","cisTopic","snapATAC","LSI"),border=NA,main="celllines wt conditions",ylim=c(0,0.8)) 
  text(p, s+0.025, round(s,3),cex=0.65)
  dev.off() 

  pdf("GSE65360_cellline_umap_comparison.pdf",width=12,height=3.25)
  par(mfrow=c(1,4))
  plot(cellline_scABC_umap[,1], cellline_scABC_umap[,2], xlab="UMAP_1", ylab="UMAP_2", col=scABC_col,main="Cellline scABC",pch=20);legend("topleft",paste0("NMI = ",round(cellline_scABC_score,3)),box.lty=0)
  text(get_orign_index(cellline_scABC_orig,cellline_scABC_umap)[,1],get_orign_index(cellline_scABC_orig,cellline_scABC_umap)[,2],unique(cellline_scABC_orig))
  plot(cellline_cisTopic_umap@dr$cell$Umap[,1], cellline_cisTopic_umap@dr$cell$Umap[,2], xlab="UMAP_1", ylab="UMAP_2", col=cisTopic_col,main="Cellline cisTopic",pch=20);legend("topleft",paste0("NMI = ",round(cellline_cisTopic_score,3)),box.lty=0)
  text(get_orign_index(cellline_cisTopic_orig,cellline_cisTopic_umap@dr$cell$Umap)[,1],get_orign_index(cellline_cisTopic_orig,cellline_cisTopic_umap@dr$cell$Umap)[,2],unique(cellline_cisTopic_orig))
  plot(cellline_snapATAC_umap[,1], cellline_snapATAC_umap[,2], xlab="UMAP_1", ylab="UMAP_2", col=snapATAC_col,main="Cellline snapATAC",pch=20);legend("topright",paste0("NMI = ",round(cellline_snapATAC_score,3)),box.lty=0)
  text(get_orign_index(cellline_snapATAC_orig,cellline_snapATAC_umap)[,1],get_orign_index(cellline_snapATAC_orig,cellline_snapATAC_umap)[,2],unique(cellline_snapATAC_orig))
  plot(cellline_MAESTRO$umap@cell.embeddings[,1], cellline_MAESTRO$umap@cell.embeddings[,2], xlab="UMAP_1", ylab="UMAP_2", col=MAESTRO_col,main="Cellline LSI",pch=20);legend("topright",paste0("NMI = ",round(cellline_MAESTRO_score,3)),box.lty=0)
  text(get_orign_index(cellline_MAESTRO_orig,cellline_MAESTRO$umap@cell.embeddings)[,1],get_orign_index(cellline_MAESTRO_orig,cellline_MAESTRO$umap@cell.embeddings)[,2],unique(cellline_MAESTRO_orig))
  dev.off()
}

# HSC data analysis
{
  HSC_scABC <- readRDS('clustering_benchmark_methods_public_data/GSE96772_HSC_scABC_cluster.rds')
  HSC_scABC_umap <- readRDS('clustering_benchmark_methods_public_data/GSE96772_HSC_scABC_umap.rds')
  HSC_cisTopic <- readRDS('clustering_benchmark_methods_public_data/GSE96772_HSC_cisTopic_cluster.rds')
  HSC_cisTopic_umap <- readRDS('clustering_benchmark_methods_public_data/GSE96772_HSC_cisTopic_umap.rds')
  HSC_snapATAC <- readRDS('clustering_benchmark_methods_public_data/GSE96772_HSC_snapATAC_cluster.rds')
  HSC_snapATAC_umap <- readRDS('clustering_benchmark_methods_public_data/GSE96772_HSC_snapATAC_umap.rds')
  HSC_MAESTRO <- readRDS('clustering_benchmark_methods_public_data/GSE96772_HSC_MAESTRO_SeuratObj.rds')


  get_orign <- function(names) {samples <- sapply(strsplit(gsub("scatac\\.","",gsub("^[0-9]*\\.","",gsub("^[0-9]*\\.scATAC\\.","",gsub("singles\\.","",names)))),"\\."), function(x) x[2])
                                samples[samples=="141017"] <- "MEP"
                                return(samples)}
  get_orign_index <- function(orig, umap) {x<-NULL;y<-NULL;
                                 for(celltype in unique(orig)) {x <- c(x,mean(umap[orig==celltype,1])); y <- c(y,mean(umap[orig==celltype,2]))}
                                 return(cbind(x,y))}
  
  HSC_scABC_orig <- get_orign(names(HSC_scABC))
  HSC_cisTopic_orig <- get_orign(names(HSC_cisTopic))
  HSC_snapATAC_orig <- toupper(get_orign(tolower(gsub("-",".",names(HSC_snapATAC)))))
  HSC_MAESTRO_orig <- get_orign(rownames(HSC_MAESTRO@meta.data))
  HSC_scABC_score <- NMI(HSC_scABC_orig, as.numeric(HSC_scABC))
  HSC_cisTopic_score <- NMI(HSC_cisTopic_orig, as.numeric(HSC_cisTopic))   
  HSC_snapATAC_score <- NMI(HSC_snapATAC_orig, as.numeric(HSC_snapATAC))
  HSC_MAESTRO_score <- NMI(HSC_MAESTRO_orig, HSC_MAESTRO@meta.data[,"seurat_clusters"])
  clusters <- rev(hue_pal()(14))
  scABC_col <- plyr::mapvalues(x = HSC_scABC_orig, from = sort(unique(HSC_scABC_orig)), to = clusters)
  cisTopic_col <- plyr::mapvalues(x = HSC_cisTopic_orig, from = sort(unique(HSC_cisTopic_orig)), to = clusters)
  snapATAC_col <- plyr::mapvalues(x = HSC_snapATAC_orig, from = sort(unique(HSC_snapATAC_orig)), to = clusters)
  MAESTRO_col <- plyr::mapvalues(x = HSC_MAESTRO_orig, from = sort(unique(HSC_MAESTRO_orig)), to = clusters)

  pdf("GSE96772_HSC_NMI_comparison.pdf",width=3.5,height=5)
  par(mar=c(8,5,3,3))
  s <- c(HSC_scABC_score,HSC_cisTopic_score,HSC_snapATAC_score, HSC_MAESTRO_score)
  p <- barplot(s,col=ncol[2],ylab="NMI",las=2, names=c("scABC","cisTopic","snapATAC","LSI"),border=NA,main="HSC",ylim=c(0,0.7)) 
  text(p, s+0.025, round(s,3),cex=0.65)
  dev.off()

  pdf("GSE96772_HSC_umap_comparison.pdf",width=12,height=3.25)
  par(mfrow=c(1,4))
  plot(HSC_scABC_umap[,1], HSC_scABC_umap[,2], xlab="UMAP_1", ylab="UMAP_2", col=scABC_col,main="HSC scABC",pch=20);legend("topleft",paste0("NMI = ",round(HSC_scABC_score,3)),box.lty=0)
  text(get_orign_index(HSC_scABC_orig,HSC_scABC_umap)[,1],get_orign_index(HSC_scABC_orig,HSC_scABC_umap)[,2],unique(HSC_scABC_orig))
  plot(HSC_cisTopic_umap@dr$cell$Umap[,1], HSC_cisTopic_umap@dr$cell$Umap[,2], xlab="UMAP_1", ylab="UMAP_2", col=cisTopic_col,main="HSC cisTopic",pch=20);legend("top",paste0("NMI = ",round(HSC_cisTopic_score,3)),box.lty=0)
  text(get_orign_index(HSC_cisTopic_orig,HSC_cisTopic_umap@dr$cell$Umap)[,1],get_orign_index(HSC_cisTopic_orig,HSC_cisTopic_umap@dr$cell$Umap)[,2],unique(HSC_cisTopic_orig))
  plot(HSC_snapATAC_umap[,1], HSC_snapATAC_umap[,2], xlab="UMAP_1", ylab="UMAP_2", col=snapATAC_col,main="HSC snapATAC",pch=20);legend("topleft",paste0("NMI = ",round(HSC_snapATAC_score,3)),box.lty=0)
  text(get_orign_index(HSC_snapATAC_orig,HSC_snapATAC_umap)[,1],get_orign_index(HSC_snapATAC_orig,HSC_snapATAC_umap)[,2],unique(HSC_snapATAC_orig))
  plot(HSC_MAESTRO$umap@cell.embeddings[,1], HSC_MAESTRO$umap@cell.embeddings[,2], xlab="UMAP_1", ylab="UMAP_2", col=MAESTRO_col,main="HSC LSI",pch=20);legend("topright",paste0("NMI = ",round(HSC_MAESTRO_score,3)),box.lty=0)
  text(get_orign_index(HSC_MAESTRO_orig,HSC_MAESTRO$umap@cell.embeddings)[,1],get_orign_index(HSC_MAESTRO_orig,HSC_MAESTRO$umap@cell.embeddings)[,2],unique(HSC_MAESTRO_orig))
  dev.off()
}

# PBMC data analysis
{
  PBMC_genescore <- read.table('clustering_benchmark_methods_public_data/10X_PBMC_gene_score.txt')
  PBMC_markerscore <- na.omit(PBMC_genescore[as.character(read.table('pbmc.maker.genes.txt')[,1]),])
  PBMC_housekeepingscore <- na.omit(PBMC_genescore[as.character(read.table('pbmc.housekeeping.genes.txt')[,1]),])
  
  PBMC_scABC <- readRDS('clustering_benchmark_methods_public_data/10X_PBMC_scABC_cluster.rds')
  PBMC_scABC_umap <- readRDS('clustering_benchmark_methods_public_data/10X_PBMC_scABC_umap.rds')
  PBMC_cisTopic <- readRDS('clustering_benchmark_methods_public_data/10X_PBMC_cisTopic_cluster.rds')
  PBMC_cisTopic_umap <- readRDS('clustering_benchmark_methods_public_data/10X_PBMC_cisTopic_umap.rds')
  PBMC_snapATAC_raw <- readRDS('clustering_benchmark_methods_public_data/10X_PBMC_snapATAC_cluster.rds') 
  PBMC_snapATAC <- as.numeric(PBMC_snapATAC_raw); names(PBMC_snapATAC) <- gsub("-",".",names(PBMC_snapATAC_raw))
  PBMC_snapATAC_umap <- readRDS('clustering_benchmark_methods_public_data/10X_PBMC_snapATAC_umap.rds')
  PBMC_MAESTRO <- readRDS('clustering_benchmark_methods_public_data/10X_PBMC_MAESTRO_SeuratObj.rds')$ATAC
   
  calculate_RAGI_score <- function(cluster)
  {
    marker_avg <- NULL
    housekeeping_avg <- NULL
    cluster <- cluster[intersect(names(cluster),colnames(PBMC_genescore))]
    for(i in unique(cluster))
    {
      if(length(names(cluster[cluster==i]))>1){
         marker_avg <- cbind(marker_avg, apply(PBMC_markerscore[,names(cluster[cluster==i])],1,mean))
         housekeeping_avg <- cbind(housekeeping_avg, apply(PBMC_housekeepingscore[,names(cluster[cluster==i])],1,mean))}
      else{
         marker_avg <- cbind(marker_avg, PBMC_markerscore[,names(cluster[cluster==i])])
         housekeeping_avg <- cbind(housekeeping_avg, PBMC_housekeepingscore[,names(cluster[cluster==i])])}      
    }
    marker_gini <- mean(na.omit(apply(marker_avg,1,function(x) gini.index(x)$statistic)))
    housekeeping_gini <- mean(na.omit(apply(housekeeping_avg,1,function(x) gini.index(x)$statistic)))
    return(marker_gini-housekeeping_gini)
  }
  PBMC_scABC_score <- calculate_RAGI_score(PBMC_scABC)
  PBMC_cisTopic_score <- calculate_RAGI_score(PBMC_cisTopic)
  PBMC_snapATAC_score <- calculate_RAGI_score(PBMC_snapATAC)
  PBMC_MAESTRO_cluster <- as.numeric(PBMC_MAESTRO@meta.data[,"seurat_clusters"])-1; names(PBMC_MAESTRO_cluster) <- rownames(PBMC_MAESTRO@meta.data)
  PBMC_MAESTRO_score <- calculate_RAGI_score(PBMC_MAESTRO_cluster)
 
  pdf("PBMC_10K_RAGI_comparison.pdf",width=3.5,height=5)
  par(mar=c(8,5,3,3))
  s <- c(PBMC_scABC_score,PBMC_cisTopic_score,PBMC_snapATAC_score, PBMC_MAESTRO_score)
  p <- barplot(s,col=ncol[2],ylab="RAGI",las=2, names=c("scABC","cisTopic","snapATAC","LSI"),border=NA,main="PBMC",ylim=c(0,0.15)) 
  text(p, s+0.01, round(s,3),cex=0.65)
  dev.off()
    
  clusters <- rev(hue_pal()(17))
  scABC_col <- plyr::mapvalues(x = PBMC_scABC, from = sort(unique(PBMC_scABC)), to = clusters[1:10])
  cisTopic_col <- plyr::mapvalues(x = PBMC_cisTopic, from = sort(unique(PBMC_cisTopic)), to = clusters[1:9])
  snapATAC_col <- plyr::mapvalues(x = PBMC_snapATAC, from = sort(unique(PBMC_snapATAC)), to = clusters[1:15])
  MAESTRO_col <- plyr::mapvalues(x = PBMC_MAESTRO_cluster, from = sort(unique(PBMC_MAESTRO_cluster)), to = clusters[1:17])
  pdf("PBMC_10K_umap_comparison.pdf",width=12,height=3.25)
  par(mfrow=c(1,4))
  plot(PBMC_scABC_umap[1:1000,1], PBMC_scABC_umap[1:1000,2], xlab="UMAP_1", ylab="UMAP_2", col=scABC_col,main="PBMC scABC",pch=20);legend("topleft",paste0("RAGI = ",round(PBMC_scABC_score,3)),box.lty=0)
  plot(PBMC_cisTopic_umap@dr$cell$Umap[1:1000,1], PBMC_cisTopic_umap@dr$cell$Umap[1:1000,2], xlab="UMAP_1", ylab="UMAP_2", col=cisTopic_col,main="PBMC cisTopic",pch=20);legend("topleft",paste0("RAGI = ",round(PBMC_cisTopic_score,3)),box.lty=0)
  plot(PBMC_snapATAC_umap[1:1000,1], PBMC_snapATAC_umap[1:1000,2], xlab="UMAP_1", ylab="UMAP_2", col=snapATAC_col,main="PBMC snapATAC",pch=20);legend("top",paste0("RAGI = ",round(PBMC_snapATAC_score,3)),box.lty=0)
  plot(PBMC_MAESTRO$umap@cell.embeddings[1:1000,1], PBMC_MAESTRO$umap@cell.embeddings[1:1000,2], xlab="UMAP_1", ylab="UMAP_2", col=MAESTRO_col,main="PBMC LSI",pch=20);legend("topright",paste0("RAGI = ",round(PBMC_MAESTRO_score,3)),box.lty=0)
  dev.off()
}





















