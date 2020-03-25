library(MAESTRO)
library(aricode)
library(cluster)
library(scales)
library(lawstat)
library(RColorBrewer)
ncol <- brewer.pal(8,"Set1")

# load cellline data from different peak count set, v1 SC only, v2 ENCODE only, v3 ENCODE+SC
MASTARObj1$ATAC <- readRDS('cluster_bentchmark_peaks/GSE65360_cellline_merged_count_v1.rds')
MASTARObj2$ATAC <- readRDS('cluster_bentchmark_peaks/GSE65360_cellline_merged_count_v2.rds')
MASTARObj3$ATAC <- readRDS('cluster_bentchmark_peaks/GSE65360_cellline_merged_count_v3.rds')

# rename the sample name and NMI analysis
get_orign_index <- function(orig, umap) {x<-NULL;y<-NULL;
                               for(celltype in unique(orig)) {x <- c(x,mean(umap[orig==celltype,1])); y <- c(y,mean(umap[orig==celltype,2]))}
                               return(cbind(x,y))}
MASTARObj1$ATAC@meta.data[,"orig.ident"] <- sapply(strsplit(rownames(MASTARObj1$ATAC@meta.data),"\\."),function(x) paste0(x[2],'.',x[3]))
cellline_peak_v1_score <- NMI(MASTARObj1$ATAC@meta.data[,"orig.ident"], MASTARObj1$ATAC@meta.data[,"ATAC_snn_res.0.6"])
MASTARObj2$ATAC@meta.data[,"orig.ident"] <- sapply(strsplit(rownames(MASTARObj2$ATAC@meta.data),"\\."),function(x) paste0(x[2],'.',x[3]))
cellline_peak_v2_score <- NMI(MASTARObj2$ATAC@meta.data[,"orig.ident"], MASTARObj2$ATAC@meta.data[,"ATAC_snn_res.0.6"])
MASTARObj3$ATAC@meta.data[,"orig.ident"] <- sapply(strsplit(rownames(MASTARObj3$ATAC@meta.data),"\\."),function(x) paste0(x[2],'.',x[3]))
cellline_peak_v3_score <- NMI(MASTARObj3$ATAC@meta.data[,"orig.ident"], MASTARObj3$ATAC@meta.data[,"ATAC_snn_res.0.6"])

# output scomparison figures
pdf("GSE65360_cellline_NMI_comparison.pdf",width=3,height=5)
par(mar=c(12,5,3,3))
s <- c(cellline_peak_v1_score,cellline_peak_v2_score,cellline_peak_v3_score)
p <- barplot(s, col=ncol[2], ylab="NMI",las=2, names=c("SC aggregated","ENCODE rDHS","SC + ENCODE"),border=NA,main="celllines", ylim=c(0,0.8)) 
text(p, s+0.025, round(s,3),cex=0.65)
dev.off()
clusters <- rev(hue_pal()(16))
pdf("GSE65360_cellline_umap_comparison.pdf",width=9,height=3.25)
par(mfrow=c(1,3))
MAESTRO_col <- plyr::mapvalues(x = MASTARObj1$ATAC@meta.data[,"orig.ident"], from = sort(unique(MASTARObj1$ATAC@meta.data[,"orig.ident"])), to = clusters)
plot(MASTARObj1$ATAC$umap@cell.embeddings[,1], MASTARObj1$ATAC$umap@cell.embeddings[,2], xlab="UMAP_1", ylab="UMAP_2", col=MAESTRO_col,main="Cellline SC aggregated",pch=20);legend("topright",paste0("NMI = ",round(cellline_peak_v1_score,3)),box.lty=0)
text(get_orign_index(MASTARObj1$ATAC@meta.data[,"orig.ident"],MASTARObj1$ATAC$umap@cell.embeddings)[,1],get_orign_index(MASTARObj1$ATAC@meta.data[,"orig.ident"],MASTARObj1$ATAC$umap@cell.embeddings)[,2],unique(MASTARObj1$ATAC@meta.data[,"orig.ident"]))
MAESTRO_col <- plyr::mapvalues(x = MASTARObj2$ATAC@meta.data[,"orig.ident"], from = sort(unique(MASTARObj2$ATAC@meta.data[,"orig.ident"])), to = clusters)
plot(MASTARObj2$ATAC$umap@cell.embeddings[,1], MASTARObj2$ATAC$umap@cell.embeddings[,2], xlab="UMAP_1", ylab="UMAP_2", col=MAESTRO_col,main="Cellline ENCODE rDHS",pch=20);legend("topright",paste0("NMI = ",round(cellline_peak_v2_score,3)),box.lty=0)
text(get_orign_index(MASTARObj2$ATAC@meta.data[,"orig.ident"],MASTARObj2$ATAC$umap@cell.embeddings)[,1],get_orign_index(MASTARObj2$ATAC@meta.data[,"orig.ident"],MASTARObj2$ATAC$umap@cell.embeddings)[,2],unique(MASTARObj2$ATAC@meta.data[,"orig.ident"]))
MAESTRO_col <- plyr::mapvalues(x = MASTARObj3$ATAC@meta.data[,"orig.ident"], from = sort(unique(MASTARObj3$ATAC@meta.data[,"orig.ident"])), to = clusters)
plot(MASTARObj3$ATAC$umap@cell.embeddings[,1], MASTARObj3$ATAC$umap@cell.embeddings[,2], xlab="UMAP_1", ylab="UMAP_2", col=MAESTRO_col,main="Cellline SC + ENCODE",pch=20);legend("topright",paste0("NMI = ",round(cellline_peak_v3_score,3)),box.lty=0)
text(get_orign_index(MASTARObj3$ATAC@meta.data[,"orig.ident"],MASTARObj3$ATAC$umap@cell.embeddings)[,1],get_orign_index(MASTARObj3$ATAC@meta.data[,"orig.ident"],MASTARObj3$ATAC$umap@cell.embeddings)[,2],unique(MASTARObj2$ATAC@meta.data[,"orig.ident"]))
dev.off()

# load HSC data from different peak count set, v1 SC only, v2 ENCODE only, v3 ENCODE+SC
MASTARObj1$ATAC <- readRDS('cluster_bentchmark_peaks/GSE96772_HSC_merged_count_v1.rds')
MASTARObj2$ATAC <- readRDS('cluster_bentchmark_peaks/GSE96772_HSC_merged_count_v2.rds')
MASTARObj3$ATAC <- readRDS('cluster_bentchmark_peaks/GSE96772_HSC_merged_count_v3.rds')

# rename the sample name and NMI analysis
get_orign <- function(names) {samples <- sapply(strsplit(gsub("scatac\\.","",gsub("^[0-9]*\\.","",gsub("^[0-9]*\\.scATAC\\.","",gsub("singles\\.","",names)))),"\\."), function(x) x[2])
                              samples[samples=="141017"] <- "MEP"
                              return(samples)}
get_orign_index <- function(orig, umap) {x<-NULL;y<-NULL;
                              for(celltype in unique(orig)) {x <- c(x,mean(umap[orig==celltype,1])); y <- c(y,mean(umap[orig==celltype,2]))}
                              return(cbind(x,y))}
MASTARObj1$ATAC@meta.data[,"orig.ident"] <- get_orign(rownames(MASTARObj1$ATAC@meta.data))
HSC_peak_v1_score <- NMI(MASTARObj1$ATAC@meta.data[,"orig.ident"], MASTARObj1$ATAC@meta.data[,"ATAC_snn_res.0.6"])
MASTARObj2$ATAC@meta.data[,"orig.ident"] <- get_orign(rownames(MASTARObj2$ATAC@meta.data))
HSC_peak_v2_score <- NMI(MASTARObj2$ATAC@meta.data[,"orig.ident"], MASTARObj2$ATAC@meta.data[,"ATAC_snn_res.0.6"])
MASTARObj3$ATAC@meta.data[,"orig.ident"] <- get_orign(rownames(MASTARObj3$ATAC@meta.data))
HSC_peak_v3_score <- NMI(MASTARObj3$ATAC@meta.data[,"orig.ident"], MASTARObj3$ATAC@meta.data[,"ATAC_snn_res.0.6"])

# ccomparison figures
pdf("GSE96772_HSC_NMI_comparison.pdf",width=3,height=5)
par(mar=c(12,5,3,3))
s <- c(HSC_peak_v1_score, HSC_peak_v2_score, HSC_peak_v3_score)
p <- barplot(s, col=ncol[2], ylab="NMI",las=2, names=c("SC aggregated","ENCODE rDHS","SC + ENCODE"),border=NA,main="HSC", ylim=c(0,0.7)) 
text(p, s+0.025, round(s,3),cex=0.65)
dev.off()  
clusters <- rev(hue_pal()(14))
pdf("GSE96772_HSC_umap_comparison.pdf",width=9,height=3.25)
par(mfrow=c(1,3))
MAESTRO_col <- plyr::mapvalues(x = MASTARObj1$ATAC@meta.data[,"orig.ident"], from = sort(unique(MASTARObj1$ATAC@meta.data[,"orig.ident"])), to = clusters)
plot(MASTARObj1$ATAC$umap@cell.embeddings[,1], MASTARObj1$ATAC$umap@cell.embeddings[,2], xlab="UMAP_1", ylab="UMAP_2", col=MAESTRO_col,main="HSC SC aggregated",pch=20);legend("topright",paste0("NMI = ",round(HSC_peak_v1_score,3)),box.lty=0)
text(get_orign_index(MASTARObj1$ATAC@meta.data[,"orig.ident"],MASTARObj1$ATAC$umap@cell.embeddings)[,1],get_orign_index(MASTARObj1$ATAC@meta.data[,"orig.ident"],MASTARObj1$ATAC$umap@cell.embeddings)[,2],unique(MASTARObj1$ATAC@meta.data[,"orig.ident"]))
MAESTRO_col <- plyr::mapvalues(x = MASTARObj2$ATAC@meta.data[,"orig.ident"], from = sort(unique(MASTARObj2$ATAC@meta.data[,"orig.ident"])), to = clusters)
plot(MASTARObj2$ATAC$umap@cell.embeddings[,1], MASTARObj2$ATAC$umap@cell.embeddings[,2], xlab="UMAP_1", ylab="UMAP_2", col=MAESTRO_col,main="HSC ENCODE rDHS",pch=20);legend("topright",paste0("NMI = ",round(HSC_peak_v2_score,3)),box.lty=0)
text(get_orign_index(MASTARObj2$ATAC@meta.data[,"orig.ident"],MASTARObj2$ATAC$umap@cell.embeddings)[,1],get_orign_index(MASTARObj2$ATAC@meta.data[,"orig.ident"],MASTARObj2$ATAC$umap@cell.embeddings)[,2],unique(MASTARObj2$ATAC@meta.data[,"orig.ident"]))
MAESTRO_col <- plyr::mapvalues(x = MASTARObj3$ATAC@meta.data[,"orig.ident"], from = sort(unique(MASTARObj3$ATAC@meta.data[,"orig.ident"])), to = clusters)
plot(MASTARObj3$ATAC$umap@cell.embeddings[,1], MASTARObj3$ATAC$umap@cell.embeddings[,2], xlab="UMAP_1", ylab="UMAP_2", col=MAESTRO_col,main="HSC SC + ENCODE",pch=20);legend("topright",paste0("NMI = ",round(HSC_peak_v3_score,3)),box.lty=0)
text(get_orign_index(MASTARObj3$ATAC@meta.data[,"orig.ident"],MASTARObj3$ATAC$umap@cell.embeddings)[,1],get_orign_index(MASTARObj3$ATAC@meta.data[,"orig.ident"],MASTARObj3$ATAC$umap@cell.embeddings)[,2],unique(MASTARObj2$ATAC@meta.data[,"orig.ident"]))
dev.off()



