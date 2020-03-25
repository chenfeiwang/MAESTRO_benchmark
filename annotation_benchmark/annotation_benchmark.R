library(MAESTRO)
library(Seurat)
options(bitmapType='cairo')

source('annotation_cross_validation.R')
source('annotation_evaluate.R')
source('annotation_run_Garnett_CV.R')

PBMC_sorted <- read.csv('Intra-dataset/Zheng_sorted/Filtered_DownSampled_SortedPBMC_data.csv',row.names = 1)
PBMC_labels <- as.matrix(read.csv('Intra-dataset/Zheng_sorted/Labels.csv'))

# Step 1 annotation
# SCINA annotation
library(SCINA)
library(preprocessCore)
Signature_Genes <- preprocess.signatures('Intra-dataset/Zheng_sorted/SCINA_Gene_signature.csv')
Data = t(as.matrix(PBMC_sorted))
Data=log(Data+1)
Data[]=normalize.quantiles(Data)
results = SCINA(Data, Signature_Genes)
write.csv(results$cell_labels,'SCINA_Zheng_sorted_Pred_Labels.csv',row.names = FALSE)

# MAESTRO annotation
PBMC_sorted_res <- RNARunSeurat(inputMat = t(PBMC_sorted), 
                                project = "PBMC_sorted", 
                                min.c = 10,
                                min.g = 0,
                                mito = TRUE,
                                dims.use = 1:15,
                                variable.genes = 2000, 
                                organism = "GRCh38",
                                cluster.res = 0.6,
                                genes.test.use = "wilcox",
                                genes.cutoff = 1e-05)
PBMC_sorted_res$RNA <- RNAAnnotateCelltype(RNA = PBMC_sorted_res$RNA, 
                                           gene = PBMC_sorted_res$gene,
                                           signatures = human.immune.CIBERSORT, 
                                           min.score = 0.1)
write.csv(PBMC_sorted_res$RNA@meta.data[,"assign.ident"],'MAESTRO_Zheng_sorted_Pred_Labels.csv',row.names = FALSE)
saveRDS(PBMC_sorted_res,'PBMC_sorted_res.rds')

# GarnettCV annotation
Cross_Validation('Intra-dataset/Zheng_sorted/Labels.csv', 1, './')
run_Garnett_CV('Intra-dataset/Zheng_sorted/Filtered_DownSampled_SortedPBMC_data.csv', 
	           'Intra-dataset/Zheng_sorted/Labels.csv', 
	           'CV_folds.RData', 
	           'Intra-dataset/Zheng_sorted/Garnett_gene_names.csv',
	           'Intra-dataset/Zheng_sorted/Garnett_hsPBMC_markers.txt', 
	           './',
	           TRUE)

# Step 2 Clustering accuracy evaluation
# Replace the labels to calculate the accuracy
PBMC_sorted_res <- readRDS("PBMC_sorted_res.rds")
PBMC_sorted_res$RNA@meta.data[,"assign.ident"] <- correct_labels(as.character(PBMC_sorted_res$RNA@meta.data[,"assign.ident"]))
p <- DimPlot(object = PBMC_sorted_res$RNA, pt.size = 0.5, group.by = "assign.ident", label = TRUE) + ggtitle("MAESTRO Zheng sorted") + theme(plot.title = element_text(size = 12, face = "bold"))
ggsave(file.path(paste0(PBMC_sorted_res$RNA@project.name, "_MAESTRO_annotated.pdf")), p, width=5.5, height=4)
write.csv(PBMC_sorted_res$RNA@meta.data[,"assign.ident"],'Zheng_sorted_Pred_Labels_MAESTRO.csv',row.names = FALSE)

Labels <- as.matrix(read.csv('SCINA_Zheng_sorted_Pred_Labels.csv'))
PBMC_sorted_res$RNA@meta.data[,"assign.ident"] <- correct_labels(Labels[,1])
p <- DimPlot(object = PBMC_sorted_res$RNA, pt.size = 0.5, group.by = "assign.ident", label = TRUE, repel = TRUE) + ggtitle("SCINA Zheng sorted") + theme(plot.title = element_text(size = 12, face = "bold"))
ggsave(file.path(paste0(PBMC_sorted_res$RNA@project.name, "_SCINA_annotated.pdf")), p, width=5.5, height=4)
write.csv(PBMC_sorted_res$RNA@meta.data[,"assign.ident"],'Zheng_sorted_Pred_Labels_SCINA.csv',row.names = FALSE)

Labels <- as.matrix(read.csv('Garnett_CV_Zheng_sorted_Pred_Labels_LM22.csv'))
PBMC_sorted_res$RNA@meta.data[,"assign.ident"] <- correct_labels(Labels[,1])
p <- DimPlot(object = PBMC_sorted_res$RNA, pt.size = 0.5, group.by = "assign.ident", label = TRUE, repel = TRUE) + ggtitle("GarnettCV LM22 Zheng sorted") + theme(plot.title = element_text(size = 12, face = "bold"))
ggsave(file.path(paste0(PBMC_sorted_res$RNA@project.name, "_Garnett_LM22_annotated.pdf")), p, width=5.5, height=4)
write.csv(PBMC_sorted_res$RNA@meta.data[,"assign.ident"],'Zheng_sorted_Pred_Labels_Garnett_LM22.csv',row.names = FALSE)

Labels <- as.matrix(read.csv('Garnett_CV_Zheng_sorted_Pred_Labels_simple.csv'))
Labels[Labels == 'B cells'] <- 'B'
Labels[Labels == 'CD34+'] <- 'Granulocyte'
Labels[Labels == 'CD4 T cells'] <- 'CD4T'
Labels[Labels == 'CD8 T cells'] <- 'CD8T'
Labels[Labels == 'Dendritic cells'] <- 'Granulocyte'
Labels[Labels == 'Monocytes'] <- 'Monocyte'
Labels[Labels == 'NK cells'] <- 'NK'
PBMC_sorted_res$RNA@meta.data[,"assign.ident"] <- Labels[,1]
p <- DimPlot(object = PBMC_sorted_res$RNA, pt.size = 0.5, group.by = "assign.ident", label = TRUE, repel = TRUE) + ggtitle("GarnettCV simple Zheng sorted") + theme(plot.title = element_text(size = 12, face = "bold"))
ggsave(file.path(paste0(PBMC_sorted_res$RNA@project.name, "_Garnett_simple_annotated.pdf")), p, width=5.5, height=4)
write.csv(PBMC_sorted_res$RNA@meta.data[,"assign.ident"],'Zheng_sorted_Pred_Labels_Garnett_simple.csv',row.names = FALSE)

Labels <- as.matrix(read.csv('Intra-dataset/Zheng_sorted/Labels.csv'))
Labels[Labels == 'CD14+ Monocyte'] <- 'Monocyte'
Labels[Labels == 'CD19+ B'] <- 'B'
Labels[Labels == 'CD34+'] <- 'Granulocyte'
Labels[Labels == 'CD4+ T Helper2'] <- 'CD4T'
Labels[Labels == 'CD4+/CD25 T Reg'] <- 'Treg'
Labels[Labels == 'CD4+/CD45RA+/CD25- Naive T'] <- 'CD4T'
Labels[Labels == 'CD4+/CD45RO+ Memory'] <- 'CD4T'
Labels[Labels == 'CD56+ NK'] <- 'NK'
Labels[Labels == 'CD8+ Cytotoxic T'] <- 'CD8T'
Labels[Labels == 'CD8+/CD45RA+ Naive Cytotoxic'] <- 'CD8T'
PBMC_sorted_res$RNA@meta.data[,"assign.ident"] <- Labels[,1]
p <- DimPlot(object = PBMC_sorted_res$RNA, pt.size = 0.5, group.by = "assign.ident", label = TRUE, repel = TRUE) + ggtitle("True Zheng sorted") + theme(plot.title = element_text(size = 12, face = "bold"))
ggsave(file.path(paste0(PBMC_sorted_res$RNA@project.name, "_True_annotated.pdf")), p, width=5.5, height=4)
write.csv(PBMC_sorted_res$RNA@meta.data[,"assign.ident"],'Zheng_sorted_Pred_Labels_TRUE.csv',row.names = FALSE)

correct_labels <- function(inlabel)
{
	inlabel[inlabel=='ActDCs'] <- 'Granulocytes'
	inlabel[inlabel=='ActMast'] <- 'Granulocytes'
	inlabel[inlabel=='ActMemCD4Tcells'] <- 'CD4T'
	inlabel[inlabel=='ActNK'] <- 'NK'
	inlabel[inlabel=='CD8Tcells'] <- 'CD8T'
	inlabel[inlabel=='Eosinophils'] <- 'Granulocytes'
	inlabel[inlabel=='MacrophagesM0'] <- 'Monocyte'
	inlabel[inlabel=='MacrophagesM1'] <- 'Monocyte'
	inlabel[inlabel=='MacrophagesM2'] <- 'Monocyte'
	inlabel[inlabel=='MemoryBcells'] <- 'B'
	inlabel[inlabel=='Monocytes'] <- 'Monocyte'
	inlabel[inlabel=='NaiveBcells'] <- 'B'
	inlabel[inlabel=='NaiveCD4Tcells'] <- 'CD4T'
	inlabel[inlabel=='Neutrophils'] <- 'Granulocytes'
	inlabel[inlabel=='PlasmaCells'] <- 'B'
	inlabel[inlabel=='RestDCs'] <- 'Granulocytes'
	inlabel[inlabel=='RestMast'] <- 'Granulocytes'
	inlabel[inlabel=='RestMemCD4Tcells'] <- 'CD4T'
	inlabel[inlabel=='RestNK'] <- 'NK'
	inlabel[inlabel=='Tfh'] <- 'CD4T'
	inlabel[inlabel=='Tgd'] <- 'CD4T'
	inlabel[inlabel=='Treg'] <- 'Treg'
	return(inlabel)
}

MAESTRO <- evaluate("Zheng_sorted_Pred_Labels_TRUE.csv","Zheng_sorted_Pred_Labels_MAESTRO.csv")
SCINA <- evaluate("Zheng_sorted_Pred_Labels_TRUE.csv","Zheng_sorted_Pred_Labels_SCINA.csv")
GarnettLM22 <- evaluate("Zheng_sorted_Pred_Labels_TRUE.csv","Zheng_sorted_Pred_Labels_Garnett_LM22.csv")
GarnettSimple <- evaluate("Zheng_sorted_Pred_Labels_TRUE.csv","Zheng_sorted_Pred_Labels_Garnett_simple.csv")

library(RColorBrewer)
ncol <- brewer.pal(8,"Set1")
pdf("annotation_bentchmark.pdf",width=3.5,height=5)
par(mar=c(8,5,3,3))
s <-  c(MAESTRO$MedF1, SCINA$MedF1, GarnettLM22$MedF1, GarnettSimple$MedF1)
p <- barplot(s,col=ncol[2], ylab="Median F1 score",las=2, names=c("MAESTRO LM22","SCINA LM22","Garnett LM22","Garnett Simple"),border=NA, main="Zheng sorted",ylim=c(0,0.8)) 
text(p, s+0.025, round(s,3),cex=0.65)
dev.off() 















































