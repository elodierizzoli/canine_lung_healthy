#############################################################################################################
#### A single cell RNA sequencing atlas of the healthy canine lung: a foundation for comparative studies ####
####                           DIFFERENTIAL GENE EXPRESSION PSEUDOBULK ANALYSIS                          ####
#############################################################################################################

# # Load library
library(Seurat)
library(SeuratObject)
library("SingleCellExperiment")
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(cowplot)
library(Glimma)
library(RColorBrewer)
library(pheatmap)
library(gplots)
library(dendextend)

#### DEFINING PATH
# set working directory

# Bring in Seurat object
Healthy_integrated <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_integrated.rds")
OutDir <- paste0(wd, "/Outs/DESeq2")
ifelse(!dir.exists(OutDir), dir.create(OutDir), "OutDir exists")

### Using AggregateExpression()
CountMatrix <- AggregateExpression(Healthy_integrated,
                                   assays ="RNA", slot='count',
                                   group.by = c("celltype","orig.ident"), return.seurat = FALSE)
write.csv(CountMatrix, paste0(OutDir, "/AggregateCounts.csv"))

# Import Counts
Counts <- read.csv(paste0(wd, "/Outs/DESeq2/AggregateCounts.csv"),
                   row.names = 1)
colnames(Counts)

# Subset counts
## FIBROBLASTS CLUSTERS ####
# Alveolar fibroblasts
OutDir <- paste0(wd, '/Outs/DESeq2/ALVFvsFib')
ifelse(!dir.exists(OutDir), dir.create(OutDir), "OutDir exists")

Dataset1 <- subset(Counts, select = grep("Alveolar.fibroblasts", colnames(Counts)) ) #CL1
colnames(Dataset1)

Dataset3 <- Counts[,c(grep(".fibroblasts", colnames(Counts))) ]
colnames(Dataset3)
Dataset4 <- subset(Dataset3, select = !grepl("myofibroblasts", colnames(Dataset3)) )
colnames(Dataset4)
Dataset5 <- subset(Dataset4, select = !grepl("RNA.Alveolar", colnames(Dataset4)) )
colnames(Dataset5)

Dataset2 <-Dataset5[,1:4]
for(i in 1:4){
  Dataset2[,i] <- rowSums(Dataset5[,grep(paste0("LUNG",i), colnames(Dataset5))])
}

# CCL19 fibroblasts
OutDir <- paste0(wd, '/Outs/DESeq2/CCL19vsFib')
ifelse(!dir.exists(OutDir), dir.create(OutDir), "OutDir exists")

Dataset1 <- subset(Counts, select = grep("CCL19", colnames(Counts)) ) #CL1
colnames(Dataset1)

Dataset3 <- Counts[,c(grep("fibroblasts", colnames(Counts))) ]
colnames(Dataset3)
Dataset4 <- subset(Dataset3, select = !grepl("myofibroblasts", colnames(Dataset3)) )
colnames(Dataset4)
Dataset5 <- subset(Dataset4, select = !grepl("CCL19", colnames(Dataset4)) )
colnames(Dataset5)

Dataset2 <-Dataset5[,1:4]
for(i in 1:4){
  Dataset2[,i] <- rowSums(Dataset5[,grep(paste0("LUNG",i), colnames(Dataset5))])
}

# CCN3 fibroblasts
OutDir <- paste0(wd, '/Outs/DESeq2/CCN3vsFib')
ifelse(!dir.exists(OutDir), dir.create(OutDir), "OutDir exists")

Dataset1 <- subset(Counts, select = grep("CCN3", colnames(Counts)) ) #CL1
colnames(Dataset1)

Dataset3 <- Counts[,c(grep("fibroblasts", colnames(Counts))) ]
colnames(Dataset3)
Dataset4 <- subset(Dataset3, select = !grepl("myofibroblasts", colnames(Dataset3)) )
colnames(Dataset4)
Dataset5 <- subset(Dataset4, select = !grepl("CCN3", colnames(Dataset4)) )
colnames(Dataset5)

Dataset2 <-Dataset5[,1:4]
for(i in 1:4){
  Dataset2[,i] <- rowSums(Dataset5[,grep(paste0("LUNG",i), colnames(Dataset5))])
}

# STMN2 fibroblasts
OutDir <- paste0(wd, '/Outs/DESeq2/STMN2vsFib')
ifelse(!dir.exists(OutDir), dir.create(OutDir), "OutDir exists")

Dataset1 <- subset(Counts, select = grep("STMN2", colnames(Counts)) ) #CL1
colnames(Dataset1)

Dataset3 <- Counts[,c(grep("fibroblasts", colnames(Counts))) ]
colnames(Dataset3)
Dataset4 <- subset(Dataset3, select = !grepl("myofibroblasts", colnames(Dataset3)) )
colnames(Dataset4)
Dataset5 <- subset(Dataset4, select = !grepl("STMN2", colnames(Dataset4)) )
colnames(Dataset5)

Dataset2 <-Dataset5[,1:4]
for(i in 1:4){
  Dataset2[,i] <- rowSums(Dataset5[,grep(paste0("LUNG",i), colnames(Dataset5))])
}

# STMN2 fibroblasts
OutDir <- paste0(wd, '/Outs/DESeq2/STMN2vsAlvFib')
ifelse(!dir.exists(OutDir), dir.create(OutDir), "OutDir exists")

Dataset1 <- subset(Counts, select = grep("STMN2", colnames(Counts)) ) #CL1
colnames(Dataset1)

Dataset2 <- Counts[,c(grep("RNA.Alveolar.fibroblasts", colnames(Counts))) ]
colnames(Dataset2)

# COL23A1 fibroblasts
OutDir <- paste0(wd, '/Outs/DESeq2/COL23A1vsFib')
ifelse(!dir.exists(OutDir), dir.create(OutDir), "OutDir exists")

Dataset1 <- subset(Counts, select = grep("COL23A1", colnames(Counts)) ) #CL1
colnames(Dataset1)

Dataset3 <- Counts[,c(grep("fibroblasts", colnames(Counts))) ]
colnames(Dataset3)
Dataset4 <- subset(Dataset3, select = !grepl("myofibroblasts", colnames(Dataset3)) )
colnames(Dataset4)
Dataset5 <- subset(Dataset4, select = !grepl("COL23A1", colnames(Dataset4)) )
colnames(Dataset5)

Dataset2 <-Dataset5[,1:4]
for(i in 1:4){
  Dataset2[,i] <- rowSums(Dataset5[,grep(paste0("LUNG",i), colnames(Dataset5))])
}

# CCBE1 fibroblasts
OutDir <- paste0(wd, '/Outs/DESeq2/CCBE1vsFib')
ifelse(!dir.exists(OutDir), dir.create(OutDir), "OutDir exists")

Dataset1 <- subset(Counts, select = grep("CCBE1", colnames(Counts)) ) #CL1
colnames(Dataset1)

Dataset3 <- Counts[,c(grep("fibroblasts", colnames(Counts))) ]
colnames(Dataset3)
Dataset4 <- subset(Dataset3, select = !grepl("myofibroblasts", colnames(Dataset3)) )
colnames(Dataset4)
Dataset5 <- subset(Dataset4, select = !grepl("CCBE1", colnames(Dataset4)) )
colnames(Dataset5)

Dataset2 <-Dataset5[,1:4]
for(i in 1:4){
  Dataset2[,i] <- rowSums(Dataset5[,grep(paste0("LUNG",i), colnames(Dataset5))])
}


## MACROPHAGE/MONOCYTE CLUSTERS ####
# CCL13+ macrophages
OutDir <- paste0(wd, '/Outs/DESeq2/CCL13MPvsAM')
ifelse(!dir.exists(OutDir), dir.create(OutDir), "OutDir exists")

Dataset1 <- subset(Counts, select = grep("CCL13", colnames(Counts)) ) #CL1
colnames(Dataset1)

Dataset3 <- Counts[,c(grep("AM", colnames(Counts))) ]
colnames(Dataset3)
Dataset2 <- subset(Dataset3, select = !grepl("C1Q", colnames(Dataset3)) )
colnames(Dataset2)

# SLC27A6 macrophages
OutDir <- paste0(wd, '/Outs/DESeq2/SLC27A6vsCCL13')
ifelse(!dir.exists(OutDir), dir.create(OutDir), "OutDir exists")

Dataset1 <- subset(Counts, select = grep("SLC27A6", colnames(Counts)) ) #CL1
colnames(Dataset1)

Dataset2 <- subset(Counts, select = grep("CCL13", colnames(Counts)) ) #CL2
colnames(Dataset2)

# C1Q AM vs AM
OutDir <- paste0(wd, '/Outs/DESeq2/C1QvsAM')
ifelse(!dir.exists(OutDir), dir.create(OutDir), "OutDir exists")

Dataset1 <- subset(Counts, select = grep("C1Q", colnames(Counts)) ) #CL1
colnames(Dataset1)

Dataset3 <- Counts[,c(grep("AM", colnames(Counts))) ]
colnames(Dataset3)
Dataset2 <- subset(Dataset3, select = !grepl("C1Q", colnames(Dataset3)) )
colnames(Dataset2)

# FN1 vs monocytes
OutDir <- paste0(wd, '/Outs/DESeq2/FN1vsMono')
ifelse(!dir.exists(OutDir), dir.create(OutDir), "OutDir exists")

Dataset1 <- subset(Counts, select = grep("FN1", colnames(Counts)) ) #CL1
colnames(Dataset1)

Dataset2 <- Counts[,c(grep("Monocytes", colnames(Counts))) ]
colnames(Dataset2)

# CD1C monocytes vs Monocytes
OutDir <- paste0(wd, '/Outs/DESeq2/CD1CvsMono')
ifelse(!dir.exists(OutDir), dir.create(OutDir), "OutDir exists")

Dataset1 <- subset(Counts, select = grep("CD1C", colnames(Counts)) ) #CL1
colnames(Dataset1)

Dataset2 <- Counts[,c(grep("Monocytes", colnames(Counts))) ]
colnames(Dataset2)

# PREPARE DATASETS ####

colnames(Dataset1) <- c("CL1_LUNG1","CL1_LUNG2","CL1_LUNG3","CL1_LUNG4")
colnames(Dataset2) <- c("CL2_LUNG1","CL2_LUNG2","CL2_LUNG3","CL2_LUNG4")

Dataset <- cbind(Dataset1, Dataset2)

name <- sapply(strsplit(colnames(Dataset),"_"), "[[",2)
celltype <- sapply(strsplit(colnames(Dataset),"_"), "[[",1)
colData <- data.frame(name=name,
                      celltype=celltype)
rownames(colData) <- colnames(Dataset)

#################################################
## DESeq2 Analysis 
#################################################;
Dataset <- round(Dataset); # dataSet for DESeq need integer count
dds   <- DESeqDataSetFromMatrix(countData = Dataset,
                                colData   = colData,
                                design    = ~  celltype)

#################################################
## SOME DATAMANAGMENT 
#################################################;
## Define REF; 
dds$celltype       <- relevel(dds$celltype, ref = "CL2"); dds$celltype

#################################################
## DESeq Analysis
#################################################;
dds <- DESeq(dds)

vsd <- varianceStabilizingTransformation(dds, blind = TRUE) 
nsel <- 500

## Selection of nsel more variable genes 
expressmatrix <- t(assay(vsd))

varvsd <- apply(expressmatrix, MARGIN = 2, var)

topnsel <- order(varvsd, decreasing = TRUE)[1:nsel]

#################################################
## CLUSTERING (HEATMAP AND TREE) Exploratory stats
#################################################;
# 1. Distances with topnsel genes
expr_dist <- dist(expressmatrix[,topnsel], method = "euclidean")   # Distance between objects

# 2. Algorithme
expr_hclust <- hclust(expr_dist, method="ward.D2")
dendro <- as.dendrogram(expr_hclust)

pdf(paste(OutDir, "dendro.pdf", sep="/"), heigh=10, width=15)
##celltype
class  <- as.factor(colData$celltype)
color  <- as.numeric(class)+1; # color
color  <- color[order.dendrogram(dendro)]
dendextend::labels_colors(dendro) <- color
dendro   <- color_branches(dendro, h=2, col = color)

plot(dendro)
rect.hclust(expr_hclust, k=2)

##name
class  <- as.factor(colData$name)
color  <- as.numeric(class)+1; # color
color  <- color[order.dendrogram(dendro)]
dendextend::labels_colors(dendro) <- color
dendro   <- color_branches(dendro, h=4, col = color)

plot(dendro)
rect.hclust(expr_hclust, k=4)
dev.off()

class_order <- class[order.dendrogram(dendro)]
id     <- rownames(colData)
id_order    <- id[order.dendrogram(dendro)]

dendro_outliers <- as.data.frame(expr_hclust$order, expr_hclust$labels)
dendro_outliers$color <- color
dendro_outliers$class <- class_order
dendro_outliers$id <- id_order

write_tsv(x = dendro_outliers, file = paste(OutDir, "dendro_outliers.txt", sep="/"), col_names=TRUE)

# 2. Heatmap managment
expr_dist_Matrix <- as.matrix(expr_dist ) 
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
cond <- data.frame(Condition=colData$celltype)
rownames(cond)=rownames(colData)

pdf(paste(OutDir, "heatmap.pdf", sep="/"), width = 10, height = 10)
ph <- pheatmap(expr_dist_Matrix,
         clustering_distance_rows=expr_dist,
         clustering_distance_cols=expr_dist,
         col=colors,
         annotation_col=cond, treeheight_col = F)
dev.off()
# pdf("output/heatmap2.pdf", width = 15, height = 15)
# ht <- heatmap.2(expr_dist_Matrix, labRow = rownames(expr_dist_Matrix), labCol = colnames(expr_dist_Matrix),
#                 col = colorRampPalette(rev(brewer.pal(9, "Blues")) )(255), density.info = "none", trace="none",
#                 colRow = as.numeric(colData$condition)+1)
# dev.off()

ggsave(
  paste(OutDir, "pheatmap.pdf", sep="/"),
  plot = ph,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 40,
  height = 30,
  units = c("cm"),
  dpi = 600,
  limitsize = TRUE,
)

expressmatrix_sample <- expressmatrix[,topnsel[1:100]]

genenames <- colnames(expressmatrix_sample)
sampleID <- rownames(expressmatrix_sample)

pdf(paste(OutDir, "heatmapgene.pdf", sep="/"), width = 15, height = 12.5)
heatmap.2(expressmatrix_sample, labRow = sampleID, labCol = genenames,
          col = rev(heat.colors(75)), density.info = "none", trace="none",
          colRow = as.numeric(class)+1)
dev.off()

#################################################
## PCA
#################################################;
## GGPLOT2
# 1. WITH DESeq data
pdf(paste(OutDir, "PlotPCA.pdf", sep="/"), width = 10, height = 10)
plotPCA(vsd, intgroup=c("celltype"))
plotPCA(vsd, intgroup=c("name"))
dev.off()

# 2. With PRCOMP on TOP500 vst
dds_PCA <- prcomp(expressmatrix[,topnsel], retx = TRUE, center = TRUE, scale. = FALSE)

eigval <- dds_PCA$sdev^2
percent_var <- eigval*100/sum(eigval)
percent_var_cum <- cumsum(percent_var)

eig_res <- cbind(eigval, percent_var, percent_var_cum)
rownames(eig_res) <- paste0("PC", 1:dim(eig_res)[1])

dds_PCA_df <- data.frame(dds_PCA$x)

pc1lim <- ceiling(max(abs(min(dds_PCA_df[,1])), max(dds_PCA_df[,1]))/10)*10
pc2lim <- ceiling(max(abs(min(dds_PCA_df[,2])), max(dds_PCA_df[,2]))/10)*10
pc3lim <- ceiling(max(abs(min(dds_PCA_df[,3])), max(dds_PCA_df[,3]))/10)*10
pc4lim <- ceiling(max(abs(min(dds_PCA_df[,4])), max(dds_PCA_df[,4]))/10)*10

# colData$group <- factor(x = colData$group, levels = c("CTL4", "ODN4", "CLT24", "TEST4h", "CTRL24h", "TCR24h", "LPS24h", "R84824h"))

pcp_manual12 <- ggplot()+
  geom_point(data=dds_PCA_df, aes(x=PC1, y=PC2, col=colData$celltype, shape=colData$name), size=3) +
  scale_colour_manual('Cond', values = c("#180091", "#008391", "#8a0091",
                                         "#910000", "#915700", "#00911f",
                                         "#4dc2d9","#93e9be")) +
  scale_shape_manual('Cond', values = c(19,18,17,16,15,20,21,22)) +
  xlim(c(-pc1lim,pc1lim)) +
  ylim(c(-pc2lim,pc2lim)) +
  theme_bw() +                # Remove ugly grey background
  coord_quickmap() +          # Prevents stretching when resizing
  xlab(paste(colnames(dds_PCA_df)[1], ": ", round(eig_res[1,2],2), "% variance", sep="")) +
  ylab(paste(colnames(dds_PCA_df)[2], ": ", round(eig_res[2,2],2), "% variance", sep="")) +
  ggtitle("PRcomp PCA TOP500 var")

pcp_manual12_lab <- ggplot()+
  geom_point(data=dds_PCA_df, aes(x=PC1, y=PC2, col=colData$name, shape=colData$celltype), size=3) +
  scale_colour_manual('Cond', values = c("#180091", "#008391", "#8a0091",
                                         "#910000", "#915700", "#00911f",
                                         "#4dc2d9","#93e9be")) +
  scale_shape_manual('Cond', values = c(19,18,17,16,15,20,21,22)) +
  xlim(c(-pc1lim,pc1lim)) +
  ylim(c(-pc2lim,pc2lim)) +
  theme_bw() +             # Remove ugly grey background
  coord_quickmap() +          # Prevents stretching when resizing
  xlab(paste(colnames(dds_PCA_df)[1], ": ", round(eig_res[1,2],2), "% variance", sep="")) +
  ylab(paste(colnames(dds_PCA_df)[2], ": ", round(eig_res[2,2],2), "% variance", sep="")) +
  geom_label_repel(data = dds_PCA_df, aes(x=PC1, y=PC2, label=rownames(dds_PCA_df), col=colData$name),
                   size=3,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey75',
                   nudge_x=2,
                   nudge_y=2)+ 
  ggtitle("PRcomp PCA TOP500 var")

pcp_manual34 <- ggplot()+
  geom_point(data=dds_PCA_df, aes(x=PC3, y=PC4, col=colData$celltype, shape=colData$name), size=3) +
  scale_colour_manual('Cond', values = c("#180091", "#008391", "#8a0091",
                                         "#910000", "#915700", "#00911f",
                                         "#4dc2d9","#93e9be")) +
  scale_shape_manual('Cond', values = c(19,18,17,16,15,20,21,22)) +
  xlim(c(-pc3lim,pc3lim)) +
  ylim(c(-pc4lim,pc4lim)) +
  theme_bw() +             # Remove ugly grey background
  coord_quickmap() +          # Prevents stretching when resizing
  xlab(paste(colnames(dds_PCA_df)[3], ": ", round(eig_res[3,2],2), "% variance", sep="")) +
  ylab(paste(colnames(dds_PCA_df)[4], ": ", round(eig_res[4,2],2), "% variance", sep="")) +
  ggtitle("PRcomp PCA TOP500 var")

pcp_manual34_lab <-  ggplot()+
  geom_point(data=dds_PCA_df, aes(x=PC3, y=PC4, col=colData$name, shape=colData$celltype), size=3) +
  scale_colour_manual('Cond', values = c("#180091", "#008391", "#8a0091",
                                         "#910000", "#915700", "#00911f",
                                         "#4dc2d9","#93e9be")) +
  scale_shape_manual('Cond', values = c(19,18,17,16,15,20,21,22)) +
  xlim(c(-pc3lim,pc3lim)) +
  ylim(c(-pc4lim,pc4lim)) +
  theme_bw() +             # Remove ugly grey background
  coord_quickmap() +          # Prevents stretching when resizing
  xlab(paste(colnames(dds_PCA_df)[3], ": ", round(eig_res[3,2],2), "% variance", sep="")) +
  ylab(paste(colnames(dds_PCA_df)[4], ": ", round(eig_res[4,2],2), "% variance", sep="")) +
  geom_label_repel(data = dds_PCA_df, aes(x=PC3, y=PC4, label=rownames(dds_PCA_df), col=colData$name),
                   size=3,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey75',
                   nudge_x=2,
                   nudge_y=2)+ 
  ggtitle("PRcomp PCA TOP500 var")

p12 <- plot_grid(pcp_manual12,pcp_manual12_lab, labels=c("A","B"), ncol = 1, nrow = 2)
p34 <- plot_grid(pcp_manual34,pcp_manual34_lab, labels=c("A","B"), ncol = 1, nrow = 2)

ggsave(
  paste(OutDir, "PCA_12.pdf", sep="/"),
  plot = p12,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 40,
  height = 30,
  units = c("cm"),
  dpi = 600,
  limitsize = TRUE,
)
ggsave(
  paste(OutDir, "PCA_34.pdf", sep="/"),
  plot = p34,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 40,
  height = 30,
  units = c("cm"),
  dpi = 600,
  limitsize = TRUE,
)
save_plot(paste(OutDir, "PCA_12.png", sep="/"), p12, dpi = 600, base_width = 12.43, base_height = 8.42)
save_plot(paste(OutDir, "PCA_34.png", sep="/"), p34, dpi = 600, base_width = 12.43, base_height = 8.42)
dev.off()


################
### PAIRWISE ###
################
### Cond1vsCond2 
################
resultsNames(dds) # Gives comparisons

res<-results(dds, name=c("celltype_CL1_vs_CL2"), cooksCutoff = T, alpha = 0.05) 
summary(res)

resShrink<-lfcShrink(dds, coef=c("celltype_CL1_vs_CL2"), res=res)
plotMA(resShrink, ylim=c(-5,5), alpha = 0.05)

#################################################
## Creata AND Export Tables
#################################################;
counts_table <- as.data.frame(counts(dds, normalized=TRUE))
counts_table.raw <- as.data.frame(counts(dds, normalized=FALSE))

### EXTRACT means
CL1  <- which(colData$celltype == "CL1") 
CL2  <- which(colData$celltype == "CL2") 

avgall    <- apply(counts_table[,c(CL1, CL2)], 1, mean)
avgCL1     <- apply(counts_table[,c(CL1)], 1, mean)
avgCL2     <- apply(counts_table[,c(CL2)], 1, mean)

avgdf <- data.frame(avgall, avgCL1, avgCL2)

avgdf$ID <- rownames(avgdf)
counts_table$ID <- rownames(counts_table)
counts_table.raw$ID <- rownames(counts_table.raw)

## 2.1 Cond1 vs Cond2
resShk_df <- as.data.frame(resShrink)
resShk_df$ID <- rownames(resShk_df)
head(resShk_df)

resShk_df <- merge(x = resShk_df, y = avgdf[,],
                   by.x = "ID", by.y = "ID")

colnames(resShk_df)[7:9] <- c("CL1.CL2_BaseMean", "CL1_BaseMean", "CL2_BaseMean")

table_CL1_CL2  <- resShk_df[order(resShk_df$padj),]; head(table_CL1_CL2)

## 2.4 TABLE of COUNTS
table_counts <- counts_table
table_counts.raw <- counts_table.raw

# tsv table
write_tsv(x = table_counts, file = paste(OutDir, "Table_counts.txt", sep="/"), col_names=TRUE)
write_tsv(x = table_counts.raw, file = paste(OutDir, "Table_counts_RAW.txt", sep="/"), col_names=TRUE)
write_tsv(x = table_CL1_CL2 , file = paste(OutDir, "Table_CL1_CL2.txt", sep="/"), col_names=TRUE)

colnames(table_CL1_CL2)
table_CL1_CL2.SubsetUpreg <- subset(table_CL1_CL2, subset= padj<=0.05 & log2FoldChange>0.58)
write_tsv(x = table_CL1_CL2.SubsetUpreg , file = paste(OutDir, "Table_CL1_CL2_Signif.txt", sep="/"), col_names=TRUE)

###################
### VolcanoPlot
###################
library(EnhancedVolcano)
de <- na.omit(table_CL1_CL2[table_CL1_CL2$padj < 0.05, ])
infosDE <- paste("(", nrow(de), " DE genes, ", nrow(de[de$log2FoldChange > 0, ]), " up, ", nrow(de[de$log2FoldChange < 0, ]), " down)", sep="")
paste("(", nrow(de), " DE genes, ", nrow(de[de$log2FoldChange > 0, ]), " up, ", nrow(de[de$log2FoldChange < 0, ]), " down)", sep="")

pdf(file = paste0(OutDir, "/Volcano_CL1_CL2.pdf"), width = 12, height = 12)
EnhancedVolcano(resShk_df,
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.05,
                lab = resShk_df$ID,
                title = paste("CL1_CL2", infosDE),
                labSize = 2,
                drawConnectors = F,
                caption = "")
dev.off()
#################################################
## GLIMMA
#################################################;
Dataset$GeneID <- rownames(Dataset)
GlimmaTable <- Dataset

dt.res <- as.numeric(resShrink$padj<0.05)
glMDPlot(resShrink, status=dt.res, counts=counts(dds, normalized=T), html = "CL1_CL2",
         groups=colData$celltype, transform = FALSE, path = paste0(OutDir),
         folder = "CL1_CL2", main = "CL1_CL2",
         samples=colnames(counts(dds, normalized=T)), anno=GlimmaTable, display.columns = c("GeneID"))

######################################################################################################
####                                   END OF DEG ANALYSIS                                        ####
######################################################################################################