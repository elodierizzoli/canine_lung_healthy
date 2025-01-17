#############################################################################################################
#### A single cell RNA sequencing atlas of the healthy canine lung: a foundation for comparative studies ####
####                                               FIGURES                                               ####
#############################################################################################################

#### LOADING REQUIRED LIBRARY AND CUSTOM FUNCTIONS ###
library(Seurat)
library(RColorBrewer)
library(cowplot)
library(viridis)
library(tidyr)
library(dplyr)
library(sctransform)
library(ggplot2)
library(scales)
library(pheatmap)
library(gridExtra)
library(ggpubr)
library(celldex)
library(tidyverse)
library(patchwork)
library(scico)
library(ggrepel)

PercentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}

### DOTPLOT CUSTOM FUNCTION ###
# Written by GIGA Genomics Platform, University of Liège
pal <- viridis(n = 10, option = "D", direction = 1)
Dotplot <-function(x,y,z){
  Data.to.plot <- FetchData(x, vars = y)
  Data.to.plot$cell <- rownames(Data.to.plot)
  Data.to.plot$cluster <- x@active.ident
  
  Data.to.plot <- Data.to.plot %>% gather(key = genes.plot, 
                                          value = expression, -c(cell, cluster))
  Data.to.plot <- Data.to.plot %>% group_by(cluster, genes.plot) %>% 
    summarize(avg.exp = ExpMean(x = expression), pct.exp = PercentAbove(x = expression, threshold = 0), Expr = "0")
  ZSCORE<-Data.to.plot %>% group_by(genes.plot) %>% 
    summarize(M = ExpMean(x=avg.exp), SD = ExpSD(x=avg.exp))
  Data.to.plot$zscore = (Data.to.plot$avg.exp - ZSCORE$M) / ZSCORE$SD
  test <- Data.to.plot %>% group_by(genes.plot) %>% summarize(gene.expr = mean(x = avg.exp), gene.max = max(x=avg.exp))
  Data.to.plot$Ratio = Data.to.plot$avg.exp / test$gene.max
  
  G1<-ggplot(Data.to.plot, mapping = aes(x = cluster, y = genes.plot)) +
    geom_point(mapping = aes(size = pct.exp, color = zscore)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_gradientn(colours = pal, guide = "colourbar", values = c(0,0.25,0.50,0.75,1)) + 
    scale_y_discrete(limits=y) + ggtitle(z) + theme(legend.text = element_text(size=8))
  
  G2<-ggplot(Data.to.plot, mapping = aes(x = cluster, y = genes.plot)) +
    geom_point(mapping = aes(size = pct.exp, color = Ratio)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_color_gradientn(colours = pal, guide = "colourbar") + 
    scale_y_discrete(limits=y) + ggtitle(z) + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),plot.title = element_blank(),legend.text = element_text(size=8))
  G3 <-plot_grid(G1,G2,ncol=2)
  return(G2)
}

# colors
pal <- viridis(n = 10, option = "D", direction = 1)
pal13 <- viridis(n = 13, option = "D", direction = 1)

#### DEFINING PATH ###
# set working directory

# FIGURE 1 ####
Healthy_integrated <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_integrated.RDS")
Cell.Proportion <- read.table("Outs/CELL_DISTRIBUTION/Cell.Proportion.txt", header = TRUE, sep = "\t")
DefaultAssay(Healthy_integrated) <- "SCT"
P1 <- FeaturePlot(Healthy_integrated, features = c("PTPRC"),cols = c("grey", "red"),label=F, order=F, pt.size = 0.2)+
  theme(axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        legend.text = element_text(size=6),
        plot.title = element_text(size=9),
        legend.key.width = unit(0.4, "cm"))
P2 <- FeaturePlot(Healthy_integrated, features = c("EPCAM"),cols = c("grey", "red"),label=F, order=F, pt.size = 0.2)+
  theme(axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        legend.text = element_text(size=6),
        plot.title = element_text(size=9),
        legend.key.width = unit(0.4, "cm"))
P3 <- FeaturePlot(Healthy_integrated, features = c("PECAM1"),cols = c("grey", "red"),label=F, order=F, pt.size = 0.2)+
  theme(axis.title.x = element_text(size=6,vjust=38),
        axis.title.y = element_text(size=6),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        legend.text = element_text(size=6),
        plot.title = element_text(size=9),
        legend.key.width = unit(0.4, "cm"))
P4 <- DimPlot(Healthy_integrated, reduction = "umap", group.by= "compartment", label = TRUE, label.size=2, repel = TRUE, pt.size = 0.2) + 
  theme(axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6))+
  NoLegend()+ ggtitle(NULL)
P5 <- DimPlot(Healthy_integrated, reduction = "umap", group.by= "orig.ident", label = F, shuffle = TRUE, pt.size = 0.05)+ 
  theme(legend.position = "right",
        legend.text = element_text(size=6),
        axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        legend.key.width = unit(0.4, "cm"))+ 
  ggtitle(NULL)
P6 <- ggplot(data=Cell.Proportion, aes(x=Compartment, y=Counts, fill=Sample, order_by(levels(Cell.Proportion$Sample))))+ 
  geom_bar(aes(), stat="identity", position="stack", width = 0.5) + 
  scale_fill_manual(limits=c("LUNG 4","LUNG 3", "LUNG 2","LUNG 1"),values = c("#C77CFF", "#00BFC4","#7CAE00", "#F8766D"))+
  theme(legend.position = "right",legend.title = element_blank(), 
        legend.text = element_text(size=6),
        plot.tag = element_text(face="bold"), 
        axis.text.x=element_text(angle = 45, size=7, vjust=0.62),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=7),
        legend.key.width = unit(0.4, "cm"),
        panel.background = element_blank(), panel.grid.major.y = element_line(color = "grey",size=0.5)) +
  ylab("Number of cells")+xlab("")+
  scale_y_continuous(expand=c(0,0))+ 
  scale_x_discrete(limits=c("Epithelial", "Endothelial", "Immune", "Mesenchymal"))+
  guides(fill = guide_legend(reverse=T))

tiff("Fig1.tiff", width=2500, height=3000, res=400)
wrap_plots(P1,P2,P3,P4,P5, P6, byrow=F, ncol = 2)+plot_annotation(tag_levels = c("A"))& 
  theme(plot.tag = element_text(size=10))
dev.off()


# FIGURE 2 ####
Mesench_integrated <- readRDS("Outs/GLOBAL_PROCESSING/Mesench_integrated.RDS")
Idents(Mesench_integrated) <- "cellgroup"
Mesench_integrated <- RenameIdents(Mesench_integrated, "Fibroblasts"="Fibroblasts","Schwann"="Schwann cells", "Muscle"= "Muscle cells")

DefaultAssay(Mesench_integrated) <- "SCT"
DP_Mesench <- DimPlot(Mesench_integrated, reduction = "umap", label = TRUE, label.size=5, repel = TRUE, pt.size = 0.5)+
  FontSize(x.title = 12, y.title = 12, x.text=10, y.text=10)+
  NoLegend()+ ggtitle(NULL)
VlnMes <- VlnPlot(Mesench_integrated, features = c("DPT", "ASPN", "SCN7A", "NRXN1", "ACTA2", "MYH11"),
                  pt.size = 0, combine = T, stack = TRUE, flip=T, fill.by="ident", sort=F)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        strip.text.y.right = element_text(size = 12),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )+NoLegend()
wrap_plots(DP_Mesench,VlnMes, ncol = 1, byrow=T, heights=1:1)+plot_annotation(tag_levels = c("A"))

# FIGURE 3 ####
Healthy_Fib_Renamed <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_Fib.RDS")
Healthy_Fib_Renamed <- RunUMAP(Healthy_Fib_Renamed, reduction = "pca", dims = 1:5, min.dist = 0.35)
DimPlot(Healthy_Fib_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

Healthy_Fib_Renamed <- RenameIdents(Healthy_Fib_Renamed,
                      "Alveolar fibroblasts"="Alveolar",
                      "STMN2+ alveolar fibroblasts"="STMN2+ alveolar",
                      "CCBE1+ adventitial fibroblasts"="CCBE1+ adventitial",
                      "CCL19+ adventitial fibroblasts"="CCL19+ adventitial",
                      "COL23A1+ adventitial fibroblasts"="COL23A1+ adventitial",
                      "CCN3+ adventitial fibroblasts"="CCN3+ adventitial")

cols<- c("#440154FF","#95D840FF", "#3CBB75FF","#39568CFF","#F4DD02FF", "#238A8DFF")

DP_Fib <-DimPlot(Healthy_Fib_Renamed, reduction = "umap", label = T, label.box=F, repel = TRUE, 
                 cols = cols, pt.size = 1)+ 
  FontSize(x.title = 10, y.title = 10, x.text=10, y.text=10)+
  theme(legend.position = "right",
        axis.title.y = element_text(vjust=-30))+
  ggtitle(NULL)+NoLegend()

DefaultAssay(Healthy_Fib_Renamed) <- "SCT"
features <- c("NAA16","MACF1", "IFT27",
              "STMN2", "PRG4", "S100A1",
              "MED13L", "HMCN1", "CCBE1",
              "CCL19", "CCL7", "SAA1",
              "NTRK2", "STRA6","COL23A1",
              "PI3", "CCN3","MATN4")
Dotplotfib <- Dotplot(Healthy_Fib_Renamed, features, "")+
  labs(size="Fraction of cells", 
       col="Expression level")
  
wrap_plots(DP_Fib, Dotplotfib, ncol=2)+
  plot_annotation(tag_levels = c("A"))+ theme(plot.tag = element_text(face="bold"))  

# if other orientation : add +coord_flip()+  
#wrap_plots(DP_Fib, Dotplotfib, nrow=2)+
#  plot_annotation(tag_levels = c("A"))+ theme(plot.tag = element_text(face="bold"))

# FIGURE 4 ####
Healthy_Muscle_Renamed <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_Muscle.RDS")
# Order idents
new_orderMu <- c("Pulmonary pericytes","Systemic pericytes","Vascular smooth muscle cells",
                "Airway smooth muscle cells","Peribronchial myofibroblasts","Myofibroblasts")
Idents(Healthy_Muscle_Renamed) <- factor(Idents(Healthy_Muscle_Renamed), levels = new_orderMu)

DefaultAssay(Healthy_Muscle_Renamed) <- "SCT"
DP_Musc <-DimPlot(Healthy_Muscle_Renamed, reduction = "umap", label = T, 
                  label.size=2.7, repel = T, pt.size = 1, cols = cols)+ 
  theme(axis.title.x = element_text(vjust=30,size=8),
        axis.title.y = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8),
        legend.text = element_text(size=8))+
  ggtitle(NULL)+NoLegend()

featuresmusc <- c("CADM1", "POSTN", "FAM162B", 
                  "APOA1", "ADRA2A", "RGS16",
                  "PLN","ADIRF","ITIH4",
                  "ACTC1", "SCG5","DGUOK",
                  "CDH4", "FGF18", "SOSTDC1",
                  "DPT","CDO1","ARHGAP45")
Dot_musc <- Dotplot(Healthy_Muscle_Renamed, featuresmusc, "")+
  labs(size="Fraction of cells", 
       col="Expression level")
wrap_plots(DP_Musc, Dot_musc, ncol = 2)+plot_annotation(tag_levels = c("A"))+ 
  theme(plot.tag = element_text(face="bold", vjust=2.5))

# FIGURE 5 ####
Healthy_Lymphoid_Renamed <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_Lymphoid.RDS")
# Order idents
new_orderLy <- c("Plasma cells","B lymphocytes", "CD4 T","CD4 T naive", 
                 "CD4 T reg", "Th17-like T", "CD8 T", "NK T", "NK", "gd T")
Idents(Healthy_Lymphoid_Renamed) <- factor(Idents(Healthy_Lymphoid_Renamed), levels = new_orderLy)

Healthy_Myeloid_Renamed <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_Myeloid.RDS")
# Order idents
new_orderMy <- c("AM","C1Q+ AM","CCL13+ macrophages","SLC27A6+ macrophages", 
                 "FN1+ monocytes","Monocytes","CD1C+ monocytes",
                "cDC1", "cDC2","Mature DC","pDC","Neutrophils", "Mast cells")
Idents(Healthy_Myeloid_Renamed) <- factor(Idents(Healthy_Myeloid_Renamed), levels = new_orderMy)

DefaultAssay(Healthy_Lymphoid_Renamed) <- "SCT"
DP_Lym <-DimPlot(Healthy_Lymphoid_Renamed, reduction = "umap", label = T, label.size=2.7, repel = TRUE, 
                  pt.size = 0.5)+ 
  theme(axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8))+
  ggtitle(NULL)+NoLegend()

DefaultAssay(Healthy_Myeloid_Renamed) <- "SCT"
DP_My <-DimPlot(Healthy_Myeloid_Renamed, reduction = "umap", label = T, label.size=2.7, repel = TRUE, 
                 pt.size = 0.5)+ 
  theme(axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8))+
  ggtitle(NULL)+NoLegend()
VlnLy <- VlnPlot(Healthy_Lymphoid_Renamed, features = c("JCHAIN","BANK1", "CD3E", 
                                                       "LEF1",  "CTLA4", "IL23R", 
                                                        "GZMK", "NCR3", "IL17RB"),
                 pt.size = 0, combine = T, stack = TRUE, flip=T, fill.by="ident",sort=F)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=8),
        axis.text.x = element_text(size=10))+NoLegend()
VlnMy <- VlnPlot(Healthy_Myeloid_Renamed, features = c("CHI3L1", "C1QB", "CCL13","SLC27A6", 
                                                       "FN1", "IL1B", "PID1","ECRG4",    
                                                       "CD1C","CCR7", "IL3RA",
                                                        "S100A12", "CPA3"),
                 pt.size = 0, combine = T, stack = TRUE, flip=T, fill.by="ident", sort=F)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=8),
        axis.text.x = element_text(size=10))+NoLegend()

wrap_plots(DP_My, DP_Lym, VlnMy, VlnLy, ncol = 2, byrow=T, heights=1:1)+plot_annotation(tag_levels = c("A"))

# FIGURE 6 ####
Healthy_Epith_Renamed <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_Epith.RDS")
DefaultAssay(Healthy_Epith_Renamed) <- "SCT"
DP_Epith <-DimPlot(Healthy_Epith_Renamed, reduction = "umap", label = T, repel = TRUE, 
                 pt.size = 1)+ theme(legend.position = "top")+ggtitle(NULL)+NoLegend()&
  theme(axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8))
FP_Epith  <-   FeaturePlot(Healthy_Epith_Renamed, features = c("RTKN2", "SFTPC", "ITPRID1", "TP63", "IQCA1"),
                         cols = c("grey", "red"),label=F, combine = TRUE, pt.size = 0.2, ncol = 3)&
  theme(axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8),
        legend.text = element_text(size=8),
        plot.title = element_text(size=10),
        legend.key.width = unit(0.4, "cm"))

Healthy_Endo_Renamed <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_Endo.RDS")
DefaultAssay(Healthy_Endo_Renamed) <- "SCT"
DP_Endo <-DimPlot(Healthy_Endo_Renamed, reduction = "umap", label = T, repel = TRUE, 
                   pt.size = 1)+ theme(legend.position = "top")+ggtitle(NULL)+NoLegend()&
  theme(axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8))
FP_Endo <- FeaturePlot(Healthy_Endo_Renamed, features = c("RELN", "LEPR", "CEMIP2", "GJA5", "VEGFC"),
                       cols = c("grey", "red"),label=F, combine = TRUE, pt.size = 0.2, ncol = 3)&
  theme(axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8),
        legend.text = element_text(size=8),
        plot.title = element_text(size=10),
        legend.key.width = unit(0.4, "cm"))

wrap_plots(DP_Epith, FP_Epith, DP_Endo, FP_Endo, ncol = 2, widths=2:3)+plot_annotation(tag_levels = c("A"))

# SUPPL FIG 1 ####
# DOWNSIZING
Healthy_Epith_Renamed <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_Epith.RDS")
cell_counts <- table(Healthy_Epith_Renamed$origin)
min_cells <- min(cell_counts)
downsampled_cells <- unlist(lapply(names(cell_counts), function(cond) {
  cells <- WhichCells(Healthy_Epith_Renamed, expression = origin == cond)
  sample(cells, min_cells)
}))
seurat_downsampled <- subset(Healthy_Epith_Renamed, cells = downsampled_cells)
seurat_downsampled <- RunUMAP(seurat_downsampled, dims = 1:10)
# verify 
DimPlot(seurat_downsampled, group.by = "origin", split.by = "origin") +
  ggtitle("UMAP Plot with Equal Cell Numbers per Condition")

#PLOT
DefaultAssay(seurat_downsampled) <- "SCT"
FeaturePlot(seurat_downsampled, features = c("EGFR","ERBB2","PCNA", "MKI67"), 
            split.by= "origin", order = T, label=T, repel=T, pt.size = 0.7,
            cols = c("grey", "red"), combine = T)&
  theme(axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key.width = unit(0.4, "cm"))

# SUPPL FIG 2 ####
Celltype.Proportion <- read.table("Outs/CELL_DISTRIBUTION/Cell.type.Proportion.txt", header = TRUE, sep = "\t")

png("SFig2.tiff", width=2480, height=3508, res=400)
ggplot(data=Celltype.Proportion, aes(x=Counts, y=Cell.type, fill=Sample))+ 
  geom_bar(aes(y=Cell.type), stat="identity", position = position_stack(reverse = TRUE), width = 1, color="white") + 
   scale_fill_manual(limits=c("LUNG1","LUNG2","LUNG3","LUNG4"),values = c("#F8766D", "#7CAE00","#00BFC4","#C77CFF"))+
  theme(legend.position = "top",legend.title = element_blank(), 
        legend.text = element_text(size=7),
        plot.tag = element_text(face="bold"), 
        axis.text.x=element_text(size=7),
        axis.title.x=element_text(size=8),
        axis.text.y=element_text(size=7),
        axis.title.y=element_text(size=8),
        legend.key.width = unit(1, "cm"),
        panel.background = element_blank()) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_discrete(limits=rev(c("Alveolar type 1","Alveolar type 2","Secretory cells","Basal cells","Ciliated cells",
                                "Lymphatic endothelial  cells","Aerocytes","General capillary endothelial cells",
                                "Arterial endothelial cells","Venous endothelial cell","AM","C1Q+ AM","CCL13+ macrophages",
                                "SLC27A6+ macrophages","Monocytes","FN1+ monocytes","CD1C+ monocytes","Mature DC","cDC1",
                                "cDC2","Plasmacytoid DC","Neutrophils","Mast cells","CD4 T cells","CD4 naive T  cells",
                                "CD4 reg T  cells","Th17-like T cells","CD8 T cells","NK T cells","NK cells","gd T cells",
                                "B lymphocytes","Plasma cells","Alveolar fibroblasts","STMN2+ alveolar fibroblasts",
                                "CCBE1+ adventitial fibroblasts","CCL19+ adventitial fibroblasts","CCN3+ adventitial fibroblasts",
                                "COL23A1+ adventitial fibroblasts","Pulmonary pericytes","Systemic pericytes",
                                "Vascular smooth muscle cells","Airway smooth muscle cells","Peribronchial myofibroblasts",
                                "Myofibroblasts","Schwann cells")))+  
  ylab("Cell types")+xlab("% of cells")
dev.off()

# SUPPL FIG 3 ####
integ.atlas <- readRDS("Outs/HUMAN_INTEGRATION/integ.atlas2.RDS")

# S FIG 3 A # 
# renaming cell types to better match humans annotations
Idents(integ.atlas) <- "ann_finest_level"
integ.atlas_R <- RenameIdents(object = integ.atlas, 
                              "Alveolar type 1" = "AT1",
                              "Alveolar type 2" = "AT2",
                              "AM" = "Alveolar macrophages",
                              "Aerocytes" = "EC aerocyte capillary",
                              "General capillary" = "EC general capillary",
                              "Venous" = "EC venous",
                              "Arterial" = "EC arterial",
                              "NK" = "NK cells",
                              "CD8 T"= "CD8 T cells",
                              "CD4 T" = "CD4 T cells",
                              "B lymphocytes"="B cells",
                              "pDC" = "Plasmacytoid DCs",
                              "cDC1"="DC1",
                              "cDC2"="DC2",
                              "Mature DC"="Migratory DCs")
DimPlot(integ.atlas_R, reduction = "umap", group.by="ident", 
        split.by="organism", pt.size = 0.1, label=T, repel=T, label.size=3) + 
  NoLegend()
integ.atlas_R$ann_finest_level <- Idents(integ.atlas_R)

UMAP <- DimPlot(integ.atlas_R, reduction = "umap", group.by="ann_finest_level", 
        split.by="organism", pt.size = 0.1, label=T, repel=T, label.size=2.5) + 
  NoLegend()+
  FontSize(x.title = 8, y.title = 8, x.text=8, y.text=8)+
  ggtitle(NULL)

png("SFIG3A.png", width=4000, height=2000, res=400)
UMAP
dev.off()

# S FIG 3 B # 
#stash new metadata with cell type name prefixed with the specices
integ.atlas_R$type <- ifelse(grepl("Canis lupus familiaris", integ.atlas_R$organism), paste0("can_",integ.atlas_R$ann_finest_level), paste0("hu_",integ.atlas_R$ann_finest_level))    
integ.atlas_R$type <- as.factor(integ.atlas_R$type)
levels(integ.atlas_R$type) <- gsub(" ","_",levels(integ.atlas_R$type))
levels(integ.atlas_R$type)

# check
DimPlot(integ.atlas_R, reduction = "umap", group.by="type", 
        split.by="organism", pt.size = 0.1, label=T, repel=T, label.size=3) + 
  NoLegend()

# hierarchical clustering
# Method adopted from Ammons et al., 2023 'A single-cell RNA sequencing atlas of circulating leukocytes from healthy and osteosarcoma affected dogs' (DOI 10.3389/fimmu.2023.1162700)
library(ape)
library(ggtree)

# extract data
metadata <- integ.atlas_R@meta.data
expression <- as.data.frame(t(integ.atlas_R@assays$integrated@data))
expression$anno_merge <- integ.atlas_R@meta.data[rownames(expression),]$type

# average expression
avg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(avg_expression) <- avg_expression$anno_merge
avg_expression$anno_merge <- NULL

# clustering
M <- (1- cor(t(avg_expression),method="pearson"))/2

png("SFIG3B.png", width=4000, height=4000, res=400)
par(mfcol=c(1,1))
hc <- hclust(as.dist(M),method="complete")
p <- plot.phylo(as.phylo(hc), type = "fan", cex = 0.8,no.margin = FALSE,tip.col = "black")
dev.off()