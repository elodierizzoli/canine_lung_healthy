#############################################################################################################
#### A single cell RNA sequencing atlas of the healthy canine lung: a foundation for comparative studies ####
####                                       DATA PROCESSING WORKFLOW                                      ####
#############################################################################################################

#### LOADING REQUIRED LIBRARY
library(Seurat)
library(cowplot)
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
library(scico)
library(ggrepel)
library(clustree)

### GENES FOR FUTURE CELL CYCLE INVESTIGATIONS ###
ccgenes.s.genes<-toupper(c("Cdc45","Pcna","Chaf1b","Pola1","Tyms", "Gmnn","Mrpl36","Mcm5","Slbp","Rrm2",
                           "Rrm1","Uhrf1","Rad51","Rad51ap1","E2f8","Mcm4","Msh2","Wdr76","Rfc2","Ccne2","Hells",
                           "Clspn","Tipin","Fen1","Dscc1","Mcm7","Casp8ap2","Ubr7","Cdc6","Nasp","Polr1b",
                           "Prim1","Usp1","Blm","Mcm6","Gins2","Cenpu", "Cdca7","Ung","Dtl","Exo1"))
ccgenes.g2m.genes<-toupper(c("Ckap2","Tpx2","Ndc80","Tacc3","Cdca2","Rangap1","Cenpa","Ube2c","Anln","Cks2",
                             "Cdk1","Gtse1","Pimreg","Nusap1","G2e3","Ncapd2","Aurka","Cdca3","Aurkb","Dlgap5",
                             "Ckap5","Ccnb2","Kif20b","Kif11","Ttk","Cdca8","Kif23","Tubb4b","Cenpe","Cdc20",
                             "Cdc25c","Kif2c","Cbx5","Top2a","Ctcf","Bub1","Ckap2l","Mki67","Hmmr","Hmgb2",
                             "Tmpo","Psrc1","Gas2l3","Smc4","Ect2","Anp32e","Cks1b","Birc5","Nuf2","Hjurp",
                             "Nek2","Cenpf","Lbr"))

#### DEFINING PATH
# set working directory

# PROCESSING SAMPLE 1 ####
### SETUP SEURAT OBJECTS ###
LUNG1.data <- Read10X(data.dir = "Raw_data/Lung1/filtered_feature_bc_matrix/")
LUNG1 <- CreateSeuratObject(counts = LUNG1.data, project = "LUNG1", min.cells = 10, min.features = 200)

### PRE-PROCESSING WORKFLOW ###
LUNG1 <- PercentageFeatureSet(LUNG1, pattern = "^MT-", col.name = "percent.mt")
LUNG1 <- PercentageFeatureSet(LUNG1, pattern = "^RPS", col.name = "percent.RPS")
LUNG1 <- PercentageFeatureSet(LUNG1, pattern = "^RPL", col.name = "percent.RPL")

VlnPlot(LUNG1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0.02, log = F, group.by = "orig.ident")
VlnPlot(LUNG1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.02, log = T, group.by = "orig.ident")
VlnPlot(LUNG1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = F, log = T, group.by = "orig.ident")
VlnPlot(LUNG1, features = c("percent.RPS", "percent.RPL", "percent.mt"),
        ncol = 3, pt.size = F, log = F, group.by = "orig.ident")
LUNG1 <- subset(LUNG1, subset = percent.mt < 20)
VlnPlot(LUNG1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =F, log = F, group.by = "orig.ident")

### PROCESSING
LUNG1 <- SCTransform(LUNG1, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG1 <- RunPCA(LUNG1, npcs = 50, verbose = T, assay = "SCT")
ElbowPlot(LUNG1)
LUNG1 <- FindNeighbors(LUNG1, dims = 1:12, assay = "SCT")
LUNG1 <- FindClusters(LUNG1, resolution = 0.8)
LUNG1 <- RunUMAP(LUNG1, reduction = "pca", dims = 1:12, min.dist = 0.25, assay = "SCT")
DimPlot(LUNG1, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.05)

### CELL CYCLE INVESTIGATIONS ### --> NO EFFECT ON CLUSTERISATION
LUNG1 <- CellCycleScoring(LUNG1, s.features = ccgenes.s.genes, g2m.features = ccgenes.g2m.genes)
DimPlot(LUNG1, label = T, group.by = "Phase")
FeaturePlot(LUNG1, features = c("S.Score","G2M.Score"))

### COMPARTMENTS
FeaturePlot(LUNG1, reduction = "umap", features = c("EPCAM", "PECAM1", "PTPRC"), order = TRUE, label = TRUE)
VlnPlot(LUNG1, c("EPCAM", "PECAM1", "PTPRC"), log=T)

#### MESENCHYMAL ###
LUNG1_Mesench <- subset(LUNG1, subset = seurat_clusters %in% c(1,4,6,11,15,17))
LUNG1_Mesench <- SCTransform(LUNG1_Mesench, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG1_Mesench <- RunPCA(LUNG1_Mesench, verbose = T)
ElbowPlot(LUNG1_Mesench)
LUNG1_Mesench <- FindNeighbors(LUNG1_Mesench, dims = 1:8)
LUNG1_Mesench <- FindClusters(LUNG1_Mesench, resolution = 0.8)
LUNG1_Mesench <- RunUMAP(LUNG1_Mesench, reduction = "pca", dims = 1:8, min.dist = 0.35)
DimPlot(LUNG1_Mesench, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)

VlnPlot(LUNG1_Mesench, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
FeaturePlot(LUNG1_Mesench, features = c("S.Score","G2M.Score"))
FeaturePlot(LUNG1_Mesench, reduction = "umap", features = c("EPCAM", "PECAM1", "PTPRC"), order = TRUE, label = TRUE)
VlnPlot(LUNG1_Mesench, c("EPCAM", "PECAM1", "PTPRC"), log=T)
# 10 are doublets endothelial

#CELL TYPE INVESTIGATIONS
#ADVENTITIAL FIB
FeaturePlot(LUNG1_Mesench, features = c("PI16", "DCN", "COL14A1", "AQP1", "LSP1", "LTBP1"), order = T, pt.size = 1)
VlnPlot(LUNG1_Mesench, c("PI16", "DCN", "COL14A1", "AQP1", "LSP1", "LTBP1"), log=T)
FeaturePlot(LUNG1_Mesench, features = c("PI3", "MATN4", "COL14A1", "DCN", "FAP", "IL33"), order = T, pt.size = 1)
VlnPlot(LUNG1_Mesench, c("PI3", "MATN4", "COL14A1", "DCN", "FAP", "IL33"), log=T)
#VSM
FeaturePlot(LUNG1_Mesench, features = c("FABP3", "MEF2C", "PLN", "ACTA2", "TAGLN", "MYH11", "SVIL", "ITIH4"), order = T, pt.size = 1)
VlnPlot(LUNG1_Mesench, c("FABP3", "MEF2C", "PLN", "ACTA2", "TAGLN", "MYH11", "SVIL", "ITIH4"), log=T)
FeaturePlot(LUNG1_Mesench, features = c("ADRA2A", "COL12A1", "CLU", "RGS5", "COL1A1", "APOA1", "HAS2", "CP", "HGF"), order = T, pt.size = 0.5)
VlnPlot(LUNG1_Mesench, c("ADRA2A", "COL12A1", "CLU", "RGS5", "COL1A1", "APOA1", "HAS2", "CP", "HGF"), log=T)
#ASM
FeaturePlot(LUNG1_Mesench, features = c("SOSTDC1", "HHIP", "ACTC1", "PRUNE2"), order = T, pt.size = 1)
VlnPlot(LUNG1_Mesench, c("SOSTDC1", "HHIP", "ACTC1", "MYLK", "DSTN", "LPP", "DES", "PRUNE2"), log=T)
#Alv fib
FeaturePlot(LUNG1_Mesench, features = c("NPNT", "COL13A1", "MACF1", "DST", "SPECC1L", "PRG4", "STMN2"), order = T, pt.size = 1)
VlnPlot(LUNG1_Mesench, c("NPNT", "COL13A1", "MACF1", "DST", "SPECC1L", "PRG4", "STMN2"), log=T)
#Pericytes
FeaturePlot(LUNG1_Mesench, features = c("FAM162B", "POSTN", "HIGD1B", "GUCY1A1", "COX4I1"), order = T, pt.size = 1)
FeaturePlot(LUNG1_Mesench, features = c("TRPC6", "PDGFRB", "CSPG4", "MCAM"), order = T, pt.size = 1)
VlnPlot(LUNG1_Mesench, c("FAM162B", "POSTN", "HIGD1B"," COX4I1", "GUCY1A1", "TRPC6", "PDGFRB", "CSPG4", "MCAM"), log=T)
#SCHWANN CELLS - not enough to cluster
FeaturePlot(LUNG1_Mesench, features = c("SCN7A", "SNCA", "CDH19", "MPZ", "MT3"), order = T, pt.size = 1)
VlnPlot(LUNG1_Mesench, c("SCN7A", "SNCA", "CDH19", "MPZ", "MT3"), log=T)
#Mesothelial - no cluster found
FeaturePlot(LUNG1_Mesench, features = c("MSLN",	"UPK3B",	"WT1", "LRRN4", "CALB2", "CXCL13", "NKAIN4", "CDH2", "GPM6A",  "RSPO1"), order = T, pt.size = 1)
VlnPlot(LUNG1_Mesench, c("MSLN",	"UPK3B",	"WT1", "LRRN4", "CALB2", "CXCL13", "NKAIN4", "CDH2", "GPM6A",  "RSPO1"), log=T)
#PERIBRONCHIAL FIBROBLASTS
FeaturePlot(LUNG1_Mesench, features = c("ASPN", "HHIP", "FGF18", "WIF1"), order = T, pt.size = 1)
VlnPlot(LUNG1_Mesench, c("ASPN", "HHIP", "FGF18", "WIF1"), log=T)
#LIPOFIBROBLASTS - no cluster found
FeaturePlot(LUNG1_Mesench, features = c("LPL", "PLIN2", "APOE"), order = T, pt.size = 1)
VlnPlot(LUNG1_Mesench, c("LPL", "PLIN2", "APOE"), log=T)
dev.off()

#####MARKERS
Markers.SubsetMesench <- FindAllMarkers(LUNG1_Mesench, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
write.table(Markers.SubsetMesench, file = "Outs/INDIV_PROCESSING/Markers.LUNG1Mesench.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

#### Rename all identities
LUNG1_Mesench_Renamed <- RenameIdents(object = LUNG1_Mesench, 
                                             "0" = "Alveolar fibroblasts",
                                             "1" = "Alveolar fibroblasts",
                                             "2" = "Alveolar fibroblasts",
                                             "3" = "VSM",
                                             "4" = "STMN2+ alveolar fibroblasts",
                                             "5" = "Alveolar fibroblasts",
                                             "6" = "Adventitial fibroblasts",
                                             "7" = "Pericytes",
                                             "8" = "VSM",
                                             "9" = "ASM+Myofib",
                                             "10" = "Doublet",
                                             "11" = "Alveolar fibroblasts")
DimPlot(LUNG1_Mesench_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
# Remove aberrant cells
LUNG1_Mesench_Renamed <- subset(LUNG1_Mesench_Renamed, idents = "Doublet", invert = TRUE)
DimPlot(LUNG1_Mesench_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
saveRDS(LUNG1_Mesench_Renamed, "Outs/INDIV_PROCESSING/LUNG1_Mesench.RDS")

###### EPITHELIAL ###
LUNG1_Epith <- subset(LUNG1, subset = seurat_clusters %in% c(3,16))
LUNG1_Epith <- SCTransform(LUNG1_Epith, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG1_Epith <- RunPCA(LUNG1_Epith, verbose = T)
ElbowPlot(LUNG1_Epith)
LUNG1_Epith <- FindNeighbors(LUNG1_Epith, dims = 1:7)
LUNG1_Epith <- FindClusters(LUNG1_Epith, resolution = 0.8)
LUNG1_Epith <- RunUMAP(LUNG1_Epith, reduction = "pca", dims = 1:7, min.dist = 0.25)
DimPlot(LUNG1_Epith, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)

VlnPlot(LUNG1_Epith, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
FeaturePlot(LUNG1_Epith, features = c("S.Score","G2M.Score"))
FeaturePlot(LUNG1_Epith, reduction = "umap", features = c("EPCAM", "DCN", "PECAM1", "PTPRC"), order = TRUE, label = TRUE)
VlnPlot(LUNG1_Epith, c("DCN", "PECAM1", "PTPRC"), log=T)
# 5 are doublet mesench, immune and endoth, 8 and doublets immune

#CELL TYPE INVESTIGATIONS
#CLUB CELLS
FeaturePlot(LUNG1_Epith, features = c("SCGB1A1", "AQP5", "SERPINA1", "ERN2", "DEPP1", "GPX2"), order = T, pt.size = 1)
VlnPlot(LUNG1_Epith, features = c("SCGB1A1", "AQP5", "SERPINA1", "ERN2", "DEPP1", "GPX2"), log=T)
#Ionocytes, Goblet
FeaturePlot(LUNG1_Epith, features = c("CFTR", "ASCL3", "MUC5B", "SPDEF"), order = T, pt.size = 1)
VlnPlot(LUNG1_Epith, features = c("CFTR", "ASCL3", "MUC5B", "SPDEF"), log=T)
#AT1
FeaturePlot(LUNG1_Epith, features = c("AGER", "AQP1", "SUSD2", "CLDN18"), order = T, pt.size = 1)
VlnPlot(LUNG1_Epith, features = c("AGER", "AQP1", "SUSD2", "CLDN18"), log=T)
#AT2
FeaturePlot(LUNG1_Epith, features = c("NAPSA", "SFTPC", "SFTPB", "MUC1"), order = T, pt.size = 1)
VlnPlot(LUNG1_Epith, features = c("NAPSA", "SFTPC", "SFTPB", "MUC1"), log=T)
#BASAL
FeaturePlot(LUNG1_Epith, features = c("KRT14", "TP63"), order = T, pt.size = 1)
VlnPlot(LUNG1_Epith, features = c("KRT14", "TP63"), log=T)
#CILIATED
FeaturePlot(LUNG1_Epith, features = c("CAPS", "FOXJ1", "CCDC78", "HYDIN"), order = T, pt.size = 1)
VlnPlot(LUNG1_Epith, features = c("CAPS", "FOXJ1", "CCDC78", "HYDIN"), log=T)

#####MARKERS
Markers.Epith <- FindAllMarkers(LUNG1_Epith, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
write.table(Markers.Epith, file = "Outs/INDIV_PROCESSING/Markers.LUNG1Epith.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

#### Rename all identities
LUNG1_Epith_Renamed <- RenameIdents(object = LUNG1_Epith,
                                                "0" = "AT2",
                                                "1" = "AT2",
                                                "2" = "AT2",
                                                "3" = "AT2",
                                                "4" = "AT2",
                                                "5" = "Doublet", 
                                                "6" = "Secretory+Basal",
                                                "7" = "AT1",
                                                "8" = "Doublet", 
                                                "9" = "Ciliated",
                                                "10" = "AT1")
DimPlot(LUNG1_Epith_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)
# Remove aberrant cells
LUNG1_Epith_Renamed <- subset(LUNG1_Epith_Renamed, idents = "Doublet", invert = TRUE)
DimPlot(LUNG1_Epith_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
saveRDS(LUNG1_Epith_Renamed, "Outs/INDIV_PROCESSING/LUNG1_Epith.RDS")

###### IMMUNE ###
LUNG1_Immune <- subset(LUNG1, subset = seurat_clusters %in% c(2,7,9,12,13,14,18,19,20,21))
LUNG1_Immune <- SCTransform(LUNG1_Immune, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG1_Immune <- RunPCA(LUNG1_Immune, verbose = T)
ElbowPlot(LUNG1_Immune)
LUNG1_Immune <- FindNeighbors(LUNG1_Immune, dims = 1:11)
LUNG1_Immune <- FindClusters(LUNG1_Immune, resolution = 1)
LUNG1_Immune <- RunUMAP(LUNG1_Immune, reduction = "pca", dims = 1:11, min.dist = 0.25)
DimPlot(LUNG1_Immune, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)

VlnPlot(LUNG1_Immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
FeaturePlot(LUNG1_Immune, features = c("S.Score","G2M.Score")) # cl16 are replicating
FeaturePlot(LUNG1_Immune, features = c("EPCAM", "PTPRC", "DCN", "PECAM1"), order = T, pt.size = 1)
VlnPlot(LUNG1_Immune, c("EPCAM", "PECAM1","COL1A1","PTPRC"), log=T)
# cl 13 are doublets endothelial, cl 15 are doublets mesenchymal

#CELL TYPE INVESTIGATIONS
####LYMPHOCYTES
FeaturePlot(LUNG1_Immune, features = c("CD3D", "CD3E", "IL7R", "ICOS", "THY1"), order = T, pt.size = 1)
VlnPlot(LUNG1_Immune, features = c("CD3D", "CD3E", "IL7R", "ICOS", "THY1"), log=T)
####CTL
FeaturePlot(LUNG1_Immune, features = c("CD3D", "CD8A", "NCR3", "GZMB", "KLRB1", "KLRK1", "CCL5"), order = T, pt.size = 1)
VlnPlot(LUNG1_Immune, features = c("CD3D", "CD8A", "NCR3", "GZMB", "KLRB1", "KLRK1", "CCL5"), log=T)
####NAIVE LT
FeaturePlot(LUNG1_Immune, features = c("CD3D", "LEF1", "TCF7", "IL7R", "CCR7", "SELL", "LTB"), order = T, pt.size = 1)
VlnPlot(LUNG1_Immune, features = c("CD3D", "LEF1", "TCF7", "IL7R", "CCR7", "SELL", "LTB"), log=T)
####LB
FeaturePlot(LUNG1_Immune, features = c("MS4A1", "CD19", "FCRLA", "CCR7"), order = T, pt.size = 1)
VlnPlot(LUNG1_Immune, features = c("MS4A1", "CD19", "FCRLA", "CCR7"), log=T)
####PLASMA CELLS
FeaturePlot(LUNG1_Immune, features = c("TNFRSF17", "JCHAIN", "POU2AF1"), order = T, pt.size = 1)
VlnPlot(LUNG1_Immune, features = c("TNFRSF17", "JCHAIN", "POU2AF1"), log=T)
####MEMORY T - no distinct cluster
FeaturePlot(LUNG1_Immune, features = c("CD3D", "IL7R", "CD7", "SAMHD1"), order = T, pt.size = 1)
VlnPlot(LUNG1_Immune, features = c("CD3D", "IL7R", "CD7", "SAMHD1"), log=T)
####T REGULATORY - no distinct cluster
FeaturePlot(LUNG1_Immune, features = c("CD3D", "CTLA4", "IKZF2", "GATA3", "CCL5", "LGALS3"), order = T, pt.size = 1)
VlnPlot(LUNG1_Immune, features = c("CD3D", "CTLA4", "IKZF2", "GATA3", "CCL5", "LGALS3"), log=T)
#T HELPER 17/1/2
FeaturePlot(LUNG1_Immune, features = c("CD3D", "IL17RB", "RORA", "TBX21", "ALDOC", "LGALS3"), order = T, pt.size = 1)
VlnPlot(LUNG1_Immune, features = c("CD3D", "IL17RB", "RORA", "TBX21", "ALDOC", "LGALS3"), log=T)
##MAST CELLS
FeaturePlot(LUNG1_Immune, features = c("KIT", "CPA3", "MS4A2", "FCER1A"), order = T, pt.size = 1)
VlnPlot(LUNG1_Immune, features = c("KIT", "CPA3", "MS4A2", "FCER1A"), log=T)
###NEUTRO
FeaturePlot(LUNG1_Immune, features = c("S100A12", "ITGAM", "CXCL8", "ADGRE1", "CD4", "SELL", "PADI3"), order = T, pt.size = 1)
VlnPlot(LUNG1_Immune, features = c("S100A12", "ITGAM", "CXCL8", "ADGRE1", "CD4", "SELL", "PADI3"), log=T)
##GENERAL MP
FeaturePlot(LUNG1_Immune, features = c("MRC1", "MSR1", "APOC1", "CD86", "CD48", "CD63"), order = T, pt.size = 1)
VlnPlot(LUNG1_Immune, features = c("MRC1", "MSR1", "APOC1", "CD86", "CD48", "CD63"), log=T)
#AM
FeaturePlot(LUNG1_Immune, features = c("MARCO", "MSR1","MRC1", "CD163", "CHI3L1", "APOC1"), order = T, pt.size = 1)
VlnPlot(LUNG1_Immune, features = c("MARCO", "MSR1","MRC1", "CD163", "CHI3L1", "APOC1"), log=T)
##C1QC MP / IM
FeaturePlot(LUNG1_Immune, features = c("C1QC", "C1QA", "CCL13", "STAB1", "NPAS1", "MRC1"), order = T, pt.size = 1)
VlnPlot(LUNG1_Immune, features = c("C1QC", "C1QA", "CCL13", "STAB1", "NPAS1", "MRC1"), log=T)
#MONOCYTES
FeaturePlot(LUNG1_Immune, features = c("LYZ", "PLEK", "VCAN", "FN1", "IL1B", "MRC1"), order = T, pt.size = 1)
VlnPlot(LUNG1_Immune, features = c("LYZ", "PLEK", "VCAN", "FN1", "IL1B", "MRC1"), log=T)
#DC
FeaturePlot(LUNG1_Immune, features = c("CD1C", "PLD4", "FCER1A",  "PPM1J", "IRF8", "DLA-DOA", "PKIB"), order = T, pt.size = 1)
VlnPlot(LUNG1_Immune, features = c("CD1C", "PLD4", "FCER1A",  "PPM1J", "IRF8", "DLA-DOA", "PKIB"), log=T)
#MATURE DC
FeaturePlot(LUNG1_Immune, features = c("CD1C","CD83", "IDO1", "IL4I1", "CCR7"), order = T, pt.size = 1)
VlnPlot(LUNG1_Immune, features = c("CD1C","CD83", "IDO1", "IL4I1", "CCR7"), log=T)

###MARKERS
Markers.Immune <- FindAllMarkers(LUNG1_Immune, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
write.table(Markers.Immune, file = "Outs/INDIV_PROCESSING/Markers.LUNG1Immune.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

## Rename all identities
LUNG1_Immune_Renamed <- RenameIdents(object = LUNG1_Immune, 
                                                 "0" = "AM", "1" = "AM", "2" = "T lymphocytes + NK",
                                                 "3" = "Mast cells", "4" = "IM", "5" = "Mast cells",
                                                 "6" = "Neutrophils", "7" = "AM", "8" = "DC",
                                                 "9" = "MMP12+ monocytes", "10" = "Monocytes", "11" = "T lymphocytes + NK", 
                                                 "12" = "T lymphocytes", "13" = "Doublet", "14" = "B+Plasma cells", 
                                                "15" = "Doublet", "16" = "Replicating immune")
DimPlot(LUNG1_Immune_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)
# Remove aberrant cells
LUNG1_Immune_Renamed <- subset(LUNG1_Immune_Renamed, idents = c("Doublet"), invert = TRUE)
DimPlot(LUNG1_Immune_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
saveRDS(LUNG1_Immune_Renamed, "Outs/INDIV_PROCESSING/LUNG1_Immune.RDS")

###### ENDOTHELIAL ###
LUNG1_Endo <- subset(LUNG1, subset = seurat_clusters %in% c(0,5,8,10))
LUNG1_Endo <- SCTransform(LUNG1_Endo, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG1_Endo <- RunPCA(LUNG1_Endo, verbose = T)
ElbowPlot(LUNG1_Endo)
LUNG1_Endo <- FindNeighbors(LUNG1_Endo, dims = 1:10)
LUNG1_Endo <- FindClusters(LUNG1_Endo, resolution = 0.4)
LUNG1_Endo <- RunUMAP(LUNG1_Endo, reduction = "pca", dims = 1:10, min.dist = 0.25)
DimPlot(LUNG1_Endo, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)

VlnPlot(LUNG1_Endo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
FeaturePlot(LUNG1_Endo, features = c("S.Score","G2M.Score"))
FeaturePlot(LUNG1_Endo, features = c("EPCAM", "PTPRC", "PECAM1", "DCN"), order = T, pt.size = 1)
VlnPlot(LUNG1_Endo, c("EPCAM", "PECAM1","DCN","PTPRC"), log=T)

#####MARKERS
Markers.Endo <- FindAllMarkers(LUNG1_Endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
write.table(Markers.Endo, file = "Outs/INDIV_PROCESSING/Markers.LUNG1Endo.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

#CELL TYPE INVESTIGATIONS
#ENDO GENERAL
FeaturePlot(LUNG1_Endo, features = c("PECAM1", "LYVE1", "CDH5", "VWF", "CD34"), order = T, pt.size = 1)
VlnPlot(LUNG1_Endo, features = c("PECAM1", "LYVE1", "CDH5", "VWF", "CD34"), log=T)
##LYMPHATIC
FeaturePlot(LUNG1_Endo, features = c("VWF", "PROX1", "PDPN", "TBX1", "THY1", "FLT4"), order = T, pt.size = 1)
VlnPlot(LUNG1_Endo, features = c("VWF", "PROX1", "PDPN", "TBX1", "THY1", "FLT4"), log=T)
#ARTERY EC
FeaturePlot(LUNG1_Endo, features = c("BMX", "MGP", "GJA5", "EFNB2"), order = T, pt.size = 1)
VlnPlot(LUNG1_Endo, features = c("BMX", "MGP", "GJA5", "EFNB2"), log=T)
##VEINOUS EC
FeaturePlot(LUNG1_Endo, features = c("SELP", "SELE", "VCAM1", "TIMP1", "ACKR1"), order = T, pt.size = 1)
VlnPlot(LUNG1_Endo, features = c("SELP", "SELE", "VCAM1", "TIMP1", "ACKR1"), log=T)
#CAPILLARY EC
FeaturePlot(LUNG1_Endo, features = c("CA4", "SPARC", "SGK1", "EMP2"), order = T, pt.size = 1)
VlnPlot(LUNG1_Endo, features = c("CA4", "SPARC", "SGK1", "EMP2"), log=T)
##AEROCYTE
FeaturePlot(LUNG1_Endo, features = c("VWF", "EMCN", "EDNRB", "TBX2"), order = T, pt.size = 1)
VlnPlot(LUNG1_Endo, features = c("VWF", "EMCN", "EDNRB", "TBX2"), log=T)
##GENERAL CAPILLARY
FeaturePlot(LUNG1_Endo, features = c("EDN1", "GPIHBP1", "CD36"), order = T, pt.size = 1)
VlnPlot(LUNG1_Endo, features = c("EDN1", "GPIHBP1", "CD36"), log=T)

#### Rename all identities
LUNG1_Endo_Renamed <- RenameIdents(object = LUNG1_Endo, 
                                               "0" = "General capillary EC",
                                               "1" = "Aerocytes",
                                               "2" = "General capillary EC",
                                               "3" = "General capillary EC",
                                               "4" = "Artery EC",
                                               "5" = "Aerocytes",
                                               "6" = "Vein EC",
                                               "7" = "Doublet",
                                               "8" = "Lymphatic EC")
DimPlot(LUNG1_Endo_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)
# Remove aberrant cells
LUNG1_Endo_Renamed <- subset(LUNG1_Endo_Renamed, idents = c("Doublet"), invert = TRUE)
DimPlot(LUNG1_Endo_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
saveRDS(LUNG1_Endo_Renamed, "Outs/INDIV_PROCESSING/LUNG1_Endo.RDS")

###FUSION OF COMPARTMENTS
LUNG1_Mesench_Renamed <- readRDS("Outs/INDIV_PROCESSING/LUNG1_Mesench.RDS")
LUNG1_Epith_Renamed <- readRDS("Outs/INDIV_PROCESSING/LUNG1_Epith.RDS")
LUNG1_Immune_Renamed <- readRDS("Outs/INDIV_PROCESSING/LUNG1_Immune.RDS")
LUNG1_Endo_Renamed <- readRDS("Outs/INDIV_PROCESSING/LUNG1_Endo.RDS")
LUNG1_integrated <- merge(LUNG1_Endo_Renamed, y = c(LUNG1_Mesench_Renamed, LUNG1_Epith_Renamed, LUNG1_Immune_Renamed), 
                                 add.cell.ids = c("Endothelial", "Mesenchymal", "Epithelial", "Immune"), project = "LUNG1")
LUNG1_integrated <- SCTransform(LUNG1_integrated, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG1_integrated <- RunPCA(LUNG1_integrated, npcs = 50, verbose = T, assay = "SCT")
LUNG1_integrated <- RunUMAP(LUNG1_integrated, reduction = "pca", dims = 1:50, min.dist = 0.25, assay = "SCT")
saveRDS(LUNG1_integrated, "Outs/INDIV_PROCESSING/LUNG1_integrated.RDS")
DimPlot(LUNG1_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.05)+NoLegend()
# end of sample 1

# PROCESSING SAMPLE 2 ####
### SETUP SEURAT OBJECTS ###
LUNG2.data <- Read10X(data.dir = "Raw_data/Lung2/filtered_feature_bc_matrix/")
LUNG2 <- CreateSeuratObject(counts = LUNG2.data, project = "LUNG2", min.cells = 10, min.features = 200)

### PRE-PROCESSING WORKFLOW ###
LUNG2 <- PercentageFeatureSet(LUNG2, pattern = "^MT-", col.name = "percent.mt")
LUNG2 <- PercentageFeatureSet(LUNG2, pattern = "^RPS", col.name = "percent.RPS")
LUNG2 <- PercentageFeatureSet(LUNG2, pattern = "^RPL", col.name = "percent.RPL")

VlnPlot(LUNG2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0.02, log = F, group.by = "orig.ident")
VlnPlot(LUNG2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.02, log = T, group.by = "orig.ident")
VlnPlot(LUNG2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = F, log = T, group.by = "orig.ident")
VlnPlot(LUNG2, features = c("percent.RPS", "percent.RPL", "percent.mt"),
        ncol = 3, pt.size = F, log = F, group.by = "orig.ident")
LUNG2 <- subset(LUNG2, subset = percent.mt < 20)
VlnPlot(LUNG2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =F, log = F, group.by = "orig.ident")

### PROCESSING
LUNG2 <- SCTransform(LUNG2, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG2 <- RunPCA(LUNG2, npcs = 50, verbose = T, assay = "SCT")
ElbowPlot(LUNG2)
LUNG2 <- FindNeighbors(LUNG2, dims = 1:10, assay = "SCT")
LUNG2 <- FindClusters(LUNG2, resolution = 0.8)
LUNG2 <- RunUMAP(LUNG2, reduction = "pca", dims = 1:10, min.dist = 0.25, assay = "SCT")
DimPlot(LUNG2, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.05)

### CELL CYCLE INVESTIGATIONS ### NO EFFECT ON CLUSTERING
LUNG2 <- CellCycleScoring(LUNG2, s.features = ccgenes.s.genes, g2m.features = ccgenes.g2m.genes)
DimPlot(LUNG2, label = T, group.by = "Phase")
FeaturePlot(LUNG2, features = c("S.Score","G2M.Score"))

####### COMPARTMENTS
FeaturePlot(LUNG2, reduction = "umap", features = c("EPCAM", "PECAM1", "PTPRC"), order = TRUE, label = TRUE)
VlnPlot(LUNG2, c("EPCAM", "PECAM1", "PTPRC"), log=T)

######## MESENCHYMAL ### 
LUNG2_Mesench <- subset(LUNG2, subset = seurat_clusters %in% c(6,7,9,10,13,14))
LUNG2_Mesench <- SCTransform(LUNG2_Mesench, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG2_Mesench <- RunPCA(LUNG2_Mesench, npcs = 50, verbose = T)
ElbowPlot(LUNG2_Mesench)
LUNG2_Mesench <- FindNeighbors(LUNG2_Mesench, dims = 1:10)
LUNG2_Mesench <- FindClusters(LUNG2_Mesench, resolution = 0.6)
LUNG2_Mesench <- RunUMAP(LUNG2_Mesench, reduction = "pca", dims = 1:10, min.dist = 0.35)
DimPlot(LUNG2_Mesench, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)

VlnPlot(LUNG2_Mesench, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
FeaturePlot(LUNG2_Mesench, features = c("S.Score","G2M.Score"))
FeaturePlot(LUNG2_Mesench, reduction = "umap", features = c("EPCAM", "PECAM1", "PTPRC"), order = TRUE, label = TRUE)
#cluster 7 are doublets with neutrophils

# CELL TYPE INVESTIGATIONS
#ADVENTITIAL FIB
FeaturePlot(LUNG2_Mesench, features = c("PI16", "DCN", "COL14A1", "AQP1", "LSP1", "LTBP1"), order = T, pt.size = 1)
VlnPlot(LUNG2_Mesench, c("PI16", "DCN", "COL14A1", "AQP1", "LSP1", "LTBP1"), log=T)
FeaturePlot(LUNG2_Mesench, features = c("PI3", "MATN4", "COL14A1", "DCN", "FAP", "IL33"), order = T, pt.size = 1)
VlnPlot(LUNG2_Mesench, c("PI3", "MATN4", "COL14A1", "DCN", "FAP", "IL33"), log=T)
#VSM
FeaturePlot(LUNG2_Mesench, features = c("FABP3", "MEF2C", "PLN", "ACTA2", "TAGLN", "MYH11", "SVIL", "ITIH4"), order = T, pt.size = 1)
VlnPlot(LUNG2_Mesench, c("FABP3", "MEF2C", "PLN", "ACTA2", "TAGLN", "MYH11", "SVIL", "ITIH4"), log=T)
FeaturePlot(LUNG2_Mesench, features = c("ADRA2A", "COL12A1", "CLU", "RGS5", "COL1A1", "APOA1", "HAS2", "CP", "HGF"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Mesench, c("ADRA2A", "COL12A1", "CLU", "RGS5", "COL1A1", "APOA1", "HAS2", "CP", "HGF"), log=T)
FeaturePlot(LUNG2_Mesench, features = c("FAM162B"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Mesench, c("ADRA2A", "COL12A1", "CLU", "RGS5", "COL1A1", "APOA1", "HAS2", "CP", "HGF"), log=T)
#ASM
FeaturePlot(LUNG2_Mesench, features = c("SOSTDC1", "HHIP", "ACTC1", "PRUNE2"), order = T, pt.size = 1)
VlnPlot(LUNG2_Mesench, c("SOSTDC1", "HHIP", "ACTC1", "MYLK", "DSTN", "LPP", "DES", "PRUNE2"), log=T)
#Alv fib
FeaturePlot(LUNG2_Mesench, features = c("NPNT", "COL13A1", "MACF1", "DST", "SPECC1L", "PRG4", "STMN2"), order = T, pt.size = 1)
VlnPlot(LUNG2_Mesench, c("NPNT", "COL13A1", "MACF1", "DST", "SPECC1L", "PRG4", "STMN2"), log=T)

#####MARKERS OF NEW CLUSTERS 
Markers.SubsetMesench <- FindAllMarkers(LUNG2_Mesench, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
write.table(Markers.SubsetMesench, file = "Outs/INDIV_PROCESSING/Markers.LUNG2Mesench.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

#### Rename all identities
LUNG2_Mesench_Renamed <- RenameIdents(object = LUNG2_Mesench, 
                                      "0" = "Alveolar fibroblasts",
                                      "1" = "Alveolar fibroblasts",
                                      "2" = "VSM+ASM+SystPericytes",
                                      "3" = "Adventitial fibroblasts",
                                      "4" = "Pericytes",
                                      "5" = "STMN2+ alveolar fibroblasts",
                                      "6" = "VSM",
                                      "7" = "Doublet",
                                      "8" = "Adventitial fibroblasts",
                                      "9" = "Adventitial fibroblasts",
                                      "10" = "Alveolar fibroblasts")
DimPlot(LUNG2_Mesench_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)
# Remove aberrant cells
LUNG2_Mesench_Renamed <- subset(LUNG2_Mesench_Renamed, idents = "Doublet", invert = TRUE)
DimPlot(LUNG2_Mesench_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
saveRDS(LUNG2_Mesench_Renamed, "Outs/INDIV_PROCESSING/LUNG2_Mesench.RDS")

###### EPITHELIAL ###
LUNG2_Epith <- subset(LUNG2, subset = seurat_clusters %in% c(3))
LUNG2_Epith <- SCTransform(LUNG2_Epith, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG2_Epith <- RunPCA(LUNG2_Epith, verbose = T)
ElbowPlot(LUNG2_Epith)
LUNG2_Epith <- FindNeighbors(LUNG2_Epith, dims = 1:9)
LUNG2_Epith <- FindClusters(LUNG2_Epith, resolution = 0.8)
LUNG2_Epith <- RunUMAP(LUNG2_Epith, reduction = "pca", dims = 1:9, min.dist = 0.25)
DimPlot(LUNG2_Epith, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)

### QC
VlnPlot(LUNG2_Epith, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
FeaturePlot(LUNG2_Epith, features = c("percent.mt"))
FeaturePlot(LUNG2_Epith, features = c("S.Score","G2M.Score"))
FeaturePlot(LUNG2_Epith, reduction = "umap", features = c("EPCAM", "PECAM1", "PTPRC", "DCN"), order = TRUE, label = TRUE)
# cl 3 is doublet with fibroblasts, 4 is doublet immune and endothelial, 5 is doublet macrophage

# CELL TYPE INVESTIGATIONS
FeaturePlot(LUNG2_Epith, features = c("SFTPC", "AGER", "SCGB1A1", "GPX2", "MUC5B", "HYDIN", "KRT14", "TP63", "CFTR"), order = T, label=T, pt.size = 1)
VlnPlot(LUNG2_Epith, features = c("SFTPC", "AGER", "SCGB1A1", "GPX2", "MUC5B", "HYDIN", "KRT14", "TP63", "CFTR"), log=T)

#####MARKERS OF NEW CLUSTERS 
Markers.Epith <- FindAllMarkers(LUNG2_Epith, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
write.table(Markers.Epith, file = "Outs/INDIV_PROCESSING/Markers.LUNG2Epith.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

#### Rename all identities
LUNG2_Epith_Renamed <- RenameIdents(object = LUNG2_Epith, 
                                           "0" = "AT1",
                                           "1" = "Secretory",
                                           "2" = "AT2",
                                           "3" = "Doublet",
                                           "4" = "Doublet",
                                           "5" = "Doublet",
                                          "6" = "Basal")
DimPlot(LUNG2_Epith_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)
# Remove aberrant cells
LUNG2_Epith_Renamed <- subset(LUNG2_Epith_Renamed, idents = "Doublet", invert = TRUE)
DimPlot(LUNG2_Epith_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)
saveRDS(LUNG2_Epith_Renamed, "Outs/INDIV_PROCESSING/LUNG2_Epith.RDS")

###### IMMUNE ###
LUNG2_Immune <- subset(LUNG2, subset = seurat_clusters %in% c(4,5,8,11,15,16,12))
LUNG2_Immune <- subset(LUNG2_Immune, subset = percent.mt < 15)
LUNG2_Immune <- SCTransform(LUNG2_Immune, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG2_Immune <- RunPCA(LUNG2_Immune, verbose = T)
ElbowPlot(LUNG2_Immune)
LUNG2_Immune <- FindNeighbors(LUNG2_Immune, dims = 1:10)
LUNG2_Immune <- FindClusters(LUNG2_Immune, resolution = 1.5)
LUNG2_Immune <- RunUMAP(LUNG2_Immune, reduction = "pca", dims = 1:10, min.dist = 0.25)
DimPlot(LUNG2_Immune, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)

#QC
VlnPlot(LUNG2_Immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
FeaturePlot(LUNG2_Immune, features = c("nFeature_RNA", "percent.mt"), ncol = 3)
VlnPlot(LUNG2_Immune, features = c("percent.RPS", "percent.RPL", "percent.mt"),ncol = 3, pt.size = F, log = F)
FeaturePlot(LUNG2_Immune, features = c("S.Score","G2M.Score")) # no effect of cell cycle
#DOUBLETS?
FeaturePlot(LUNG2_Immune, features = c("EPCAM", "DCN", "PECAM1"), order = T, pt.size = 1, label=T)
VlnPlot(LUNG2_Immune, c("EPCAM", "PECAM1", "DCN"), log=T)
# cl 4,7,10,11,16,19 are doublets or low quality cells. 14 are stressed/dying cells + mast cells

# CELL TYPE INVESTIGATIONS
####LYMPHOCYTES
FeaturePlot(LUNG2_Immune, features = c("CD3D", "CD3E", "IL7R", "ICOS", "THY1"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Immune, features = c("CD3D", "CD3E", "IL7R", "ICOS", "THY1"), log=T)
####CTL
FeaturePlot(LUNG2_Immune, features = c("CD3D", "CD8A", "NCR3", "GZMB", "KLRB1", "KLRK1", "CCL5"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Immune, features = c("CD3D", "CD8A", "NCR3", "GZMB", "KLRB1", "KLRK1", "CCL5"), log=T)
####NAIVE LT
FeaturePlot(LUNG2_Immune, features = c("CD3D", "LEF1", "TCF7", "IL7R", "CCR7", "SELL", "LTB"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Immune, features = c("CD3D", "LEF1", "TCF7", "IL7R", "CCR7", "SELL", "LTB"), log=T)
####LB
FeaturePlot(LUNG2_Immune, features = c("MS4A1", "CD19", "FCRLA", "CCR7"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Immune, features = c("MS4A1", "CD19", "FCRLA", "CCR7"), log=T)
####PLASMA CELLS
FeaturePlot(LUNG2_Immune, features = c("TNFRSF17", "JCHAIN", "POU2AF1"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Immune, features = c("TNFRSF17", "JCHAIN", "POU2AF1"), log=T)
####MEMORY T
FeaturePlot(LUNG2_Immune, features = c("CD3D", "IL7R", "CD7", "SAMHD1"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Immune, features = c("CD3D", "IL7R", "CD7", "SAMHD1"), log=T)
####T REGULATORY
FeaturePlot(LUNG2_Immune, features = c("CD3D", "CTLA4", "IKZF2", "GATA3", "CCL5", "LGALS3"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Immune, features = c("CD3D", "CTLA4", "IKZF2", "GATA3", "CCL5", "LGALS3"), log=T)
#T HELPER 17/1/2
FeaturePlot(LUNG2_Immune, features = c("CD3D", "IL17RB", "RORA", "TBX21", "ALDOC", "LGALS3"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Immune, features = c("CD3D", "IL17RB", "RORA", "TBX21", "ALDOC", "LGALS3"), log=T)
##MAST CELLS
FeaturePlot(LUNG2_Immune, features = c("KIT", "CPA3", "MS4A2", "FCER1A"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Immune, features = c("KIT", "CPA3", "MS4A2", "FCER1A"), log=T)
###NEUTRO
FeaturePlot(LUNG2_Immune, features = c("S100A12", "ITGAM", "CXCL8", "ADGRE1", "CD4", "SELL", "PADI3"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Immune, features = c("S100A12", "ITGAM", "CXCL8", "ADGRE1", "CD4", "SELL", "PADI3"), log=T)
##GENERAL MP
FeaturePlot(LUNG2_Immune, features = c("MRC1", "MSR1", "APOC1", "CD86", "CD48", "CD63"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Immune, features = c("MRC1", "MSR1", "APOC1", "CD86", "CD48", "CD63"), log=T)
#AM
FeaturePlot(LUNG2_Immune, features = c("MARCO", "MSR1","MRC1", "CD163", "CHI3L1", "APOC1"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Immune, features = c("MARCO", "MSR1","MRC1", "CD163", "CHI3L1", "APOC1"), log=T)
##C1QC MP
FeaturePlot(LUNG2_Immune, features = c("C1QC", "C1QA", "CCL13", "STAB1", "NPAS1", "MRC1"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Immune, features = c("C1QC", "C1QA", "CCL13", "STAB1", "NPAS1", "MRC1"), log=T)
#MONOCYTES
FeaturePlot(LUNG2_Immune, features = c("LYZ", "PLEK", "VCAN", "FN1", "IL1B", "MRC1"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Immune, features = c("LYZ", "PLEK", "VCAN", "FN1", "IL1B", "MRC1"), log=T)
#DC
FeaturePlot(LUNG2_Immune, features = c("CD1C", "PLD4", "FCER1A",  "PPM1J", "IRF8", "DLA-DOA", "PKIB"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Immune, features = c("CD1C", "PLD4", "FCER1A",  "PPM1J", "IRF8", "DLA-DOA", "PKIB"), log=T)
#MATURE DC
FeaturePlot(LUNG2_Immune, features = c("CD1C","CD83", "IDO1", "IL4I1", "CCR7"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Immune, features = c("CD1C","CD83", "IDO1", "IL4I1", "CCR7"), log=T)
#REPLICATING
FeaturePlot(LUNG2_Immune, features = c("CENPF", "MKI67", "TOP2A"), order = T, pt.size = 0.5)
VlnPlot(LUNG2_Immune, features = c("CENPF", "MKI67", "TOP2A"), log=T)
dev.off()

###MARKERS OF NEW CLUSTERS 
Markers.Immune <- FindAllMarkers(LUNG2_Immune, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
write.table(Markers.Immune, file = "Outs/INDIV_PROCESSING/Markers.LUNG2Immune.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

## Rename all identities
LUNG2_Immune_Renamed <- RenameIdents(object = LUNG2_Immune, 
                                     "0" = "Neutrophils", "1" = "CD4 T", "2" = "CD4 T", "3" = "AM + IM", 
                                     "4" = "Doublet", "5" = "Cytotoxic T", "6" = "Monocytes", 
                                     "7" = "Doublet", "8" = "Monocytes", "9" = "FN1+ monocytes",
                                     "10" = "Doublet", "11" = "Doublet", "12" = "DC",
                                     "13" = "CD4 T", "14" = "Mast cells + Dying cells", "15" = "Neutrophils",
                                     "16" = "Doublet", "17" = "B+Plasma cells", "18" = "DC", 
                                     "19" = "Doublet")
DimPlot(LUNG2_Immune_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)
# Remove aberrant cells
LUNG2_Immune_Renamed <- subset(LUNG2_Immune_Renamed, idents = c("Doublet"), invert = TRUE)
DimPlot(LUNG2_Immune_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1.5)
saveRDS(LUNG2_Immune_Renamed, "Outs/INDIV_PROCESSING/LUNG2_Immune.RDS")

###### ENDOTHELIAL ###
LUNG2_Endo <- subset(LUNG2, subset = seurat_clusters %in% c(0,1,2,17))
LUNG2_Endo <- SCTransform(LUNG2_Endo, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG2_Endo <- RunPCA(LUNG2_Endo, verbose = T)
ElbowPlot(LUNG2_Endo)
LUNG2_Endo <- FindNeighbors(LUNG2_Endo, dims = 1:8)
LUNG2_Endo <- FindClusters(LUNG2_Endo, resolution = 0.4)
LUNG2_Endo <- RunUMAP(LUNG2_Endo, reduction = "pca", dims = 1:8, min.dist = 0.25)
DimPlot(LUNG2_Endo, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)

#QC
VlnPlot(LUNG2_Endo, features = c("percent.RPS", "percent.RPL", "percent.mt"), ncol = 3, pt.size = F, log = T)
VlnPlot(LUNG2_Endo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
FeaturePlot(LUNG2_Endo, features = c("EPCAM", "DCN", "PTPRC"), order = T, label=T, pt.size = 1)
# 0 and 4 are doublets

# CELL TYPE INVESTIGATION
#ENDO GENERAL
FeaturePlot(LUNG2_Endo, features = c("PECAM1", "LYVE1", "CDH5", "VWF", "CD34"), order = T, pt.size = 1)
VlnPlot(LUNG2_Endo, features = c("PECAM1", "LYVE1", "CDH5", "VWF", "CD34"), log=T)
##LYMPHATIC
FeaturePlot(LUNG2_Endo, features = c("VWF", "PROX1", "PDPN", "TBX1", "THY1", "FLT4"), order = T, pt.size = 1)
VlnPlot(LUNG2_Endo, features = c("VWF", "PROX1", "PDPN", "TBX1", "THY1", "FLT4"), log=T)
#ARTERY EC
FeaturePlot(LUNG2_Endo, features = c("BMX", "MGP", "GJA5", "EFNB2"), order = T, pt.size = 1)
VlnPlot(LUNG2_Endo, features = c("BMX", "MGP", "GJA5", "EFNB2"), log=T)
##VEINOUS EC
FeaturePlot(LUNG2_Endo, features = c("SELP", "SELE", "VCAM1", "TIMP1", "ACKR1"), order = T, pt.size = 1)
VlnPlot(LUNG2_Endo, features = c("SELP", "SELE", "VCAM1", "TIMP1", "ACKR1"), log=T)
#CAPILLARY EC
FeaturePlot(LUNG2_Endo, features = c("CA4", "SPARC", "SGK1", "EMP2"), order = T, pt.size = 1)
VlnPlot(LUNG2_Endo, features = c("CA4", "SPARC", "SGK1", "EMP2"), log=T)
##AEROCYTE
FeaturePlot(LUNG2_Endo, features = c("VWF", "EMCN", "EDNRB", "TBX2"), order = T, pt.size = 1)
VlnPlot(LUNG2_Endo, features = c("VWF", "EMCN", "EDNRB", "TBX2"), log=T)
##GENERAL CAPILLARY
FeaturePlot(LUNG2_Endo, features = c("EDN1", "GPIHBP1", "CD36"), order = T, pt.size = 1)
VlnPlot(LUNG2_Endo, features = c("EDN1", "GPIHBP1", "CD36"), log=T)

#####MARKERS OF NEW CLUSTERS 
Markers.Endo <- FindAllMarkers(LUNG2_Endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
write.table(Markers.Endo, file = "Outs/INDIV_PROCESSING/Markers.LUNG2Endo.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

#### Rename all identities
LUNG2_Endo_Renamed <- RenameIdents(object = LUNG2_Endo, 
                                          "0" = "Doublet",
                                          "1" = "Aerocytes",
                                          "2" = "GC + Artery EC",
                                          "3" = "General capillary EC",
                                          "4" = "Doublet", 
                                          "5" = "Vein EC",
                                          "6" = "Lymphatic EC")
DimPlot(LUNG2_Endo_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)
# Remove aberrant cells
LUNG2_Endo_Renamed <- subset(LUNG2_Endo_Renamed, idents = c("Doublet"), invert = TRUE)
DimPlot(LUNG2_Endo_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
saveRDS(LUNG2_Endo_Renamed, "Outs/INDIV_PROCESSING/LUNG2_Endo.RDS")

###FUSION OF COMPARTMENTS
LUNG2_integrated <- merge(LUNG2_Endo_Renamed, y = c(LUNG2_Mesench_Renamed, LUNG2_Epith_Renamed, LUNG2_Immune_Renamed), 
                            add.cell.ids = c("Endothelial", "Mesenchymal", "Epithelial", "Immune"), project = "LUNG2")
LUNG2_integrated <- SCTransform(LUNG2_integrated, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG2_integrated <- RunPCA(LUNG2_integrated, npcs = 50, verbose = T, assay = "SCT")
LUNG2_integrated <- RunUMAP(LUNG2_integrated, reduction = "pca", dims = 1:10, min.dist = 0.25, assay = "SCT")
saveRDS(LUNG2_integrated, "Outs/INDIV_PROCESSING/LUNG2_integrated.RDS")
DimPlot(LUNG2_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.05)+NoLegend()
# end of sample 2

# PROCESSING SAMPLE 3 ####
### SETUP SEURAT OBJECTS ###
LUNG3.data <- Read10X(data.dir = "Raw_data/Lung3/filtered_feature_bc_matrix/")
LUNG3 <- CreateSeuratObject(counts = LUNG3.data, project = "LUNG3",min.cells = 10, min.features = 200)

### PRE-PROCESSING WORKFLOW ###
LUNG3 <- PercentageFeatureSet(LUNG3, pattern = "^MT-", col.name = "percent.mt")
LUNG3 <- PercentageFeatureSet(LUNG3, pattern = "^RPS", col.name = "percent.RPS")
LUNG3 <- PercentageFeatureSet(LUNG3, pattern = "^RPL", col.name = "percent.RPL")

VlnPlot(LUNG3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0.02, log = F, group.by = "orig.ident")
VlnPlot(LUNG3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.02, log = T, group.by = "orig.ident")
VlnPlot(LUNG3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = F, log = T, group.by = "orig.ident")
VlnPlot(LUNG3, features = c("percent.RPS", "percent.RPL", "percent.mt"),
        ncol = 3, pt.size = F, log = F, group.by = "orig.ident")
LUNG3 <- subset(LUNG3, subset = percent.mt < 20)
VlnPlot(LUNG3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =F, log = F, group.by = "orig.ident")

#PROCESSING
LUNG3 <- SCTransform(LUNG3, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG3 <- RunPCA(LUNG3, npcs = 50, verbose = T, assay = "SCT")
ElbowPlot(LUNG3)
LUNG3 <- FindNeighbors(LUNG3, dims = 1:17, assay = "SCT")
LUNG3 <- FindClusters(LUNG3, resolution = 0.8)
LUNG3 <- RunUMAP(LUNG3, reduction = "pca", dims = 1:17, min.dist = 0.25, assay = "SCT")
DimPlot(LUNG3, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.05)

### CELL CYCLE INVESTIGATIONS ### 
LUNG3 <- CellCycleScoring(LUNG3, s.features = ccgenes.s.genes, g2m.features = ccgenes.g2m.genes)
DimPlot(LUNG3, label = T, group.by = "Phase")
FeaturePlot(LUNG3, features = c("S.Score","G2M.Score"))
# NO EFFECT OF CC ON CLUSTERING

####### COMPARTMENTS
FeaturePlot(LUNG3, reduction = "umap", features = c("EPCAM", "PECAM1", "PTPRC"), order = TRUE, label = TRUE)
VlnPlot(LUNG3, c("EPCAM", "PECAM1", "PTPRC"), log=T)

################ MESENCHYMAL ###
LUNG3_Mesench <- subset(LUNG3, subset = seurat_clusters %in% c(4,5,6,7,8,16,18))
LUNG3_Mesench <- SCTransform(LUNG3_Mesench, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG3_Mesench <- RunPCA(LUNG3_Mesench, verbose = T)
ElbowPlot(LUNG3_Mesench)
LUNG3_Mesench <- FindNeighbors(LUNG3_Mesench, dims = 1:20)
LUNG3_Mesench <- FindClusters(LUNG3_Mesench, resolution = 0.8)
LUNG3_Mesench <- RunUMAP(LUNG3_Mesench, reduction = "pca", dims = 1:20, min.dist = 0.35)
DimPlot(LUNG3_Mesench, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

FeaturePlot(LUNG3_Mesench, features = c("S.Score","G2M.Score")) # no cluster 
VlnPlot(LUNG3_Mesench, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
FeaturePlot(LUNG3_Mesench, features = c("EPCAM", "PTPRC", "JCHAIN", "PECAM1"), order = T, pt.size = 1)

# CELL TYPE INVESTIGATIONS
#GENERAL MARKERS
FeaturePlot(LUNG3_Mesench, features = c("ACTA2", "TAGLN", "MYH11", "NOTCH3"), order = T, pt.size = 1)
VlnPlot(LUNG3_Mesench, c("ACTA2", "TAGLN", "MYH11", "NOTCH3"), log=T)
FeaturePlot(LUNG3_Mesench, features = c("PDGFRA", "COL1A1", "PDGFRB"), order = T, pt.size = 1)
VlnPlot(LUNG3_Mesench, c("PDGFRA", "COL1A1", "PDGFRB"), log=T)
#VSM
FeaturePlot(LUNG3_Mesench, features = c("FABP3", "MEF2C", "PLN", "SVIL", "ITIH4"), order = T, pt.size = 1)
VlnPlot(LUNG3_Mesench, c("FABP3", "MEF2C", "PLN", "SVIL", "ITIH4"), log=T)
#APOA1+ cells
FeaturePlot(LUNG3_Mesench, features = c("ADRA2A", "COL12A1", "CLU", "RGS5", "APOA1", "HAS2", "CP", "HGF"), order = T, pt.size = 0.5)
VlnPlot(LUNG3_Mesench, c("ADRA2A", "COL12A1", "CLU", "RGS5", "APOA1", "HAS2", "CP", "HGF"), log=T)
#Pericytes
FeaturePlot(LUNG3_Mesench, features = c("FAM162B", "POSTN", "HIGD1B", "GUCY1A1"), order = T, pt.size = 1)
FeaturePlot(LUNG3_Mesench, features = c("TRPC6", "CSPG4", "MCAM", "COX4I1"), order = T, pt.size = 1)
VlnPlot(LUNG3_Mesench, c("FAM162B", "POSTN", "HIGD1B"," COX4I1", "GUCY1A1", "TRPC6", "CSPG4", "MCAM"), log=T)
#ASM
FeaturePlot(LUNG3_Mesench, features = c("ACTC1", "MYLK", "DSTN", "LPP", "DES", "PRUNE2"), order = T, pt.size = 1)
VlnPlot(LUNG3_Mesench, c("ACTC1", "MYLK", "DSTN", "LPP", "DES", "PRUNE2"), log=T)
#DUCTAL MYOFIBROBLASTS
FeaturePlot(LUNG3_Mesench, features = c("SOSTDC1", "ASPN", "HHIP", "FGF18", "CDH4"), order = T, pt.size = 1)
VlnPlot(LUNG3_Mesench, c("SOSTDC1", "ASPN", "HHIP", "FGF18", "CDH4"), log=T)
###ADV vs ALV FIBROBLASTS
FeaturePlot(LUNG3_Mesench, features = c("COL14A1", "GLI1", "COL13A1", "WNT2"), order = T, pt.size = 1)
VlnPlot(LUNG3_Mesench, c("COL14A1", "GLI1", "COL13A1", "WNT2"), log=T)
#ADVENTITIAL FIB
FeaturePlot(LUNG3_Mesench, features = c("PI16", "DCN", "AQP1", "LSP1", "LTBP1"), order = T, pt.size = 1)
VlnPlot(LUNG3_Mesench, c("PI16", "DCN", "AQP1", "LSP1", "LTBP1"), log=T)
#MATN4 ADV FIB
FeaturePlot(LUNG3_Mesench, features = c("PI3", "MATN4", "FAP", "IL33"), order = T, pt.size = 1)
VlnPlot(LUNG3_Mesench, c("PI3", "MATN4", "FAP", "IL33"), log=T)
#COL23A1 ADV FIB
FeaturePlot(LUNG3_Mesench, features = c("COL23A1", "TMEM132C", "COL15A1", "ENTPD1", "PLCL1"), order = T, pt.size = 1)
VlnPlot(LUNG3_Mesench, c("COL23A1", "TMEM132C", "COL15A1", "ENTPD1", "PLCL1"), log=T)
#Alv fib
FeaturePlot(LUNG3_Mesench, features = c("NPNT", "MACF1", "DST", "SPECC1L", "PRG4"), order = T, pt.size = 1)
VlnPlot(LUNG3_Mesench, c("NPNT", "MACF1", "DST", "SPECC1L", "PRG4"), log=T)
#IMMUNE-RECRUITING ?
FeaturePlot(LUNG3_Mesench, features = c("SAA1", "CCL19", "COL14A1", "CXCL12"), order = T, pt.size = 1)
VlnPlot(LUNG3_Mesench, c("SAA1", "CCL19", "COL14A1", "CXCL12"), log=T)
#SCHWANN CELLS - no cluster
FeaturePlot(LUNG3_Mesench, features = c("SCN7A", "SNCA", "CDH19", "MPZ"), order = T, pt.size = 1)
VlnPlot(LUNG3_Mesench, c("SCN7A", "SNCA", "CDH19", "MPZ"), log=T)
#LIPOFIBROBLASTS - No lipofibroblasts cluster found
FeaturePlot(LUNG3_Mesench, features = c("LPL", "PLIN2", "APOE"), order = T, pt.size = 1)
VlnPlot(LUNG3_Mesench, c("LPL", "PLIN2", "APOE"), log=T)
#Mesothelial cells ? no distinct cluster
FeaturePlot(LUNG3_Mesench, features = c("MSLN",	"UPK3B",	"WT1", "LRRN4", "CALB2", "CXCL13", "NKAIN4", "CDH2", "GPM6A",  "RSPO1"), order = T, pt.size = 1)
VlnPlot(LUNG3_Mesench, c("MSLN",	"UPK3B",	"WT1", "LRRN4", "CALB2", "CXCL13", "NKAIN4", "CDH2", "GPM6A",  "RSPO1"), log=T)

#####MARKERS OF NEW CLUSTERS 
Markers.SubsetMesench <- FindAllMarkers(LUNG3_Mesench, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
write.table(Markers.SubsetMesench, file = "Outs/INDIV_PROCESSING/Markers.LUNG3Mesench.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

#### Rename all identities
LUNG3_Mesench_Renamed <- RenameIdents(object = LUNG3_Mesench, 
                                      "0" = "Alveolar fibroblasts",
                                      "1" = "Pericytes",
                                      "2" = "VSM",
                                      "3" = "Adventitial fibroblasts",
                                      "4" = "Alveolar fibroblasts",
                                      "5" = "STMN2+ alveolar fibroblasts",
                                      "6" = "Adventitial fibroblasts",
                                      "7" = "Alveolar fibroblasts",
                                      "8" = "APOA1+ systemic pericytes",
                                      "9" = "Doublet",
                                      "10" = "Alveolar fibroblasts",
                                      "11" = "APOA1+ systemic pericytes",
                                      "12" = "ASM+Myofib")
LUNG3_Mesench_Renamed <- subset(LUNG3_Mesench_Renamed, idents = "Doublet", invert = TRUE)
DimPlot(LUNG3_Mesench_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
saveRDS(LUNG3_Mesench_Renamed, "Outs/INDIV_PROCESSING/LUNG3_Mesench.RDS")

############### EPITHELIAL ###
LUNG3_Epith <- subset(LUNG3, subset = seurat_clusters %in% c(14,19))
LUNG3_Epith <- SCTransform(LUNG3_Epith, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG3_Epith <- RunPCA(LUNG3_Epith, verbose = T)
ElbowPlot(LUNG3_Epith)
LUNG3_Epith <- FindNeighbors(LUNG3_Epith, dims = 1:10)
LUNG3_Epith <- FindClusters(LUNG3_Epith, resolution = 0.8)
LUNG3_Epith <- RunUMAP(LUNG3_Epith, reduction = "pca", dims = 1:10, min.dist = 0.25)
DimPlot(LUNG3_Epith, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

VlnPlot(LUNG3_Epith, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
FeaturePlot(LUNG3_Epith, features = c("nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 1)
FeaturePlot(LUNG3_Epith, features = c("EPCAM", "PTPRC", "PECAM1"), order = T, pt.size = 1)
VlnPlot(LUNG3_Epith, features = c("EPCAM", "PTPRC", "PECAM1"), ncol = 3, pt.size = F, log = T)
# cluster 4 = low quality cluster

# CELL TYPE INVESTIGATIONS
FeaturePlot(LUNG3_Epith, features = c("SFTPC", "AGER", "SCGB1A1", "GPX2", "MUC5B", "HYDIN", "KRT14", "TP63", "CFTR"), order = T, label=T, pt.size = 1)
VlnPlot(LUNG3_Epith, features = c("SFTPC", "AGER", "SCGB1A1", "GPX2", "MUC5B", "HYDIN", "KRT14", "TP63", "CFTR"), log=T)

#####MARKERS OF NEW CLUSTERS 
Markers.Epith <- FindAllMarkers(LUNG3_Epith, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
write.table(Markers.Epith, file = "Outs/INDIV_PROCESSING/Markers.LUNG3Epith.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

#### Rename all identities
LUNG3_Epith_Renamed <- RenameIdents(object = LUNG3_Epith, 
                                           "0" = "Secretory",
                                           "1" = "AT1+AT2",
                                           "2" = "Basal",
                                           "3" = "Ciliated",
                                           "4" = "LQ")
LUNG3_Epith_Renamed <- subset(LUNG3_Epith_Renamed, idents = "LQ", invert = TRUE)
DimPlot(LUNG3_Epith_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)
saveRDS(LUNG3_Epith_Renamed, "Outs/INDIV_PROCESSING/LUNG3_Epith.RDS")

################### IMMUNE ###
LUNG3_Immune <- subset(LUNG3, subset = seurat_clusters %in% c(2,9,10,11,13,15))
LUNG3_Immune <- SCTransform(LUNG3_Immune, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG3_Immune <- RunPCA(LUNG3_Immune, verbose = T)
ElbowPlot(LUNG3_Immune)
LUNG3_Immune <- FindNeighbors(LUNG3_Immune, dims = 1:7)
LUNG3_Immune <- FindClusters(LUNG3_Immune, resolution = 1)
LUNG3_Immune <- RunUMAP(LUNG3_Immune, reduction = "pca", dims = 1:7, min.dist = 0.25)
DimPlot(LUNG3_Immune, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

FeaturePlot(LUNG3_Immune, features = c("S.Score","G2M.Score")) # no cluster of replicating cells
VlnPlot(LUNG3_Immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
FeaturePlot(LUNG3_Immune, features = c("nFeature_RNA", "percent.mt"), ncol = 3)
FeaturePlot(LUNG3_Immune, features = c("EPCAM", "PTPRC", "DCN", "PECAM1"), order = T, pt.size = 1) 
# cluster 9 are doublets, cl 10 have high MT

# CELL TYPE INVESTIGATIONS
####T LYMPHOCYTES
FeaturePlot(LUNG3_Immune, features = c("CD3D", "CD3E", "IL7R", "ICOS", "THY1"), order = T, pt.size = 0.5)
VlnPlot(LUNG3_Immune, features = c("CD3D", "CD3E", "IL7R", "ICOS", "THY1"), log=T)
####CTL
FeaturePlot(LUNG3_Immune, features = c("CD3D", "CD8A", "NCR3", "GZMB", "KLRB1", "KLRK1", "CCL5"), order = T, pt.size = 0.5)
VlnPlot(LUNG3_Immune, features = c("CD3D", "CD8A", "NCR3", "GZMB", "KLRB1", "KLRK1", "CCL5"), log=T)
####NAIVE LT
FeaturePlot(LUNG3_Immune, features = c("CD3D", "LEF1", "TCF7", "IL7R", "CCR7", "SELL", "LTB"), order = T, pt.size = 0.5)
VlnPlot(LUNG3_Immune, features = c("CD3D", "LEF1", "TCF7", "IL7R", "CCR7", "SELL", "LTB"), log=T)
####LB
FeaturePlot(LUNG3_Immune, features = c("MS4A1", "CD19", "FCRLA", "CCR7"), order = T, pt.size = 0.5)
VlnPlot(LUNG3_Immune, features = c("MS4A1", "CD19", "FCRLA", "CCR7"), log=T)
####PLASMA CELLS
FeaturePlot(LUNG3_Immune, features = c("TNFRSF17", "JCHAIN", "POU2AF1"), order = T, pt.size = 0.5)
VlnPlot(LUNG3_Immune, features = c("TNFRSF17", "JCHAIN", "POU2AF1"), log=T)
####MEMORY T
FeaturePlot(LUNG3_Immune, features = c("CD3D", "IL7R", "CD7", "SAMHD1"), order = T, pt.size = 0.5)
VlnPlot(LUNG3_Immune, features = c("CD3D", "IL7R", "CD7", "SAMHD1"), log=T)
####T REGULATORY
FeaturePlot(LUNG3_Immune, features = c("CD3D", "CTLA4", "IKZF2", "GATA3", "CCL5", "LGALS3"), order = T, pt.size = 0.5)
VlnPlot(LUNG3_Immune, features = c("CD3D", "CTLA4", "IKZF2", "GATA3", "CCL5", "LGALS3"), log=T)
#T HELPER 17/1/2
FeaturePlot(LUNG3_Immune, features = c("CD3D", "IL17RB", "RORA", "TBX21", "ALDOC", "LGALS3"), order = T, pt.size = 0.5)
VlnPlot(LUNG3_Immune, features = c("CD3D", "IL17RB", "RORA", "TBX21", "ALDOC", "LGALS3"), log=T)
##MAST CELLS
FeaturePlot(LUNG3_Immune, features = c("KIT", "CPA3", "MS4A2", "FCER1A"), order = T, pt.size = 0.5)
VlnPlot(LUNG3_Immune, features = c("KIT", "CPA3", "MS4A2", "FCER1A"), log=T)
###NEUTRO
FeaturePlot(LUNG3_Immune, features = c("S100A12", "ITGAM", "CXCL8", "ADGRE1", "CD4", "SELL", "PADI3"), order = T, pt.size = 0.5)
VlnPlot(LUNG3_Immune, features = c("S100A12", "ITGAM", "CXCL8", "ADGRE1", "CD4", "SELL", "PADI3"), log=T)
##GENERAL MP
FeaturePlot(LUNG3_Immune, features = c("MRC1", "MSR1", "APOC1", "CD86", "CD48", "CD63"), order = T, pt.size = 0.5)
VlnPlot(LUNG3_Immune, features = c("MRC1", "MSR1", "APOC1", "CD86", "CD48", "CD63"), log=T)
#AM
FeaturePlot(LUNG3_Immune, features = c("MARCO", "MSR1","MRC1", "CD163", "CHI3L1", "APOC1"), order = T, pt.size = 0.5)
VlnPlot(LUNG3_Immune, features = c("MARCO", "MSR1","MRC1", "CD163", "CHI3L1", "APOC1"), log=T)
##C1QC MP
FeaturePlot(LUNG3_Immune, features = c("C1QC", "C1QA", "CCL13", "STAB1", "NPAS1", "MRC1"), order = T, pt.size = 0.5)
VlnPlot(LUNG3_Immune, features = c("C1QC", "C1QA", "CCL13", "STAB1", "NPAS1", "MRC1"), log=T)
#MONOCYTES
FeaturePlot(LUNG3_Immune, features = c("LYZ", "PLEK", "VCAN", "FN1", "IL1B", "MRC1"), order = T, pt.size = 0.5)
VlnPlot(LUNG3_Immune, features = c("LYZ", "PLEK", "VCAN", "FN1", "IL1B", "MRC1"), log=T)
#DC
FeaturePlot(LUNG3_Immune, features = c("CD1C", "PLD4", "FCER1A",  "PPM1J", "IRF8", "DLA-DOA", "PKIB"), order = T, pt.size = 0.5)
VlnPlot(LUNG3_Immune, features = c("CD1C", "PLD4", "FCER1A",  "PPM1J", "IRF8", "DLA-DOA", "PKIB"), log=T)
#MATURE DC
FeaturePlot(LUNG3_Immune, features = c("CD1C","CD83", "IDO1", "IL4I1", "CCR7"), order = T, pt.size = 0.5)
VlnPlot(LUNG3_Immune, features = c("CD1C","CD83", "IDO1", "IL4I1", "CCR7"), log=T)
#REPLICATING
FeaturePlot(LUNG3_Immune, features = c("CENPF", "MKI67", "TOP2A"), order = T, pt.size = 0.5)
VlnPlot(LUNG3_Immune, features = c("CENPF", "MKI67", "TOP2A"), log=T)

###MARKERS OF NEW CLUSTERS 
Markers.Immune <- FindAllMarkers(LUNG3_Immune, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
write.table(Markers.Immune, file = "Outs/INDIV_PROCESSING/Markers.LUNG3Immune.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

## Rename all identities
LUNG3_Immune_Renamed <- RenameIdents(object = LUNG3_Immune, 
                                    "0" = "DC", "1" = "Neutrophils", "2" = "CD4 T",
                                    "3" = "CD4 T", "4" = "Monocytes", "5" = "Neutrophils",
                                    "6" = "Cytotoxic T", "7" = "IM", "8" = "B+Plasma cells",
                                    "9" = "Doublet", "10" = "CD4 T high MT", "11" = "AM")
DimPlot(LUNG3_Immune_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)
# Remove aberrant cells
LUNG3_Immune_Renamed <- subset(LUNG3_Immune_Renamed, idents = "Doublet", invert = TRUE)
DimPlot(LUNG3_Immune_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
saveRDS(LUNG3_Immune_Renamed, "Outs/INDIV_PROCESSING/LUNG3_Immune.RDS")

###### ENDOTHELIAL ###
LUNG3_Endo <- subset(LUNG3, subset = seurat_clusters %in% c(0,1,3,12,17))
LUNG3_Endo <- SCTransform(LUNG3_Endo, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG3_Endo <- RunPCA(LUNG3_Endo, verbose = T)
ElbowPlot(LUNG3_Endo)
LUNG3_Endo <- FindNeighbors(LUNG3_Endo, dims = 1:10)
LUNG3_Endo <- FindClusters(LUNG3_Endo, resolution = 0.7)
LUNG3_Endo <- RunUMAP(LUNG3_Endo, reduction = "pca", dims = 1:10, min.dist = 0.25)
DimPlot(LUNG3_Endo, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

# QC
VlnPlot(LUNG3_Endo, features = c("percent.RPS", "percent.RPL", "percent.mt"), ncol = 3, pt.size = F, log = T)
VlnPlot(LUNG3_Endo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3, pt.size = F, log = T)
FeaturePlot(LUNG3_Endo, features = c("nFeature_RNA", "percent.mt"), ncol = 3)
FeaturePlot(LUNG3_Endo, features = c("S.Score","G2M.Score")) # no cluster of replicating cells
FeaturePlot(LUNG3_Endo, features = c("EPCAM", "PTPRC", "DCN", "PECAM1", "ACTA2"), order = T, pt.size = 1)
#9=mesenchymal doublets, 11=immune, 10=schwann cells

# CELL TYPE INVESTIGATIONS
FeaturePlot(LUNG3_Endo, features = c("VWF", "CD34", "BMX", "PROX1", "PDPN", "MGP", "SELP"), label=T, order = T, pt.size = 1)
VlnPlot(LUNG3_Endo, features = c("PECAM1","LYVE1","CDH5","VWF", "CD34", "BMX", "PROX1", "PDPN", "MGP", "SELP"), log=T)

#####MARKERS OF NEW CLUSTERS 
Markers.Endo <- FindAllMarkers(LUNG3_Endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
write.table(Markers.Endo, file = "Outs/Markers.LUNG3Endo.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

#### Rename all identities
LUNG3_Endo_Renamed <- RenameIdents(object = LUNG3_Endo, 
                                          "0" = "General capillary EC",
                                          "1" = "Aerocytes",
                                          "2" = "General capillary EC",
                                          "3" = "General capillary EC",
                                          "4" = "Lymphatic EC",
                                          "5" = "General capillary EC",
                                          "6" = "Artery EC",
                                          "7" = "Vein EC",
                                          "8" = "General capillary EC",
                                          "9" = "Doublet",
                                          "10" = "Schwann cells",
                                          "11" = "Doublet")
DimPlot(LUNG3_Endo_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)
# Remove aberrant cells
LUNG3_Endo_Renamed <- subset(LUNG3_Endo_Renamed, idents = c("Doublet"), invert = TRUE)
DimPlot(LUNG3_Endo_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
saveRDS(LUNG3_Endo_Renamed, "Outs/INDIV_PROCESSING/LUNG3_Endo.RDS")

###FUSION OF COMPARTMENTS
LUNG3_integrated <- merge(LUNG3_Endo_Renamed, y = c(LUNG3_Mesench_Renamed, LUNG3_Epith_Renamed, LUNG3_Immune_Renamed), 
                            add.cell.ids = c("Endothelial", "Mesenchymal", "Epithelial", "Immune"), project = "LUNG3")
LUNG3_integrated <- SCTransform(LUNG3_integrated, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG3_integrated <- RunPCA(LUNG3_integrated, npcs = 50, verbose = T, assay = "SCT")
LUNG3_integrated <- RunUMAP(LUNG3_integrated, reduction = "pca", dims = 1:50, min.dist = 0.25, assay = "SCT")
DimPlot(LUNG3_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.05)+NoLegend()
saveRDS(LUNG3_integrated, "Outs/INDIV_PROCESSING/LUNG3_integrated.RDS")
# end of sample 3

# PROCESSING SAMPLE 4 ####
### SETUP SEURAT OBJECTS ###
LUNG4.data <- Read10X(data.dir = "Raw_data/Lung4/filtered_feature_bc_matrix/")
LUNG4 <- CreateSeuratObject(counts = LUNG4.data, project = "LUNG4", min.cells = 10, min.features = 200)

### PRE-PROCESSING WORKFLOW ###
LUNG4 <- PercentageFeatureSet(LUNG4, pattern = "^MT-", col.name = "percent.mt")
LUNG4 <- PercentageFeatureSet(LUNG4, pattern = "^RPS", col.name = "percent.RPS")
LUNG4 <- PercentageFeatureSet(LUNG4, pattern = "^RPL", col.name = "percent.RPL")

VlnPlot(LUNG4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0.02, log = F, group.by = "orig.ident")
VlnPlot(LUNG4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.02, log = T, group.by = "orig.ident")
VlnPlot(LUNG4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = F, log = T, group.by = "orig.ident")
VlnPlot(LUNG4, features = c("percent.RPS", "percent.RPL", "percent.mt"),
        ncol = 3, pt.size = F, log = F, group.by = "orig.ident")
LUNG4 <- subset(LUNG4, subset = percent.mt < 20)
VlnPlot(LUNG4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =F, log = F, group.by = "orig.ident")

### PROCESSING
LUNG4 <- SCTransform(LUNG4, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG4 <- RunPCA(LUNG4, npcs = 50, verbose = T, assay = "SCT")
ElbowPlot(LUNG4)
LUNG4 <- FindNeighbors(LUNG4, dims = 1:50, assay = "SCT")
LUNG4 <- FindClusters(LUNG4, resolution = 0.8)
LUNG4 <- RunUMAP(LUNG4, reduction = "pca", dims = 1:50, min.dist = 0.25, assay = "SCT")
DimPlot(LUNG4, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.05)

### CELL CYCLE INVESTIGATIONS ###
LUNG4 <- CellCycleScoring(LUNG4, s.features = ccgenes.s.genes, g2m.features = ccgenes.g2m.genes)
DimPlot(LUNG4, label = T, group.by = "Phase")
FeaturePlot(LUNG4, features = c("S.Score","G2M.Score"))

####### COMPARTMENTS
FeaturePlot(LUNG4, reduction = "umap", 
            features = c("EPCAM", "PECAM1", "PTPRC"), order = TRUE, label = TRUE)
VlnPlot(LUNG4, c("EPCAM", "PECAM1", "PTPRC"), log=T)

############# MESENCHYMAL ###
LUNG4_Mesench <- subset(LUNG4, subset = seurat_clusters %in% c(6,12,2,5,7,9))
LUNG4_Mesench <- SCTransform(LUNG4_Mesench, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG4_Mesench <- RunPCA(LUNG4_Mesench, verbose = T)
LUNG4_Mesench <- FindNeighbors(LUNG4_Mesench, dims = 1:50)
LUNG4_Mesench <- FindClusters(LUNG4_Mesench, resolution = 1)
LUNG4_Mesench <- RunUMAP(LUNG4_Mesench, reduction = "pca", dims = 1:50, min.dist = 0.25)
DimPlot(LUNG4_Mesench, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

LUNG4_Mesench <- CellCycleScoring(LUNG4_Mesench, s.features = ccgenes.s.genes, g2m.features = ccgenes.g2m.genes)
FeaturePlot(LUNG4_Mesench, features = c("S.Score","G2M.Score"))
VlnPlot(LUNG4_Mesench, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
# CLUSTER 7 IS LOW QUALITY

# DOUBLETS ?
FeaturePlot(LUNG4_Mesench, features = c("EPCAM", "PTPRC", "PECAM1"), order = T, pt.size = 1) # CLUSTER 7 = immune doublets
FeaturePlot(LUNG4_Mesench, features = c("percent.mt"), order = T, pt.size = 1)

# 1st CELL TYPE INVESTIGATION
#GENERAL
FeaturePlot(LUNG4_Mesench, features = c("ACTA2", "TAGLN", "MYH11", "NOTCH3"), order = T, pt.size = 1)
VlnPlot(LUNG4_Mesench, c("ACTA2", "TAGLN", "MYH11", "NOTCH3"), log=T)
FeaturePlot(LUNG4_Mesench, features = c("PDGFRA", "COL1A1", "PDGFRB"), order = T, pt.size = 1)
VlnPlot(LUNG4_Mesench, c("PDGFRA", "COL1A1", "PDGFRB"), log=T)
#VSM
FeaturePlot(LUNG4_Mesench, features = c("FABP3", "MEF2C", "PLN", "SVIL", "ITIH4"), order = T, pt.size = 1)
VlnPlot(LUNG4_Mesench, c("FABP3", "MEF2C", "PLN", "SVIL", "ITIH4"), log=T)
#APOA1 cells
FeaturePlot(LUNG4_Mesench, features = c("ADRA2A", "COL12A1", "CLU", "RGS5", "APOA1", "HAS2", "CP", "HGF"), order = T, pt.size = 0.5)
VlnPlot(LUNG4_Mesench, c("ADRA2A", "COL12A1", "CLU", "RGS5", "APOA1", "HAS2", "CP", "HGF"), log=T)
#Pericytes
FeaturePlot(LUNG4_Mesench, features = c("FAM162B", "POSTN", "HIGD1B", "GUCY1A1"), order = T, pt.size = 1)
VlnPlot(LUNG4_Mesench, c("FAM162B", "POSTN", "HIGD1B"," COX4I1", "GUCY1A1", "TRPC6", "CSPG4", "MCAM"), log=T)
#ASM
FeaturePlot(LUNG4_Mesench, features = c("ACTC1", "MYLK", "DSTN", "LPP", "DES", "PRUNE2"), order = T, pt.size = 1)
VlnPlot(LUNG4_Mesench, c("ACTC1", "MYLK", "DSTN", "LPP", "DES", "PRUNE2"), log=T)
#DUCTAL MYOFIBROBLASTS
FeaturePlot(LUNG4_Mesench, features = c("SOSTDC1", "ASPN", "HHIP", "FGF18", "CDH4"), order = T, pt.size = 1)
VlnPlot(LUNG4_Mesench, c("SOSTDC1", "ASPN", "HHIP", "FGF18", "CDH4"), log=T)
###ADV vs ALV FIBROBLASTS
FeaturePlot(LUNG4_Mesench, features = c("COL14A1", "GLI1", "COL13A1", "WNT2"), order = T, pt.size = 1)
VlnPlot(LUNG4_Mesench, c("COL14A1", "GLI1", "COL13A1", "WNT2"), log=T)
#Alv fib
FeaturePlot(LUNG4_Mesench, features = c("NPNT", "MACF1", "DST", "SPECC1L", "PRG4"), order = T, pt.size = 1)
VlnPlot(LUNG4_Mesench, c("NPNT", "MACF1", "DST", "SPECC1L", "PRG4"), log=T)
#SCHWANN CELLS
FeaturePlot(LUNG4_Mesench, features = c("SCN7A", "SNCA", "CDH19", "MPZ"), order = T, pt.size = 1)
VlnPlot(LUNG4_Mesench, c("SCN7A", "SNCA", "CDH19", "MPZ"), log=T)
#LIPOFIBROBLASTS - no defined cluster
FeaturePlot(LUNG4_Mesench, features = c("LPL", "PLIN2", "APOE"), order = T, pt.size = 1)
VlnPlot(LUNG4_Mesench, c("LPL", "PLIN2", "APOE"), log=T)
#Mesothelial - no cluster
FeaturePlot(LUNG4_Mesench, features = c("MSLN",	"UPK3B",	"WT1", "LRRN4", "CALB2", "CXCL13", "NKAIN4", "CDH2", "GPM6A",  "RSPO1"), order = T, pt.size = 1)
VlnPlot(LUNG4_Mesench, c("MSLN",	"UPK3B",	"WT1", "LRRN4", "CALB2", "CXCL13", "NKAIN4", "CDH2", "GPM6A",  "RSPO1"), log=T)

#####MARKERS OF NEW CLUSTERS 
Markers.SubsetMesench <- FindAllMarkers(LUNG4_Mesench, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
write.table(Markers.SubsetMesench, file = "Outs/INDIV_PROCESSING/Markers.LUNG4Mesench.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

#### Rename all identities
LUNG4_Mesench_Renamed <- RenameIdents(object = LUNG4_Mesench, 
                                              "0" = "Alveolar fibroblasts",
                                              "1" = "Alveolar fibroblasts",
                                              "2" = "Adventitial fibroblasts",
                                              "3" = "Pericytes",
                                              "4" = "Alveolar fibroblasts", 
                                              "5" = "Alveolar fibroblasts",
                                              "6" = "VSM+ASM+Myofib",
                                              "7" = "Doublet",
                                              "8" = "APOA1+ systemic pericytes")
DimPlot(LUNG4_Mesench_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
# Remove aberrant cells
LUNG4_Mesench_Renamed <- subset(LUNG4_Mesench_Renamed, idents = "Doublet", invert = TRUE)
DimPlot(LUNG4_Mesench_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
saveRDS(LUNG4_Mesench_Renamed, "Outs/INDIV_PROCESSING/LUNG4_Mesench.RDS")

#LUNG4_Mesench_Renamed <- readRDS("Outs/INDIV_PROCESSING/LUNG4_Mesench.RDS")

################# EPITHELIAL ###
LUNG4_Epith <- subset(LUNG4, subset = seurat_clusters %in% c(3,10,13))
LUNG4_Epith <- SCTransform(LUNG4_Epith, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG4_Epith <- RunPCA(LUNG4_Epith, npcs = 50, verbose = T)
ElbowPlot(LUNG4_Epith)
LUNG4_Epith <- FindNeighbors(LUNG4_Epith, dims = 1:5)
LUNG4_Epith <- FindClusters(LUNG4_Epith, resolution = 0.6)
LUNG4_Epith <- RunUMAP(LUNG4_Epith, reduction = "pca", dims = 1:5, min.dist = 0.15)
DimPlot(LUNG4_Epith, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

LUNG4_Epith <- CellCycleScoring(LUNG4_Epith, s.features = ccgenes.s.genes, g2m.features = ccgenes.g2m.genes)
FeaturePlot(LUNG4_Epith, features = c("S.Score","G2M.Score")) # NO CLUSTER IS EFFECT OF CELL CYCLE
VlnPlot(LUNG4_Epith, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
FeaturePlot(LUNG4_Epith, features = c("EPCAM", "PTPRC", "PECAM1", "JCHAIN"), order = T, label=T,pt.size = 1)
#CL5= DOUBLET ENDOTHELIAL / CL6= DOUBLET IMMUNE

# CELL TYPE INVESTIGATION
#CLUB CELLS
FeaturePlot(LUNG4_Epith, features = c("SCGB1A1", "AQP5", "SERPINA1", "ERN2", "DEPP1", "GPX2"), order = T, label=T, pt.size = 1)
VlnPlot(LUNG4_Epith, features = c("SCGB1A1", "AQP5", "SERPINA1", "ERN2", "DEPP1", "GPX2"), log=T)
#Ionocytes, Goblet
FeaturePlot(LUNG4_Epith, features = c("CFTR", "ASCL3", "MUC5B", "SPDEF"), order = T, label=T, pt.size = 1)
VlnPlot(LUNG4_Epith, features = c("CFTR", "ASCL3", "MUC5B", "SPDEF"), log=T)
#AT1
FeaturePlot(LUNG4_Epith, features = c("AGER", "AQP1", "SUSD2", "CLDN18"), order = T, label=T,pt.size = 1)
VlnPlot(LUNG4_Epith, features = c("AGER", "AQP1", "SUSD2", "CLDN18"), log=T)
#AT2
FeaturePlot(LUNG4_Epith, features = c("NAPSA", "SFTPC", "SFTPB", "MUC1"), order = T, label=T,pt.size = 1)
VlnPlot(LUNG4_Epith, features = c("NAPSA", "SFTPC", "SFTPB", "MUC1"), log=T)
#BASAL
FeaturePlot(LUNG4_Epith, features = c("KRT14", "TP63"), order = T, label=T,pt.size = 1)
VlnPlot(LUNG4_Epith, features = c("KRT14", "TP63"), log=T)
#CILIATED
FeaturePlot(LUNG4_Epith, features = c("CAPS", "FOXJ1", "CCDC78", "HYDIN"), order = T, label=T,pt.size = 1)
VlnPlot(LUNG4_Epith, features = c("CAPS", "FOXJ1", "CCDC78", "HYDIN"), log=T)
#REPLICATING
FeaturePlot(LUNG4_Epith, features = c("MKI67","CENPF"), order = T, label=T,pt.size = 1)
VlnPlot(LUNG4_Epith, c("MKI67","CENPF"), log=T)

#####MARKERS OF NEW CLUSTERS 
Markers.Epith <- FindAllMarkers(LUNG4_Epith, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
write.table(Markers.Epith, file = "Outs/INDIV_PROCESSING/Markers.LUNG4_Epith.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

#### Rename all identities
LUNG4_Epith_Renamed <- RenameIdents(object = LUNG4_Epith,
                                             "0" = "AT2",
                                             "1" = "AT2",
                                             "2" = "AT2",
                                             "3" = "AT2",
                                             "4" = "AT1",
                                             "5" = "Doublet",
                                             "6" = "Doublet",
                                             "7" = "Secretory+Basal+Ciliated")
DimPlot(LUNG4_Epith_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)
# Remove aberrant cells
LUNG4_Epith_Renamed <- subset(LUNG4_Epith_Renamed, idents = "Doublet", invert = TRUE)
DimPlot(LUNG4_Epith_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
saveRDS(LUNG4_Epith_Renamed, "Outs/INDIV_PROCESSING/LUNG4_Epith.RDS")

###### IMMUNE ###
LUNG4_Immune <- subset(LUNG4, subset = seurat_clusters %in% c(0,4,8,15,11,14,16))
LUNG4_Immune <- SCTransform(LUNG4_Immune, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG4_Immune <- RunPCA(LUNG4_Immune, npcs = 20, verbose = T)
ElbowPlot(LUNG4_Immune)
LUNG4_Immune <- FindNeighbors(LUNG4_Immune, dims = 1:7)
LUNG4_Immune <- FindClusters(LUNG4_Immune, resolution = 0.6)
LUNG4_Immune <- RunUMAP(LUNG4_Immune, reduction = "pca", dims = 1:7, min.dist = 0.25)
DimPlot(LUNG4_Immune, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

LUNG4_Immune <- CellCycleScoring(LUNG4_Immune, s.features = ccgenes.s.genes, g2m.features = ccgenes.g2m.genes)
FeaturePlot(LUNG4_Immune, features = c("S.Score","G2M.Score"))
VlnPlot(LUNG4_Immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
FeaturePlot(LUNG4_Immune, features = c("EPCAM", "COL1A1", "PECAM1"), order = T, pt.size = 1)#no cluster of doublets

# CELL TYPE INVESTIGATION
####LYMPHOCYTES
FeaturePlot(LUNG4_Immune, features = c("CD3D", "CD3E", "IL7R", "ICOS", "THY1"), order = T, pt.size = 0.5)
VlnPlot(LUNG4_Immune, features = c("CD3D", "CD3E", "IL7R", "ICOS", "THY1"), log=T)
####CTL
FeaturePlot(LUNG4_Immune, features = c("CD3D", "CD8A", "NCR3", "GZMB", "KLRB1", "KLRK1", "CCL5"), order = T, pt.size = 0.5)
VlnPlot(LUNG4_Immune, features = c("CD3D", "CD8A", "NCR3", "GZMB", "KLRB1", "KLRK1", "CCL5"), log=T)
####NAIVE LT
FeaturePlot(LUNG4_Immune, features = c("CD3D", "LEF1", "TCF7", "IL7R", "CCR7", "SELL", "LTB"), order = T, pt.size = 0.5)
VlnPlot(LUNG4_Immune, features = c("CD3D", "LEF1", "TCF7", "IL7R", "CCR7", "SELL", "LTB"), log=T)
####LB
FeaturePlot(LUNG4_Immune, features = c("MS4A1", "CD19", "FCRLA", "CCR7"), order = T, pt.size = 0.5)
VlnPlot(LUNG4_Immune, features = c("MS4A1", "CD19", "FCRLA", "CCR7"), log=T)
####PLASMA CELLS
FeaturePlot(LUNG4_Immune, features = c("TNFRSF17", "JCHAIN", "POU2AF1"), order = T, pt.size = 0.5)
VlnPlot(LUNG4_Immune, features = c("TNFRSF17", "JCHAIN", "POU2AF1"), log=T)
####MEMORY T
FeaturePlot(LUNG4_Immune, features = c("CD3D", "IL7R", "CD7", "SAMHD1"), order = T, pt.size = 0.5)
VlnPlot(LUNG4_Immune, features = c("CD3D", "IL7R", "CD7", "SAMHD1"), log=T)
####T REGULATORY
FeaturePlot(LUNG4_Immune, features = c("CD3D", "CTLA4", "IKZF2", "GATA3", "CCL5", "LGALS3"), order = T, pt.size = 0.5)
VlnPlot(LUNG4_Immune, features = c("CD3D", "CTLA4", "IKZF2", "GATA3", "CCL5", "LGALS3"), log=T)
#T HELPER 17/1/2
FeaturePlot(LUNG4_Immune, features = c("CD3D", "IL17RB", "RORA", "TBX21", "ALDOC", "LGALS3"), order = T, pt.size = 0.5)
VlnPlot(LUNG4_Immune, features = c("CD3D", "IL17RB", "RORA", "TBX21", "ALDOC", "LGALS3"), log=T)
##MAST CELLS
FeaturePlot(LUNG4_Immune, features = c("KIT", "CPA3", "MS4A2", "FCER1A"), order = T, pt.size = 0.5)
VlnPlot(LUNG4_Immune, features = c("KIT", "CPA3", "MS4A2", "FCER1A"), log=T)
###NEUTRO
FeaturePlot(LUNG4_Immune, features = c("S100A12", "ITGAM", "CXCL8", "ADGRE1", "CD4", "SELL", "PADI3"), order = T, pt.size = 0.5)
VlnPlot(LUNG4_Immune, features = c("S100A12", "ITGAM", "CXCL8", "ADGRE1", "CD4", "SELL", "PADI3"), log=T)
##GENERAL MP
FeaturePlot(LUNG4_Immune, features = c("MRC1", "MSR1", "APOC1", "CD86", "CD48", "CD63"), order = T, pt.size = 0.5)
VlnPlot(LUNG4_Immune, features = c("MRC1", "MSR1", "APOC1", "CD86", "CD48", "CD63"), log=T)
#AM
FeaturePlot(LUNG4_Immune, features = c("MARCO", "MSR1","MRC1", "CD163", "CHI3L1", "APOC1"), order = T, pt.size = 0.5)
VlnPlot(LUNG4_Immune, features = c("MARCO", "MSR1","MRC1", "CD163", "CHI3L1", "APOC1"), log=T)
##C1QC MP
FeaturePlot(LUNG4_Immune, features = c("C1QC", "C1QA", "CCL13", "STAB1", "NPAS1", "MRC1"), order = T, pt.size = 0.5)
VlnPlot(LUNG4_Immune, features = c("C1QC", "C1QA", "CCL13", "STAB1", "NPAS1", "MRC1"), log=T)
#MONOCYTES
FeaturePlot(LUNG4_Immune, features = c("LYZ", "PLEK", "VCAN", "FN1", "IL1B", "MRC1"), order = T, pt.size = 0.5)
VlnPlot(LUNG4_Immune, features = c("LYZ", "PLEK", "VCAN", "FN1", "IL1B", "MRC1"), log=T)
#DC
FeaturePlot(LUNG4_Immune, features = c("CD1C", "PLD4", "FCER1A",  "PPM1J", "IRF8", "DLA-DOA", "PKIB"), order = T, pt.size = 0.5)
VlnPlot(LUNG4_Immune, features = c("CD1C", "PLD4", "FCER1A",  "PPM1J", "IRF8", "DLA-DOA", "PKIB"), log=T)
#MATURE DC
FeaturePlot(LUNG4_Immune, features = c("CD1C","CD83", "IDO1", "IL4I1", "CCR7"), order = T, pt.size = 0.5)
VlnPlot(LUNG4_Immune, features = c("CD1C","CD83", "IDO1", "IL4I1", "CCR7"), log=T)
#REPLICATING
FeaturePlot(LUNG4_Immune, features = c("CENPF", "MKI67", "TOP2A"), order = T, pt.size = 0.5)
VlnPlot(LUNG4_Immune, features = c("CENPF", "MKI67", "TOP2A"), log=T)
## Rename all identities
LUNG4_Immune_Renamed <- RenameIdents(object = LUNG4_Immune,
                                    "0" = "AM",
                                    "1" = "AM",
                                    "2" = "AM",
                                    "3" = "IM",
                                    "4" = "IM",
                                    "5" = "T Lymphocytes + NK",
                                    "6" = "DC",
                                    "7" = "Neutrophils",
                                    "8" = "Monocytes",
                                    "9" = "Plasma cells")
DimPlot(LUNG4_Immune_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
saveRDS(LUNG4_Immune_Renamed, "Outs/INDIV_PROCESSING/LUNG4_Immune.RDS")

###################### ENDOTHELIAL ###
LUNG4_Endo <- subset(LUNG4, subset = seurat_clusters %in% c(1,17))
LUNG4_Endo <- SCTransform(LUNG4_Endo, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG4_Endo <- RunPCA(LUNG4_Endo, npcs = 20, verbose = T)
ElbowPlot(LUNG4_Endo)
LUNG4_Endo <- FindNeighbors(LUNG4_Endo, dims = 1:5)
LUNG4_Endo <- FindClusters(LUNG4_Endo, resolution = 0.6)
LUNG4_Endo <- RunUMAP(LUNG4_Endo, reduction = "pca", dims = 1:5, min.dist = 0.25)
DimPlot(LUNG4_Endo, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

LUNG4_Endo <- CellCycleScoring(LUNG4_Endo, s.features = ccgenes.s.genes, g2m.features = ccgenes.g2m.genes)
FeaturePlot(LUNG4_Endo, features = c("S.Score","G2M.Score"))
VlnPlot(LUNG4_Endo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.RPS"), ncol = 3, pt.size = F, log = T)
FeaturePlot(LUNG4_Endo, features = c("EPCAM", "COL1A1", "PTPRC"), order = T, pt.size = 1)
VlnPlot(LUNG4_Endo, features = c("EPCAM", "COL1A1", "PTPRC"), ncol = 3, pt.size = F, log = T) 
# cluster 6 is doublet immune

# CELL TYPE INVESTIGATION
#ENDO GENERAL
FeaturePlot(LUNG4_Endo, features = c("PECAM1", "LYVE1", "CDH5", "VWF", "CD34"), order = T, pt.size = 1)
VlnPlot(LUNG4_Endo, features = c("PECAM1", "LYVE1", "CDH5", "VWF", "CD34"), log=T)
##LYMPHATIC
FeaturePlot(LUNG4_Endo, features = c("VWF", "PROX1", "PDPN", "TBX1", "THY1", "FLT4"), order = T, pt.size = 1)
VlnPlot(LUNG4_Endo, features = c("VWF", "PROX1", "PDPN", "TBX1", "THY1", "FLT4"), log=T)
#ARTERY EC
FeaturePlot(LUNG4_Endo, features = c("BMX", "MGP", "GJA5", "EFNB2"), order = T, pt.size = 1)
VlnPlot(LUNG4_Endo, features = c("BMX", "MGP", "GJA5", "EFNB2"), log=T)
##VEINOUS EC
FeaturePlot(LUNG4_Endo, features = c("SELP", "SELE", "VCAM1", "TIMP1", "ACKR1"), order = T, pt.size = 1)
VlnPlot(LUNG4_Endo, features = c("SELP", "SELE", "VCAM1", "TIMP1", "ACKR1"), log=T)
#CAPILLARY EC
FeaturePlot(LUNG4_Endo, features = c("CA4", "SPARC", "SGK1", "EMP2"), order = T, pt.size = 1)
VlnPlot(LUNG4_Endo, features = c("CA4", "SPARC", "SGK1", "EMP2"), log=T)
##AEROCYTE
FeaturePlot(LUNG4_Endo, features = c("VWF", "EMCN", "EDNRB", "TBX2"), order = T, pt.size = 1)
VlnPlot(LUNG4_Endo, features = c("VWF", "EMCN", "EDNRB", "TBX2"), log=T)
##GENERAL CAPILLARY
FeaturePlot(LUNG4_Endo, features = c("EDN1", "GPIHBP1", "CD36"), order = T, pt.size = 1)
VlnPlot(LUNG4_Endo, features = c("EDN1", "GPIHBP1", "CD36"), log=T)

# #####MARKERS OF NEW CLUSTERS 
# Markers.Endo <- FindAllMarkers(LUNG4_Endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
# write.table(Markers.Endo, file = "Outs/INDIV_PROCESSING/Markers.LUNG4_Endo.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

#### Rename all identities
LUNG4_Endo_Renamed <- RenameIdents(object = LUNG4_Endo, 
                                           "0" = "Endothelial cells","1" = "Endothelial cells",
                                           "2" = "Aerocytes", "3" = "Aerocytes",
                                           "4" = "Endothelial cells","5" = "Lymphatic EC","6" = "Doublet")
DimPlot(LUNG4_Endo_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)
# Remove aberrant cells
LUNG4_Endo_Renamed <- subset(LUNG4_Endo_Renamed, idents = c("Doublet"), invert = TRUE)
# Re-visualize the clusters
DimPlot(LUNG4_Endo_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1.5)
saveRDS(LUNG4_Endo_Renamed, "Outs/INDIV_PROCESSING/LUNG4_Endo.RDS")

###FUSION OF COMPARTMENTS
LUNG4_integrated <- merge(LUNG4_Endo_Renamed, y = c(LUNG4_Mesench_Renamed, LUNG4_Epith_Renamed, LUNG4_Immune_Renamed), 
                                   add.cell.ids = c("Endothelial", "Mesenchymal", "Epithelial", "Immune"), project = "LUNG4")
LUNG4_integrated <- SCTransform(LUNG4_integrated, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
LUNG4_integrated <- RunPCA(LUNG4_integrated, npcs = 20, verbose = T, assay = "SCT")
LUNG4_integrated <- RunUMAP(LUNG4_integrated, reduction = "pca", dims = 1:20, min.dist = 0.25, assay = "SCT")
DimPlot(LUNG4_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.05)+NoLegend()
saveRDS(LUNG4_integrated, "Outs/INDIV_PROCESSING/LUNG4_integrated.RDS")
#end of sample 4

# INTEGRATION OF 4 SAMPLES ####

# LOADING DATA
cells1 <- readRDS("Outs/INDIV_PROCESSING/LUNG1_integrated.RDS")
cells2 <- readRDS("Outs/INDIV_PROCESSING/LUNG2_integrated.RDS")
cells3 <- readRDS("Outs/INDIV_PROCESSING/LUNG3_integrated.RDS")
cells4 <- readRDS("Outs/INDIV_PROCESSING/LUNG4_integrated.RDS")

#MERGING ALL in one
Healthymerged <- merge(cells1, y = c(cells2, cells3, cells4), 
                         add.cell.ids = c("LUNG1", "LUNG2", "LUNG3", "LUNG4"), project = "HEALTHY")

# save previous annotation to metadata
Healthymerged$initial_annotation <- Idents(Healthymerged)

### Investigate merging before integration
# Healthymerged_SCT <- SCTransform(Healthymerged, vst.flavor = "v2", vars.to.regress = c("percent.mt"), verbose = FALSE)
# Healthymerged_SCT <- CellCycleScoring(Healthymerged_SCT, s.features = ccgenes.s.genes, g2m.features = ccgenes.g2m.genes)
# DefaultAssay(Healthymerged_SCT) <- "RNA"
# Healthymerged_SCT <- SCTransform(Healthymerged_SCT, vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)
# Healthymerged_SCT <- RunPCA(Healthymerged_SCT, npcs = 50, verbose = T, assay = "SCT")
# ElbowPlot(Healthymerged_SCT)
# Healthymerged_SCT <- FindNeighbors(Healthymerged_SCT, dims = 1:20, assay = "SCT")
# Healthymerged_SCT <- FindClusters(Healthymerged_SCT, resolution = 0.1)
# Healthymerged_SCT <- RunUMAP(Healthymerged_SCT, reduction = "pca", dims = 1:20, min.dist = 0.25, assay = "SCT")
# DimPlot(Healthymerged_SCT, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.05)
# DimPlot(Healthymerged_SCT, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.05, group.by = "orig.ident")
# Conclusion:
# Merge dataset shows bias by dog, batch effect, integration is needed

# Split object for integration
Healthysplit <- SplitObject(Healthymerged, split.by = "orig.ident")
#NORMALIZATION
for (i in 1:length(Healthysplit)) {
  Healthysplit[[i]] <- SCTransform(Healthysplit[[i]], vst.flavor = "v2", vars.to.regress = c("percent.mt"), verbose = FALSE)
}
for (i in 1:length(Healthysplit)) {
  Healthysplit[[i]] <- CellCycleScoring(Healthysplit[[i]], s.features = ccgenes.s.genes, g2m.features = ccgenes.g2m.genes)
}
for (i in 1:length(Healthysplit)) {
  Healthysplit[[i]] <- SCTransform(Healthysplit[[i]], vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)
}

#INTEGRATION
integ_features <- SelectIntegrationFeatures(object.list = Healthysplit, nfeatures = 3000) 
Healthysplit <- PrepSCTIntegration(object.list = Healthysplit, anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = Healthysplit, normalization.method = "SCT", anchor.features = integ_features)
seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")

####50 NPCS
DefaultAssay(seurat_integrated) <- "integrated"
seurat_integrated <- RunPCA(seurat_integrated, npcs = 100, verbose = FALSE)
ElbowPlot(seurat_integrated)
seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "pca", dims = 1:100)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 2.8)
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:100, min.dist = 0.3)
saveRDS(seurat_integrated, "Outs/Unannotated/Integration_Healthy.RDS")

seurat_integrated <- readRDS("Outs/Unannotated/Integration_Healthy.RDS")

pdf("Outs/GLOBAL_PROCESSING/Integration_clusters.pdf")
DimPlot(seurat_integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
DimPlot(seurat_integrated, reduction = "umap", group.by = "orig.ident")
DimPlot(seurat_integrated, reduction = "umap", split.by = "orig.ident")
DimPlot(seurat_integrated, reduction = "umap", group.by = "initial_annotation", label=T)+NoLegend()
dev.off()

### CELL CYCLE INVESTIGATIONS ###
DimPlot(seurat_integrated, label = T, group.by = "Phase")
FeaturePlot(seurat_integrated, features = c("S.Score","G2M.Score"))

########## IDENTIFYING CELL COMPARTMENTS ###
DefaultAssay(seurat_integrated) <- "SCT"
FeaturePlot(seurat_integrated, reduction = "umap", features = c("EPCAM", "PECAM1", "PTPRC"), order = TRUE, label = TRUE)
FeaturePlot(seurat_integrated, reduction = "umap", features = c("EPCAM"), order = TRUE, label = TRUE)
FeaturePlot(seurat_integrated, reduction = "umap", features = c("PECAM1"), order = TRUE, label = TRUE)
FeaturePlot(seurat_integrated, reduction = "umap", features = c("PTPRC"), order = TRUE, label = TRUE)
VlnPlot(seurat_integrated, c("EPCAM", "PECAM1", "PTPRC"), log=T)

#### MESENCHYMAL COMPARTMENT ####
DefaultAssay(seurat_integrated) <- "integrated"
Mesench_merged <- subset(seurat_integrated, subset = seurat_clusters %in% c(2,3,9,11,14,21,29,32,34,44,7,22,25,33,38,43,49))
Mesench_split <- SplitObject(Mesench_merged, split.by = "orig.ident")
for (i in 1:length(Mesench_split)) {
  Mesench_split[[i]] <- SCTransform(Mesench_split[[i]], vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)
}
integ_features <- SelectIntegrationFeatures(object.list = Mesench_split, nfeatures = 3000) 
Mesench_split <- PrepSCTIntegration(object.list = Mesench_split, anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = Mesench_split, normalization.method = "SCT", anchor.features = integ_features)
Healthy_Mesench <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
Healthy_Mesench <- RunPCA(Healthy_Mesench, npcs = 50, verbose = T)
ElbowPlot(Healthy_Mesench)
Healthy_Mesench <- FindNeighbors(Healthy_Mesench, dims = 1:15)
Healthy_Mesench <- FindClusters(Healthy_Mesench, resolution = 0.6)
Healthy_Mesench <- RunUMAP(Healthy_Mesench, reduction = "pca", dims = 1:15, min.dist = 0.35)
DimPlot(Healthy_Mesench, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

DefaultAssay(Healthy_Mesench) <- "SCT"
pdf("Outs/GLOBAL_PROCESSING/Healthy_MesenchymalClusters.pdf")
DimPlot(Healthy_Mesench, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(Healthy_Mesench, reduction = "umap", label = TRUE, repel = TRUE, split.by = "orig.ident", pt.size = 1)
DimPlot(Healthy_Mesench, reduction = "umap", label = TRUE, repel = TRUE, group.by = "orig.ident", pt.size = 1)
VlnPlot(Healthy_Mesench, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
DimPlot(Healthy_Mesench, reduction = "umap", group.by = "initial_annotation", label=T)+NoLegend()
FeaturePlot(Healthy_Mesench, reduction = "umap", features = c("EPCAM"), order = TRUE, label = TRUE)
FeaturePlot(Healthy_Mesench, reduction = "umap", features = c("PECAM1"), order = TRUE, label = TRUE)
FeaturePlot(Healthy_Mesench, reduction = "umap", features = c("PTPRC"), order = TRUE, label = TRUE)
VlnPlot(Healthy_Mesench, c("EPCAM"), log=T)
VlnPlot(Healthy_Mesench, c("PECAM1"), log=T)
VlnPlot(Healthy_Mesench, c("PTPRC"), log=T)
dev.off()
# cluster 8 are doublets with immune cells

# CELL TYPE INVESTIGATION
DefaultAssay(Healthy_Mesench) <- "SCT"
pdf("Outs/GLOBAL_PROCESSING/Healthy_MesenchymalCellTypes.pdf")
#ADVENTITIAL FIB
FeaturePlot(Healthy_Mesench, features = c("PI16", "DCN", "COL14A1", "AQP1", "LSP1", "LTBP1"), order = T, pt.size = 1)
VlnPlot(Healthy_Mesench, c("PI16", "DCN", "COL14A1", "AQP1", "LSP1", "LTBP1"), log=T)
FeaturePlot(Healthy_Mesench, features = c("PI3", "MATN4", "COL14A1", "DCN", "FAP", "IL33"), order = T, pt.size = 1)
VlnPlot(Healthy_Mesench, c("PI3", "MATN4", "COL14A1", "DCN", "FAP", "IL33"), log=T)
#Alv fib
FeaturePlot(Healthy_Mesench, features = c("NPNT", "COL13A1", "MACF1", "DST", "SPECC1L", "PRG4", "STMN2"), order = T, pt.size = 1)
VlnPlot(Healthy_Mesench, c("NPNT", "COL13A1", "MACF1", "DST", "SPECC1L", "PRG4", "STMN2"), log=T)
#PERIBRONCHIAL FIBROBLASTS
FeaturePlot(Healthy_Mesench, features = c("ASPN", "HHIP", "FGF18", "WIF1"), order = T, pt.size = 1)
VlnPlot(Healthy_Mesench, c("ASPN", "HHIP", "FGF18", "WIF1"), log=T)
#VSM
FeaturePlot(Healthy_Mesench, features = c("FABP3", "MEF2C", "PLN", "ACTA2", "TAGLN", "MYH11", "SVIL", "ITIH4"), order = T, pt.size = 1)
VlnPlot(Healthy_Mesench, c("FABP3", "MEF2C", "PLN", "ACTA2", "TAGLN", "MYH11", "SVIL", "ITIH4"), log=T)
FeaturePlot(Healthy_Mesench, features = c("ADRA2A", "COL12A1", "CLU", "RGS5", "COL1A1", "APOA1", "HAS2", "CP", "HGF"), order = T, pt.size = 0.5)
VlnPlot(Healthy_Mesench, c("ADRA2A", "COL12A1", "CLU", "RGS5", "COL1A1", "APOA1", "HAS2", "CP", "HGF"), log=T)
#ASM
FeaturePlot(Healthy_Mesench, features = c("SOSTDC1", "HHIP", "ACTC1", "PRUNE2"), order = T, pt.size = 1)
VlnPlot(Healthy_Mesench, c("SOSTDC1", "HHIP", "ACTC1", "MYLK", "DSTN", "LPP", "DES", "PRUNE2"), log=T)
#Pericytes
FeaturePlot(Healthy_Mesench, features = c("FAM162B", "POSTN", "HIGD1B", "GUCY1A1", "COX4I1"), order = T, pt.size = 1)
FeaturePlot(Healthy_Mesench, features = c("TRPC6", "PDGFRB", "CSPG4", "MCAM"), order = T, pt.size = 1)
VlnPlot(Healthy_Mesench, c("FAM162B", "POSTN", "HIGD1B", "GUCY1A1", "TRPC6", "PDGFRB", "CSPG4", "MCAM"), log=T)
#SCHWANN CELLS
FeaturePlot(Healthy_Mesench, features = c("SCN7A", "SNCA", "CDH19", "MPZ", "MT3"), order = T, pt.size = 1)
VlnPlot(Healthy_Mesench, c("SCN7A", "SNCA", "CDH19", "MPZ", "MT3"), log=T)
#Mesothelial cells - no distinct cluster
FeaturePlot(Healthy_Mesench, features = c("CDH2", "GPM6A",  "RSPO1"), order = T, pt.size = 1)
VlnPlot(Healthy_Mesench, c("CDH2", "GPM6A",  "RSPO1"), log=T)
#LIPOFIBROBLASTS - no distinct cluster
FeaturePlot(Healthy_Mesench, features = c("LPL", "PLIN2", "APOE"), order = T, pt.size = 1)
VlnPlot(Healthy_Mesench, c("LPL", "PLIN2", "APOE"), log=T)
dev.off()

#####MARKERS OF NEW CLUSTERS 
Markers.SubsetMesench <- FindAllMarkers(Healthy_Mesench, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
write.table(Markers.SubsetMesench, file = "Outs/GLOBAL_PROCESSING/Markers.Healthy_Mesench.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

### Rename all identities
Healthy_Mesench_Renamed <- RenameIdents(object = Healthy_Mesench,
                                                  "0" = "Fibroblasts",
                                                  "1" = "Fibroblasts",
                                                  "2" = "Fibroblasts",
                                                  "3" = "Fibroblasts",
                                                  "4" = "Muscle cells",
                                                  "5" = "Muscle cells",
                                                  "6" = "Fibroblasts",
                                                  "7" = "Fibroblasts",
                                                  "8" = "Doublet",
                                                  "9" = "Fibroblasts",
                                                  "10" = "Muscle cells",
                                                  "11" = "Muscle cells",
                                                  "12" = "Muscle cells",
                                                  "13" = "Schwann cells + Doublets")
DimPlot(Healthy_Mesench_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
Healthy_Mesench_Renamed <- subset(Healthy_Mesench_Renamed, idents = "Doublet", invert = TRUE)
DimPlot(Healthy_Mesench_Renamed, reduction = "umap", label = TRUE, group.by="seurat_clusters", repel = TRUE, pt.size = 1)
saveRDS(Healthy_Mesench_Renamed, "Outs/GLOBAL_PROCESSING/Healthy_Mesench.RDS")

############### SUBSET CONTRACTILE CELLS ###
Healthy_Mesench <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_Mesench.RDS")
DimPlot(Healthy_Mesench, reduction = "umap", label = TRUE, group.by="seurat_clusters", repel = TRUE, pt.size = 1)
Muscle_merged <- subset(Healthy_Mesench, subset = seurat_clusters %in% c(4,5,10,11,12))
Muscle_split <- SplitObject(Muscle_merged, split.by = "orig.ident")
for (i in 1:length(Muscle_split)) {
  Muscle_split[[i]] <- SCTransform(Muscle_split[[i]], vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)
}
integ_features <- SelectIntegrationFeatures(object.list = Muscle_split, nfeatures = 3000) 
Muscle_split <- PrepSCTIntegration(object.list = Muscle_split, anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = Muscle_split, normalization.method = "SCT", anchor.features = integ_features)
Healthy_Muscle <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
Healthy_Muscle <- RunPCA(Healthy_Muscle, npcs = 50, verbose = T)
ElbowPlot(Healthy_Muscle)
Healthy_Muscle <- FindNeighbors(Healthy_Muscle, dims = 1:15)
clustree(Healthy_Muscle, prefix = "integrated_snn_res.")
Healthy_Muscle <- FindClusters(Healthy_Muscle, resolution = 0.9)
Healthy_Muscle <- RunUMAP(Healthy_Muscle, reduction = "pca", dims = 1:15, min.dist = 0.35)
DimPlot(Healthy_Muscle, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)

DefaultAssay(Healthy_Muscle) <- "SCT"
pdf("Outs/GLOBAL_PROCESSING/Healthy_MuscleClusters.pdf")
DimPlot(Healthy_Muscle, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(Healthy_Muscle, reduction = "umap", label = TRUE, repel = TRUE, split.by = "orig.ident", pt.size = 1)
DimPlot(Healthy_Muscle, reduction = "umap", label = TRUE, repel = TRUE, group.by = "orig.ident", pt.size = 1)
VlnPlot(Healthy_Muscle, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
FeaturePlot(Healthy_Muscle, features = c("S.Score","G2M.Score"))
DimPlot(Healthy_Muscle, reduction = "umap", group.by = "initial_annotation", label=T)+NoLegend()
FeaturePlot(Healthy_Muscle, features = c("PECAM1", "EPCAM", "PTPRC"), order = T, pt.size = 1)
dev.off()

# CELL TYPE INVESTIGATION
DefaultAssay(Healthy_Muscle) <- "SCT"
pdf("Outs/GLOBAL_PROCESSING/Healthy_MuscleCellTypes.pdf")
#GENERAL
FeaturePlot(Healthy_Muscle, features = c("ACTA2", "TAGLN", "MYH11", "NOTCH3"), order = T, pt.size = 1)
VlnPlot(Healthy_Muscle, c("ACTA2", "TAGLN", "MYH11", "NOTCH3"), log=T)
FeaturePlot(Healthy_Muscle, features = c("PDGFRA", "COL1A1", "PDGFRB"), order = T, pt.size = 1)
VlnPlot(Healthy_Muscle, c("PDGFRA", "COL1A1", "PDGFRB"), log=T)
#VSM
FeaturePlot(Healthy_Muscle, features = c("FABP3", "MEF2C", "PLN", "SVIL", "ITIH4"), order = T, pt.size = 1)
VlnPlot(Healthy_Muscle, c("FABP3", "MEF2C", "PLN", "SVIL", "ITIH4"), log=T)
#SYSTEMIC PERICYTES
FeaturePlot(Healthy_Muscle, features = c("ADRA2A", "COL12A1", "CLU", "RGS5", "APOA1", "HAS2", "CP", "HGF"), order = T, pt.size = 0.5)
VlnPlot(Healthy_Muscle, c("ADRA2A", "COL12A1", "CLU", "RGS5", "APOA1", "HAS2", "CP", "HGF"), log=T)
#Pericytes
FeaturePlot(Healthy_Muscle, features = c("FAM162B", "POSTN", "HIGD1B", "GUCY1A1"), order = T, pt.size = 1)
FeaturePlot(Healthy_Muscle, features = c("TRPC6", "CSPG4", "MCAM"), order = T, pt.size = 1)
VlnPlot(Healthy_Muscle, c("FAM162B", "POSTN", "HIGD1B","GUCY1A1", "TRPC6", "CSPG4", "MCAM"), log=T)
#ASM
FeaturePlot(Healthy_Muscle, features = c("ACTC1", "MYLK", "DSTN", "LPP", "DES", "PRUNE2"), order = T, pt.size = 1)
VlnPlot(Healthy_Muscle, c("ACTC1", "MYLK", "DSTN", "LPP", "DES", "PRUNE2"), log=T)
#PERIBRONCHIAL FIBROBLASTS/MYOFIBROBLASTS
FeaturePlot(Healthy_Muscle, features = c("SOSTDC1", "ASPN", "HHIP", "FGF18", "CDH4"), order = T, pt.size = 1)
VlnPlot(Healthy_Muscle, c("SOSTDC1", "ASPN", "HHIP", "FGF18", "CDH4", "WIF1"), log=T)
dev.off()

#####MARKERS OF NEW CLUSTERS 
Markers.SubsetMuscle <- FindAllMarkers(Healthy_Muscle, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
write.table(Markers.SubsetMuscle, file = "Outs/GLOBAL_PROCESSING/Markers.Healthy_Muscle.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

#### Rename all identities
Idents(Healthy_Muscle) <- "seurat_clusters"
Healthy_Muscle_Renamed <- RenameIdents(object = Healthy_Muscle,
                                       "0" = "Vascular smooth muscle cells",
                                       "1" = "Pulmonary pericytes",
                                       "2" = "Pulmonary pericytes",
                                       "3" = "Vascular smooth muscle cells",
                                       "4" = "Vascular smooth muscle cells",
                                       "5" = "Systemic pericytes",
                                       "6" = "Pulmonary pericytes",
                                       "7" = "Doublet",
                                       "8" = "Pulmonary pericytes",
                                       "9" = "Peribronchial myofibroblasts",
                                       "10" = "Airway smooth muscle cells",
                                       "11" = "Myofibroblasts")
DimPlot(Healthy_Muscle_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
# Remove aberrant cells
Healthy_Muscle_Renamed <- subset(Healthy_Muscle_Renamed, idents = "Doublet", invert = TRUE)
DimPlot(Healthy_Muscle_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
# RE-RUN UMAP AFTER DOUBLET EXCLUSION
Healthy_Muscle_Renamed <- RunUMAP(Healthy_Muscle_Renamed, reduction = "pca", dims = 1:15, min.dist = 0.35)
DimPlot(Healthy_Muscle_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)
# store final annotation in metadata
Healthy_Muscle_Renamed$celltype <- Idents(Healthy_Muscle_Renamed)
saveRDS(Healthy_Muscle_Renamed, "Outs/GLOBAL_PROCESSING/Healthy_Muscle.RDS")

###FIND MARKERS OF CELL TYPES
Markers.SubsetMuscleR <- FindAllMarkers(Healthy_Muscle_Renamed, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
write.table(Markers.SubsetMuscleR, file = "Outs/FINAL_MARKERS/Muscle_Cells.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

##################### SUBSET FIBROBLASTES ###
Healthy_Mesench <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_Mesench.RDS")
Fib_merged <- subset(Healthy_Mesench, subset = seurat_clusters %in% c(0,1,2,3,6,7,9))
Fib_split <- SplitObject(Fib_merged, split.by = "orig.ident")
for (i in 1:length(Fib_split)) {
  Fib_split[[i]] <- SCTransform(Fib_split[[i]], vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)
}
integ_features <- SelectIntegrationFeatures(object.list = Fib_split, nfeatures = 3000) 
Fib_split <- PrepSCTIntegration(object.list = Fib_split, anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = Fib_split, normalization.method = "SCT", anchor.features = integ_features)
Healthy_Fib <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
Healthy_Fib <- RunPCA(Healthy_Fib, npcs = 50, verbose = T)
ElbowPlot(Healthy_Fib)
Healthy_Fib <- FindNeighbors(Healthy_Fib, dims = 1:5)
clustree(Healthy_Fib, prefix = "integrated_snn_res.")
Healthy_Fib <- FindClusters(Healthy_Fib, resolution = 0.7)
Healthy_Fib <- RunUMAP(Healthy_Fib, reduction = "pca", dims = 1:5, min.dist = 0.35)
DimPlot(Healthy_Fib, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

DefaultAssay(Healthy_Fib) <- "SCT"
pdf("Outs/GLOBAL_PROCESSING/Healthy_FibClusters.pdf")
DimPlot(Healthy_Fib, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(Healthy_Fib, reduction = "umap", label = TRUE, repel = TRUE, split.by = "orig.ident", pt.size = 1)
DimPlot(Healthy_Fib, reduction = "umap", label = TRUE, repel = TRUE, group.by = "orig.ident", pt.size = 1)
VlnPlot(Healthy_Fib, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
FeaturePlot(Healthy_Fib, features = c("S.Score","G2M.Score")) # no replicating cluster
DimPlot(Healthy_Fib, reduction = "umap", group.by = "initial_annotation", label=T)+NoLegend()
FeaturePlot(Healthy_Fib, features = c("PECAM1", "EPCAM", "PTPRC"), order = T, pt.size = 1) # no doublets
dev.off()

# CELL TYPE INVESTIGATION
DefaultAssay(Healthy_Fib) <- "SCT"
pdf("Outs/GLOBAL_PROCESSING/Healthy_FibCellTypes.pdf")
FeaturePlot(Healthy_Fib, features = c("PDGFRA", "COL1A1", "PDGFRB"), order = T, pt.size = 1)
VlnPlot(Healthy_Fib, c("PDGFRA", "COL1A1", "PDGFRB"), log=T)
###ADV vs ALV
FeaturePlot(Healthy_Fib, features = c("COL14A1", "GLI1", "COL13A1", "WNT2"), order = T, pt.size = 1)
VlnPlot(Healthy_Fib, c("COL14A1", "GLI1", "COL13A1", "WNT2"), log=T)
#ADVENTITIAL FIB
FeaturePlot(Healthy_Fib, features = c("PI16", "DCN", "AQP1", "LSP1", "LTBP1"), order = T, pt.size = 1)
VlnPlot(Healthy_Fib, c("PI16", "DCN", "AQP1", "LSP1", "LTBP1"), log=T)
#ALVEOLAR fib
FeaturePlot(Healthy_Fib, features = c("NPNT", "MACF1", "DST", "SPECC1L", "PRG4"), order = T, pt.size = 1)
VlnPlot(Healthy_Fib, c("NPNT", "MACF1", "DST", "SPECC1L", "PRG4"), log=T)
# 'IMMUNE-RECRUITING' FIBROBLASTS ?
FeaturePlot(Healthy_Fib, features = c("CCL19", "CXCL12"), label=T, order = T, pt.size = 1)
VlnPlot(Healthy_Fib, c("CCL19", "CXCL12"), log=T)
# 'PERIBRONCHIAL FIBROBLASTS'?
FeaturePlot(Healthy_Fib, features = c("TMEM132C", "COL15A1", "ENTPD1", "PLCL1"), order = T, pt.size = 1)
VlnPlot(Healthy_Fib, c("TMEM132C", "COL15A1", "ENTPD1", "PLCL1", "ATRNL1","PLPPR4"), log=T)
VlnPlot(Healthy_Fib, c("F13A1", "FGF14", "CHN1", "NTRK3", "SUGCT", "PRDM6", "ENOX1"), log=T)
dev.off()

###FIND MARKERS OF CELL TYPES
Markers.SubsetFib <- FindAllMarkers(Healthy_Fib, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
write.table(Markers.SubsetFib, file = "Outs/GLOBAL_PROCESSING/Markers.Fib.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

#### Rename all identities
Healthy_Fib_Renamed <- RenameIdents(object = Healthy_Fib,
                                    "0" = "Alveolar fibroblasts",
                                    "1" = "Alveolar fibroblasts",
                                    "2" = "Alveolar fibroblasts",
                                    "3" = "Alveolar fibroblasts",
                                    "4" = "STMN2+ alveolar fibroblasts",
                                    "5" = "Alveolar fibroblasts",
                                    "6" = "COL23A1+ adventitial fibroblasts",
                                    "7" = "CCBE1+ adventitial fibroblasts",
                                    "8" = "CCL19+ adventitial fibroblasts",
                                    "9" = "CCN3+ adventitial fibroblasts")
DimPlot(Healthy_Fib_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1.5)+NoLegend()
# store final annotation in metadata
Healthy_Fib_Renamed$celltype <- Idents(Healthy_Fib_Renamed)
saveRDS(Healthy_Fib_Renamed, "Outs/GLOBAL_PROCESSING/Healthy_Fib.RDS")

###FIND MARKERS OF CELL TYPES
Markers.SubsetFibR <- FindAllMarkers(Healthy_Fib_Renamed, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
write.table(Markers.SubsetFibR, file = "Outs/FINAL_MARKERS/Fib_Cells.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

##################### SUBSET SCHWANN CELLS ###
Healthy_Mesench <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_Mesench.RDS")
Schwann_isolated <- subset(Healthy_Mesench, subset = seurat_clusters %in% c(13))
DimPlot(Schwann_isolated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DefaultAssay(Schwann_isolated) <- "integrated"
Schwann_isolated <- RunPCA(Schwann_isolated, npcs = 20, verbose = T)
Schwann_isolated <- FindNeighbors(Schwann_isolated, dims = 1:8)
Schwann_isolated <- FindClusters(Schwann_isolated, resolution = 1.5)
Schwann_isolated <- RunUMAP(Schwann_isolated, reduction = "pca", dims = 1:8, min.dist = 0.15)
DimPlot(Schwann_isolated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

# separating schwann cells (cl 2) from doublets (cl 1,3,4)
DefaultAssay(Schwann_isolated) <- "SCT"
VlnPlot(Schwann_isolated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
FeaturePlot(Schwann_isolated, features = c("EPCAM", "PTPRC", "PECAM1", "AGER"), order = T, pt.size = 2)
FeaturePlot(Schwann_isolated, features = c("SCN7A", "SNCA", "CDH19", "MPZ"), order = T, pt.size = 2) # "SNCA", "CDH19", "MPZ"

Schwann_isolated_renamed <- RenameIdents(object = Schwann_isolated, 
                                         "0" = "Doublet",
                                         "1" = "Schwann cells",
                                         "2" = "Doublet",
                                         "3" = "Doublet")
DimPlot(Schwann_isolated_renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
# keeping doublets for know to avoid issues with future merging
# store final annotation in metadata
Schwann_isolated_renamed$celltype <- Idents(Schwann_isolated_renamed)
saveRDS(Schwann_isolated_renamed, "Outs/GLOBAL_PROCESSING/Schwann_Doublets.RDS")

##################### RE-INTEGRATION ALL MESENCHYMAL CELLS ###
Healthy_Muscle_Renamed <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_Muscle.RDS")
Healthy_Fib_Renamed <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_Fib.RDS")
Schwann_isolated_renamed <- readRDS("Outs/GLOBAL_PROCESSING/Schwann_Doublets.RDS")

Mesench_merged <- merge(Healthy_Fib_Renamed, y = c(Healthy_Muscle_Renamed, Schwann_isolated_renamed), 
                       add.cell.ids = c("Fibroblasts", "Muscle", "Schwann"), project = "MESENCH")
Mesench_integrated_split <- SplitObject(Mesench_merged, split.by = "orig.ident")
for (i in 1:length(Mesench_integrated_split)) {
  Mesench_integrated_split[[i]] <- SCTransform(Mesench_integrated_split[[i]], vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)
}
integ_features <- SelectIntegrationFeatures(object.list = Mesench_integrated_split, nfeatures = 3000)
Mesench_integrated_split <- PrepSCTIntegration(object.list = Mesench_integrated_split, anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = Mesench_integrated_split, normalization.method = "SCT", anchor.features = integ_features)
Mesench_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
Mesench_integrated <- RunPCA(Mesench_integrated, npcs = 50, verbose = T)
Mesench_integrated <- RunUMAP(Mesench_integrated, reduction = "pca", dims = 1:15, min.dist = 0.35)
DimPlot(Mesench_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.05)+ NoLegend()
# exclusion of previously identified remaining doublets
Mesench_integrated <- subset(Mesench_integrated, idents = "Doublet", invert = TRUE)
Mesench_integrated <- RunUMAP(Mesench_integrated, reduction = "pca", dims = 1:15, min.dist = 0.4)
DimPlot(Mesench_integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.05)+ NoLegend()
saveRDS(Mesench_integrated, "Outs/GLOBAL_PROCESSING/Mesench_integrated.RDS")

# Adding group info to metadata
head(colnames(Mesench_integrated))
cellgroup <- sapply(X = strsplit(colnames(Mesench_integrated), split = "_"), FUN = "[", 1)
Mesench_integrated$cellgroup <- cellgroup
DimPlot(Mesench_integrated, reduction = "umap", label = TRUE, repel = TRUE, group.by= "cellgroup", pt.size = 0.05)+ NoLegend()
saveRDS(Mesench_integrated, "Outs/GLOBAL_PROCESSING/Mesench_integrated.RDS")

#### EPITHELIAL COMPARTMENT ####
Epith_merged <- subset(seurat_integrated, subset = seurat_clusters %in% c(16,17,37,42,47,53))
Epith_split <- SplitObject(Epith_merged, split.by = "orig.ident")
for (i in 1:length(Epith_split)) {
  Epith_split[[i]] <- SCTransform(Epith_split[[i]], vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)
}
integ_features <- SelectIntegrationFeatures(object.list = Epith_split, nfeatures = 3000)
Epith_split <- PrepSCTIntegration(object.list = Epith_split, anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = Epith_split, normalization.method = "SCT", anchor.features = integ_features)
Healthy_Epith <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
Healthy_Epith <- RunPCA(Healthy_Epith, npcs = 50, verbose = T)
ElbowPlot(Healthy_Epith)
Healthy_Epith <- FindNeighbors(Healthy_Epith, dims = 1:8)
clustree(Healthy_Epith, prefix = "integrated_snn_res.")
Healthy_Epith <- FindClusters(Healthy_Epith, resolution = 0.7)
Healthy_Epith <- RunUMAP(Healthy_Epith, reduction = "pca", dims = 1:8, min.dist = 0.35)
DimPlot(Healthy_Epith, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

DefaultAssay(Healthy_Epith) <- "SCT"
pdf("Outs/GLOBAL_PROCESSING/Healthy_EpithClusters.pdf")
DimPlot(Healthy_Epith, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(Healthy_Epith, reduction = "umap", label = TRUE, repel = TRUE, split.by = "orig.ident", pt.size = 1.5)
DimPlot(Healthy_Epith, reduction = "umap", label = TRUE, repel = TRUE, group.by = "orig.ident", pt.size = 1)
VlnPlot(Healthy_Epith, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
DimPlot(Healthy_Epith, reduction = "umap", group.by = "initial_annotation", label=T)+NoLegend()
FeaturePlot(Healthy_Epith, features = c("S.Score","G2M.Score"))
FeaturePlot(Healthy_Epith, features = c("PECAM1", "EPCAM", "PTPRC", "DCN"), order = T, pt.size = 1)
VlnPlot(Healthy_Epith, c("EPCAM"), log=T)
VlnPlot(Healthy_Epith, c("PECAM1"), log=T)
VlnPlot(Healthy_Epith, c("PTPRC"), log=T)
VlnPlot(Healthy_Epith, c("DCN"), log=T)
dev.off()

# cl 10 are doublets, exclusion of cluster 10 and re-run processing
Healthy_Epith <- subset(Healthy_Epith, subset = seurat_clusters %in% c(10), invert=T)
DimPlot(Healthy_Epith, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

# CELL TYPE INVESTIGATION
DefaultAssay(Healthy_Epith) <- "SCT"
pdf("Outs/GLOBAL_PROCESSING/Healthy_EpithCellTypes.pdf")
#CLUB CELLS
FeaturePlot(Healthy_Epith, features = c("SCGB1A1", "AQP5", "SERPINA1", "ERN2", "DEPP1", "GPX2"), label=T, order = T, pt.size = 1)
VlnPlot(Healthy_Epith, features = c("SCGB1A1", "AQP5", "SERPINA1", "ERN2", "DEPP1", "GPX2"), log=T)
#Ionocytes, Goblet
FeaturePlot(Healthy_Epith, features = c("CFTR", "MUC5B", "SPDEF"), order = T, pt.size = 1)
VlnPlot(Healthy_Epith, features = c("CFTR", "MUC5B", "SPDEF"), log=T)
#AT1
FeaturePlot(Healthy_Epith, features = c("AGER", "AQP1", "SUSD2", "CLDN18"), order = T, pt.size = 1)
VlnPlot(Healthy_Epith, features = c("AGER", "AQP1", "SUSD2", "CLDN18"), log=T)
#AT2
FeaturePlot(Healthy_Epith, features = c("NAPSA", "SFTPC", "SFTPB", "MUC1"), order = T, pt.size = 1)
VlnPlot(Healthy_Epith, features = c("NAPSA", "SFTPC", "SFTPB", "MUC1"), log=T)
#BASAL
FeaturePlot(Healthy_Epith, features = c("KRT14", "TP63"), order = T, pt.size = 1)
VlnPlot(Healthy_Epith, features = c("KRT14", "TP63"), log=T)
#CILIATED
FeaturePlot(Healthy_Epith, features = c("CAPS", "FOXJ1", "CCDC78", "HYDIN"), order = T, pt.size = 1)
VlnPlot(Healthy_Epith, features = c("CAPS", "FOXJ1", "CCDC78", "HYDIN"), log=T)
dev.off()

#####MARKERS OF NEW CLUSTERS 
Markers.Epith <- FindAllMarkers(Healthy_Epith, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
write.table(Markers.Epith, file = "Outs/GLOBAL_PROCESSING/Markers.Healthy_Epith.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

#### Rename all identities
Healthy_Epith_Renamed <- RenameIdents(object = Healthy_Epith,
                                                "0" = "Alveolar type 2", "1" = "Alveolar type 2", "2" = "Alveolar type 2",
                                                "3" = "Secretory","4" = "Alveolar type 2", "5" = "Basal",
                                                "6" = "Alveolar type 2", "7" = "Alveolar type 2", "8" = "Alveolar type 1",
                                                "9" = "Ciliated")
DimPlot(Healthy_Epith_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)
# store final annotation in metadata
Healthy_Epith_Renamed$celltype <- Idents(Healthy_Epith_Renamed)
saveRDS(Healthy_Epith_Renamed, "Outs/GLOBAL_PROCESSING/Healthy_Epith.RDS")

###FIND MARKERS OF CELL TYPES
Markers.EpithR <- FindAllMarkers(Healthy_Epith_Renamed, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
write.table(Markers.EpithR, file = "Outs/FINAL_MARKERS/Epithelial_Cells.CSV", sep = "\t", dec = ".", col.names = T, row.names = F)

### INVESTIGATION OF CANCER MARKERS ###
#rename origin
Healthy_Epith_Renamed <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_Epith.RDS")
new.origin.ids <- c("ADJACENT", "HEALTHY", "HEALTHY", "ADJACENT")
Healthy_Epith_Renamed$origin <- Healthy_Epith_Renamed$orig.ident
Healthy_Epith_Renamed$origin <- as.factor(Healthy_Epith_Renamed$origin)
levels(Healthy_Epith_Renamed@meta.data$origin)<-new.origin.ids
DimPlot(Healthy_Epith_Renamed, reduction = "umap", label = TRUE, split.by= "origin", group.by="orig.ident", pt.size = 2)
saveRDS(Healthy_Epith_Renamed, "Outs/GLOBAL_PROCESSING/Healthy_Epith.RDS")

DefaultAssay(Healthy_Epith_Renamed) <- "RNA"
Idents(Healthy_Epith_Renamed) <- "origin"
ADJTUMOR.effect <- FindMarkers(Healthy_Epith_Renamed, ident.1 = "ADJACENT", ident.2 = "HEALTHY", only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
write.table(ADJTUMOR.effect, file = "Outs/FINAL_MARKERS/DEG_ADJTUMOR_EPITH.CSV", sep = "\t", dec = ".", col.names = T, row.names = T)

#### IMMUNE COMPARTMENT ####
Immune_merged <- subset(seurat_integrated, subset = seurat_clusters %in% c(0,4,15,19,26,28,50,51,46,48,10,40,8,12,35,36,20,52))
Immune_split <- SplitObject(Immune_merged, split.by = "orig.ident")
for (i in 1:length(Immune_split)) {
  Immune_split[[i]] <- SCTransform(Immune_split[[i]], vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)
}
integ_features <- SelectIntegrationFeatures(object.list = Immune_split, nfeatures = 3000)
Immune_split <- PrepSCTIntegration(object.list = Immune_split, anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = Immune_split, normalization.method = "SCT", anchor.features = integ_features)
Healthy_Immune <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
Healthy_Immune <- RunPCA(Healthy_Immune, npcs = 50, verbose = T)
ElbowPlot(Healthy_Epith)
Healthy_Immune <- FindNeighbors(Healthy_Immune, dims = 1:8)
Healthy_Immune <- FindClusters(Healthy_Immune, resolution = 1)
Healthy_Immune <- RunUMAP(Healthy_Immune, reduction = "pca", dims = 1:8, min.dist = 0.25)
DimPlot(Healthy_Immune, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

DefaultAssay(Healthy_Immune) <- "SCT"
pdf("Outs/GLOBAL_PROCESSING/Healthy_ImmuneClusters.pdf")
DimPlot(Healthy_Immune, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(Healthy_Immune, reduction = "umap", label = TRUE, repel = TRUE, split.by = "orig.ident", pt.size = 1)
DimPlot(Healthy_Immune, reduction = "umap", label = TRUE, repel = TRUE, group.by = "orig.ident", pt.size = 1)
VlnPlot(Healthy_Immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
VlnPlot(Healthy_Immune, features = c("percent.RPS", "percent.RPL", "percent.mt"), ncol = 3, pt.size = F, log = F)
DimPlot(Healthy_Immune, reduction = "umap", group.by = "initial_annotation", label=T)+NoLegend()
FeaturePlot(Healthy_Immune, features = c("S.Score","G2M.Score"))
FeaturePlot(Healthy_Immune, features = c("PECAM1", "EPCAM", "PTPRC", "COL1A1"), order = T, pt.size = 1) 
VlnPlot(Healthy_Immune, c("PECAM1", "EPCAM", "PTPRC", "DCN"), log=T)
dev.off()

# Global cell type identification 
DefaultAssay(Healthy_Immune) <- "SCT"
pdf("Outs/GLOBAL_PROCESSING/Healthy_ImmuneCellTypes.pdf")
#LYMPHOCYTES
FeaturePlot(Healthy_Immune, features = c("CD3D", "CD3E", "IL7R", "ICOS", "THY1"), order = T, pt.size = 1)
VlnPlot(Healthy_Immune, features = c("CD3D", "CD3E", "IL7R", "ICOS", "THY1"), log=T)
#CTL
FeaturePlot(Healthy_Immune, features = c("CD3D", "CD8A", "NCR3", "GZMB", "KLRB1", "KLRK1", "CCL5"), order = T, pt.size = 1)
VlnPlot(Healthy_Immune, features = c("CD3D", "CD8A", "NCR3", "GZMB", "KLRB1", "KLRK1", "CCL5"), log=T)
#NAIVE LT
FeaturePlot(Healthy_Immune, features = c("CD3D", "LEF1", "TCF7", "IL7R", "CCR7", "SELL", "LTB"), order = T, pt.size = 1)
VlnPlot(Healthy_Immune, features = c("CD3D", "LEF1", "TCF7", "IL7R", "CCR7", "SELL", "LTB"), log=T)
#LB
FeaturePlot(Healthy_Immune, features = c("MS4A1", "CD19", "FCRLA", "CCR7"), order = T, pt.size = 1)
VlnPlot(Healthy_Immune, features = c("MS4A1", "CD19", "FCRLA", "CCR7"), log=T)
#PLASMA CELLS
FeaturePlot(Healthy_Immune, features = c("TNFRSF17", "JCHAIN", "POU2AF1"), order = T, pt.size = 1)
VlnPlot(Healthy_Immune, features = c("TNFRSF17", "JCHAIN", "POU2AF1"), log=T)
#T REGULATORY
FeaturePlot(Healthy_Immune, features = c("CD3D", "CTLA4", "IKZF2", "GATA3", "CCL5", "LGALS3"), order = T, pt.size = 1)
VlnPlot(Healthy_Immune, features = c("CD3D", "CTLA4", "IKZF2", "GATA3", "CCL5", "LGALS3"), log=T)
#MAST CELLS
FeaturePlot(Healthy_Immune, features = c("KIT", "CPA3", "MS4A2", "FCER1A"), order = T, pt.size = 1)
VlnPlot(Healthy_Immune, features = c("KIT", "CPA3", "MS4A2", "FCER1A"), log=T)
#NEUTRO
FeaturePlot(Healthy_Immune, features = c("S100A12", "ITGAM", "CXCL8", "ADGRE1", "CD4", "SELL", "PADI3"), order = T, pt.size = 1)
VlnPlot(Healthy_Immune, features = c("S100A12", "ITGAM", "CXCL8", "ADGRE1", "CD4", "SELL", "PADI3"), log=T)
#GENERAL MP
FeaturePlot(Healthy_Immune, features = c("MRC1", "MSR1", "APOC1", "CD86", "CD48", "CD63"), order = T, pt.size = 1)
VlnPlot(Healthy_Immune, features = c("MRC1", "MSR1", "APOC1", "CD86", "CD48", "CD63"), log=T)
#AM
FeaturePlot(Healthy_Immune, features = c("MARCO", "MSR1","MRC1", "CD163", "CHI3L1", "APOC1"), order = T, pt.size = 1)
VlnPlot(Healthy_Immune, features = c("MARCO", "MSR1","MRC1", "CD163", "CHI3L1", "APOC1"), log=T)
##CCL13 MP
FeaturePlot(Healthy_Immune, features = c("C1QC", "C1QA", "CCL13", "STAB1", "NPAS1", "MRC1"), order = T, pt.size = 1)
VlnPlot(Healthy_Immune, features = c("C1QC", "C1QA", "CCL13", "STAB1", "NPAS1", "MRC1"), log=T)
#MONOCYTES
FeaturePlot(Healthy_Immune, features = c("LYZ", "PLEK", "VCAN", "FN1", "IL1B", "MRC1"), order = T, pt.size = 1)
VlnPlot(Healthy_Immune, features = c("LYZ", "PLEK", "VCAN", "FN1", "IL1B", "MRC1"), log=T)
#DC
FeaturePlot(Healthy_Immune, features = c("CD1C", "PLD4", "FCER1A",  "PPM1J", "IRF8", "DLA-DOA", "PKIB"), order = T, pt.size = 1)
VlnPlot(Healthy_Immune, features = c("CD1C", "PLD4", "FCER1A",  "PPM1J", "IRF8", "DLA-DOA", "PKIB"), log=T)
#MATURE DC
FeaturePlot(Healthy_Immune, features = c("CD1C","CD83", "IDO1", "IL4I1", "CCR7"), order = T, pt.size = 1)
VlnPlot(Healthy_Immune, features = c("CD1C","CD83", "IDO1", "IL4I1", "CCR7"), log=T)
dev.off()

##### Rename all identities
Healthy_Immune_Renamed <- RenameIdents(object = Healthy_Immune,
                                                 "0" = "Lymphoid", "1" = "Myeloid", "2" = "Myeloid",
                                                 "3" = "Myeloid", "4" = "Myeloid",
                                                 "5" = "Lymphoid", "6" = "Myeloid",
                                                 "7" = "Myeloid", "8" = "Myeloid",
                                                 "9" = "Myeloid", "10" = "Myeloid",
                                                 "11" = "Myeloid", "12" = "Lymphoid",
                                                 "13" = "Lymphoid", "14" = "Myeloid",
                                                 "15" = "Lymphoid", "16" = "Lymphoid",
                                                 "17" = "Doublet")
DimPlot(Healthy_Immune_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)
# Remove aberrant cells
Healthy_Immune_Renamed <- subset(Healthy_Immune_Renamed, idents = c("Doublet"), invert = TRUE)
DimPlot(Healthy_Immune_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
saveRDS(Healthy_Immune_Renamed, "Outs/GLOBAL_PROCESSING/Healthy_Immune.RDS")

##################### MYELOID CELLS ###
Subset_Myeloid <- subset(Healthy_Immune_Renamed, subset = seurat_clusters %in% c(1,2,3,4,6,7,8,9,10,11,14))
Myeloid_split <- SplitObject(Subset_Myeloid, split.by = "orig.ident")
for (i in 1:length(Myeloid_split)) {
  Myeloid_split[[i]] <- SCTransform(Myeloid_split[[i]], vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)
}
integ_features <- SelectIntegrationFeatures(object.list = Myeloid_split, nfeatures = 3000)
Myeloid_split <- PrepSCTIntegration(object.list = Myeloid_split, anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = Myeloid_split, normalization.method = "SCT", anchor.features = integ_features)
Healthy_Myeloid <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
Healthy_Myeloid <- RunPCA(Healthy_Myeloid, npcs = 60, verbose = T)
ElbowPlot(Healthy_Myeloid)
clustree(Healthy_Myeloid, prefix = "integrated_snn_res.")
Healthy_Myeloid <- FindNeighbors(Healthy_Myeloid, dims = 1:60)
Healthy_Myeloid <- FindClusters(Healthy_Myeloid, resolution = 1.8)
Healthy_Myeloid <- RunUMAP(Healthy_Myeloid, reduction = "pca", dims = 1:60, min.dist = 0.35)
DimPlot(Healthy_Myeloid, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)

DefaultAssay(Healthy_Myeloid) <- "SCT"
pdf("Outs/GLOBAL_PROCESSING/Healthy_MyeloidClusters.pdf")
DimPlot(Healthy_Myeloid, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(Healthy_Myeloid, reduction = "umap", label = TRUE, repel = TRUE, split.by = "orig.ident", pt.size = 1)
DimPlot(Healthy_Myeloid, reduction = "umap", label = TRUE, repel = TRUE, group.by = "orig.ident", pt.size = 1)
VlnPlot(Healthy_Myeloid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
VlnPlot(Healthy_Myeloid, features = c("percent.RPS", "percent.RPL", "percent.mt"), ncol = 3, pt.size = F, log = F)
FeaturePlot(Healthy_Myeloid, features = c("percent.mt"), order = T, pt.size = 1) 
DimPlot(Healthy_Myeloid, reduction = "umap", group.by = "initial_annotation", label=T)+NoLegend()
FeaturePlot(Healthy_Myeloid, features = c("S.Score","G2M.Score"))
FeaturePlot(Healthy_Myeloid, features = c("PECAM1", "EPCAM", "DCN"), order = T, pt.size = 1) 
VlnPlot(Healthy_Myeloid, c("PECAM1", "EPCAM", "PTPRC", "DCN"), log=T)
dev.off()

# cluster 5 and 11 are low quality, filtering out before re-clustering
Healthy_Myeloid <- subset(Healthy_Myeloid, subset = seurat_clusters %in% c(5,11), invert=T)
DimPlot(Healthy_Myeloid, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

DefaultAssay(Healthy_Myeloid) <- "SCT"
pdf("Outs/GLOBAL_PROCESSING/Healthy_MyeloidCellTypes.pdf")
##MAST CELLS
FeaturePlot(Healthy_Myeloid, features = c("KIT", "CPA3", "MS4A2", "FCER1A"), order = T, pt.size = 1)
VlnPlot(Healthy_Myeloid, features = c("KIT", "CPA3", "MS4A2", "FCER1A"), log=T)
###NEUTRO
FeaturePlot(Healthy_Myeloid, features = c("S100A12", "ITGAM", "CXCL8", "ADGRE1", "CD4", "SELL", "PADI3"), order = T, pt.size = 1)
VlnPlot(Healthy_Myeloid, features = c("S100A12", "ITGAM", "CXCL8", "ADGRE1", "CD4", "SELL", "PADI3"), log=T)
##GENERAL MP
FeaturePlot(Healthy_Myeloid, features = c("MRC1", "MSR1", "APOC1", "CD86", "CD48", "CD63"), order = T, pt.size = 1)
VlnPlot(Healthy_Myeloid, features = c("MRC1", "MSR1", "APOC1", "CD86", "CD48", "CD63"), log=T)
#AM
FeaturePlot(Healthy_Myeloid, features = c("MARCO", "MSR1","MRC1", "CD163", "CHI3L1", "APOC1"), order = T, pt.size = 1)
VlnPlot(Healthy_Myeloid, features = c("MARCO", "MSR1","MRC1", "CD163", "CHI3L1", "APOC1"), log=T)
#CCL13 MP
FeaturePlot(Healthy_Myeloid, features = c("C1QC", "C1QA", "CCL13", "STAB1", "NPAS1", "MRC1"), order = T, pt.size = 1)
VlnPlot(Healthy_Myeloid, features = c("C1QC", "C1QA", "CCL13", "STAB1", "NPAS1", "MRC1"), log=T)
#MONOCYTES
FeaturePlot(Healthy_Myeloid, features = c("LYZ", "PLEK", "VCAN", "FN1", "IL1B", "MRC1"), order = T, pt.size = 1)
VlnPlot(Healthy_Myeloid, features = c("LYZ", "PLEK", "VCAN", "FN1", "IL1B", "MRC1"), log=T)
#DC
FeaturePlot(Healthy_Myeloid, features = c("PLD4", "FCER1A",  "PPM1J", "DLA-DOA", "PKIB"), order = T, pt.size = 1)
VlnPlot(Healthy_Myeloid, features = c("PLD4", "FCER1A",  "PPM1J", "DLA-DOA", "PKIB"), log=T)
#mcDC1
FeaturePlot(Healthy_Myeloid, features = c("ECRG4", "IRF8", "BATF3", "CADM1"), order = T, pt.size = 1)
VlnPlot(Healthy_Myeloid, features = c("ECRG4", "IRF8", "BATF3", "CADM1"), log=T)
#mcDC2
FeaturePlot(Healthy_Myeloid, features = c("CD1C", "PKIB", "IRF4", "ITGAM", "IRF2"), order = T, pt.size = 1)
VlnPlot(Healthy_Myeloid, features = c("CD1C", "PKIB", "IRF4", "ITGAM", "IRF2"), log=T)
#pDC
FeaturePlot(Healthy_Myeloid, features = c("IL3RA", "NRP1", "IRF8", "FCER1A"), order = T, pt.size = 1)
VlnPlot(Healthy_Myeloid, features = c("IL3RA", "NRP1", "IRF8", "FCER1A"), log=T)
#MATURE DC
FeaturePlot(Healthy_Myeloid, features = c("CD1C","CD83", "IDO1", "IL4I1", "CCR7"), order = T, pt.size = 1)
VlnPlot(Healthy_Myeloid, features = c("CD1C","CD83", "IDO1", "IL4I1", "CCR7"), log=T)
dev.off()

###MARKERS MYELOID
Markers.SubsetMyeloid <- FindAllMarkers(Healthy_Myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
write.table(Markers.SubsetMyeloid, file = "Outs/GLOBAL_PROCESSING/Markers.Healthy_Myeloid.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

##### Rename all identities
Idents(Healthy_Myeloid) <- "seurat_clusters"
Healthy_Myeloid_Renamed <- RenameIdents(object = Healthy_Myeloid,
                                       "0" = "Neutrophils", "1" = "AM", "2" = "FN1+ monocytes",
                                       "3" = "Mast cells","4" = "AM","5" = "cDC2",
                                       "6" = "Monocytes","7" = "cDC2","8" = "CD1C+ monocytes",
                                       "9" = "Neutrophils", "10" = "CCL13+ macrophages","11" = "CCL13+ macrophages",
                                       "12" = "C1Q+ AM","13" = "FN1+ monocytes","14" = "CCL13+ macrophages",
                                       "15" = "Neutrophils","16" = "Monocytes","17" = "Mast cells",
                                       "18" = "SLC27A6+ macrophages","19" = "Mature DC","20" = "cDC1", 
                                       "21"="Neutrophils", "22"="pDC")
DimPlot(Healthy_Myeloid_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)+NoLegend()
# store final annotation in metadata
Healthy_Myeloid_Renamed$celltype <- Idents(Healthy_Myeloid_Renamed)
saveRDS(Healthy_Myeloid_Renamed, "Outs/GLOBAL_PROCESSING/Healthy_Myeloid.RDS")

###FIND MARKERS OF CELL TYPES
Markers.SubsetMyeloidR <- FindAllMarkers(Healthy_Myeloid_Renamed, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
write.table(Markers.SubsetMyeloidR, file = "Outs/FINAL_MARKERS/Myeloid_Cellsres1.8.CSV", sep = "\t", dec = ".", col.names = T, row.names = F)

##################### LYMPHOID CELLS ###
Healthy_Immune <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_Immune.RDS")
Subset_Lymphoid <- subset(Healthy_Immune, subset = seurat_clusters %in% c(0,5,12,13,15,16))
Lymphoid_split <- SplitObject(Subset_Lymphoid, split.by = "orig.ident")
for (i in 1:length(Lymphoid_split)) {
  Lymphoid_split[[i]] <- SCTransform(Lymphoid_split[[i]], vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)
}
integ_features <- SelectIntegrationFeatures(object.list = Lymphoid_split, nfeatures = 3000)
Lymphoid_split <- PrepSCTIntegration(object.list = Lymphoid_split, anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = Lymphoid_split, normalization.method = "SCT", anchor.features = integ_features)
Healthy_Lymphoid <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
Healthy_Lymphoid <- RunPCA(Healthy_Lymphoid, npcs = 60, verbose = T)
ElbowPlot(Healthy_Lymphoid)
Healthy_Lymphoid <- FindNeighbors(Healthy_Lymphoid, dims = 1:50)
clustree(Healthy_Lymphoid, prefix = "integrated_snn_res.")
Healthy_Lymphoid <- FindClusters(Healthy_Lymphoid, resolution = 1.5)
Healthy_Lymphoid <- RunUMAP(Healthy_Lymphoid, reduction = "pca", dims = 1:50, min.dist = 0.35)
DimPlot(Healthy_Lymphoid, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

DefaultAssay(Healthy_Lymphoid) <- "SCT"
pdf("Outs/GLOBAL_PROCESSING/Healthy_LymphoidClusters.pdf")
DimPlot(Healthy_Lymphoid, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(Healthy_Lymphoid, reduction = "umap", label = TRUE, repel = TRUE, split.by = "orig.ident", pt.size = 1)
DimPlot(Healthy_Lymphoid, reduction = "umap", label = TRUE, repel = TRUE, group.by = "orig.ident", pt.size = 1)
VlnPlot(Healthy_Lymphoid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
VlnPlot(Healthy_Lymphoid, features = c("percent.RPS", "percent.RPL", "percent.mt"), ncol = 3, pt.size = F, log = F)
DimPlot(Healthy_Lymphoid, reduction = "umap", group.by = "initial_annotation", label=T)+NoLegend()
FeaturePlot(Healthy_Lymphoid, features = c("S.Score","G2M.Score"))
FeaturePlot(Healthy_Lymphoid, features = c("PECAM1", "EPCAM", "DCN"), order = T, pt.size = 1)
VlnPlot(Healthy_Lymphoid, c("PECAM1", "EPCAM", "DCN"), log=T)
dev.off()

# CELL TYPE INVESTIGATION
DefaultAssay(Healthy_Lymphoid) <- "SCT"
pdf("Outs/GLOBAL_PROCESSING/Healthy_LymphoidCellTypes.pdf")
#LYMPHOCYTES
FeaturePlot(Healthy_Lymphoid, features = c("CD3D", "CD3E", "IL7R", "ICOS", "THY1"), order = T, pt.size = 1)
VlnPlot(Healthy_Lymphoid, features = c("CD3D", "CD3E", "IL7R", "ICOS", "THY1"), log=T)
#CTL
FeaturePlot(Healthy_Lymphoid, features = c("CD3D", "CD8A", "NCR3", "GZMB", "KLRB1", "KLRK1", "CCL5"), order = T, pt.size = 1)
VlnPlot(Healthy_Lymphoid, features = c("CD3D", "CD8A", "NCR3", "GZMB", "KLRB1", "KLRK1", "CCL5"), log=T)
#NAIVE LT
FeaturePlot(Healthy_Lymphoid, features = c("CD3D", "LEF1", "TCF7", "IL7R", "CCR7", "SELL", "LTB"), order = T, pt.size = 1)
VlnPlot(Healthy_Lymphoid, features = c("CD3D", "LEF1", "TCF7", "IL7R", "CCR7", "SELL", "LTB"), log=T)
#LB
FeaturePlot(Healthy_Lymphoid, features = c("MS4A1", "CD19", "FCRLA", "CCR7"), order = T, pt.size = 1)
VlnPlot(Healthy_Lymphoid, features = c("MS4A1", "CD19", "FCRLA", "CCR7"), log=T)
#PLASMA CELLS
FeaturePlot(Healthy_Lymphoid, features = c("TNFRSF17", "JCHAIN", "POU2AF1"), order = T, pt.size = 1)
VlnPlot(Healthy_Lymphoid, features = c("TNFRSF17", "JCHAIN", "POU2AF1"), log=T)
#MEMORY T ?
FeaturePlot(Healthy_Lymphoid, features = c("CD3D", "IL7R", "CD7", "SAMHD1"), order = T, pt.size = 1)
VlnPlot(Healthy_Lymphoid, features = c("CD3D", "IL7R", "CD7", "SAMHD1"), log=T)
#T REGULATORY
FeaturePlot(Healthy_Lymphoid, features = c("CD3D", "CTLA4", "IKZF2", "GATA3", "CCL5", "LGALS3"), order = T, pt.size = 1)
VlnPlot(Healthy_Lymphoid, features = c("CD3D", "CTLA4", "IKZF2", "GATA3", "CCL5", "LGALS3"), log=T)
#T HELPER 17/1/2
FeaturePlot(Healthy_Lymphoid, features = c("CD3D", "IL17RB", "RORA", "TBX21", "ALDOC", "LGALS3"), order = T, pt.size = 1)
VlnPlot(Healthy_Lymphoid, features = c("CD3D", "IL17RB", "RORA", "TBX21", "ALDOC", "LGALS3"), log=T)
dev.off()

###MARKERS Lymphoid
Markers.SubsetLymphoid <- FindAllMarkers(Healthy_Lymphoid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
write.table(Markers.SubsetLymphoid, file = "Outs/GLOBAL_PROCESSING/Markers.Healthy_Lymphoid.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

##### Rename all identities
Healthy_Lymphoid_Renamed <- RenameIdents(object = Healthy_Lymphoid,
                                         "0" = "CD4 T","1" = "NK T","2" = "CD4 T",
                                         "3" = "CD8 T","4" = "CD4 T","5" = "CD4 T",
                                         "6" = "CD4 T reg","7" = "CD4 T","8" = "Plasma cells",
                                         "9" = "CD4 T naive","10" = "B lymphocytes","11" = "Th17-like T",
                                         "12" = "Doublet","13" = "gd T","14" = "NK")
DimPlot(Healthy_Lymphoid_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)
# Remove aberrant cells
Healthy_Lymphoid_Renamed <- subset(Healthy_Lymphoid_Renamed, idents = c("Doublet"), invert = TRUE)
# RE-RUN UMAP AFTER DOUBLET EXCLUSION
Healthy_Lymphoid_Renamed <- RunUMAP(Healthy_Lymphoid_Renamed, reduction = "pca", dims = 1:50, min.dist = 0.35)
# Re-visualize the clusters
DimPlot(Healthy_Lymphoid_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
# store final annotation in metadata
Healthy_Lymphoid_Renamed$celltype <- Idents(Healthy_Lymphoid_Renamed)
saveRDS(Healthy_Lymphoid_Renamed, "Outs/GLOBAL_PROCESSING/Healthy_Lymphoid.RDS")

###FIND MARKERS OF CELL TYPES
Markers.SubsetLymphoidR <- FindAllMarkers(Healthy_Lymphoid_Renamed, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
write.table(Markers.SubsetLymphoidR, file = "Outs/FINAL_MARKERS/Lymphoid_Cells.CSV", sep = "\t", dec = ".", col.names = T, row.names = F)

#### ENDOTHELIAL COMPARTMENT ####
DefaultAssay(seurat_integrated) <- "integrated"
Endo_merged <- subset(seurat_integrated, subset = seurat_clusters %in% c(1,5,6,13,18,23,24,27,30,31,39,41,45))
Endo_split <- SplitObject(Endo_merged, split.by = "orig.ident")
for (i in 1:length(Endo_split)) {
  Endo_split[[i]] <- SCTransform(Endo_split[[i]], vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)
}
integ_features <- SelectIntegrationFeatures(object.list = Endo_split, nfeatures = 3000)
Endo_split <- PrepSCTIntegration(object.list = Endo_split, anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = Endo_split, normalization.method = "SCT", anchor.features = integ_features)
Healthy_Endo <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
Healthy_Endo <- RunPCA(Healthy_Endo, npcs = 50, verbose = T)
ElbowPlot(Healthy_Endo)
Healthy_Endo <- FindNeighbors(Healthy_Endo, dims = 1:12)
clustree(Healthy_Endo, prefix = "integrated_snn_res.")
Healthy_Endo <- FindClusters(Healthy_Endo, resolution = 0.8)
Healthy_Endo <- RunUMAP(Healthy_Endo, reduction = "pca", dims = 1:12, min.dist = 0.35)
DimPlot(Healthy_Endo, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

# cluster 9 are doublets, filter out and re-process after integration
Healthy_Endo <- subset(Healthy_Endo, subset = seurat_clusters %in% c(9), invert=T)
DimPlot(Healthy_Endo, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)

DefaultAssay(Healthy_Endo) <- "SCT"
pdf("Outs/GLOBAL_PROCESSING/Healthy_EndoClusters.pdf")
DimPlot(Healthy_Endo, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)
DimPlot(Healthy_Endo, reduction = "umap", label = TRUE, repel = TRUE, split.by = "orig.ident", pt.size = 1)
DimPlot(Healthy_Endo, reduction = "umap", label = TRUE, repel = TRUE, group.by = "orig.ident", pt.size = 1)
VlnPlot(Healthy_Endo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = F, log = T)
VlnPlot(Healthy_Endo, features = c("percent.RPS", "percent.RPL", "percent.mt"), ncol = 3, pt.size = F, log = T)
DimPlot(Healthy_Endo, reduction = "umap", group.by = "initial_annotation", label=T)+NoLegend()
FeaturePlot(Healthy_Endo, features = c("S.Score","G2M.Score"))
FeaturePlot(Healthy_Endo, features = c("PTPRC", "EPCAM","ACTA2", "COL1A1"), order = T, pt.size = 1)
VlnPlot(Healthy_Endo, features = c("PTPRC", "EPCAM", "ACTA2", "COL1A1"), log=T)
dev.off()

DefaultAssay(Healthy_Endo) <- "SCT"
pdf("Outs/GLOBAL_PROCESSING/Healthy_EndoCellTypes.pdf")
#ENDO GENERAL
FeaturePlot(Healthy_Endo, features = c("PECAM1", "LYVE1", "CDH5", "VWF", "CD34"), order = T, pt.size = 1)
VlnPlot(Healthy_Endo, features = c("PECAM1", "LYVE1", "CDH5", "VWF", "CD34"), log=T)
##LYMPHATIC
FeaturePlot(Healthy_Endo, features = c("VWF", "PROX1", "PDPN", "TBX1", "THY1", "FLT4"), order = T, pt.size = 1)
VlnPlot(Healthy_Endo, features = c("VWF", "PROX1", "PDPN", "TBX1", "THY1", "FLT4"), log=T)
#ARTERY EC
FeaturePlot(Healthy_Endo, features = c("BMX", "MGP", "GJA5", "EFNB2"), order = T, pt.size = 1)
VlnPlot(Healthy_Endo, features = c("BMX", "MGP", "GJA5", "EFNB2"), log=T)
##VENOUS EC
FeaturePlot(Healthy_Endo, features = c("SELP", "SELE", "VCAM1", "TIMP1", "ACKR1"), order = T, pt.size = 1)
VlnPlot(Healthy_Endo, features = c("SELP", "SELE", "VCAM1", "TIMP1", "ACKR1"), log=T)
#CAPILLARY EC
FeaturePlot(Healthy_Endo, features = c("CA4", "SPARC", "SGK1", "EMP2"), order = T, pt.size = 1)
VlnPlot(Healthy_Endo, features = c("CA4", "SPARC", "SGK1", "EMP2"), log=T)
##AEROCYTE
FeaturePlot(Healthy_Endo, features = c("VWF", "EMCN", "EDNRB", "TBX2"), order = T, pt.size = 1)
VlnPlot(Healthy_Endo, features = c("VWF", "EMCN", "EDNRB", "TBX2"), log=T)
##GENERAL CAPILLARY
FeaturePlot(Healthy_Endo, features = c("EDN1", "GPIHBP1", "CD36"), order = T, pt.size = 1)
VlnPlot(Healthy_Endo, features = c("EDN1", "GPIHBP1", "CD36"), log=T)
dev.off()

#####MARKERS OF NEW CLUSTERS 
Markers.Endo <- FindAllMarkers(Healthy_Endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
write.table(Markers.Endo, file = "Outs/GLOBAL_PROCESSING/Markers.Healthy_Endo.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

#### Rename all identities
Healthy_Endo_Renamed <- RenameIdents(object = Healthy_Endo,
                                     "0" = "General capillary","1" = "Aerocytes",
                                     "2" = "General capillary","3" = "General capillary",
                                     "4" = "General capillary","5" = "Aerocytes",
                                     "6" = "Venous","7" = "Arterial", "8" = "Aerocytes",
                                     "9" = "General capillary","10" = "General capillary",
                                     "11" = "Lymphatic")
DimPlot(Healthy_Endo_Renamed, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2)
# store final annotation in metadata
Healthy_Endo_Renamed$celltype <- Idents(Healthy_Endo_Renamed)
saveRDS(Healthy_Endo_Renamed, "Outs/GLOBAL_PROCESSING/Healthy_Endo.RDS")

#####FIND MARKERS OF CELL TYPES
Markers.EndoR <- FindAllMarkers(Healthy_Endo_Renamed, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
write.table(Markers.EndoR, file = "Outs/FINAL_MARKERS/Endo_Cells.CSV", sep = "\t", dec = ".", col.names = T, row.names = F)

#### FUSION OF ALL COMPARTMENTS ####
Healthy_Epith_Renamed <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_Epith.RDS")
Healthy_Myeloid_Renamed <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_Myeloid.RDS")
Healthy_Lymphoid_Renamed <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_Lymphoid.RDS")
Healthy_Endo_Renamed <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_Endo.RDS")
Healthy_Mesench_integrated <- readRDS("Outs/GLOBAL_PROCESSING/Mesench_integrated.RDS")

Healthy_integrated_merged <- merge(Healthy_Epith_Renamed, y = c(Healthy_Endo_Renamed, Healthy_Myeloid_Renamed, Healthy_Lymphoid_Renamed, Healthy_Mesench_integrated),
                                   add.cell.ids = c("Epithelial", "Endothelial", "Immune", "Immune", "Mesenchymal"), project = "HEALTHY")
Healthy_integrated_split <- SplitObject(Healthy_integrated_merged, split.by = "orig.ident")
for (i in 1:length(Healthy_integrated_split)) {
  Healthy_integrated_split[[i]] <- SCTransform(Healthy_integrated_split[[i]], vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)
}
integ_features <- SelectIntegrationFeatures(object.list = Healthy_integrated_split, nfeatures = 3000)
Healthy_integrated_split <- PrepSCTIntegration(object.list = Healthy_integrated_split, anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = Healthy_integrated_split, normalization.method = "SCT", anchor.features = integ_features)
Healthy_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
Healthy_integrated <- RunPCA(Healthy_integrated, npcs = 50, verbose = T)
Healthy_integrated <- RunUMAP(Healthy_integrated, reduction = "pca", dims = 1:50, min.dist = 0.25)
saveRDS(Healthy_integrated, "Outs/GLOBAL_PROCESSING/Healthy_integrated.RDS")
DimPlot(Healthy_integrated, reduction = "umap", label = TRUE, shuffle=T, 
        group.by="celltype", repel = TRUE, pt.size = 0.05)+ NoLegend()

# Store compartment information in metadata
head(colnames(Healthy_integrated))
compartment <- sapply(X = strsplit(colnames(Healthy_integrated), split = "_"), FUN = "[", 1)
Healthy_integrated$compartment <- compartment
DimPlot(Healthy_integrated, reduction = "umap", label = TRUE, shuffle=T, 
        group.by="compartment", repel = TRUE, pt.size = 0.05)+ NoLegend()
saveRDS(Healthy_integrated, "Outs/GLOBAL_PROCESSING/Healthy_integrated.RDS")

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(Healthy_integrated, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)
View(n_cells)
write.csv2(n_cells, file="Outs/CELL_DISTRIBUTION/CellRepartition.csv")

# Extract compartment and sample information 
n_cells_cpt <- FetchData(Healthy_integrated, vars = c("compartment", "orig.ident")) %>%
  dplyr::count(compartment, orig.ident) %>%
  tidyr::spread(compartment, n)
View(n_cells_cpt)
write.csv2(n_cells_cpt, file="Outs/CELL_DISTRIBUTION/CellRepartition_bycompartment.csv")

##GENERAL MARKERS
Markers.AllCells <- FindAllMarkers(Healthy_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
write.table(Markers.AllCells, file = "Outs/FINAL_MARKERS/All_Cells.CSV", sep = "\t", dec = ".", col.names = T, row.names = F)

######################################################################################################
####                                   END OF DATA PROCESSING                                     ####
######################################################################################################