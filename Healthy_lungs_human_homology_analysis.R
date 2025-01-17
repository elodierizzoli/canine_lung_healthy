#############################################################################################################
#### A single cell RNA sequencing atlas of the healthy canine lung: a foundation for comparative studies ####
####                                    CANINE/HUMAN HOMOLOGY ANALYSIS                                   ####
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
library(orthogene)

### GENES FOR CELL CYCLE INVESTIGATIONS ###
ccgenes.s.genes<-toupper(c("Cdc45","Pcna","Chaf1b","Pola1","Tyms",
                           "Gmnn","Mrpl36","Mcm5","Slbp","Rrm2",
                           "Rrm1","Uhrf1","Rad51","Rad51ap1","E2f8",
                           "Mcm4","Msh2","Wdr76","Rfc2","Ccne2","Hells",
                           "Clspn","Tipin","Fen1","Dscc1","Mcm7",
                           "Casp8ap2","Ubr7","Cdc6","Nasp","Polr1b",
                           "Prim1","Usp1","Blm","Mcm6","Gins2","Cenpu",
                           "Cdca7","Ung","Dtl","Exo1"))
ccgenes.g2m.genes<-toupper(c("Ckap2","Tpx2","Ndc80","Tacc3","Cdca2",
                             "Rangap1","Cenpa","Ube2c","Anln","Cks2",
                             "Cdk1","Gtse1","Pimreg","Nusap1","G2e3",
                             "Ncapd2","Aurka","Cdca3","Aurkb","Dlgap5",
                             "Ckap5","Ccnb2","Kif20b","Kif11","Ttk",
                             "Cdca8","Kif23","Tubb4b","Cenpe","Cdc20",
                             "Cdc25c","Kif2c","Cbx5","Top2a","Ctcf",
                             "Bub1","Ckap2l","Mki67","Hmmr","Hmgb2",
                             "Tmpo","Psrc1","Gas2l3","Smc4","Ect2",
                             "Anp32e","Cks1b","Birc5","Nuf2","Hjurp",
                             "Nek2","Cenpf","Lbr"))

#### DEFINING PATH
# set working directory

# 1 # LOAD HUMAN REFERENCE DATASET ####
# human reference
lung.ref <- readRDS("Outs/HUMAN_INTEGRATION/humanlung.RDS")

# transformation into v4 object
counts <- GetAssayData(lung.ref, slot = "counts")
data <- GetAssayData(lung.ref, slot = "data")
meta_data <- lung.ref@meta.data
reductions <- lung.ref@reductions
seurat_v4 <- CreateSeuratObject(counts = counts, meta.data = meta_data)
DefaultAssay(seurat_v4) <- "RNA"
seurat_v4 <- SetAssayData(seurat_v4, slot = "data", new.data = data)
seurat_v4@reductions <- reductions
saveRDS(seurat_v4, file = "Outs/HUMAN_INTEGRATION/seurat_v4_object.rds")

# interspecies gene mapping (DOI: 10.18129/B9.bioc.orthogene)
convert <- human_lung@assays$RNA@counts
convert <- orthogene::convert_orthologs(gene_df = cnts2,
                                      gene_input = "rownames", 
                                      gene_output = "rownames", 
                                      input_species = "human",
                                      output_species = "dog",
                                      non121_strategy = "drop_both_species") 
rownames(convert) <- unname(rownames(convert))

seu.obj.human <- CreateSeuratObject(convert, project = "humanConvert", assay = "RNA",
                                    min.cells = 0, min.features = 0, names.field = 1,
                                    names.delim = "_", meta.data = human_lung@meta.data)

# metadata "percent.mt"
seu.obj.human <- PercentageFeatureSet(seu.obj.human, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(seu.obj.human, features = c("percent.mt"),
        ncol = 3, pt.size = F, log = F)
# dataset does not contain mitochondrial reads
# vars.to.regress will be only c("S.Score", "G2M.Score")

# decreasing size of dataset
subset_hu <- subset(seu.obj.human, subset = tissue == "lung parenchyma")
table(subset_hu$tissue)

# cleaning metadata
colnames(subset_hu@meta.data)
subset_hu@meta.data <- subset_hu@meta.data[, !grepl("ontology", colnames(subset_hu@meta.data))]
subset_hu@meta.data <- subset_hu@meta.data[, !grepl("leiden", colnames(subset_hu@meta.data))]
subset_hu@meta.data <- subset_hu@meta.data[, !(colnames(subset_hu@meta.data) %in% c("tissue", "tissue_type", 
                                                                                    "observation_joinid",
                                                                                    "fresh_or_frozen", "percent.mt",
                                                                                    "development_stage",
                                                                                    "cause_of_death",
                                                                                    "donor_id", "BMI",
                                                                                    "age_or_mean_of_age_range",
                                                                                    "age_range","scanvi_label",
                                                                                    "size_factors","reference_genome",
                                                                                    "sex","self_reported_ethnicity",
                                                                                    "subject_type", "tissue_sampling_method", 
                                                                                    "disease", "suspension_type",
                                                                                    "is_primary_data", "anatomical_region_ccf_score",
                                                                                    "mixed_ancestry"
                                                                                    ))]
print(object.size(subset_hu), units = "MB")

# downsampling
subset_cells <- sample(Cells(subset_hu), size = 50000) 
subset_hu50K <- subset(subset_hu, cells = subset_cells)
table(subset_hu50K$ann_finest_level)
print(object.size(subset_hu50K), units = "MB")
saveRDS(subset_hu50K, file = "Outs/HUMAN_INTEGRATION/subset_hu50K.rds")

# 2 # LOAD CANINE DATASET ####
Healthy_integrated <- readRDS("Outs/GLOBAL_PROCESSING/Healthy_integrated.RDS")

# add useful metadata
Healthy_integrated$ann_finest_level <- as.factor(Healthy_integrated$celltype)
Healthy_integrated$organism <- as.factor("Canis lupus familiaris")

# clean environment
rm(subset_hu)
gc()

# 3 # INTEGRATION ####
options(future.globals.maxSize = 7000 * 1024^2)
atlas_merged <- merge(Healthy_integrated, y = c(subset_hu50K),
                                   add.cell.ids = c("Canine", "Human"), project = "INTEG")
atlas_split <- SplitObject(atlas_merged, split.by = "orig.ident")
for (i in 1:length(atlas_split)) {
  atlas_split[[i]] <- CellCycleScoring(atlas_split[[i]], s.features = ccgenes.s.genes, g2m.features = ccgenes.g2m.genes)
}
for (i in 1:length(atlas_split)) {
  atlas_split[[i]] <- SCTransform(atlas_split[[i]], vst.flavor = "v2", vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
}
integ_features <- SelectIntegrationFeatures(atlas_split, nfeatures = 3000)
atlas_split <- PrepSCTIntegration(object.list = atlas_split, anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = atlas_split, normalization.method = "SCT", anchor.features = integ_features)
integ.atlas <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
ElbowPlot(integ.atlas)
integ.atlas <- RunPCA(integ.atlas, npcs = 10, verbose = T)
integ.atlas <- FindNeighbors(integ.atlas, reduction = "pca", dims = 1:10)
clustree(integ.atlas, prefix = "integrated_snn_res.")
integ.atlas <- FindClusters(integ.atlas, resolution = 0.5)
integ.atlas <- RunUMAP(integ.atlas, reduction = "pca", dims = 1:10, min.dist = 0.25)
saveRDS(integ.atlas, "Outs/HUMAN_INTEGRATION/integ.atlas.RDS")

# 4 # VISUALIZATION ####
# dim plots
integ.atlas <- readRDS("Outs/HUMAN_INTEGRATION/integ.atlas.RDS")
DimPlot(integ.atlas, reduction = "umap", label = TRUE, group.by="seurat_clusters", 
        split.by="organism", pt.size = 0.05)

DimPlot(integ.atlas, reduction = "umap", group.by="ann_finest_level", 
        split.by="organism", pt.size = 0.1, label=T, repel=T, label.size=3) + 
  NoLegend()

DimPlot(integ.atlas, reduction = "umap", group.by="organism", 
       pt.size = 0.05, repel=T, label=T, label.size = 3)

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

#stash new metadata with cell type name prefixed with the specices
integ.atlas_R$type <- ifelse(grepl("Canis lupus familiaris", integ.atlas_R$organism), paste0("can_",integ.atlas_R$ann_finest_level), paste0("hu_",integ.atlas_R$ann_finest_level))    
integ.atlas_R$type <- as.factor(integ.atlas_R$type)
levels(integ.atlas_R$type) <- gsub(" ","_",levels(integ.atlas_R$type))
levels(integ.atlas_R$type)

DimPlot(integ.atlas_R, reduction = "umap", group.by="type", 
        split.by="organism", pt.size = 0.1, label=T, repel=T, label.size=3) + 
  NoLegend()

# 5 # ANALYSIS OF INTEGRATED DATA ####
# hierarchical clustering
# Method adopted from Ammons et al., 2023 'A single-cell RNA sequencing atlas of circulating leukocytes from healthy and osteosarcoma affected dogs' (DOI 10.3389/fimmu.2023.1162700)
#https://github.com/dyammons/canine_osteosarcoma_atlas/
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

png("integ_type.png", width=4000, height=4000, res=400)
par(mfcol=c(1,1))
hc <- hclust(as.dist(M),method="complete")
p <- plot.phylo(as.phylo(hc), type = "fan", cex = 0.8,no.margin = FALSE,tip.col = "black")
dev.off()

######################################################################################################
####                                 END OF HOMOLOGY ANALYSIS                                     ####
######################################################################################################