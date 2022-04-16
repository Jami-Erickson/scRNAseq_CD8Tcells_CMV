## Script for integrating single cell data in Seurat v4
## CD3+ sorted T population as well as Lin-HLADR+ population from various tissues
## JRE, 2020-01-28

## Load required packages, including MAST
library(Seurat)
library(tidyverse)
library(MAST)
library(RColorBrewer)
library(doMC)
library(data.table)
library(circlize)
library(reticulate)
library(future)
library(ggplot2)
library(sctransform)
library(dplyr)
library(scater)
library(loomR)
library(scran)
library(SeuratWrappers)
library(batchelor)
library(devtools)
library(harmony)
library(liger)



## Enable parallelization in Seurat and MAST
plan("multisession", workers = 7)

options(mc.cores = 7)
knitr::opts_chunk$set(message = F, error = F, warning = F, cache = F, fig.width = 2, fig.height = 6, auto.dep = T)

## Deal with the RAM issue
options(future.globals.maxSize = 12000 * 1024^2)

## Set working directory
setwd('~/Desktop/R/CMV/')

######## Process all CD8 data sets
## Load datasets
CD8_Patient1.data <- Read10X(data.dir = 'Patient1_CD8_Aggr')

CD8_Patient2.data <- Read10X(data.dir = 'Patient2_CD8_Aggr')

CD8_Patient3.data <- Read10X(data.dir = 'Patient3_CD8_Aggr')

CD8_Patient4.data <- Read10X(data.dir = 'Patient4_CD8_Aggr')


## Create all Seurat objects
#Patient1

Genotype1 <- read.csv("Patient1DvR.csv", header = T)
Genotype1 <- as.matrix(Genotype4)
head(Genotype1)
tail(Genotype1)

AnnotationsPatient1 <- read.csv("Patient1annotations.csv", header = T)
AnnotationsPatient1 <- as.matrix(AnnotationsPatient1)
head(AnnotationsPatient1)
tail(AnnotationsPatient1)

GenoAnnoPatient1 <- merge(Genotype4, AnnotationsPatient1, by.x = 'barcode', by.y = 'barcode', all.x = T, all.y = T, sort = F)
head(GenoAnnoPatient1)
tail(GenoAnnoPatient1)
head(CD8_Patient1.data[1,])

cellids10x1 <- as.data.frame(CD8_Patient1.data[1,])
head(cellids10x1)
setDT(cellids10x1, keep.rownames = T)[]
colnames(cellids10x1) <- c('barcodes', 'NA')

GenoAnnoPatient1b <- merge(cellids10x1, GenoAnnoPatient1, by.x = 'barcodes', by.y = 'barcode', all.x = T, all.y = F, sort = F)
head(GenoAnnoPatient1b)
tail(GenoAnnoPatient1b)
GenoAnnoPatient1b <- as.matrix(GenoAnnoPatient1b)

cellids1 <- c(GenoAnnoPatient1b[,1])
genotype1 <- c(GenoAnnoPatient1b[,3])
annotations1 <- c(GenoAnnoPatient1b[,4])
clonotype_id1 <- c(GenoAnnoPatient1b[,5])
clustermap_row_number1 <- c(GenoAnnoPatient1b[,6])

write.csv(GenoAnnoPatient1b, "GenoAnnoPatient1b.csv")


CD8_Patient1_Seurat <- CreateSeuratObject(counts = CD8_Patient1.data, 
                                          assay = 'RNA',
                                          project = 'CMV_CD8_Patient1')
CD8_Patient1_Seurat # 33694 genes across 15600 samples
CD8_Patient1_Seurat@meta.data$CellID <- cellids1
CD8_Patient1_Seurat@meta.data$Genotype <- genotype1
CD8_Patient1_Seurat@meta.data$Annotations <- annotations1
CD8_Patient1_Seurat@meta.data$ClonotypeID <- clonotype_id1
CD8_Patient1_Seurat@meta.data$ClustermapRowNumber <- clustermap_row_number1
CD8_Patient1_Seurat@meta.data$SampleID <- sapply(rownames(CD8_Patient1_Seurat@meta.data), function(x) str_split(x, '-')[[1]][2])
CD8_Patient1_Seurat@meta.data$Tissue <- plyr::mapvalues(x = CD8_Patient1_Seurat@meta.data$SampleID,
                                                        from = c('1', '2', '3'),
                                                        to = c('CD8_30_Patient1', 'CD8_60_Patient1', 'CD8_90_Patient1'))
CD8_Patient1_Seurat@meta.data$Time <- plyr::mapvalues(x = CD8_Patient1_Seurat@meta.data$SampleID,
                                                      from = c('1', '2', '3'),
                                                      to = c('30', '60', '90'))

CD8_Patient1_Seurat <- RenameCells(CD8_Patient1_Seurat, add.cell.id = "Patient1")
CD8_Patient1_Seurat@meta.data$sample <- "Patient1"
CD8_Patient1_Seurat[["percent.mt"]] <- PercentageFeatureSet(CD8_Patient1_Seurat, pattern = "^MT-")
CD8_Patient1_Seurat[["percent.ribo"]] <- PercentageFeatureSet(CD8_Patient1_Seurat, pattern = "^RP")
VlnPlot(CD8_Patient1_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
CD8_Patient1_Seurat_2 <- subset(CD8_Patient1_Seurat, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 17500)

saveRDS(CD8_Patient1_Seurat_2, "CD8_Patient1_Seurat_2.RDS")

#Patient2

Genotype2 <- read.csv("Patient2DP2.csv", header = T)
Genotype2 <- as.matrix(Genotype2)

AnnotationsPatient2 <- read.csv("Patient2annotations.csv", header = T)
AnnotationsPatient2 <- as.matrix(AnnotationsPatient2)

GenoAnnoPatient2 <- merge(Genotype2, AnnotationsPatient2, by.x = 'barcode', by.y = 'barcode', all.x = T, all.y = T)
head(GenoAnnoPatient2)
head(CD8_Patient2.data[1,])

cellids10x2 <- as.data.frame(CD8_Patient2.data[1,])
head(cellids10x2)
library(data.table)
setDT(cellids10x2, keep.rownames = T)[]
colnames(cellids10x2) <- c('barcodes', 'NA')

GenoAnnoPatient2b <- merge(cellids10x2, GenoAnnoPatient2, by.x = 'barcodes', by.y = 'barcode', all.x = T, all.y = F, sort = F)
head(GenoAnnoPatient2b)
GenoAnnoPatient2b <- as.matrix(GenoAnnoPatient2b)

cellids2 <- c(GenoAnnoPatient2b[,1]) 
genotype2 <- c(GenoAnnoPatient2b[,3])
annotations2 <- c(GenoAnnoPatient2b[,4])
clonotype_id2 <- c(GenoAnnoPatient2b[,5])
clustermap_row_number2 <- c(GenoAnnoPatient2b[,6])

write.csv(GenoAnnoPatient2b, "GenoAnnoPatient2b.csv")


#- Setup seurat object class
CD8_Patient2_Seurat <- CreateSeuratObject(counts = CD8_Patient2.data, 
                                       assay = 'RNA',
                                       project = 'CMV_CD8_Patient2')
CD8_Patient2_Seurat # 33694 genes across 15600 samples


#- Sample info
CD8_Patient2_Seurat@meta.data$CellID <- cellids2
CD8_Patient2_Seurat@meta.data$Genotype <- genotype2
CD8_Patient2_Seurat@meta.data$Annotations <- annotations2
CD8_Patient2_Seurat@meta.data$ClonotypeID <- clonotype_id2
CD8_Patient2_Seurat@meta.data$ClustermapRowNumber <- clustermap_row_number
CD8_Patient2_Seurat@meta.data$SampleID <- sapply(rownames(CD8_Patient2_Seurat@meta.data), function(x) str_split(x, '-')[[1]][2])
CD8_Patient2_Seurat@meta.data$Tissue <- plyr::mapvalues(x = CD8_Patient2_Seurat@meta.data$SampleID,
                                                     from = c('1', '2', '3'),
                                                     to = c('CD8_30_Patient2', 'CD8_60_Patient2', 'CD8_90_Patient2'))
CD8_Patient2_Seurat@meta.data$Time <- plyr::mapvalues(x = CD8_Patient2_Seurat@meta.data$SampleID,
                                                   from = c('1', '2', '3'),
                                                   to = c('30', '60', '90'))

table(CD8_Patient2_Seurat@meta.data$Time)

CD8_Patient2_Seurat <- RenameCells(CD8_Patient2_Seurat, add.cell.id = "Patient2")
CD8_Patient2_Seurat@meta.data$sample <- "Patient2"
CD8_Patient2_Seurat[["percent.mt"]] <- PercentageFeatureSet(CD8_Patient2_Seurat, pattern = "^MT-")
CD8_Patient2_Seurat[["percent.ribo"]] <- PercentageFeatureSet(CD8_Patient2_Seurat, pattern = "^RP")
VlnPlot(CD8_Patient2_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)

plot1 <- FeatureScatter(object = CD8_Patient2_Seurat, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
plot2 <- FeatureScatter(object = CD8_Patient2_Seurat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
CombinePlots(plots = list(plot1, plot2))

##-- Filtering (removing the outlier cells with percent.mito being higher than 10%
CD8_Patient2_Seurat_2 <- subset(CD8_Patient2_Seurat, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 2600 & nCount_RNA < 10000)



# -Inf and Inf should be used if you don't want a lower or upper threshold
CD8_Patient2_Seurat_2 ##-- 33694 features across 17903 samples.

#Re-plot after filtering
plot1 <- FeatureScatter(object = CD8_Patient2_Seurat_2, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
plot2 <- FeatureScatter(object = CD8_Patient2_Seurat_2, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
CombinePlots(plots = list(plot1, plot2))

saveRDS(CD8_Patient2_Seurat_2, "CD8_Patient2_Seurat_2.RDS")

#Patient3

Genotype3 <- read.csv("Patient3DvR_time.csv", header = T)
Genotype3 <- as.matrix(Genotype3)

AnnotationsPatient3 <- read.csv("Patient3annotations.csv", header = T)
AnnotationsPatient3 <- as.matrix(AnnotationsPatient3)

GenoAnnoPatient3 <- merge(Genotype3, AnnotationsPatient3, by.x = 'barcode', by.y = 'barcode', all.x = T, all.y = T, sort = F)
head(GenoAnnoPatient3)
head(CD8_Patient3.data[1,])

cellids10x3 <- as.data.frame(CD8_Patient3.data[1,])
head(cellids10x3)
library(data.table)
setDT(cellids10x3, keep.rownames = T)[]
colnames(cellids10x3) <- c('barcodes', 'NA')

GenoAnnoPatient3b <- merge(cellids10x3, GenoAnnoPatient3, by.x = 'barcodes', by.y = 'barcode', all.x = T, all.y = F, sort = F)
head(GenoAnnoPatient3b)
GenoAnnoPatient3b <- as.matrix(GenoAnnoPatient3b)

cellids3 <- c(GenoAnnoPatient3b[,1])
genotype3 <- c(GenoAnnoPatient3b[,3])
annotations3 <- c(GenoAnnoPatient3b[,4])
clonotype_id3 <- c(GenoAnnoPatient3b[,5])
clustermap_row_number3 <- c(GenoAnnoPatient3b[,6])

write.csv(GenoAnnoPatient3b, "GenoAnnoPatient3b.csv")

#- Setup seurat object class
CD8_Patient3_Seurat <- CreateSeuratObject(counts = CD8_Patient3.data, 
                                       assay = 'RNA',
                                       project = 'CMV_CD8_Patient3')
CD8_Patient3_Seurat # 33694 genes across 15600 samples


#- Sample info
CD8_Patient3_Seurat@meta.data$CellID <- cellids3
CD8_Patient3_Seurat@meta.data$Genotype <- genotype3
CD8_Patient3_Seurat@meta.data$Annotations <- annotations3
CD8_Patient3_Seurat@meta.data$ClonotypeID <- clonotype_id3
CD8_Patient3_Seurat@meta.data$ClustermapRowNumber <- clustermap_row_number3
CD8_Patient3_Seurat@meta.data$SampleID <- sapply(rownames(CD8_Patient3_Seurat@meta.data), function(x) str_split(x, '-')[[1]][2])
CD8_Patient3_Seurat@meta.data$Tissue <- plyr::mapvalues(x = CD8_Patient3_Seurat@meta.data$SampleID,
                                                     from = c('1', '2', '3'),
                                                     to = c('CD8_30_Patient3', 'CD8_60_Patient3', 'CD8_90_Patient3'))
CD8_Patient3_Seurat@meta.data$Time <- plyr::mapvalues(x = CD8_Patient3_Seurat@meta.data$SampleID,
                                                   from = c('1', '2', '3'),
                                                   to = c('30', '60', '90'))

CD8_Patient3_Seurat <- RenameCells(CD8_Patient3_Seurat, add.cell.id = "Patient3")
CD8_Patient3_Seurat@meta.data$sample <- "Patient3"
CD8_Patient3_Seurat[["percent.mt"]] <- PercentageFeatureSet(CD8_Patient3_Seurat, pattern = "^MT-")
CD8_Patient3_Seurat[["percent.ribo"]] <- PercentageFeatureSet(CD8_Patient3_Seurat, pattern = "^RP")
VlnPlot(CD8_Patient3_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
CD8_Patient3_Seurat_2 <- subset(CD8_Patient3_Seurat, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 4200 & nCount_RNA < 21000)

table(CD8_Patient3_Seurat@meta.data$Time)

saveRDS(CD8_Patient3_Seurat_2, "CD8_Patient3_Seurat_2.RDS")

#Patient4
Genotype4 <- read.csv("Patient4DP2.csv", header = T)
Genotype4 <- as.matrix(Genotype3)

AnnotationsPatient4 <- read.csv("Patient4annotations.csv", header = T)
AnnotationsPatient4 <- as.matrix(AnnotationsPatient4)

GenoAnnoPatient4 <- merge(Genotype4, AnnotationsPatient4, by.x = 'barcode', by.y = 'barcode', all.x = T, all.y = T, sort = F)
head(GenoAnnoPatient4)
head(CD8_Patient4.data[1,])

cellids10x4 <- as.data.frame(CD8_Patient4.data[1,])
head(cellids10x4)
library(data.table)
setDT(cellids10x4, keep.rownames = T)[]
colnames(cellids10x4) <- c('barcodes', 'NA')

GenoAnnoPatient4b <- merge(cellids10x4, GenoAnnoPatient4, by.x = 'barcodes', by.y = 'barcode', all.x = T, all.y = F, sort = F)
head(GenoAnnoPatient4b)
GenoAnnoPatient4b <- as.matrix(GenoAnnoPatient4b)

cellids4 <- c(GenoAnnoPatient4b[,1])
genotype4 <- c(GenoAnnoPatient4b[,3])
annotations4 <- c(GenoAnnoPatient4b[,4])
clonotype_id4 <- c(GenoAnnoPatient4b[,5])
clustermap_row_number4 <- c(GenoAnnoPatient4b[,6])

write.csv(GenoAnnoPatient4b, "GenoAnnoPatient4b.csv")

#- Setup seurat object class
CD8_Patient4_Seurat <- CreateSeuratObject(counts = CD8_Patient4.data, 
                                       assay = 'RNA',
                                       project = 'CMV_CD8_Patient4')
CD8_Patient4_Seurat # 33694 genes across 15600 samples


#- Sample info
CD8_Patient4_Seurat@meta.data$CellID <- cellids4
CD8_Patient4_Seurat@meta.data$Genotype <- genotype4
CD8_Patient4_Seurat@meta.data$Annotations <- annotations4
CD8_Patient4_Seurat@meta.data$ClonotypeID <- clonotype_id4
CD8_Patient4_Seurat@meta.data$ClustermapRowNumber <- clustermap_row_number4
CD8_Patient4_Seurat@meta.data$SampleID <- sapply(rownames(CD8_Patient4_Seurat@meta.data), function(x) str_split(x, '-')[[1]][2])
CD8_Patient4_Seurat@meta.data$Tissue <- plyr::mapvalues(x = CD8_Patient4_Seurat@meta.data$SampleID,
                                                     from = c('1', '2', '3'),
                                                     to = c('CD8_30_Patient4', 'CD8_60_Patient4', 'CD8_90_Patient4'))
CD8_Patient4_Seurat@meta.data$Time <- plyr::mapvalues(x = CD8_Patient4_Seurat@meta.data$SampleID,
                                                   from = c('1', '2', '3'),
                                                   to = c('30', '60', '90'))

CD8_Patient4_Seurat <- RenameCells(CD8_Patient4_Seurat, add.cell.id = "Patient4")
CD8_Patient4_Seurat@meta.data$sample <- "Patient4"
CD8_Patient4_Seurat[["percent.mt"]] <- PercentageFeatureSet(CD8_Patient4_Seurat, pattern = "^MT-")
CD8_Patient4_Seurat[["percent.ribo"]] <- PercentageFeatureSet(CD8_Patient4_Seurat, pattern = "^RP")
VlnPlot(CD8_Patient4_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
CD8_Patient4_Seurat_2 <- subset(CD8_Patient4_Seurat, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 3200 & nCount_RNA < 12500)

table(CD8_Patient4_Seurat_2@meta.data$Time)

saveRDS(CD8_Patient4_Seurat_2, "CD8_Patient4_Seurat_2.RDS")



CD8_Patient1_Seurat_2 <- readRDS("CD8_Patient1_Seurat_2.RDS")
CD8_Patient2_Seurat_2 <- readRDS("CD8_Patient2_Seurat_2.RDS")
CD8_Patient3_Seurat_2 <- readRDS("CD8_Patient3_Seurat_2.RDS")
CD8_Patient4_Seurat_2 <- readRDS("CD8_Patient4_Seurat_2.RDS")

## Generate merged CD8 object
CD8 <- merge(CD8_Patient2_Seurat_2, c(CD8_Patient3_Seurat_2, CD8_Patient4_Seurat_2, CD8_Patient1_Seurat_2))
CD8 #33694 features across 65100 samples within 1 assay  

table(CD8@meta.data$Tissue)
table(CD8_Patient1_Seurat_2@meta.data$Tissue)
table(CD8_Patient2_Seurat_2@meta.data$Tissue)
table(CD8_Patient3_Seurat_2@meta.data$Tissue)
table(CD8_Patient4_Seurat_2@meta.data$Tissue)



############# Harmony workflow
## Generate merged object with all data for HARMONY integration
CD8 <- NormalizeData(CD8, verbose = T)
CD8 <- FindVariableFeatures(CD8, selection.method = "vst", nfeatures = 10000, verbose = T)
CD8 <- ScaleData(CD8, verbose = T)

## take sex-linked genes out of var.features

var.genes <- as.vector(CD8@assays$RNA@var.features)
genes.filtered.xist <- grep(pattern = "XIST", var.genes, value = TRUE, invert = TRUE)
vargenes.filtered <- as.vector(grep(pattern = "RPS4Y1", genes.filtered.xist, value = TRUE, invert = TRUE))
vargenes.filtered2 <- as.vector(grep(pattern = "^TTTY", vargenes.filtered, value = TRUE, invert = TRUE))
vargenes.filtered3 <- as.vector(grep(pattern = "UTY", vargenes.filtered2, value = TRUE, invert = TRUE))
vargenes.filtered4 <- as.vector(grep(pattern = "ZFY", vargenes.filtered3, value = TRUE, invert = TRUE))
vargenes.filtered5 <- as.vector(grep(pattern = "EIF1AY", vargenes.filtered4, value = TRUE, invert = TRUE))




write.csv(vargenes.filtered, "vargenes.csv")

CD8 <- RunPCA(CD8, pc.genes = vargenes.filtered5, npcs = 50, verbose = T)

CD8 <- RunHarmony(CD8, group.by.vars = "sample",
                  lambda = 1, theta = 2, plot_convergence = T, verbose = T, max.iter.harmony = 20)
CD8
table(CD8$sample)
table(CD8$Time)
table(CD8@meta.data$Genotype)

ElbowPlot(object = CD8,
          ndims = 50)

## Perform UMAP calculation and graph-based clustering #switch "pca" to "harmony" if running harmony. 
CD8 <- RunUMAP(CD8, reduction = "harmony", dims = 1:40, seed.use = 34)
CD8 <- FindNeighbors(CD8, reduction = "harmony", dims = 1:40, verbose = T)
CD8 <- FindClusters(CD8, resolution = 0.8, random.seed = 34, verbose = T)



## Plotting
color.vector <- c(brewer.pal(8, "Set1"), brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))
DimPlot(CD8, reduction = "umap", pt.size = 0.1, label = F, group.by="sample", cols = color.vector)
DimPlot(CD8, reduction = "umap", pt.size = 0.1, label = F, group.by = "Tissue", cols = color.vector)
DimPlot(CD8, reduction = "umap", pt.size = 0.1, label = T, group.by = "seurat_clusters", cols = color.vector)
DimPlot(CD8, reduction = "umap", pt.size = 0.1, label = F, group.by="orig.ident", cols = brewer.pal(4, "Set1"))
DimPlot(CD8, reduction = "umap", pt.size = 0.1, label = F, group.by="Genotype", cols = brewer.pal(3, "Set1"))
DimPlot(CD8, reduction = "umap", pt.size = 0.1, label = F, group.by="sample", cols = color.vector)



FeaturePlot(CD8, features = c("CD14", "ITGAX", "IL3RA", "CD79A", "IGHM", "CD3E", "FOXP3", "CD4", "CD8B"), 
            pt.size = 0.1, ncol = 3)

CD8 <- SetIdent(CD8, value = "RNA_snn_res.0.8")
CD8.clean <- subset(CD8, idents = "22", invert = T)
plot <- DimPlot(CD8.clean, reduction = "umap", cols = color.vector, pt.size = 0.8, label = T)
HoverLocator(plot = plot, information = FetchData(CD8.clean, vars = c("CellID", "Tissue", 'Genotype')))

all.genes <- rownames(CD8.clean2@assays$RNA@data)
genes.filtered.mt <- grep(pattern = "^MT-", all.genes, value = TRUE, invert = TRUE)
genes.filtered <- as.vector(grep(pattern = "^RP", genes.filtered.mt, value = TRUE, invert = TRUE))
length(all.genes)
length(genes.filtered)



CD8.clean <- SetIdent(CD8.clean, value = "CellID")
CD8.clean2 <- subset(CD8.clean, ident = c('GAAATGACATTATCTC-1'), invert = T)

DimPlot(CD8.clean2, reduction = "umap", pt.size = 0.1, label = F, group.by="Genotype", cols = brewer.pal(3, "Set1"))


CD8.clean2 <- SetIdent(CD8.clean2, value = "Genotype")
CD8.clean3 <- subset(CD8.clean2, ident = 'Unknown', invert = T)

GenoColor <- c("#C7C2E3", "#276419")
DimPlot(CD8.clean3, reduction = "umap", pt.size = 0.1, label = F, group.by="Genotype", cols = GenoColor)
SampleColor <- c("#e7298a", "#1b9e77", "#d95f02", "#7570b3")
DimPlot(CD8.clean3, reduction = "umap", pt.size = 0.1, label = F, group.by="sample", cols = SampleColor)
TimeColor <- c("#a6cee3", "#b2df8a", "#fb9a99")
DimPlot(CD8.clean3, reduction = "umap", pt.size = 0.1, label = F, group.by="Time", cols = TimeColor)


## Save object for later
saveRDS(data.integrated, file = "~/Desktop/CD3_and_HLADR_integrated_HARMONY.RDS")


# Making Stacked Bar Plots for Genotype Distributions (Figure 2E)
counts <- with(CD8.clean3@meta.data, table(orig.ident,Genotype))

ggplot(as.data.frame(counts), aes(x = Tissue, y = Freq, fill = Genotype)) + 
  geom_col(
    position = "fill" #Include if you want each bar to sum to 1 for relative abundance
  ) +
  scale_fill_manual(values=GenoColor) +
  ylab("Frequency of Cells") +
  theme(
    #legend.position = "none", #removes the legend
    axis.title.y = element_text(size=12),
    axis.text.y = element_text(size=10),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(angle=90, vjust=0.5, size=10))



#Use SingleR to classify cell types
library(SingleCellExperiment)
library(Seurat)
library(scater)
library(SingleR)
library(scRNAseq)



### Workflow based on vignette
# either read in a Seurat object as RDS, or take any Seurat object and convert it to a singlecell experiment
sce <- as.SingleCellExperiment(CD8.clean3)

# option for subsampling to speed things up (only used for testing)
#sce <- sce[,c(1:10000,90000:100000)]

# reference data sets included in SingleR - the monaco data set has worked ok for me (details about the reference see vignette)
monaco.se <- MonacoImmuneData()
#monaco.se <- NovershternHematopoieticData()

# Grab the common genes between the query and reference data set
common <- intersect(rownames(sce), rownames(monaco.se))
monaco.se <- monaco.se[common,]
sce <- sce[common,]

# Test and reference sets should always be log normalized. The included reference sets are already normalized. 
sce <- logNormCounts(sce)

# Perform the actual SingleR annotation (one is fine labels, one is main labels, havent figured out how to do it in one go)
annotated.sce.fine <- SingleR(test = sce, ref = monaco.se, 
                              labels = monaco.se$label.fine, assay.type.ref = "logcounts")
annotated.sce <- SingleR(test = sce, ref = monaco.se, 
                         labels = monaco.se$label.main, assay.type.ref = "logcounts")
table(annotated.sce$labels)
table(annotated.sce.fine$labels)

# just take the labels from the sce object and turn them into a vector. Use that Vector to put it into the meta-data of your Seurat object
monaco.labels <- annotated.sce$labels
monaco.labels.fine <- annotated.sce.fine$labels

CD8.clean3@meta.data$monaco.labels <- monaco.labels
CD8.clean3@meta.data$monaco.labels.fine <- monaco.labels.fine

DimPlot(CD8.clean3, reduction = "umap", pt.size = 0.1, label = T, group.by="monaco.labels.fine", cols = color.vector)

###LOOK AT MAITS

CD8.clean3 <- SetIdent(CD8.clean3, value = "monaco.labels.fine")
MAITs <- subset(CD8.clean3, idents = 'MAIT cells')

DimPlot(MAITs, reduction = "umap", pt.size = 0.1, label = F, group.by="Genotype", cols = GenoColor)
DimPlot(MAITs, reduction = "umap", pt.size = 0.1, label = F, group.by="Time", cols = TimeColor)
DimPlot(MAITs, reduction = "umap", pt.size = 0.01, label = F, group.by="sample", cols = SampleColor)
FeaturePlot(MAITs, features = c("TRAV1-2", "KLRB1"), 
            pt.size = 0.1, ncol = 2)


counts <- with(MAITs@meta.data, table(Tissue,Genotype))

ggplot(as.data.frame(counts), aes(x = sample, y = Freq, fill = Genotype)) + 
  geom_col(
    position = "fill" #Include if you want each bar to sum to 1 for relative abundance
  ) +
  scale_fill_manual(values=GenoColor) +
  ylab("Frequency of Cells") +
  theme(
    #legend.position = "none", #removes the legend
    axis.title.y = element_text(size=12),
    axis.text.y = element_text(size=10),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(angle=90, vjust=0.5, size=10))

MAITcells <- SetIdent(MAITcells, value = "Time")
markers.MAITs <- FindAllMarkers(MAITcells,
                                features = vargenes.filtered5,
                                only.pos = TRUE, 
                                min.pct = 0.1, 
                                logfc.threshold = 0.25,
                                test.use = "MAST",
                                verbose = T,
                                latent.vars = "nCount_RNA")


MAITcells <- SetIdent(MAITcells, value = "sample")
MAITcells2 <- subset(MAITcells, ident = "Patient2", invert = F)

MAITcells2 <- SetIdent(MAITcells2, value = "Time")
markers.MAITs2 <- FindAllMarkers(MAITcells2,
                                 features = vargenes.filtered5,
                                 only.pos = TRUE, 
                                 min.pct = 0.1, 
                                 logfc.threshold = 0.25,
                                 test.use = "MAST",
                                 verbose = T,
                                 latent.vars = "nCount_RNA")

#write.csv(markers.MAITs2, "markers.MAITs2.csv")



MAITcells <- SetIdent(MAITcells, value = "sample")
MAITs3 <- subset(MAITcells, ident = "Patient3", invert = F)

MAITs3 <- SetIdent(MAITs3, value = "Time")
markers.MAITs3 <- FindAllMarkers(MAITs3,
                                 features =vargenes.filtered5,
                                 only.pos = TRUE, 
                                 min.pct = 0.1, 
                                 logfc.threshold = 0.25,
                                 test.use = "MAST",
                                 verbose = T,
                                 latent.vars = "nCount_RNA")

#write.csv(markers.MAITs3, "markers.MAITs3.csv")




MAITcells <- SetIdent(MAITcells, value = "sample")
MAITs4 <- subset(MAITcells, ident = "Patient4", invert = F)

MAITs4 <- SetIdent(MAITs4, value = "Time")
markers.MAITs4 <- FindAllMarkers(MAITs4,
                                 features = vargenes.filtered5,
                                 only.pos = TRUE, 
                                 min.pct = 0.1, 
                                 logfc.threshold = 0.25,
                                 test.use = "MAST",
                                 verbose = T,
                                 latent.vars = "nCount_RNA")

#write.csv(markers.MAITs4, "markers.MAITs4.csv")


MAITcells <- SetIdent(MAITcells, value = "sample")
MAITs1 <- subset(MAITcells, ident = "Patient1", invert = F)

MAITs1 <- SetIdent(MAITs1, value = "Time")
markers.MAITs1 <- FindAllMarkers(MAITs1,
                                 features = vargenes.filtered5,
                                 only.pos = TRUE, 
                                 min.pct = 0.1, 
                                 logfc.threshold = 0.25,
                                 test.use = "MAST",
                                 verbose = T,
                                 latent.vars = "nCount_RNA")





#### make .csv with barcode, umap coordinates, and monaco lables

barcodeumap <- as.data.frame(sce@int_colData@listData$reducedDims@listData$UMAP)
head(barcodeumap)
tail(barcodeumap)

setDT(barcodeumap, keep.rownames = T)[]
colnames(barcodeumap) <- c('barcode', 'UMAP_1', 'UMAP_2')

head(barcodeumap)
tail(barcodeumap)


table(CD8.clean3@meta.data$monaco.labels.fine)


metadata <- as.data.frame(CD8.clean3@meta.data)
head(metadata)
tail(metadata)

#merge sample and cellIDs

metadata$barcode <- paste(metadata$sample, metadata$CellID, sep = "_")
head(metadata)

setDT(metadata, keep.rownames = T)[]
colnames(metadata) <- c('Number', 'orig.ident', 'nCount_RNA', "nFeature_RNA", "CellID", "Genotype", "Annotations", "ClonotypeID", "ClustermapRowNumber", "SampleID", "Tissue", "Time", "sample", "percent.mt", "prcent.ribo", "RNA_snn_res.0.8", "seurat_clusters", "monaco.labels", "monaco.labels.fine", "barcode")

head(metadata)
tail(metadata)

UMAPMeta <- merge(metadata, barcodeumap, by.x = 'barcode', by.y = 'barcode', all.x = T, all.y = F, sort = F)
head(UMAPMeta)

write.csv(UMAPMeta, "CMV_CD8_UMAP_MonacoLabels.csv")

#Make figure 2C D with filtered cells

CD8.clean3@meta.data$XIST <- CD8.clean3@assays$RNA@data["XIST",]
CD8.clean3@meta.data$RPS4Y1 <- CD8.clean3@assays$RNA@data["RPS4Y1",]

write.csv(CD8.clean3@meta.data, "Meta.csv")

#sexlinked gene abundance 
write.csv(t(CD8_Patient2_Seurat_2@assays$RNA@data[c("RPS4Y1", "XIST"),]), "Patient2.csv")
write.csv(t(CD8_Patient1_Seurat_2@assays$RNA@data[c("RPS4Y1", "XIST"),]), "Patient1.csv")
write.csv(t(CD8_Patient3_Seurat_2@assays$RNA@data[c("RPS4Y1", "XIST"),]), "Patient3.csv")
write.csv(t(CD8_Patient4_Seurat_2@assays$RNA@data[c("RPS4Y1", "XIST"),]), "Patient4.csv")

#Load CMV-tetramer+ CD8+ T cell data
#CMV-tetramer+ CD8+ T cells were added to a sample that also contained NK cells. Once samples are loaded, we will isolate CD8+ T cells from NK cells and remove duplicate barcodes across samples. .

#Patient1
NKCMV_Patient1.data <- Read10X(data.dir = 'Patient1_NKCMV_Aggr')


Genotype1 <- read.csv("Patient1NKCMVDR.csv", header = T)
#read in a CSV files with timepoint, barcode, and genotype columns. 
Genotype1 <- as.matrix(Genotype1)
head(Genotype1)
tail(Genotype1)


cellids10x1 <- as.data.frame(NKCMV_Patient1.data[1,])
head(cellids10x1)
setDT(cellids10x1, keep.rownames = T)[]
colnames(cellids10x1) <- c('barcodes', 'NA')

GenoAnnoPatient1b <- merge(cellids10x1, Genotype1, by.x = 'barcodes', by.y = 'barcode', all.x = T, all.y = F, sort = F)
head(GenoAnnoPatient1b)
tail(GenoAnnoPatient1b)
GenoAnnoPatient1b <- as.matrix(GenoAnnoPatient1b)

cellids1 <- c(GenoAnnoPatient1b[,1])
genotype1 <- c(GenoAnnoPatient1b[,4])



NKCMV_Patient1_Seurat <- CreateSeuratObject(counts = NKCMV_Patient1.data, 
                                         assay = 'RNA',
                                         project = 'CMV_NKCMV_Patient1')
NKCMV_Patient1_Seurat # 33694 features across 15954 samples within 1 assay 
NKCMV_Patient1_Seurat@meta.data$CellID <- cellids1
NKCMV_Patient1_Seurat@meta.data$Genotype <- genotype1
NKCMV_Patient1_Seurat@meta.data$SampleID <- sapply(rownames(NKCMV_Patient1_Seurat@meta.data), function(x) str_split(x, '-')[[1]][2])
NKCMV_Patient1_Seurat@meta.data$Tissue <- plyr::mapvalues(x = NKCMV_Patient1_Seurat@meta.data$SampleID,
                                                       from = c('1', '2', '3'),
                                                       to = c('CD8_30_Patient1', 'CD8_60_Patient1', 'CD8_90_Patient1'))
NKCMV_Patient1_Seurat@meta.data$Time <- plyr::mapvalues(x = NKCMV_Patient1_Seurat@meta.data$SampleID,
                                                     from = c('1', '2', '3'),
                                                     to = c('30', '60', '90'))

NKCMV_Patient1_Seurat <- RenameCells(NKCMV_Patient1_Seurat, add.cell.id = "Patient1")
NKCMV_Patient1_Seurat@meta.data$sample <- "Patient1"
NKCMV_Patient1_Seurat[["percent.mt"]] <- PercentageFeatureSet(NKCMV_Patient1_Seurat, pattern = "^MT-")
NKCMV_Patient1_Seurat[["percent.ribo"]] <- PercentageFeatureSet(NKCMV_Patient1_Seurat, pattern = "^RP")
VlnPlot(NKCMV_Patient1_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
NKCMV_Patient1_Seurat_2 <- subset(NKCMV_Patient1_Seurat, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 3200 & nCount_RNA < 15000)
NKCMV_Patient1_Seurat_2 ##-- 33694 features across 15809 samples within 1 assay.

saveRDS(NKCMV_Patient1_Seurat_2, "NKCMV_Patient1_Seurat_2.RDS")

#Patient 2

NKCMV_Patient2.data <- Read10X(data.dir = 'Patient2_NKCMV_Aggr')

Genotype2 <- read.csv("Patient2NKCMVDR.csv", header = T)
Genotype2 <- as.data.frame(Genotype2)

cellids10x2 <- as.data.frame(NKCMV_Patient2.data[1,])
head(cellids10x2)
library(data.table)
setDT(cellids10x2, keep.rownames = T)[]
colnames(cellids10x2) <- c('barcodes', 'NA')

GenoCellIDPatient2 <- merge(cellids10x2, Genotype2, by.x = 'barcodes', by.y = 'barcode', all.x = T, all.y = F, sort = F)
head(GenoCellIDPatient2)
GenoCellIDPatient2 <- as.matrix(GenoCellIDPatient2)

length(GenoCellIDPatient2)

cellids1 <- c(GenoCellIDPatient2[,1])
genotype1 <- c(GenoCellIDPatient2[,3])

#- Setup seurat object class
NKCMV_Patient2_Seurat <- CreateSeuratObject(counts = NKCMV_Patient2.data,
                                         assay = 'RNA',
                                         project = 'CMV_NKCMV_Patient2')
NKCMV_Patient2_Seurat # 33694 features across 18303 samples within 1 assay 


#- Sample info
NKCMV_Patient2_Seurat@meta.data$CellID <- cellids2
NKCMV_Patient2_Seurat@meta.data$Genotype <- genotype2
NKCMV_Patient2_Seurat@meta.data$SampleID <- sapply(rownames(NKCMV_Patient2_Seurat@meta.data), function(x) str_split(x, '-')[[1]][2])
NKCMV_Patient2_Seurat@meta.data$Tissue <- plyr::mapvalues(x = NKCMV_Patient2_Seurat@meta.data$SampleID,
                                                       from = c('1', '2', '3'),
                                                       to = c('NKCMV_30_Patient2', 'NKCMV_60_Patient2', 'NKCMV_90_Patient2'))
NKCMV_Patient2_Seurat[["percent.mt"]] <- PercentageFeatureSet(NKCMV_Patient2_Seurat, pattern = "^MT-")
NKCMV_Patient2_Seurat[["percent.ribo"]] <- PercentageFeatureSet(NKCMV_Patient2_Seurat, pattern = "^RP")
VlnPlot(NKCMV_Patient2_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
NKCMV_Patient2_Seurat_2 <- subset(NKCMV_Patient2_Seurat, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA < 8000)

# -Inf and Inf should be used if you don't want a lower or upper threshold
NKCMV_Patient2_Seurat_2 ##-- 33694 features across 18137 samples within 1 assay .

saveRDS(NKCMV_Patient2_Seurat_2, "NKCMV_Patient2_Seurat_2.RDS")


#Patient 3

NKCMV_Patient3.data <- Read10X(data.dir = 'Patient3_NKCMV_Aggr')

Genotype3 <- read.csv("Patient3NKCMVDR.csv", header = T)
Genotype3 <- as.data.frame(Genotype3)

cellids10x3 <- as.data.frame(NKCMV_Patient3.data[1,])
head(cellids10x3)
library(data.table)
setDT(cellids10x3, keep.rownames = T)[]
colnames(cellids10x3) <- c('barcodes', 'NA')

GenoCellIDPatient3 <- merge(cellids10x3, Genotype3, by.x = 'barcodes', by.y = 'barcode', all.x = T, all.y = F, sort = F)
head(GenoCellIDPatient3)
GenoCellIDPatient3 <- as.matrix(GenoCellIDPatient3)

length(GenoCellIDPatient3)

cellids3 <- c(GenoCellIDPatient3[,1])
genotype3 <- c(GenoCellIDPatient3[,3])

#- Setup seurat object class
NKCMV_Patient3_Seurat <- CreateSeuratObject(counts = NKCMV_Patient3.data,
                                            assay = 'RNA',
                                            project = 'CMV_NKCMV_Patient3')
NKCMV_Patient3_Seurat 


#- Sample info
NKCMV_Patient3_Seurat@meta.data$CellID <- cellids3
NKCMV_Patient3_Seurat@meta.data$Genotype <- genotype3
NKCMV_Patient3_Seurat@meta.data$SampleID <- sapply(rownames(NKCMV_Patient3_Seurat@meta.data), function(x) str_split(x, '-')[[1]][2])
NKCMV_Patient3_Seurat@meta.data$Tissue <- plyr::mapvalues(x = NKCMV_Patient3_Seurat@meta.data$SampleID,
                                                          from = c('1', '2', '3'),
                                                          to = c('NKCMV_30_Patient3', 'NKCMV_60_Patient3', 'NKCMV_90_Patient3'))
NKCMV_Patient3_Seurat[["percent.mt"]] <- PercentageFeatureSet(NKCMV_Patient3_Seurat, pattern = "^MT-")
NKCMV_Patient3_Seurat[["percent.ribo"]] <- PercentageFeatureSet(NKCMV_Patient3_Seurat, pattern = "^RP")
VlnPlot(NKCMV_Patient3_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
NKCMV_Patient3_Seurat_2 <- subset(NKCMV_Patient3_Seurat, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA < 8000)


NKCMV_Patient3_Seurat_2

saveRDS(NKCMV_Patient3_Seurat_2, "NKCMV_Patient3_Seurat_2.RDS")



#Patient 4

NKCMV_Patient4.data <- Read10X(data.dir = 'Patient4_NKCMV_Aggr')

Genotype4 <- read.csv("Patient4NKCMVDR.csv", header = T)
Genotype4 <- as.data.frame(Genotype4)

cellids10x4 <- as.data.frame(NKCMV_Patient4.data[1,])
head(cellids10x4)
library(data.table)
setDT(cellids10x4, keep.rownames = T)[]
colnames(cellids10x4) <- c('barcodes', 'NA')

GenoCellIDPatient4 <- merge(cellids10x4, Genotype4, by.x = 'barcodes', by.y = 'barcode', all.x = T, all.y = F, sort = F)
head(GenoCellIDPatient4)
GenoCellIDPatient4 <- as.matrix(GenoCellIDPatient4)

length(GenoCellIDPatient4)

cellids4 <- c(GenoCellIDPatient4[,1])
genotype4 <- c(GenoCellIDPatient4[,3])

#- Setup seurat object class
NKCMV_Patient4_Seurat <- CreateSeuratObject(counts = NKCMV_Patient4.data,
                                            assay = 'RNA',
                                            project = 'CMV_NKCMV_Patient4')
NKCMV_Patient4_Seurat 


#- Sample info
NKCMV_Patient4_Seurat@meta.data$CellID <- cellids4
NKCMV_Patient4_Seurat@meta.data$Genotype <- genotype4
NKCMV_Patient4_Seurat@meta.data$SampleID <- sapply(rownames(NKCMV_Patient4_Seurat@meta.data), function(x) str_split(x, '-')[[1]][2])
NKCMV_Patient4_Seurat@meta.data$Tissue <- plyr::mapvalues(x = NKCMV_Patient4_Seurat@meta.data$SampleID,
                                                          from = c('1', '2', '3'),
                                                          to = c('NKCMV_30_Patient4', 'NKCMV_60_Patient4', 'NKCMV_90_Patient4'))
NKCMV_Patient4_Seurat[["percent.mt"]] <- PercentageFeatureSet(NKCMV_Patient4_Seurat, pattern = "^MT-")
NKCMV_Patient4_Seurat[["percent.ribo"]] <- PercentageFeatureSet(NKCMV_Patient4_Seurat, pattern = "^RP")
VlnPlot(NKCMV_Patient4_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
NKCMV_Patient3_Seurat_4 <- subset(NKCMV_Patient4_Seurat, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA < 8000)


NKCMV_Patient4_Seurat_2

saveRDS(NKCMV_Patient4_Seurat_2, "NKCMV_Patient4_Seurat_2.RDS")

#Merge NKCMV
NKCMV <- merge(NKCMV_Patient2_Seurat_2, c(NKCMV_Patient3_Seurat_2, NKCMV_Patient4_Seurat_2, NKCMV_Patient1_Seurat_2))
NKCMV #33694 features across 65100 samples within 1 assay  

table(NKCMV@meta.data$Tissue)
table(NKCMV_Patient1_Seurat_2@meta.data$Tissue)
table(NKCMV_Patient2_Seurat_2@meta.data$Tissue)
table(NKCMV_Patient3_Seurat_2@meta.data$Tissue)
table(NKCMV_Patient4_Seurat_2@meta.data$Tissue)



############# Harmony workflow
## Generate merged object with all data for HARMONY integration
NKCMV <- NormalizeData(NKCMV, verbose = T)
NKCMV <- FindVariableFeatures(NKCMV, selection.method = "vst", nfeatures = 10000, verbose = T)
NKCMV <- ScaleData(NKCMV, verbose = T)

## take sex-linked genes out of var.features

var.genes <- as.vector(NKCMV@assays$RNA@var.features)
genes.filtered.xist <- grep(pattern = "XIST", var.genes, value = TRUE, invert = TRUE)
vargenes.filtered <- as.vector(grep(pattern = "RPS4Y1", genes.filtered.xist, value = TRUE, invert = TRUE))
vargenes.filtered2 <- as.vector(grep(pattern = "^TTTY", vargenes.filtered, value = TRUE, invert = TRUE))
vargenes.filtered3 <- as.vector(grep(pattern = "UTY", vargenes.filtered2, value = TRUE, invert = TRUE))
vargenes.filtered4 <- as.vector(grep(pattern = "ZFY", vargenes.filtered3, value = TRUE, invert = TRUE))
vargenes.filtered5 <- as.vector(grep(pattern = "EIF1AY", vargenes.filtered4, value = TRUE, invert = TRUE))




write.csv(vargenes.filtered, "vargenes.csv")

NKCMV <- RunPCA(NKCMV, pc.genes = vargenes.filtered5, npcs = 50, verbose = T)

NKCMV <- RunHarmony(NKCMV, group.by.vars = "sample",
                  lambda = 1, theta = 2, plot_convergence = T, verbose = T, max.iter.harmony = 20)
NKCMV
table(NKCMV$sample)
table(NKCMV$Time)
table(NKCMV@meta.data$Genotype)

ElbowPlot(object = NKCMV,
          ndims = 50)

## Perform UMAP calculation and graph-based clustering #switch "pca" to "harmony" if running harmony. 
NKCMV <- RunUMAP(NKCMV, reduction = "harmony", dims = 1:40, seed.use = 34)
NKCMV <- FindNeighbors(NKCMV, reduction = "harmony", dims = 1:40, verbose = T)
NKCMV <- FindClusters(NKCMV, resolution = 0.8, random.seed = 34, verbose = T)



## Plotting
color.vector <- c(brewer.pal(8, "Set1"), brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))
DimPlot(NKCMV, reduction = "umap", pt.size = 0.1, label = F, group.by="sample", cols = color.vector)
DimPlot(NKCMV, reduction = "umap", pt.size = 0.1, label = F, group.by = "Tissue", cols = color.vector)
DimPlot(NKCMV, reduction = "umap", pt.size = 0.1, label = T, group.by = "seurat_clusters", cols = color.vector)
DimPlot(NKCMV, reduction = "umap", pt.size = 0.1, label = F, group.by="orig.ident", cols = brewer.pal(4, "Set1"))
DimPlot(NKCMV, reduction = "umap", pt.size = 0.1, label = F, group.by="Genotype", cols = brewer.pal(3, "Set1"))
DimPlot(NKCMV, reduction = "umap", pt.size = 0.1, label = F, group.by="sample", cols = color.vector)


#Filter to select T cells (CMV-Tet+ T Cells)
FeaturePlot(NKCMV, features = c("CD3D", "CD3E", "CD8B", "CD8A", "NCAM1", "KLRK1", "NCR1", "FCGR3A"), 
            pt.size = 0.1, ncol = 3)

CMV <- subset(NKCMV, idents = c(3, 4, 8), invert = F)
CMV <- RunUMAP(CMV, reduction = "pca", dims = 1:20, seed.use = 34)
CMV <- FindNeighbors(CMV, reduction = "pca", dims = 1:20, verbose = T)
CMV <- FindClusters(CMV, resolution = 0.8, random.seed = 34, verbose = T)

color.vector <- c(brewer.pal(8, "Set1"), brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))
DimPlot(CMV, reduction = "umap", pt.size = 0.1, label = T, group.by="seurat_clusters", cols = color.vector)

FeaturePlot(CMV, features = c("CD3D", "CD3E", "CD8B", "CD8A", "NCAM1", "KLRK1", "NCR1", "FCGR3A"), 
            pt.size = 0.1, ncol = 3)

CMV.clean <- subset(CMV, idents = c(3, 6), invert = T) #cluster number may be different -- ensure taking T cell clusters
CMV.clean <- RunUMAP(CMV.clean, reduction = "pca", dims = 1:20, seed.use = 34)
CMV.clean <- FindNeighbors(CMV.clean, reduction = "pca", dims = 1:20, verbose = T)
CMV.clean <- FindClusters(CMV.clean, resolution = 3, random.seed = 34, verbose = T)

DimPlot(CMV.clean, reduction = "umap", pt.size = 0.8, label = T, group.by="seurat_clusters", cols = color.vector)

FeaturePlot(CMV2.clean, features = c("CD3D", "CD3E", "CD8B", "CD8A", "NCAM1", "KLRK1", "NCR1", "FCGR3A"), 
            pt.size = 0.1, ncol = 3)


CMV.clean2 <- subset(CMV.clean, idents = c(13), invert = T)
CMV.clean2 <- RunUMAP(CMV.clean2, reduction = "pca", dims = 1:20, seed.use = 34)
CMV.clean2 <- FindNeighbors(CMV.clean2, reduction = "pca", dims = 1:20, verbose = T)
CMV.clean2 <- FindClusters(CMV.clean2, resolution = 0.8, random.seed = 34, verbose = T)

DimPlot(CMV.clean2, reduction = "umap", pt.size = 0.1, label = T, group.by="seurat_clusters", cols = color.vector)

FeaturePlot(CMV2.clean2, features = c("CD3D", "CD3E", "CD8B", "CD8A", "NCAM1", "KLRK1", "NCR1", "FCGR3A"), 
            pt.size = 0.1, ncol = 3)

DimPlot(CMV.clean2, reduction = "umap", pt.size = 1, label = F, group.by="Genotype", cols = color.vector)

#make table
table(CMV.clean2@meta.data$Genotype)

table(CMV.clean2@meta.data$Genotype, CMV2.clean@meta.data$Tissue)

GenColors <- c("#c7c2e3", "#276419")

counts <- with(CMV.clean2@meta.data, table(Tissue,Genotype))

ggplot(as.data.frame(counts), aes(x = Tissue, y = Freq, fill = Genotype)) + 
  geom_col(
    position = "fill" #Include if you want each bar to sum to 1 for relative abundance
  ) +
  scale_fill_manual(values=GenColors) +
  ylab("Frequency of Cells") +
  theme(
    #legend.position = "none", #removes the legend
    axis.title.y = element_text(size=12),
    axis.text.y = element_text(size=10),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(angle=90, vjust=0.5, size=10)) #390x770


write.csv(CMV2.clean2@meta.data$CellID, "CellID2.csv")
write.csv(CMV2.clean2@assays$RNA@scale.data[1,], "scaledatachr2.csv")


CMV@meta.data$sample <- "Tet"

saveRDS(CMV2.clean2, "Patient3_CMVTet.RDS")

CMV.clean <- subset(NKCMV_Patient3_Seurat_3, idents = c(9), invert = F)
CMV.clean@meta.data$sample <- "Tet"

DimPlot(CMV.clean, reduction = "umap", pt.size = 1, label = F, group.by="Genotype", cols = color.vector)





##-- Load the Demo PBMC dataset, required format is the unzipped file from the CellRanger workflow
DMSOCD137neg.data <- Read10X(data.dir = 'DMSO_CD137neg')
IE1CD137pos.data <- Read10X(data.dir = 'IE1_CD137pos')
pp65CD137pos.data <- Read10X(data.dir = 'pp65_CD137pos')

Genotype <- read.csv("Patient3stimDP.csv", header = T)
Genotype <- as.matrix(Genotype)

cellids <- c(Genotype[,1])
genotype <- c(Genotype[,2])

#- Setup seurat object class
DMSOCD137neg <- CreateSeuratObject(DMSOCD137neg.data, names.field = 1,
                                   project = 'Patient3 Peptide Stim')
DMSOCD137neg # 33694 features across 2111 samples.
DMSOCD137neg@meta.data$SampleID <- "1"
DMSOCD137neg@meta.data$Condition <- "Control"

IE1CD137pos <- CreateSeuratObject(IE1CD137pos.data, names.field = 1,
                                  project = 'Patient3 Peptide Stim')
IE1CD137pos # 33694 features across 2111 samples.
IE1CD137pos@meta.data$SampleID <- "2"
IE1CD137pos@meta.data$Condition <- "Stim"

pp65CD137pos <- CreateSeuratObject(pp65CD137pos.data, names.field = 1,
                                   project = 'Patient3 Peptide Stim')
pp65CD137pos # 33694 features across 2111 samples.
pp65CD137pos@meta.data$SampleID <- "3"
pp65CD137pos@meta.data$Condition <- "Stim"

#merge
Pepstim <- merge(DMSOCD137neg, y = c(IE1CD137pos, pp65CD137pos), add.cell.ids = c("DMSOCD137neg", "IE1CD137pos", "pp65CD137pos"), project = "Patient3 Peptide Stim")

#- Sample info
Pepstim@meta.data$CellID <- cellids
Pepstim@meta.data$Genotype <- genotype
#Pepstim@meta.data$SampleID <- sapply(rownames(Pepstim@meta.data), function(x) str_split(x, '-')[[1]][2])
Pepstim@meta.data$Tissue <- plyr::mapvalues(x = Pepstim@meta.data$SampleID,
                                                 from = c('1', '2', '3'),
                                                 to = c('DMSOCD137neg', 'IE1CD137pos', 'pp65CD137pos'))

table(Pepstim@meta.data$SampleID)

with(Pepstim@meta.data, table(Tissue,genotype))


##-- QC metrics - Percentage of mitochondrial genes
#- percent.mito variable
mito.genes <- grep(pattern = '^MT-', x = rownames(x = Pepstim@assays$RNA@data), value = TRUE) # 13 mitochondrial genes
percent.mito <- Matrix::colSums(Pepstim@assays$RNA@counts[mito.genes, ]) / Matrix::colSums(Pepstim@assays$RNA@counts) # Matrix package is needed


#- Use of AddMetaData to add column to object@meta.data (great place to store QC metrics)
Pepstim <- AddMetaData(object = Pepstim, 
                            metadata = percent.mito, 
                            col.name = 'percent.mito')



##-- QC metrics - draw a VlnPlot()
VlnPlot(object = Pepstim, 
        features = c('nFeature_RNA', 'nCount_RNA', 'percent.mito'), 
        ncol = 3, pt.size =0.1, group.by = "Tissue")

##-- QC metrics - another QC metric using GenePlot(), correlation displayed on top
par(mfrow = c(1, 2))
FeatureScatter(object = Pepstim, feature1 = 'nCount_RNA', feature2 = 'percent.mito', group.by = "Tissue")
FeatureScatter(object = Pepstim, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', group.by = "Tissue")

##-- Filtering (removing the outlier cells with percent.mito being higher than 10%
Pepstim_2 <- subset(Pepstim, subset = percent.mito < 0.2)
# -Inf and Inf should be used if you don't want a lower or upper threshold
Pepstim_2 ##-- 33694 features across 2311 samples

#Re-plot after filtering
par(mfrow = c(1, 2))
FeatureScatter(object = Pepstim_2, feature1 = 'nCount_RNA', feature2 = 'percent.mito', group.by = "Tissue")
FeatureScatter(object = Pepstim_2, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', group.by = "Tissue")


##-- Normalization, see R Markdownsheet
Pepstim_2 <- NormalizeData(object = Pepstim_2, 
                                normalization.method = 'LogNormalize', 
                                scale.factor = 10000)

##-- Highly variable genes, y is the dispersion/variance, cutoffs need to be set based on dataset
Pepstim_2 <- FindVariableFeatures(object = Pepstim_2, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Pepstim_2), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Pepstim_2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
##-- Highly variable genes, number
length(Pepstim_2@assays$RNA@var.features) #5000

##-- Scaling the data and removing unwanted sources of variation - this function takes about 1-2 minutes, might be longer
Pepstim_2 <- ScaleData(object = Pepstim_2, 
                            model.use = 'linear',
                            vars.to.regress = c('nUMI', 'percent.mito')) # nUMI is used as proxy for Cellular Detection Rate
##-- add batch if you have different batches of data. Usually regressing for number of UMI sufficient

##-- perform PCA using the scale-corrected matrix, pcs.compute is the number of principal components
Pepstim_2a <- RunPCA(object = Pepstim_2, 
                          pc.genes = Pepstim_2@assays$RNA@var.features, 
                          pcs.compute = 50,
                          do.print = FALSE)

##-- PCA - PrintPCA
print(Pepstim_2a[["pca"]], dims = 1:5, nfeatures = 5) # nfeatures is number of genes to print for each PC

##-- PCA - VizPCA, plot top 30 genes for PC1 and PC2, 
##-- check whether one of the PCs is driven my mitochondrial genes or cell cycle-genes
VizDimLoadings(Pepstim_2a, dims = 1:2, reduction = "pca", nfeatures = 30)


Pepstim_2a <- SetAllIdent(Pepstim_2a, id = "Genotype")
##-- PCA - PCAPlot
PCAPlot(object = Pepstim_2a, 
        dim.1 = 1, 
        dim.2 = 2)

##-- PCElbowPlot, plots the standard deviation for each PC, look at the "elbow" of the curve (JackStraw procedure would take much longer)
ElbowPlot(object = Pepstim_2a,
          ndims = 50)



##-- UMAP
Pepstim_2a <- FindNeighbors(object = Pepstim_2a, reduction = 'pca', dims = 1:33, k.param = 30, force.recalc = T)
Pepstim_2a <- FindClusters(object = Pepstim_2a, verbose = TRUE, n.start = 100)
Pepstim_2d <- RunUMAP(Pepstim_2a, reduction.use = "pca", dims = 1:33, seed.use = 34)
Pepstim_2d <- SetIdent(Pepstim_2d, value = "Genotype")
DimPlot(Pepstim_2d, reduction = "umap", cols = brewer.pal(3, "Paired"), pt.size = 0.8)

with(Pepstim_2d@meta.data, table(Tissue,Genotype))
with(Pepstim_2d@meta.data, table(Condition,Genotype))

Pepstim_2d <- SetIdent(Pepstim_2d, value = "Tissue")
DimPlot(Pepstim_2d, reduction.use = "umap", cols = brewer.pal(9, "Set1"), pt.size = 0.8)

Pepstim_2d <- SetIdent(Pepstim_2d, value = "RNA_snn_res.0.8")
DimPlot(Pepstim_2d, reduction.use = "umap", cols = c('#cb181d', '#f16913', '#ffff33', '#006837', '#66c2a4', '#081d58', '#41b6c4', '#810f7c', '#807dba'), pt.size = 0.8)

FeaturePlot(Pepstim_2d, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3)


# Example of Making Stacked Bar Plots for Cluster Distributions
colors <- c(brewer.pal(10, "Paired"), brewer.pal(8, "Dark2"))
counts <- with(Pepstim_2d@meta.data, table(Tissue,Genotype))

ggplot(as.data.frame(counts), aes(x = Tissue, y = Freq, fill = Genotype)) + 
  geom_col(
    position = "fill" #Include if you want each bar to sum to 1 for relative abundance
  ) +
  scale_fill_manual(values=colors) +
  ylab("Frequency of Cells") +
  theme(
    #legend.position = "none", #removes the legend
    axis.title.y = element_text(size=12),
    axis.text.y = element_text(size=10),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(angle=90, vjust=0.5, size=10))





#compare differentiation states between CMV-tetramer+ cells and CD8s

#read in NKCMV

Patient1Tet <- readRDS("Patient1_CMVTet.RDS")

Patient2Tet <- readRDS("Patient2_CMVTet.RDS")

Patient3Tet <- readRDS("Patient3_CMVTet.RDS")

Patient4Tet <- readRDS("Patient4_CMVTet.RDS")

#read in CD8s

CD8_Patient1_Seurat_2 <- readRDS("CD8_Patient1_Seurat_2.RDS")
CD8_Patient2_Seurat_2 <- readRDS("CD8_Patient2_Seurat_2.RDS")
CD8_Patient3_Seurat_2 <- readRDS("CD8_Patient3_Seurat_2.RDS")
CD8_Patient4_Seurat_2 <- readRDS("CD8_Patient4_Seurat_2.RDS")

Patient1Tet@meta.data$CellIdent <- 'CMVTet'
Patient1Tet@meta.data$Tissue <- plyr::mapvalues(x = Patient1Tet@meta.data$Tissue,
                                                          from = c('CD8_30_Patient1', 'CD8_60_Patient1', 'CD8_90_Patient1'),
                                                          to = c('NKCMV_30_Patient1', 'NKCMV_60_Patient1', 'NKCMV_90_Patient1'))

Patient2Tet@meta.data$CellIdent <- 'CMVTet'
Patient3Tet@meta.data$CellIdent <- 'CMVTet'
Patient4Tet@meta.data$CellIdent <- 'CMVTet'

CD8_Patient1_Seurat_2@meta.data$CellIdent <- 'CD8'
CD8_Patient2_Seurat_2@meta.data$CellIdent <- 'CD8'
CD8_Patient3_Seurat_2@meta.data$CellIdent <- 'CD8'
CD8_Patient4_Seurat_2@meta.data$CellIdent <- 'CD8'

## Generate merged CD3 object
CD8wCMVtet <- merge(CD8_Patient1_Seurat_2, c(CD8_Patient2_Seurat_2, CD8_Patient3_Seurat_2, CD8_Patient4_Seurat_2, Patient1Tet, Patient2Tet, Patient3Tet, Patient4Tet))
CD8wCMVtet #33694 features across 68368 samples within 1 assay  

table(CD8wCMVtet@meta.data$Tissue)


############# Harmony workflow
## Generate merged object with all data for HARMONY integration
CD8wCMVtet <- NormalizeData(CD8wCMVtet, verbose = T)
CD8wCMVtet <- FindVariableFeatures(CD8wCMVtet, selection.method = "vst", nfeatures = 10000, verbose = T)
CD8wCMVtet <- ScaleData(CD8wCMVtet, verbose = T)

## take sex-linked genes out of var.features

var.genes <- as.vector(CD8wCMVtet@assays$RNA@var.features)
genes.filtered.xist <- grep(pattern = "XIST", var.genes, value = TRUE, invert = TRUE)
vargenes.filtered <- as.vector(grep(pattern = "RPS4Y1", genes.filtered.xist, value = TRUE, invert = TRUE))
vargenes.filtered2 <- as.vector(grep(pattern = "^TTTY", vargenes.filtered, value = TRUE, invert = TRUE))
vargenes.filtered3 <- as.vector(grep(pattern = "UTY", vargenes.filtered2, value = TRUE, invert = TRUE))
vargenes.filtered4 <- as.vector(grep(pattern = "ZFY", vargenes.filtered3, value = TRUE, invert = TRUE))
vargenes.filtered5 <- as.vector(grep(pattern = "EIF1AY", vargenes.filtered4, value = TRUE, invert = TRUE))




#write.csv(vargenes.filtered, "vargenes.csv")

CD8wCMVtet <- RunPCA(CD8wCMVtet, pc.genes = vargenes.filtered5, npcs = 50, verbose = T)

CD8wCMVtet <- RunHarmony(CD8wCMVtet, group.by.vars = "sample",
                  lambda = 1, theta = 2, plot_convergence = T, verbose = T, max.iter.harmony = 20)


ElbowPlot(object = CD8wCMVtet,
          ndims = 50)

## Perform UMAP calculation and graph-based clustering #switch "pca" to "harmony" if running harmony. 
CD8wCMVtet <- RunUMAP(CD8wCMVtet, reduction = "harmony", dims = 1:30, seed.use = 34)
CD8wCMVtet <- FindNeighbors(CD8wCMVtet, reduction = "harmony", dims = 1:30, verbose = T)
CD8wCMVtet <- FindClusters(CD8wCMVtet, resolution = 0.8, random.seed = 34, future.seed = T, verbose = T)



## Plotting
color.vector <- c(brewer.pal(8, "Set1"), brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))
DimPlot(CD8wCMVtet, reduction = "umap", pt.size = 0.1, label = F, group.by="sample", cols = color.vector)
DimPlot(CD8wCMVtet, reduction = "umap", pt.size = 0.1, label = F, group.by = "Tissue", cols = color.vector)
DimPlot(CD8wCMVtet, reduction = "umap", pt.size = 0.1, label = T, group.by = "seurat_clusters", cols = color.vector)
DimPlot(CD8wCMVtet, reduction = "umap", pt.size = 0.1, label = F, group.by="orig.ident", cols = color.vector)
DimPlot(CD8wCMVtet, reduction = "umap", pt.size = 0.1, label = F, group.by="Genotype", cols = brewer.pal(3, "Set1"))
DimPlot(CD8wCMVtet, reduction = "umap", pt.size = 0.1, label = F, group.by="sample", cols = color.vector)



FeaturePlot(CD8wCMVtet, features = c("CD14", "ITGAX", "IL3RA", "CD79A", "IGHM", "CD3E", "FOXP3", "CD4", "CD8B"), 
            pt.size = 0.1, ncol = 3)

CD8wCMVtet <- SetIdent(CD8wCMVtet, value = "RNA_snn_res.0.8")
CD8wCMVtet.clean <- subset(CD8wCMVtet, idents = "16", invert = T)
plot <- DimPlot(CD8wCMVtet.clean, reduction = "umap", cols = color.vector, pt.size = 0.8, label = T)
HoverLocator(plot = plot, information = FetchData(CD8wCMVtet.clean, vars = c("CellID", "Tissue", 'Genotype')))


CD8wCMVtet.clean <- SetIdent(CD8wCMVtet.clean, value = "CellID")
CD8wCMVtet.clean2 <- subset(CD8wCMVtet.clean, ident = c('AGAATAGTCTCTAGGA-3', 'TACGGTAAGGAGTCTG-1', 'TCAGGTATCCAGTAGT-3'), invert = T)

DimPlot(CD8wCMVtet.clean2, reduction = "umap", pt.size = 0.1, label = F, group.by="Genotype", cols = brewer.pal(3, "Set1"))


CD8wCMVtet.clean2 <- SetIdent(CD8wCMVtet.clean2, value = "Genotype")
CD8wCMVtet.clean3 <- subset(CD8wCMVtet.clean2, ident = 'Unknown', invert = T)

GenoColor <- c("#C7C2E3", "#276419")
DimPlot(CD8wCMVtet.clean3, reduction = "umap", pt.size = 0.1, label = F, group.by="Genotype", cols = GenoColor)
SampleColor <- c("#e7298a", "#1b9e77", "#d95f02", "#7570b3")
DimPlot(CD8wCMVtet.clean3, reduction = "umap", pt.size = 0.1, label = F, group.by="sample", cols = SampleColor)
TimeColor <- c("#a6cee3", "#b2df8a", "#fb9a99")
DimPlot(CD8wCMVtet.clean3, reduction = "umap", pt.size = 0.1, label = F, group.by="Time", cols = TimeColor)

DimPlot(CD8wCMVtet.clean3, reduction = "umap", pt.size = 0.1, label = F, group.by="monaco.labels.fine", cols = color.vector)
DimPlot(CD8wCMVtet.clean3, reduction = "umap", pt.size = 0.1, label = F, group.by="CellIdent", cols = c("Dark Gray", "Blue"), order="CMVTet")


## Save object for later
saveRDS(CD8wCMVtet.clean3, file = "CD8wCMVtet.clean3.RDS")

### Workflow based on vignette
# either read in a Seurat object as RDS, or take any Seurat object and convert it to a singlecell experiment
sce <- as.SingleCellExperiment(CD8wCMVtet.clean3)

# option for subsampling to speed things up (only used for testing)
#sce <- sce[,c(1:10000,90000:100000)]

# reference data sets included in SingleR - the monaco data set has worked ok for me (details about the reference see vignette)
monaco.se <- MonacoImmuneData()
#monaco.se <- NovershternHematopoieticData()

# Grab the common genes between the query and reference data set
common <- intersect(rownames(sce), rownames(monaco.se))
monaco.se <- monaco.se[common,]
sce <- sce[common,]

# Test and reference sets should always be log normalized. The included reference sets are already normalized. 
sce <- logNormCounts(sce)

# Perform the actual SingleR annotation (one is fine labels, one is main labels, havent figured out how to do it in one go)
annotated.sce.fine <- SingleR(test = sce, ref = monaco.se, 
                              labels = monaco.se$label.fine, assay.type.ref = "logcounts")
annotated.sce <- SingleR(test = sce, ref = monaco.se, 
                         labels = monaco.se$label.main, assay.type.ref = "logcounts")
table(annotated.sce$labels)
table(annotated.sce.fine$labels)

# just take the labels from the sce object and turn them into a vector. Use that Vector to put it into the meta-data of your Seurat object
monaco.labels <- annotated.sce$labels
monaco.labels.fine <- annotated.sce.fine$labels

CD8wCMVtet.clean3@meta.data$monaco.labels <- monaco.labels
CD8wCMVtet.clean3@meta.data$monaco.labels.fine <- monaco.labels.fine

DimPlot(CD8wCMVtet.clean3, reduction = "umap", pt.size = 0.1, label = F, group.by="monaco.labels.fine", cols = color.vector)

CD8wCMVtet.clean3 <- SetIdent(CD8wCMVtet.clean3, value = "Tissue")
CMVTet1 <- subset(CD8wCMVtet.clean3, idents = c("NKCMV_30_Patient1", "NKCMV_60_Patient1", "NKCMV_90_Patient1"))

DimPlot(CMVTet1, reduction = "umap", pt.size = 0.1, label = F, group.by="monaco.labels.fine", cols = c("#e41a1c", '#377eb8', '#984ea3', '#ffff33', '#a65628', '#8dd3c7', '#80b1d3', '#fbd462', '#b3de69', '#fccde5', '#d9d9d9', "#bc80bd", "#ccebc5"))
DimPlot(CMVTet1, reduction = "umap", pt.size = 0.1, label = F, group.by="Time", cols = TimeColor)


CD8wCMVtet.clean3 <- SetIdent(CD8wCMVtet.clean3, value = "Tissue")
CMVTet3 <- subset(CD8wCMVtet.clean3, idents = c("NKCMV_30_Patient3", "NKCMV_60_Patient3", "NKCMV_90_Patient3"))

DimPlot(CMVTet3, reduction = "umap", pt.size = 0.1, label = F, group.by="monaco.labels.fine", cols = c("#e41a1c", '#377eb8', '#984ea3', '#f781bf', '#8dd3c7', '#fb8072', '#80b1d3', '#fbd462', "#ccebc5"))
DimPlot(CMVTet3, reduction = "umap", pt.size = 0.1, label = F, group.by="Tissue", cols = TimeColor)

CD8wCMVtet.clean3 <- SetIdent(CD8wCMVtet.clean3, value = "Tissue")
CMVTet2 <- subset(CD8wCMVtet.clean3, idents = c("NKCMV_30_Patient2", "NKCMV_60_Patient2", "NKCMV_90_Patient2"))

DimPlot(CMVTet2, reduction = "umap", pt.size = 0.1, label = F, group.by="monaco.labels.fine", cols = c("#e41a1c", '#377eb8', '#984ea3', '#f781bf', '#8dd3c7', '#80b1d3', '#fbd462', "#ccebc5"))
DimPlot(CMVTet2, reduction = "umap", pt.size = 0.1, label = F, group.by="Tissue", cols = TimeColor)

CD8wCMVtet.clean3 <- SetIdent(CD8wCMVtet.clean3, value = "Tissue")
CMVTet4 <- subset(CD8wCMVtet.clean3, idents = c("NKCMV_30_Patient4", "NKCMV_90_Patient4"))

DimPlot(CMVTet4, reduction = "umap", pt.size = 0.1, label = F, group.by="monaco.labels.fine", cols = c('#377eb8', '#984ea3', '#f781bf', '#8dd3c7', '#fb8072', '#80b1d3', '#fbd462', "#ccebc5"))
DimPlot(CMVTet4, reduction = "umap", pt.size = 0.1, label = F, group.by="Tissue", cols =  c("#a6cee3", "#fb9a99"))

#DEGs for similar clones

Recip1 =  c('ACTATCTTCTCTGCTG-2', 'AGTGAGGCAAGTTGTC-2', 'GAAACTCGTCCATGAT-2', 'GAACATCCATCAGTCA-2', 'GACGTTAGTCGCGAAA-2', 'GTGCTTCCATGCCTTC-2', 'AAAGTAGAGACTAGAT-3', 'ACTGTCCGTCCGTCAG-3', 'CAGGTGCAGAGGACGG-3', 'CCGTTCACAGTATAAG-3', 'TGAGAGGAGCGCCTTG-3', 'TTTGGTTAGCGTTGCC-3')
Don1 =  c('AACCATGTCTGGCGAC-3', 'ACTGCTCAGCTGATAA-3', 'AGAGTGGCATAGGATA-3', 'ATCGAGTTCCGGCACA-3', 'CAGCATAGTGAAGGCT-3', 'CAGCTGGTCACATGCA-3', 'CCGGGATCAGCTGCTG-3', 'GCTGCAGTCCGCGCAA-3', 'TCAACGAGTGTAAGTA-3', 'TCCACACGTGAACCTT-3', 'TCGTACCAGGGTTTCT-3', 'TCGTAGAAGGGCTTGA-3', 'TGGACGCGTACGACCC-3', 'TTCGAAGTCATTCACT-3', 'TCAGATGTCAAGGTAA-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
markers.clone1 <- FindAllMarkers(Clone1,
                                 features = vargenes.filtered5,
                                 only.pos = FALSE, 
                                 min.pct = 0.01, 
                                 logfc.threshold = 0.25,
                                 test.use = "MAST",
                                 verbose = T,
                                 latent.vars = "nCount_RNA")

setDT(markers.clone1, keep.rownames = T)[]
colnames(markers.clone1) <- c('gene', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj')


markers.clone1top20 <- markers.clone1 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


Clone1 <- subset(CD8.clean3, ident = c('ACTATCTTCTCTGCTG-2', 'AGTGAGGCAAGTTGTC-2', 'GAAACTCGTCCATGAT-2', 'GAACATCCATCAGTCA-2', 'GACGTTAGTCGCGAAA-2', 'GTGCTTCCATGCCTTC-2', 'AAAGTAGAGACTAGAT-3', 'ACTGTCCGTCCGTCAG-3', 'CAGGTGCAGAGGACGG-3', 'CCGTTCACAGTATAAG-3', 'TGAGAGGAGCGCCTTG-3', 'TTTGGTTAGCGTTGCC-3', 'AACCATGTCTGGCGAC-3', 'ACTGCTCAGCTGATAA-3', 'AGAGTGGCATAGGATA-3', 'ATCGAGTTCCGGCACA-3', 'CAGCATAGTGAAGGCT-3', 'CAGCTGGTCACATGCA-3', 'CCGGGATCAGCTGCTG-3', 'GCTGCAGTCCGCGCAA-3', 'TCAACGAGTGTAAGTA-3', 'TCCACACGTGAACCTT-3', 'TCGTACCAGGGTTTCT-3', 'TCGTAGAAGGGCTTGA-3', 'TGGACGCGTACGACCC-3', 'TTCGAAGTCATTCACT-3', 'TCAGATGTCAAGGTAA-3'), invert = F)

Clone1 <- SetIdent(Clone1, value = "Genotype")
heatmap <- DoHeatmap(Clone1, 
                     features = markers.clone1top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g1_recip <- WhichCells(CD8.clean3, idents = c('ACTATCTTCTCTGCTG-2', 'AGTGAGGCAAGTTGTC-2', 'GAAACTCGTCCATGAT-2', 'GAACATCCATCAGTCA-2', 'GACGTTAGTCGCGAAA-2', 'GTGCTTCCATGCCTTC-2', 'AAAGTAGAGACTAGAT-3', 'ACTGTCCGTCCGTCAG-3', 'CAGGTGCAGAGGACGG-3', 'CCGTTCACAGTATAAG-3', 'TGAGAGGAGCGCCTTG-3', 'TTTGGTTAGCGTTGCC-3'))
g1_don <- WhichCells(CD8.clean3, idents = c('AACCATGTCTGGCGAC-3', 'ACTGCTCAGCTGATAA-3', 'AGAGTGGCATAGGATA-3', 'ATCGAGTTCCGGCACA-3', 'CAGCATAGTGAAGGCT-3', 'CAGCTGGTCACATGCA-3', 'CCGGGATCAGCTGCTG-3', 'GCTGCAGTCCGCGCAA-3', 'TCAACGAGTGTAAGTA-3', 'TCCACACGTGAACCTT-3', 'TCGTACCAGGGTTTCT-3', 'TCGTAGAAGGGCTTGA-3', 'TGGACGCGTACGACCC-3', 'TTCGAAGTCATTCACT-3', 'TCAGATGTCAAGGTAA-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g1_recip, g1_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")

DimPlot(CD8.clean3, reduction = "umap", pt.size = 0.1, label = F, group.by="CellID", cols = c("#e0e0e0"), 
        cells.highlight = c('ACTATCTTCTCTGCTG-2', 'AGTGAGGCAAGTTGTC-2', 'GAAACTCGTCCATGAT-2', 'GAACATCCATCAGTCA-2', 'GACGTTAGTCGCGAAA-2', 'GTGCTTCCATGCCTTC-2', 'AAAGTAGAGACTAGAT-3', 'ACTGTCCGTCCGTCAG-3', 'CAGGTGCAGAGGACGG-3', 'CCGTTCACAGTATAAG-3', 'TGAGAGGAGCGCCTTG-3', 'TTTGGTTAGCGTTGCC-3', 'AACCATGTCTGGCGAC-3', 'ACTGCTCAGCTGATAA-3', 'AGAGTGGCATAGGATA-3', 'ATCGAGTTCCGGCACA-3', 'CAGCATAGTGAAGGCT-3', 'CAGCTGGTCACATGCA-3', 'CCGGGATCAGCTGCTG-3', 'GCTGCAGTCCGCGCAA-3', 'TCAACGAGTGTAAGTA-3', 'TCCACACGTGAACCTT-3', 'TCGTACCAGGGTTTCT-3', 'TCGTAGAAGGGCTTGA-3', 'TGGACGCGTACGACCC-3', 'TTCGAAGTCATTCACT-3', 'TCAGATGTCAAGGTAA-3'),
        cols.highlight = c('#e41a1c'), sizes.highlight = 0.8)

DimPlot(CD8.clean3, reduction = "umap", pt.size = 0.1, label = F, 
        cells.highlight = "Don1", cols.highlight = c('#e41a1c'), sizes.highlight = 0.8)

#2
Recip2 =  c('CGCGGTAGTTAAAGTG-1','TGACTTTCATTTCACT-1', 'CAGAATCCACAGTCGC-2', 'CTCTACGGTCTCGTTC-2', 'TCATTTGCACGGACAA-2', 'TGCTGCTCATTCCTGC-2', 'AAAGATGCAAAGGAAG-3', 'GCAAACTAGTGGGTTG-3', 'GCCAAATAGAGTCTGG-3', 'GCGACCAAGCTTTGGT-3', 'TCTTCGGAGAAAGTGG-3', 'TGCGCAGTCTTTACGT-3')
Don2 =  c('AACCATGTCTGGCGAC-3', 'ACTGCTCAGCTGATAA-3', 'AGAGTGGCATAGGATA-3', 'ATCGAGTTCCGGCACA-3', 'CAGCATAGTGAAGGCT-3', 'CAGCTGGTCACATGCA-3', 'CCGGGATCAGCTGCTG-3', 'GCTGCAGTCCGCGCAA-3',  'TCAACGAGTGTAAGTA-3', 'TCCACACGTGAACCTT-3', 'TCGTACCAGGGTTTCT-3', 'TCGTAGAAGGGCTTGA-3', 'TGGACGCGTACGACCC-3', 'TTCGAAGTCATTCACT-3', 'TCAGATGTCAAGGTAA-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
markers.clone2 <- FindAllMarkers(Clone2,
                                 features = vargenes.filtered5,
                                 only.pos = FALSE, 
                                 min.pct = 0.01, 
                                 logfc.threshold = 0.25,
                                 test.use = "MAST",
                                 verbose = T,
                                 latent.vars = "nCount_RNA")

setDT(markers.clone2, keep.rownames = T)[]
colnames(markers.clone2) <- c('gene', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj')


markers.clone2top20 <- markers.clone2 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


Clone2 <- subset(CD8.clean3, ident = c('CGCGGTAGTTAAAGTG-1','TGACTTTCATTTCACT-1', 'CAGAATCCACAGTCGC-2', 'CTCTACGGTCTCGTTC-2', 'TCATTTGCACGGACAA-2', 'TGCTGCTCATTCCTGC-2', 'AAAGATGCAAAGGAAG-3', 'GCAAACTAGTGGGTTG-3', 'GCCAAATAGAGTCTGG-3', 'GCGACCAAGCTTTGGT-3', 'TCTTCGGAGAAAGTGG-3', 'TGCGCAGTCTTTACGT-3', 'AACCATGTCTGGCGAC-3', 'ACTGCTCAGCTGATAA-3', 'AGAGTGGCATAGGATA-3', 'ATCGAGTTCCGGCACA-3', 'CAGCATAGTGAAGGCT-3', 'CAGCTGGTCACATGCA-3', 'CCGGGATCAGCTGCTG-3', 'GCTGCAGTCCGCGCAA-3',  'TCAACGAGTGTAAGTA-3', 'TCCACACGTGAACCTT-3', 'TCGTACCAGGGTTTCT-3', 'TCGTAGAAGGGCTTGA-3', 'TGGACGCGTACGACCC-3', 'TTCGAAGTCATTCACT-3', 'TCAGATGTCAAGGTAA-3'), invert = F)

Clone2 <- SetIdent(Clone2, value = "Genotype")
heatmap <- DoHeatmap(Clone2, 
                     features = markers.clone2top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g2_recip <- WhichCells(CD8.clean3, idents = c('CGCGGTAGTTAAAGTG-1','TGACTTTCATTTCACT-1', 'CAGAATCCACAGTCGC-2', 'CTCTACGGTCTCGTTC-2', 'TCATTTGCACGGACAA-2', 'TGCTGCTCATTCCTGC-2', 'AAAGATGCAAAGGAAG-3', 'GCAAACTAGTGGGTTG-3', 'GCCAAATAGAGTCTGG-3', 'GCGACCAAGCTTTGGT-3', 'TCTTCGGAGAAAGTGG-3', 'TGCGCAGTCTTTACGT-3'))
g2_don <- WhichCells(CD8.clean3, idents = c('AACCATGTCTGGCGAC-3', 'ACTGCTCAGCTGATAA-3', 'AGAGTGGCATAGGATA-3', 'ATCGAGTTCCGGCACA-3', 'CAGCATAGTGAAGGCT-3', 'CAGCTGGTCACATGCA-3', 'CCGGGATCAGCTGCTG-3', 'GCTGCAGTCCGCGCAA-3',  'TCAACGAGTGTAAGTA-3', 'TCCACACGTGAACCTT-3', 'TCGTACCAGGGTTTCT-3', 'TCGTAGAAGGGCTTGA-3', 'TGGACGCGTACGACCC-3', 'TTCGAAGTCATTCACT-3', 'TCAGATGTCAAGGTAA-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g2_recip, g2_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")

#3
Recip3 =  c('ATTTCTGAGGAATGGA-1', 'CGACCTTAGGAATCGC-1', 'GGAACTTGTAGAGTGC-1', 'GGGCACTTCATCTGTT-1', 'GTATCTTAGCCCAGCT-1', 'AGAATAGAGTCGTACT-2', 'CACATAGGTGCAACGA-2', 'CATGCCTTCGGCGCAT-2', 'GAATAAGAGAGCTGCA-2', 'GACGGCTAGGAGCGTT-2', 'TCGTAGAAGTCCGGTC-2', 'TGCCCTAAGGATGGTC-2', 'TGCTACCGTTTCCACC-2', 'TTAGGCACAGTGACAG-2', 'TTTATGCGTTACGTCA-2', 'TTAGGACTCAAACCAC-2', 'ATAACGCTCTCGCTTG-3', 'GGCTGGTGTACTCTCC-3', 'GGGAGATGTAGCTGCC-3', 'TAGCCGGGTTTGGCGC-3', 'TCGCGTTGTGCAACTT-3', 'TGAGAGGCAAGCTGGA-3')
Don3 =  c('AAAGCAACATCGGGTC-3', 'CGTGAGCTCCAGTATG-3', 'GAGTCCGCAAGGTTCT-3', 'GTATTCTCAGCTGTTA-3', 'GTCACAAGTGTAAGTA-3', 'TCAGGTAGTCATACTG-3', 'TGTGGTAGTATCTGCA-3', 'TTCTCAAGTTCCACGG-3', 'GTGTTAGGTATATGAG-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone3 <- subset(CD8.clean3, ident = c('ATTTCTGAGGAATGGA-1', 'CGACCTTAGGAATCGC-1', 'GGAACTTGTAGAGTGC-1', 'GGGCACTTCATCTGTT-1', 'GTATCTTAGCCCAGCT-1', 'AGAATAGAGTCGTACT-2', 'CACATAGGTGCAACGA-2', 'CATGCCTTCGGCGCAT-2', 'GAATAAGAGAGCTGCA-2', 'GACGGCTAGGAGCGTT-2', 'TCGTAGAAGTCCGGTC-2', 'TGCCCTAAGGATGGTC-2', 'TGCTACCGTTTCCACC-2', 'TTAGGCACAGTGACAG-2', 'TTTATGCGTTACGTCA-2', 'TTAGGACTCAAACCAC-2', 'ATAACGCTCTCGCTTG-3', 'GGCTGGTGTACTCTCC-3', 'GGGAGATGTAGCTGCC-3', 'TAGCCGGGTTTGGCGC-3', 'TCGCGTTGTGCAACTT-3', 'TGAGAGGCAAGCTGGA-3', 'AAAGCAACATCGGGTC-3', 'CGTGAGCTCCAGTATG-3', 'GAGTCCGCAAGGTTCT-3', 'GTATTCTCAGCTGTTA-3', 'GTCACAAGTGTAAGTA-3', 'TCAGGTAGTCATACTG-3', 'TGTGGTAGTATCTGCA-3', 'TTCTCAAGTTCCACGG-3', 'GTGTTAGGTATATGAG-3'), invert = F)

Clone3 <- SetIdent(Clone3, value = "Genotype")
markers.clone3 <- FindAllMarkers(Clone3,
                                 features = vargenes.filtered5,
                                 only.pos = FALSE, 
                                 min.pct = 0.01, 
                                 logfc.threshold = 0.25,
                                 test.use = "MAST",
                                 verbose = T,
                                 latent.vars = "nCount_RNA")


markers.clone3top20 <- markers.clone3 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone3, 
                     features = markers.clone3top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g3_recip <- WhichCells(CD8.clean3, idents = c('ATTTCTGAGGAATGGA-1', 'CGACCTTAGGAATCGC-1', 'GGAACTTGTAGAGTGC-1', 'GGGCACTTCATCTGTT-1', 'GTATCTTAGCCCAGCT-1', 'AGAATAGAGTCGTACT-2', 'CACATAGGTGCAACGA-2', 'CATGCCTTCGGCGCAT-2', 'GAATAAGAGAGCTGCA-2', 'GACGGCTAGGAGCGTT-2', 'TCGTAGAAGTCCGGTC-2', 'TGCCCTAAGGATGGTC-2', 'TGCTACCGTTTCCACC-2', 'TTAGGCACAGTGACAG-2', 'TTTATGCGTTACGTCA-2', 'TTAGGACTCAAACCAC-2', 'ATAACGCTCTCGCTTG-3', 'GGCTGGTGTACTCTCC-3', 'GGGAGATGTAGCTGCC-3', 'TAGCCGGGTTTGGCGC-3', 'TCGCGTTGTGCAACTT-3', 'TGAGAGGCAAGCTGGA-3'))
g3_don <- WhichCells(CD8.clean3, idents = c('AAAGCAACATCGGGTC-3', 'CGTGAGCTCCAGTATG-3', 'GAGTCCGCAAGGTTCT-3', 'GTATTCTCAGCTGTTA-3', 'GTCACAAGTGTAAGTA-3', 'TCAGGTAGTCATACTG-3', 'TGTGGTAGTATCTGCA-3', 'TTCTCAAGTTCCACGG-3', 'GTGTTAGGTATATGAG-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g3_recip, g3_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")


#4
Recip4 =  c('AAAGATGTCAGATAAG-2', 'AATCCAGCATCAGTAC-2', 'ACGCCAGCACATCCGG-2', 'AGACGTTGTGAGGGAG-2', 'AGGGTGATCAACCATG-2', 'ATTTCTGCAATGAATG-2', 'CACACAAAGTAACCCT-2', 'CATCCACTCTTGGGTA-2', 'CCACCTAGTAACGTTC-2', 'CGCGGTAGTCTAGTCA-2', 'CGGACACGTTCCCTTG-2', 'CTGAAGTTCGTGACAT-2', 'GACGCGTAGTGATCGG-2', 'GACTAACCATGAACCT-2', 'GCACATATCGTATCAG-2', 'GCAGTTATCAGCTTAG-2', 'GTGCAGCAGCTGTTCA-2', 'TACACGACACAAGACG-2', 'TACAGTGCAAGGTTTC-2', 'TCGGGACCAGGATTGG-2', 'TCTGAGAAGCTAGTCT-2', 'TGAGCATTCGTCTGAA-2', 'CTACACCTCTGGCGTG-2', 'ACGGCCATCCAGAAGG-3', 'AGCCTAAAGGGTTCCC-3', 'AGGTCATCAAAGGTGC-3', 'AGTAGTCAGAGTCTGG-3', 'ATAACGCGTGTTTGTG-3', 'CGAACATAGGATGCGT-3', 'CGGAGCTTCTACTATC-3', 'CTCGAAATCCACGACG-3', 'GAATGAAGTAAGTGTA-3', 'GACGGCTTCTTTAGGG-3', 'GATCTAGTCAGGTTCA-3', 'GCGCCAAAGAGTCTGG-3', 'GGCGACTAGAGAACAG-3', 'TAAACCGGTGTTTGTG-3', 'TCATTACAGAGTAAGG-3', 'TCTTCGGAGAATCTCC-3', 'TGCCAAAGTATAATGG-3', 'TTCTCCTGTTGACGTT-3', 'GTAACGTAGGGATGGG-3')
Don4 =  c('CGCTTCAGTCTAGCGC-2', 'GCGCCAAAGGCATGGT-2', 'TACTCGCCAGATGGGT-2', 'TGGTTAGAGAGTTGGC-2', 'AGATCTGTCGCTGATA-3', 'CGAGCACGTGGCGAAT-3', 'CTCGGGACAGGCAGTA-3', 'CTGCCTACAAGCGTAG-3', 'TGACTTTGTTCACCTC-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone4 <- subset(CD8.clean3, ident = c('AAAGATGTCAGATAAG-2', 'AATCCAGCATCAGTAC-2', 'ACGCCAGCACATCCGG-2', 'AGACGTTGTGAGGGAG-2', 'AGGGTGATCAACCATG-2', 'ATTTCTGCAATGAATG-2', 'CACACAAAGTAACCCT-2', 'CATCCACTCTTGGGTA-2', 'CCACCTAGTAACGTTC-2', 'CGCGGTAGTCTAGTCA-2', 'CGGACACGTTCCCTTG-2', 'CTGAAGTTCGTGACAT-2', 'GACGCGTAGTGATCGG-2', 'GACTAACCATGAACCT-2', 'GCACATATCGTATCAG-2', 'GCAGTTATCAGCTTAG-2', 'GTGCAGCAGCTGTTCA-2', 'TACACGACACAAGACG-2', 'TACAGTGCAAGGTTTC-2', 'TCGGGACCAGGATTGG-2', 'TCTGAGAAGCTAGTCT-2', 'TGAGCATTCGTCTGAA-2', 'CTACACCTCTGGCGTG-2', 'ACGGCCATCCAGAAGG-3', 'AGCCTAAAGGGTTCCC-3', 'AGGTCATCAAAGGTGC-3', 'AGTAGTCAGAGTCTGG-3', 'ATAACGCGTGTTTGTG-3', 'CGAACATAGGATGCGT-3', 'CGGAGCTTCTACTATC-3', 'CTCGAAATCCACGACG-3', 'GAATGAAGTAAGTGTA-3', 'GACGGCTTCTTTAGGG-3', 'GATCTAGTCAGGTTCA-3', 'GCGCCAAAGAGTCTGG-3', 'GGCGACTAGAGAACAG-3', 'TAAACCGGTGTTTGTG-3', 'TCATTACAGAGTAAGG-3', 'TCTTCGGAGAATCTCC-3', 'TGCCAAAGTATAATGG-3', 'TTCTCCTGTTGACGTT-3', 'GTAACGTAGGGATGGG-3', 'CGCTTCAGTCTAGCGC-2', 'GCGCCAAAGGCATGGT-2', 'TACTCGCCAGATGGGT-2', 'TGGTTAGAGAGTTGGC-2', 'AGATCTGTCGCTGATA-3', 'CGAGCACGTGGCGAAT-3', 'CTCGGGACAGGCAGTA-3', 'CTGCCTACAAGCGTAG-3', 'TGACTTTGTTCACCTC-3'), invert = F)

Clone4 <- SetIdent(Clone4, value = "Genotype")
markers.clone4 <- FindAllMarkers(Clone4,
                                 features = vargenes.filtered5,
                                 only.pos = FALSE, 
                                 min.pct = 0.01, 
                                 logfc.threshold = 0.25,
                                 test.use = "MAST",
                                 verbose = T,
                                 latent.vars = "nCount_RNA")


markers.clone4top20 <- markers.clone4 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone4, 
                     features = markers.clone4top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g4_recip <- WhichCells(CD8.clean3, idents = c('AAAGATGTCAGATAAG-2', 'AATCCAGCATCAGTAC-2', 'ACGCCAGCACATCCGG-2', 'AGACGTTGTGAGGGAG-2', 'AGGGTGATCAACCATG-2', 'ATTTCTGCAATGAATG-2', 'CACACAAAGTAACCCT-2', 'CATCCACTCTTGGGTA-2', 'CCACCTAGTAACGTTC-2', 'CGCGGTAGTCTAGTCA-2', 'CGGACACGTTCCCTTG-2', 'CTGAAGTTCGTGACAT-2', 'GACGCGTAGTGATCGG-2', 'GACTAACCATGAACCT-2', 'GCACATATCGTATCAG-2', 'GCAGTTATCAGCTTAG-2', 'GTGCAGCAGCTGTTCA-2', 'TACACGACACAAGACG-2', 'TACAGTGCAAGGTTTC-2', 'TCGGGACCAGGATTGG-2', 'TCTGAGAAGCTAGTCT-2', 'TGAGCATTCGTCTGAA-2', 'CTACACCTCTGGCGTG-2', 'ACGGCCATCCAGAAGG-3', 'AGCCTAAAGGGTTCCC-3', 'AGGTCATCAAAGGTGC-3', 'AGTAGTCAGAGTCTGG-3', 'ATAACGCGTGTTTGTG-3', 'CGAACATAGGATGCGT-3', 'CGGAGCTTCTACTATC-3', 'CTCGAAATCCACGACG-3', 'GAATGAAGTAAGTGTA-3', 'GACGGCTTCTTTAGGG-3', 'GATCTAGTCAGGTTCA-3', 'GCGCCAAAGAGTCTGG-3', 'GGCGACTAGAGAACAG-3', 'TAAACCGGTGTTTGTG-3', 'TCATTACAGAGTAAGG-3', 'TCTTCGGAGAATCTCC-3', 'TGCCAAAGTATAATGG-3', 'TTCTCCTGTTGACGTT-3', 'GTAACGTAGGGATGGG-3'))
g4_don <- WhichCells(CD8.clean3, idents = c('CGCTTCAGTCTAGCGC-2', 'GCGCCAAAGGCATGGT-2', 'TACTCGCCAGATGGGT-2', 'TGGTTAGAGAGTTGGC-2', 'AGATCTGTCGCTGATA-3', 'CGAGCACGTGGCGAAT-3', 'CTCGGGACAGGCAGTA-3', 'CTGCCTACAAGCGTAG-3', 'TGACTTTGTTCACCTC-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g4_recip, g4_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")



#5
Recip5 =  c('AAAGATGTCAGATAAG-2', 'AATCCAGCATCAGTAC-2', 'ACGCCAGCACATCCGG-2', 'AGACGTTGTGAGGGAG-2', 'AGGGTGATCAACCATG-2', 'ATTTCTGCAATGAATG-2', 'CACACAAAGTAACCCT-2', 'CATCCACTCTTGGGTA-2', 'CCACCTAGTAACGTTC-2', 'CGCGGTAGTCTAGTCA-2', 'CGGACACGTTCCCTTG-2', 'CTGAAGTTCGTGACAT-2', 'GACGCGTAGTGATCGG-2', 'GACTAACCATGAACCT-2', 'GCACATATCGTATCAG-2', 'GCAGTTATCAGCTTAG-2', 'GTGCAGCAGCTGTTCA-2', 'TACACGACACAAGACG-2', 'TACAGTGCAAGGTTTC-2', 'TCGGGACCAGGATTGG-2', 'TCTGAGAAGCTAGTCT-2', 'TGAGCATTCGTCTGAA-2', 'CTACACCTCTGGCGTG-2', 'ACGGCCATCCAGAAGG-3', 'AGCCTAAAGGGTTCCC-3', 'AGGTCATCAAAGGTGC-3', 'AGTAGTCAGAGTCTGG-3', 'ATAACGCGTGTTTGTG-3', 'CGAACATAGGATGCGT-3', 'CGGAGCTTCTACTATC-3', 'CTCGAAATCCACGACG-3', 'GAATGAAGTAAGTGTA-3', 'GACGGCTTCTTTAGGG-3', 'GATCTAGTCAGGTTCA-3', 'GCGCCAAAGAGTCTGG-3', 'GGCGACTAGAGAACAG-3', 'TAAACCGGTGTTTGTG-3', 'TCATTACAGAGTAAGG-3', 'TCTTCGGAGAATCTCC-3', 'TGCCAAAGTATAATGG-3', 'TTCTCCTGTTGACGTT-3', 'GTAACGTAGGGATGGG-3')
Don5 =  c('AACACGTTCACTATTC-2', 'AATCCAGGTCGCGAAA-2', 'ACTGATGTCGCATGAT-2', 'ACTGTCCGTGACTCAT-2', 'ACTTGTTGTCATTAGC-2', 'AGCATACAGATCTGCT-2', 'AGCATACCACCGGAAA-2', 'ATGGGAGCATGGTTGT-2', 'CAGATCACAGCTGCAC-2', 'CAGGTGCAGCGACGTA-2', 'CATCAAGGTGAGCGAT-2', 'CCAATCCCAGCTGTTA-2', 'CCACTACGTTGGGACA-2', 'CTCGTACCACTAGTAC-2', 'CTCGTCATCGCGGATC-2', 'CTGTTTATCCTATTCA-2', 'CTTCTCTTCCCGACTT-2', 'GATCAGTCACACGCTG-2', 'GCATGCGAGTGGACGT-2', 'GCCTCTAAGTTAGGTA-2', 'GCCTCTATCGGAAACG-2', 'GCGCAGTAGATGAGAG-2', 'GGACAGATCACGCGGT-2', 'GGACATTGTTTGGGCC-2', 'GTTCATTAGACGCAAC-2', 'TACCTATTCGCGTAGC-2', 'TATCTCAGTACATGTC-2', 'TGACAACAGGACCACA-2', 'TGACTTTGTACGAAAT-2', 'TGACTTTGTCTCACCT-2', 'TGGCTGGTCCTAGGGC-2', 'CTCTGGTGTCCAGTAT-2', 'AACTCCCGTCTTGCGG-3', 'AACTTTCCAGGTGGAT-3', 'AAGGTTCGTGAGGCTA-3', 'ACACTGATCAAACAAG-3', 'ACCCACTAGGATGGTC-3', 'AGCGTATAGTGCGTGA-3', 'ATTTCTGCATAGGATA-3', 'CAACCAACAGATCCAT-3', 'CACAGGCGTAATAGCA-3', 'CACTCCAGTATAGGTA-3', 'CAGAGAGCAGTCACTA-3', 'CATTATCAGCCCTAAT-3', 'CCTTACGGTGGTCTCG-3', 'CCTTACGTCTACTTAC-3', 'CGGAGTCTCAGGCAAG-3', 'CGTAGCGTCGAATCCA-3', 'CGTGTAACAGTAAGCG-3', 'GAAATGAGTAAATGAC-3', 'GACCTGGCAAGCCCAC-3', 'GACGCGTCATCGGAAG-3', 'GGGAATGAGTTTCCTT-3', 'TAAGTGCTCGGCGGTT-3', 'TACACGAGTTTGTTTC-3', 'TACCTTACAAGCCGTC-3', 'TACCTTACACCAACCG-3', 'TACTCATTCCTGTAGA-3', 'TACTCGCGTTAAGACA-3', 'TCGCGTTAGACTAAGT-3', 'TCTGGAAGTGTAACGG-3', 'TGAGGGACATACTCTT-3', 'TTGACTTGTGAGGCTA-3', 'ATCATGGCATCCTTGC-3', 'GCCAAATCAACTTGAC-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone5 <- subset(CD8.clean3, ident = c('AAAGATGTCAGATAAG-2', 'AATCCAGCATCAGTAC-2', 'ACGCCAGCACATCCGG-2', 'AGACGTTGTGAGGGAG-2', 'AGGGTGATCAACCATG-2', 'ATTTCTGCAATGAATG-2', 'CACACAAAGTAACCCT-2', 'CATCCACTCTTGGGTA-2', 'CCACCTAGTAACGTTC-2', 'CGCGGTAGTCTAGTCA-2', 'CGGACACGTTCCCTTG-2', 'CTGAAGTTCGTGACAT-2', 'GACGCGTAGTGATCGG-2', 'GACTAACCATGAACCT-2', 'GCACATATCGTATCAG-2', 'GCAGTTATCAGCTTAG-2', 'GTGCAGCAGCTGTTCA-2', 'TACACGACACAAGACG-2', 'TACAGTGCAAGGTTTC-2', 'TCGGGACCAGGATTGG-2', 'TCTGAGAAGCTAGTCT-2', 'TGAGCATTCGTCTGAA-2', 'CTACACCTCTGGCGTG-2', 'ACGGCCATCCAGAAGG-3', 'AGCCTAAAGGGTTCCC-3', 'AGGTCATCAAAGGTGC-3', 'AGTAGTCAGAGTCTGG-3', 'ATAACGCGTGTTTGTG-3', 'CGAACATAGGATGCGT-3', 'CGGAGCTTCTACTATC-3', 'CTCGAAATCCACGACG-3', 'GAATGAAGTAAGTGTA-3', 'GACGGCTTCTTTAGGG-3', 'GATCTAGTCAGGTTCA-3', 'GCGCCAAAGAGTCTGG-3', 'GGCGACTAGAGAACAG-3', 'TAAACCGGTGTTTGTG-3', 'TCATTACAGAGTAAGG-3', 'TCTTCGGAGAATCTCC-3', 'TGCCAAAGTATAATGG-3', 'TTCTCCTGTTGACGTT-3', 'GTAACGTAGGGATGGG-3', 'AACACGTTCACTATTC-2', 'AATCCAGGTCGCGAAA-2', 'ACTGATGTCGCATGAT-2', 'ACTGTCCGTGACTCAT-2', 'ACTTGTTGTCATTAGC-2', 'AGCATACAGATCTGCT-2', 'AGCATACCACCGGAAA-2', 'ATGGGAGCATGGTTGT-2', 'CAGATCACAGCTGCAC-2', 'CAGGTGCAGCGACGTA-2', 'CATCAAGGTGAGCGAT-2', 'CCAATCCCAGCTGTTA-2', 'CCACTACGTTGGGACA-2', 'CTCGTACCACTAGTAC-2', 'CTCGTCATCGCGGATC-2', 'CTGTTTATCCTATTCA-2', 'CTTCTCTTCCCGACTT-2', 'GATCAGTCACACGCTG-2', 'GCATGCGAGTGGACGT-2', 'GCCTCTAAGTTAGGTA-2', 'GCCTCTATCGGAAACG-2', 'GCGCAGTAGATGAGAG-2', 'GGACAGATCACGCGGT-2', 'GGACATTGTTTGGGCC-2', 'GTTCATTAGACGCAAC-2', 'TACCTATTCGCGTAGC-2', 'TATCTCAGTACATGTC-2', 'TGACAACAGGACCACA-2', 'TGACTTTGTACGAAAT-2', 'TGACTTTGTCTCACCT-2', 'TGGCTGGTCCTAGGGC-2', 'CTCTGGTGTCCAGTAT-2', 'AACTCCCGTCTTGCGG-3', 'AACTTTCCAGGTGGAT-3', 'AAGGTTCGTGAGGCTA-3', 'ACACTGATCAAACAAG-3', 'ACCCACTAGGATGGTC-3', 'AGCGTATAGTGCGTGA-3', 'ATTTCTGCATAGGATA-3', 'CAACCAACAGATCCAT-3', 'CACAGGCGTAATAGCA-3', 'CACTCCAGTATAGGTA-3', 'CAGAGAGCAGTCACTA-3', 'CATTATCAGCCCTAAT-3', 'CCTTACGGTGGTCTCG-3', 'CCTTACGTCTACTTAC-3', 'CGGAGTCTCAGGCAAG-3', 'CGTAGCGTCGAATCCA-3', 'CGTGTAACAGTAAGCG-3', 'GAAATGAGTAAATGAC-3', 'GACCTGGCAAGCCCAC-3', 'GACGCGTCATCGGAAG-3', 'GGGAATGAGTTTCCTT-3', 'TAAGTGCTCGGCGGTT-3', 'TACACGAGTTTGTTTC-3', 'TACCTTACAAGCCGTC-3', 'TACCTTACACCAACCG-3', 'TACTCATTCCTGTAGA-3', 'TACTCGCGTTAAGACA-3', 'TCGCGTTAGACTAAGT-3', 'TCTGGAAGTGTAACGG-3', 'TGAGGGACATACTCTT-3', 'TTGACTTGTGAGGCTA-3', 'ATCATGGCATCCTTGC-3', 'GCCAAATCAACTTGAC-3'), invert = F)

Clone5 <- SetIdent(Clone5, value = "Genotype")
markers.clone5 <- FindAllMarkers(Clone5,
                                 features = vargenes.filtered5,
                                 only.pos = FALSE, 
                                 min.pct = 0.01, 
                                 logfc.threshold = 0.25,
                                 test.use = "MAST",
                                 verbose = T,
                                 latent.vars = "nCount_RNA")


markers.clone5top20 <- markers.clone5 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone5, 
                     features = markers.clone5top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g5_recip <- WhichCells(CD8.clean3, idents = c('AAAGATGTCAGATAAG-2', 'AATCCAGCATCAGTAC-2', 'ACGCCAGCACATCCGG-2', 'AGACGTTGTGAGGGAG-2', 'AGGGTGATCAACCATG-2', 'ATTTCTGCAATGAATG-2', 'CACACAAAGTAACCCT-2', 'CATCCACTCTTGGGTA-2', 'CCACCTAGTAACGTTC-2', 'CGCGGTAGTCTAGTCA-2', 'CGGACACGTTCCCTTG-2', 'CTGAAGTTCGTGACAT-2', 'GACGCGTAGTGATCGG-2', 'GACTAACCATGAACCT-2', 'GCACATATCGTATCAG-2', 'GCAGTTATCAGCTTAG-2', 'GTGCAGCAGCTGTTCA-2', 'TACACGACACAAGACG-2', 'TACAGTGCAAGGTTTC-2', 'TCGGGACCAGGATTGG-2', 'TCTGAGAAGCTAGTCT-2', 'TGAGCATTCGTCTGAA-2', 'CTACACCTCTGGCGTG-2', 'ACGGCCATCCAGAAGG-3', 'AGCCTAAAGGGTTCCC-3', 'AGGTCATCAAAGGTGC-3', 'AGTAGTCAGAGTCTGG-3', 'ATAACGCGTGTTTGTG-3', 'CGAACATAGGATGCGT-3', 'CGGAGCTTCTACTATC-3', 'CTCGAAATCCACGACG-3', 'GAATGAAGTAAGTGTA-3', 'GACGGCTTCTTTAGGG-3', 'GATCTAGTCAGGTTCA-3', 'GCGCCAAAGAGTCTGG-3', 'GGCGACTAGAGAACAG-3', 'TAAACCGGTGTTTGTG-3', 'TCATTACAGAGTAAGG-3', 'TCTTCGGAGAATCTCC-3', 'TGCCAAAGTATAATGG-3', 'TTCTCCTGTTGACGTT-3', 'GTAACGTAGGGATGGG-3'))
g5_don <- WhichCells(CD8.clean3, idents = c('AACACGTTCACTATTC-2', 'AATCCAGGTCGCGAAA-2', 'ACTGATGTCGCATGAT-2', 'ACTGTCCGTGACTCAT-2', 'ACTTGTTGTCATTAGC-2', 'AGCATACAGATCTGCT-2', 'AGCATACCACCGGAAA-2', 'ATGGGAGCATGGTTGT-2', 'CAGATCACAGCTGCAC-2', 'CAGGTGCAGCGACGTA-2', 'CATCAAGGTGAGCGAT-2', 'CCAATCCCAGCTGTTA-2', 'CCACTACGTTGGGACA-2', 'CTCGTACCACTAGTAC-2', 'CTCGTCATCGCGGATC-2', 'CTGTTTATCCTATTCA-2', 'CTTCTCTTCCCGACTT-2', 'GATCAGTCACACGCTG-2', 'GCATGCGAGTGGACGT-2', 'GCCTCTAAGTTAGGTA-2', 'GCCTCTATCGGAAACG-2', 'GCGCAGTAGATGAGAG-2', 'GGACAGATCACGCGGT-2', 'GGACATTGTTTGGGCC-2', 'GTTCATTAGACGCAAC-2', 'TACCTATTCGCGTAGC-2', 'TATCTCAGTACATGTC-2', 'TGACAACAGGACCACA-2', 'TGACTTTGTACGAAAT-2', 'TGACTTTGTCTCACCT-2', 'TGGCTGGTCCTAGGGC-2', 'CTCTGGTGTCCAGTAT-2', 'AACTCCCGTCTTGCGG-3', 'AACTTTCCAGGTGGAT-3', 'AAGGTTCGTGAGGCTA-3', 'ACACTGATCAAACAAG-3', 'ACCCACTAGGATGGTC-3', 'AGCGTATAGTGCGTGA-3', 'ATTTCTGCATAGGATA-3', 'CAACCAACAGATCCAT-3', 'CACAGGCGTAATAGCA-3', 'CACTCCAGTATAGGTA-3', 'CAGAGAGCAGTCACTA-3', 'CATTATCAGCCCTAAT-3', 'CCTTACGGTGGTCTCG-3', 'CCTTACGTCTACTTAC-3', 'CGGAGTCTCAGGCAAG-3', 'CGTAGCGTCGAATCCA-3', 'CGTGTAACAGTAAGCG-3', 'GAAATGAGTAAATGAC-3', 'GACCTGGCAAGCCCAC-3', 'GACGCGTCATCGGAAG-3', 'GGGAATGAGTTTCCTT-3', 'TAAGTGCTCGGCGGTT-3', 'TACACGAGTTTGTTTC-3', 'TACCTTACAAGCCGTC-3', 'TACCTTACACCAACCG-3', 'TACTCATTCCTGTAGA-3', 'TACTCGCGTTAAGACA-3', 'TCGCGTTAGACTAAGT-3', 'TCTGGAAGTGTAACGG-3', 'TGAGGGACATACTCTT-3', 'TTGACTTGTGAGGCTA-3', 'ATCATGGCATCCTTGC-3', 'GCCAAATCAACTTGAC-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g5_recip, g5_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")


#6
Recip6 =  c('CCGGGATTCGACGGAA-2', 'GATCGTAGTAGCAAAT-2', 'AAGTCTGCAACTGGCC-3', 'ATCCACCGTAACGCGA-3', 'GCACATACAGCTGCTG-3', 'TAAGCGTAGAGTACCG-3', 'TGAGGGACAGTCGATT-3')
Don6 =  c('AAATGCCTCAGAGCTT-2', 'CGATGTACATCCGCGA-2', 'CGCTATCCAGCTTAAC-2', 'GAACGGATCACGATGT-2', 'GTACTTTTCAAACCAC-2', 'TACACGAAGGGTATCG-2', 
          'TACTTACAGTCAATAG-2', 'TACTTACGTAGCACGA-2', 'TGACTAGAGGTAGCTG-2', 'TGGTTCCTCTTAACCT-2', 'CGTCACTAGTCGAGTG-2', 'AAAGCAAGTCGCGGTT-3', 
          'AACCGCGGTGATAAAC-3', 'ACACCCTCATTGCGGC-3', 'ACACCGGGTTCACGGC-3', 'ACGCCAGCAGACGCTC-3', 'ACTGAACTCAAGGCTT-3', 'AGAGCGAGTGGAAAGA-3', 
          'AGATCTGTCCGCAGTG-3', 'AGCGGTCGTTATCCGA-3', 'AGCTTGAAGGCTCAGA-3', 'AGGGTGAGTTGGTTTG-3', 'AGGTCCGTCGTATCAG-3', 'ATCATCTGTCCGTCAG-3', 
          'ATCATGGCATGCCACG-3', 'CACATTTTCACTATTC-3', 'CAGAGAGCAGCCACCA-3', 'CAGCTAATCACATAGC-3', 'CAGTCCTGTAGGGTAC-3', 'CATCGGGAGACCTTTG-3', 
          'CCGGTAGAGCGGCTTC-3', 'CCGTACTTCTATCGCC-3', 'CCGTGGAGTACCTACA-3', 'CCTTACGCAAGCGTAG-3', 'CCTTTCTGTACTTGAC-3', 'CGAGCACGTAGCGTGA-3', 
          'CGAGCCAAGTATGACA-3', 'CGCGTTTGTGATAAGT-3', 'CGGAGTCCATGGTAGG-3', 'CGTAGCGAGCGGCTTC-3', 'CGTGAGCGTTCGTTGA-3', 'CGTGAGCTCTTGTTTG-3', 
          'CTAAGACGTTTGTTTC-3', 'CTACACCGTTAAGAAC-3', 'CTCACACTCTCCAGGG-3', 'CTCATTAGTGTTCGAT-3', 'CTCCTAGGTCGACTAT-3', 'CTCGTACTCTCGAGTA-3', 
          'CTCTAATAGGGCATGT-3', 'CTTAACTAGTGGGCTA-3', 'CTTACCGTCGGGAGTA-3', 'CTTAGGATCCGCGGTA-3', 'GAAACTCGTCTAGCCG-3', 'GAACGGACACAGCCCA-3', 
          'GACAGAGCAGCCACCA-3', 'GATGAGGGTAAGTAGT-3', 'GCAATCATCACGGTTA-3', 'GCAGTTAAGCGTAGTG-3', 'GCATGATGTGCGAAAC-3', 'GCGAGAAGTGGCAAAC-3', 
          'GCGCCAAAGAGGTAGA-3', 'GCGCCAACAGGTGGAT-3', 'GGAAAGCCATCACGTA-3', 'GGTATTGAGGCTAGGT-3', 'GTACGTATCGGCGCTA-3', 'GTAGGCCCACATTTCT-3', 
          'GTCGTAAGTGATGTGG-3', 'GTGTGCGAGACCTTTG-3', 'GTGTGCGTCCACGACG-3', 'GTTAAGCCAGCGATCC-3', 'GTTCATTAGTACTTGC-3', 'TAAGTGCGTTACTGAC-3', 
          'TACGGTAGTTGACGTT-3', 'TACTCGCCACAAGACG-3', 'TAGTGGTAGAAGGGTA-3', 'TATTACCTCAACGAAA-3', 'TCAACGAGTTCTGTTT-3', 'TCAATCTCAAAGGTGC-3', 
          'TCCCGATAGGAGTTGC-3', 'TCGCGAGGTTAAGAAC-3', 'TCTCATAAGTCTCAAC-3', 'TCTCATAGTTATGCGT-3', 'TCTGAGACATACCATG-3', 'TGCCCATCACGTCTCT-3', 
          'TGCGGGTGTGACGGTA-3', 'TGGACGCCATCCAACA-3', 'TGTGTTTAGTACGACG-3', 'TTAGGACTCGGTCCGA-3', 'TTATGCTGTTCGTCTC-3', 'TTGGCAAAGTCTCCTC-3', 
          'GTTCGGGGTAGGAGTC-3', 'TAGACCAGTAGAGCTG-3', 'CAGCGACTCATCGATG-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone6 <- subset(CD8.clean3, ident = c('CCGGGATTCGACGGAA-2', 'GATCGTAGTAGCAAAT-2', 'AAGTCTGCAACTGGCC-3', 'ATCCACCGTAACGCGA-3', 'GCACATACAGCTGCTG-3', 'TAAGCGTAGAGTACCG-3', 'TGAGGGACAGTCGATT-3', 'AAATGCCTCAGAGCTT-2', 'CGATGTACATCCGCGA-2', 'CGCTATCCAGCTTAAC-2', 'GAACGGATCACGATGT-2', 'GTACTTTTCAAACCAC-2', 'TACACGAAGGGTATCG-2', 
                                       'TACTTACAGTCAATAG-2', 'TACTTACGTAGCACGA-2', 'TGGTTCCTCTTAACCT-2', 'CGTCACTAGTCGAGTG-2', 'AAAGCAAGTCGCGGTT-3', 
                                       'AACCGCGGTGATAAAC-3', 'ACACCCTCATTGCGGC-3', 'ACACCGGGTTCACGGC-3', 'ACGCCAGCAGACGCTC-3', 'ACTGAACTCAAGGCTT-3', 'AGAGCGAGTGGAAAGA-3', 
                                       'AGATCTGTCCGCAGTG-3', 'AGCGGTCGTTATCCGA-3', 'AGCTTGAAGGCTCAGA-3', 'AGGGTGAGTTGGTTTG-3', 'AGGTCCGTCGTATCAG-3', 'ATCATCTGTCCGTCAG-3', 
                                       'ATCATGGCATGCCACG-3', 'CACATTTTCACTATTC-3', 'CAGAGAGCAGCCACCA-3', 'CAGCTAATCACATAGC-3', 'CAGTCCTGTAGGGTAC-3', 'CATCGGGAGACCTTTG-3', 
                                       'CCGGTAGAGCGGCTTC-3', 'CCGTACTTCTATCGCC-3', 'CCGTGGAGTACCTACA-3', 'CCTTACGCAAGCGTAG-3', 'CCTTTCTGTACTTGAC-3', 'CGAGCACGTAGCGTGA-3', 
                                       'CGAGCCAAGTATGACA-3', 'CGCGTTTGTGATAAGT-3', 'CGGAGTCCATGGTAGG-3', 'CGTAGCGAGCGGCTTC-3', 'CGTGAGCGTTCGTTGA-3', 'CGTGAGCTCTTGTTTG-3', 
                                       'CTAAGACGTTTGTTTC-3', 'CTACACCGTTAAGAAC-3', 'CTCACACTCTCCAGGG-3', 'CTCATTAGTGTTCGAT-3', 'CTCCTAGGTCGACTAT-3', 'CTCGTACTCTCGAGTA-3', 
                                       'CTCTAATAGGGCATGT-3', 'CTTAACTAGTGGGCTA-3', 'CTTACCGTCGGGAGTA-3', 'CTTAGGATCCGCGGTA-3', 'GAAACTCGTCTAGCCG-3', 'GAACGGACACAGCCCA-3', 
                                       'GACAGAGCAGCCACCA-3', 'GATGAGGGTAAGTAGT-3', 'GCAATCATCACGGTTA-3', 'GCAGTTAAGCGTAGTG-3', 'GCATGATGTGCGAAAC-3', 'GCGAGAAGTGGCAAAC-3', 
                                       'GCGCCAAAGAGGTAGA-3', 'GCGCCAACAGGTGGAT-3', 'GGAAAGCCATCACGTA-3', 'GGTATTGAGGCTAGGT-3', 'GTACGTATCGGCGCTA-3', 'GTAGGCCCACATTTCT-3', 
                                       'GTCGTAAGTGATGTGG-3', 'GTGTGCGAGACCTTTG-3', 'GTGTGCGTCCACGACG-3', 'GTTAAGCCAGCGATCC-3', 'GTTCATTAGTACTTGC-3', 'TAAGTGCGTTACTGAC-3', 
                                       'TACGGTAGTTGACGTT-3', 'TACTCGCCACAAGACG-3', 'TAGTGGTAGAAGGGTA-3', 'TATTACCTCAACGAAA-3', 'TCAACGAGTTCTGTTT-3', 'TCAATCTCAAAGGTGC-3', 
                                       'TCCCGATAGGAGTTGC-3', 'TCGCGAGGTTAAGAAC-3', 'TCTCATAAGTCTCAAC-3', 'TCTCATAGTTATGCGT-3', 'TCTGAGACATACCATG-3', 'TGCCCATCACGTCTCT-3', 
                                       'TGCGGGTGTGACGGTA-3', 'TGGACGCCATCCAACA-3', 'TGTGTTTAGTACGACG-3', 'TTAGGACTCGGTCCGA-3', 'TTATGCTGTTCGTCTC-3', 'TTGGCAAAGTCTCCTC-3', 
                                       'GTTCGGGGTAGGAGTC-3', 'TAGACCAGTAGAGCTG-3', 'CAGCGACTCATCGATG-3'), invert = F)

Clone6 <- SetIdent(Clone6, value = "Genotype")
markers.clone6 <- FindAllMarkers(Clone6,
                                 features = vargenes.filtered5,
                                 only.pos = FALSE, 
                                 min.pct = 0.01, 
                                 logfc.threshold = 0.25,
                                 test.use = "MAST",
                                 verbose = T,
                                 latent.vars = "nCount_RNA")


markers.clone6top20 <- markers.clone6 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone6, 
                     features = markers.clone6top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g6_recip <- WhichCells(CD8.clean3, idents = c('CCGGGATTCGACGGAA-2', 'GATCGTAGTAGCAAAT-2', 'AAGTCTGCAACTGGCC-3', 'ATCCACCGTAACGCGA-3', 'GCACATACAGCTGCTG-3', 'TAAGCGTAGAGTACCG-3', 'TGAGGGACAGTCGATT-3'))
g6_don <- WhichCells(CD8.clean3, idents = c('AAATGCCTCAGAGCTT-2', 'CGATGTACATCCGCGA-2', 'CGCTATCCAGCTTAAC-2', 'GAACGGATCACGATGT-2', 'GTACTTTTCAAACCAC-2', 'TACACGAAGGGTATCG-2', 
                                            'TACTTACAGTCAATAG-2', 'TACTTACGTAGCACGA-2', 'TGGTTCCTCTTAACCT-2', 'CGTCACTAGTCGAGTG-2', 'AAAGCAAGTCGCGGTT-3', 
                                            'AACCGCGGTGATAAAC-3', 'ACACCCTCATTGCGGC-3', 'ACACCGGGTTCACGGC-3', 'ACGCCAGCAGACGCTC-3', 'ACTGAACTCAAGGCTT-3', 'AGAGCGAGTGGAAAGA-3', 
                                            'AGATCTGTCCGCAGTG-3', 'AGCGGTCGTTATCCGA-3', 'AGCTTGAAGGCTCAGA-3', 'AGGGTGAGTTGGTTTG-3', 'AGGTCCGTCGTATCAG-3', 'ATCATCTGTCCGTCAG-3', 
                                            'ATCATGGCATGCCACG-3', 'CACATTTTCACTATTC-3', 'CAGAGAGCAGCCACCA-3', 'CAGCTAATCACATAGC-3', 'CAGTCCTGTAGGGTAC-3', 'CATCGGGAGACCTTTG-3', 
                                            'CCGGTAGAGCGGCTTC-3', 'CCGTACTTCTATCGCC-3', 'CCGTGGAGTACCTACA-3', 'CCTTACGCAAGCGTAG-3', 'CCTTTCTGTACTTGAC-3', 'CGAGCACGTAGCGTGA-3', 
                                            'CGAGCCAAGTATGACA-3', 'CGCGTTTGTGATAAGT-3', 'CGGAGTCCATGGTAGG-3', 'CGTAGCGAGCGGCTTC-3', 'CGTGAGCGTTCGTTGA-3', 'CGTGAGCTCTTGTTTG-3', 
                                            'CTAAGACGTTTGTTTC-3', 'CTACACCGTTAAGAAC-3', 'CTCACACTCTCCAGGG-3', 'CTCATTAGTGTTCGAT-3', 'CTCCTAGGTCGACTAT-3', 'CTCGTACTCTCGAGTA-3', 
                                            'CTCTAATAGGGCATGT-3', 'CTTAACTAGTGGGCTA-3', 'CTTACCGTCGGGAGTA-3', 'CTTAGGATCCGCGGTA-3', 'GAAACTCGTCTAGCCG-3', 'GAACGGACACAGCCCA-3', 
                                            'GACAGAGCAGCCACCA-3', 'GATGAGGGTAAGTAGT-3', 'GCAATCATCACGGTTA-3', 'GCAGTTAAGCGTAGTG-3', 'GCATGATGTGCGAAAC-3', 'GCGAGAAGTGGCAAAC-3', 
                                            'GCGCCAAAGAGGTAGA-3', 'GCGCCAACAGGTGGAT-3', 'GGAAAGCCATCACGTA-3', 'GGTATTGAGGCTAGGT-3', 'GTACGTATCGGCGCTA-3', 'GTAGGCCCACATTTCT-3', 
                                            'GTCGTAAGTGATGTGG-3', 'GTGTGCGAGACCTTTG-3', 'GTGTGCGTCCACGACG-3', 'GTTAAGCCAGCGATCC-3', 'GTTCATTAGTACTTGC-3', 'TAAGTGCGTTACTGAC-3', 
                                            'TACGGTAGTTGACGTT-3', 'TACTCGCCACAAGACG-3', 'TAGTGGTAGAAGGGTA-3', 'TATTACCTCAACGAAA-3', 'TCAACGAGTTCTGTTT-3', 'TCAATCTCAAAGGTGC-3', 
                                            'TCCCGATAGGAGTTGC-3', 'TCGCGAGGTTAAGAAC-3', 'TCTCATAAGTCTCAAC-3', 'TCTCATAGTTATGCGT-3', 'TCTGAGACATACCATG-3', 'TGCCCATCACGTCTCT-3', 
                                            'TGCGGGTGTGACGGTA-3', 'TGGACGCCATCCAACA-3', 'TGTGTTTAGTACGACG-3', 'TTAGGACTCGGTCCGA-3', 'TTATGCTGTTCGTCTC-3', 'TTGGCAAAGTCTCCTC-3', 
                                            'GTTCGGGGTAGGAGTC-3', 'TAGACCAGTAGAGCTG-3', 'CAGCGACTCATCGATG-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g6_recip, g6_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")



#7
Recip7 =  c('ACGGGTCGTTTGACAC-2', 'CTCAGAATCTGCTGTC-2', 'CTCGTCATCTCGTATT-2', 'GTCGTAAGTAATCGTC-2', 'TTGAACGTCTGGCGAC-2', 'GTGGGTCCAGGCTGAA-3')
Don7 =  c('AAGGCAGAGGTAGCCA-2', 'AGGTCATCACCAGGCT-2', 'CGCTGGATCGTTTATC-2', 'CTAATGGAGTGTTAGA-2', 'GAATGAATCAGGCAAG-2', 'GACTAACAGCACGCCT-2', 'GACTGCGCACTTCTGC-2', 'GGGCATCAGTGCGATG-2', 'CATGGCGAGTCTCCTC-2', 'AGCCTAAGTCGATTGT-3', 'ATTCTACTCCATTCTA-3', 'GAATGAATCTTGTATC-3', 'GACTAACAGCAATCTC-3', 'TTCTACAGTTTGACAC-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone7 <- subset(CD8.clean3, ident = c('CTCAGAATCTGCTGTC-2', 'CTCGTCATCTCGTATT-2', 'GTCGTAAGTAATCGTC-2', 'TTGAACGTCTGGCGAC-2', 'GTGGGTCCAGGCTGAA-3', 'AAGGCAGAGGTAGCCA-2', 'AGGTCATCACCAGGCT-2', 'CGCTGGATCGTTTATC-2', 'CTAATGGAGTGTTAGA-2', 'GAATGAATCAGGCAAG-2', 'GACTAACAGCACGCCT-2', 'GACTGCGCACTTCTGC-2', 'GGGCATCAGTGCGATG-2', 'CATGGCGAGTCTCCTC-2', 'AGCCTAAGTCGATTGT-3', 'ATTCTACTCCATTCTA-3', 'GAATGAATCTTGTATC-3', 'GACTAACAGCAATCTC-3', 'TTCTACAGTTTGACAC-3'), invert = F)

Clone7 <- SetIdent(Clone7, value = "Genotype")
markers.clone7 <- FindAllMarkers(Clone7,
                                 features = vargenes.filtered5,
                                 only.pos = FALSE, 
                                 min.pct = 0.01, 
                                 logfc.threshold = 0.25,
                                 test.use = "MAST",
                                 verbose = T,
                                 latent.vars = "nCount_RNA")


markers.clone7top20 <- markers.clone7 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone7, 
                     features = markers.clone7top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g7_recip <- WhichCells(CD8.clean3, idents = c('CTCAGAATCTGCTGTC-2', 'CTCGTCATCTCGTATT-2', 'GTCGTAAGTAATCGTC-2', 'TTGAACGTCTGGCGAC-2', 'GTGGGTCCAGGCTGAA-3'))
g7_don <- WhichCells(CD8.clean3, idents = c('AAGGCAGAGGTAGCCA-2', 'AGGTCATCACCAGGCT-2', 'CGCTGGATCGTTTATC-2', 'CTAATGGAGTGTTAGA-2', 'GAATGAATCAGGCAAG-2', 'GACTAACAGCACGCCT-2', 'GACTGCGCACTTCTGC-2', 'GGGCATCAGTGCGATG-2', 'CATGGCGAGTCTCCTC-2', 'AGCCTAAGTCGATTGT-3', 'ATTCTACTCCATTCTA-3', 'GAATGAATCTTGTATC-3', 'GACTAACAGCAATCTC-3', 'TTCTACAGTTTGACAC-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g7_recip, g7_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")



#8
Recip4 =  c('ACGGGTCGTTTGACAC-2', 'CTCAGAATCTGCTGTC-2', 'CTCGTCATCTCGTATT-2', 'GTCGTAAGTAATCGTC-2', 'TTGAACGTCTGGCGAC-2', 'GTGGGTCCAGGCTGAA-3')
Don4 =  c('CATCAGATCTGTTGAG-2', 'CCGTACTTCCTCTAGC-2', 'GACGTGCGTAACGTTC-2', 'GTGAAGGAGTGGAGAA-2', 'TAAGTGCAGGGAGTAA-2', 'GTACTTTGTCAAACTC-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone8 <- subset(CD8.clean3, ident = c('CTCAGAATCTGCTGTC-2', 'CTCGTCATCTCGTATT-2', 'GTCGTAAGTAATCGTC-2', 'TTGAACGTCTGGCGAC-2', 'GTGGGTCCAGGCTGAA-3', 'CATCAGATCTGTTGAG-2', 'CCGTACTTCCTCTAGC-2', 'GACGTGCGTAACGTTC-2', 'GTGAAGGAGTGGAGAA-2', 'TAAGTGCAGGGAGTAA-2', 'GTACTTTGTCAAACTC-3'), invert = F)

Clone8 <- SetIdent(Clone8, value = "Genotype")
markers.clone8 <- FindAllMarkers(Clone8,
                                 features = vargenes.filtered5,
                                 only.pos = FALSE, 
                                 min.pct = 0.01, 
                                 logfc.threshold = 0.25,
                                 test.use = "MAST",
                                 verbose = T,
                                 latent.vars = "nCount_RNA")


markers.clone8top20 <- markers.clone8 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone8, 
                     features = markers.clone8top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g8_recip <- WhichCells(CD8.clean3, idents = c('CTCAGAATCTGCTGTC-2', 'CTCGTCATCTCGTATT-2', 'GTCGTAAGTAATCGTC-2', 'TTGAACGTCTGGCGAC-2', 'GTGGGTCCAGGCTGAA-3'))
g8_don <- WhichCells(CD8.clean3, idents = c('CATCAGATCTGTTGAG-2', 'CCGTACTTCCTCTAGC-2', 'GACGTGCGTAACGTTC-2', 'GTGAAGGAGTGGAGAA-2', 'TAAGTGCAGGGAGTAA-2', 'GTACTTTGTCAAACTC-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g8_recip, g8_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")



#9
Recip4 =  c('ATTGGACGTCTCCATC-2', 'CGAGAAGAGATGTCGG-2', 'CGCGGTATCTCTGTCG-2', 'AAGGTTCTCAGGTAAA-3', 'ACCAGTATCCTTCAAT-3', 'ACCGTAACAATCAGAA-3', 'AGAGCGAGTCCAACTA-3', 'AGGCCACCAGTGGAGT-3', 'CCTCTGATCTCAAACG-3', 'CGTTCTGGTCAAGCGA-3', 'CTAAGACTCGTAGATC-3', 'CTCTAATGTGTTCGAT-3', 'GGCCGATCAGATGGCA-3', 'GGCGACTGTCATGCAT-3', 'TCACGAACATTGAGCT-3', 'TCAGGATAGTGTCTCA-3', 'TCTCTAACAGCTTCGG-3', 'TGAGCATCAATTCCTT-3', 'TGCCCTAAGGGCACTA-3', 'TTTGGTTTCCCTTGCA-3')
Don4 =  c('AACCGCGTCACTTCAT-3', 'AACTCTTGTCGGGTCT-3', 'AAGGCAGAGGTAGCCA-3', 'AAGGTTCGTCGCCATG-3', 'AAGGTTCTCCACGTGG-3', 'AATCCAGAGTCTCGGC-3', 
          'ACACCAATCCACTGGG-3', 'ACACCGGGTCTAGTCA-3', 'ACACTGAAGAGTCTGG-3', 'ACAGCTATCTCTGCTG-3', 'ACATACGAGATCACGG-3', 'ACCAGTAAGTTACCCA-3', 
          'ACCGTAAGTGGTTTCA-3', 'ACGAGGACAAGTACCT-3', 'ACGATGTAGCTGGAAC-3', 'ACGGGCTAGCGTTTAC-3', 'ACGTCAATCCTAAGTG-3', 
          'ACGTCAATCTGCTTGC-3', 'ACTGAGTAGAGTAATC-3', 'ACTGCTCTCTACTTAC-3', 'AGAATAGAGGTTACCT-3', 'AGAGCGAGTCTTGATG-3', 'AGATCTGTCTTCAACT-3', 
          'AGCGTCGTCGCTAGCG-3', 'AGCTTGAAGGTGCACA-3', 'AGCTTGAGTGAGCGAT-3', 'AGGCCACGTGGGTCAA-3', 'AGGCCACTCGGAGCAA-3', 'AGGCCGTGTCTGCAAT-3', 
          'AGGGAGTTCCACGCAG-3', 'AGGGATGAGCGATTCT-3', 'AGGTCATGTGTGCGTC-3', 'AGGTCATGTTCGCGAC-3', 'AGTCTTTTCAGGTAAA-3', 'ATAACGCGTAAGGATT-3', 
          'ATAAGAGAGCAGCGTA-3', 'ATAAGAGGTGCACGAA-3', 'ATCATCTGTGTCAATC-3', 'ATCCGAAGTAGAAAGG-3', 'ATCTGCCTCACCGGGT-3', 'ATGTGTGGTATAGGGC-3', 
          'ATTGGACAGCTCCCAG-3', 'ATTGGACCAGGGAGAG-3', 'ATTGGTGCATCGGGTC-3', 'ATTTCTGAGACAAGCC-3', 'CAACCAAAGACAGGCT-3', 'CAAGATCGTAAGTAGT-3', 
          'CAAGATCGTATAGGGC-3', 'CAAGTTGTCCTGCCAT-3', 'CAAGTTGTCGTACGGC-3', 'CACACCTCACAGAGGT-3', 'CACACTCCATGAAGTA-3', 'CACATAGTCCTAGGGC-3', 
          'CACATTTAGGCGCTCT-3', 'CACCTTGAGATAGCAT-3', 'CAGATCACAGTTCCCT-3', 'CAGATCATCACCAGGC-3', 'CAGCATATCATCGATG-3', 'CATCAAGAGTAGTGCG-3', 
          'CATCCACCAGGCTGAA-3', 'CATCCACGTAGCGCAA-3', 'CATCGGGCAGCTTAAC-3', 'CCAATCCAGTTAACGA-3', 'CCACCTAGTGGTTTCA-3', 'CCAGCGAAGTGCAAGC-3', 
          'CCAGCGAGTAGCGCAA-3', 'CCCAATCGTTGAGTTC-3', 'CCCTCCTAGTATCTCG-3', 'CCCTCCTTCCCATTAT-3', 'CCCTCCTTCGCGGATC-3', 'CCGGGATCAGGGTACA-3', 
          'CCGTGGAAGCTAGGCA-3', 'CCGTGGATCCATGCTC-3', 'CCTAAAGCACCATCCT-3', 'CCTACCACATCAGTCA-3', 'CGAACATTCAACGGCC-3', 'CGAGAAGGTGGTAACG-3', 
          'CGAGCACGTGCTCTTC-3', 'CGAGCACTCAGCTCGG-3', 'CGATTGAAGCGTAATA-3', 'CGCCAAGAGACCTAGG-3', 'CGCCAAGGTCTGCAAT-3', 'CGCGTTTTCAACCAAC-3', 
          'CGCGTTTTCGAACGGA-3', 'CGCTGGAAGTGATCGG-3', 'CGCTTCAAGACAAGCC-3', 'CGGACACGTCATCCCT-3', 'CGGGTCACAAGAAAGG-3', 'CGTCAGGCACCACCAG-3', 
          'CGTCCATTCTGCGTAA-3', 'CGTGAGCAGCGTTTAC-3', 'CGTGAGCAGTTGAGTA-3', 'CGTGTAACATGTCTCC-3', 'CGTGTAAGTCTGCGGT-3', 'CGTTCTGCAGACGTAG-3', 
          'CGTTCTGTCTGGAGCC-3', 'CTAATGGGTGATGATA-3', 'CTACGTCGTTCTCATT-3', 'CTAGTGAAGCAGATCG-3', 'CTAGTGAAGGAGTACC-3', 'CTCAGAACAAAGAATC-3', 
          'CTCGAAAAGCGATATA-3', 'CTCGTCAGTTGGTGGA-3', 'CTGAAACAGGATGCGT-3', 'CTGAAGTAGGGTTCCC-3', 'CTGCCTAAGACACTAA-3', 'CTGTTTATCCCTAACC-3', 
          'CTTAGGACAAAGGCGT-3', 'CTTCTCTTCCTTAATC-3', 'CTTGGCTCAGGATTGG-3', 'GAAACTCAGGTAGCTG-3', 'GAAACTCCAGATCGGA-3', 'GAAACTCGTACAGTGG-3', 
          'GAAACTCGTAGGGACT-3', 'GAAACTCTCCCTTGCA-3', 'GAACCTAGTACCGTAT-3', 'GAATAAGGTATAATGG-3', 'GAATAAGTCAACGGGA-3', 'GAATAAGTCGAGCCCA-3', 
          'GACAGAGGTTGGTTTG-3', 'GACGCGTTCGATCCCT-3', 'GACGGCTGTAGCGCTC-3', 'GAGCAGAAGAAACGAG-3', 'GATCTAGCAGTTTACG-3', 'GATGAGGCAGCTGTTA-3', 
          'GCAAACTTCCGTAGTA-3', 'GCATGTACAGGAATCG-3', 'GCGCAGTTCCAACCAA-3', 'GCGCCAACATATGCTG-3', 'GCGGGTTTCAACGAAA-3', 'GCTTCCAGTTCGTGAT-3', 
          'GCTTGAAGTCTAGAGG-3', 'GGAAAGCAGTTATCGC-3', 'GGAATAATCGTTGCCT-3', 'GGACAGAGTCGCTTTC-3', 'GGACAGATCAGGTTCA-3', 'GGACATTAGTAATCCC-3', 
          'GGACATTCACGGTAAG-3', 'GGATGTTTCCACGCAG-3', 'GGATTACCATAGGATA-3', 'GGATTACGTACTCGCG-3', 'GGCTGGTAGCCGGTAA-3', 'GGGAGATTCTTTAGGG-3', 
          'GGGATGAAGATCCTGT-3', 'GGGCACTAGCCGGTAA-3', 'GGGCACTAGTTACCCA-3', 'GGGCATCCAATCCAAC-3', 'GGGTTGCCATGGGAAC-3', 'GGTATTGAGCAACGGT-3', 
          'GGTATTGTCGGTGTTA-3', 'GTATTCTTCGTGGACC-3', 'GTCAAGTCAGTCGTGC-3', 'GTCAAGTCATTTGCCC-3', 'GTCACGGAGGCTAGGT-3', 'GTCATTTCAAGACACG-3', 
          'GTCATTTTCGTATCAG-3', 'GTCGGGTAGCCCAACC-3', 'GTGCAGCTCACTCTTA-3', 'GTGCAGCTCTTACCGC-3', 'GTGGGTCAGCACCGTC-3', 'GTGGGTCAGCGTGTCC-3', 
          'GTTACAGTCATTGCGA-3', 'GTTCATTGTGAGCGAT-3', 'GTTCTCGCATTAACCG-3', 'TAAACCGCAAGCGAGT-3', 'TAAGAGAGTTCTGTTT-3', 'TAAGCGTCATTCTTAC-3', 
          'TAAGTGCGTGCAACTT-3', 'TACACGATCTTTACGT-3', 'TACGGTATCCATTCTA-3', 'TACTCATGTCGCTTCT-3', 'TACTTACGTACTCAAC-3', 'TACTTACTCACCCTCA-3', 
          'TATGCCCTCACAGGCC-3', 'TCAATCTAGTGCTGCC-3', 'TCAGATGCAGCTTCGG-3', 'TCAGATGGTTCCACAA-3', 'TCAGCAAAGCCCAATT-3', 'TCAGGATTCTGGCGTG-3', 
          'TCAGGTACAATAGCAA-3', 'TCATTTGTCACCACCT-3', 'TCCCGATAGATGCCTT-3', 'TCCCGATTCGCCATAA-3', 'TCGGTAATCACCAGGC-3', 'TCTGAGACAAGCTGGA-3', 
          'TCTGGAAAGTCTCAAC-3', 'TGAGCCGGTTGAACTC-3', 'TGAGGGAGTGCAGACA-3', 'TGCCCATCATACAGCT-3', 'TGCCCTACACTTAACG-3', 'TGGACGCTCCCGACTT-3', 
          'TGGCTGGCACGAAATA-3', 'TGGTTAGCACCTCGTT-3', 'TGGTTCCTCCCTTGTG-3', 'TGGTTCCTCGCCAGCA-3', 'TGTATTCCACTAAGTC-3', 'TGTGTTTAGCCAGGAT-3', 
          'TGTGTTTGTGTGACCC-3', 'TTAACTCAGCTACCTA-3', 'TTCCCAGAGCTAACTC-3', 'TTCTACAGTCTAGTCA-3', 'TTGACTTAGACTTGAA-3', 'TTGACTTGTCATCGGC-3', 
          'TTGCCGTAGACAATAC-3', 'TTGCCGTTCTATCGCC-3', 'TTGCGTCTCGTGGTCG-3', 'TTGTAGGCAATGGACG-3', 'TTTACTGAGCAGGCTA-3', 'TTTATGCAGTGACTCT-3', 
          'TTTATGCGTTCGAATC-3', 'TTTGGTTAGAATTGTG-3', 'CGACCTTGTCTCAACA-3', 'CAAGAAATCATGCATG-3', 'ATTCTACGTAGGGACT-3', 'ATCACGATCGGACAAG-3', 
          'ACACTGAAGAATGTGT-3', 'AATCGGTCAATGACCT-3', 'TTCCCAGCAAACGCGA-3', 'TACCTTATCCTAGAAC-3', 'GATTCAGTCATCTGCC-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone9 <- subset(CD8.clean3, ident = c('ATTGGACGTCTCCATC-2', 'CGAGAAGAGATGTCGG-2', 'CGCGGTATCTCTGTCG-2', 'AAGGTTCTCAGGTAAA-3', 'ACCAGTATCCTTCAAT-3', 'ACCGTAACAATCAGAA-3', 'AGAGCGAGTCCAACTA-3', 'AGGCCACCAGTGGAGT-3', 'CCTCTGATCTCAAACG-3', 'CGTTCTGGTCAAGCGA-3', 'CTAAGACTCGTAGATC-3', 'CTCTAATGTGTTCGAT-3', 'GGCCGATCAGATGGCA-3', 'GGCGACTGTCATGCAT-3', 'TCACGAACATTGAGCT-3', 'TCAGGATAGTGTCTCA-3', 'TCTCTAACAGCTTCGG-3', 'TGAGCATCAATTCCTT-3', 'TGCCCTAAGGGCACTA-3', 'TTTGGTTTCCCTTGCA-3', 'AACCGCGTCACTTCAT-3', 'AACTCTTGTCGGGTCT-3', 'AAGGCAGAGGTAGCCA-3', 'AAGGTTCGTCGCCATG-3', 'AAGGTTCTCCACGTGG-3', 'AATCCAGAGTCTCGGC-3', 
                                       'ACACCAATCCACTGGG-3', 'ACACCGGGTCTAGTCA-3', 'ACACTGAAGAGTCTGG-3', 'ACAGCTATCTCTGCTG-3', 'ACATACGAGATCACGG-3', 'ACCAGTAAGTTACCCA-3', 
                                       'ACCGTAAGTGGTTTCA-3', 'ACGAGGACAAGTACCT-3', 'ACGATGTAGCTGGAAC-3', 'ACGGGCTAGCGTTTAC-3', 'ACGTCAATCCTAAGTG-3', 
                                       'ACGTCAATCTGCTTGC-3', 'ACTGAGTAGAGTAATC-3', 'ACTGCTCTCTACTTAC-3', 'AGAATAGAGGTTACCT-3', 'AGAGCGAGTCTTGATG-3', 'AGATCTGTCTTCAACT-3', 
                                       'AGCGTCGTCGCTAGCG-3', 'AGCTTGAAGGTGCACA-3', 'AGCTTGAGTGAGCGAT-3', 'AGGCCACGTGGGTCAA-3', 'AGGCCACTCGGAGCAA-3', 'AGGCCGTGTCTGCAAT-3', 
                                       'AGGGAGTTCCACGCAG-3', 'AGGGATGAGCGATTCT-3', 'AGGTCATGTGTGCGTC-3', 'AGGTCATGTTCGCGAC-3', 'AGTCTTTTCAGGTAAA-3', 'ATAACGCGTAAGGATT-3', 
                                       'ATAAGAGAGCAGCGTA-3', 'ATAAGAGGTGCACGAA-3', 'ATCATCTGTGTCAATC-3', 'ATCCGAAGTAGAAAGG-3', 'ATCTGCCTCACCGGGT-3', 'ATGTGTGGTATAGGGC-3', 
                                       'ATTGGACAGCTCCCAG-3', 'ATTGGACCAGGGAGAG-3', 'ATTGGTGCATCGGGTC-3', 'ATTTCTGAGACAAGCC-3', 'CAACCAAAGACAGGCT-3', 'CAAGATCGTAAGTAGT-3', 
                                       'CAAGATCGTATAGGGC-3', 'CAAGTTGTCCTGCCAT-3', 'CAAGTTGTCGTACGGC-3', 'CACACCTCACAGAGGT-3', 'CACACTCCATGAAGTA-3', 'CACATAGTCCTAGGGC-3', 
                                       'CACATTTAGGCGCTCT-3', 'CACCTTGAGATAGCAT-3', 'CAGATCACAGTTCCCT-3', 'CAGATCATCACCAGGC-3', 'CAGCATATCATCGATG-3', 'CATCAAGAGTAGTGCG-3', 
                                       'CATCCACCAGGCTGAA-3', 'CATCCACGTAGCGCAA-3', 'CATCGGGCAGCTTAAC-3', 'CCAATCCAGTTAACGA-3', 'CCACCTAGTGGTTTCA-3', 'CCAGCGAAGTGCAAGC-3', 
                                       'CCAGCGAGTAGCGCAA-3', 'CCCAATCGTTGAGTTC-3', 'CCCTCCTAGTATCTCG-3', 'CCCTCCTTCCCATTAT-3', 'CCCTCCTTCGCGGATC-3', 'CCGGGATCAGGGTACA-3', 
                                       'CCGTGGAAGCTAGGCA-3', 'CCGTGGATCCATGCTC-3', 'CCTAAAGCACCATCCT-3', 'CCTACCACATCAGTCA-3', 'CGAACATTCAACGGCC-3', 'CGAGAAGGTGGTAACG-3', 
                                       'CGAGCACGTGCTCTTC-3', 'CGAGCACTCAGCTCGG-3', 'CGATTGAAGCGTAATA-3', 'CGCCAAGAGACCTAGG-3', 'CGCCAAGGTCTGCAAT-3', 'CGCGTTTTCAACCAAC-3', 
                                       'CGCGTTTTCGAACGGA-3', 'CGCTGGAAGTGATCGG-3', 'CGCTTCAAGACAAGCC-3', 'CGGACACGTCATCCCT-3', 'CGGGTCACAAGAAAGG-3', 'CGTCAGGCACCACCAG-3', 
                                       'CGTCCATTCTGCGTAA-3', 'CGTGAGCAGCGTTTAC-3', 'CGTGAGCAGTTGAGTA-3', 'CGTGTAACATGTCTCC-3', 'CGTGTAAGTCTGCGGT-3', 'CGTTCTGCAGACGTAG-3', 
                                       'CGTTCTGTCTGGAGCC-3', 'CTAATGGGTGATGATA-3', 'CTACGTCGTTCTCATT-3', 'CTAGTGAAGCAGATCG-3', 'CTAGTGAAGGAGTACC-3', 'CTCAGAACAAAGAATC-3', 
                                       'CTCGAAAAGCGATATA-3', 'CTCGTCAGTTGGTGGA-3', 'CTGAAACAGGATGCGT-3', 'CTGAAGTAGGGTTCCC-3', 'CTGCCTAAGACACTAA-3', 'CTGTTTATCCCTAACC-3', 
                                       'CTTAGGACAAAGGCGT-3', 'CTTCTCTTCCTTAATC-3', 'CTTGGCTCAGGATTGG-3', 'GAAACTCAGGTAGCTG-3', 'GAAACTCCAGATCGGA-3', 'GAAACTCGTACAGTGG-3', 
                                       'GAAACTCGTAGGGACT-3', 'GAAACTCTCCCTTGCA-3', 'GAACCTAGTACCGTAT-3', 'GAATAAGGTATAATGG-3', 'GAATAAGTCAACGGGA-3', 'GAATAAGTCGAGCCCA-3', 
                                       'GACAGAGGTTGGTTTG-3', 'GACGCGTTCGATCCCT-3', 'GACGGCTGTAGCGCTC-3', 'GAGCAGAAGAAACGAG-3', 'GATCTAGCAGTTTACG-3', 'GATGAGGCAGCTGTTA-3', 
                                       'GCAAACTTCCGTAGTA-3', 'GCATGTACAGGAATCG-3', 'GCGCAGTTCCAACCAA-3', 'GCGCCAACATATGCTG-3', 'GCGGGTTTCAACGAAA-3', 'GCTTCCAGTTCGTGAT-3', 
                                       'GCTTGAAGTCTAGAGG-3', 'GGAAAGCAGTTATCGC-3', 'GGAATAATCGTTGCCT-3', 'GGACAGAGTCGCTTTC-3', 'GGACAGATCAGGTTCA-3', 'GGACATTAGTAATCCC-3', 
                                       'GGACATTCACGGTAAG-3', 'GGATGTTTCCACGCAG-3', 'GGATTACCATAGGATA-3', 'GGATTACGTACTCGCG-3', 'GGCTGGTAGCCGGTAA-3', 'GGGAGATTCTTTAGGG-3', 
                                       'GGGATGAAGATCCTGT-3', 'GGGCACTAGCCGGTAA-3', 'GGGCACTAGTTACCCA-3', 'GGGCATCCAATCCAAC-3', 'GGGTTGCCATGGGAAC-3', 'GGTATTGAGCAACGGT-3', 
                                       'GGTATTGTCGGTGTTA-3', 'GTATTCTTCGTGGACC-3', 'GTCAAGTCAGTCGTGC-3', 'GTCAAGTCATTTGCCC-3', 'GTCACGGAGGCTAGGT-3', 'GTCATTTCAAGACACG-3', 
                                       'GTCATTTTCGTATCAG-3', 'GTCGGGTAGCCCAACC-3', 'GTGCAGCTCACTCTTA-3', 'GTGCAGCTCTTACCGC-3', 'GTGGGTCAGCACCGTC-3', 'GTGGGTCAGCGTGTCC-3', 
                                       'GTTACAGTCATTGCGA-3', 'GTTCATTGTGAGCGAT-3', 'GTTCTCGCATTAACCG-3', 'TAAACCGCAAGCGAGT-3', 'TAAGAGAGTTCTGTTT-3', 'TAAGCGTCATTCTTAC-3', 
                                       'TAAGTGCGTGCAACTT-3', 'TACACGATCTTTACGT-3', 'TACGGTATCCATTCTA-3', 'TACTCATGTCGCTTCT-3', 'TACTTACGTACTCAAC-3', 'TACTTACTCACCCTCA-3', 
                                       'TATGCCCTCACAGGCC-3', 'TCAATCTAGTGCTGCC-3', 'TCAGATGCAGCTTCGG-3', 'TCAGATGGTTCCACAA-3', 'TCAGCAAAGCCCAATT-3', 'TCAGGATTCTGGCGTG-3', 
                                       'TCAGGTACAATAGCAA-3', 'TCATTTGTCACCACCT-3', 'TCCCGATAGATGCCTT-3', 'TCCCGATTCGCCATAA-3', 'TCGGTAATCACCAGGC-3', 'TCTGAGACAAGCTGGA-3', 
                                       'TCTGGAAAGTCTCAAC-3', 'TGAGCCGGTTGAACTC-3', 'TGAGGGAGTGCAGACA-3', 'TGCCCATCATACAGCT-3', 'TGCCCTACACTTAACG-3', 'TGGACGCTCCCGACTT-3', 
                                       'TGGCTGGCACGAAATA-3', 'TGGTTAGCACCTCGTT-3', 'TGGTTCCTCCCTTGTG-3', 'TGGTTCCTCGCCAGCA-3', 'TGTATTCCACTAAGTC-3', 'TGTGTTTAGCCAGGAT-3', 
                                       'TGTGTTTGTGTGACCC-3', 'TTAACTCAGCTACCTA-3', 'TTCCCAGAGCTAACTC-3', 'TTCTACAGTCTAGTCA-3', 'TTGACTTAGACTTGAA-3', 'TTGACTTGTCATCGGC-3', 
                                       'TTGCCGTAGACAATAC-3', 'TTGCCGTTCTATCGCC-3', 'TTGCGTCTCGTGGTCG-3', 'TTGTAGGCAATGGACG-3', 'TTTACTGAGCAGGCTA-3', 'TTTATGCAGTGACTCT-3', 
                                       'TTTATGCGTTCGAATC-3', 'TTTGGTTAGAATTGTG-3', 'CGACCTTGTCTCAACA-3', 'CAAGAAATCATGCATG-3', 'ATTCTACGTAGGGACT-3', 'ATCACGATCGGACAAG-3', 
                                       'ACACTGAAGAATGTGT-3', 'AATCGGTCAATGACCT-3', 'TTCCCAGCAAACGCGA-3', 'TACCTTATCCTAGAAC-3', 'GATTCAGTCATCTGCC-3'), invert = F)

Clone9 <- SetIdent(Clone9, value = "Genotype")
markers.clone9 <- FindAllMarkers(Clone9,
                                 features = vargenes.filtered5,
                                 only.pos = FALSE, 
                                 min.pct = 0.01, 
                                 logfc.threshold = 0.25,
                                 test.use = "MAST",
                                 verbose = T,
                                 latent.vars = "nCount_RNA")


markers.clone9top20 <- markers.clone9 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone9, 
                     features = markers.clone9top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g9_recip <- WhichCells(CD8.clean3, idents = c('ATTGGACGTCTCCATC-2', 'CGAGAAGAGATGTCGG-2', 'CGCGGTATCTCTGTCG-2', 'AAGGTTCTCAGGTAAA-3', 'ACCAGTATCCTTCAAT-3', 'ACCGTAACAATCAGAA-3', 'AGAGCGAGTCCAACTA-3', 'AGGCCACCAGTGGAGT-3', 'CCTCTGATCTCAAACG-3', 'CGTTCTGGTCAAGCGA-3', 'CTAAGACTCGTAGATC-3', 'CTCTAATGTGTTCGAT-3', 'GGCCGATCAGATGGCA-3', 'GGCGACTGTCATGCAT-3', 'TCACGAACATTGAGCT-3', 'TCAGGATAGTGTCTCA-3', 'TCTCTAACAGCTTCGG-3', 'TGAGCATCAATTCCTT-3', 'TGCCCTAAGGGCACTA-3', 'TTTGGTTTCCCTTGCA-3'))
g9_don <- WhichCells(CD8.clean3, idents = c('AACCGCGTCACTTCAT-3', 'AACTCTTGTCGGGTCT-3', 'AAGGCAGAGGTAGCCA-3', 'AAGGTTCGTCGCCATG-3', 'AAGGTTCTCCACGTGG-3', 'AATCCAGAGTCTCGGC-3', 
                                            'ACACCAATCCACTGGG-3', 'ACACCGGGTCTAGTCA-3', 'ACACTGAAGAGTCTGG-3', 'ACAGCTATCTCTGCTG-3', 'ACATACGAGATCACGG-3', 'ACCAGTAAGTTACCCA-3', 
                                            'ACCGTAAGTGGTTTCA-3', 'ACGAGGACAAGTACCT-3', 'ACGATGTAGCTGGAAC-3', 'ACGGGCTAGCGTTTAC-3', 'ACGTCAATCCTAAGTG-3', 
                                            'ACGTCAATCTGCTTGC-3', 'ACTGAGTAGAGTAATC-3', 'ACTGCTCTCTACTTAC-3', 'AGAATAGAGGTTACCT-3', 'AGAGCGAGTCTTGATG-3', 'AGATCTGTCTTCAACT-3', 
                                            'AGCGTCGTCGCTAGCG-3', 'AGCTTGAAGGTGCACA-3', 'AGCTTGAGTGAGCGAT-3', 'AGGCCACGTGGGTCAA-3', 'AGGCCACTCGGAGCAA-3', 'AGGCCGTGTCTGCAAT-3', 
                                            'AGGGAGTTCCACGCAG-3', 'AGGGATGAGCGATTCT-3', 'AGGTCATGTGTGCGTC-3', 'AGGTCATGTTCGCGAC-3', 'AGTCTTTTCAGGTAAA-3', 'ATAACGCGTAAGGATT-3', 
                                            'ATAAGAGAGCAGCGTA-3', 'ATAAGAGGTGCACGAA-3', 'ATCATCTGTGTCAATC-3', 'ATCCGAAGTAGAAAGG-3', 'ATCTGCCTCACCGGGT-3', 'ATGTGTGGTATAGGGC-3', 
                                            'ATTGGACAGCTCCCAG-3', 'ATTGGACCAGGGAGAG-3', 'ATTGGTGCATCGGGTC-3', 'ATTTCTGAGACAAGCC-3', 'CAACCAAAGACAGGCT-3', 'CAAGATCGTAAGTAGT-3', 
                                            'CAAGATCGTATAGGGC-3', 'CAAGTTGTCCTGCCAT-3', 'CAAGTTGTCGTACGGC-3', 'CACACCTCACAGAGGT-3', 'CACACTCCATGAAGTA-3', 'CACATAGTCCTAGGGC-3', 
                                            'CACATTTAGGCGCTCT-3', 'CACCTTGAGATAGCAT-3', 'CAGATCACAGTTCCCT-3', 'CAGATCATCACCAGGC-3', 'CAGCATATCATCGATG-3', 'CATCAAGAGTAGTGCG-3', 
                                            'CATCCACCAGGCTGAA-3', 'CATCCACGTAGCGCAA-3', 'CATCGGGCAGCTTAAC-3', 'CCAATCCAGTTAACGA-3', 'CCACCTAGTGGTTTCA-3', 'CCAGCGAAGTGCAAGC-3', 
                                            'CCAGCGAGTAGCGCAA-3', 'CCCAATCGTTGAGTTC-3', 'CCCTCCTAGTATCTCG-3', 'CCCTCCTTCCCATTAT-3', 'CCCTCCTTCGCGGATC-3', 'CCGGGATCAGGGTACA-3', 
                                            'CCGTGGAAGCTAGGCA-3', 'CCGTGGATCCATGCTC-3', 'CCTAAAGCACCATCCT-3', 'CCTACCACATCAGTCA-3', 'CGAACATTCAACGGCC-3', 'CGAGAAGGTGGTAACG-3', 
                                            'CGAGCACGTGCTCTTC-3', 'CGAGCACTCAGCTCGG-3', 'CGATTGAAGCGTAATA-3', 'CGCCAAGAGACCTAGG-3', 'CGCCAAGGTCTGCAAT-3', 'CGCGTTTTCAACCAAC-3', 
                                            'CGCGTTTTCGAACGGA-3', 'CGCTGGAAGTGATCGG-3', 'CGCTTCAAGACAAGCC-3', 'CGGACACGTCATCCCT-3', 'CGGGTCACAAGAAAGG-3', 'CGTCAGGCACCACCAG-3', 
                                            'CGTCCATTCTGCGTAA-3', 'CGTGAGCAGCGTTTAC-3', 'CGTGAGCAGTTGAGTA-3', 'CGTGTAACATGTCTCC-3', 'CGTGTAAGTCTGCGGT-3', 'CGTTCTGCAGACGTAG-3', 
                                            'CGTTCTGTCTGGAGCC-3', 'CTAATGGGTGATGATA-3', 'CTACGTCGTTCTCATT-3', 'CTAGTGAAGCAGATCG-3', 'CTAGTGAAGGAGTACC-3', 'CTCAGAACAAAGAATC-3', 
                                            'CTCGAAAAGCGATATA-3', 'CTCGTCAGTTGGTGGA-3', 'CTGAAACAGGATGCGT-3', 'CTGAAGTAGGGTTCCC-3', 'CTGCCTAAGACACTAA-3', 'CTGTTTATCCCTAACC-3', 
                                            'CTTAGGACAAAGGCGT-3', 'CTTCTCTTCCTTAATC-3', 'CTTGGCTCAGGATTGG-3', 'GAAACTCAGGTAGCTG-3', 'GAAACTCCAGATCGGA-3', 'GAAACTCGTACAGTGG-3', 
                                            'GAAACTCGTAGGGACT-3', 'GAAACTCTCCCTTGCA-3', 'GAACCTAGTACCGTAT-3', 'GAATAAGGTATAATGG-3', 'GAATAAGTCAACGGGA-3', 'GAATAAGTCGAGCCCA-3', 
                                            'GACAGAGGTTGGTTTG-3', 'GACGCGTTCGATCCCT-3', 'GACGGCTGTAGCGCTC-3', 'GAGCAGAAGAAACGAG-3', 'GATCTAGCAGTTTACG-3', 'GATGAGGCAGCTGTTA-3', 
                                            'GCAAACTTCCGTAGTA-3', 'GCATGTACAGGAATCG-3', 'GCGCAGTTCCAACCAA-3', 'GCGCCAACATATGCTG-3', 'GCGGGTTTCAACGAAA-3', 'GCTTCCAGTTCGTGAT-3', 
                                            'GCTTGAAGTCTAGAGG-3', 'GGAAAGCAGTTATCGC-3', 'GGAATAATCGTTGCCT-3', 'GGACAGAGTCGCTTTC-3', 'GGACAGATCAGGTTCA-3', 'GGACATTAGTAATCCC-3', 
                                            'GGACATTCACGGTAAG-3', 'GGATGTTTCCACGCAG-3', 'GGATTACCATAGGATA-3', 'GGATTACGTACTCGCG-3', 'GGCTGGTAGCCGGTAA-3', 'GGGAGATTCTTTAGGG-3', 
                                            'GGGATGAAGATCCTGT-3', 'GGGCACTAGCCGGTAA-3', 'GGGCACTAGTTACCCA-3', 'GGGCATCCAATCCAAC-3', 'GGGTTGCCATGGGAAC-3', 'GGTATTGAGCAACGGT-3', 
                                            'GGTATTGTCGGTGTTA-3', 'GTATTCTTCGTGGACC-3', 'GTCAAGTCAGTCGTGC-3', 'GTCAAGTCATTTGCCC-3', 'GTCACGGAGGCTAGGT-3', 'GTCATTTCAAGACACG-3', 
                                            'GTCATTTTCGTATCAG-3', 'GTCGGGTAGCCCAACC-3', 'GTGCAGCTCACTCTTA-3', 'GTGCAGCTCTTACCGC-3', 'GTGGGTCAGCACCGTC-3', 'GTGGGTCAGCGTGTCC-3', 
                                            'GTTACAGTCATTGCGA-3', 'GTTCATTGTGAGCGAT-3', 'GTTCTCGCATTAACCG-3', 'TAAACCGCAAGCGAGT-3', 'TAAGAGAGTTCTGTTT-3', 'TAAGCGTCATTCTTAC-3', 
                                            'TAAGTGCGTGCAACTT-3', 'TACACGATCTTTACGT-3', 'TACGGTATCCATTCTA-3', 'TACTCATGTCGCTTCT-3', 'TACTTACGTACTCAAC-3', 'TACTTACTCACCCTCA-3', 
                                            'TATGCCCTCACAGGCC-3', 'TCAATCTAGTGCTGCC-3', 'TCAGATGCAGCTTCGG-3', 'TCAGATGGTTCCACAA-3', 'TCAGCAAAGCCCAATT-3', 'TCAGGATTCTGGCGTG-3', 
                                            'TCAGGTACAATAGCAA-3', 'TCATTTGTCACCACCT-3', 'TCCCGATAGATGCCTT-3', 'TCCCGATTCGCCATAA-3', 'TCGGTAATCACCAGGC-3', 'TCTGAGACAAGCTGGA-3', 
                                            'TCTGGAAAGTCTCAAC-3', 'TGAGCCGGTTGAACTC-3', 'TGAGGGAGTGCAGACA-3', 'TGCCCATCATACAGCT-3', 'TGCCCTACACTTAACG-3', 'TGGACGCTCCCGACTT-3', 
                                            'TGGCTGGCACGAAATA-3', 'TGGTTAGCACCTCGTT-3', 'TGGTTCCTCCCTTGTG-3', 'TGGTTCCTCGCCAGCA-3', 'TGTATTCCACTAAGTC-3', 'TGTGTTTAGCCAGGAT-3', 
                                            'TGTGTTTGTGTGACCC-3', 'TTAACTCAGCTACCTA-3', 'TTCCCAGAGCTAACTC-3', 'TTCTACAGTCTAGTCA-3', 'TTGACTTAGACTTGAA-3', 'TTGACTTGTCATCGGC-3', 
                                            'TTGCCGTAGACAATAC-3', 'TTGCCGTTCTATCGCC-3', 'TTGCGTCTCGTGGTCG-3', 'TTGTAGGCAATGGACG-3', 'TTTACTGAGCAGGCTA-3', 'TTTATGCAGTGACTCT-3', 
                                            'TTTATGCGTTCGAATC-3', 'TTTGGTTAGAATTGTG-3', 'CGACCTTGTCTCAACA-3', 'CAAGAAATCATGCATG-3', 'ATTCTACGTAGGGACT-3', 'ATCACGATCGGACAAG-3', 
                                            'ACACTGAAGAATGTGT-3', 'AATCGGTCAATGACCT-3', 'TTCCCAGCAAACGCGA-3', 'TACCTTATCCTAGAAC-3', 'GATTCAGTCATCTGCC-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g9_recip, g9_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")


#10
Recip10 =  c('CCGGGATTCGACGGAA-2', 'GATCGTAGTAGCAAAT-2', 'TTCTACACACATGGGA-2', 'AAGTCTGCAACTGGCC-3', 'ATCCACCGTAACGCGA-3', 'GCACATACAGCTGCTG-3', 'TAAGCGTAGAGTACCG-3', 'TGAGGGACAGTCGATT-3')
Don10 =  c('ACGATGTGTAGCGTAG-2', 'GAACATCCACCACGTG-2', 'TGAGCATGTACCAGTT-2', 'ATGTGTGTCACCGTAA-3', 'ATTGGTGTCGTAGGAG-3', 'GCTCCTATCGTTGCCT-3', 'AACCGCGCAGTATGCT-3', 'AGGTCCGCAAGCTGAG-3', 'ATGTGTGAGAGGTTGC-3', 'CACACTCCAGATGAGC-3', 'CCTACACGTGGTTTCA-3', 'CGAATGTGTGAGTGAC-3', 'CGGCTAGCACGGTAGA-3', 'GAAATGACAATCTACG-3', 'TAAGCGTCAGCGTTCG-3', 'TGCGTGGTCTGTCCGT-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone10 <- subset(CD8.clean3, ident = c('CCGGGATTCGACGGAA-2', 'GATCGTAGTAGCAAAT-2', 'AAGTCTGCAACTGGCC-3', 'ATCCACCGTAACGCGA-3', 'GCACATACAGCTGCTG-3', 'TAAGCGTAGAGTACCG-3', 'TGAGGGACAGTCGATT-3', 'ACGATGTGTAGCGTAG-2', 'GAACATCCACCACGTG-2', 'TGAGCATGTACCAGTT-2', 'ATGTGTGTCACCGTAA-3', 'ATTGGTGTCGTAGGAG-3', 'GCTCCTATCGTTGCCT-3', 'AACCGCGCAGTATGCT-3', 'AGGTCCGCAAGCTGAG-3', 'ATGTGTGAGAGGTTGC-3', 'CACACTCCAGATGAGC-3', 'CCTACACGTGGTTTCA-3', 'CGAATGTGTGAGTGAC-3', 'CGGCTAGCACGGTAGA-3', 'GAAATGACAATCTACG-3', 'TAAGCGTCAGCGTTCG-3', 'TGCGTGGTCTGTCCGT-3'), invert = F)

Clone10 <- SetIdent(Clone10, value = "Genotype")
markers.clone10 <- FindAllMarkers(Clone10,
                                  features = vargenes.filtered5,
                                  only.pos = FALSE, 
                                  min.pct = 0.01, 
                                  logfc.threshold = 0.25,
                                  test.use = "MAST",
                                  verbose = T,
                                  latent.vars = "nCount_RNA")


markers.clone10top20 <- markers.clone10 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone10, 
                     features = markers.clone10top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g10_recip <- WhichCells(CD8.clean3, idents = c('CCGGGATTCGACGGAA-2', 'GATCGTAGTAGCAAAT-2', 'AAGTCTGCAACTGGCC-3', 'ATCCACCGTAACGCGA-3', 'GCACATACAGCTGCTG-3', 'TAAGCGTAGAGTACCG-3', 'TGAGGGACAGTCGATT-3'))
g10_don <- WhichCells(CD8.clean3, idents = c('ACGATGTGTAGCGTAG-2', 'GAACATCCACCACGTG-2', 'TGAGCATGTACCAGTT-2', 'ATGTGTGTCACCGTAA-3', 'ATTGGTGTCGTAGGAG-3', 'GCTCCTATCGTTGCCT-3', 'AACCGCGCAGTATGCT-3', 'AGGTCCGCAAGCTGAG-3', 'ATGTGTGAGAGGTTGC-3', 'CACACTCCAGATGAGC-3', 'CCTACACGTGGTTTCA-3', 'CGAATGTGTGAGTGAC-3', 'CGGCTAGCACGGTAGA-3', 'GAAATGACAATCTACG-3', 'TAAGCGTCAGCGTTCG-3', 'TGCGTGGTCTGTCCGT-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g10_recip, g10_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")



#11
Recip11 =  c('ACAGCCGGTCGAATCT-2', 'AGGTCATGTTGGTTTG-2', 'CAGCTAAGTTGGGACA-2', 'CCTAGCTAGAGACGAA-2', 'CTGCCTACAAACGCGA-2', 'CTGTGCTAGTCACGCC-2', 'GAACCTAAGTCCTCCT-2', 'GACAGAGGTCCATGAT-2', 'GCTCCTATCCTTTACA-2', 'TGAGGGAGTCCGAATT-2', 'AGCGTATAGTTTCCTT-2', 'ACGAGGATCGTGACAT-3', 'CCATTCGTCTCATTCA-3', 'CGTGTAATCTACCAGA-3', 'GCCTCTAGTCGTGGCT-3', 'TCAGCAAAGACCTAGG-3', 'TCTCTAACACGGCCAT-3', 'TTCGGTCCAGCTCCGA-3', 'TGAGGGAAGAGCCTAG-3')
Don11 =  c('AAATGCCTCAGAGCTT-2', 'CGATGTACATCCGCGA-2', 'CGCTATCCAGCTTAAC-2', 'GAACGGATCACGATGT-2', 'GTACTTTTCAAACCAC-2', 'TACACGAAGGGTATCG-2',
           'TACTTACAGTCAATAG-2', 'TACTTACGTAGCACGA-2', 'TGACTAGAGGTAGCTG-2', 'TGGTTCCTCTTAACCT-2', 'CGTCACTAGTCGAGTG-2', 'AAAGCAAGTCGCGGTT-3',
           'AACCGCGGTGATAAAC-3', 'ACACCCTCATTGCGGC-3', 'ACACCGGGTTCACGGC-3', 'ACGCCAGCAGACGCTC-3', 'ACTGAACTCAAGGCTT-3', 'AGAGCGAGTGGAAAGA-3',
           'AGATCTGTCCGCAGTG-3', 'AGCGGTCGTTATCCGA-3', 'AGCTTGAAGGCTCAGA-3', 'AGGGTGAGTTGGTTTG-3', 'AGGTCCGTCGTATCAG-3', 'ATCATCTGTCCGTCAG-3', 
           'ATCATGGCATGCCACG-3', 'CACATTTTCACTATTC-3', 'CAGAGAGCAGCCACCA-3', 'CAGCTAATCACATAGC-3', 'CAGTCCTGTAGGGTAC-3', 'CATCGGGAGACCTTTG-3', 
           'CCGGTAGAGCGGCTTC-3', 'CCGTACTTCTATCGCC-3', 'CCGTGGAGTACCTACA-3', 'CCTTACGCAAGCGTAG-3', 'CCTTTCTGTACTTGAC-3', 'CGAGCACGTAGCGTGA-3',
           'CGAGCCAAGTATGACA-3', 'CGCGTTTGTGATAAGT-3', 'CGGAGTCCATGGTAGG-3', 'CGTAGCGAGCGGCTTC-3', 'CGTGAGCGTTCGTTGA-3', 'CGTGAGCTCTTGTTTG-3',
           'CTAAGACGTTTGTTTC-3', 'CTACACCGTTAAGAAC-3', 'CTCACACTCTCCAGGG-3', 'CTCATTAGTGTTCGAT-3', 'CTCCTAGGTCGACTAT-3', 'CTCGTACTCTCGAGTA-3',
           'CTCTAATAGGGCATGT-3', 'CTTAACTAGTGGGCTA-3', 'CTTACCGTCGGGAGTA-3', 'CTTAGGATCCGCGGTA-3', 'GAAACTCGTCTAGCCG-3', 'GAACGGACACAGCCCA-3',
           'GACAGAGCAGCCACCA-3', 'GATGAGGGTAAGTAGT-3', 'GCAATCATCACGGTTA-3', 'GCAGTTAAGCGTAGTG-3', 'GCATGATGTGCGAAAC-3', 'GCGAGAAGTGGCAAAC-3', 
           'GCGCCAAAGAGGTAGA-3', 'GCGCCAACAGGTGGAT-3', 'GGAAAGCCATCACGTA-3', 'GGTATTGAGGCTAGGT-3', 'GTACGTATCGGCGCTA-3', 'GTAGGCCCACATTTCT-3', 
           'GTCGTAAGTGATGTGG-3', 'GTGTGCGAGACCTTTG-3', 'GTGTGCGTCCACGACG-3', 'GTTAAGCCAGCGATCC-3', 'GTTCATTAGTACTTGC-3', 'TAAGTGCGTTACTGAC-3',
           'TACGGTAGTTGACGTT-3', 'TACTCGCCACAAGACG-3', 'TAGTGGTAGAAGGGTA-3', 'TATTACCTCAACGAAA-3', 'TCAACGAGTTCTGTTT-3', 'TCAATCTCAAAGGTGC-3',
           'TCCCGATAGGAGTTGC-3', 'TCGCGAGGTTAAGAAC-3', 'TCTCATAAGTCTCAAC-3', 'TCTCATAGTTATGCGT-3', 'TCTGAGACATACCATG-3', 'TGCCCATCACGTCTCT-3', 
           'TGCGGGTGTGACGGTA-3', 'TGGACGCCATCCAACA-3', 'TGTGTTTAGTACGACG-3', 'TTAGGACTCGGTCCGA-3', 'TTATGCTGTTCGTCTC-3', 'TTGGCAAAGTCTCCTC-3',
           'GTTCGGGGTAGGAGTC-3', 'TAGACCAGTAGAGCTG-3', 'CAGCGACTCATCGATG-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone11 <- subset(CD8.clean3, ident = c('ACAGCCGGTCGAATCT-2', 'AGGTCATGTTGGTTTG-2', 'CAGCTAAGTTGGGACA-2', 'CCTAGCTAGAGACGAA-2', 'CTGCCTACAAACGCGA-2', 
                                        'CTGTGCTAGTCACGCC-2', 'GAACCTAAGTCCTCCT-2', 'GACAGAGGTCCATGAT-2', 'GCTCCTATCCTTTACA-2', 'TGAGGGAGTCCGAATT-2', 
                                        'AGCGTATAGTTTCCTT-2', 'ACGAGGATCGTGACAT-3', 'CCATTCGTCTCATTCA-3', 'CGTGTAATCTACCAGA-3', 'GCCTCTAGTCGTGGCT-3', 
                                        'TCAGCAAAGACCTAGG-3', 'TCTCTAACACGGCCAT-3', 'TTCGGTCCAGCTCCGA-3', 'TGAGGGAAGAGCCTAG-3', 'AAATGCCTCAGAGCTT-2', 'CGATGTACATCCGCGA-2', 'CGCTATCCAGCTTAAC-2', 'GAACGGATCACGATGT-2', 'GTACTTTTCAAACCAC-2', 'TACACGAAGGGTATCG-2',
                                        'TACTTACAGTCAATAG-2', 'TACTTACGTAGCACGA-2', 'TGGTTCCTCTTAACCT-2', 'CGTCACTAGTCGAGTG-2', 'AAAGCAAGTCGCGGTT-3',
                                        'AACCGCGGTGATAAAC-3', 'ACACCCTCATTGCGGC-3', 'ACACCGGGTTCACGGC-3', 'ACGCCAGCAGACGCTC-3', 'ACTGAACTCAAGGCTT-3', 'AGAGCGAGTGGAAAGA-3',
                                        'AGATCTGTCCGCAGTG-3', 'AGCGGTCGTTATCCGA-3', 'AGCTTGAAGGCTCAGA-3', 'AGGGTGAGTTGGTTTG-3', 'AGGTCCGTCGTATCAG-3', 'ATCATCTGTCCGTCAG-3', 
                                        'ATCATGGCATGCCACG-3', 'CACATTTTCACTATTC-3', 'CAGAGAGCAGCCACCA-3', 'CAGCTAATCACATAGC-3', 'CAGTCCTGTAGGGTAC-3', 'CATCGGGAGACCTTTG-3', 
                                        'CCGGTAGAGCGGCTTC-3', 'CCGTACTTCTATCGCC-3', 'CCGTGGAGTACCTACA-3', 'CCTTACGCAAGCGTAG-3', 'CCTTTCTGTACTTGAC-3', 'CGAGCACGTAGCGTGA-3',
                                        'CGAGCCAAGTATGACA-3', 'CGCGTTTGTGATAAGT-3', 'CGGAGTCCATGGTAGG-3', 'CGTAGCGAGCGGCTTC-3', 'CGTGAGCGTTCGTTGA-3', 'CGTGAGCTCTTGTTTG-3',
                                        'CTAAGACGTTTGTTTC-3', 'CTACACCGTTAAGAAC-3', 'CTCACACTCTCCAGGG-3', 'CTCATTAGTGTTCGAT-3', 'CTCCTAGGTCGACTAT-3', 'CTCGTACTCTCGAGTA-3',
                                        'CTCTAATAGGGCATGT-3', 'CTTAACTAGTGGGCTA-3', 'CTTACCGTCGGGAGTA-3', 'CTTAGGATCCGCGGTA-3', 'GAAACTCGTCTAGCCG-3', 'GAACGGACACAGCCCA-3',
                                        'GACAGAGCAGCCACCA-3', 'GATGAGGGTAAGTAGT-3', 'GCAATCATCACGGTTA-3', 'GCAGTTAAGCGTAGTG-3', 'GCATGATGTGCGAAAC-3', 'GCGAGAAGTGGCAAAC-3', 
                                        'GCGCCAAAGAGGTAGA-3', 'GCGCCAACAGGTGGAT-3', 'GGAAAGCCATCACGTA-3', 'GGTATTGAGGCTAGGT-3', 'GTACGTATCGGCGCTA-3', 'GTAGGCCCACATTTCT-3', 
                                        'GTCGTAAGTGATGTGG-3', 'GTGTGCGAGACCTTTG-3', 'GTGTGCGTCCACGACG-3', 'GTTAAGCCAGCGATCC-3', 'GTTCATTAGTACTTGC-3', 'TAAGTGCGTTACTGAC-3',
                                        'TACGGTAGTTGACGTT-3', 'TACTCGCCACAAGACG-3', 'TAGTGGTAGAAGGGTA-3', 'TATTACCTCAACGAAA-3', 'TCAACGAGTTCTGTTT-3', 'TCAATCTCAAAGGTGC-3',
                                        'TCCCGATAGGAGTTGC-3', 'TCGCGAGGTTAAGAAC-3', 'TCTCATAAGTCTCAAC-3', 'TCTCATAGTTATGCGT-3', 'TCTGAGACATACCATG-3', 'TGCCCATCACGTCTCT-3', 
                                        'TGCGGGTGTGACGGTA-3', 'TGGACGCCATCCAACA-3', 'TGTGTTTAGTACGACG-3', 'TTAGGACTCGGTCCGA-3', 'TTATGCTGTTCGTCTC-3', 'TTGGCAAAGTCTCCTC-3',
                                        'GTTCGGGGTAGGAGTC-3', 'TAGACCAGTAGAGCTG-3', 'CAGCGACTCATCGATG-3'), invert = F)

Clone11 <- SetIdent(Clone11, value = "Genotype")
markers.clone11 <- FindAllMarkers(Clone11,
                                  features = vargenes.filtered5,
                                  only.pos = FALSE, 
                                  min.pct = 0.01, 
                                  logfc.threshold = 0.25,
                                  test.use = "MAST",
                                  verbose = T,
                                  latent.vars = "nCount_RNA")


markers.clone11top20 <- markers.clone11 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone11, 
                     features = markers.clone11top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g11_recip <- WhichCells(CD8.clean3, idents = c('ACAGCCGGTCGAATCT-2', 'AGGTCATGTTGGTTTG-2', 'CAGCTAAGTTGGGACA-2', 'CCTAGCTAGAGACGAA-2', 'CTGCCTACAAACGCGA-2', 'CTGTGCTAGTCACGCC-2', 'GAACCTAAGTCCTCCT-2', 'GACAGAGGTCCATGAT-2', 'GCTCCTATCCTTTACA-2', 'TGAGGGAGTCCGAATT-2', 'AGCGTATAGTTTCCTT-2', 'ACGAGGATCGTGACAT-3', 'CCATTCGTCTCATTCA-3', 'CGTGTAATCTACCAGA-3', 'GCCTCTAGTCGTGGCT-3', 'TCAGCAAAGACCTAGG-3', 'TCTCTAACACGGCCAT-3', 'TTCGGTCCAGCTCCGA-3', 'TGAGGGAAGAGCCTAG-3'))
g11_don <- WhichCells(CD8.clean3, idents = c('AAATGCCTCAGAGCTT-2', 'CGATGTACATCCGCGA-2', 'CGCTATCCAGCTTAAC-2', 'GAACGGATCACGATGT-2', 'GTACTTTTCAAACCAC-2', 'TACACGAAGGGTATCG-2',
                                             'TACTTACAGTCAATAG-2', 'TACTTACGTAGCACGA-2',  'TGGTTCCTCTTAACCT-2', 'CGTCACTAGTCGAGTG-2', 'AAAGCAAGTCGCGGTT-3',
                                             'AACCGCGGTGATAAAC-3', 'ACACCCTCATTGCGGC-3', 'ACACCGGGTTCACGGC-3', 'ACGCCAGCAGACGCTC-3', 'ACTGAACTCAAGGCTT-3', 'AGAGCGAGTGGAAAGA-3',
                                             'AGATCTGTCCGCAGTG-3', 'AGCGGTCGTTATCCGA-3', 'AGCTTGAAGGCTCAGA-3', 'AGGGTGAGTTGGTTTG-3', 'AGGTCCGTCGTATCAG-3', 'ATCATCTGTCCGTCAG-3', 
                                             'ATCATGGCATGCCACG-3', 'CACATTTTCACTATTC-3', 'CAGAGAGCAGCCACCA-3', 'CAGCTAATCACATAGC-3', 'CAGTCCTGTAGGGTAC-3', 'CATCGGGAGACCTTTG-3', 
                                             'CCGGTAGAGCGGCTTC-3', 'CCGTACTTCTATCGCC-3', 'CCGTGGAGTACCTACA-3', 'CCTTACGCAAGCGTAG-3', 'CCTTTCTGTACTTGAC-3', 'CGAGCACGTAGCGTGA-3',
                                             'CGAGCCAAGTATGACA-3', 'CGCGTTTGTGATAAGT-3', 'CGGAGTCCATGGTAGG-3', 'CGTAGCGAGCGGCTTC-3', 'CGTGAGCGTTCGTTGA-3', 'CGTGAGCTCTTGTTTG-3',
                                             'CTAAGACGTTTGTTTC-3', 'CTACACCGTTAAGAAC-3', 'CTCACACTCTCCAGGG-3', 'CTCATTAGTGTTCGAT-3', 'CTCCTAGGTCGACTAT-3', 'CTCGTACTCTCGAGTA-3',
                                             'CTCTAATAGGGCATGT-3', 'CTTAACTAGTGGGCTA-3', 'CTTACCGTCGGGAGTA-3', 'CTTAGGATCCGCGGTA-3', 'GAAACTCGTCTAGCCG-3', 'GAACGGACACAGCCCA-3',
                                             'GACAGAGCAGCCACCA-3', 'GATGAGGGTAAGTAGT-3', 'GCAATCATCACGGTTA-3', 'GCAGTTAAGCGTAGTG-3', 'GCATGATGTGCGAAAC-3', 'GCGAGAAGTGGCAAAC-3', 
                                             'GCGCCAAAGAGGTAGA-3', 'GCGCCAACAGGTGGAT-3', 'GGAAAGCCATCACGTA-3', 'GGTATTGAGGCTAGGT-3', 'GTACGTATCGGCGCTA-3', 'GTAGGCCCACATTTCT-3', 
                                             'GTCGTAAGTGATGTGG-3', 'GTGTGCGAGACCTTTG-3', 'GTGTGCGTCCACGACG-3', 'GTTAAGCCAGCGATCC-3', 'GTTCATTAGTACTTGC-3', 'TAAGTGCGTTACTGAC-3',
                                             'TACGGTAGTTGACGTT-3', 'TACTCGCCACAAGACG-3', 'TAGTGGTAGAAGGGTA-3', 'TATTACCTCAACGAAA-3', 'TCAACGAGTTCTGTTT-3', 'TCAATCTCAAAGGTGC-3',
                                             'TCCCGATAGGAGTTGC-3', 'TCGCGAGGTTAAGAAC-3', 'TCTCATAAGTCTCAAC-3', 'TCTCATAGTTATGCGT-3', 'TCTGAGACATACCATG-3', 'TGCCCATCACGTCTCT-3', 
                                             'TGCGGGTGTGACGGTA-3', 'TGGACGCCATCCAACA-3', 'TGTGTTTAGTACGACG-3', 'TTAGGACTCGGTCCGA-3', 'TTATGCTGTTCGTCTC-3', 'TTGGCAAAGTCTCCTC-3',
                                             'GTTCGGGGTAGGAGTC-3', 'TAGACCAGTAGAGCTG-3', 'CAGCGACTCATCGATG-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g11_recip, g11_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")

#12
Recip12 =  c('ACAGCCGGTCGAATCT-2', 'AGGTCATGTTGGTTTG-2', 'CAGCTAAGTTGGGACA-2', 'CCTAGCTAGAGACGAA-2', 'CTGCCTACAAACGCGA-2', 'CTGTGCTAGTCACGCC-2',
             'GAACCTAAGTCCTCCT-2', 'GACAGAGGTCCATGAT-2', 'GCTCCTATCCTTTACA-2', 'TGAGGGAGTCCGAATT-2', 'AGCGTATAGTTTCCTT-2', 'ACGAGGATCGTGACAT-3',
             'CCATTCGTCTCATTCA-3', 'CGTGTAATCTACCAGA-3', 'GCCTCTAGTCGTGGCT-3', 'TCAGCAAAGACCTAGG-3', 'TCTCTAACACGGCCAT-3', 'TTCGGTCCAGCTCCGA-3',
             'TGAGGGAAGAGCCTAG-3')
Don12 =  c('ACGATGTGTAGCGTAG-2', 'GAACATCCACCACGTG-2', 'TGAGCATGTACCAGTT-2', 'ATGTGTGTCACCGTAA-3', 'ATTGGTGTCGTAGGAG-3', 'GCTCCTATCGTTGCCT-3',
           'AACCGCGCAGTATGCT-3', 'AGGTCCGCAAGCTGAG-3', 'ATGTGTGAGAGGTTGC-3', 'CACACTCCAGATGAGC-3', 'CCTACACGTGGTTTCA-3', 'CGAATGTGTGAGTGAC-3',
           'CGGCTAGCACGGTAGA-3', 'GAAATGACAATCTACG-3', 'TAAGCGTCAGCGTTCG-3', 'TGCGTGGTCTGTCCGT-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone12 <- subset(CD8.clean3, ident = c('ACAGCCGGTCGAATCT-2', 'AGGTCATGTTGGTTTG-2', 'CAGCTAAGTTGGGACA-2', 'CCTAGCTAGAGACGAA-2', 'CTGCCTACAAACGCGA-2', 'CTGTGCTAGTCACGCC-2',
                                        'GAACCTAAGTCCTCCT-2', 'GACAGAGGTCCATGAT-2', 'GCTCCTATCCTTTACA-2', 'TGAGGGAGTCCGAATT-2', 'AGCGTATAGTTTCCTT-2', 'ACGAGGATCGTGACAT-3',
                                        'CCATTCGTCTCATTCA-3', 'CGTGTAATCTACCAGA-3', 'GCCTCTAGTCGTGGCT-3', 'TCAGCAAAGACCTAGG-3', 'TCTCTAACACGGCCAT-3', 'TTCGGTCCAGCTCCGA-3',
                                        'TGAGGGAAGAGCCTAG-3', 'ACGATGTGTAGCGTAG-2', 'GAACATCCACCACGTG-2', 'TGAGCATGTACCAGTT-2', 'ATGTGTGTCACCGTAA-3', 'ATTGGTGTCGTAGGAG-3', 'GCTCCTATCGTTGCCT-3',
                                        'AACCGCGCAGTATGCT-3', 'AGGTCCGCAAGCTGAG-3', 'ATGTGTGAGAGGTTGC-3', 'CACACTCCAGATGAGC-3', 'CCTACACGTGGTTTCA-3', 'CGAATGTGTGAGTGAC-3',
                                        'CGGCTAGCACGGTAGA-3', 'GAAATGACAATCTACG-3', 'TAAGCGTCAGCGTTCG-3', 'TGCGTGGTCTGTCCGT-3'), invert = F)

Clone12 <- SetIdent(Clone12, value = "Genotype")
markers.clone12 <- FindAllMarkers(Clone12,
                                  features = vargenes.filtered5,
                                  only.pos = FALSE, 
                                  min.pct = 0.01, 
                                  logfc.threshold = 0.25,
                                  test.use = "MAST",
                                  verbose = T,
                                  latent.vars = "nCount_RNA")


markers.clone12top20 <- markers.clone12 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone12, 
                     features = markers.clone12top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g12_recip <- WhichCells(CD8.clean3, idents = c('ACAGCCGGTCGAATCT-2', 'AGGTCATGTTGGTTTG-2', 'CAGCTAAGTTGGGACA-2', 'CCTAGCTAGAGACGAA-2', 'CTGCCTACAAACGCGA-2', 'CTGTGCTAGTCACGCC-2',
                                               'GAACCTAAGTCCTCCT-2', 'GACAGAGGTCCATGAT-2', 'GCTCCTATCCTTTACA-2', 'TGAGGGAGTCCGAATT-2', 'AGCGTATAGTTTCCTT-2', 'ACGAGGATCGTGACAT-3',
                                               'CCATTCGTCTCATTCA-3', 'CGTGTAATCTACCAGA-3', 'GCCTCTAGTCGTGGCT-3', 'TCAGCAAAGACCTAGG-3', 'TCTCTAACACGGCCAT-3', 'TTCGGTCCAGCTCCGA-3',
                                               'TGAGGGAAGAGCCTAG-3'))
g12_don <- WhichCells(CD8.clean3, idents = c('ACGATGTGTAGCGTAG-2', 'GAACATCCACCACGTG-2', 'TGAGCATGTACCAGTT-2', 'ATGTGTGTCACCGTAA-3', 'ATTGGTGTCGTAGGAG-3', 'GCTCCTATCGTTGCCT-3',
                                             'AACCGCGCAGTATGCT-3', 'AGGTCCGCAAGCTGAG-3', 'ATGTGTGAGAGGTTGC-3', 'CACACTCCAGATGAGC-3', 'CCTACACGTGGTTTCA-3', 'CGAATGTGTGAGTGAC-3',
                                             'CGGCTAGCACGGTAGA-3', 'GAAATGACAATCTACG-3', 'TAAGCGTCAGCGTTCG-3', 'TGCGTGGTCTGTCCGT-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g12_recip, g12_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")

#13
Recip13 =  c('ACGAGGAGTGTGAATA-1', 'AGGGATGGTTTGTTTC-1', 'GATCGATTCCACGTGG-1', 'GGGATGAGTACATCCA-1', 'TCGAGGCCAGGGCATA-1', 'AAACCTGAGGGATACC-2',
             'AACTCCCAGTCGATAA-2', 'ATAGACCCAGGACCCT-2', 'CAGAGAGTCCGATATG-2', 'CGCTGGATCGTCGTTC-2', 'CTGTGCTCAGGAATGC-2', 'TGGCCAGCATCGGGTC-2', 
             'ACGCAGCAGGCTCTTA-3', 'CCCTCCTGTCTCCCTA-3', 'GTCTTCGCAGTCAGAG-3', 'TGCTACCGTTGCGTTA-3')
Don13 =  c('AACCGCGTCCTGCTTG-3', 'ACGGGCTGTCTCACCT-3', 'ATGGGAGAGGTGCAAC-3', 'ATTGGACTCTTACCGC-3', 'CAGATCAAGAAGGTTT-3', 'CAGCAGCTCTCTTATG-3',
           'CAGCGACCAGGATTGG-3', 'CATCCACCAATAGAGT-3', 'CGTTGGGTCGTAGATC-3', 'CTGAAACGTCGTGGCT-3', 'CTGAAACTCTGTCTAT-3', 'CTTGGCTCAAGCCCAC-3',
           'GAAATGACAAGCCTAT-3', 'GCACTCTTCGCAAGCC-3', 'GCAGTTAGTTGAGTTC-3', 'GCATGCGTCCTTTACA-3', 'GCGAGAACAGATCCAT-3', 'GTCCTCAGTCGGGTCT-3',
           'GTCGGGTCAGGGTATG-3', 'GTGCATACAGCATACT-3', 'GTTACAGGTATCGCAT-3', 'TAAGCGTCAAGGTTCT-3', 'TGGACGCAGAGCAATT-3', 'TTTGGTTTCAGATAAG-3',
           'CAACTAGCACATCCAA-3', 'ACTGATGTCTTAGCCC-3', 'TAGTTGGGTTAAAGAC-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone13 <- subset(CD8.clean3, ident = c('ACGAGGAGTGTGAATA-1', 'AGGGATGGTTTGTTTC-1', 'GATCGATTCCACGTGG-1', 'GGGATGAGTACATCCA-1', 'TCGAGGCCAGGGCATA-1', 'AAACCTGAGGGATACC-2',
                                        'AACTCCCAGTCGATAA-2', 'ATAGACCCAGGACCCT-2', 'CAGAGAGTCCGATATG-2', 'CGCTGGATCGTCGTTC-2', 'CTGTGCTCAGGAATGC-2', 'TGGCCAGCATCGGGTC-2', 
                                        'ACGCAGCAGGCTCTTA-3', 'CCCTCCTGTCTCCCTA-3', 'GTCTTCGCAGTCAGAG-3', 'TGCTACCGTTGCGTTA-3', 'AACCGCGTCCTGCTTG-3', 'ACGGGCTGTCTCACCT-3', 'ATGGGAGAGGTGCAAC-3', 'ATTGGACTCTTACCGC-3', 'CAGATCAAGAAGGTTT-3', 'CAGCAGCTCTCTTATG-3',
                                        'CAGCGACCAGGATTGG-3', 'CATCCACCAATAGAGT-3', 'CGTTGGGTCGTAGATC-3', 'CTGAAACGTCGTGGCT-3', 'CTGAAACTCTGTCTAT-3', 'CTTGGCTCAAGCCCAC-3',
                                        'GAAATGACAAGCCTAT-3', 'GCACTCTTCGCAAGCC-3', 'GCAGTTAGTTGAGTTC-3', 'GCATGCGTCCTTTACA-3', 'GCGAGAACAGATCCAT-3', 'GTCCTCAGTCGGGTCT-3',
                                        'GTCGGGTCAGGGTATG-3', 'GTGCATACAGCATACT-3', 'GTTACAGGTATCGCAT-3', 'TAAGCGTCAAGGTTCT-3', 'TGGACGCAGAGCAATT-3', 'TTTGGTTTCAGATAAG-3',
                                        'CAACTAGCACATCCAA-3', 'ACTGATGTCTTAGCCC-3', 'TAGTTGGGTTAAAGAC-3'), invert = F)

Clone13 <- SetIdent(Clone13, value = "Genotype")
markers.clone13 <- FindAllMarkers(Clone13,
                                  features = vargenes.filtered5,
                                  only.pos = FALSE, 
                                  min.pct = 0.01, 
                                  logfc.threshold = 0.25,
                                  test.use = "MAST",
                                  verbose = T,
                                  latent.vars = "nCount_RNA")


markers.clone13top20 <- markers.clone13 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone13, 
                     features = markers.clone13top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g13_recip <- WhichCells(CD8.clean3, idents = c('ACGAGGAGTGTGAATA-1', 'AGGGATGGTTTGTTTC-1', 'GATCGATTCCACGTGG-1', 'GGGATGAGTACATCCA-1', 'TCGAGGCCAGGGCATA-1', 'AAACCTGAGGGATACC-2',
                                               'AACTCCCAGTCGATAA-2', 'ATAGACCCAGGACCCT-2', 'CAGAGAGTCCGATATG-2', 'CGCTGGATCGTCGTTC-2', 'CTGTGCTCAGGAATGC-2', 'TGGCCAGCATCGGGTC-2', 
                                               'ACGCAGCAGGCTCTTA-3', 'CCCTCCTGTCTCCCTA-3', 'GTCTTCGCAGTCAGAG-3', 'TGCTACCGTTGCGTTA-3'))
g13_don <- WhichCells(CD8.clean3, idents = c('AACCGCGTCCTGCTTG-3', 'ACGGGCTGTCTCACCT-3', 'ATGGGAGAGGTGCAAC-3', 'ATTGGACTCTTACCGC-3', 'CAGATCAAGAAGGTTT-3', 'CAGCAGCTCTCTTATG-3',
                                             'CAGCGACCAGGATTGG-3', 'CATCCACCAATAGAGT-3', 'CGTTGGGTCGTAGATC-3', 'CTGAAACGTCGTGGCT-3', 'CTGAAACTCTGTCTAT-3', 'CTTGGCTCAAGCCCAC-3',
                                             'GAAATGACAAGCCTAT-3', 'GCACTCTTCGCAAGCC-3', 'GCAGTTAGTTGAGTTC-3', 'GCATGCGTCCTTTACA-3', 'GCGAGAACAGATCCAT-3', 'GTCCTCAGTCGGGTCT-3',
                                             'GTCGGGTCAGGGTATG-3', 'GTGCATACAGCATACT-3', 'GTTACAGGTATCGCAT-3', 'TAAGCGTCAAGGTTCT-3', 'TGGACGCAGAGCAATT-3', 'TTTGGTTTCAGATAAG-3',
                                             'CAACTAGCACATCCAA-3', 'ACTGATGTCTTAGCCC-3', 'TAGTTGGGTTAAAGAC-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g13_recip, g13_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")

#14
Recip14 =  c('AAATGCCGTCGGGTCT-1', 'AACTGGTAGTTACGGG-1', 'AAGTCTGGTCCGAACC-1', 'ACGCAGCCAACTGGCC-1', 'AGGTCCGCAGGTGGAT-1', 'AGTGTCAGTACGACCC-1', 'ATCATGGTCAACCAAC-1', 'ATTCTACTCTTTAGGG-1', 'CCTTACGCATCTGGTA-1', 'CTCTGGTCAAGGTTTC-1', 'GAATGAAGTCTCACCT-1', 'GCAGTTAGTCTTGCGG-1', 'GCTTGAATCGACCAGC-1', 'GTATTCTCATGGTCTA-1', 'TAAACCGGTTGGTGGA-1', 'TACAGTGAGATAGTCA-1', 'AAATGCCAGAAGGGTA-2', 'ACTGAACCAAGACGTG-2', 'AGATTGCGTACCGTAT-2', 'CATGCCTAGGAGTAGA-2', 'CCTATTAAGCGCTTAT-2', 'CGCTGGATCCTTGACC-2', 'CGGACACTCGGCGGTT-2', 'CGTGAGCAGTCAATAG-2', 'CTACACCTCGTCACGG-2', 'GAAATGACAGTATGCT-2', 'GATGAAAAGAGATGAG-2', 'GCACATATCTCTAAGG-2', 'GCACTCTCAGACGTAG-2', 'GCTTCCAGTTGTCTTT-2', 'GCTTGAAAGTGGCACA-2', 'GGATGTTTCTTCATGT-2', 'TACACGATCTACTTAC-2', 'TAGGCATTCTCTTGAT-2', 'TATCTCATCGGAATCT-2', 'TCACAAGCACCACCAG-2', 'TCATTACAGATACACA-2', 'TCGCGTTAGACTAGAT-2', 'TCGGTAACACAGACAG-2', 'TCTTCGGAGCCTCGTG-2', 'TGCCAAACACACAGAG-2', 'TTAACTCTCCGATATG-2', 'ACACCGGGTGTTGAGG-3', 'ACACTGACACCCATGG-3', 'ACAGCCGAGATAGGAG-3', 'ACATCAGGTAAAGGAG-3', 'ATGGGAGCAATGGTCT-3', 'CTACGTCGTACATCCA-3', 'CTGAAGTGTTATCGGT-3', 'GACTACAGTGATGTCT-3', 'GATGCTAAGTACACCT-3', 'GCGACCATCAGTCAGT-3', 'GGACATTTCCCAGGTG-3', 'GGGTTGCAGTCTCCTC-3', 'TATGCCCGTCTCTTAT-3', 'TGGTTAGCACGGTAGA-3', 'TTGGAACCAATCCAAC-3', 'TTGGCAAGTGACTACT-3', 'TGAGAGGTCCCTTGCA-3')
Don14 =  c('CTAAGACTCATCGATG-2', 'GACGGCTTCCACGCAG-2',  'AAACGGGTCGCAGGCT-3', 'ACCGTAACAAGGTGTG-3', 'CATTATCTCCAACCAA-3', 'CATTATCTCGCCATAA-3', 'CGAACATGTAAACACA-3', 'CGCCAAGGTAGCGCTC-3', 'CGTCTACCACCTGGTG-3', 'CTCTACGAGAGCAATT-3', 'GATCGCGCATCCTTGC-3', 'TCAACGAAGTGAAGAG-3', 'TGAGGGAGTCTCAACA-3', 'GTCGGGTCAGGTGGAT-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone14 <- subset(CD8.clean3, ident = c('AAATGCCGTCGGGTCT-1', 'AACTGGTAGTTACGGG-1', 'AAGTCTGGTCCGAACC-1', 'ACGCAGCCAACTGGCC-1', 'AGGTCCGCAGGTGGAT-1', 'AGTGTCAGTACGACCC-1', 'ATCATGGTCAACCAAC-1', 'ATTCTACTCTTTAGGG-1', 'CCTTACGCATCTGGTA-1', 'CTCTGGTCAAGGTTTC-1', 'GAATGAAGTCTCACCT-1', 'GCAGTTAGTCTTGCGG-1', 'GCTTGAATCGACCAGC-1', 'GTATTCTCATGGTCTA-1', 'TAAACCGGTTGGTGGA-1', 'TACAGTGAGATAGTCA-1', 'AAATGCCAGAAGGGTA-2', 'ACTGAACCAAGACGTG-2', 'AGATTGCGTACCGTAT-2', 'CATGCCTAGGAGTAGA-2', 'CCTATTAAGCGCTTAT-2', 'CGCTGGATCCTTGACC-2', 'CGGACACTCGGCGGTT-2', 'CGTGAGCAGTCAATAG-2', 'CTACACCTCGTCACGG-2', 'GAAATGACAGTATGCT-2', 'GATGAAAAGAGATGAG-2', 'GCACATATCTCTAAGG-2', 'GCACTCTCAGACGTAG-2', 'GCTTCCAGTTGTCTTT-2', 'GCTTGAAAGTGGCACA-2', 'GGATGTTTCTTCATGT-2', 'TACACGATCTACTTAC-2', 'TAGGCATTCTCTTGAT-2', 'TATCTCATCGGAATCT-2', 'TCACAAGCACCACCAG-2', 'TCATTACAGATACACA-2', 'TCGCGTTAGACTAGAT-2', 'TCGGTAACACAGACAG-2', 'TCTTCGGAGCCTCGTG-2', 'TGCCAAACACACAGAG-2', 'TTAACTCTCCGATATG-2', 'ACACCGGGTGTTGAGG-3', 'ACACTGACACCCATGG-3', 'ACAGCCGAGATAGGAG-3', 'ACATCAGGTAAAGGAG-3', 'ATGGGAGCAATGGTCT-3', 'CTACGTCGTACATCCA-3', 'CTGAAGTGTTATCGGT-3', 'GACTACAGTGATGTCT-3', 'GATGCTAAGTACACCT-3', 'GCGACCATCAGTCAGT-3', 'GGACATTTCCCAGGTG-3', 'GGGTTGCAGTCTCCTC-3', 'TATGCCCGTCTCTTAT-3', 'TGGTTAGCACGGTAGA-3', 'TTGGAACCAATCCAAC-3', 'TTGGCAAGTGACTACT-3', 'TGAGAGGTCCCTTGCA-3', 'AAACGGGTCGCAGGCT-3', 'ACCGTAACAAGGTGTG-3', 'CATTATCTCCAACCAA-3', 'CATTATCTCGCCATAA-3', 'CGAACATGTAAACACA-3', 'CGCCAAGGTAGCGCTC-3', 'CGTCTACCACCTGGTG-3', 'CTCTACGAGAGCAATT-3', 'GATCGCGCATCCTTGC-3', 'TCAACGAAGTGAAGAG-3', 'TGAGGGAGTCTCAACA-3', 'GTCGGGTCAGGTGGAT-3'), invert = F)

Clone14 <- SetIdent(Clone14, value = "Genotype")
markers.clone14 <- FindAllMarkers(Clone14,
                                  features = vargenes.filtered5,
                                  only.pos = FALSE, 
                                  min.pct = 0.01, 
                                  logfc.threshold = 0.25,
                                  test.use = "MAST",
                                  verbose = T,
                                  latent.vars = "nCount_RNA")


markers.clone14top20 <- markers.clone14 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone14, 
                     features = markers.clone14top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g14_recip <- WhichCells(CD8.clean3, idents = c('AAATGCCGTCGGGTCT-1', 'AACTGGTAGTTACGGG-1', 'AAGTCTGGTCCGAACC-1', 'ACGCAGCCAACTGGCC-1', 'AGGTCCGCAGGTGGAT-1', 'AGTGTCAGTACGACCC-1', 'ATCATGGTCAACCAAC-1', 'ATTCTACTCTTTAGGG-1', 'CCTTACGCATCTGGTA-1', 'CTCTGGTCAAGGTTTC-1', 'GAATGAAGTCTCACCT-1', 'GCAGTTAGTCTTGCGG-1', 'GCTTGAATCGACCAGC-1', 'GTATTCTCATGGTCTA-1', 'TAAACCGGTTGGTGGA-1', 'TACAGTGAGATAGTCA-1', 'AAATGCCAGAAGGGTA-2', 'ACTGAACCAAGACGTG-2', 'AGATTGCGTACCGTAT-2', 'CATGCCTAGGAGTAGA-2', 'CCTATTAAGCGCTTAT-2', 'CGCTGGATCCTTGACC-2', 'CGGACACTCGGCGGTT-2', 'CGTGAGCAGTCAATAG-2', 'CTACACCTCGTCACGG-2', 'GAAATGACAGTATGCT-2', 'GATGAAAAGAGATGAG-2', 'GCACATATCTCTAAGG-2', 'GCACTCTCAGACGTAG-2', 'GCTTCCAGTTGTCTTT-2', 'GCTTGAAAGTGGCACA-2', 'GGATGTTTCTTCATGT-2', 'TACACGATCTACTTAC-2', 'TAGGCATTCTCTTGAT-2', 'TATCTCATCGGAATCT-2', 'TCACAAGCACCACCAG-2', 'TCATTACAGATACACA-2', 'TCGCGTTAGACTAGAT-2', 'TCGGTAACACAGACAG-2', 'TCTTCGGAGCCTCGTG-2', 'TGCCAAACACACAGAG-2', 'TTAACTCTCCGATATG-2', 'ACACCGGGTGTTGAGG-3', 'ACACTGACACCCATGG-3', 'ACAGCCGAGATAGGAG-3', 'ACATCAGGTAAAGGAG-3', 'ATGGGAGCAATGGTCT-3', 'CTACGTCGTACATCCA-3', 'CTGAAGTGTTATCGGT-3', 'GACTACAGTGATGTCT-3', 'GATGCTAAGTACACCT-3', 'GCGACCATCAGTCAGT-3', 'GGACATTTCCCAGGTG-3', 'GGGTTGCAGTCTCCTC-3', 'TATGCCCGTCTCTTAT-3', 'TGGTTAGCACGGTAGA-3', 'TTGGAACCAATCCAAC-3', 'TTGGCAAGTGACTACT-3', 'TGAGAGGTCCCTTGCA-3'))
g14_don <- WhichCells(CD8.clean3, idents = c('AAACGGGTCGCAGGCT-3', 'ACCGTAACAAGGTGTG-3', 'CATTATCTCCAACCAA-3', 'CATTATCTCGCCATAA-3', 'CGAACATGTAAACACA-3', 'CGCCAAGGTAGCGCTC-3', 'CGTCTACCACCTGGTG-3', 'CTCTACGAGAGCAATT-3', 'GATCGCGCATCCTTGC-3', 'TCAACGAAGTGAAGAG-3', 'TGAGGGAGTCTCAACA-3', 'GTCGGGTCAGGTGGAT-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g14_recip, g14_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")


#15
Recip15 =  c('ACGGCCAAGAGAACAG-1', 'AGCATACGTGTTCTTT-1', 'ATAAGAGAGCTGAACG-1', 'CAGAATCCATATACCG-1', 'CATTATCCACTAAGTC-1', 'GGGCACTCAGCTGTGC-1', 'GTTCATTCAAGGACTG-1', 'TCACAAGTCAGAGCTT-1', 'TCATTACCAGACACTT-1', 'TCATTACGTCACTGGC-1', 'TCGGGACAGCCGTCGT-1', 'GGGTTGCTCGCCATAA-1', 'ACCGTAAGTTACGCGC-2', 'ACGTCAAAGGTGACCA-2', 'CGGCTAGAGGCATGTG-2', 'GCGCAGTGTACAAGTA-2', 'GGAATAACAAATCCGT-2', 'GTGTGCGGTTGTTTGG-2', 'TAGCCGGCAAGACACG-2', 'TTGTAGGGTCAGAATA-2', 'AAGGCAGCAGAAGCAC-3', 'ACGATACTCATACGGT-3', 'CATCAAGAGTATGACA-3', 'CGTGAGCAGGGTATCG-3', 'GGACAGAAGTGAAGTT-3', 'TGGCCAGGTGACTCAT-3')
Don15 =  c('AACCGCGTCCTGCTTG-3', 'ACGGGCTGTCTCACCT-3', 'ATGGGAGAGGTGCAAC-3', 'ATTGGACTCTTACCGC-3', 'CAGATCAAGAAGGTTT-3', 'CAGCAGCTCTCTTATG-3', 'CAGCGACCAGGATTGG-3', 'CATCCACCAATAGAGT-3', 'CGTTGGGTCGTAGATC-3', 'CTGAAACGTCGTGGCT-3', 'CTGAAACTCTGTCTAT-3', 'CTTGGCTCAAGCCCAC-3', 'GAAATGACAAGCCTAT-3', 'GCACTCTTCGCAAGCC-3', 'GCAGTTAGTTGAGTTC-3', 'GCATGCGTCCTTTACA-3', 'GCGAGAACAGATCCAT-3', 'GTCCTCAGTCGGGTCT-3', 'GTCGGGTCAGGGTATG-3', 'GTGCATACAGCATACT-3', 'GTTACAGGTATCGCAT-3', 'TAAGCGTCAAGGTTCT-3', 'TGGACGCAGAGCAATT-3', 'TTTGGTTTCAGATAAG-3', 'CAACTAGCACATCCAA-3', 'ACTGATGTCTTAGCCC-3', 'TAGTTGGGTTAAAGAC-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone15 <- subset(CD8.clean3, ident = c('ACGGCCAAGAGAACAG-1', 'AGCATACGTGTTCTTT-1', 'ATAAGAGAGCTGAACG-1', 'CAGAATCCATATACCG-1', 'CATTATCCACTAAGTC-1', 'GGGCACTCAGCTGTGC-1', 'GTTCATTCAAGGACTG-1', 'TCACAAGTCAGAGCTT-1', 'TCATTACCAGACACTT-1', 'TCATTACGTCACTGGC-1', 'TCGGGACAGCCGTCGT-1', 'GGGTTGCTCGCCATAA-1', 'ACCGTAAGTTACGCGC-2', 'ACGTCAAAGGTGACCA-2', 'CGGCTAGAGGCATGTG-2', 'GCGCAGTGTACAAGTA-2', 'GGAATAACAAATCCGT-2', 'GTGTGCGGTTGTTTGG-2', 'TAGCCGGCAAGACACG-2', 'TTGTAGGGTCAGAATA-2', 'AAGGCAGCAGAAGCAC-3', 'ACGATACTCATACGGT-3', 'CATCAAGAGTATGACA-3', 'CGTGAGCAGGGTATCG-3', 'GGACAGAAGTGAAGTT-3', 'TGGCCAGGTGACTCAT-3', 'AACCGCGTCCTGCTTG-3', 'ACGGGCTGTCTCACCT-3', 'ATGGGAGAGGTGCAAC-3', 'ATTGGACTCTTACCGC-3', 'CAGATCAAGAAGGTTT-3', 'CAGCAGCTCTCTTATG-3', 'CAGCGACCAGGATTGG-3', 'CATCCACCAATAGAGT-3', 'CGTTGGGTCGTAGATC-3', 'CTGAAACGTCGTGGCT-3', 'CTGAAACTCTGTCTAT-3', 'CTTGGCTCAAGCCCAC-3', 'GAAATGACAAGCCTAT-3', 'GCACTCTTCGCAAGCC-3', 'GCAGTTAGTTGAGTTC-3', 'GCATGCGTCCTTTACA-3', 'GCGAGAACAGATCCAT-3', 'GTCCTCAGTCGGGTCT-3', 'GTCGGGTCAGGGTATG-3', 'GTGCATACAGCATACT-3', 'GTTACAGGTATCGCAT-3', 'TAAGCGTCAAGGTTCT-3', 'TGGACGCAGAGCAATT-3', 'TTTGGTTTCAGATAAG-3', 'CAACTAGCACATCCAA-3', 'ACTGATGTCTTAGCCC-3', 'TAGTTGGGTTAAAGAC-3'), invert = F)

Clone15 <- SetIdent(Clone15, value = "Genotype")
markers.clone15 <- FindAllMarkers(Clone15,
                                  features = vargenes.filtered5,
                                  only.pos = FALSE, 
                                  min.pct = 0.01, 
                                  logfc.threshold = 0.25,
                                  test.use = "MAST",
                                  verbose = T,
                                  latent.vars = "nCount_RNA")


markers.clone15top20 <- markers.clone15 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone15, 
                     features = markers.clone2top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g15_recip <- WhichCells(CD8.clean3, idents = c('ACGGCCAAGAGAACAG-1', 'AGCATACGTGTTCTTT-1', 'ATAAGAGAGCTGAACG-1', 'CAGAATCCATATACCG-1', 'CATTATCCACTAAGTC-1', 'GGGCACTCAGCTGTGC-1', 'GTTCATTCAAGGACTG-1', 'TCACAAGTCAGAGCTT-1', 'TCATTACCAGACACTT-1', 'TCATTACGTCACTGGC-1', 'TCGGGACAGCCGTCGT-1', 'GGGTTGCTCGCCATAA-1', 'ACCGTAAGTTACGCGC-2', 'ACGTCAAAGGTGACCA-2', 'CGGCTAGAGGCATGTG-2', 'GCGCAGTGTACAAGTA-2', 'GGAATAACAAATCCGT-2', 'GTGTGCGGTTGTTTGG-2', 'TAGCCGGCAAGACACG-2', 'TTGTAGGGTCAGAATA-2', 'AAGGCAGCAGAAGCAC-3', 'ACGATACTCATACGGT-3', 'CATCAAGAGTATGACA-3', 'CGTGAGCAGGGTATCG-3', 'GGACAGAAGTGAAGTT-3', 'TGGCCAGGTGACTCAT-3'))
g15_don <- WhichCells(CD8.clean3, idents = c('AACCGCGTCCTGCTTG-3', 'ACGGGCTGTCTCACCT-3', 'ATGGGAGAGGTGCAAC-3', 'ATTGGACTCTTACCGC-3', 'CAGATCAAGAAGGTTT-3', 'CAGCAGCTCTCTTATG-3', 'CAGCGACCAGGATTGG-3', 'CATCCACCAATAGAGT-3', 'CGTTGGGTCGTAGATC-3', 'CTGAAACGTCGTGGCT-3', 'CTGAAACTCTGTCTAT-3', 'CTTGGCTCAAGCCCAC-3', 'GAAATGACAAGCCTAT-3', 'GCACTCTTCGCAAGCC-3', 'GCAGTTAGTTGAGTTC-3', 'GCATGCGTCCTTTACA-3', 'GCGAGAACAGATCCAT-3', 'GTCCTCAGTCGGGTCT-3', 'GTCGGGTCAGGGTATG-3', 'GTGCATACAGCATACT-3', 'GTTACAGGTATCGCAT-3', 'TAAGCGTCAAGGTTCT-3', 'TGGACGCAGAGCAATT-3', 'TTTGGTTTCAGATAAG-3', 'CAACTAGCACATCCAA-3', 'ACTGATGTCTTAGCCC-3', 'TAGTTGGGTTAAAGAC-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g15_recip, g15_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")

#16
Recip16 =  c('ACGGCCAAGAGAACAG-1', 'AGCATACGTGTTCTTT-1', 'ATAAGAGAGCTGAACG-1', 'CAGAATCCATATACCG-1', 'CATTATCCACTAAGTC-1', 'GGGCACTCAGCTGTGC-1', 'GTTCATTCAAGGACTG-1', 'TCACAAGTCAGAGCTT-1', 'TCATTACCAGACACTT-1', 'TCATTACGTCACTGGC-1', 'TCGGGACAGCCGTCGT-1', 'GGGTTGCTCGCCATAA-1',  'ACCGTAAGTTACGCGC-2', 'ACGTCAAAGGTGACCA-2', 'CGGCTAGAGGCATGTG-2', 'GCGCAGTGTACAAGTA-2', 'GGAATAACAAATCCGT-2', 'GTGTGCGGTTGTTTGG-2', 'TAGCCGGCAAGACACG-2', 'TTGTAGGGTCAGAATA-2', 'ACCGTAAGTTACGCGC-2', 'ACGTCAAAGGTGACCA-2', 'CGGCTAGAGGCATGTG-2', 'GCGCAGTGTACAAGTA-2', 'GGAATAACAAATCCGT-2', 'GTGTGCGGTTGTTTGG-2', 'TAGCCGGCAAGACACG-2', 'TTGTAGGGTCAGAATA-2',  'AAGGCAGCAGAAGCAC-3', 'ACGATACTCATACGGT-3', 'CATCAAGAGTATGACA-3', 'CGTGAGCAGGGTATCG-3', 'GGACAGAAGTGAAGTT-3', 'TGGCCAGGTGACTCAT-3')
Don16 =  c('AAAGCAAGTCATTAGC-3', 'ACGGGCTAGGACTGGT-3', 'ACTATCTGTTGAACTC-3', 'AGCATACCAGGATTGG-3', 'AGCGGTCAGGCAAAGA-3', 'AGGCCACAGCTGTTCA-3', 'ATCATCTGTTCAGTAC-3', 'ATCTGCCCATGGGACA-3', 'ATTATCCTCTGATTCT-3', 'ATTGGTGTCAAAGACA-3', 'ATTTCTGAGACCCACC-3', 'CAAGGCCAGTAGTGCG-3', 'CAAGGCCCAGCTGCTG-3', 'CACACCTCATGGAATA-3', 'CACTCCAAGAAGAAGC-3', 'CCCAGTTTCTAACGGT-3', 'CCTAAAGGTGTCTGAT-3', 'CGATGGCGTGAGTATA-3', 'CGTCCATCAAGCCTAT-3', 'CGTTCTGAGCTAAACA-3', 'CTGAAGTGTGACGCCT-3', 'CTGATAGAGGGCTTGA-3', 'CTGGTCTTCTGCCCTA-3', 'GAACCTAAGGACACCA-3', 'GATGAGGTCTTATCTG-3', 'GATGCTATCACTATTC-3', 'GCCTCTAAGCTCCTCT-3', 'GCTCTGTAGACATAAC-3', 'GGACAGACAAGCGATG-3', 'GGCGTGTAGGTGCTTT-3', 'GTAGTCAGTTCCACAA-3', 'GTGGGTCTCAGTACGT-3', 'GTTACAGAGCGCTTAT-3', 'TAGACCAAGGAACTGC-3', 'TCGTAGAGTGAAGGCT-3', 'TGATTTCCACCGAAAG-3', 'TTAGGCAGTCCTCTTG-3', 'TTCTCAAGTCATACTG-3', 'TTGGAACGTTGCGTTA-3', 'TTTGTCAAGGCTAGGT-3', 'CACCACTGTAAAGGAG-3', 'CAACCTCGTTCGTTGA-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone16 <- subset(CD8.clean3, ident = c('ACGGCCAAGAGAACAG-1', 'AGCATACGTGTTCTTT-1', 'ATAAGAGAGCTGAACG-1', 'CAGAATCCATATACCG-1', 'CATTATCCACTAAGTC-1', 'GGGCACTCAGCTGTGC-1', 'GTTCATTCAAGGACTG-1', 'TCACAAGTCAGAGCTT-1', 'TCATTACCAGACACTT-1', 'TCATTACGTCACTGGC-1', 'TCGGGACAGCCGTCGT-1', 'GGGTTGCTCGCCATAA-1',  'ACCGTAAGTTACGCGC-2', 'ACGTCAAAGGTGACCA-2', 'CGGCTAGAGGCATGTG-2', 'GCGCAGTGTACAAGTA-2', 'GGAATAACAAATCCGT-2', 'GTGTGCGGTTGTTTGG-2', 'TAGCCGGCAAGACACG-2', 'TTGTAGGGTCAGAATA-2', 'ACCGTAAGTTACGCGC-2', 'ACGTCAAAGGTGACCA-2', 'CGGCTAGAGGCATGTG-2', 'GCGCAGTGTACAAGTA-2', 'GGAATAACAAATCCGT-2', 'GTGTGCGGTTGTTTGG-2', 'TAGCCGGCAAGACACG-2', 'TTGTAGGGTCAGAATA-2',  'AAGGCAGCAGAAGCAC-3', 'ACGATACTCATACGGT-3', 'CATCAAGAGTATGACA-3', 'CGTGAGCAGGGTATCG-3', 'GGACAGAAGTGAAGTT-3', 'TGGCCAGGTGACTCAT-3', 'AAAGCAAGTCATTAGC-3', 'ACGGGCTAGGACTGGT-3', 'ACTATCTGTTGAACTC-3', 'AGCATACCAGGATTGG-3', 'AGCGGTCAGGCAAAGA-3', 'AGGCCACAGCTGTTCA-3', 'ATCATCTGTTCAGTAC-3', 'ATCTGCCCATGGGACA-3', 'ATTATCCTCTGATTCT-3', 'ATTGGTGTCAAAGACA-3', 'ATTTCTGAGACCCACC-3', 'CAAGGCCAGTAGTGCG-3', 'CAAGGCCCAGCTGCTG-3', 'CACACCTCATGGAATA-3', 'CACTCCAAGAAGAAGC-3', 'CCCAGTTTCTAACGGT-3', 'CCTAAAGGTGTCTGAT-3', 'CGATGGCGTGAGTATA-3', 'CGTCCATCAAGCCTAT-3', 'CGTTCTGAGCTAAACA-3', 'CTGAAGTGTGACGCCT-3', 'CTGATAGAGGGCTTGA-3', 'CTGGTCTTCTGCCCTA-3', 'GAACCTAAGGACACCA-3', 'GATGAGGTCTTATCTG-3', 'GATGCTATCACTATTC-3', 'GCCTCTAAGCTCCTCT-3', 'GCTCTGTAGACATAAC-3', 'GGACAGACAAGCGATG-3', 'GGCGTGTAGGTGCTTT-3', 'GTAGTCAGTTCCACAA-3', 'GTGGGTCTCAGTACGT-3', 'GTTACAGAGCGCTTAT-3', 'TAGACCAAGGAACTGC-3', 'TCGTAGAGTGAAGGCT-3', 'TGATTTCCACCGAAAG-3', 'TTAGGCAGTCCTCTTG-3', 'TTCTCAAGTCATACTG-3', 'TTGGAACGTTGCGTTA-3', 'TTTGTCAAGGCTAGGT-3', 'CACCACTGTAAAGGAG-3', 'CAACCTCGTTCGTTGA-3'), invert = F)

Clone16 <- SetIdent(Clone16, value = "Genotype")
markers.clone16 <- FindAllMarkers(Clone16,
                                  features = vargenes.filtered5,
                                  only.pos = FALSE, 
                                  min.pct = 0.01, 
                                  logfc.threshold = 0.25,
                                  test.use = "MAST",
                                  verbose = T,
                                  latent.vars = "nCount_RNA")


markers.clone16top20 <- markers.clone16 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone16, 
                     features = markers.clone16top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g16_recip <- WhichCells(CD8.clean3, idents = c('ACGGCCAAGAGAACAG-1', 'AGCATACGTGTTCTTT-1', 'ATAAGAGAGCTGAACG-1', 'CAGAATCCATATACCG-1', 'CATTATCCACTAAGTC-1', 'GGGCACTCAGCTGTGC-1', 'GTTCATTCAAGGACTG-1', 'TCACAAGTCAGAGCTT-1', 'TCATTACCAGACACTT-1', 'TCATTACGTCACTGGC-1', 'TCGGGACAGCCGTCGT-1', 'GGGTTGCTCGCCATAA-1',  'ACCGTAAGTTACGCGC-2', 'ACGTCAAAGGTGACCA-2', 'CGGCTAGAGGCATGTG-2', 'GCGCAGTGTACAAGTA-2', 'GGAATAACAAATCCGT-2', 'GTGTGCGGTTGTTTGG-2', 'TAGCCGGCAAGACACG-2', 'TTGTAGGGTCAGAATA-2', 'ACCGTAAGTTACGCGC-2', 'ACGTCAAAGGTGACCA-2', 'CGGCTAGAGGCATGTG-2', 'GCGCAGTGTACAAGTA-2', 'GGAATAACAAATCCGT-2', 'GTGTGCGGTTGTTTGG-2', 'TAGCCGGCAAGACACG-2', 'TTGTAGGGTCAGAATA-2',  'AAGGCAGCAGAAGCAC-3', 'ACGATACTCATACGGT-3', 'CATCAAGAGTATGACA-3', 'CGTGAGCAGGGTATCG-3', 'GGACAGAAGTGAAGTT-3', 'TGGCCAGGTGACTCAT-3'))
g16_don <- WhichCells(CD8.clean3, idents = c('AAAGCAAGTCATTAGC-3', 'ACGGGCTAGGACTGGT-3', 'ACTATCTGTTGAACTC-3', 'AGCATACCAGGATTGG-3', 'AGCGGTCAGGCAAAGA-3', 'AGGCCACAGCTGTTCA-3', 'ATCATCTGTTCAGTAC-3', 'ATCTGCCCATGGGACA-3', 'ATTATCCTCTGATTCT-3', 'ATTGGTGTCAAAGACA-3', 'ATTTCTGAGACCCACC-3', 'CAAGGCCAGTAGTGCG-3', 'CAAGGCCCAGCTGCTG-3', 'CACACCTCATGGAATA-3', 'CACTCCAAGAAGAAGC-3', 'CCCAGTTTCTAACGGT-3', 'CCTAAAGGTGTCTGAT-3', 'CGATGGCGTGAGTATA-3', 'CGTCCATCAAGCCTAT-3', 'CGTTCTGAGCTAAACA-3', 'CTGAAGTGTGACGCCT-3', 'CTGATAGAGGGCTTGA-3', 'CTGGTCTTCTGCCCTA-3', 'GAACCTAAGGACACCA-3', 'GATGAGGTCTTATCTG-3', 'GATGCTATCACTATTC-3', 'GCCTCTAAGCTCCTCT-3', 'GCTCTGTAGACATAAC-3', 'GGACAGACAAGCGATG-3', 'GGCGTGTAGGTGCTTT-3', 'GTAGTCAGTTCCACAA-3', 'GTGGGTCTCAGTACGT-3', 'GTTACAGAGCGCTTAT-3', 'TAGACCAAGGAACTGC-3', 'TCGTAGAGTGAAGGCT-3', 'TGATTTCCACCGAAAG-3', 'TTAGGCAGTCCTCTTG-3', 'TTCTCAAGTCATACTG-3', 'TTGGAACGTTGCGTTA-3', 'TTTGTCAAGGCTAGGT-3', 'CACCACTGTAAAGGAG-3', 'CAACCTCGTTCGTTGA-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g16_recip, g16_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")

#17
Recip17 =  c('AAGGAGCCACTCTGTC-1', 'ATGAGGGCATTGGGCC-1', 'ATGCGATTCGGTCTAA-1', 'ATTCTACTCTTGTACT-1', 'ATTTCTGCATCACAAC-1', 'CATCAAGTCCACTCCA-1', 'GGACAGAAGCAACGGT-1', 'GGGTCTGCACAGTCGC-1', 'TCTCTAACACGAAGCA-1', 'TTGGAACAGCTAGTCT-1', 'AAGCCGCGTTTAGGAA-2', 'ACACCCTCATCGATTG-2', 'ACACTGATCGGCGGTT-2', 'CCTCTGACAATAGCAA-2', 'CGGTTAAAGTGCGTGA-2', 'GTTTCTAAGATCTGAA-2', 'TACTTACCAGCAGTTT-2',  'AAAGATGTCCTCTAGC-3', 'AGTTGGTCAAGAAGAG-3', 'CAGCAGCAGACGCACA-3', 'CTCTACGCATGCCCGA-3', 'GGTGCGTGTGTTCTTT-3', 'GTAGTCATCCTTGGTC-3', 'GTTCGGGTCCGAGCCA-3', 'TGGGAAGAGACAGACC-3')
Don17 =  c('AACCGCGTCCTGCTTG-3', 'ACGGGCTGTCTCACCT-3', 'ATGGGAGAGGTGCAAC-3', 'ATTGGACTCTTACCGC-3', 'CAGATCAAGAAGGTTT-3', 'CAGCAGCTCTCTTATG-3', 'CAGCGACCAGGATTGG-3', 'CATCCACCAATAGAGT-3', 'CGTTGGGTCGTAGATC-3', 'CTGAAACGTCGTGGCT-3', 'CTGAAACTCTGTCTAT-3', 'CTTGGCTCAAGCCCAC-3', 'GAAATGACAAGCCTAT-3', 'GCACTCTTCGCAAGCC-3', 'GCAGTTAGTTGAGTTC-3', 'GCATGCGTCCTTTACA-3', 'GCGAGAACAGATCCAT-3', 'GTCCTCAGTCGGGTCT-3', 'GTCGGGTCAGGGTATG-3', 'GTGCATACAGCATACT-3', 'GTTACAGGTATCGCAT-3', 'TAAGCGTCAAGGTTCT-3', 'TGGACGCAGAGCAATT-3', 'TTTGGTTTCAGATAAG-3', 'CAACTAGCACATCCAA-3', 'ACTGATGTCTTAGCCC-3', 'TAGTTGGGTTAAAGAC-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone17 <- subset(CD8.clean3, ident = c('AAGGAGCCACTCTGTC-1', 'ATGAGGGCATTGGGCC-1', 'ATGCGATTCGGTCTAA-1', 'ATTCTACTCTTGTACT-1', 'ATTTCTGCATCACAAC-1', 'CATCAAGTCCACTCCA-1', 'GGACAGAAGCAACGGT-1', 'GGGTCTGCACAGTCGC-1', 'TCTCTAACACGAAGCA-1', 'TTGGAACAGCTAGTCT-1', 'AAGCCGCGTTTAGGAA-2', 'ACACCCTCATCGATTG-2', 'ACACTGATCGGCGGTT-2', 'CCTCTGACAATAGCAA-2', 'CGGTTAAAGTGCGTGA-2', 'GTTTCTAAGATCTGAA-2', 'TACTTACCAGCAGTTT-2',  'AAAGATGTCCTCTAGC-3', 'AGTTGGTCAAGAAGAG-3', 'CAGCAGCAGACGCACA-3', 'CTCTACGCATGCCCGA-3', 'GGTGCGTGTGTTCTTT-3', 'GTAGTCATCCTTGGTC-3', 'GTTCGGGTCCGAGCCA-3', 'TGGGAAGAGACAGACC-3', 'AACCGCGTCCTGCTTG-3', 'ACGGGCTGTCTCACCT-3', 'ATGGGAGAGGTGCAAC-3', 'ATTGGACTCTTACCGC-3', 'CAGATCAAGAAGGTTT-3', 'CAGCAGCTCTCTTATG-3', 'CAGCGACCAGGATTGG-3', 'CATCCACCAATAGAGT-3', 'CGTTGGGTCGTAGATC-3', 'CTGAAACGTCGTGGCT-3', 'CTGAAACTCTGTCTAT-3', 'CTTGGCTCAAGCCCAC-3', 'GAAATGACAAGCCTAT-3', 'GCACTCTTCGCAAGCC-3', 'GCAGTTAGTTGAGTTC-3', 'GCATGCGTCCTTTACA-3', 'GCGAGAACAGATCCAT-3', 'GTCCTCAGTCGGGTCT-3', 'GTCGGGTCAGGGTATG-3', 'GTGCATACAGCATACT-3', 'GTTACAGGTATCGCAT-3', 'TAAGCGTCAAGGTTCT-3', 'TGGACGCAGAGCAATT-3', 'TTTGGTTTCAGATAAG-3', 'CAACTAGCACATCCAA-3', 'ACTGATGTCTTAGCCC-3', 'TAGTTGGGTTAAAGAC-3'), invert = F)

Clone17 <- SetIdent(Clone17, value = "Genotype")
markers.clone17 <- FindAllMarkers(Clone17,
                                  features = vargenes.filtered5,
                                  only.pos = FALSE, 
                                  min.pct = 0.01, 
                                  logfc.threshold = 0.25,
                                  test.use = "MAST",
                                  verbose = T,
                                  latent.vars = "nCount_RNA")


markers.clone17top20 <- markers.clone17 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone17, 
                     features = markers.clone17top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g17_recip <- WhichCells(CD8.clean3, idents = c('AAGGAGCCACTCTGTC-1', 'ATGAGGGCATTGGGCC-1', 'ATGCGATTCGGTCTAA-1', 'ATTCTACTCTTGTACT-1', 'ATTTCTGCATCACAAC-1', 'CATCAAGTCCACTCCA-1', 'GGACAGAAGCAACGGT-1', 'GGGTCTGCACAGTCGC-1', 'TCTCTAACACGAAGCA-1', 'TTGGAACAGCTAGTCT-1', 'AAGCCGCGTTTAGGAA-2', 'ACACCCTCATCGATTG-2', 'ACACTGATCGGCGGTT-2', 'CCTCTGACAATAGCAA-2', 'CGGTTAAAGTGCGTGA-2', 'GTTTCTAAGATCTGAA-2', 'TACTTACCAGCAGTTT-2',  'AAAGATGTCCTCTAGC-3', 'AGTTGGTCAAGAAGAG-3', 'CAGCAGCAGACGCACA-3', 'CTCTACGCATGCCCGA-3', 'GGTGCGTGTGTTCTTT-3', 'GTAGTCATCCTTGGTC-3', 'GTTCGGGTCCGAGCCA-3', 'TGGGAAGAGACAGACC-3'))
g17_don <- WhichCells(CD8.clean3, idents = c('AACCGCGTCCTGCTTG-3', 'ACGGGCTGTCTCACCT-3', 'ATGGGAGAGGTGCAAC-3', 'ATTGGACTCTTACCGC-3', 'CAGATCAAGAAGGTTT-3', 'CAGCAGCTCTCTTATG-3', 'CAGCGACCAGGATTGG-3', 'CATCCACCAATAGAGT-3', 'CGTTGGGTCGTAGATC-3', 'CTGAAACGTCGTGGCT-3', 'CTGAAACTCTGTCTAT-3', 'CTTGGCTCAAGCCCAC-3', 'GAAATGACAAGCCTAT-3', 'GCACTCTTCGCAAGCC-3', 'GCAGTTAGTTGAGTTC-3', 'GCATGCGTCCTTTACA-3', 'GCGAGAACAGATCCAT-3', 'GTCCTCAGTCGGGTCT-3', 'GTCGGGTCAGGGTATG-3', 'GTGCATACAGCATACT-3', 'GTTACAGGTATCGCAT-3', 'TAAGCGTCAAGGTTCT-3', 'TGGACGCAGAGCAATT-3', 'TTTGGTTTCAGATAAG-3', 'CAACTAGCACATCCAA-3', 'ACTGATGTCTTAGCCC-3', 'TAGTTGGGTTAAAGAC-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g17_recip, g17_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")

#18
Recip18 =  c('AAGGAGCCACTCTGTC-1', 'ATGAGGGCATTGGGCC-1', 'ATGCGATTCGGTCTAA-1', 'ATTCTACTCTTGTACT-1', 'ATTTCTGCATCACAAC-1', 'CATCAAGTCCACTCCA-1', 'GGACAGAAGCAACGGT-1', 'GGGTCTGCACAGTCGC-1', 'TCTCTAACACGAAGCA-1', 'TTGGAACAGCTAGTCT-1', 'AAGCCGCGTTTAGGAA-2', 'ACACCCTCATCGATTG-2', 'ACACTGATCGGCGGTT-2', 'CCTCTGACAATAGCAA-2', 'CGGTTAAAGTGCGTGA-2', 'GTTTCTAAGATCTGAA-2', 'TACTTACCAGCAGTTT-2', 'AAAGATGTCCTCTAGC-3', 'AGTTGGTCAAGAAGAG-3', 'CAGCAGCAGACGCACA-3', 'CTCTACGCATGCCCGA-3', 'GGTGCGTGTGTTCTTT-3', 'GTAGTCATCCTTGGTC-3', 'GTTCGGGTCCGAGCCA-3', 'TGGGAAGAGACAGACC-3')
Don18 =  c('AAAGCAAGTCATTAGC-3', 'ACGGGCTAGGACTGGT-3', 'ACTATCTGTTGAACTC-3', 'AGCATACCAGGATTGG-3', 'AGCGGTCAGGCAAAGA-3', 'AGGCCACAGCTGTTCA-3', 'ATCATCTGTTCAGTAC-3', 'ATCTGCCCATGGGACA-3', 'ATTATCCTCTGATTCT-3', 'ATTGGTGTCAAAGACA-3', 'ATTTCTGAGACCCACC-3', 'CAAGGCCAGTAGTGCG-3', 'CAAGGCCCAGCTGCTG-3', 'CACACCTCATGGAATA-3', 'CACTCCAAGAAGAAGC-3', 'CCCAGTTTCTAACGGT-3', 'CCTAAAGGTGTCTGAT-3', 'CGATGGCGTGAGTATA-3', 'CGTCCATCAAGCCTAT-3', 'CGTTCTGAGCTAAACA-3', 'CTGAAGTGTGACGCCT-3', 'CTGATAGAGGGCTTGA-3', 'CTGGTCTTCTGCCCTA-3', 'GAACCTAAGGACACCA-3', 'GATGAGGTCTTATCTG-3', 'GATGCTATCACTATTC-3', 'GCCTCTAAGCTCCTCT-3', 'GCTCTGTAGACATAAC-3', 'GGACAGACAAGCGATG-3', 'GGCGTGTAGGTGCTTT-3', 'GTAGTCAGTTCCACAA-3', 'GTGGGTCTCAGTACGT-3', 'GTTACAGAGCGCTTAT-3', 'TAGACCAAGGAACTGC-3', 'TCGTAGAGTGAAGGCT-3', 'TGATTTCCACCGAAAG-3', 'TTAGGCAGTCCTCTTG-3', 'TTCTCAAGTCATACTG-3', 'TTGGAACGTTGCGTTA-3', 'TTTGTCAAGGCTAGGT-3', 'CACCACTGTAAAGGAG-3', 'CAACCTCGTTCGTTGA-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone18 <- subset(CD8.clean3, ident = c('AAGGAGCCACTCTGTC-1', 'ATGAGGGCATTGGGCC-1', 'ATGCGATTCGGTCTAA-1', 'ATTCTACTCTTGTACT-1', 'ATTTCTGCATCACAAC-1', 'CATCAAGTCCACTCCA-1', 'GGACAGAAGCAACGGT-1', 'GGGTCTGCACAGTCGC-1', 'TCTCTAACACGAAGCA-1', 'TTGGAACAGCTAGTCT-1', 'AAGCCGCGTTTAGGAA-2', 'ACACCCTCATCGATTG-2', 'ACACTGATCGGCGGTT-2', 'CCTCTGACAATAGCAA-2', 'CGGTTAAAGTGCGTGA-2', 'GTTTCTAAGATCTGAA-2', 'TACTTACCAGCAGTTT-2', 'AAAGATGTCCTCTAGC-3', 'AGTTGGTCAAGAAGAG-3', 'CAGCAGCAGACGCACA-3', 'CTCTACGCATGCCCGA-3', 'GGTGCGTGTGTTCTTT-3', 'GTAGTCATCCTTGGTC-3', 'GTTCGGGTCCGAGCCA-3', 'TGGGAAGAGACAGACC-3', 'AAAGCAAGTCATTAGC-3', 'ACGGGCTAGGACTGGT-3', 'ACTATCTGTTGAACTC-3', 'AGCATACCAGGATTGG-3', 'AGCGGTCAGGCAAAGA-3', 'AGGCCACAGCTGTTCA-3', 'ATCATCTGTTCAGTAC-3', 'ATCTGCCCATGGGACA-3', 'ATTATCCTCTGATTCT-3', 'ATTGGTGTCAAAGACA-3', 'ATTTCTGAGACCCACC-3', 'CAAGGCCAGTAGTGCG-3', 'CAAGGCCCAGCTGCTG-3', 'CACACCTCATGGAATA-3', 'CACTCCAAGAAGAAGC-3', 'CCCAGTTTCTAACGGT-3', 'CCTAAAGGTGTCTGAT-3', 'CGATGGCGTGAGTATA-3', 'CGTCCATCAAGCCTAT-3', 'CGTTCTGAGCTAAACA-3', 'CTGAAGTGTGACGCCT-3', 'CTGATAGAGGGCTTGA-3', 'CTGGTCTTCTGCCCTA-3', 'GAACCTAAGGACACCA-3', 'GATGAGGTCTTATCTG-3', 'GATGCTATCACTATTC-3', 'GCCTCTAAGCTCCTCT-3', 'GCTCTGTAGACATAAC-3', 'GGACAGACAAGCGATG-3', 'GGCGTGTAGGTGCTTT-3', 'GTAGTCAGTTCCACAA-3', 'GTGGGTCTCAGTACGT-3', 'GTTACAGAGCGCTTAT-3', 'TAGACCAAGGAACTGC-3', 'TCGTAGAGTGAAGGCT-3', 'TGATTTCCACCGAAAG-3', 'TTAGGCAGTCCTCTTG-3', 'TTCTCAAGTCATACTG-3', 'TTGGAACGTTGCGTTA-3', 'TTTGTCAAGGCTAGGT-3', 'CACCACTGTAAAGGAG-3', 'CAACCTCGTTCGTTGA-3'), invert = F)

Clone18 <- SetIdent(Clone18, value = "Genotype")
markers.clone18 <- FindAllMarkers(Clone18,
                                  features = vargenes.filtered5,
                                  only.pos = FALSE, 
                                  min.pct = 0.01, 
                                  logfc.threshold = 0.25,
                                  test.use = "MAST",
                                  verbose = T,
                                  latent.vars = "nCount_RNA")


markers.clone18top20 <- markers.clone18 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone18, 
                     features = markers.clone18top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g18_recip <- WhichCells(CD8.clean3, idents = c('AAGGAGCCACTCTGTC-1', 'ATGAGGGCATTGGGCC-1', 'ATGCGATTCGGTCTAA-1', 'ATTCTACTCTTGTACT-1', 'ATTTCTGCATCACAAC-1', 'CATCAAGTCCACTCCA-1', 'GGACAGAAGCAACGGT-1', 'GGGTCTGCACAGTCGC-1', 'TCTCTAACACGAAGCA-1', 'TTGGAACAGCTAGTCT-1', 'AAGCCGCGTTTAGGAA-2', 'ACACCCTCATCGATTG-2', 'ACACTGATCGGCGGTT-2', 'CCTCTGACAATAGCAA-2', 'CGGTTAAAGTGCGTGA-2', 'GTTTCTAAGATCTGAA-2', 'TACTTACCAGCAGTTT-2', 'AAAGATGTCCTCTAGC-3', 'AGTTGGTCAAGAAGAG-3', 'CAGCAGCAGACGCACA-3', 'CTCTACGCATGCCCGA-3', 'GGTGCGTGTGTTCTTT-3', 'GTAGTCATCCTTGGTC-3', 'GTTCGGGTCCGAGCCA-3', 'TGGGAAGAGACAGACC-3'))
g18_don <- WhichCells(CD8.clean3, idents = c('AAAGCAAGTCATTAGC-3', 'ACGGGCTAGGACTGGT-3', 'ACTATCTGTTGAACTC-3', 'AGCATACCAGGATTGG-3', 'AGCGGTCAGGCAAAGA-3', 'AGGCCACAGCTGTTCA-3', 'ATCATCTGTTCAGTAC-3', 'ATCTGCCCATGGGACA-3', 'ATTATCCTCTGATTCT-3', 'ATTGGTGTCAAAGACA-3', 'ATTTCTGAGACCCACC-3', 'CAAGGCCAGTAGTGCG-3', 'CAAGGCCCAGCTGCTG-3', 'CACACCTCATGGAATA-3', 'CACTCCAAGAAGAAGC-3', 'CCCAGTTTCTAACGGT-3', 'CCTAAAGGTGTCTGAT-3', 'CGATGGCGTGAGTATA-3', 'CGTCCATCAAGCCTAT-3', 'CGTTCTGAGCTAAACA-3', 'CTGAAGTGTGACGCCT-3', 'CTGATAGAGGGCTTGA-3', 'CTGGTCTTCTGCCCTA-3', 'GAACCTAAGGACACCA-3', 'GATGAGGTCTTATCTG-3', 'GATGCTATCACTATTC-3', 'GCCTCTAAGCTCCTCT-3', 'GCTCTGTAGACATAAC-3', 'GGACAGACAAGCGATG-3', 'GGCGTGTAGGTGCTTT-3', 'GTAGTCAGTTCCACAA-3', 'GTGGGTCTCAGTACGT-3', 'GTTACAGAGCGCTTAT-3', 'TAGACCAAGGAACTGC-3', 'TCGTAGAGTGAAGGCT-3', 'TGATTTCCACCGAAAG-3', 'TTAGGCAGTCCTCTTG-3', 'TTCTCAAGTCATACTG-3', 'TTGGAACGTTGCGTTA-3', 'TTTGTCAAGGCTAGGT-3', 'CACCACTGTAAAGGAG-3', 'CAACCTCGTTCGTTGA-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g18_recip, g18_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")


#19
Recip19 =  c('CTGATCCTCACCCTCA-1', 'CCCAATCCAATCAGAA-2', 'GAGGTGACACATCCAA-2', 'GGATGTTTCGCGTAGC-2', 'GTGCATACAATGTAAG-2', 'TCGCGAGCAGGTGCCT-2', 'TCTGAGAAGGCTAGAC-2', 'TGTATTCAGTTATCGC-2', 'ATTCTACCATCCCATC-2', 'ATCGAGTGTCTAGCCG-3', 'CAGCATACACTTGGAT-3', 'GCGCAGTTCTGACCTC-3')
Don19 =  c('AACCGCGTCCTGCTTG-3', 'ACGGGCTGTCTCACCT-3', 'ATGGGAGAGGTGCAAC-3', 'ATTGGACTCTTACCGC-3', 'CAGATCAAGAAGGTTT-3', 'CAGCAGCTCTCTTATG-3', 'CAGCGACCAGGATTGG-3', 'CATCCACCAATAGAGT-3', 'CGTTGGGTCGTAGATC-3', 'CTGAAACGTCGTGGCT-3', 'CTGAAACTCTGTCTAT-3', 'CTTGGCTCAAGCCCAC-3', 'GAAATGACAAGCCTAT-3', 'GCACTCTTCGCAAGCC-3', 'GCAGTTAGTTGAGTTC-3', 'GCATGCGTCCTTTACA-3', 'GCGAGAACAGATCCAT-3', 'GTCCTCAGTCGGGTCT-3', 'GTCGGGTCAGGGTATG-3', 'GTGCATACAGCATACT-3', 'GTTACAGGTATCGCAT-3', 'TAAGCGTCAAGGTTCT-3', 'TGGACGCAGAGCAATT-3', 'TTTGGTTTCAGATAAG-3', 'CAACTAGCACATCCAA-3', 'ACTGATGTCTTAGCCC-3', 'TAGTTGGGTTAAAGAC-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone19 <- subset(CD8.clean3, ident = c('CTGATCCTCACCCTCA-1', 'CCCAATCCAATCAGAA-2', 'GAGGTGACACATCCAA-2', 'GGATGTTTCGCGTAGC-2', 'GTGCATACAATGTAAG-2', 'TCGCGAGCAGGTGCCT-2', 'TCTGAGAAGGCTAGAC-2', 'ATTCTACCATCCCATC-2', 'ATCGAGTGTCTAGCCG-3', 'CAGCATACACTTGGAT-3', 'GCGCAGTTCTGACCTC-3', 'AACCGCGTCCTGCTTG-3', 'ACGGGCTGTCTCACCT-3', 'ATGGGAGAGGTGCAAC-3', 'ATTGGACTCTTACCGC-3', 'CAGATCAAGAAGGTTT-3', 'CAGCAGCTCTCTTATG-3', 'CAGCGACCAGGATTGG-3', 'CATCCACCAATAGAGT-3', 'CGTTGGGTCGTAGATC-3', 'CTGAAACGTCGTGGCT-3', 'CTGAAACTCTGTCTAT-3', 'CTTGGCTCAAGCCCAC-3', 'GAAATGACAAGCCTAT-3', 'GCACTCTTCGCAAGCC-3', 'GCAGTTAGTTGAGTTC-3', 'GCATGCGTCCTTTACA-3', 'GCGAGAACAGATCCAT-3', 'GTCCTCAGTCGGGTCT-3', 'GTCGGGTCAGGGTATG-3', 'GTGCATACAGCATACT-3', 'GTTACAGGTATCGCAT-3', 'TAAGCGTCAAGGTTCT-3', 'TGGACGCAGAGCAATT-3', 'TTTGGTTTCAGATAAG-3', 'CAACTAGCACATCCAA-3', 'ACTGATGTCTTAGCCC-3', 'TAGTTGGGTTAAAGAC-3'), invert = F)

Clone19 <- SetIdent(Clone19, value = "Genotype")
markers.clone19 <- FindAllMarkers(Clone19,
                                  features = vargenes.filtered5,
                                  only.pos = FALSE, 
                                  min.pct = 0.01, 
                                  logfc.threshold = 0.25,
                                  test.use = "MAST",
                                  verbose = T,
                                  latent.vars = "nCount_RNA")


markers.clone19top20 <- markers.clone19 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone19, 
                     features = markers.clone19top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g19_recip <- WhichCells(CD8.clean3, idents = c('CTGATCCTCACCCTCA-1', 'CCCAATCCAATCAGAA-2', 'GAGGTGACACATCCAA-2', 'GGATGTTTCGCGTAGC-2', 'GTGCATACAATGTAAG-2', 'TCGCGAGCAGGTGCCT-2', 'TCTGAGAAGGCTAGAC-2',  'ATTCTACCATCCCATC-2', 'ATCGAGTGTCTAGCCG-3', 'CAGCATACACTTGGAT-3', 'GCGCAGTTCTGACCTC-3'))
g19_don <- WhichCells(CD8.clean3, idents = c('AACCGCGTCCTGCTTG-3', 'ACGGGCTGTCTCACCT-3', 'ATGGGAGAGGTGCAAC-3', 'ATTGGACTCTTACCGC-3', 'CAGATCAAGAAGGTTT-3', 'CAGCAGCTCTCTTATG-3', 'CAGCGACCAGGATTGG-3', 'CATCCACCAATAGAGT-3', 'CGTTGGGTCGTAGATC-3', 'CTGAAACGTCGTGGCT-3', 'CTGAAACTCTGTCTAT-3', 'CTTGGCTCAAGCCCAC-3', 'GAAATGACAAGCCTAT-3', 'GCACTCTTCGCAAGCC-3', 'GCAGTTAGTTGAGTTC-3', 'GCATGCGTCCTTTACA-3', 'GCGAGAACAGATCCAT-3', 'GTCCTCAGTCGGGTCT-3', 'GTCGGGTCAGGGTATG-3', 'GTGCATACAGCATACT-3', 'GTTACAGGTATCGCAT-3', 'TAAGCGTCAAGGTTCT-3', 'TGGACGCAGAGCAATT-3', 'TTTGGTTTCAGATAAG-3', 'CAACTAGCACATCCAA-3', 'ACTGATGTCTTAGCCC-3', 'TAGTTGGGTTAAAGAC-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g19_recip, g19_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")

#20
Recip20 =  c('ACGAGGAGTGTGAATA-1', 'AGGGATGGTTTGTTTC-1', 'GATCGATTCCACGTGG-1', 'GGGATGAGTACATCCA-1', 'TCGAGGCCAGGGCATA-1', 'AAACCTGAGGGATACC-2', 'AACTCCCAGTCGATAA-2', 'ATAGACCCAGGACCCT-2', 'CAGAGAGTCCGATATG-2', 'CGCTGGATCGTCGTTC-2', 'CTGTGCTCAGGAATGC-2', 'TGGCCAGCATCGGGTC-2', 'ACGCAGCAGGCTCTTA-3', 'CCCTCCTGTCTCCCTA-3', 'GTCTTCGCAGTCAGAG-3', 'TGCTACCGTTGCGTTA-3')
Don20 =  c('AAAGCAAGTCATTAGC-3', 'ACGGGCTAGGACTGGT-3', 'ACTATCTGTTGAACTC-3', 'AGCATACCAGGATTGG-3', 'AGCGGTCAGGCAAAGA-3', 'AGGCCACAGCTGTTCA-3', 'ATCATCTGTTCAGTAC-3', 'ATCTGCCCATGGGACA-3', 'ATTATCCTCTGATTCT-3', 'ATTGGTGTCAAAGACA-3', 'ATTTCTGAGACCCACC-3', 'CAAGGCCAGTAGTGCG-3', 'CAAGGCCCAGCTGCTG-3', 'CACACCTCATGGAATA-3', 'CACTCCAAGAAGAAGC-3', 'CCCAGTTTCTAACGGT-3', 'CCTAAAGGTGTCTGAT-3', 'CGATGGCGTGAGTATA-3', 'CGTCCATCAAGCCTAT-3', 'CGTTCTGAGCTAAACA-3', 'CTGAAGTGTGACGCCT-3', 'CTGATAGAGGGCTTGA-3', 'CTGGTCTTCTGCCCTA-3', 'GAACCTAAGGACACCA-3', 'GATGAGGTCTTATCTG-3', 'GATGCTATCACTATTC-3', 'GCCTCTAAGCTCCTCT-3', 'GCTCTGTAGACATAAC-3', 'GGACAGACAAGCGATG-3', 'GGCGTGTAGGTGCTTT-3', 'GTAGTCAGTTCCACAA-3', 'GTGGGTCTCAGTACGT-3', 'GTTACAGAGCGCTTAT-3', 'TAGACCAAGGAACTGC-3', 'TCGTAGAGTGAAGGCT-3', 'TGATTTCCACCGAAAG-3', 'TTAGGCAGTCCTCTTG-3', 'TTCTCAAGTCATACTG-3', 'TTGGAACGTTGCGTTA-3', 'TTTGTCAAGGCTAGGT-3', 'CACCACTGTAAAGGAG-3', 'CAACCTCGTTCGTTGA-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone20 <- subset(CD8.clean3, ident = c('ACGAGGAGTGTGAATA-1', 'AGGGATGGTTTGTTTC-1', 'GATCGATTCCACGTGG-1', 'GGGATGAGTACATCCA-1', 'TCGAGGCCAGGGCATA-1', 'AAACCTGAGGGATACC-2', 'AACTCCCAGTCGATAA-2', 'ATAGACCCAGGACCCT-2', 'CAGAGAGTCCGATATG-2', 'CGCTGGATCGTCGTTC-2', 'CTGTGCTCAGGAATGC-2', 'TGGCCAGCATCGGGTC-2', 'ACGCAGCAGGCTCTTA-3', 'CCCTCCTGTCTCCCTA-3', 'GTCTTCGCAGTCAGAG-3', 'TGCTACCGTTGCGTTA-3', 'AAAGCAAGTCATTAGC-3', 'ACGGGCTAGGACTGGT-3', 'ACTATCTGTTGAACTC-3', 'AGCATACCAGGATTGG-3', 'AGCGGTCAGGCAAAGA-3', 'AGGCCACAGCTGTTCA-3', 'ATCATCTGTTCAGTAC-3', 'ATCTGCCCATGGGACA-3', 'ATTATCCTCTGATTCT-3', 'ATTGGTGTCAAAGACA-3', 'ATTTCTGAGACCCACC-3', 'CAAGGCCAGTAGTGCG-3', 'CAAGGCCCAGCTGCTG-3', 'CACACCTCATGGAATA-3', 'CACTCCAAGAAGAAGC-3', 'CCCAGTTTCTAACGGT-3', 'CCTAAAGGTGTCTGAT-3', 'CGATGGCGTGAGTATA-3', 'CGTCCATCAAGCCTAT-3', 'CGTTCTGAGCTAAACA-3', 'CTGAAGTGTGACGCCT-3', 'CTGATAGAGGGCTTGA-3', 'CTGGTCTTCTGCCCTA-3', 'GAACCTAAGGACACCA-3', 'GATGAGGTCTTATCTG-3', 'GATGCTATCACTATTC-3', 'GCCTCTAAGCTCCTCT-3', 'GCTCTGTAGACATAAC-3', 'GGACAGACAAGCGATG-3', 'GGCGTGTAGGTGCTTT-3', 'GTAGTCAGTTCCACAA-3', 'GTGGGTCTCAGTACGT-3', 'GTTACAGAGCGCTTAT-3', 'TAGACCAAGGAACTGC-3', 'TCGTAGAGTGAAGGCT-3', 'TGATTTCCACCGAAAG-3', 'TTAGGCAGTCCTCTTG-3', 'TTCTCAAGTCATACTG-3', 'TTGGAACGTTGCGTTA-3', 'TTTGTCAAGGCTAGGT-3', 'CACCACTGTAAAGGAG-3', 'CAACCTCGTTCGTTGA-3'), invert = F)

Clone20 <- SetIdent(Clone20, value = "Genotype")
markers.clone20 <- FindAllMarkers(Clone20,
                                  features = vargenes.filtered5,
                                  only.pos = FALSE, 
                                  min.pct = 0.01, 
                                  logfc.threshold = 0.25,
                                  test.use = "MAST",
                                  verbose = T,
                                  latent.vars = "nCount_RNA")


markers.clone20top20 <- markers.clone20 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone20, 
                     features = markers.clone20top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g20_recip <- WhichCells(CD8.clean3, idents = c('ACGAGGAGTGTGAATA-1', 'AGGGATGGTTTGTTTC-1', 'GATCGATTCCACGTGG-1', 'GGGATGAGTACATCCA-1', 'TCGAGGCCAGGGCATA-1', 'AAACCTGAGGGATACC-2', 'AACTCCCAGTCGATAA-2', 'ATAGACCCAGGACCCT-2', 'CAGAGAGTCCGATATG-2', 'CGCTGGATCGTCGTTC-2', 'CTGTGCTCAGGAATGC-2', 'TGGCCAGCATCGGGTC-2', 'ACGCAGCAGGCTCTTA-3', 'CCCTCCTGTCTCCCTA-3', 'GTCTTCGCAGTCAGAG-3', 'TGCTACCGTTGCGTTA-3'))
g20_don <- WhichCells(CD8.clean3, idents = c('AAAGCAAGTCATTAGC-3', 'ACGGGCTAGGACTGGT-3', 'ACTATCTGTTGAACTC-3', 'AGCATACCAGGATTGG-3', 'AGCGGTCAGGCAAAGA-3', 'AGGCCACAGCTGTTCA-3', 'ATCATCTGTTCAGTAC-3', 'ATCTGCCCATGGGACA-3', 'ATTATCCTCTGATTCT-3', 'ATTGGTGTCAAAGACA-3', 'ATTTCTGAGACCCACC-3', 'CAAGGCCAGTAGTGCG-3', 'CAAGGCCCAGCTGCTG-3', 'CACACCTCATGGAATA-3', 'CACTCCAAGAAGAAGC-3', 'CCCAGTTTCTAACGGT-3', 'CCTAAAGGTGTCTGAT-3', 'CGATGGCGTGAGTATA-3', 'CGTCCATCAAGCCTAT-3', 'CGTTCTGAGCTAAACA-3', 'CTGAAGTGTGACGCCT-3', 'CTGATAGAGGGCTTGA-3', 'CTGGTCTTCTGCCCTA-3', 'GAACCTAAGGACACCA-3', 'GATGAGGTCTTATCTG-3', 'GATGCTATCACTATTC-3', 'GCCTCTAAGCTCCTCT-3', 'GCTCTGTAGACATAAC-3', 'GGACAGACAAGCGATG-3', 'GGCGTGTAGGTGCTTT-3', 'GTAGTCAGTTCCACAA-3', 'GTGGGTCTCAGTACGT-3', 'GTTACAGAGCGCTTAT-3', 'TAGACCAAGGAACTGC-3', 'TCGTAGAGTGAAGGCT-3', 'TGATTTCCACCGAAAG-3', 'TTAGGCAGTCCTCTTG-3', 'TTCTCAAGTCATACTG-3', 'TTGGAACGTTGCGTTA-3', 'TTTGTCAAGGCTAGGT-3', 'CACCACTGTAAAGGAG-3', 'CAACCTCGTTCGTTGA-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g20_recip, g20_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")

#21
Recip21 =  c('CTGATCCTCACCCTCA-1', 'GAGGTGACACATCCAA-2', 'GGATGTTTCGCGTAGC-2', 'GTGCATACAATGTAAG-2', 'TCGCGAGCAGGTGCCT-2', 'TCTGAGAAGGCTAGAC-2', 'ATTCTACCATCCCATC-2', 'ATCGAGTGTCTAGCCG-3', 'CAGCATACACTTGGAT-3', 'GCGCAGTTCTGACCTC-3')
Don21 =  c('AAAGCAAGTCATTAGC-3', 'ACGGGCTAGGACTGGT-3', 'ACTATCTGTTGAACTC-3', 'AGCATACCAGGATTGG-3', 'AGCGGTCAGGCAAAGA-3', 'AGGCCACAGCTGTTCA-3', 'ATCATCTGTTCAGTAC-3', 'ATCTGCCCATGGGACA-3', 'ATTATCCTCTGATTCT-3', 'ATTGGTGTCAAAGACA-3', 'ATTTCTGAGACCCACC-3', 'CAAGGCCAGTAGTGCG-3', 'CAAGGCCCAGCTGCTG-3', 'CACACCTCATGGAATA-3', 'CACTCCAAGAAGAAGC-3', 'CCCAGTTTCTAACGGT-3', 'CCTAAAGGTGTCTGAT-3', 'CGATGGCGTGAGTATA-3', 'CGTCCATCAAGCCTAT-3', 'CGTTCTGAGCTAAACA-3', 'CTGAAGTGTGACGCCT-3', 'CTGATAGAGGGCTTGA-3', 'CTGGTCTTCTGCCCTA-3', 'GAACCTAAGGACACCA-3', 'GATGAGGTCTTATCTG-3', 'GATGCTATCACTATTC-3', 'GCCTCTAAGCTCCTCT-3', 'GCTCTGTAGACATAAC-3', 'GGACAGACAAGCGATG-3', 'GGCGTGTAGGTGCTTT-3', 'GTAGTCAGTTCCACAA-3', 'GTGGGTCTCAGTACGT-3', 'GTTACAGAGCGCTTAT-3', 'TAGACCAAGGAACTGC-3', 'TCGTAGAGTGAAGGCT-3', 'TGATTTCCACCGAAAG-3', 'TTAGGCAGTCCTCTTG-3', 'TTCTCAAGTCATACTG-3', 'TTGGAACGTTGCGTTA-3', 'TTTGTCAAGGCTAGGT-3', 'CACCACTGTAAAGGAG-3', 'CAACCTCGTTCGTTGA-3')

CD8.clean3 <- SetIdent(CD8.clean3, value = "CellID")
Clone21 <- subset(CD8.clean3, ident = c('CTGATCCTCACCCTCA-1', 'GAGGTGACACATCCAA-2', 'GGATGTTTCGCGTAGC-2', 'GTGCATACAATGTAAG-2', 'TCGCGAGCAGGTGCCT-2', 'TCTGAGAAGGCTAGAC-2', 'ATTCTACCATCCCATC-2', 'ATCGAGTGTCTAGCCG-3', 'CAGCATACACTTGGAT-3', 'GCGCAGTTCTGACCTC-3', 'AAAGCAAGTCATTAGC-3', 'ACGGGCTAGGACTGGT-3', 'ACTATCTGTTGAACTC-3', 'AGCATACCAGGATTGG-3', 'AGCGGTCAGGCAAAGA-3', 'AGGCCACAGCTGTTCA-3', 'ATCATCTGTTCAGTAC-3', 'ATCTGCCCATGGGACA-3', 'ATTATCCTCTGATTCT-3', 'ATTGGTGTCAAAGACA-3', 'ATTTCTGAGACCCACC-3', 'CAAGGCCAGTAGTGCG-3', 'CAAGGCCCAGCTGCTG-3', 'CACACCTCATGGAATA-3', 'CACTCCAAGAAGAAGC-3', 'CCCAGTTTCTAACGGT-3', 'CCTAAAGGTGTCTGAT-3', 'CGATGGCGTGAGTATA-3', 'CGTCCATCAAGCCTAT-3', 'CGTTCTGAGCTAAACA-3', 'CTGAAGTGTGACGCCT-3', 'CTGATAGAGGGCTTGA-3', 'CTGGTCTTCTGCCCTA-3', 'GAACCTAAGGACACCA-3', 'GATGAGGTCTTATCTG-3', 'GATGCTATCACTATTC-3', 'GCCTCTAAGCTCCTCT-3', 'GCTCTGTAGACATAAC-3', 'GGACAGACAAGCGATG-3', 'GGCGTGTAGGTGCTTT-3', 'GTAGTCAGTTCCACAA-3', 'GTGGGTCTCAGTACGT-3', 'GTTACAGAGCGCTTAT-3', 'TAGACCAAGGAACTGC-3', 'TCGTAGAGTGAAGGCT-3', 'TGATTTCCACCGAAAG-3', 'TTAGGCAGTCCTCTTG-3', 'TTCTCAAGTCATACTG-3', 'TTGGAACGTTGCGTTA-3', 'TTTGTCAAGGCTAGGT-3', 'CACCACTGTAAAGGAG-3', 'CAACCTCGTTCGTTGA-3'), invert = F)

Clone21 <- SetIdent(Clone21, value = "Genotype")
markers.clone21 <- FindAllMarkers(Clone21,
                                  features = vargenes.filtered5,
                                  only.pos = FALSE, 
                                  min.pct = 0.01, 
                                  logfc.threshold = 0.25,
                                  test.use = "MAST",
                                  verbose = T,
                                  latent.vars = "nCount_RNA")


markers.clone21top20 <- markers.clone21 %>% 
  group_by("Genotype") %>% 
  top_n(n = 20, wt = avg_log2FC)


heatmap <- DoHeatmap(Clone21, 
                     features = markers.clone21top20$gene, 
                     group.bar = T,
                     draw.lines = F,
                     disp.min = -2.5, disp.max = 2.5,
                     raster = F)

heatmap + scale_fill_gradientn(colors = c("#bdbbff", "#FFFFFF", "#c20000"))

g21_recip <- WhichCells(CD8.clean3, idents = c('CTGATCCTCACCCTCA-1', 'GAGGTGACACATCCAA-2', 'GGATGTTTCGCGTAGC-2', 'GTGCATACAATGTAAG-2', 'TCGCGAGCAGGTGCCT-2', 'TCTGAGAAGGCTAGAC-2', 'ATTCTACCATCCCATC-2', 'ATCGAGTGTCTAGCCG-3', 'CAGCATACACTTGGAT-3', 'GCGCAGTTCTGACCTC-3'))
g21_don <- WhichCells(CD8.clean3, idents = c('AAAGCAAGTCATTAGC-3', 'ACGGGCTAGGACTGGT-3', 'ACTATCTGTTGAACTC-3', 'AGCATACCAGGATTGG-3', 'AGCGGTCAGGCAAAGA-3', 'AGGCCACAGCTGTTCA-3', 'ATCATCTGTTCAGTAC-3', 'ATCTGCCCATGGGACA-3', 'ATTATCCTCTGATTCT-3', 'ATTGGTGTCAAAGACA-3', 'ATTTCTGAGACCCACC-3', 'CAAGGCCAGTAGTGCG-3', 'CAAGGCCCAGCTGCTG-3', 'CACACCTCATGGAATA-3', 'CACTCCAAGAAGAAGC-3', 'CCCAGTTTCTAACGGT-3', 'CCTAAAGGTGTCTGAT-3', 'CGATGGCGTGAGTATA-3', 'CGTCCATCAAGCCTAT-3', 'CGTTCTGAGCTAAACA-3', 'CTGAAGTGTGACGCCT-3', 'CTGATAGAGGGCTTGA-3', 'CTGGTCTTCTGCCCTA-3', 'GAACCTAAGGACACCA-3', 'GATGAGGTCTTATCTG-3', 'GATGCTATCACTATTC-3', 'GCCTCTAAGCTCCTCT-3', 'GCTCTGTAGACATAAC-3', 'GGACAGACAAGCGATG-3', 'GGCGTGTAGGTGCTTT-3', 'GTAGTCAGTTCCACAA-3', 'GTGGGTCTCAGTACGT-3', 'GTTACAGAGCGCTTAT-3', 'TAGACCAAGGAACTGC-3', 'TCGTAGAGTGAAGGCT-3', 'TGATTTCCACCGAAAG-3', 'TTAGGCAGTCCTCTTG-3', 'TTCTCAAGTCATACTG-3', 'TTGGAACGTTGCGTTA-3', 'TTTGTCAAGGCTAGGT-3', 'CACCACTGTAAAGGAG-3', 'CAACCTCGTTCGTTGA-3'))
DimPlot(CD8.clean3, label=F, group.by="Genotype", cells.highlight= list(g21_recip, g21_don), cols.highlight = c("darkred", "darkblue"), cols= "grey")



#Filter expanding clones

CD8.clean3 <- SetIdent(CD8.clean3, value = "ClonotypeID")
CD8Expanding60 <- subset(CD8.clean3, ident = c("60_clonotype10",
                                               '60_clonotype11',
                                               '60_clonotype12',
                                               '60_clonotype13',
                                               '60_clonotype17',
                                               '60_clonotype19',
                                               '60_clonotype2',
                                               '60_clonotype20',
                                               '60_clonotype26',
                                               '60_clonotype3',
                                               '60_clonotype31',
                                               '60_clonotype34',
                                               '60_clonotype35',
                                               '60_clonotype38',
                                               '60_clonotype41',
                                               '60_clonotype43',
                                               '60_clonotype5',
                                               '60_clonotype6',
                                               '60_clonotype60',
                                               '60_clonotype7',
                                               '60_clonotype8',
                                               '90_clonotype1',
                                               '90_clonotype19',
                                               '90_clonotype24',
                                               '90_clonotype3',
                                               '90_clonotype53',
                                               '90_clonotype6',
                                               '90_clonotype7',
                                               '60_clonotype104',
                                               '60_clonotype16',
                                               '60_clonotype23',
                                               '60_clonotype24',
                                               '60_clonotype29',
                                               '60_clonotype33',
                                               '60_clonotype36',
                                               '60_clonotype37',
                                               '60_clonotype39',
                                               '60_clonotype50',
                                               '60_clonotype53',
                                               '60_clonotype54',
                                               '60_clonotype55',
                                               '60_clonotype57',
                                               '60_clonotype58',
                                               '60_clonotype61',
                                               '60_clonotype66',
                                               '60_clonotype75',
                                               '60_clonotype82',
                                               '60_clonotype85',
                                               '60_clonotype93',
                                               "60_clonotype94"))

counts <- with(CD8Expanding60@meta.data, table(Time,ClonotypeID,Genotype,sample))

ggplot(as.data.frame(counts), aes(x = Time, y = Freq, group = ClonotypeID)) + 
  geom_line(aes(group. = Genotype)) +
  scale_y_log10()


































#### make .csv with barcode, umap coordinates, and monaco lables

barcodeumap <- as.data.frame(sce@int_colData@listData$reducedDims@listData$UMAP)
head(barcodeumap)
tail(barcodeumap)

setDT(barcodeumap, keep.rownames = T)[]
colnames(barcodeumap) <- c('barcode', 'UMAP_1', 'UMAP_2')

head(barcodeumap)
tail(barcodeumap)


table(CD8.clean3@meta.data$monaco.labels.fine)


metadata <- as.data.frame(CD8.clean3@meta.data)
head(metadata)
tail(metadata)

#merge sample and cellIDs

metadata$barcode <- paste(metadata$sample, metadata$CellID, sep = "_")
head(metadata)

setDT(metadata, keep.rownames = T)[]
colnames(metadata) <- c('Number', 'orig.ident', 'nCount_RNA', "nFeature_RNA", "CellID", "Genotype", "Annotations", "ClonotypeID", "ClustermapRowNumber", "SampleID", "Tissue", "Time", "sample", "percent.mt", "prcent.ribo", "RNA_snn_res.0.8", "seurat_clusters", "monaco.labels", "monaco.labels.fine", "barcode")

head(metadata)
tail(metadata)

UMAPMeta <- merge(metadata, barcodeumap, by.x = 'barcode', by.y = 'barcode', all.x = T, all.y = F, sort = F)
head(UMAPMeta)

write.csv(UMAPMeta, "CMV_CD8_UMAP_MonacoLabels.csv")

#Make figure 2C D with filtered cells

CD8.clean3@meta.data$XIST <- CD8.clean3@assays$RNA@data["XIST",]
CD8.clean3@meta.data$RPS4Y1 <- CD8.clean3@assays$RNA@data["RPS4Y1",]

write.csv(CD8.clean3@meta.data, "Meta.csv")


write.csv(t(Patient2@assays$RNA@data[c("RPS4Y1", "XIST"),]), "Patient2.csv")
write.csv(t(Patient1@assays$RNA@data[c("RPS4Y1", "XIST"),]), "Patient1.csv")
write.csv(t(Patient3@assays$RNA@data[c("RPS4Y1", "XIST"),]), "Patient3.csv")
write.csv(t(Patient4@assays$RNA@data[c("RPS4Y1", "XIST"),]), "Patient4.csv")

