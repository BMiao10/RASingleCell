#install.packages("pheatmap")
#BiocManager::install("SingleCellExperiment")
#library("scater")
#library("scRNAseq")
#library("dittoSeq")
library("Seurat")
library("SingleCellExperiment")
library("dplyr")
#library("ggsci")
library("Matrix")

# --Load data  ---------------------------------------------------------------------------------------------
setwd("~/PBMC")

### params
inputFolder = './dataOutput/toR/dim-reduc'
inputPCs = "./dataOutput/harmony_PCs.csv"
outputPath = './dataOutput/batch1/'

# load data
inputMatrix = paste(inputFolder, ".mtx", sep="")
inputGenes = paste(inputFolder, "_genes.csv", sep="")
inputMeta = paste(inputFolder, "_metadata.csv", sep="")

adataX = readMM(inputMatrix)
adataX = as.sparse(adataX)
genes = read.csv(file = inputGenes)
metadata = read.csv(file = inputMeta)
PCs = read.csv(file=inputPCs)

# create seurat object
dimnames(adataX) = list(genes$genes,metadata$barcodes)
UMI_RA <- CreateSeuratObject(adataX, assay = "RNA")

# all character columns to factor:
#metadata[] <- lapply(metadata, as.factor)

# --PCA projection  ----------------------------------------------------------------------------------------------------

#all.gene list
all.genes = rownames(UMI_RA)

#scaling on all.genes 
#PCA
UMI_RA = ScaleData(UMI_RA, features = all.genes)
UMI_RA = RunPCA(UMI_RA, features = VariableFeatures(UMI_RA)) 
variable_feature = print(UMI_RA[["pca"]], dims = 1:5, nfeatures = 5)

#variable feature exportation
cat("Variable feature\n", file = "variables_feature_PCA.txt" , append = TRUE )
capture.output(variable_feature, file = "variables_feature_PCA.txt", append = TRUE)

#PCA plot in dimension reduction and heatmap
VizDimLoadings(UMI_RA, dims = 1:2, reduction = "pca")
DimPlot(UMI_RA, reduction="pca", group.by = 'sample', label = FALSE)
DimHeatmap(UMI_RA, dims = 1:5, cells = 500, balanced=TRUE)
FeaturePlot(UMI_RA, features=c("nCount_RNA"))











