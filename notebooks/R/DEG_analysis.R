
BiocManager::install("cole-trapnell-lab/monocle3")
library("SingleCellExperiment")
library("Matrix")
library("Matrix.utils")
library("dplyr")
library("tidyverse")
library("Seurat")
library("scater")
library("magrittr")
library("DESeq2")
library("RColorBrewer")
library("ashr")
library("tibble")
#library("muscat")
library("edgeR")
library("scran")
library("presto")
library("monocle3")


# libraries used chosen based on:
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2599-6

#####################################
############## EdgeR ################
#####################################

# read in aggregated data all batches
inputFolder = './dataOutput/all-batches/toRfiles/pseudobulk'

inputMatrix = paste(inputFolder, ".mtx", sep="")
inputGenes = paste(inputFolder, "_genes.csv", sep="")
inputMeta = paste(inputFolder, "_metadata.csv", sep="")

# load datasets
counts = readMM(inputMatrix)
counts = as.sparse(counts)
genes = read.csv(file = inputGenes)
metadata = read.csv(file = inputMeta)

# label dimensions
dimnames(counts) = list(metadata$label, genes$genes)

# create design
metadata$"sampleType" = as.factor(metadata$"sampleType")
metadata$"sex" = as.factor(metadata$"sex")
metadata$"Lane" = as.factor(metadata$"Lane")

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(counts), pattern = "_"),  `[`, 1)

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
counts_df <- split.data.frame(counts,  factor(splitf)) 
counts_df = lapply(counts_df, function(u) t(u))

for (celltype in counts_df) {
  curr_dimnames = celltype@Dimnames[[2]]
  curr_type = str_split(curr_dimnames[1], "_")[[1]][1]
  
  new_dimnames = lapply(strsplit(curr_dimnames, "_"), function(u) u[[2]])
  new_dimnames = unlist(new_dimnames)
  
  dim1names = unlist(counts_df[[curr_type]]@Dimnames[1])
  
  dimnames(counts_df[[curr_type]]) = list(dim1names,new_dimnames)
}

for (k in 1:7) {
  # get current cell type
  clusters <- unique(metadata$singleR_blueprint)
  curr_celltype = clusters[k]
  print(curr_celltype)
  
  # Subset the metadata to only the B cells
  cluster_metadata <- metadata[which(metadata$singleR_blueprint == curr_celltype), ]
  
  # Assign the rownames of the metadata to be the sample IDs
  rownames(cluster_metadata) <- cluster_metadata$Sample
  
  # Subset the counts to only the B cells
  counts_sub <- counts_df[[curr_celltype]]
  
  # add cell detection rate as factor
  #cdr <- scale(colMeans(counts_sub) > 0)
  
  ### Make single cell experiment
  sceAdata <- SingleCellExperiment(assays = counts_sub, colData = cluster_metadata)
  
  ## convert to edgR object
  dge <- convertTo(sceAdata, type="edgeR", assay.type = 1, )
  
  # correct cell count
  dge$samples$lib.size = dge$samples$cellcount
  
  # convert design variables to factor
  dge$samples$sampleType = as.factor(dge$samples$sampleType)
  dge$samples$batch = as.factor(dge$samples$batch)
  
  # initialize norm factors and design
  dge <- calcNormFactors(dge)
  
  # construct design & contrast matrix
  design <- model.matrix(~ 0 + sampleType + sex + batch, data=dge$samples)
  
  # make contrasts
  (contrast <- makeContrasts("sampleTypeRA-sampleTypeControl", levels = design))
  
  # estimate dispersion
  dge <- estimateDisp(dge, design = design)
  
  # fit model with glmQLFit (provides more 
  # accurate type I error rate control than glmFit)
  fit <- glmQLFit(dge, design = design)
  qlf_as.is <- glmQLFTest(fit, contrast=contrast)
  
  # top genes
  qlf_as.adj = topTags(qlf_as.is, n=nrow(counts_sub))
  
  # Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- qlf_as.adj %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  
  # create outputPath
  outputPath = paste("./dataOutput/all-batches/pvalsOut/", curr_celltype, sep="")
  
  # Write all results to file
  write.csv(res_tbl,
            paste(outputPath,"_RAvControl_edgeR.csv", sep=""),
            quote = FALSE, 
            row.names = FALSE)
}




#####################################
############## PRESTO ###############
#####################################
# all batches
inputFolder = './dataOutput/all-batches/toRfiles/all-batches_subset'
outputPath = './dataOutput/all-batches/'

inputMatrix = paste(inputFolder, ".mtx", sep="")
inputGenes = paste(inputFolder, "_genes.csv", sep="")
inputMeta = paste(inputFolder, "_metadata.csv", sep="")

adataX = readMM(inputMatrix)
adataX = as.sparse(adataX)
genes = read.csv(file = inputGenes)
metadata = read.csv(file = inputMeta)

dimnames(adataX) = list(genes$genes,metadata$X)
curr_celltype = unique(metadata$singleR_blueprint)

res <- wilcoxauc(adataX, metadata$sampleType)

# top genes
top_markers(res, n = nrow(adataX))

# create outputPath
outputPath = paste("./dataOutput/all-batches/pvalsOut/", curr_celltype, sep="")

# Write all results to file
write.csv(res,
          paste(outputPath,"_RAvControl_presto.csv", sep=""),
          quote = FALSE, 
          row.names = FALSE)



#####################################
############# Monocle ##############
#####################################


#http://cole-trapnell-lab.github.io/monocle-release/monocle3/#step-3-partition-the-cells-into-supergroups
# all batches
inputFolder = './dataOutput/all-batches/toRfiles/all-batches_normalized_genes.csv'
outputPath = './dataOutput/all-batches/'

inputMatrix = paste(inputFolder, ".mtx", sep="")
inputGenes = paste(inputFolder, "_genes.csv", sep="")
inputMeta = paste(inputFolder, "_metadata.csv", sep="")

expression_matrix = readMM(inputMatrix)
expression_matrix = as.sparse(adataX)
gene_annotation = read.csv(file = inputGenes)
cell_metadata = read.csv(file = inputMeta)

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cell.type")

cds <- learn_graph(cds)

plot_cells(cds,
           color_cells_by = "cell.type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

  


#####################################
############# DEsingle ##############
#####################################

#https://miaozhun.github.io/DEsingle/
BiocManager::install("DEsingle")


#####################################
############## muscat ############### ARCHIVE
#####################################


# all batches
inputFolder = './dataOutput/all-batches/toRfiles/averageCounts'

inputMatrix = paste(inputFolder, ".mtx", sep="")
inputGenes = paste(inputFolder, "_genes.csv", sep="")
inputMeta = paste(inputFolder, "_metadata.csv", sep="")

# load datasets
counts = readMM(inputMatrix)
counts = as.sparse(counts)
genes = read.csv(file = inputGenes)
metadata = read.csv(file = inputMeta)

# label dimensions
dimnames(counts) = list(metadata$label, genes$genes)

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(counts), pattern = "_"),  `[`, 1)

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
counts_df <- split.data.frame(counts,  factor(splitf)) 
counts_df = lapply(counts_df, function(u) t(u))

for (celltype in counts_df) {
  curr_dimnames = celltype@Dimnames[[2]]
  curr_type = str_split(curr_dimnames[1], "_")[[1]][1]
  
  new_dimnames = lapply(strsplit(curr_dimnames, "_"), function(u) u[[2]])
  new_dimnames = unlist(new_dimnames)
  
  dim1names = unlist(counts_df[[curr_type]]@Dimnames[1])
  
  dimnames(counts_df[[curr_type]]) = list(dim1names,new_dimnames)
}

# pseudobulk calculations per cell type  

# get current cell type
clusters <- unique(metadata$singleR_blueprint)
curr_celltype = clusters[7]

# Subset the metadata to only the B cells
cluster_metadata <- metadata[which(metadata$singleR_blueprint == curr_celltype), ]

# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$Sample

# Subset the counts to only the B cells
counts_sub <- counts_df[[curr_celltype]]

# make new data frame
cluster_counts <- data.frame(counts_sub[, which(colnames(counts_sub) %in% rownames(cluster_metadata))], check.names=FALSE)

adataSeurat <- CreateSeuratObject(cluster_counts, assay = "RNA")

# all character columns to factor:
#metadata[] <- lapply(metadata, as.factor)

# add leiden column to seurat object
adataSeurat@meta.data = cluster_metadata

sceAdata = as.SingleCellExperiment(x = adataSeurat)

#contrast <- c("sampleType", levels(cluster_metadata$sampleType)[2], levels(cluster_metadata$sampleType)[1])
#model.matrix(formula(~ group))

# add attributes to sce
adding_meta = list()
adding_meta$"experiment_info" = cluster_metadata
adding_meta$"agg_pars" = select(cluster_metadata, "assay", "by", "fun", "scaled")
adding_meta$"n_cells" = select(cluster_metadata, "cellcount")
sceAdata@metadata = adding_meta

# run pseudobulk (aggregation-based) DS analysis
ds_pb <- pbDS(sceAdata, method = "DESeq2", design = ~ batch + sampleType)



# Write all results to file
write.csv(res_tbl,outputPath, quote = FALSE, row.names = FALSE)



#####################################
####### DESeq2 - pseudobulk #########
#####################################

#https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html

# read in aggregated data all batches
inputFolder = './dataOutput/all-batches/toRfiles/pseudobulk'

inputMatrix = paste(inputFolder, ".mtx", sep="")
inputGenes = paste(inputFolder, "_genes.csv", sep="")
inputMeta = paste(inputFolder, "_metadata.csv", sep="")

# load datasets
counts = readMM(inputMatrix)
counts = as.sparse(counts)
genes = read.csv(file = inputGenes)
metadata = read.csv(file = inputMeta)

# label dimensions
dimnames(counts) = list(metadata$label, genes$genes)

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(counts), pattern = "_"),  `[`, 1)

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
counts_df <- split.data.frame(counts,  factor(splitf)) 
counts_df = lapply(counts_df, function(u) t(u))

for (celltype in counts_df) {
  curr_dimnames = celltype@Dimnames[[2]]
  curr_type = str_split(curr_dimnames[1], "_")[[1]][1]
  
  new_dimnames = lapply(strsplit(curr_dimnames, "_"), function(u) u[[2]])
  new_dimnames = unlist(new_dimnames)
  
  dim1names = unlist(counts_df[[curr_type]]@Dimnames[1])
  
  dimnames(counts_df[[curr_type]]) = list(dim1names,new_dimnames)
}

# pseudobulk calculations per cell type  
for (i in 1:7){
  
  # get current cell type
  clusters <- unique(metadata$singleR_blueprint)
  curr_celltype = clusters[i]
  
  # Subset the metadata to only the B cells
  cluster_metadata <- metadata[which(metadata$singleR_blueprint == curr_celltype), ]
  
  # Assign the rownames of the metadata to be the sample IDs
  rownames(cluster_metadata) <- cluster_metadata$Sample
  
  # Subset the counts to only the B cells
  counts_sub <- counts_df[[curr_celltype]]
  
  # make new data frame
  cluster_counts <- data.frame(counts_sub[, which(colnames(counts_sub) %in% rownames(cluster_metadata))], check.names=FALSE)
  
  # Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
  all(rownames(cluster_metadata) == colnames(cluster_counts)) 
  
  cluster_metadata['batch'] = as.factor(cluster_metadata$'batch')
  cluster_metadata['sampleType'] = as.factor(cluster_metadata$'sampleType')
  cluster_metadata['Lane'] = as.factor(cluster_metadata$"Lane")
  cluster_metadata['sex'] = as.factor(cluster_metadata$"sex")
  
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ sampleType + sex + Lane)
  
  # Transform counts for data visualization
  #rld <- vst(dds, blind=TRUE, fitType="local")
  rld <- varianceStabilizingTransformation(dds, blind=TRUE, fitType="parametric")
  
  # Plot PCA
  DESeq2::plotPCA(rld, intgroup = "sampleType")
  
  # Run DESeq2 differential expression analysis
  dds <- DESeq(dds)
  
  # Plot dispersion estimates
  plotDispEsts(dds)
  
  contrast <- c("sampleType", levels(cluster_metadata$sampleType)[2], levels(cluster_metadata$sampleType)[1])
  
  # resultsNames(dds)
  res <- results(dds, 
                 contrast = contrast,
                 alpha = 0.05)
  
  '
  # instead of shrinking
  
  res <- lfcShrink(dds, 
                   contrast =  contrast,
                   res=res, type="ashr")
  '
  
  # just find significant genes
  
  # Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  
  # Check results output
  arrange(res_tbl, res_tbl$padj)
  
  # create outputPath
  outputPath = paste('./dataOutput/all-batches/pvalsOut/', curr_celltype, sep="")
  
  # Write all results to file
  write.csv(res_tbl,
            paste(outputPath,"_RAvControl_DEseq2.csv", sep=""),
            quote = FALSE, 
            row.names = FALSE)
}






