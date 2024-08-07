
#####################################
############# EdgeR ###############
#####################################
#BiocManager::install("scater")

# libraries
library(scran)
library(tidyverse)
library(Matrix.utils)
library(edgeR)
library(magrittr)
library(purrr)
library(reshape2)
library(S4Vectors)
#library(tibble)
library(apeglm)
library(DESeq2)
library("pheatmap")

# set working directory
setwd("~/PBMC/dataOutput/DEGs/")

# params
inputpath = 'DEG_subsets_raw.csv'
inputMeta = 'DEG_metadata.csv'

# read in pseudobulk matrix
counts_df = read.csv(inputpath, row.names=1)
counts_df = as.matrix(counts_df)

# read in metadata for linear model design
metadata = read.csv(file = inputMeta)
metadata$celltype_id <- as.factor(metadata$celltype_id)
metadata$sample_id <- as.factor(metadata$sample_id)
metadata$cohort <- as.factor(metadata$cohort)

# Create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(counts_df), 
                                    pattern = "_",  
                                    n = 2), `[`, 1)

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
counts_df <- split.data.frame(counts_df, 
                              factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

###
n = length(counts_df)
full_tbl = data.frame()

# Run for loop to detect DEGs in each cell type
for (k in 1:n) {
  # get current cell type
  clusters <- levels(metadata$celltype_id)
  curr_celltype = clusters[k]
  print(curr_celltype)
  
  # Subset the metadata to only the B cells
  cluster_metadata <- metadata[which(metadata$celltype_id == curr_celltype), ]
  
  # Assign the rownames of the metadata to be the sample IDs
  rownames(cluster_metadata) <- cluster_metadata$sample_id
  
  # Subset the counts to only the B cells
  cluster_counts <- counts_df[[curr_celltype]]
  
  # Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
  all(rownames(cluster_metadata) == colnames(cluster_counts))        
  
  # add cell detection rate as factor
  #cdr <- scale(colMeans(cluster_counts) > 0)
  
  # create edgeR object
  dge <- DGEList(cluster_counts, samples=cluster_metadata)
  
  # correct cell count by normalization factor category
  #dge$samples$lib.size = 
  
  # convert design variables to factor
  dge$samples$batch = as.factor(dge$samples$batch)
  dge$samples$batch = relevel(dge$samples$batch, ref="3")
  
  # construct design & contrast matrix
  design <- model.matrix(~dx+batch, data=dge$samples)
  
  # initialize norm factors and design
  dge <- calcNormFactors(dge, )

  # make contrasts
  #(contrast <- makeContrasts(RA_vs_HC = dxRA-dxcontrol, levels=design))
  
  # estimate dispersion
  dge <- estimateDisp(dge, design = design)
  
  # fit model with glmQLFit provides more 
  # accurate type I error rate control than glmFit
  fit <- glmQLFit(dge, design = design, robust=TRUE)
  #fit = glmFit(dge, design, robust = TRUE)
  res <- glmLRT(fit, coef="dxRA")
  #res <- glmQLFTest(fit, contrast = contrast)
  #res <- glmLRT(fit, contrast = contrast)
  
  topTags(res)
  
  # top genes
  qlf_as.adj = topTags(res, n=nrow(cluster_counts))
  
  # Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- qlf_as.adj %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  
  # Check results output
  arrange(res_tbl, res_tbl$FDR)
  
  # concatenate all results
  res_tbl$subset = curr_celltype
  
  full_tbl <- rbind(full_tbl, res_tbl)
  
}

# Write all results to file
write.csv(full_tbl,"RAvControl_edgeR.csv",
          quote = FALSE, 
          row.names = FALSE)


