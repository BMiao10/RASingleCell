#install.packages("pheatmap")
#BiocManager::install("scRNAseq")
#library("scater")
#library("scRNAseq")
#library("dittoSeq")
library("Seurat")
library("SingleCellExperiment")
library("dplyr")
#library("ggsci")
library("Matrix")
#BiocManager::install("SingleR", build_vignettes = TRUE)
library("BiocParallel")
library("SingleR")
#library("hdf5r")
#library("ensembldb")

#####################################
############# SingleR ###############
#####################################

#######################
### Coarse annotation
### on full dataset
#######################

### params
# set working directory
setwd(dir = '~/PBMC/batch3/')

# set folder paths
inputFolder = './dataOutput/singleR/all-batches_normalized'
outputPath = './dataOutput/singleR/'

# get reference datasets
refs = list(BlueprintEncodeData(), 
            #celldex::DatabaseImmuneCellExpressionData(), 
            celldex::MonacoImmuneData())

datasets = list('blueprint','monaco') #'immunecell',
leidens = list('leiden_r3.0') #'leiden_r0.5', 'leiden_r2.5', 
textLabels = list('label.main', 'label.fine')

### Load data
adataX = readMM(paste0(inputFolder, ".mtx"))
adataX = as.sparse(adataX)
genes = read.csv(file = paste0(inputFolder, "_genes.csv"))
metadata = read.csv(file = paste0(inputFolder, "_metadata.csv"))

dimnames(adataX) = list(genes$genes,metadata$barcodes)
adataSeurat <- CreateSeuratObject(adataX, assay = "RNA")

# all character columns to factor:
metadata[] <- lapply(metadata, as.factor)

# add leiden column to seurat object
adataSeurat <- AddMetaData(adataSeurat, metadata$leiden_r0.5, col.name = 'leiden_r0.5')
adataSeurat <- AddMetaData(adataSeurat, metadata$leiden_r3.0, col.name = 'leiden_r3.0')

# convert to SingleCellExperiment for SingleR
sceAdata= as.SingleCellExperiment(adataSeurat)

### run SingleR for all datasets
for (l in leidens) {
  print(l)
  for (idx in 1:length(refs)) {
    for (t in textLabels) {
      r=refs[[idx]]
      d=datasets[idx]
      
      # singleR

      pred.cluster <- SingleR(test = sceAdata, 
                              ref = r, 
                              labels = r[[t]],
                              #labels = r[['label.fine']], 
                              #method = "cluster",
                              de.method = "classic",
                              fine.tune=TRUE,
                              prune=TRUE,
                              clusters = unlist(adataSeurat[[l]]),
                              BPPARAM = SerialParam())#MulticoreParam(8))
      plotScoreHeatmap(pred.cluster)
      table(pred.cluster$labels)
      
      adataSeurat[[paste0("singleR_",l,"_",d,"_",t)]]<- 
        pred.cluster$labels[match(adataSeurat[[]][[l]], rownames(pred.cluster))]
    }
  }
}


labeledMeta <- adataSeurat@meta.data
write.csv(labeledMeta, paste0(outputPath, 'singleR_labels.csv'))

##### ARCHIVE #####

# ref <- MonacoImmuneData()
# dataset = 'monaco'

# 
# #######################
# ### Fine annotation
# ### on subset
# #######################
# 
### Load values

ref <- celldex::MonacoImmuneData()
dataset = 'monaco'

# set working directory
setwd("~/Desktop/ButteLab/RASingleCell/PBMC/")

cells = cbind('CD4Tcells','Monocytes','NKcells','Bcells','CD8Tcells')
leiden = cbind('0.6', '0.6',  '0.6', '0.4', '0.6')

for (n in (1:length(cells))) {

  # all batches
  celltype =  cells[n]
  inputFolder = paste('./dataOutput/singleR/', celltype, "_normalized", sep="")
  outputPath = paste('./dataOutput/singleR/', celltype, sep="")
  
  inputMatrix = paste(inputFolder, ".mtx", sep="")
  inputGenes = paste(inputFolder, "_genes.csv", sep="")
  inputMeta = paste(inputFolder, "_metadata.csv", sep="")
  
  adataX = readMM(inputMatrix)
  adataX = as.sparse(adataX)
  genes = read.csv(file = inputGenes)
  metadata = read.csv(file = inputMeta)
  
  dimnames(adataX) = list(genes$genes,metadata$barcodes)
  adataSeurat <- CreateSeuratObject(adataX, assay = "RNA")
  
  # all character columns to factor:
  metadata[] <- lapply(metadata, as.factor)
  
  # add leiden column to seurat object
  currleiden = paste(celltype, "_r", leiden[n], sep="")
  adataSeurat <- AddMetaData(adataSeurat, metadata[, currleiden], col.name = 'leiden')
  
  # cluster
  sceAdata= as.SingleCellExperiment(adataSeurat)
  
  pred.cluster <- SingleR(test = sceAdata,
                          ref = ref,
                          #labels = ref$label.main,
                          labels = ref$label.fine,
                          #de.n = 50,
                          #labels = ref$label.ont,
                          #ref = list(ref1=ref1, ref2=ref2),
                          #labels = list(ref1$label.ont, ref2$label.ont),
                          #method="cluster",
                          clusters = adataSeurat$leiden,
                          BPPARAM = SerialParam())#MulticoreParam(16))
  
  plotScoreHeatmap(pred.cluster)
  table(pred.cluster$labels)
  
  adataSeurat[['SingleR_cluster_labels']]<-
    pred.cluster$labels[match(adataSeurat[[]][["leiden"]], rownames(pred.cluster))]
  
  labeledMeta <- adataSeurat@meta.data
  write.csv(labeledMeta, paste(outputPath, '_singleR_', dataset, '.csv', sep=""))
}

#####################################
######## Batch correction ###########
#####################################
# 
# BiocManager::install("immunogenomics/harmony")
# library("harmony")
# 
# # all-batches analysis
# inputFolder = './dataOutput/all-batches/all-batches_toR'
# outputPath = './dataOutput/all-batches/'
# 
# inputMatrix = paste(inputFolder, ".mtx", sep="")
# inputGenes = paste(inputFolder, "_genes.csv", sep="")
# inputMeta = paste(inputFolder, "_metadata.csv", sep="")
# 
# adataX = readMM(inputMatrix)
# genes = read.csv(file = inputGenes)
# metadata = read.csv(file = inputMeta)
# 
# dimnames(adataX) = list(genes$genes,metadata$barcodes)
# adataSeurat <- CreateSeuratObject(adataX, assay = "RNA")
# 
# # all character columns to factor:
# metadata[] <- lapply(metadata, as.factor)
# 
# # add mito percent column to seurat object
# adataSeurat <- AddMetaData(adataSeurat, metadata$batch, col.name = 'batch')
# adataSeurat <- AddMetaData(adataSeurat, metadata$Lane, col.name = 'Lane')
# 
# seurat.filtered <-  FindVariableFeatures(adataSeurat) %>% ScaleData(vars.to.regress = c("nCount_RNA")) %>% RunPCA(verbose = TRUE)
# 
# seurat.filtered[["percent.ribo"]] <- PercentageFeatureSet(seurat.filtered,
#                                                           pattern = "^RP[SL]\\d+")
# 
# seurat.filtered <- AddMetaData(seurat.filtered, metadata$Lane, col.name = 'Lane')
# 
# ###### Adding the metadata
# 
# # Run Harmony
# seurat.filtered <- RunHarmony(seurat.filtered, 
#                               group.by.vars = "Lane",
#                               # group.by.vars = "batch",
#                               kmeans_init_nstart = 20, 
#                               kmeans_init_iter_max = 100, 
#                               # epsilon.harmony = -Inf,
#                               max.iter.harmony = 20,
#                               plot_convergence = TRUE)
# 
# 
# harmony_embeds <- seurat.filtered@reductions$harmony@cell.embeddings
# write.csv(harmony_embeds, paste(outputPath, 'singleRHarmony', '.csv', sep=""))


#' @title edgeR data preparation
#'
#' @description Performs a data normalization using edgeR on the dataset. Then generates a quant matrix, and computes Pvalues and Fold Changes for each transcript
#' 
#' @param dge_list  a EdgER object containing the DGE list
#' @param samples_info a character vector containing the list of all sample names and metadata
#' @param DE_method a string giving the type of analysis :- "ET" Fisher Exact Test
#'                                                        - "QLF" Quasi likelihood ratio F test
#'                                                        - "LR" Likelihood Ratio
#'' @param padj_method a string giving the method of padj :- "holm"
#'                                                         - "hochberg"
#'                                                         - "bonferroni"
#'                                                         - "BH"
#'                                                         - "BY"
#'                                                         - "fdr"
#'                                                         - "none"
#' @param biomart a data frame containing the Biomart table containing Transcript IDs, Gene IDs and Gene names#' 
#' @param transc_level a boolean indicating whether the data must be loaded at transcript- or gene-level
#' @param treat_batch_effect a boolean indicating whether a batch effect must be corrected or not
#' @export
#'
#' @return expression_matrix, a list containing a quant matrix, list stat (Pvalue, FC, gene/transcript name and ID) amd cpm,
#'
#'
#'
edgeR_stat = function(dge_list, samples_info,DE_method, padj_method, trans_level, biomart, treat_batch_effect){
  
  conditions = samples_info$dx
  batchs = samples_info$batch
  
  print("Exporatoty MDS")
  plotMDS(dge_list)
  
  print("Design Matrix")
  if(treat_batch_effect == "TRUE"){
    batchs = samples_info$batch
    design = model.matrix(~conditions+batchs)
  } else {
    design = model.matrix(~conditions)
  }
  
  print("Dispersion Estimation")
  dge_list = estimateDisp(dge_list, design)
  print(dge_list$common.dispersion)
  plotBCV(dge_list)
  
  print("Quant Matrix")
  quant_matrix = dge_list$counts
  colnames(quant_matrix) = dge_list$samples$ID
  quant_matrix = data.frame(quant_matrix)
  quant_matrix = cbind(rownames(quant_matrix), quant_matrix)
  
  print("cpms")
  cpms  = cpm(dge_list, dge_list = y$offset, log = FALSE)
  cpms = cbind(rownames(cpms), cpms)
  print("Differential Expression")
  if(DE_method == "ET"){
    print("Exact Test Fisher")
    et = exactTest(dge_list)
    toptags_samples = topTags(et, n=nrow(dge_list), sort="none")
    
  }else if(DE_method == "QLF"){
    print("Generalized Linear Model quasi-likelihood (QL) F-test")
    
    fit = glmQLFit(dge_list, design, robust = TRUE)
    lft = glmQLFTest(fit, coef = 2)
    toptags_samples = topTags(lft, n=nrow(dge_list), sort="none")
    
  }else if(DE_method == "LR"){
    print("Generalized Linear Model Likelihood ratio")
    fit = glmFit(dge_list, design, robust = TRUE)
    lrt = glmLRT(fit, coef = 2)
    toptags_samples = topTags(lrt, n=nrow(dge_list), sort="Pvalue")
  }
  
  print("statistical analysis")
  if(trans_level == TRUE){
    
    list_stat = data.frame("Transcript" = rownames(toptags_samples$table),
                           "Pvalue" = toptags_samples$table$PValue,
                           "FC" = toptags_samples$table$logFC)
    
    list_stat = list_stat[order(list_stat$Transcript),]
    list_stat$Transcript = sub("\\..*", "", list_stat$Transcript)
    biomart = biomart[order(biomart$Transcript_stable_ID),]
    list_stat$Gene_name = biomart[biomart$Transcript_stable_ID %in% list_stat$Transcript,]$Gene_name
    list_stat = list_stat[order(list_stat$Pvalue),]
    colnames(quant_matrix)[1] = "Transcript"
    colnames(cpms)[1] = "Transcript"  
    
  } else {
    list_stat = data.frame("Gene_ID" = rownames(toptags_samples$table),
                           "Pvalue" = toptags_samples$table$PValue,
                           "FC" = toptags_samples$table$logFC)
    
    # list_stat = list_stat[order(list_stat$Gene_ID),]
    # biomart = biomart[order(biomart$Gene_stable_ID),]
    # biomart_unique_genes = biomart[match(unique(biomart$Gene_stable_ID), biomart$Gene_stable_ID),]
    # list_stat$Gene_name = biomart_unique_genes[biomart_unique_genes$Gene_stable_ID %in% list_stat$Gene_ID,]$Gene_name
    # list_stat = list_stat[order(list_stat$Pvalue),]
    
    colnames(quant_matrix)[1] = "Gene"
    colnames(cpms)[1] = "Gene"  
  }
  
  print("calculate FDR")
  hist(list_stat$Pvalue, breaks=50)
  list_stat$FDR = rep(NA, times = length(list_stat$Pvalue))
  significant_results = list_stat[which(list_stat$Pvalue <= 0.05),]
  significant_results$FDR = p.adjust(significant_results$Pvalue, method = padj_method)
  
  list_stat$FDR  = replace(list_stat$FDR,which(list_stat$Pvalue <= 0.05),significant_results$FDR )
  
  expression_matrix = list(quant_matrix,list_stat, cpms)
  
  return(expression_matrix)
  
}

edgeR_stat(dge_list = dge, samples_info=cluster_metadata, DE_method="QLF",
           padj_method="BH", trans_level=FALSE, biomart=NaN, treat_batch_effect=TRUE)







