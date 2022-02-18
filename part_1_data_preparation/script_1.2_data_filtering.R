################################################################################
### Name:        script_1.2_data_filtering.R
### Description: Perform data filtering steps
### Author:      Clemens Pieter Dick Hendrickx
################################################################################

################################################################################
# Load Packages
################################################################################

LoadPackages <- function(){
  if(!require(dplyr))     {install.packages("dplyr")};        library(dplyr)
  if(!require(cowplot))   {install.packages("cowplot")};      library(cowplot)
  
  # Install procedure via BioManager
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if(!require(scater))    {BiocManager::install("scater")};   library(scater)
}
LoadPackages()

################################################################################
# Setup of hyperparameters
################################################################################
dir <- "~/Thesis/Thesis Data/Selected Data/Intermediate data/"

selected_healthStatus <- c('PersA', 'Ctrl') # Options: 'PersA', 'Ctrl', 'ClinR', 'ComR'
filt_crit <- list(# Cell filtering criteria
                  minGenesDetected = 1000,   # CA
                  maxMitoPercentage = 20,    # CB

                  # Gene filtering criteria
                  minCellsExpressed = 10,    # GA
                  minHighExpressedCells = 5  # GB  - Used in Cole et al. (2019)
                  )

################################################################################
# Load data
################################################################################
gc()

counts <- readRDS(paste0(dir, 'expression_counts.rds'))
metadata <- readRDS(paste0(dir, 'expression_metadata.rds'))

# Filter on Asthmatic and Control group (from 4693 obs. to 3094 obs.)
metadata <- metadata %>% filter(healthStatus %in% selected_healthStatus)
counts <- counts[,metadata$seqID]

################################################################################
# Quality control
################################################################################
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  rowData = data.frame(gene_names = rownames(counts)),
  colData = data.frame(metadata)
)

# Identify mitochondrial genes
mito_genes <- rownames(sce)[grep("^MT\\.", rownames(sce))]
sce <- addPerCellQC(sce, subsets = list(Mito = mito_genes))

#----------- Filtering ---------------------------------------------------------
# Two types of filtering:
#   + (I)   Cell filtering: remove low quality cells
#             Consists of filters based on:
#             - Minimum number of genes detected per cell
#             - Percentage of Mitochondrial mRNA
#             N.B.: Has to be done before gene filtering because low quality genes may only be expressed in low quality cells
# 
#   + (II)  Gene filtering: remove low quality genes
#             Consists of filters based on:
#             - Minimum number of cells to be detected in
#             - Minimum number of cells to be highly expressed in

#---Cell filtering---

# Number of genes detected in each cell
range(sce$detected) # For every cell, between 1001 and 7970 genes detected
selected_cells_detection <- sce$detected >= filt_crit[['minGenesDetected']]

# Mitochondrial filtering
percentage_Mito <- sce@colData$subsets_Mito_percent

#   Max percentage Mitochondrial mRNA
print(paste(sum(percentage_Mito < filt_crit[['maxMitoPercentage']])/length(percentage_Mito)*100, 'percent of cells retained',
            'with maxMitoPercentage =', filt_crit[['maxMitoPercentage']]))
selected_cells_mito <- sce$subsets_Mito_percent < filt_crit[['maxMitoPercentage']]
cell_filter <- selected_cells_detection & selected_cells_mito

sce.filt <- sce[, cell_filter]
dim(sce.filt) %>% print

#---Gene filtering---

#   Minimum number of genes detected per cell 
selected_genes_min_expr <- rowSums(counts(sce.filt) != 0) >= filt_crit[['minCellsExpressed']]

# Minimum number of highly expressed cells
num_reads  <-  quantile(assay(sce.filt)[assay(sce.filt) > 0])[4] # Define "High" expression as upper quartile
selected_genes_highly_expressed <-  rowSums(assay(sce.filt) >= num_reads ) >= filt_crit[['minHighExpressedCells']]

gene_filter <- selected_genes_min_expr & selected_genes_highly_expressed

sce.filt <- sce.filt[gene_filter,]
dim(sce.filt) %>% print # 14548 genes by 3094 cells

################################################################################
# Set-up for Scone normalisation comparison
################################################################################
print('Metadata features')
sce.filt@colData %>% names %>% print

run <- factor(colData(sce.filt)$run)
seqRun <- factor(colData(sce.filt)$seqRun)
bio  <- factor(colData(sce.filt)$healthStatus)

# Remove intermediate objects
remove(counts, metadata, sce, cell_filter, gene_filter, mito_genes,
       num_reads, percentage_Mito, selected_cells_detection, selected_cells_mito,
       selected_genes_highly_expressed, selected_genes_min_expr, 
       selected_healthStatus, LoadPackages)

# Save metadata and count objects
filt_used <- paste0('CA',  filt_crit[['minGenesDetected']],
                    '_CB', filt_crit[['maxMitoPercentage']],
                    '_GA', filt_crit[['minCellsExpressed']],
                    '_GB', filt_crit[['minHighExpressedCells']])
folder_name <- paste0('QC_OUTPUT_minG.',filt_used)
dir.create(paste0(dir, folder_name))

saveRDS(sce.filt, paste0(dir, folder_name, '/', 'scExperiment_QC',filt_used,'.rds'))
saveRDS(bio,      paste0(dir, folder_name, '/', 'scone_meta_bio',filt_used,'.rds'))
saveRDS(run,      paste0(dir, folder_name, '/', 'scone_meta_run',filt_used,'.rds'))
saveRDS(seqRun,   paste0(dir, folder_name, '/', 'scone_meta_seqRun',filt_used,'.rds'))

