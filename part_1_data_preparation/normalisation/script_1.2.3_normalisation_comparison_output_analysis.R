################################################################################
### Name:        script_1.2.X_normalisation_comparison_output_analysis.R
### Description: Analysis of the output of the normalisation comparison
### Author:      Clemens Pieter Dick Hendrickx
################################################################################

################################################################################
# Load Packages
################################################################################

LoadPackages <- function(){
  
  # sconeReport requirements
  # if(!require(NMF))             {install.packages("NMF")};          library(NMF)
  # if(!require(plotly))          {install.packages("plotly")};       library(plotly)
  # if(!require(visNetwork))      {install.packages("visNetwork")};   library(visNetwork)
  if(!require(dplyr))           {install.packages("dplyr")};        library(dplyr)
  # Install procedure via BioManager
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if(!require(scater))          {BiocManager::install("scater")};         library(scater)
  if(!require(scone))           {BiocManager::install("scone")};          library(scone)
}
LoadPackages()

################################################################################
# Setup of hyperparameters
################################################################################
dir <- "~/Thesis/Thesis Data/Selected Data/Intermediate data/"

filt_used <- 'CA1000_CB20_GA10_GB5'
folder_name <- paste0('QC_OUTPUT_minG.', filt_used)

# batchChoice <- 'run' #seqRun'

################################################################################
# Auxilliary functions
################################################################################

######## Read in data from Scone normalisation comparison
sconeExp <-    readRDS(paste0(dir, folder_name, '/', 'scExperiment_QC_normalised_2.0_',filt_used,'.rds'))

# View Metric Scores
get_scores(sconeExp)

# View Mean Score Rank
get_score_ranks(sconeExp)


# Extract normalized data from top method
out_norm  <-  get_normalized(sconeExp,
                             method = rownames(get_params(sconeExp))[1])

# Save RDS
saveRDS(out_norm, paste0(dir, folder_name, '/', 'normalised_counts_scone',filt_used,'.rds'))
sce <-    readRDS(paste0(dir, folder_name, '/', 'scExperiment_QC',filt_used,'.rds'))
normcounts(sce) <- out_norm
saveRDS(sce, paste0(dir, folder_name, '/', 'sce_normalised_incl',filt_used,'.rds'))

# Compute descriptive statistics

mean(counts(sce) == 0) # Fraction of zero's (sparsity) before normalisation
mean(normcounts(sce) == 0) # Fraction of zero's (sparsity) normalised

counts_PersA <- counts(sce)[,sce$healthStatus == 'PersA']
normcounts_PersA <- normcounts(sce)[,sce$healthStatus == 'PersA']

counts_Ctrl <- counts(sce)[,sce$healthStatus == 'Ctrl']
normcounts_Ctrl <- normcounts(sce)[,sce$healthStatus == 'Ctrl']

all(rownames(counts_PersA) %in% rownames(normcounts_PersA))
all(colnames(counts_PersA) %in% colnames(normcounts_PersA))
all(rownames(counts_Ctrl) %in% rownames(normcounts_Ctrl))
all(colnames(counts_Ctrl) %in% colnames(normcounts_Ctrl))


mean(counts_PersA == 0) # Fraction of zero's (sparsity) after gene filtering asthma
mean(normcounts_PersA == 0) # Fraction of zero's (sparsity) after gene filtering asthma

mean(counts_Ctrl == 0) # Fraction of zero's (sparsity) after gene filtering control
mean(normcounts_Ctrl == 0) # Fraction of zero's (sparsity) after gene filtering control

mean(colMeans(counts_PersA)) # Average read count per Asthma cell
mean(colMeans(counts_Ctrl)) # Average read count per Control cell
mean(colMeans(counts(sce))) # Average read count total

mean(colMeans(normcounts_PersA)) # Average read count per Asthma cell
mean(colMeans(normcounts_Ctrl)) # Average read count per Control cell
mean(colMeans(normcounts(sce))) # Average read count total

sd(colMeans(counts_PersA)) # Average read count per Asthma cell
sd(colMeans(counts_Ctrl)) # Average read count per Control cell
sd(colMeans(counts(sce))) # Average read count total

sd(colMeans(normcounts_PersA)) # Average read count per Asthma cell
sd(colMeans(normcounts_Ctrl)) # Average read count per Control cell
sd(colMeans(normcounts(sce))) # Average read count total


metadata$run %>% table # Cells per run
metadata %>% select(healthStatus, run) %>% table # Cells per run per healthStatus

# Fraction of zero entries
(sum(out_norm == 0)/(nrow(out_norm)*ncol(out_norm))) %>% print

exclude_ionocytes_and_mixed_immune <- TRUE
if(exclude_ionocytes_and_mixed_immune) {
  
  subset <- sce@colData %>% data.frame %>% filter(cluster != 'Ionocytes' & cluster != 'Mixed_Immune') %>% select(seqID) %>% unlist
  
  # Save as CSV
  write.csv(normcounts(sce[,subset]), file = paste0(dir, folder_name, '/', 'normalised_counts_scone',filt_used,'_excluding.csv'))
  
  write.csv(colData(sce[,subset])[,c('index', 'donor', 'healthStatus', 'sex', 'cluster', 'seqRun', 'seqID', 'run', 'lane')],
            file = paste0(dir, folder_name, '/', 'metadata_normalised_counts_scone',filt_used,'_excluding.csv'))
  
} else {
  
  # Save as CSV
  write.csv(out_norm, file = paste0(dir, folder_name, '/', 'normalised_counts_scone',filt_used,'.csv'))
  
  write.csv(colData(sce)[,c('index', 'donor', 'healthStatus', 'sex', 'cluster', 'seqRun', 'seqID', 'run', 'lane')],
            file = paste0(dir, folder_name, '/', 'metadata_normalised_counts_scone',filt_used,'.csv'))
  
}

bio <-    readRDS(paste0(dir, folder_name, '/', 'scone_meta_bio',filt_used,'.rds'))
run <-    readRDS(paste0(dir, folder_name, '/', 'scone_meta_run',filt_used,'.rds'))
seqRun <- readRDS(paste0(dir, folder_name, '/', 'scone_meta_seqRun',filt_used,'.rds'))

seqRun %>% unique
run %>% unique


