################################################################################
### Name:        script_1.2.X_normalisation_hpc.R
### Description: performs comparison of normalisation procedures using 
###              scone workflow on High Performance Computing Cluster
### Author:      Clemens Pieter Dick Hendrickx
################################################################################

################################################################################
# Load Packages
################################################################################

LoadPackages <- function(){
  
  # sconeReport requirements
  if(!require(NMF))             {install.packages("NMF")};          library(NMF)
  if(!require(plotly))          {install.packages("plotly")};       library(plotly)
  if(!require(visNetwork))      {install.packages("visNetwork")};   library(visNetwork)
  
  if(!require(dplyr))           {install.packages("dplyr")};        library(dplyr)
  
  
  # Install procedure via BioManager
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if(!require(scater))          {BiocManager::install("scater")};         library(scater)
  if(!require(scone))           {BiocManager::install("scone")};          library(scone)
  if(!require(SCnorm))          {BiocManager::install("SCnorm")};         library(SCnorm)
  if(!require(scran))           {BiocManager::install("scran")};          library(scran)
  if(!require(BiocParallel))    {BiocManager::install("BiocParallel")};   library(BiocParallel)
  
}
LoadPackages()

################################################################################
# Setup of hyperparameters
################################################################################

# Get arguments if running on HPC cluster
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1] 
filt_used <- args[2] 
batchChoice <- args[3] 
threads <- as.numeric(args[4])

# Decomment if running locally
# dir <- "~/Thesis/Thesis Data/Selected Data/Intermediate data/"
# filt_used <- 'CA1000_CB20_GA10_GB5'
# batchChoice <- 'run' #seqRun'

################################################################################
# Auxilliary functions
################################################################################

################################################################################
# Load data
################################################################################
gc()

folder_name <- paste0('QC_OUTPUT_minG.', filt_used)

sce <-    readRDS(paste0(dir, folder_name, '/', 'scExperiment_QC',filt_used,'.rds'))
bio <-    readRDS(paste0(dir, folder_name, '/', 'scone_meta_bio',filt_used,'.rds'))
run <-    readRDS(paste0(dir, folder_name, '/', 'scone_meta_run',filt_used,'.rds'))
seqRun <- readRDS(paste0(dir, folder_name, '/', 'scone_meta_seqRun',filt_used,'.rds'))

################################################################################
# Scone set-up
################################################################################

# Creating a SconeExperiment Object (batch options: run or seqRun)
sconeExp <- SconeExperiment(as.matrix(assay(sce)), 
                            bio = bio, 
                            batch = eval(parse(text = batchChoice)), 
                            is_log = FALSE
                            )
#----------- Definining normalisation methods ----------------------------------

# User-defined normalisation: Divide cell count by number of detected genes
# Does not really make sense for imputed data if done later, but now not imputed yet so we can use?
EFF_FN <- function (ei) {
  sums <- colSums(ei > 0)
  eo <- t(t(ei)*sums/mean(sums))
  return(eo)
}

# SCnorm normalisation
# SCNORM_FN <- function (ei) {
#   SCnorm(ei,
#          Conditions = rep(c(1), each = ncol(ei)),
#          FilterCellNum = 10,
#          PrintProgressPlots = FALSE,
#          NCores = 1#11#threads - 1
#          )
# }

# Define the list of functions used for normalisation
scaling_list  <- list(none = identity # Identity - do nothing
                      ,eff = EFF_FN # User-defined function
                      ,sum = SUM_FN # SCONE library wrappers...
                      ,scran = SCRAN_FN
                      ,tmm = TMM_FN 
                      ,uq = UQ_FN
                      ,fq = FQT_FN
                      ,deseq = DESEQ_FN
                      ,clr = CLR_FN
                      # Scater SCNorm
                      #,scnorm = SCNORM_FN
)

sconeExp <- scone(sconeExp,
                  scaling = scaling_list,
                  k_qc = 0,  # Do not perform QC factor adjustment (We have no qc argument)
                  k_ruv = 0, # Do not perform RUV adjustment (We have no negative controls)
                  adjust_bio = "no", # Biological origin is only used as for evaluation purposes
                  adjust_batch = "yes", # Batch information only used for evaluation purposes
                  run = FALSE # Do not run normalisation and evaluation yet
                  )

# Inspect scone experiment
get_params(sconeExp) %>% print
apply(get_params(sconeExp),2,unique) %>% print

################################################################################
# Running scone
################################################################################


#----Performance metrics:
# https://www.bioconductor.org/packages/devel/bioc/vignettes/scone/inst/doc/sconeTutorial.html

# BIO_SIL:    +                                              (Positive signature)
# (PCA)       Preservation of Biological Difference. The average silhouette width 
#             of clusters defined by bio, defined with respect to a Euclidean distance 
#             metric over the first 3 expression PCs.

# BATCH_SIL:  -                                              (Negative signature)
# (PCA)       Removal of Batch Structure. The average silhouette width of clusters 
#             defined by batch, defined with respect to a Euclidean distance 
#             metric over the first 3 expression PCs.

# PAM_SIL:    +                                              (Positive signature)
# (PCA)       Preservation of Single-Cell Heterogeneity. The maximum average silhouette 
#             width of clusters defined by PAM clustering, defined with respect to a 
#             Euclidean distance metric over the first 3 expression PCs.

# RLE_MED:    -                                              (Negative signature)
#             Reduction of Global Differential Expression. The mean squared-median 
#             Relative Log Expression (RLE).

# RLE_IQR:    -                                              (Negative signature)
#             Reduction of Global Differential Variability. The variance of 
#             the inter-quartile range (IQR) of the RLE.

# Not Used
# X EXP_QC_COR: -                                              (Negative signature)
#             Removal of Alignment Artifacts. R^2 measure for regression of first 
#             3 expression PCs on first k_qc QPCs. 

# X EXP_UV_COR: -                                              (Negative signature)
#             Removal of Expression Artifacts. R^2 measure for regression of first 3 
#             expression PCs on first 3 PCs of the negative control (specified by 
#             eval_negcon or ruv_negcon by default) sub-matrix of the original 
#             (raw) data.

# X EXP_WV_COR: +                                              (Positive signature)
#             Preservation of Biological Variance. R^2 measure for regression of 
#             first 3 expression PCs on first 3 PCs of the positive control (specified 
#             by eval_poscon) sub-matrix of the original (raw) data.

BiocParallel::register(
  #BiocParallel::SerialParam()
  #BiocParallel::SnowParam() # Use when on Windows
  BiocParallel::MulticoreParam(workers= threads - 1)
) # Register BiocParallel

registered() %>% print

sconeExp <- scone(sconeExp,
                  scaling = scaling_list,
                  k_qc = 0,  # Do not perform QC factor adjustment (We have no qc argument)
                  k_ruv = 0, # Do not perform RUV adjustment (We have no negative controls)
                  adjust_bio = "no", # Biological origin is only used as for evaluation purposes
                  adjust_batch = "no", # Batch information only used for evaluation purposes
                  run = TRUE,
                  eval_kclust = 2:10, # run: max = 74, seqRun: max = 4
                  stratified_pam = TRUE,
                  stratified_rle = TRUE,
                  return_norm = "in_memory",
                  zero = "postadjust",
                  verbose = TRUE
                  )

print('Finished running scone!')

# View Metric Scores
get_scores(sconeExp) %>% print

# View Mean Score Rank
get_score_ranks(sconeExp) %>% print

saveRDS(sconeExp, paste0(dir, folder_name, '/', 'scExperiment_QC_normalised_2.0_',filt_used,'.rds'))
