################################################################################
### Name:        script_1.3.3_scGNN_output_selection.R
### Description: Selects and saves the relevant subset of the imputed data based 
##               on the output of scGNN 
### Author:      Clemens Pieter Dick Hendrickx
################################################################################

################################################################################
# Load Packages
################################################################################

LoadPackages <- function(){
  
  if(!require(dplyr))           {install.packages("dplyr")};        library(dplyr)
  if(!require(readr))           {install.packages("readr")};        library(readr)
  if(!require(qlcMatrix))       {install.packages("qlcMatrix")};    library(qlcMatrix)
  if(!require(Rfast))           {install.packages("Rfast")};        library(Rfast)
  if(!require(janitor))         {install.packages("janitor")};      library(janitor)
  
  # if(!require(moments))         {install.packages("moments")};      library(moments)
  # if(!require(tseries))         {install.packages("tseries")};      library(tseries)
  # if(!require(pdist))           {install.packages("pdist")};        library(pdist)
  # if(!require(philentropy))     {install.packages("philentropy")};  library(philentropy)
  # if(!require(Rcpp))            {install.packages("Rcpp")};         library(Rcpp)
  }
LoadPackages()

gc()

dir <- "~/Thesis/Thesis Data/Selected Data/Intermediate data/QC_OUTPUT_minG.CA1000_CB20_GA10_GB5/"
num_genes <- 2000
excluding <- '_excl'
n_clusters <- '_clusters_10'

Use_expression <- read_csv(paste0(dir, 'scGNN_pre_output/ng_', num_genes, '_log', excluding, '/Use_expression.csv')) %>% as.data.frame
rownames(Use_expression) <- Use_expression[,'X1']
Use_expression[, 'X1'] <- NULL

scGNN_clustering <- read_csv(paste0(dir, 'scGNN_output/ng_', num_genes, '_log', n_clusters, excluding, '/ng_', num_genes, '_log', excluding, '_results.txt')) %>% as.data.frame
rownames(scGNN_clustering) <- scGNN_clustering[,'X1']

metadata <- read_csv(paste0(dir, 'metadata_normalised_counts_sconeCA1000_CB20_GA10_GB5.csv')) %>% as.data.frame
rownames(metadata) <- metadata[,'X1']

metadata_scGNN <- merge(x = scGNN_clustering, y = metadata, by = 'X1', all.x = TRUE, all.y = FALSE)

Celltypes_clusters <- tabyl(metadata_scGNN, Celltype, cluster) %>% print

Celltypes_healthStatus <- tabyl(metadata_scGNN, Celltype, healthStatus) %>% print

# old_clustering <- metadata %>% select(cluster, X1)
# clusterings <- merge(x = scGNN_clustering, y = old_clustering, by = 'X1', all.x = TRUE, all.y = FALSE)


# 
# compare_clustering <- function(clusters) {
#   old_types <- clusters %>% select(cluster) %>% table %>% names
#   k_range <- 0:(length(table(select(clusters, Celltype))) - 1)
#   
#   clustering_compare <- data.frame(NULL)
#   
#   for (k in k_range){
#     counts <- rep(0,length(old_types))
#     names(counts) <- old_types
#     
#     this_counts <- clusters %>% filter(Celltype == k) %>% select(cluster) %>% table
#     counts[names(this_counts)] <- this_counts 
#     clustering_compare <- rbind(clustering_compare, counts)
#   }
#   colnames(clustering_compare) <- old_types
#   rownames(clustering_compare) <- k_range
#   
#   return(clustering_compare)
# }
# comparison <- compare_clustering(clusterings) %>% print


Celltypes_clusters <- Celltypes_clusters[,-1]
desired_cols <- match(c('Basal', 'Secretory'), colnames(Celltypes_clusters))
clusters_selected <- c()
for(i in 1:nrow(Celltypes_clusters)){
  this_row <- Celltypes_clusters[i,]
  if (match(max(this_row), this_row) %in% desired_cols) {
    clusters_selected <- c(clusters_selected, i-1)
  }
}
clusters_selected <- unique(clusters_selected) %>% print


imputed_data <- read_csv(paste0(dir, '/scGNN_output/ng_', num_genes, '_log', n_clusters, excluding, '/ng_', num_genes, '_log', excluding, '_recon.csv')) %>% as.data.frame()
rownames(imputed_data) <- imputed_data[,'X1']
imputed_data[,'X1'] <- NULL

input_data <- Use_expression[rownames(imputed_data),colnames(imputed_data)]

cells_selected <- metadata[metadata_scGNN %>% filter(Celltype %in% clusters_selected) %>% select(X1) %>% unlist,]
control_group <- cells_selected %>% filter(healthStatus == 'Ctrl') %>% select(seqID) %>% unlist
asthma_group <- cells_selected %>% filter(healthStatus == 'PersA') %>% select(seqID) %>% unlist

expression_asthma_original  <- input_data[,asthma_group]
expression_asthma_imputed   <- imputed_data[,asthma_group]
expression_control_original <- input_data[,control_group]
expression_control_imputed  <- imputed_data[,control_group]

remove(cells_selected, metadata, scGNN_clustering, this_row, Use_expression,
       asthma_group, clusters_selected, control_group, desired_cols, i)


# Distribution statistics
describe_data <- function(data) {
  data <- data %>% as.matrix %>% as.vector
  par(mfrow=c(1,2))
  data %>% hist
  data %>% density %>% plot
  print(summary(data))
  print(paste('Variance:',var(data)))
  print(paste('Skewness:',skewness(data)))
  print(paste('Kurtosis:',kurtosis(data)))
  # print(jarque.bera.test(data))
}


# Compare fractions of zero --> Conclusion: similar fractions before and after, 
fraction_zeros <- function(data) {
  zeros <- data == 0
  print(sum(zeros)/ (nrow(data)*ncol(data)))
}


# Complete data
describe_data(expression_asthma_original  )
describe_data(expression_control_original )
describe_data(expression_asthma_imputed   )
describe_data(expression_control_imputed  )

expression_asthma_original %>% fraction_zeros * 100
expression_control_original %>% fraction_zeros * 100
expression_asthma_imputed %>% fraction_zeros * 100
expression_control_imputed %>% fraction_zeros * 100

# Totals of each gene
describe_data(rowSums(expression_asthma_original))
describe_data(rowSums(expression_control_original))
describe_data(rowSums(expression_asthma_imputed))
describe_data(rowSums(expression_control_imputed))


# Made zero inspection
rowSums(expression_asthma_original)[rowSums(expression_asthma_imputed) == 0] %>% mean
rowSums(expression_control_original)[rowSums(expression_control_imputed) == 0] %>% mean

intersect(names(rowSums(expression_asthma_original)[rowSums(expression_asthma_imputed) == 0]),
      names(rowSums(expression_control_original)[rowSums(expression_control_imputed) == 0])) %>% length

madeZero_asthma_names <- names(rowSums(expression_asthma_original)[rowSums(expression_asthma_imputed) == 0])
madeZero_control_names <- names(rowSums(expression_control_original)[rowSums(expression_control_imputed) == 0])

madeZero_control_rest_names <- setdiff(madeZero_control_names, madeZero_asthma_names)


madeZero_asthma <- expression_asthma_original[madeZero_asthma_names,]
madeZero_control <- expression_control_original[madeZero_asthma_names,]

madeZero_control_rest <- expression_control_original[madeZero_control_rest_names,]

describe_data(madeZero_asthma)
describe_data(madeZero_control)


describe_data(expression_asthma_original[madeZero_control_rest_names,])
describe_data(madeZero_control_rest)

describe_data(madeZero_asthma[madeZero_asthma != 0])
describe_data(madeZero_control[madeZero_control != 0])


describe_data(madeZero_control_rest[madeZero_control_rest != 0])




# Save as CSV


folder_name <- paste0('scGNN_output/ng_', num_genes, '_log', n_clusters, excluding, '/')

paste0(dir, folder_name, deparse(substitute(expression_asthma_original)),'.csv')

save_as_csv <- function(data) {
  write.csv(data, file = paste0(dir, folder_name, deparse(substitute(data)),'.csv'))
}

save_as_csv(expression_asthma_original)
save_as_csv(expression_control_original)
save_as_csv(expression_asthma_imputed)
save_as_csv(expression_control_imputed)

