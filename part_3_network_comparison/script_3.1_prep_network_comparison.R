################################################################################
### Name:        script_3.1_prep_network_comparison.R
### Description: Preparation of network comparison
###              - Convert GRNBoost2 output (edge list) to weight matrix
###              - Create .csv edge list for use in Cytoscape
###                (used for creating network visualisations)
### Author:      Clemens Pieter Dick Hendrickx
################################################################################

################################################################################
# Load Packages
################################################################################
LoadPackages <- function(){
  if(!require(tidyverse))       {install.packages("tidyverse")};    library(tidyverse)
  if(!require(readr))           {install.packages("readr")};        library(readr)
  if(!require(foreach))         {install.packages("foreach")};      library(foreach)
  if(!require(doParallel))      {install.packages("doParallel")};   library(doParallel)
  if(!require(pracma))          {install.packages("pracma")};       library(pracma)
  if(!require(maotai))          {install.packages("maotai")};       library(maotai)
}
LoadPackages()

gc()

################################################################################
# Hyperparameters
################################################################################
dir <- "~/Thesis/Thesis Data/Selected Data/Intermediate data/QC_OUTPUT_minG.CA1000_CB20_GA10_GB5/"
num_genes <- 2000
excluding <- '_excl'
n_clusters <- '_clusters_10'

folder_name <- paste0('scGNN_output/ng_', num_genes, '_log', n_clusters, excluding, '/')

seedA <- 1913
seedC <- 1913
impute_status <- 'imputed'

################################################################################
# Load data
################################################################################
df_Asthma <- read_csv(paste0(dir, folder_name, 'network_asthma_',impute_status,'_',seedA, '.csv')) %>% as.data.frame
colnames(df_Asthma)[colnames(df_Asthma) == 'importance'] <- 'weight'

df_Control <- read_csv(paste0(dir, folder_name, 'network_control_',impute_status,'_',seedC, '.csv')) %>% as.data.frame
colnames(df_Control)[colnames(df_Control) == 'importance'] <- 'weight'

################################################################################
# Auxilliary functions
################################################################################

makeComparable <- function(A, C) {
  # Function that makes edgelists of networks comparable by ensuring that both
  #   networks consist of the same sets of nodes.
  includedGenes <- intersect(unique(c(A$TF, A$target)),
                             unique(c(C$TF, C$target)))
  
  A <- A %>% filter(TF %in% includedGenes & target %in% includedGenes)
  C <- C %>% filter(TF %in% includedGenes & target %in% includedGenes)
  return(list(Asthma = A,
              Control = C
  )
  )
}

edges2adjacency <- function(edge_df, parallel = FALSE) {
  # Function that converts edge DataFrame into adjacency matrix
  nodes <- union(edge_df$TF, edge_df$target) %>% sort
  base_row <- rep(0, length(nodes))
  names(base_row) <- nodes
  
  if (!parallel){
    result <- foreach(node =  nodes, .combine = rbind) %do% {
      this_row <- base_row
      this_row[(filter(edge_df, TF == node))$target] <- 1
      return(this_row)
    }
  } else {
    n.cores <- parallel::detectCores() - 1
    my.cluster <- parallel::makeCluster(
      n.cores,
      type = "PSOCK"
    )
    doParallel::registerDoParallel(cl = my.cluster)
    foreach::getDoParRegistered()
    result <- foreach(node =  nodes, .combine = rbind) %dopar% {
      this_row <- base_row
      this_row[(filter(edge_df, TF == node))$target] <- 1
      return(this_row)
    }
    parallel::stopCluster(cl = my.cluster)
  }
  rownames(result) <- nodes
  return(result)
}


edges2weight <- function(edge_df) {
  # Function that converts edge DataFrame into weight matrix
  nodes <- union(edge_df$TF, edge_df$target) %>% sort
  base_row <- rep(0, length(nodes))
  names(base_row) <- nodes
  result <- foreach(node =  nodes, .combine = rbind) %do% {
    this_row <- base_row
    this_filt <- (filter(edge_df, TF == node))
    this_row[this_filt$target] <- this_filt$weight
    return(this_row)
  }
  
  rownames(result) <- nodes
  return(result)
}

reduceNetworksByCutOff <- function(Asthma_weights, Control_weights, CutOff = 0) {
  # Function that reduces the size of the networks by eliminating edges with weight
  # below CutOff value. Nodes without edges in one or both networks are automatically
  # removed.
  Asthma_weights[Asthma_weights < CutOff] <- 0
  Control_weights[Control_weights < CutOff] <- 0
  
  not_exclude <- !((rowSums(Asthma_weights) == 0 & colSums(Asthma_weights) == 0) | 
                     (rowSums(Control_weights) == 0 & colSums(Control_weights) == 0)
  )
  
  Asthma_weights <- Asthma_weights[not_exclude, not_exclude]
  Control_weights <- Control_weights[not_exclude, not_exclude]
  
  list(Asthma = Asthma_weights,
       Control = Control_weights)
}

Edges2MeanWeightPlusCutOff <- function(edges_df, cutOff) {
  # Function that removes directionality from network and removes edges under
  # CutOff threshold through manipulating the edgelist
  # N.B.: Only used to construct edgelist as input for CytoScape (visualisation)
  edges_toMerge <- edges_df
  colnames(edges_toMerge) <- c("index", "target", "source", "weight", "interaction")
  edges_merged <- merge(x = edges_df, y = edges_toMerge, by = c("source", "target"), all.x = TRUE, all.y = FALSE)
  edges_merged[is.na(edges_merged$weight.y), ]$weight.y <- 0
  edges_merged$weight <- (edges_merged$weight.x + edges_merged$weight.y)/2
  edges_merged <- edges_merged %>% filter(weight > cutOff) %>% {.[, c("source", "interaction.x", "weight", "target")]}
  colnames(edges_merged) <- c("source", "interaction", "weight", "target")
  
  return(edges_merged)
}

################################################################################
# Prepare network data for network comparison
################################################################################

networks <- makeComparable(df_Asthma, df_Control)
edges_Asthma <- networks$Asthma
edges_Control <- networks$Control

edges_Asthma$weight %>% hist

colnames(edges_Asthma) <- c('index', 'source', 'target', 'weight')
colnames(edges_Control) <- c('index', 'source', 'target', 'weight')

edges_Asthma$interaction <- 'pred' 
edges_Control$interaction <- 'pred' 


# Create weight matrices
Weight_Asthma_0 <- edges2weight(edges_Asthma)
Weight_Control_0 <- edges2weight(edges_Control)

# Save Weight matrices
saveRDS(Weight_Asthma_0, paste0(dir, folder_name, 'Weight_matrix_','Asthma_0','.rds'))
saveRDS(Weight_Control_0, paste0(dir, folder_name, 'Weight_matrix_','Control_0','.rds'))



############ Save edge list for use in Cytoscape for visualisation #############
edges_Asthma_10_mean <- Edges2MeanWeightPlusCutOff(edges_Asthma, 10)
edges_Control_10_mean <- Edges2MeanWeightPlusCutOff(edges_Control, 10)

write_csv(edges_Asthma_10_mean, paste0(dir, folder_name, 'edges_','Asthma_10_mean','.csv'))
write_csv(edges_Control_10_mean, paste0(dir, folder_name, 'edges_','Control_10_mean','.csv'))












