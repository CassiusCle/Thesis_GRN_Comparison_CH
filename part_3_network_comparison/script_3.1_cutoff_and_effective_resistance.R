################################################################################
### Name:        script_3.1_cutoff_and_effective_resistance.R
### Description: First script for network comparison.
###              - CutOff edges at thresholds
###              - Compute effective resistance distance
### Author:      Clemens Pieter Dick Hendrickx
################################################################################

################################################################################
# Load Packages
################################################################################
LoadPackages <- function(){
  if(!require(tidyverse))       {install.packages("tidyverse")};    library(tidyverse)
  if(!require(readr))           {install.packages("readr")};        library(readr)
  if(!require(foreach))         {install.packages("foreach")};      library(foreach)
  if(!require(doSNOW))         {install.packages("doSNOW")};      library(doSNOW)
  if(!require(doParallel))      {install.packages("doParallel")};   library(doParallel)
  if(!require(pracma))          {install.packages("pracma")};       library(pracma)
  if(!require(matrixcalc))        {install.packages("matrixcalc")};   library(matrixcalc)
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
cutOff <- 0

results_directory <- "C:\\Users\\clemenshendrickx\\Documents\\Thesis\\Thesis Data\\Results v1.0\\"

################################################################################
# Load data
################################################################################
Weight_Asthma_0 <- readRDS(paste0(dir, folder_name, 'Weight_matrix_','Asthma_',cutOff,'.rds'))
Weight_Control_0 <- readRDS(paste0(dir, folder_name, 'Weight_matrix_','Control_',cutOff,'.rds'))

if( all((rownames(Weight_Asthma_0) == colnames(Weight_Asthma_0))&
        all(rownames(Weight_Control_0) == colnames(Weight_Control_0))&
        all(rownames(Weight_Control_0) == rownames(Weight_Asthma_0)))){
  gene_names <- rownames(Weight_Asthma_0)
}

################################################################################
# Auxilliary functions
################################################################################
GetLaplacian <- function(matrix, validate = FALSE) {
  # Function that computes the (weighted) Laplacian matrix of the network
  #     N.B.: The adjacency matrix is a special case of the input weight matrix
  laplacian <- -matrix
  diag(laplacian) <- rowSums(matrix)
  
  if (validate){
    print('Checking if Laplacian satisfies properties outlined in Ellens (2011)')
    # Laplacian Matrix should be positive semi-definite
    cat('\tLaplacian is Positive Semi-Definite:', is.positive.semi.definite(laplacian),'\n')
    
    # Laplacian Matrix should be symmetric
    cat('\tLaplacian is symmetric:', isSymmetric(laplacian),'\n')
    
    # Laplacian Matrix rows should sum to zero
    cat('\tSum of absolute rowSums of Laplacian is effectively zero:', sum(abs(rowSums(laplacian))),'\n')
    
    # Eigenvalues are real, non-negative and smallest is 0
    L_eigen <- eigs_sym(laplacian, 2, sigma = 0)$values
    cat('\tEigenvalues of Laplacian are real:', sum(Im(eig(laplacian))) == 0,'\n')
    cat('\tEigenvalues of Laplacian are non-negative and smallest is 0, Minimum Eigenvalue =', min(L_eigen),'\n')
    cat('\tNetwork is connected, i.e. Algebraic Connectivity is not zero:', L_eigen[1],'\n')
  }
  return(laplacian)
}


weight2adjacency <- function(weights) {
  # Function that converts weight matrix to adjacency matrix
  (1.0*(weights != 0))
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
  
  list(asthma = Asthma_weights,
       control = Control_weights)
}

AnalyseLaplacian <- function(laplacian, tolerance = 1e6*.Machine$double.eps) {
  # Function that checks a network for certain properties by analysing its
  # laplacian matrix.
  
  number_nodes <- nrow(laplacian)
  number_edges <- sum(laplacian[lower.tri(laplacian)] != 0)
  is_psd <- is.positive.semi.definite(laplacian, tol = tolerance)
  is_symmetric <- isSymmetric(laplacian, tol = tolerance)
  is_rowsum_zero <- sum(abs(rowSums(laplacian))) < tolerance
  
  # Spectral measures
  eigen_vals <- eigen(laplacian, symmetric = is_symmetric, only.values = TRUE)$values
  is_eigen_real <- sum(Im(eigen_vals)) == 0
  min_eigen <- eigen_vals[length(eigen_vals)]
  algebraic_connectivity <- eigen_vals[length(eigen_vals)-1]
  total_effective_resistance <- length(eigen_vals)*sum(1/eigen_vals[1:(length(eigen_vals)-1)])
  
  output <- data.frame(numb_nodes = number_nodes,
                       numb_edges = number_edges,
                       is_PSD = is_psd,
                       is_symm = is_symmetric,
                       is_rowsum_zero = is_rowsum_zero,
                       is_eigen_real = is_eigen_real,
                       min_eigenvalue = min_eigen,
                       alg_conn = algebraic_connectivity,
                       total_ER = total_effective_resistance
                       )
  return(output)
}

SetUpParallelCluster <- function(configurations) {
  # Function that sets up cluster for parallel computation
  
  myCluster <- makeCluster(min(nrow(configurations), detectCores()-1))
  registerDoSNOW(myCluster)
  
  pb <- txtProgressBar(max = nrow(configurations), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  return(list(cluster = myCluster,
              opts = opts))
}

ComputeEffectiveResistance <- function(laplacian, is_symmetric = TRUE) {
  # Function that computes the effective resistance distances over a network 
  # using from the (weighted) laplacian matrix.
  # (Method described in Vos (2016))
  eigen_decomp <- eigen(laplacian, symmetric = is_symmetric)
  n <- nrow(laplacian)
  
  U <- eigen_decomp$vectors
  D_plus <- diag(c(1/eigen_decomp$values[1:(length(eigen_decomp$values)-1)],0))
  laplacian_pseudoinverse <- U %*% D_plus %*% t(U)
  effective_resistance <- laplacian_pseudoinverse %>% {(matrix(rep(diag(.), n), nrow = n, byrow = FALSE) +
                                                          matrix(rep(diag(.), n), nrow = n, byrow = TRUE) - 2*.)}
  rownames(effective_resistance) <- colnames(effective_resistance) <- rownames(laplacian)
  return(effective_resistance)
}

################################################################################
# Network comparison
################################################################################


########### Test properties of networks for different CutOff values  ###########
cutoff_options <- seq(1, 30, by = 0.5)
network_options <- c("mean") # Options: c(dir, mean, max, min)
weighted_options <- c(TRUE)

Weight_Asthma <- Weight_Asthma_0
Weight_Control <- Weight_Control_0

configurations <- data.frame(cutoff = cutoff_options, net = network_options, is_weighted = weighted_options)
configurations

network_weights <- list()
if("dir" %in% unique(configurations$net)){
  network_weights[["dir"]] <- list(asthma = Weight_Asthma,
                                   control = Weight_Control
                                   )
}
if("mean" %in% unique(configurations$net)){
  network_weights[["mean"]] <- list(asthma = (Weight_Asthma + t(Weight_Asthma))/2,
                                    control = (Weight_Control + t(Weight_Control))/2
                                    )
}
if("max" %in% unique(configurations$net)){
  network_weights[["max"]] <- list(asthma = pmax(Weight_Asthma, t(Weight_Asthma)),
                                   control = pmax(Weight_Control, t(Weight_Control))
                                   )
}
if("min" %in% unique(configurations$net)){
  network_weights[["min"]] <- list(asthma = pmin(Weight_Asthma_0, t(Weight_Asthma_0)),
                                   control = pmin(Weight_Control_0, t(Weight_Control_0))
                                   )
}

# Set up parallel cluster
cluster_list <- SetUpParallelCluster(configurations)

cutoff_analysis <- foreach(i = 1:nrow(configurations), 
                           .combine = "rbind", 
                           .packages = c("dplyr", "matrixcalc"),
                           .options.snow = cluster_list$opts) %dopar% {
  this_config <- configurations[i, ]
  
  this_asthma <- network_weights[[this_config$net]]$asthma
  this_control <- network_weights[[this_config$net]]$control
  
  this_networks <- reduceNetworksByCutOff(this_asthma, this_control, this_config$cutoff)
  this_asthma <- this_networks$asthma
  this_control <- this_networks$control
  
  if(!this_config$is_weighted) {
    this_asthma <- weight2adjacency(this_asthma)
    this_control <- weight2adjacency(this_control)
  }
  
  out_asthma <- AnalyseLaplacian(GetLaplacian(this_asthma))
  out_control <- AnalyseLaplacian(GetLaplacian(this_control))
  
  output <- cbind(this_config, out_asthma, out_control)
  colnames(output) <- c(colnames(this_config),
                        paste0("asthma_",colnames(out_asthma)),
                        paste0("control_",colnames(out_control))
  )
  return(output)
}
stopCluster(cluster_list$cluster)

cutoff_results_weighted <- cutoff_analysis %>% filter(net == "mean" & is_weighted)
cutoff_results_unweighted <- cutoff_analysis %>% filter(net == "mean" & !is_weighted)
config_vars <- colnames(cutoff_results_weighted) %in% c("cutoff", "net", "is_weighted")
cutoff_analysis_asthma <- cutoff_analysis[, config_vars | str_detect(colnames(cutoff_analysis), ".*asthma.*")]
cutoff_analysis_control <- cutoff_analysis[, config_vars | str_detect(colnames(cutoff_analysis), ".*control.*")]
colnames(cutoff_analysis_asthma) <- colnames(cutoff_analysis_control) <- str_replace(colnames(cutoff_analysis_asthma), "asthma_", "")
cutoff_analysis_asthma$healthStatus <- "Asthma" 
cutoff_analysis_control$healthStatus <- "Control" 
cutoff_analysis <- rbind(cutoff_analysis_asthma, cutoff_analysis_control)

plot(x = cutoff_analysis$cutoff, y = cutoff_analysis$numb_edges)

saveRDS(cutoff_analysis, paste0(results_directory, 'cutoff_analysis','_8_12_2021','.rds'))
remove(cutoff_results_weighted, 
       cutoff_results_unweighted, 
       config_vars, 
       cutoff_analysis_asthma, 
       cutoff_analysis_control)
cutoff_analysis <- readRDS(paste0(results_directory, 'cutoff_analysis','_8_12_2021','.rds'))


# Edges plots
ggplot(filter(cutoff_analysis, is_weighted == TRUE), 
       aes(x=cutoff, y = numb_edges, color = healthStatus)) + 
       geom_line(linetype = "dashed") + geom_point()

# Nodes plots
ggplot(filter(cutoff_analysis, is_weighted == TRUE), 
       aes(x=cutoff, y = numb_nodes, color = healthStatus)) + 
       geom_line(linetype = "dashed") + geom_point()

# Algebraic connectivity plots
ggplot(filter(cutoff_analysis, is_weighted == TRUE), 
       aes(x=cutoff, y = alg_conn, color = healthStatus)) + 
  geom_line(linetype = "dashed") + geom_point()

# Total effective resistance plots
ggplot(filter(cutoff_analysis, is_weighted == TRUE, cutoff <= 22), 
       aes(x=cutoff, y = total_ER, color = healthStatus)) + 
  geom_line(linetype = "dashed") + geom_point()


########### Calculate effective resistances for different CutOff values ########
cutoff_options <- seq(from = 5, to = 22, by = 0.5)
weighted_options <- c(TRUE, FALSE)

configurations <- crossing(cutoff = cutoff_options, is_weighted = weighted_options)
configurations

# Create Symmetric Matrices
Weight_Asthma_Mean <- (Weight_Asthma_0 + t(Weight_Asthma_0))/2
Weight_Control_Mean <- (Weight_Control_0 + t(Weight_Control_0))/2

network_weights <- list(asthma = Weight_Asthma_Mean,
                       control = Weight_Control_Mean
                       )


# Set up parallel cluster
cluster_list <- SetUpParallelCluster(configurations)

cutoff_deep_analysis <- foreach(i = 1:nrow(configurations), 
                           .packages = c("dplyr", "matrixcalc"),
                           .options.snow = cluster_list$opts) %dopar% {
                             this_config <- configurations[i, ]
                             
                             this_asthma <- network_weights$asthma
                             this_control <- network_weights$control
                             
                             this_networks <- reduceNetworksByCutOff(this_asthma, this_control, this_config$cutoff)
                             this_asthma <- this_networks$asthma
                             this_control <- this_networks$control
                             
                             if(!this_config$is_weighted) {
                               this_asthma <- weight2adjacency(this_asthma)
                               this_control <- weight2adjacency(this_control)
                             }
                             
                             out_asthma <- ComputeEffectiveResistance(GetLaplacian(this_asthma))
                             out_control <- ComputeEffectiveResistance(GetLaplacian(this_control))
                             
                             output <- list(cutoff = this_config$cutoff, 
                                            is_weighted = this_config$is_weighted, 
                                            ER_matrices = list(
                                              asthma = out_asthma,
                                              control = out_control
                                              ))
                             return(output)
                           }
stopCluster(cluster_list$cluster)

configurations$weight_situation <- sapply(configurations$is_weighted, FUN = function(x){if(x) {"weighted"} else {"unweighted"}})
names(cutoff_deep_analysis) <- paste0(configurations$weight_situation, "_", configurations$cutoff)

saveRDS(cutoff_deep_analysis, paste0(results_directory, 'cutoff_deep_analysis_8_12_2021.rds'))
