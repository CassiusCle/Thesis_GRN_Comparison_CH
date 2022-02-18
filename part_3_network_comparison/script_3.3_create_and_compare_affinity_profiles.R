

################################################################################
### Name:        script_3.3_create_and_compare_affinity_profiles.R
### Description: Create and compare gene affinity profiles
###              - Convert effective resistance distances to relative
###              - Transform using RBF kernel
###              - Analyse results and visualise using boxplot
### Author:      Clemens Pieter Dick Hendrickx
################################################################################

LoadPackages <- function(){
  if(!require(tidyverse))       {install.packages("tidyverse")};    library(tidyverse)
  if(!require(readr))           {install.packages("readr")};        library(readr)
  if(!require(rootSolve))       {install.packages("rootSolve")};    library(rootSolve)
  if(!require(extraDistr))      {install.packages("extraDistr")};   library(extraDistr)
  if(!require(foreach))         {install.packages("foreach")};      library(foreach)
  if(!require(doSNOW))          {install.packages("doSNOW")};       library(doSNOW)
  if(!require(doParallel))      {install.packages("doParallel")};   library(doParallel)
  if(!require(matrixStats))     {install.packages("matrixStats")};  library(matrixStats)
  if(!require(numDeriv))        {install.packages("numDeriv")};     library(numDeriv)
}
LoadPackages()

gc()

################################################################################
# Hyperparameters
################################################################################
results_directory <- "C:\\Users\\clemenshendrickx\\Documents\\Thesis\\Thesis Data\\Results v1.0\\"
cutoff_options <- seq(from = 5, to = 22, by = 0.5)
weighted_options <- TRUE #c(TRUE, FALSE)


################################################################################
# Load data
################################################################################
cutoff_deep_analysis <- readRDS(paste0(results_directory, 'cutoff_deep_analysis_8_12_2021.rds'))

################################################################################
# Auxilliary functions
################################################################################

Convert2Relative <- function(ER_matrices){
  # Function that converts the effective resistance distances to relative on range [0, 1]
  R_Asthma_normalised <- ER_matrices$asthma/rowSums(as.matrix(ER_matrices$asthma))
  R_Control_normalised <- ER_matrices$control/rowSums(as.matrix(ER_matrices$control))
  return(list(R_asthma_normalised = R_Asthma_normalised, R_control_normalised = R_Control_normalised))
}

RBF <- function(x, s = 1) {
  # Radial basis funciton kernel function
  exp(-(x^2)/(2*s^2))
}

IsMonotonicSeries <- function(series) {
  # Function that checks if a series is strictly increasing or decreasing
  is_strictly_decreasing <- all(diff(series) >= -.Machine$double.eps) 
  is_strictly_increasing <- all(diff(series) <= .Machine$double.eps)
  return(is_strictly_increasing | is_strictly_decreasing)
}

transformUsingRBFKernel <- function(matrix, sigma, 
                                    minVal = 2.2250738585072E-308){
  # Function that transforms input matrix using Radial Basis Function kernel
  
  transformedMatrix <- RBF(matrix, s = sigma) %>% pmax(., minVal)
  rownames(transformedMatrix) <- colnames(transformedMatrix) <- rownames(matrix)
  return(transformedMatrix)
}

computeSimilarityRBF <- function(R_Asthma, R_Control, 
                                 sigma = NULL,
                                 nbhd_size = 100,
                                 factor_increase = 5,
                                 minVal = 2.2250738585072E-308,
                                 verbose = FALSE
) {
  # Function that transforms ER distances for asthma and control group using 
  # an Radial Basis Function kernel
  
  if(is.null(sigma)){
    complete_data <- rbind(R_Asthma, R_Control)
    nbhd_target <- mean(apply(complete_data, MARGIN = 1, FUN = function(this_row) {
      sort(this_row, decreasing = FALSE)[nbhd_size]
    }))
    complete_data <- c(R_Asthma[upper.tri(R_Asthma)], 
                       R_Control[upper.tri(R_Control)]
    )
    hist_ER <- hist(complete_data, breaks = 10000, plot = FALSE)
    sigma <- uniroot(function(s){
      prob_density <- RBF(hist_ER$mids, s)
      distribution <- prob_density*hist_ER$counts
      density_percentage <- sum(distribution[hist_ER$mids <= nbhd_target])/sum(distribution)
      populations_percentage <- sum(hist_ER$counts[hist_ER$mids <= nbhd_target])/sum(hist_ER$counts)
      return(density_percentage/populations_percentage - factor_increase)
    }, c(.Machine$double.eps^(1/4), 5), tol = .Machine$double.eps)$root # / 10 for mean

    is_gradient_monotonic <- IsMonotonicSeries(grad(function(x){RBF(x, s=sigma)}, hist_ER$mids))
    
    if (verbose){
      print(paste("Gradient is monotonic:", is_gradient_monotonic))
      print(paste("Sigma found:", sigma))
    }
    
  }
  
  if(is_gradient_monotonic){
    P_Asthma <- transformUsingRBFKernel(R_Asthma, sigma, minVal)
    P_Control <- transformUsingRBFKernel(R_Control, sigma, minVal)
    return(list(Similarities_asthma = P_Asthma, Similarities_control = P_Control))
  } else {
    return(list(Similarities_asthma = 0, Similarities_control = 0))
  }
}

AnalyseTS <- function(values, 
                      name=NULL, 
                      plot = TRUE,
                      out_length = 10,
                      roundOutput = TRUE) {
  # Function that takes as input a vector of Test Statistic scores and outputs
  # the highest scores and plots the distribution
  if(values[1] == "Gradient not monotonic"){
    return("Gradient not monotonic")
  } else {
    values <- values %>% sort(decreasing = TRUE)
    out <- values[1:out_length]
    if(roundOutput){
      out <- round(out, digits=2)
    }
    if(!plot){
      return(out)
    } else if(!is.null(name)) {
      plot(values, main = paste0("Gene network profile analysis - ",name))
    } else {
      plot(values)
    }
    return(out)
  }
}


GetL1Distribution <- function(delta_matrix) {
  # Function that computes the L1 norm over all rows for an input matrix
  # of differences.
  if(sum(delta_matrix) == 0) {
    return("Gradient not monotonic")
  } else {
    values <- rowSums(abs(delta_matrix))
    return(values)
  }
}

GetL2Distribution <- function(delta_matrix) {
  # Function that computes the L2 norm over all rows for an input matrix
  # of differences.
  if(sum(delta_matrix) == 0) {
    return("Gradient not monotonic")
  } else {
    values <- rowSums(delta_matrix^2) %>% sqrt()
    return(values)
  }
}


################################################################################
# Analyse results
################################################################################

cutoff_deep_weighted <- cutoff_deep_analysis[paste0("weighted_", cutoff_options)]
cutoff_deep_unweighted <- cutoff_deep_analysis[paste0("unweighted_", cutoff_options)]

# Convert ER distances to weight distribution for all configurations
for (i in 1:length(cutoff_deep_weighted)){
  cutoff_deep_weighted[[i]]$ER_matrices_DIST <- Convert2Relative(cutoff_deep_weighted[[i]]$ER_matrices)
  cutoff_deep_unweighted[[i]]$ER_matrices_DIST <- Convert2Relative(cutoff_deep_unweighted[[i]]$ER_matrices)
  print(round(i/length(cutoff_deep_weighted)*100))
}
remove(cutoff_deep_analysis)
gc()


###########              Visualise kernel for cutoff = 10            ###########

nbhd_size <- 50
factor_increase <- 5

R_Asthma <- cutoff_deep_weighted$weighted_10$ER_matrices_DIST$R_asthma_normalised
R_Control <- cutoff_deep_weighted$weighted_10$ER_matrices_DIST$R_control_normalised

complete_data <- rbind(R_Asthma, R_Control)
nbhd_count <- apply(complete_data, MARGIN = 1, FUN = function(this_row) {
  sort(this_row, decreasing = FALSE)[nbhd_size]
})
nbhd_border <- mean(nbhd_count) %>% print # check 

x <- seq(range(complete_data)[1], range(complete_data)[2], length.out = 1000)
minimum_data <- min(complete_data[complete_data != 0])

complete_data <- c(R_Asthma[upper.tri(R_Asthma)], 
                   R_Control[upper.tri(R_Control)]
)

hist_ER <- hist(complete_data, breaks = 10000, plot = FALSE)
sigma_RBF <- uniroot(function(s){
  prob_density <- RBF(hist_ER$mids, s)
  distribution <- prob_density*hist_ER$counts
  density_percentage <- sum(distribution[hist_ER$mids <= nbhd_border])/sum(distribution)
  populations_percentage <- sum(hist_ER$counts[hist_ER$mids <= nbhd_border])/sum(hist_ER$counts)
  return(density_percentage/populations_percentage - factor_increase)
}, c(.Machine$double.eps^(1/4), 1), tol = .Machine$double.eps)$root

print(paste("Gradient is monotonic:", IsMonotonicSeries(grad(function(x){RBF(x, s=sigma_RBF)}, hist_ER$mids))))
print(paste("Sigma found:", sigma_RBF))

P_RBF_kernel <- RBF(x, s=sigma_RBF)

perc_of_data <- 0.95
x_limits <- c(0, quantile(complete_data, probs = c(perc_of_data)))

complete_data <- c(R_Asthma[upper.tri(R_Asthma)], 
                   R_Control[upper.tri(R_Control)]
)
complete_hist <- hist(complete_data[complete_data < x_limits[2]], breaks = 160, plot = FALSE)

P_RBF_kernel_adj <- P_RBF_kernel*max(complete_hist$counts)/RBF(0, s= sigma_RBF)


# RBF only
plot(x = x, y = P_RBF_kernel_adj, type = "l", col = "red",# xaxt = "n", yaxt = "n",
     xlim = x_limits,
     main=paste0("RBF kernel"),
     xlab="Relative resistance distance",
     ylab="Affinity / Similarity"
)
grad_RBF <- grad(function(x){RBF(x, s=sigma_RBF)}, x)
points(x = x[which.min(grad_RBF)], y = P_RBF_kernel_adj[which.min(grad_RBF)], type = "p", col = "red")
abline(v = min(complete_data), col="blue", lwd=1, lty=2)
abline(v = nbhd_border, col="black", lwd=1, lty=2)


### Perform network comparison for different hyperparameter configurations #####

neighbourhood_size <- 50
factor_increase <- 5

selection_factors <- c(3,5,8)
selection_sizes <- c(5, 10, 15, 25, 35, 50, 75)

length(selection_factors)*length(selection_sizes)

iter <- 1
for(neighbourhood_size in selection_sizes){
  for(factor_increase in selection_factors){
    
    print(paste('Neighbourhood size:', neighbourhood_size))
    print(paste('Factor increase:', factor_increase))
    print(paste("Iteration", iter, "v/d", length(selection_factors)*length(selection_sizes)))
    iter <- iter + 1
    
    
    # Compute delta similarities normalised
    for (i in 1:length(cutoff_deep_weighted)){
      print(paste0(i, ' v/d ', length(cutoff_deep_weighted)))
      if (max(weighted_options) == 1){
        cutoff_deep_weighted[[i]]$delta_similarities <- computeSimilarityRBF(cutoff_deep_weighted[[i]]$ER_matrices_DIST$R_asthma_normalised,
                                                                             cutoff_deep_weighted[[i]]$ER_matrices_DIST$R_control_normalised,
                                                                             factor_increase = factor_increase,
                                                                             nbhd_size = neighbourhood_size,
                                                                             verbose = TRUE
        ) %>% {.$Similarities_asthma - .$Similarities_control}
      }
      
      if (min(weighted_options) == 0){
        cutoff_deep_unweighted[[i]]$delta_similarities <- computeSimilarityRBF(cutoff_deep_unweighted[[i]]$ER_matrices_DIST$R_asthma_normalised,
                                                                               cutoff_deep_unweighted[[i]]$ER_matrices_DIST$R_control_normalised,
                                                                               factor_increase = factor_increase,
                                                                               nbhd_size = neighbourhood_size,
                                                                               verbose = TRUE
        ) %>% {.$Similarities_asthma - .$Similarities_control}
      }
    }
    
    # Analyse  Test Statistics
    outlength <- 10
    # Weighted
    weighted_results_L1 <- foreach(i = 1:length(cutoff_deep_weighted), .combine = "rbind") %do% {
      out <- AnalyseTS(GetL1Distribution(cutoff_deep_weighted[[i]]$delta_similarities),
                       name=paste("L1", names(cutoff_deep_weighted)[i]),
                       plot = FALSE,
                       out_length = 15,
                       roundOutput = TRUE)
      if(out[1] == "Gradient not monotonic"){
        return(c(names(cutoff_deep_weighted)[i], "Gradient not monotonic", rep("", 2*outlength)))
      } else {
        return(c(names(cutoff_deep_weighted)[i],
                 which.min(diff(out)),#diff(out) %>% {(1:length(.))[. == min(.)]},
                 c(rbind(names(out)[1:outlength], out[1:outlength]))))
      }
    } %>% as.data.frame %>% print
    
    # write.csv(weighted_results_L1, paste0(results_directory,"wghtd_L1_nbhd(", 
    # neighbourhood_size,")_dof(", dof,")_frac(",fraction_in_neighbourhood, ").csv"), row.names = FALSE)
    write.csv(weighted_results_L1, paste0(results_directory,"F_L1_", 
                                          neighbourhood_size,"_F",factor_increase, ".csv"), row.names = FALSE)
    
    weighted_results_L2 <- foreach(i = 1:length(cutoff_deep_weighted), .combine = "rbind") %do% {
      out <- AnalyseTS(GetL2Distribution(cutoff_deep_weighted[[i]]$delta_similarities),
                       name=paste("L1", names(cutoff_deep_weighted)[i]),
                       plot = FALSE,
                       out_length = 15)
      if(out[1] == "Gradient not monotonic"){
        return(c(names(cutoff_deep_weighted)[i], "Gradient not monotonic", rep("", 2*outlength)))
      }
      return(c(names(cutoff_deep_weighted)[i],
               which.min(diff(out)),#diff(out) %>% {(1:length(.))[. == min(.)]},
               c(rbind(names(out)[1:outlength], out[1:outlength]))))
    } %>% as.data.frame %>% print
    
    write.csv(weighted_results_L2, paste0(results_directory,"F_L2_", 
                                          neighbourhood_size,"_F",factor_increase, ".csv"), row.names = FALSE)
    # write.csv(weighted_results_L2, paste0(results_directory,"wghtd_L2_nbhd(", 
    # neighbourhood_size,")_dof(", dof,")_frac(",fraction_in_neighbourhood, ").csv"), row.names = FALSE)
  }
}



###########     Output comparison results for cutoff = 10            ###########
cutoff_deep_weighted[["weighted_10"]]$delta_similarities <- computeSimilarityRBF(cutoff_deep_weighted[["weighted_10"]]$ER_matrices_DIST$R_asthma_normalised,
                                                                     cutoff_deep_weighted[["weighted_10"]]$ER_matrices_DIST$R_control_normalised,
                                                                     factor_increase = factor_increase,
                                                                     nbhd_size = neighbourhood_size,
                                                                     verbose = TRUE
) %>% {.$Similarities_asthma - .$Similarities_control}

out <- AnalyseTS(GetL2Distribution(cutoff_deep_weighted[["weighted_10"]]$delta_similarities),
                 name=paste("L2", "weighted_10"),
                 plot = TRUE,
                 out_length = nrow(cutoff_deep_weighted[["weighted_10"]]$delta_similarities),
                 roundOutput = FALSE
)


temp_L2 <- GetL2Distribution(cutoff_deep_weighted[["weighted_10"]]$delta_similarities)
temp_L2_out <- AnalyseTS(temp_L2,
                         name=paste("L2", "weighted_10"),
                         plot = TRUE,
                         out_length = length(temp_L2),
                         roundOutput = TRUE)

temp <- boxplot(temp_L2_out)

# Print first ten results
temp_L2_out[1:10]


# Visualise results in boxplot
identified_genes <- c("NDUFA5", "SNX3", "ALDOA", "FAM3D") 
par(mar = c(1, 5, 2, 2))  
boxplot(temp_L2_out, border = "black", medcol = "darkblue", col = "lightblue1",
        main = "Distances between gene affinity profiles",
        ylab = "Euclidean distance",
        xlab = NULL)
points(x = temp_L2_out["FAM3D"], pch = 21, col = "black", bg = "red", lwd = 1.5)
points(x = temp_L2_out["ALDOA"], pch = 21, col = "black", bg = "red", lwd = 1.5)
points(x = temp_L2_out["SNX3"], pch = 21, col = "black", bg = "red", lwd = 1.5)
points(x = temp_L2_out["NDUFA5"], pch = 21, col = "black", bg = "red", lwd = 1.5)

quantile(temp_L2_out, probs = c(0.25, 0.75)) %>% {.[2]+1.5*(.[2] - .[1])} # Interquartile range * 1.5 tukey citation




