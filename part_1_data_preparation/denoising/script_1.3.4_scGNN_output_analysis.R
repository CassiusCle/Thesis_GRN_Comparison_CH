################################################################################
### Name:        script_1.3.4_scGNN_output_analysis.R
### Description: Analyses output of scGNN
### Author:      Clemens Pieter Dick Hendrickx
################################################################################

################################################################################
# Load Packages
################################################################################

LoadPackages <- function(){
  if(!require(tidyverse))       {install.packages("tidyverse")};    library(tidyverse)
  if(!require(foreach))         {install.packages("foreach")};      library(foreach)
  
  if(!require(moments))         {install.packages("moments")};      library(moments)
  if(!require(tseries))         {install.packages("tseries")};      library(tseries)
  if(!require(pdist))           {install.packages("pdist")};        library(pdist)
  if(!require(philentropy))     {install.packages("philentropy")};  library(philentropy)
  if(!require(Rcpp))            {install.packages("Rcpp")};         library(Rcpp)
  if(!require(lavaan))          {install.packages("lavaan")};       library(lavaan)
  }
LoadPackages()

gc()

# Use expression inspectation
imputation_pre_data <- read_csv("~/Thesis/Thesis Data/Selected Data/Intermediate data/QC_OUTPUT_minG.CA1000_CB20_GA10_GB5/scGNN_pre_output/ng_2000_log_excl/Use_expression.csv"
                                ) %>% as.data.frame
colnames(imputation_pre_data)[1] <- 'Gene'
rownames(imputation_pre_data) <- imputation_pre_data[,'Gene']
imputation_pre_data[,'Gene'] <- NULL

# rownames(imputation_pre_data) <- imputation_pre_data[,'X1']
# imputation_pre_data[,'X1'] <- NULL

# Verified that Use_expression (for 2000) is just a subset of Use_expression (for 14458) (done in v1.3)

# Transformation is just log1p, normalised counts and use_expression is the same
# Input is the same, only very small difference, but log1p() transformed, confirmed. (done in v1.3)

# 2000 output, 10 clusters excl
ng_2000_log_results <- read_csv("~/Thesis/Thesis Data/Selected Data/Intermediate data/QC_OUTPUT_minG.CA1000_CB20_GA10_GB5/scGNN_output/ng_2000_log_clusters_10_excl/ng_2000_log_excl_results.txt") %>% as.data.frame()
colnames(ng_2000_log_results)[1] <- 'Cell'
rownames(ng_2000_log_results) <- ng_2000_log_results[,'Cell']

metadata_normalised <- read_csv("~/Thesis/Thesis Data/Selected Data/Intermediate data/QC_OUTPUT_minG.CA1000_CB20_GA10_GB5/metadata_normalised_counts_sconeCA1000_CB20_GA10_GB5.csv") %>% as.data.frame()
# rownames(metadata_normalised) <- metadata_normalised[,'X1']
colnames(metadata_normalised)[1] <- 'Cell'
rownames(metadata_normalised) <- metadata_normalised[,'Cell']
#metadata_normalised[,'X1'] <- NULL
old_clustering <- metadata_normalised %>% select(cluster, healthStatus, Cell)

clusters_2000 <- merge(x = ng_2000_log_results, y = old_clustering, by = 'Cell', all.x = TRUE, all.y = FALSE)

compare_clustering <- function(clusters) {
  old_types <- clusters %>% select(cluster) %>% table %>% names
  k_range <- 0:(length(table(select(clusters, Celltype))) - 1)
  
  clustering_compare <- data.frame(NULL)
  
  for (k in k_range){
    counts <- rep(0,length(old_types))
    names(counts) <- old_types
    
    this_counts <- clusters %>% filter(Celltype == k) %>% select(cluster) %>% table
    counts[names(this_counts)] <- this_counts 
    clustering_compare <- rbind(clustering_compare, counts)
  }
  colnames(clustering_compare) <- old_types
  rownames(clustering_compare) <- k_range
  
  return(clustering_compare)
}
compare_clustering(clusters_2000)

compare_healthStatus <- function(clusters) {
  old_types <- clusters %>% select(healthStatus) %>% table %>% names
  k_range <- 0:(length(table(select(clusters, Celltype))) - 1)
  
  clustering_compare <- data.frame(NULL)
  
  for (k in k_range){
    counts <- rep(0,length(old_types))
    names(counts) <- old_types
    
    this_counts <- clusters %>% filter(Celltype == k) %>% select(healthStatus) %>% table
    counts[names(this_counts)] <- this_counts 
    clustering_compare <- rbind(clustering_compare, counts)
  }
  colnames(clustering_compare) <- old_types
  rownames(clustering_compare) <- k_range
  
  return(clustering_compare)
}
compare_healthStatus(clusters_2000)

(clusters_2000$healthStatus %>% table)*100/nrow(clusters_2000)

imputation_post_data <- read_csv("~/Thesis/Thesis Data/Selected Data/Intermediate data/QC_OUTPUT_minG.CA1000_CB20_GA10_GB5/scGNN_output/ng_2000_log_clusters_10_excl/ng_2000_log_excl_recon.csv") %>% as.data.frame()
colnames(imputation_post_data)[1] <- 'Gene'
rownames(imputation_post_data) <- imputation_post_data[,'Gene']
imputation_post_data[,'Gene'] <- NULL
imputation_post_data <- imputation_post_data[rownames(imputation_pre_data), colnames(imputation_pre_data)]

selection_basal_secr <- clusters_2000 %>% filter(Celltype %in% c(0, 2, 3, 6, 7, 9)) %>% select(Cell) %>% unlist

imp_pre_sel_data <- imputation_pre_data[,selection_basal_secr]
imp_post_sel_data <- imputation_post_data[,selection_basal_secr]



get_statistics <- function(df){
  result <- foreach(gene = rownames(df), .combine = rbind ) %do% {
    this_row <- df[gene,] %>% unlist
    c(name = gene,
      sum = sum(this_row),
      mean = mean(this_row),
      var = var(this_row),
      frac_zero = sum(this_row == 0)/length(this_row)
      )
  } %>% as.data.frame
  result$sum <- as.numeric(result$sum)
  result$mean <- as.numeric(result$mean)
  result$var <- as.numeric(result$var)
  result$frac_zero <- as.numeric(result$frac_zero)
  return(result)
}

# Check lowest variance and stuff
metadata_2000 <- merge(x = ng_2000_log_results, y = metadata_normalised, by = 'Cell', all.x = TRUE, all.y = FALSE)
ids_asthma <- metadata_2000 %>% filter(healthStatus == 'PersA') %>% select(seqID) %>% unlist
ids_control <- metadata_2000 %>% filter(healthStatus == 'Ctrl') %>% select(seqID) %>% unlist

quickStats <- imp_post_sel_data %>% get_statistics
quickStats_ast <- imp_post_sel_data[, colnames(imp_post_sel_data) %in% ids_asthma] %>% get_statistics
quickStats_contr <- imp_post_sel_data[, colnames(imp_post_sel_data) %in% ids_control] %>% get_statistics

quickStats_ast %>% filter(var == 0) %>% nrow
quickStats_contr %>% filter(var == 0) %>% nrow



pre_data_choice <- imputation_pre_data#imp_pre_sel_data
post_data_choice <- imputation_post_data#imp_post_sel_data

ids_asthma <- metadata_normalised %>% filter(healthStatus == 'PersA') %>% select(seqID) %>% unlist
ids_asthma <- ids_asthma[ids_asthma %in% colnames(pre_data_choice)]
ids_control <- metadata_normalised %>% filter(healthStatus == 'Ctrl') %>% select(seqID) %>% unlist
ids_control <- ids_control[ids_control %in% colnames(post_data_choice)]


imputation_pre_asthma <- pre_data_choice[,ids_asthma]
imputation_pre_control <- pre_data_choice[,ids_control]

imputation_post_asthma <- post_data_choice[,ids_asthma]
imputation_post_control <- post_data_choice[,ids_control]

##### Stats

describe_data <- function(data) {
  data <- data %>% as.matrix 
  par(mfrow=c(1,2))
  data %>% hist
  data %>% density %>% plot
  print(summary(data))
  print(paste('Variance:',var(data)))
  print(paste('Skewness:',skewness(data)))
  print(paste('Kurtosis:',kurtosis(data)))
  # print(jarque.bera.test(data))
  # dev.off() 
}

describe_data(as.vector(as.matrix(imputation_pre_data)))
describe_data(as.vector(as.matrix(imputation_post_data)))

par(mfrow=c(1,2))
hist(as.vector(as.matrix(imputation_pre_data)), main = 'Pre-imputation',
     xlab = 'Log read counts')
hist(as.vector(as.matrix(imputation_post_data)), main = 'Post-imputation',
     xlab = 'Log read counts')

describe_data(as.vector(as.matrix(imputation_pre_asthma)))
describe_data(as.vector(as.matrix(imputation_post_asthma)))

describe_data(as.vector(as.matrix(imputation_pre_control)))
describe_data(as.vector(as.matrix(imputation_post_control)))
#####


#### More than 10% change check

check_change <- function(pre, post, percentage = 10, dim = 'cells') {
  if (dim == 'cells') 
    df <- as.data.frame(list(pre = colMeans(pre), post = colMeans(post))) 
  else 
    df <- as.data.frame(list(pre = rowMeans(pre), post = rowMeans(post))) 
  df$delta_perc <- (df$post-df$pre)/df$pre*100
  
  if (dim == 'genes') {
    to_zero <- df %>% filter(post == 0)
    print(paste('Out of', nrow(df), 'total genes,', nrow(to_zero), 'genes put to zero'))
    print(paste('Average pre expression zero\'s:', mean(to_zero$pre)))
    df <- (gene_averages %>% filter(abs(delta_perc) != 100))
  }
  
  print(paste('Average expression change:', mean(df$delta_perc)))
  print(paste('Number of',dim, 'more than',percentage, 'change:', nrow(filter(df, abs(delta_perc) > 10))))
  
  if(dim == 'genes')
    return(rownames(filter(df, abs(delta_perc) > 10)))
}

df <- as.data.frame(list(pre = rowMeans(imputation_pre_data), post = rowMeans(imputation_post_data))) 
df$delta_perc <- (df$post-df$pre)/df$pre*100
df$delta_perc %>% abs %>% max
df$delta_perc %>% hist

rownames(filter(df, delta_perc == -100))

imputation_pre_stats$rank_zero <- rank(imputation_pre_stats$frac_zero)

temp <- imputation_pre_stats %>% filter(name %in% rownames(filter(df, delta_perc == -100)))
temp <- imputation_pre_stats %>% filter(name %in% rownames(filter(df, post < 0.0001)))

temp <- filter(df, post < 0.0006 )

temp$frac_zero %>% mean
temp2 <- imputation_pre_stats %>% filter( !(name %in% rownames(filter(df, delta_perc == -100))))


filter(df, delta_perc != -100 & abs(delta_perc) > 10 & abs(delta_perc) < 90)$delta_perc %>% abs %>% mean


imputation_pre_stats$frac_zero %>% mean
imputation_pre_stats$mean %>% mean

check_change(pre_data_choice, post_data_choice, dim = 'cells')
check_change(imputation_pre_asthma, imputation_post_asthma, dim = 'cells')
check_change(imputation_pre_control, imputation_post_control, dim = 'cells')

check_change(pre_data_choice, post_data_choice, dim = 'genes')
changed_asthma <- check_change(imputation_pre_asthma, imputation_post_asthma, dim = 'genes')
changed_control <- check_change(imputation_pre_control, imputation_post_control, dim = 'genes')

all(changed_asthma %in% changed_control)
all(changed_control %in% changed_asthma)

