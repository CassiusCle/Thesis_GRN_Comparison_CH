################################################################################
### Name:        script_3.4_visualise_connectivity_differences.R
### Description: Create visualisations of connectivity differences
### Author:      Clemens Pieter Dick Hendrickx
################################################################################

################################################################################
# Load Packages
################################################################################

LoadPackages <- function(){
  if(!require(tidyverse))       {install.packages("tidyverse")};    library(tidyverse)
  if(!require(readr))           {install.packages("readr")};        library(readr)
  if(!require(foreach))         {install.packages("foreach")};      library(foreach)
  if(!require(doSNOW))          {install.packages("doSNOW")};       library(doSNOW)
  if(!require(doParallel))      {install.packages("doParallel")};   library(doParallel)
  if(!require(pracma))          {install.packages("pracma")};       library(pracma)
  if(!require(ggraph))          {install.packages("ggraph")};       library(ggraph)
  if(!require(igraph))          {install.packages("igraph")};       library(igraph)
  if(!require(circlize))        {install.packages("circlize")};     library(circlize)
  if(!require(tidyverse))       {install.packages("tidyverse")};    library(tidyverse)
  if(!require(matrixcalc))      {install.packages("matrixcalc")};   library(matrixcalc)
  if(!require(hash))            {install.packages("hash")};         library(hash)
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
cutOff <- 10

results_directory <- "C:\\Users\\clemenshendrickx\\Documents\\Thesis\\Thesis Data\\Results v1.0\\"

cutoff_options <- seq(from = 10, to = 10, by = 0.5)
weighted_options <- TRUE

################################################################################
# Load data
################################################################################
Weight_Asthma_10_mean <- readRDS(paste0(dir, folder_name, 'Weight_matrix_','Asthma_',cutOff,'_mean.rds'))
Weight_Control_10_mean <- readRDS(paste0(dir, folder_name, 'Weight_matrix_','Control_',cutOff,'_mean.rds'))

cutoff_deep_analysis <- readRDS(paste0(results_directory, 'cutoff_deep_analysis_8_12_2021.rds'))

cutoff_deep_weighted <- cutoff_deep_analysis[paste0("weighted_", cutoff_options)]
remove(cutoff_deep_analysis)


################################################################################
# Auxilliary functions
################################################################################
makeNetworksComparable <- function(Asthma_weights, Control_weights) {
  # Function that makes sure that both networks contain the same genes
  included_genes <- intersect(rownames(Asthma_weights), rownames(Control_weights))
  
  Asthma_weights <- Asthma_weights[included_genes, included_genes]
  Control_weights <- Control_weights[included_genes, included_genes]
  
  list(asthma = Asthma_weights,
       control = Control_weights)
}

################################################################################
# Network comparison
################################################################################

networks <- makeNetworksComparable(Weight_Asthma_10_mean, Weight_Control_10_mean)

Weight_Asthma_Sym10 <- networks$asthma
Weight_Control_Sym10 <- networks$control

Adjacency_Asthma <- Weight_Asthma_Sym10 %>% {1.0*(. > 0)}
Adjacency_Control <- Weight_Control_Sym10 %>% {1.0*(. > 0)}

if( all((rownames(Weight_Asthma_Sym10) == colnames(Weight_Asthma_Sym10))&
        all(rownames(Weight_Control_Sym10) == colnames(Weight_Control_Sym10))&
        all(rownames(Weight_Control_Sym10) == rownames(Weight_Asthma_Sym10)))){
  gene_names <- rownames(Weight_Asthma_Sym10)
}

# Distance matrices
Rel_ER_Asthma_10 <- cutoff_deep_weighted$weighted_10$ER_matrices$asthma %>% {./rowSums(.)}
Rel_ER_Control_10 <- cutoff_deep_weighted$weighted_10$ER_matrices$control %>% {./rowSums(.)}

ER_Asthma_10 <- cutoff_deep_weighted$weighted_10$ER_matrices$asthma
ER_Control_10 <- cutoff_deep_weighted$weighted_10$ER_matrices$control

####### Create circular graphs : regions of candidate genes ####################
# NDUFA5	ALDOA SNX3 FAM3D (82.25)
# gene <- "NDUFA5" # 82
# gene <- "ALDOA" # 85
# gene <- "SNX3" # 84
gene <- "FAM3D" # 78


nbhd_size <- 50

asthma_FirstFifty <- Rel_ER_Asthma_10[gene, ] %>% sort %>% {.[1:(nbhd_size+1)]} %>% names %>% print
asthma_FirstFifty_dist <- Rel_ER_Asthma_10[gene, ] %>% sort %>% {.[1:(nbhd_size+1)]} %>% print

control_FirstFifty <- Rel_ER_Control_10[gene, ] %>% sort %>% {.[1:(nbhd_size+1)]} %>% names %>% print
control_FirstFifty_dist <- Rel_ER_Control_10[gene, ] %>% sort %>% {.[1:(nbhd_size+1)]}

FirstFifty_in_common <- intersect(asthma_FirstFifty, control_FirstFifty) %>% print
FirstFifty_asthma_only <- setdiff(asthma_FirstFifty, control_FirstFifty) %>% print
FirstFifty_control_only <- setdiff(control_FirstFifty, asthma_FirstFifty) %>% print

FirstFifty_union <- union(asthma_FirstFifty, control_FirstFifty) %>% print
FirstFifty_union %>% length

Asthma_dist_matrix_FF_union <- ER_Asthma_10[FirstFifty_union, FirstFifty_union]
Control_dist_matrix_FF_union <- ER_Control_10[FirstFifty_union, FirstFifty_union]
Delta_dist_matrix_FF_union <- Asthma_dist_matrix_FF_union - Control_dist_matrix_FF_union # Positive means distance increased

groups <- c("Asthma", "Control", "Both")
group_lengths <- c(length(FirstFifty_asthma_only), length(FirstFifty_control_only), length(FirstFifty_in_common))
group_members <- c(FirstFifty_asthma_only, FirstFifty_control_only, FirstFifty_in_common)

d1 <- data.frame(from="origin", to= groups)
d2 <- data.frame(from=rep(d1$to, times = group_lengths), to = group_members)
hierarchy <- rbind(d1, d2)

all((d2 %>% filter(from == "Asthma") %>% select(to) %>% unlist) == FirstFifty_asthma_only)
all((d2 %>% filter(from == "Control") %>% select(to) %>% unlist) == FirstFifty_control_only)
all((d2 %>% filter(from == "Both") %>% select(to) %>% unlist) == FirstFifty_in_common)


# Edge table
percentile_change <- 0.5
Delta_dist_edges <-  get.data.frame(graph.adjacency(Delta_dist_matrix_FF_union, weighted=TRUE, mode = "undirected"))
Delta_dist_edges$value <- Delta_dist_edges$weight

h <- hash() 
h[d2$to] <- d2$from
Delta_dist_edges$from_group <- values(h, keys=Delta_dist_edges$from)
Delta_dist_edges$to_group <- values(h, keys=Delta_dist_edges$to)

Delta_dist_edges <- Delta_dist_edges %>% filter(abs(value) > quantile(as.numeric(abs(Delta_dist_matrix_FF_union)), probs = percentile_change))



connect <- Delta_dist_edges
col <- connect$value

vertices  <-  data.frame(
  name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to))) 
) 
vertices$group  <-  hierarchy$from[ match( vertices$name, hierarchy$to ) ]


mygraph <- graph_from_data_frame( hierarchy, vertices=vertices )

from  <-  match( connect$from, vertices$name)
to  <-  match( connect$to, vertices$name)

p_graphic <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  # Edges between nodes
  geom_conn_bundle(data = get_con(from = from, to = to, col = col), 
                   tension = 0.78
                   , aes(alpha = abs(col), colour=col)
                   #, show.legend = FALSE
  ) +
  scale_edge_color_gradient2(low = "red", mid = "lightgrey", high = "blue", name = "Change in \nER distance", labels = c("", "Decrease", "No change", "Increase", ""))+#, labels = NULL) +
  scale_edge_alpha(labels = NULL, guide = "none") +
  
  # Node points
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, fill = group, shape = group), colour = "black", stroke = 1, size = 4)+
  scale_shape_manual(values= c(19, 10, 1), name = "Subject group")+
  scale_fill_manual(values= c("red", "darkgrey", "blue"), name = "Subject group") +
  
  # Gene labels
  geom_node_text(aes(x = x*1.1, y=y*1.1, filter = leaf
                     , label=name 
                     , angle = ifelse(abs(node_angle(x, y) - 180) < 90, node_angle(x, y)+180, node_angle(x, y) )
                     , hjust = ifelse(abs(node_angle(x, y) - 180) < 90, 1, 0 )
                     ),
                 size = 3,
                 alpha=1
                 ) +
  # Theme stuff
  theme(panel.background = element_rect(fill='transparent')) +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) +
  labs(title = "Connectivity difference between subject groups",
       subtitle = paste0("Union of neighbourhoods containing the ", nbhd_size, " most strongly connected genes to ", gene))
  
p_graphic

ggsave(paste0(results_directory, "circular_network_", gene,"_", nbhd_size, ".png"), width = 9.26, height = 8.51, units = "in")


####### Create circular graph : Overlapping set of genes    ####################
genes <- c("NDUFA5",	"ALDOA", "SNX3", "FAM3D")
asthma_total <- list()
control_total <- list()
for (i in genes){
  asthma_total[[i]] <- ER_Asthma_10[i, ] %>% sort %>% {.[1:(nbhd_size+1)]} %>% names 
  control_total[[i]] <- ER_Control_10[i, ] %>% sort %>% {.[1:(nbhd_size+1)]} %>% names
}

consistently_asthma_full <- intersect(asthma_total[[1]], asthma_total[[2]]) %>% intersect(., asthma_total[[3]]) %>% intersect(., asthma_total[[4]]) %>% sort
consistently_control_full <- intersect(control_total[[1]], control_total[[2]]) %>% intersect(., control_total[[3]]) %>% intersect(., control_total[[4]]) %>% sort

consistently_asthma_full
consistently_control_full

consistently_both <- intersect(consistently_asthma_full, consistently_control_full) %>% sort
consistently_asthma <- setdiff(consistently_asthma_full, consistently_both)
consistently_control <- setdiff(consistently_control_full, consistently_both)

consistently_union <- union(union(consistently_asthma, consistently_control), consistently_both)

Asthma_dist_matrix_FF_union <- ER_Asthma_10[consistently_union, consistently_union]
Control_dist_matrix_FF_union <- ER_Control_10[consistently_union, consistently_union]
Delta_dist_matrix_FF_union <- Asthma_dist_matrix_FF_union - Control_dist_matrix_FF_union # Positive means distance increased

groups <- c("Asthma", "Control", "Both")
group_lengths <- c(length(consistently_asthma), length(consistently_control), length(consistently_both))
group_members <- c(consistently_asthma, consistently_control, consistently_both)

d1 <- data.frame(from="origin", to= groups)
d2 <- data.frame(from=rep(d1$to, times = group_lengths), to = group_members)
hierarchy <- rbind(d1, d2)

all((d2 %>% filter(from == "Asthma") %>% select(to) %>% unlist) == consistently_asthma)
all((d2 %>% filter(from == "Control") %>% select(to) %>% unlist) == consistently_control)
all((d2 %>% filter(from == "Both") %>% select(to) %>% unlist) == consistently_both)


# Edge table
percentile_change <- 0
Delta_dist_edges <-  get.data.frame(graph.adjacency(Delta_dist_matrix_FF_union, weighted=TRUE, mode = "undirected"))
Delta_dist_edges$value <- Delta_dist_edges$weight

h <- hash() 
h[d2$to] <- d2$from
Delta_dist_edges$from_group <- values(h, keys=Delta_dist_edges$from)
Delta_dist_edges$to_group <- values(h, keys=Delta_dist_edges$to)

Delta_dist_edges <- Delta_dist_edges %>% filter(abs(value) > quantile(as.numeric(abs(Delta_dist_matrix_FF_union)), probs = percentile_change))

connect <- Delta_dist_edges
col <- connect$value

vertices  <-  data.frame(
  name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to))) 
) 

vertices$group  <-  hierarchy$from[ match( vertices$name, hierarchy$to ) ]

# Create a graph object
mygraph <- graph_from_data_frame( hierarchy, vertices=vertices )

from  <-  match( connect$from, vertices$name)
to  <-  match( connect$to, vertices$name)

p_graphic <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  # Edges between nodes
  geom_conn_bundle(data = get_con(from = from, to = to, col = col), 
                   tension = 0.78
                   , aes(alpha = abs(col), colour=col)
                   , n = 200
  ) +
  scale_edge_color_gradient2(low = "red", mid ="white", high = "blue", name = "Change in \nER distance", labels = c("", "Decrease", "No change", "Increase", ""))+#, labels = NULL) +
  scale_edge_alpha(labels = NULL, guide = "none") +
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, fill = group, shape = group), colour = "black", stroke = 1, size = 4)+
  scale_shape_manual(values= c(19, 10, 1), name = "Subject group")+
  scale_fill_manual(values= c("red", "darkgrey", "blue"), name = "Subject group") +
  # Gene labels
  geom_node_text(aes(x = x*1.1, y=y*1.1, filter = leaf
                     , label=name 
                     , angle = ifelse(abs(node_angle(x, y) - 180) < 90, node_angle(x, y)+180, node_angle(x, y) )
                     , hjust = ifelse(abs(node_angle(x, y) - 180) < 90, 1, 0 )
  ),
  size = 4,
  alpha=1
  ) +
  theme(panel.background = element_rect(fill='transparent')) +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) +
  labs(title = "Connectivity difference between subject groups",
       subtitle = paste0("Genes consistently found in the neighbourhoods of NDUFA5, ALDOA, SNX3 and FAM3D")
       )


p_graphic

ggsave(paste0(results_directory, "circular_network_union_consistent","_", nbhd_size, ".png"), width = 9.26, height = 8.51, units = "in")


