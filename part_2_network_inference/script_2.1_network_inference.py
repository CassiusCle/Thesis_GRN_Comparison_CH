#!/usr/bin/env python3
# coding: utf-8

'''
    File name: network_inference_v1.0_hpc.py
    Author: Clemens Pieter Dick Hendrickx
    Python Version: 3.9.5
'''

# Load packages
import pandas as pd
import sys
from distributed import LocalCluster, Client

from arboreto.algo import grnboost2 #, genie3
#from arboreto.utils import load_tf_names

# Run on HPC cluster?
hpc_cluster = True

if (hpc_cluster): # Run cluster
    argument_list = sys.argv
    print('argument_list:', argument_list)
    directory = argument_list[1]
    print('directory:', directory)
    num_genes = argument_list[2]
    print('num_genes:', num_genes)
    excluding = argument_list[3]
    print('excluding:', excluding)
    n_clusters = argument_list[4]
    print('n_clusters:', n_clusters)
    imputeStatus = argument_list[5]
    print('imputeStatus:', imputeStatus)
    cores_used = argument_list[6]
    print('cores_used:', cores_used)
    seedA = argument_list[7]
    print('seedA:', seedA)
    seedC = argument_list[8]
    print('seedC:', seedC)
else: # Run local
    directory = "C:/Users/clemenshendrickx/Documents/Thesis/Thesis Data/Selected Data/Intermediate data/QC_OUTPUT_minG.CA1000_CB20_GA10_GB5/"
    num_genes = 2000
    excluding = '_excl'
    n_clusters = '_clusters_10'
    imputeStatus = 'imputed'
    cores_used = 4

folder_name = 'scGNN_output/ng_'+ str(num_genes) + '_log'+ n_clusters + excluding + '/'
data_dir = directory + folder_name

path_asthma_imputed  = data_dir + 'expression_asthma_'+ imputeStatus + '.csv'
path_control_imputed = data_dir + 'expression_control_'+ imputeStatus + '.csv'

# Load input data
# N.B.: arboreto requires expression matrix to have cells as rows and genes as columns --> Need to transpose data
df_A_imp = pd.read_csv(path_asthma_imputed)
df_C_imp = pd.read_csv(path_control_imputed)

df_A_imp = df_A_imp.rename(columns={'Unnamed: 0': 'Gene'})
df_C_imp = df_C_imp.rename(columns={'Unnamed: 0': 'Gene'})

df_A_imp = df_A_imp.set_index('Gene').transpose()
df_C_imp = df_C_imp.set_index('Gene').transpose()

# Check dimensions
print('Input Asthma data dimensions:', df_A_imp.shape)
#df_A_imp.head()

print('Input Control data dimensions:', df_C_imp.shape)
#df_C_imp.head()

# Setup parallel computing
if (hpc_cluster):
    local_cluster = LocalCluster(n_workers= 1,#int(cores_used), 
                             threads_per_worker= int(cores_used), #1,
                             processes = False

                             )
else:
    local_cluster = LocalCluster()

custom_client = Client(local_cluster)
print(custom_client)


# Infer networks
network_A = grnboost2(expression_data= df_A_imp,
                     gene_names = list(df_A_imp),
                     seed = int(seedA),
                     verbose = True,
                     client_or_address=custom_client
                     )


network_C = grnboost2(expression_data= df_C_imp,
                     gene_names = list(df_C_imp),
                     seed = int(seedC),
                     verbose = True,
                     client_or_address=custom_client
                     )

print('network_A - Head')
print(network_A.head())

print('network_C - Head')
print(network_C.head())

print('Length Asthma network:', len(network_A))
print('Length Control network:', len(network_C))

# Save networks to CSV
network_A.to_csv(data_dir + 'network_asthma_imputed_' + seedA + '.csv', header=True, index=True)
network_C.to_csv(data_dir + 'network_control_imputed_' + seedC + '.csv', header=True, index=True)

# Close local cluster
custom_client.close()
local_cluster.close()
