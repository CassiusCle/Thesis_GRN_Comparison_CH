#!/bin/bash
#SBATCH --job-name=network_inference_job
#SBATCH --time=04:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=21
#SBATCH --output=/groups/umcg-micompany/tmp01/scRNA_Networks_CPDH/code/network_inference/infer_networks_hpc_imputed_v3.0.out
#SBATCH --mem=10gb

version=3.0

# Arguments
directory=/groups/umcg-micompany/tmp01/scRNA_Networks_CPDH/data/intermediate_data/QC_OUTPUT_minG.CA1000_CB20_GA10_GB5/ 
num_genes=2000 # Either 2000 or 4000
excluding=_excl
n_clusters=_clusters_10
imputeStatus=original
cores_used=21
seedA=1913
seedC=1913

# folder_name=scGNN_output/ng_${num_genes}_log${n_clusters}${excluding}/

# Modules
module purge
module load Anaconda3/5.3.0 
module list

# Setup
#conda config --add envs_dirs /groups/umcg-micompany/tmp01/scRNA_Networks_CPDH/environments/ # Setup
#conda config --add pkgs_dirs /groups/umcg-micompany/tmp01/scRNA_Networks_CPDH/environments/pkgs # Setup
#conda create -n arboreto-env python=3.9.5 pip # Setup

source activate arboreto-env

#conda install -c bioconda arboreto

# Full command
#python -W ignore network_inference_v${version}_hpc.py \
#$directory \
#$num_genes \
#$excluding \
#$n_clusters \
#$imputeStatus \
#$cores_used \
#$seedA \ 
#$seedC \
#>& python_out_network_inference_hpc_v${version}.out

python -W ignore -m sklearnex network_inference_v${version}_hpc.py \
$directory \
$num_genes \
$excluding \
$n_clusters \
$imputeStatus \
$cores_used \
$seedA \ 
$seedC \
>& python_out_network_inference_hpc_v${version}.out
