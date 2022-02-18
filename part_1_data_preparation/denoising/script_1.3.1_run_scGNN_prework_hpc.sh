#!/bin/bash
#SBATCH --job-name=scGNN_preprocessing_job
#SBATCH --time=06:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/groups/umcg-micompany/tmp01/scRNA_Networks_CPDH/code/imputation/scGNN_pre_2000_log_v4.0_used.out
#SBATCH --mem=25gb

# This script does the preprocessing for scGNN

# scGNN Hyperparameters
numGenes=2000
transform=log

# Define directory on Gearshift cluster
directory=/groups/umcg-micompany/tmp01/scRNA_Networks_CPDH/ 
scGNN_dir=${directory}environments/scGNN/
data_dir=${directory}data/intermediate_data/QC_OUTPUT_minG.CA1000_CB20_GA10_GB5/
data_name=normalised_counts_sconeCA1000_CB20_GA10_GB5_excluding.csv

cd $scGNN_dir

# Modules
module purge

module load Anaconda3/5.3.0 
module list

source activate scgnnEnv

python3 -W ignore ${scGNN_dir}PreprocessingscGNN.py \
--datasetName $data_name \
--datasetDir $data_dir \
--LTMGDir $data_dir \
--filetype CSV \
--transform $transform \
--geneSelectnum $numGenes \
--inferLTMGTag

conda deactivate

this_data_dir=${data_dir}scGNN_pre_output/ng_${numGenes}_${transform}_excl/

mv ${data_dir}Use_expression.csv ${this_data_dir}Use_expression.csv
mv ${data_dir}LTMG_sparse.mtx ${this_data_dir}LTMG_sparse.mtx
