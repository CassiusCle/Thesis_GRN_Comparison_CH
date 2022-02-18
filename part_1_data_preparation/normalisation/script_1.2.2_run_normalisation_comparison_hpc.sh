#!/bin/bash
#SBATCH --job-name=scone_norm
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --output=scone_norm_hpc_v2.0_all_options.out
#SBATCH --mem=80gb
#SBATCH --partition= regular


directory='/groups/umcg-micompany/tmp01/scRNA_Networks_CPDH/'

cd $directory
pwd

# Modules
module purge
module load R/4.0.3-foss-2018b-bare  
module load RPlus/4.0.3-foss-2018b-v21.05.1
module list

# Set parameters for Rscript
#data_dir="${directory}data/intermediate_data/"
data_dir="data/intermediate_data/"
batch_choice="run"
threads=22

# Run Rscript
filter_used="CA1000_CB20_GA10_GB5"
Rscript code/preprocessing/scone_normalisation_hpc_v2.0.R $data_dir $filter_used $batch_choice $threads >& ${directory}code/preprocessing/scone_norm_R_console_v2.0_opt_${filter_used}.out

# Run Rscript
filter_used="CA1000_CB15_GA10_GB5"
Rscript code/preprocessing/scone_normalisation_hpc_v2.0.R $data_dir $filter_used $batch_choice $threads >& ${directory}code/preprocessing/scone_norm_R_console_v2.0_opt_${filter_used}.out

# Run Rscript
filter_used="CA1000_CB10_GA10_GB5"
Rscript code/preprocessing/scone_normalisation_hpc_v2.0.R $data_dir $filter_used $batch_choice $threads >& ${directory}code/preprocessing/scone_norm_R_console_v2.0_opt_${filter_used}.out

# Run Rscript
filter_used="CA1000_CB5_GA10_GB5"
Rscript code/preprocessing/scone_normalisation_hpc_v2.0.R $data_dir $filter_used $batch_choice $threads >& ${directory}code/preprocessing/scone_norm_R_console_v2.0_opt_${filter_used}.out

