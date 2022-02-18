#!/bin/bash
#SBATCH --job-name=scGNN_job_mC
#SBATCH --time=04:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=/groups/umcg-micompany/tmp01/scRNA_Networks_CPDH/code/imputation/scGNN_hpc_2000_7.0_excl.out
#SBATCH --mem=10gb

# This script does the actual scGNN

# Hyperparameters from preprocessing
numGenes=2000
transform=log

# Define directory on Gearshift cluster
directory=/groups/umcg-micompany/tmp01/scRNA_Networks_CPDH/ 
scGNN_dir=${directory}environments/scGNN/

datasetDir=${directory}data/intermediate_data/QC_OUTPUT_minG.CA1000_CB20_GA10_GB5/scGNN_pre_output/
datasetName=ng_${numGenes}_${transform}_excl

#data_dir_high=${directory}data/intermediate_data/QC_OUTPUT_minG.CA1000_CB20_GA10_GB5/
#data_input_dir=${data_dir_high}scGNN_pre_output/ng_${numGenes}_${transform}/
#datasetName=''


# Parameters for the script
# 	Required parameters:
		#datasetName
			# defines the folder of scRNA-Seq
		LTMGDir=$datasetDir$datasetName
			# defines folder of the preprocessed LTMG output
		
		# Output dir moved down!
		#outputDir=scGNN_output/ng_${numGenes}_${transform}_clusters_${n_clusters}/
			#${data_dir_high}scGNN_output/ng_${numGenes}_${transform}/
			# output folder of the results

# 	Clustering related
		clustering_method=KMeans
		    #LouvainK 	
			# clustering method on identifying celltypes from the embedding. Default is LouvainK
			#Supporting clustering type: Louvain/KMeans/SpectralClustering/AffinityPropagation/AgglomerativeClustering/AgglomerativeClusteringK/Birch/BirchN/MeanShift/OPTICS/LouvainK/LouvainB
		n_clusters=10				
			# predefines the number of clusters, it only used for clustering methods need a number of clusters input as KNN
		maxClusterNumber=30
			# defines the maximum number of cluster allowed, default is 30. This parameter prevents extreme cases that too many clusters identified by Louvian clustering
		minMemberinCluster=5
			# defines the minimum number of cells in a cluster, default is 5. This parameter prevents extreme cases that too many clusters identified by Louvain clustering.
		resolution=auto
			# controls the number of clusters identified by Louvain clustering. This parameter can be set between 0.4 and 1.2 in most cases. According to results on benchmarks, we set default 'auto'.


outputDir=scGNN_output/ng_${numGenes}_${transform}_clusters_${n_clusters}_excl/


#	Optional: Hyperparameters
		EM_iteration=10
			# defines the number of iteration, default is 10
		Regu_epochs=500
			# defines epochs in Feature Autoencoder initially, default is 500
		EM_epochs=200
			# defines epochs in Feature Autoencoder in the iteration, default is 200
		cluster_epochs=200
			# defines epochs in the Cluster Autoencoder, default is 200
		#k=	
			# is k of the K-Nearest-Neighour Graph
		#knn-distance=euclidean
			# distance type of building K-Nearest-Neighour Graph, supported type: euclidean/cosine/correlation (default: euclidean)
		#GAEepochs=
			# number of epochs to train in Graph Autoencoder

#	Optional: Performance
		quickmode=							
			# whether or not to bypass the Cluster Autoencoder. Opt: (blanc) or --quickmode
		useGAEembedding=--useGAEembedding	
			# whether use Graph Autoencoder
		regulized_type=LTMG					
			# is the regularized type: noregu/LTMG, default is to use LTMG
		#alphaRegularizePara= 				
			# alpha in the manuscript, the intensity of the regularizer
		EMregulized_type=Celltype			
			# defines the imputation regularizer type:noregu/Graph/Celltype, default: Celltype
		#gammaImputePara= 					
			# defines the intensity of LTMG regularizer in Imputation
		#graphImputePara= 					
			# defines the intensity of graph regularizer in Imputation
		#celltypeImputePara=				
			# defines the intensity of celltype regularizer in Imputation
		L1Para=1.0 							
			# defines the intensity of L1 regularizer, default: 1.0
		L2Para=0.0 							
			# defines the intensity of L2 regularizer, defualt: 0.0
		#saveinternal= 						
			# whether output internal results for debug usage

#	Optional: Speed
		no_cuda=--no-cuda	
			# defines devices in usage. Default is using GPU, add --no-cuda in command line if you only have CPU.
		coresUsage=12		
			# defines how many cores can be used. default: 1. Change this value if you want to use more.


#cd $scGNN_dir
cd $datasetDir
cd ..


# Modules
module purge
module load Anaconda3/5.3.0 
module list

source activate scgnnEnv

# Full command
python3 -W ignore ${scGNN_dir}scGNN.py \
--datasetName $datasetName \
--datasetDir $datasetDir \
--LTMGDir $LTMGDir \
--outputDir $outputDir \
--nonsparseMode \
--n-clusters $n_clusters \
--clustering-method $clustering_method \
--coresUsage $coresUsage \
$no_cuda


# Full command
#python3 -W ignore ${scGNN_dir}scGNN.py \
#--datasetName $datasetName \
#--datasetDir $datasetDir \
#--LTMGDir $LTMGDir \
#--outputDir $outputDir \
#--nonsparseMode \
#--coresUsage $coresUsage \
#$no_cuda \
#--clustering-method $clustering_method \
#--maxClusterNumber $maxClusterNumber \
#--minMemberinCluster $minMemberinCluster \
#--resolution $resolution \
#--EM-iteration $EM_iteration \ 
#--Regu-epochs $Regu_epochs \
#--EM-epochs $EM_epochs \
#--cluster-epochs $cluster_epochs \
#$quickmode \
#$useGAEembedding \ 
#--regulized-type $regulized_type \
#--EMregulized-type $EMregulized_type
#--L1Para $L1Para \
#--L2Para $L2Para 


