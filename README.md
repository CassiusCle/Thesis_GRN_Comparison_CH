# Thesis_GRN_Comparison_CH
The public repository of the code used for the master's thesis in Econometrics and Management Science by Clemens Hendrickx. 
The repository follows the same structure as this thesis, with folders for the different parts and steps of the research. 

Scripts whose name ends in "hpc" are written to be run on the HPC cluster.

## part\_1\_data\_preparation

### filtering
#### script\_1.1.1\_raw\_data\_processing.R
R script that imports the raw data sources and merges them into the right format. The output is one table for the gene expression data and one for the metadata.

#### script\_1.1.2\_data\_filtering.R
R script that implements the data filtering steps described in Section \ref{sec:data}. The output is the filtered expression data set.

### normalisation
#### script\_1.2.1\_normalisation\_hpc.R
R script that implements the scone workflow from \cite{cole2019performance} on the HPC cluster. This workflow implements a variety of normalisation configurations and scores these configurations using different metrics. This comparison is described in detail in Appendix \ref{app:intermediate_results:normalisation}. It outputs a data object containing both the scores from the comparison as the different normalised data sets.

#### script\_1.2.2\_run\_normalisation\_comparison\_hpc.sh
Shell script used for running script 1.2.1 on the HPC cluster. It does not have any output itself.

#### script\_1.2.3\_normalisation\_comparison\_output\_analysis.R
R script in which the output of the normalisation comparison is analysed and the best performer is chosen. It outputs the data set normalised using FQ-normalisation. 

### denoising
#### script\_1.3.1\_run\_scGNN\_prework\_hpc.sh
Shell script that runs the prework routine of scGNN on the HPC cluster. It prepares the data for being used in scGNN by for example training the LTGM model. It outputs the data in the format necessary for running scGNN in script 1.3.2.

#### script\_1.3.2\_run\_scGNN\_clustering\_and\_imputation\_hpc.sh
Shell script that runs scGNN on the normalised (and prepared) data on the HPC cluster. Its relevant outputs are the cell-type clustering and the imputed expression data set.

#### script\_1.3.3\_scGNN\_output\_selection.R
R script that creates the relevant subset of the imputed data set by filtering out cell types based on the cell-type clustering of scGNN. The output is the imputed expression data set.

#### script\_1.3.3\_scGNN\_output\_analysis.R
R script that analyses the effects of the imputation step on the data. It outputs the results of several tests and statistics which are described in Appendix \ref{app:intermediate_results:denoising}.

## part\_2\_network\_inference
#### script\_2.1\_network\_inference\_hpc.py
Python script that builds the networks from the imputed data set by running GRNBoost2 on the HPC cluster.

#### script\_2.2\_run\_network\_inference\_hpc.sh
Shell script used for running script 2.1 on the HPC cluster. It does not have any output itself.

## part\_3\_network\_comparison
#### script\_3.1\_prep\_network\_comparison.R
R script that prepares the networks outputted by script 2.1 for the network comparison step. It converts the edge lists that GRNBoost2 outputted to weight matrices and creates \texttt{.csv} files for visualising the networks in Cytoscape.
The outputs are these weight matrices and the \texttt{.csv} files.

#### script\_3.2\_cutoff\_and\_effective\_resistance.R
R script that implements part of the network comparison. It implements the edge weight cutoff threshold and computes the effective resistance distance between all genes in the networks. It outputs the effective resistance distance matrices of both networks.

#### script\_3.3\_create\_and\_compare\_affinity\_profiles.R
R script that finishes the creation of the gene affinity profiles, implements the comparison of the profiles and reports the results. It outputs the results from the comparison and the boxplot presented in Section \ref{sec:results:NetworkComparison}.

#### script\_3.4\_visualise\_connectivity\_differences.R
R script that creates the visualisations of the connectivity differences for the neighbourhoods of the different candidate genes (see Appendix \ref{app:four_regions}) and of the overlapping genes reported in Section \ref{res:exploration_regions}.
