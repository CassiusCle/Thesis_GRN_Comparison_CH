################################################################################
### Name:        script_1.1_process_raw_data.R
### Description: Imports raw data and merges into an expression count dataframe
###              and metadata metaframe
### Author:      Clemens Pieter Dick Hendrickx
################################################################################

################################################################################
# Load Packages
################################################################################

LoadPackages <- function(){
  
  if(!require(readr))     {install.packages("readr")};        library(readr)
  if(!require(plyr))      {install.packages("plyr")};         library(plyr)
  if(!require(dplyr))     {install.packages("dplyr")};        library(dplyr)
  if(!require(tidyverse)) {install.packages("tidyverse")};    library(tidyverse)
  if(!require(foreach))   {install.packages("foreach")};      library(foreach)
}
LoadPackages()

################################################################################
# Hyperparameters
################################################################################
dir <- "~/Thesis/Thesis Data/Selected Data/"

################################################################################
# Auxilliary functions
################################################################################

# Function that reads in data from .csv customised to the format of my raw data
read_csv_custom <- function(directory, fileName, colType_default = 'c') {
  colNames <- read_delim(paste0(directory, fileName),
                         delim = ' ', n_max = 1, col_names = FALSE) %>% t %>% as.vector
  colNames <- c('index', colNames)
  read_delim(paste0(directory, fileName),
             delim = ' ', skip = 1, col_names = colNames, col_types = cols(.default = colType_default)) %>% as.data.frame
}

# Function that counts the number of NA and unique entries of each column
count_NAs_Uniques <- function(df) foreach(col = colnames(df), .combine = 'rbind') %do% 
  list(colName = col, missingValues = sum(is.na(df[,col])), 
       uniqueValues = length(unique(df[,col]))) %>% {data.frame(.[,-1], row.names = .[,1])}

# Function that removes NA values from a vector
noNAs <- function(vector) vector[!is.na(vector)]

# Function that finds the columns in a dataframe
find_cols_with_pattern <- function(df, pattern, ignoreCase = TRUE) {
  foreach(col = colnames(df), .combine = 'rbind') %do% {
    list(colName = col, numMatches = sum(grepl(pattern, df[,col], ignore.case = ignoreCase)))
  } %>% data.frame %>%  filter(numMatches > 0) %>% {.[,'colName']} %>% unlist
}

################################################################################
# Load data
################################################################################
gc()

# Read in count data
colNames <- read_csv(paste0(dir, "ACR_epithelial_brushes_counts.csv"),
                     n_max = 1, col_names = FALSE) %>% t %>% as.vector

counts_temp <- read_csv(paste0(dir, "ACR_epithelial_brushes_counts.csv"),
                       skip = 1, col_names = FALSE)

counts <- data.frame(counts_temp[,-1])
rownames(counts) <- unlist(counts_temp[,1])
colnames(counts) <- colNames
remove(counts_temp, colNames)

# Read in metadata
metadata <- read_csv_custom(dir, "ACR_epithelial_brushes_metadata_w_clusterIdents.csv")
EpithelialPlus_sampleInfo <- read_csv_custom(dir, "EpithelialPlus_sampleInfo.csv")
rownames(EpithelialPlus_sampleInfo) <- EpithelialPlus_sampleInfo$seqID

################################################################################
# Merge metadatasets
################################################################################

#----------- Prepare metadata dataframe ----------------------------------------
count_NAs_Uniques(metadata) # No missing values

# Extract features from index variable
metadata$seqID <- str_extract(metadata$index, '[0-9]{5}_[0-9]#[0-9]{1,3}') 
metadata$location <-  str_extract(metadata$index, '^[A-Z][0-9]{1,2}')      
metadata$plateThisDonor <- str_extract(str_extract(metadata$index, '_P[0-9]{1,2}_?$'), 'P[0-9]{1,2}')
count_NAs_Uniques(metadata) %>% print   # 2407 obs. with seqID, 2346 obs. with location and plateThisDonor

# Format location to always contain two digits
loc_toFix <- metadata$location[grepl('^[A-Z][0-9]{1}$', metadata$location)] 
metadata$location[grepl('^[A-Z][0-9]{1}$', metadata$location)] <- paste0(str_extract(loc_toFix, '^[A-Z]'),'0', str_extract(loc_toFix, '[0-9]$'))
metadata$location %>% unique %>% print # Always two digits

# Format plateThisDonor to always contain two digits
plate_toFix <- metadata$plateThisDonor[grepl('P[0-9]{1}$', metadata$plateThisDonor)] 
metadata$plateThisDonor[grepl('P[0-9]{1}$', metadata$plateThisDonor)] <- paste0('P0', str_extract(plate_toFix, '[0-9]$'))
metadata$plateThisDonor %>% unique %>% print # Always two digits

# Construct MergeID key using donor, location and plateThisDonor
unique(metadata$donor) # Donors all in format: ARMS[0-9]{3} <-- Note all have three digits
metadata$MergeID <- paste0(metadata$donor, metadata$location, metadata$plateThisDonor)

# Remove incomplete MergeIDs 
mergeID_incomplete <- !grepl('ARM?S[0-9]{3}[A-Z][0-9]{2}P[0-9]{2}', metadata$MergeID) # MergeID not good format (due to NAs)
mergeID_incomplete %>% sum %>% print # Out: 2406
metadata$MergeID[mergeID_incomplete] <- NA
count_NAs_Uniques(metadata) # MergeID has 2347 unique values (1 is NA), other 2406 obs. are identified by seqID


#----------- Prepare EpithelialPlus_sampleInfo dataframe -----------------------
count_NAs_Uniques(EpithelialPlus_sampleInfo) 

# Rename index to avoid confusion
EpithelialPlus_sampleInfo <- EpithelialPlus_sampleInfo %>% dplyr::rename(index_epith = index)

# Extract donorID feature
pattern_donorID <- 'ARM?S[0-9]{2,3}'
cols_donorID <- find_cols_with_pattern(EpithelialPlus_sampleInfo, pattern_donorID) %>% print # Output: "V2"
EpithelialPlus_sampleInfo$donorID <- str_extract(EpithelialPlus_sampleInfo$V2, pattern_donorID)
EpithelialPlus_sampleInfo$donorID %>% unique # Some donorID have two instead of three digits, one has ARS instead of ARMS

# Format donorID to always have three digits and to have prefix ARMS
EpithelialPlus_sampleInfo <- EpithelialPlus_sampleInfo %>% mutate(donorID = case_when(
  grepl('ARM?S[0-9]{2}$', donorID, ignore.case = TRUE) ~ paste0(str_extract(donorID, 'ARM?S'), '0', str_extract(donorID, '[0-9]{2}$')),
  grepl('ARM?S[0-9]{3}$', donorID, ignore.case = TRUE) ~ donorID
))
donorID_toFix <- EpithelialPlus_sampleInfo$donorID[grepl('ARS', EpithelialPlus_sampleInfo$donorID)]
EpithelialPlus_sampleInfo$donorID[grepl('ARS', EpithelialPlus_sampleInfo$donorID)] <- paste0('ARMS',str_extract(donorID_toFix, '[0-9]{3}$'))
EpithelialPlus_sampleInfo$donorID %>% unique # All donorID have three digits
EpithelialPlus_sampleInfo$donorID %>% unique %>% length # 30 - 1 donors (one is NA)

# Format location to always contain two digits
EpithelialPlus_sampleInfo$location %>% unique # Some have one digit instead of two
loc_toFix <- EpithelialPlus_sampleInfo$location[grepl('^[A-Z][0-9]{1}$', EpithelialPlus_sampleInfo$location)] 
EpithelialPlus_sampleInfo$location[grepl('^[A-Z][0-9]{1}$', EpithelialPlus_sampleInfo$location)] <- paste0(str_extract(loc_toFix, '^[A-Z]'),
                                                                         '0', str_extract(loc_toFix, '[0-9]$'))
EpithelialPlus_sampleInfo$location %>% unique # Always two digits

# Extract plateThisDonor feature
find_cols_with_pattern(EpithelialPlus_sampleInfo, '_P[0-9]{1,2}_?$') %>% print # Out: "V2"
EpithelialPlus_sampleInfo$plateThisDonor <- str_extract(str_extract(EpithelialPlus_sampleInfo$V2, '_P[0-9]{1,2}_?$'), 'P[0-9]{1,2}')

# Format plateThisDonor to always contain two digits
plate_toFix <- EpithelialPlus_sampleInfo$plateThisDonor[grepl('P[0-9]{1}$', EpithelialPlus_sampleInfo$plateThisDonor)] 
EpithelialPlus_sampleInfo$plateThisDonor[grepl('P[0-9]{1}$', EpithelialPlus_sampleInfo$plateThisDonor)] <- paste0('P0', str_extract(plate_toFix, '[0-9]$'))
EpithelialPlus_sampleInfo$plateThisDonor %>% unique # Always two digits

# Construct MergeID key using donor, location and plateThisDonor
EpithelialPlus_sampleInfo$MergeID <- paste0(EpithelialPlus_sampleInfo$donorID, EpithelialPlus_sampleInfo$location, EpithelialPlus_sampleInfo$plateThisDonor) 
count_NAs_Uniques(EpithelialPlus_sampleInfo) # MergeID has 7170 unique values (1 is NA)

# Remove incomplete MergeIDs 
mergeID_incomplete <- !grepl('ARMS[0-9]{3}[A-Z][0-9]{2}P[0-9]{2}', EpithelialPlus_sampleInfo$MergeID) # MergeID not good format (due to NAs)
mergeID_incomplete %>% sum %>% print # 1962
EpithelialPlus_sampleInfo$MergeID[mergeID_incomplete] <- NA
count_NAs_Uniques(EpithelialPlus_sampleInfo) # MergeID has 6776 unique values (1 is NA)

EpithelialPlus_sampleInfo$MID_is_duplicate[!is.na(EpithelialPlus_sampleInfo$MergeID)] <- duplicated(EpithelialPlus_sampleInfo$MergeID[!is.na(EpithelialPlus_sampleInfo$MergeID)])
EpithelialPlus_sampleInfo$MID_is_duplicate[is.na(EpithelialPlus_sampleInfo$MID_is_duplicate)] <- FALSE
EpithelialPlus_sampleInfo$MID_is_duplicate %>% sum %>% print # 479 duplicates
unique(EpithelialPlus_sampleInfo$MergeID[EpithelialPlus_sampleInfo$MID_is_duplicate]) %>% length # Real duplicates 384

noNAs(unique(EpithelialPlus_sampleInfo$MergeID)) %>% length # Unique 6775 (excl the NA in this case), some are duplicated 
non_duplicated_uniques <- setdiff(noNAs(unique(EpithelialPlus_sampleInfo$MergeID)), unique(EpithelialPlus_sampleInfo$MergeID[EpithelialPlus_sampleInfo$MID_is_duplicate]))
non_duplicated_uniques %>% length %>% print  # 6391 total
EpithelialPlus_sampleInfo$MID_not_duplicated <- FALSE
EpithelialPlus_sampleInfo$MID_not_duplicated[EpithelialPlus_sampleInfo$MergeID %in% non_duplicated_uniques] <- TRUE

#----------- Merge dataframes --------------------------------------------------
metadata_complete <- rbind.fill(merge(x = filter(metadata, !is.na(MergeID)), 
                                      y = filter(EpithelialPlus_sampleInfo, !is.na(MergeID) & MID_not_duplicated), 
                                      all.x = TRUE, all.y = FALSE,
                                      by.x = 'MergeID', by.y = 'MergeID'),
                                filter(metadata, is.na(MergeID))
                                )

# Reduce number of rows and columns
metadata_complete <- metadata_complete %>% filter(!is.na(seqID) | !is.na(seqID.y)) # Remove observations without seqID
metadata_complete <- metadata_complete %>% select_if(~!(all(is.na(.)) | all(. == ""))) # Remove empty columns (28 -> 23)

metadata_complete$seqID[is.na(metadata_complete$seqID)] <- metadata_complete$seqID.y[is.na(metadata_complete$seqID)] # Combine seqID's
count_NAs_Uniques(metadata_complete) # --> seqID now available for all observations

# Extract features seqRun, run and lane (batch information)
metadata_complete$seqRun <- str_extract(metadata_complete$seqID, '[0-9]+_[0-9]')
metadata_complete$run <-  str_extract(metadata_complete$seqRun, '^[0-9]+') 
metadata_complete$lane <-  str_extract(metadata_complete$seqRun, '[0-9]+$') 

metadata_complete$seqRun %>% unique %>% length %>% print # 19 run and lane combinations ==> batches
metadata_complete$run %>% unique %>% length %>% print # 4 runs

# Construct definitive metadata and count dataframes
metadata <- metadata_complete %>% select(index, donor, healthStatus, sex, cluster, seqRun, seqID, run, lane)
rownames(metadata) <- metadata$index
count_NAs_Uniques(metadata) # No missing values

counts <- counts[metadata$index] # Select only cells present in metadata dataframe

# Make both use the same key
colnames(counts) <- metadata[colnames(counts),]$seqID 
rownames(metadata) <- metadata$seqID

# Remove intermediate objects
remove(EpithelialPlus_sampleInfo, metadata_complete, 
       cols_donorID, donorID_toFix, loc_toFix, plate_toFix,
       mergeID_incomplete, non_duplicated_uniques, pattern_donorID,
       count_NAs_Uniques, find_cols_with_pattern, noNAs, read_csv_custom)

# Save metadata and count objects
saveRDS(counts,   paste0(dir,'/Intermediate data/', 'expression_counts.rds'))
saveRDS(metadata, paste0(dir,'/Intermediate data/', 'expression_metadata.rds'))

