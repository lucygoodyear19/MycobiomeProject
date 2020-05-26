######################################################################
########## Ecuador: HPC docs vs Phil's docs comparison ###############
######################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1
# Date: 03/02/20

# clear workspace
rm(list=ls())


############ initial set up, library and data loading ################


# set working directory
setwd("~/Documents/MRes/MycobiomeProject/References/FromPhil/Ecuador_pipelineITS/filtN_Phil/")
# load Phil's docs
filt_phil <- read.csv("filtering_output8.csv", stringsAsFactors = FALSE)
abun_phil <- read.table("Ecuador Abundance Table new.txt", stringsAsFactors = FALSE)
metadata_phil <- read.csv("Ecuador Metadata new.txt", stringsAsFactors = FALSE, header = FALSE)
tax_phil <- read.csv("Ecuador Tax Table new.txt", stringsAsFactors = FALSE)

# set working directory
setwd("~/Documents/MRes/MycobiomeProject/Analysis/Countries_Runs/Ecuador_CR/Ecuador/HPC_Results_200401/")
# load my docs
filt_lg <- read.csv("filtering_output8.csv", stringsAsFactors = FALSE)
abun_lg <- read.table("Abundance_Table.txt", stringsAsFactors = FALSE)
metadata_lg <- read.csv("Metadata.txt", stringsAsFactors = FALSE, header = FALSE)
tax_lg <- read.table("Tax_Table.txt", stringsAsFactors = FALSE)


####################### function for comparison ###########################


# function to compare two corresponding files
compare <- function(file1, file2){
  # reorder data for comparison
  file1 <- as.character(file1[,1])
  file2 <- as.character(file2[,1])
  file1 <- plyr::arrange(file1, file1[,1])
  file2 <- plyr::arrange(file2, file2[,1])
  # generate differences between file1 and file2
  differences_1to2 <- setdiff(file1, file2)
  # subset file1 by differences
  file1_diff <- subset(file1, file1[,1] %in% differences_1to2[,1])
  # generate differences bewteen file2 and file1
  differences_2to1 <- setdiff(file2, file1)
  # subset file2 by differences
  file2_diff <- subset(file2, file2[,1] %in% differences_2to1[,1])
  # create combined dataframe with both differences
  df_compare <- cbind(file1_diff, file2_diff)
  return(df_compare)
}


########################### compare files ####################################


# compare all correpsonding files pairwise using above function
filt_diff_compare <- compare(filt_lg, filt_phil)
abun_compare <- compare(abun_lg, abun_phil)
metadata_compare <- compare(metadata_lg, metadata_phil)
tax_compare <- compare(tax_lg, tax_phil)


## end of script
