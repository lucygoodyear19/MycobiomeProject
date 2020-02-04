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
setwd("~/Documents/CMEECourseWork/Project/Ecuador_pipelineITS/Plates/filtN_Phil/")
# load Phil's docs
abun_phil <- read.csv("Ecuador Abundance Table new.txt", stringsAsFactors = FALSE)
metadata_phil <- read.csv("Ecuador Metadata new.txt", stringsAsFactors = FALSE, header = FALSE)
tax_phil <- read.csv("Ecuador Tax Table new.txt", stringsAsFactors = FALSE)
filt_phil <- read.csv("filtering_output8.csv", stringsAsFactors = FALSE)

# set working directory
setwd("~/Documents/CMEECourseWork/Project/Ecuador_pipelineITS/HPC_results/Ecuador_results_lg/")
# load my docs
abun_lg <- read.csv("Ecuador_Abundance_Table_lg.txt", stringsAsFactors = FALSE)
metadata_lg <- read.csv("Ecuador_Metadata_lg.txt", stringsAsFactors = FALSE, header = FALSE)
tax_lg <- read.csv("Ecuador_Tax_Table_lg.txt", stringsAsFactors = FALSE)
filt_lg <- read.csv("filtering_output8_lg.csv", stringsAsFactors = FALSE)


####################### functions for comparison ###########################


# reorder data for comparison
filt <- plyr::arrange(filt, X, decreasing = TRUE)
filt_phil <- plyr::arrange(filt_phil, X, decreasing = TRUE)

filt_diff_l_to_p <- setdiff(filt,filt_phil)
filt_diff_lg <- subset(filt, X %in% filt_diff_l_to_p[,1])
colnames(filt_diff_lg) = c("X_lg", "reads.in_lg", "reads.out_lg")
filt_diff_p_to_l <- setdiff(filt_phil,filt)
filt_diff_phil <- subset(filt_phil, X %in% filt_diff_p_to_l[,1])
filt_diff_compare <- cbind(filt_diff_lg, filt_diff_phil)

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

filt_diff_compare <- compare(filt_lg, filt_phil)
tax_compare <- compare(tax_lg, tax_phil)
metadata_compare <- compare(metadata_lg, metadata_phil)

Sample = as.data.frame(c(1,2,3,4,4,3,2,6,7,1))
for (i in 1:nrow(Sample)) {
    if (duplicated(Sample)[i] == TRUE){
    Sample[i,1] = paste0(Sample[i,1], "a")
    }
}
