##############################################################################
######################### DADA2 and Esto Cooccurence #########################
##############################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


##############################################################################
################################## Set up ####################################

# load packages
library(phyloseq)
library(cooccur)

# import arguments to run script on specific country data
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) # setup to accept arguments from command line
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied in the form of an R-script containing the following arguments: 
       1) dada2_data_path - path to directory containing input data (phyloseq object)")
}
# load arguments into script
source(args)
# print arguments as check
print(dada2_data_path)
print(results_path)

# load phyloseq object for DADA2 pipeline
dada2 <- readRDS(paste0(dada2_data_path, "physeqob_DADA2.rds"))

print("Data loaded")


##############################################################################
################################ Co-occurence ################################


# extract Bd column from sample data
samps <- as(sample_data(dada2), "data.frame")
bd <- samps$Bd

# add Bd binary column to otu_table
cooccur_object <- cbind(as.data.frame(otu_table(dada2)), bd)

# remove samples with NA Bd qPCR
cooccur_object <- cooccur_object[!is.na(cooccur_object$bd),]

# set Bd column to numeric
cooccur_object[,ncol(cooccur_object)] <- sapply(cooccur_object[,ncol(cooccur_object)], as.character)
cooccur_object[,ncol(cooccur_object)] <- sapply(cooccur_object[,ncol(cooccur_object)], as.numeric)

# set all abundances > 0 to 1 to get a presence absence table
cooccur_object[cooccur_object > 0] <- 1

# remove all ASVs found only in one sample
cooccur_object_t <- as.data.frame(t(cooccur_object))
cooccur_object_t$sum <- rowSums(cooccur_object_t)
cooccur_object_t <- subset(cooccur_object_t, sum >= 2)
cooccur_object_t <- subset(cooccur_object_t, select=-c(sum))
# save object for use in analysis
saveRDS(cooccur_object_t, paste0(results_path,"cooccur_full_dataframe.rds"))

# main cooccur function to calculate pairwise co-occurrence patterns
print("Cooccur function start:")
Sys.time()
cooccur_full <- cooccur(cooccur_object_t, 
                        type="spp_site",
                        thresh=TRUE, 
                        spp_names=TRUE)
print("Cooccur function finished:")
Sys.time()

# save co-occur object
saveRDS(cooccur_full, paste0(results_path,"cooccur_full_dada2.rds"))


## end of script
