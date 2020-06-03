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
       1) path containing results from DADA2 pipeline")
}
# load arguments into script
source(args)
# print arguments as check
print(path_out)

# load phyloseq object for DADA2 pipeline
dada2 <- readRDS(paste0(path_out, "physeqob_DADA2.rds"))
print("Data loaded")


##############################################################################
########################### Calculate cooccurence ############################


# extract Bd column from sample data
bd <- sample_data(dada2)$Bd

# add Bd binary column to otu_table
cooccur_object <- cbind(as.data.frame(otu_table(dada2)), bd)

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

# main coocur function to calculate pairwise cooccurence patterns
print("Cooccur function start:")
Sys.time()
cooccur_full <- cooccur(cooccur_object_t, 
                        type="spp_site",
                        thresh=TRUE, 
                        spp_names=TRUE)
print("Cooccur function finished:")
Sys.time()


cooccur_full <- readRDS("/Users/lucy/Documents/MRes/MycobiomeProject/Analysis/Results/Taiwan_2016/Cooccur/cooccur_full_dada2.rds")


##############################################################################
################################# Bd Analysis ################################


# look at cooccurence patterns just for Bd
cooccur_bd <- pair(cooccur_full, spp="bd")

# view all significant results for bd
print(cooccur_bd)
summary(cooccur_bd)
pair.profile()


sort(cooccur_bd, )

# plot as a heatmap
heatmap_bd <- plot(cooccur_bd)


##############################################################################
################################## Save outputs ##############################


saveRDS(cooccur_full, paste0(path_out,"cooccur_full_dada2.rds"))


## end of script
