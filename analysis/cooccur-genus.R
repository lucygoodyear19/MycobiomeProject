##############################################################################
######################### Cooccurence - genus level ##########################
##############################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


##############################################################################
################################## Set up ####################################

# load packages
library("phyloseq")
library("cooccur")

# import arguments to run script on specific country data
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) # setup to accept arguments from command line
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied in the form of an R-script containing the following arguments: 
       1) dada2_data_path - path to directory containing input data (phyloseq object)
       2) results_path - path to results directory")
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


# add genus to abundance table
abun <- as.data.frame(otu_table(dada2))
abun_t <- as.data.frame(t(abun))
taxa <- as.data.frame(tax_table(dada2))
genus <- taxa$Genus
genus_abun <- cbind(genus, abun_t)

# merge by samples by genus
genus_grid <- plyr::ddply(genus_abun, "genus", plyr::numcolwise(sum))

# remove NA row
genus_grid <- subset(genus_grid, !is.na(genus))

# rename genus to remove prefix
genus_grid$genus <- gsub("g__", "", genus_grid$genus)

# rename rows to genus 
rownames(genus_grid) <- genus_grid$genus

#remove genus column
genus_grid <- genus_grid[,-1]

# transpose to add Bd column
genus_grid_t <- as.data.frame(t(genus_grid))

# extract Bd column from sample data
samps <- as(sample_data(dada2), "data.frame")
bd <- samps$Bd

# add Bd binary column to otu_table
cooccur_object <- cbind(genus_grid_t, bd)

# remove samples with NA Bd qPCR
cooccur_object <- cooccur_object[!is.na(cooccur_object$bd),]

# set Bd column to numeric
cooccur_object[,ncol(cooccur_object)] <- sapply(cooccur_object[,ncol(cooccur_object)], as.character)
cooccur_object[,ncol(cooccur_object)] <- sapply(cooccur_object[,ncol(cooccur_object)], as.numeric)

# set all abundances > 0 to 1 to get a presence absence table
cooccur_object[cooccur_object > 0] <- 1

# remove all taxa found only in one sample
cooccur_object_t <- as.data.frame(t(cooccur_object))
cooccur_object_t$sum <- rowSums(cooccur_object_t)
cooccur_object_t <- subset(cooccur_object_t, sum >= 2)
cooccur_object_t <- subset(cooccur_object_t, select=-c(sum))

# main coocur function to calculate pairwise cooccurence patterns
print("Cooccur function start:")
Sys.time()
cooccur_genus <- cooccur(cooccur_object_t, 
                        type="spp_site",
                        thresh=TRUE, 
                        spp_names=TRUE)
print("Cooccur function finished:")
Sys.time()

# save output
saveRDS(cooccur_genus, paste0(results_path,"cooccur_genus_dada2.rds"))


## end of script
