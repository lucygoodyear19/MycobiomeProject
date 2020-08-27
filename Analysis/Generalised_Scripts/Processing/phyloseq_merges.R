######################################################################################
################################### Phyloseq merges ##################################
######################################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear worksapce
rm(list=ls())


###################################### Set up #########################################


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) # setup to accept arguments from command line
# test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("Two or more phyloseq objects must be supplied as arguments for merging")
}

# load phyloseq objects and save to list
physeqs <- c()
print("Arguments supplied:")
for (obj in 1:length(args)) {
  print(args[obj])
  assign(paste0("arg", obj), readRDS(args[obj]))
  physeqs <- c(physeqs, get(paste0("arg", obj)))
}

#load phyloseq package
library("phyloseq")

# list the sample data categories wanted for merge
metadata_global <- c("MiSeqCode", 
                     "A_Order",
                     "A_Family",
                     "A_Genus_Species",
                     "Country", 
                     "Lifestage", 
                     "Elevation_m",
                     "Latitude",
                     "Longitude",
                     "Bd",
                     "Alpha_Shannon")


########################## Prepare sample data for merging ############################


# for loop to ensure sample_data contains all the required columns and subset by them
for (obj in 1:length(physeqs)) {
  physeqob <- physeqs[[obj]] # unlist phyloseq object
  samp <- as.data.frame(as.matrix((sample_data(physeqob)))) # extract sample_data as dataframe
  print(paste0("Number of samples in phyloseq object ", obj, " : ", nrow(samp)))
  for (data in 2:length(metadata_global)) {
    if (is.null(samp[[metadata_global[data]]])) {
      samp[[metadata_global[data]]] <- NA # if metadata column doesn't exist, create and fill with NA
    }
  }
  samp <- subset(samp, select=metadata_global) # subset by required metadata
  sample_data(physeqob) <- sample_data(samp) # set subsetted metadata as sample data
  physeqs[[obj]] <- physeqob # update list with updated phyloseq object
}

# merge first two phyloseq objects
dada2 <- merge_phyloseq(physeqs[[1]], physeqs[[2]])
# merge all phyloseq objects
if (length(physeqs) > 2) {
  for (obj in 3:length(physeqs)) {
    dada2 <- merge_phyloseq(dada2, physeqs[[obj]])
  }
}
print(paste0("(Check) Number of samples in combined phyloseq object: ", nrow(sample_data(dada2))))

# save with raw sequences for phylogenetic tree script
saveRDS(dada2, "physeqob_DADA2_tree.rds")

# save sequences to refseqs slot and rename ASVs for convenience
asv_seqs <- Biostrings::DNAStringSet(taxa_names(dada2))
names(asv_seqs) <- taxa_names(dada2)
dada2 <- merge_phyloseq(dada2, asv_seqs)
taxa_names(dada2) <- paste0("ASV", seq(ntaxa(dada2)))

# print out phyloseq object to screen
print("View merged phyloseq object: ")
dada2

# save phyloseq object to be imported into analysis scripts
saveRDS(dada2, "physeqob_DADA2.rds")

print("Script completed")

## end of script
