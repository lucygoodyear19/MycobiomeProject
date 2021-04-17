########################################################################################
#################################### Subset FastTree ###################################
########################################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


########################################################################################
###################################### Set up ##########################################


# load packages
library(phyloseq)
library(Biostrings)
#library(ggtree)
library(ggplot2)
#library(gridExtra)
library(phytools)
library(adephylo) # for listTips() function

# load filtered phyloseq object for both DADA2 pipeline
dada2 <- readRDS(paste0(dada2_data_path, "physeqob_DADA2.rds"))

# load FastTree tree
tree <- read_tree(paste0(dada2_data_path, "fasttree_tree"))

# merge fasttree with dada2 phyloseq object
dada2 <- merge_phyloseq(dada2, 
                        phy_tree(tree))


###################################################################################
######################### Subset phylogenetic tree ################################


# subset tree via (mostly) taxonomically independent method

# add Genus/species name columns for tip labels
temp <- as(tax_table(dada2),"matrix")
tempdf <- as.data.frame(temp)
tempdf$name<-paste0(tempdf[,"Genus"], "_",tempdf[,"Species"])
tempdf$name <- gsub("g__", "", tempdf$name)
tempdf$name <- gsub("s__", "", tempdf$name)
tempdf$name <- gsub("NA_NA", "", tempdf$name)
tempdf <- as.matrix(tempdf)
dada2 <- merge_phyloseq(dada2, tax_table(tempdf))

# create tree subsetting by order
chytrids = subset_taxa(dada2, Order == "o__Rhizophydiales")
chytrids_pruned = prune_samples(sample_sums(chytrids) > 10, chytrids)

# generate previously subsetted chytrid tree to find ASVs to use as ends
#plotTree(phy_tree(chytrids_pruned), node.numbers=TRUE)
#nodelabels() 

# find the type labels at each end
one <- rownames(as.data.frame(listTips(phy_tree(chytrids_pruned))[[1]]))[1]
two <- rownames(as.data.frame(listTips(phy_tree(chytrids_pruned))[[length(listTips(phy_tree(chytrids_pruned)))]]))[1]

# use ASV ends to find MRCA node number
mrca_node <- getMRCA(as.phylo(phy_tree(dada2)), c(one, two))
# extract clade from node number
clade_chytrids <- extract.clade(as.phylo(phy_tree(dada2)), mrca_node)
# subset dada2 by clade
chytrids_final = subset_taxa(dada2, taxa_names(dada2) %in% clade1$tip.label)

# save chytrids to separate phyloseq object
saveRDS(chytrids, paste0(results_path, "physeqob_DADA2_chytrids.rds"))


## end of script
