##############################################################################
########################### Plot Phylogenetic Trees ##########################
##############################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


##############################################################################
################################## Set up ####################################


# set wd if needed
setwd("~/Documents/MRes/MycobiomeProject/Analysis/Runs_Countries/Taiwan_Vietnam_2016/Taiwan")

# load packages
library(phyloseq)
library(Biostrings)
library(ggtree)
library(ggplot2)
library(gridExtra)

# load filtered phyloseq object for both Estonian pipeline and DADA2 pipeline
esto <- readRDS("physeqob_jen.rds")
#dada2 <- readRDS("physeqob_dada2.rds")
dada2 <- readRDS("DADA2_Results/physeqob_DADA2_complete.rds")

asv_seqs <- Biostrings::DNAStringSet(taxa_names(dada2))
names(asv_seqs) <- taxa_names(dada2)
dada2 <- merge_phyloseq(dada2, asv_seqs)
taxa_names(dada2) <- paste0("ASV", seq(ntaxa(dada2)))
dada2


###################################################################################
######################### Create phylogenetic tree ################################


# create tree subsetting by order
chytrids = subset_taxa(dada2, Order == "o__Rhizophydiales")
chytrids_pruned = prune_samples(sample_sums(chytrids) > 4, chytrids)

pdf("chytrid_tree.pdf")
plot_tree(chytrids_pruned, 
          shape = "Family", 
          label.tips = "Species", 
          #size = "abundance", 
          plot.margin = 0.5, 
          ladderize = TRUE,
          text.size = 3,
          color = "Genus")
#coord_polar(theta="y")
dev.off()

# create tree via taxonomically independent method

#library(TreeTools)
#tree <- Preorder(phy_tree(dada2))
#plot(tree)
#ape::nodelabels()
#ape::nodelabels(40, 40, bg='yellow')

#plot(Subtree(tree, 40))

# add Genus/species name columns for tip labels
temp <- as(tax_table(dada2),"matrix")
tempdf <- as.data.frame(temp)
tempdf$name<-paste0(tempdf[,"Genus"], "_",tempdf[,"Species"])
tempdf$name <- gsub("g__", "", tempdf$name)
tempdf$name <- gsub("s__", "", tempdf$name)
tempdf$name <- gsub("NA_NA", "", tempdf$name)
tempdf <- as.matrix(tempdf)
dada2 <- merge_phyloseq(dada2, tax_table(tempdf))

# generate previously subsetted chytrid tree to find ASVs to use as ends
plotTree(phy_tree(chytrids_pruned))
nodelabels() 

library(adephylo)
one <- rownames(as.data.frame(listTips(phy_tree(chytrids_pruned))[[1]]))[1]
two <- rownames(as.data.frame(listTips(phy_tree(chytrids_pruned))[[length(listTips(phy_tree(chytrids_pruned)))]]))[1]

# use ASV ends to find MRCA node number
mrca_node <- getMRCA(as.phylo(phy_tree(dada2)), c("ASV2776", "ASV2365")) #one, two))
# extract clade from node number
clade1 <- extract.clade(as.phylo(phy_tree(dada2)), mrca_node)
# subset dada2 by clade
chytrids1 = subset_taxa(dada2, taxa_names(dada2) %in% clade1$tip.label)

# plot tree
pdf("Analysis_Results/chytrid_tree_tax_ind.pdf", 8, 25)
plot_tree(chytrids1,           
          shape = "Family", 
          label.tips = "name", 
          #size = "abundance", 
          plot.margin = 0.5, 
          ladderize = TRUE,
          text.size = 3,
          color = "Genus")
dev.off()

#p <- ggtree(chytrids, ladderize = FALSE) +
#  geom_tiplab(aes(label=Genus), hjust=-.3) +
#  geom_point(aes(x=x+hjust, color=Species, shape=Family, size=Abundance),na.rm=TRUE) +
#  theme(legend.position="right") + ggtitle("reproduce phyloseq by ggtree")
#print(p)

#p2 <- open_tree(p, 180)
#print(p2)

#plot_heatmap(chytrids1, na.value = "black")

#x <- plot_heatmap(chytrids_pruned)


###############################################################################
#################### BLAST chytrids against NCBI database #####################

#ASV2292, ASV2580
seqs <- as.data.frame(refseq(chytrids1))
seq1 <- seqs[rownames(seqs) == "ASV2580",]
seq2 <- seqs[rownames(seqs) == "ASV2292",]

refseq(chytrids1)

taxa <- assignTaxonomy(refseq(chytrids1), unite.ref, multithread = TRUE, tryRC = TRUE)

# add new assignments to phyloseq object
chytrids1 <- merge_phyloseq(chytrids1, tax_table(taxa))


################################################################################
########### Create phylogenetic tree including new assignments #################


pdf("Analysis_Results/chytrid_tree_tax_ind_new.pdf", 8, 25)
plot_tree(chytrids1,           
          shape = "Family", 
          label.tips = "name", 
          #size = "abundance", 
          plot.margin = 0.5, 
          ladderize = TRUE,
          text.size = 3,
          color = "Genus")
dev.off()


## end of script

