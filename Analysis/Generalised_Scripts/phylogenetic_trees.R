######################################################################
################### Phylogenetic Tree Generator ######################
######################################################################

# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


############ initial set up, library and data loading ################


# load packages
library(DECIPHER)
library(phyloseq)
library(phangorn)
library(dada2)
library(ggplot2)
library(gridExtra)

# personal laptop
setwd("~/Documents/MRes/MycobiomeProject/Analysis/Countries_Runs/Ecuador_CR/Ecuador/HPC_Results_200401/")
# hpc
#setwd(path_out)

# load my docs
abun <- read.csv("Abundance_Table.txt", stringsAsFactors = FALSE)
tax <- read.table("Taxa_Table.txt", stringsAsFactors = FALSE)
seqtab <- read.table("Seq_Abun_Table.txt", stringsAsFactors = FALSE)
info <- read.csv("Metadata.csv", stringsAsFactors = FALSE)


############################# Build tree ##################################


# sum abundances for each sequence and store as vector
seqtab[nrow(seqtab)+1,] = apply(seqtab, 2, sum)
abun_sums <- seqtab[nrow(seqtab),]
abundance <- as.numeric(abun_sums)

# extract seqeunces
sequence <- colnames(seqtab)

# save sequences as total abundances as new dataframe
seqsfortree <- cbind(sequence,abundance)
seqsfortree <- as.data.frame(seqsfortree)

# source for following code: https://f1000research.com/articles/5-1492

# get sequences
seqs <- getSequences(seqsfortree)
names(seqs) <- seqs # this propagates to the tip labels of the tree
# perform an allignment of all the sequences
alignment <- AlignSeqs(DNAStringSet(seqss), anchor=NA)

# transform alignment data into correct format for phangorn
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
# 
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
# compute likelihood of tree
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
# optimise tree fit parameters
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
# unload phangorn package
detach("package:phangorn", unload=TRUE)


######################## Create phyloseq object #############################


physeqob <- phyloseq(tax_table(tax), otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data, phy_tree(fitGTR))
saveRDS(physeqob,paste(path_out,"physeqob.rds",sep=''))


## end of script

