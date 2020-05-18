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
#library(ggplot2)
#library(gridExtra)

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

# load phyloseq object
dada2 <- readRDS(paste0(path_out, "physeqob_DADA2.rds"))


############################# Build tree ##################################


# sum abundances for each sequence and store as vector
seqtab <- as.data.frame(rbind(otu_table(dada2), apply(otu_table(dada2), 2, sum)))
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
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

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


############################ create phyloseq object ##############################


# add to phyloseq object
dada2 <- merge_phyloseq(dada2, 
                     phy_tree(fitGTR$tree))

# save phyloseq object to be imported into analysis scripts
saveRDS(dada2,paste0(path_out,"physeqob_DADA2_complete.rds"))


## end of script


######################## Extracting chytrids ##################################


## end of script

