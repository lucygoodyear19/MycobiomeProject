###########################################################################################
############################### Phylogenetic Tree Prep ####################################
###########################################################################################

# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


###################### initial set up, library and data loading ###########################


# load packages
library(Biostrings)
library(DECIPHER)
library(phyloseq)
library(phangorn)
library(dada2)
#library(seqinr) # for saving to fasta
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
print(paste0("Data path: ", dada2_data_path))
print(paste0("Results path: ", results_path))

# load phyloseq object
dada2 <- readRDS(paste0(dada2_data_path, "physeqob_DADA2_tree.rds"))


##################################### Build tree prep #####################################


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
names(seqs) <- paste0("ASV", 1:length(seqs))
#seqs_ls <- as.list(seqs)
#names(seqs_ls) <- paste0("ASV", 1:length(seqs_ls)) # this propagates to the tip labels of the tree

# order sequences by length
#seqs_ord <- seqs_ls[order(sapply(seqs_ls, function(x) nchar(x)))]

# perform an allignment of all the sequences
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)


#################################### Write to fasta #######################################


# write to fasta file to run fast_tree
#write.fasta(sequences = seqs_ord, names = names(seqs_ord), file.out = paste0(results_path, "seqs_ord_for_tree.fasta"))

writeXStringSet(alignment, paste0(results_path, "alignment.fasta"), append=FALSE,
                compress=FALSE, format="fasta")

## end of script
