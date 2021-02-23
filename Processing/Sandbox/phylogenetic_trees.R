###########################################################################################
############################# Phylogenetic Tree Generator #################################
###########################################################################################

# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


###################### initial set up, library and data loading ###########################


# load packages
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

#readDNAStringSet(paste0(dada2_data_path, "alignment.fasta"), format="fasta")


##################################### Build tree prep #####################################


# create tree subsetting by order
chytrids_na <- subset_taxa(dada2, is.na(Order))
chytrids_rhi <- subset_taxa(dada2, Order == "o__Rhizophydiales")
chytrids <- merge_phyloseq(chytrids_rhi, chytrids_na)

# sum abundances for each sequence and store as vector
seqtab <- as.data.frame(rbind(otu_table(chytrids), apply(otu_table(chytrids), 2, sum)))
abun_sums <- seqtab[nrow(seqtab),]
abundance <- as.numeric(abun_sums)

# extract seqeunces
sequence <- colnames(seqtab)

# save sequences as total abundances as new dataframe
seqsfortree <- cbind(sequence,abundance)
seqsfortree <- as.data.frame(seqsfortree)
seqsfortree$abundance <- as.numeric(as.character(seqsfortree$abundance))

# source for following code: https://f1000research.com/articles/5-1492

# get sequences
seqs <- getSequences(seqsfortree)
names(seqs) <- seqs

#seqs_ls <- as.list(seqs)
#names(seqs_ls) <- paste0("ASV", 1:length(seqs_ls)) # this propagates to the tip labels of the tree

# order sequences by length
#seqs_ord <- seqs_ls[order(sapply(seqs_ls, function(x) nchar(x)))]


###################################### Fasttree ###########################################


# write to fasta file to run fast_tree
#write.fasta(sequences = seqs_ord, names = names(seqs_ord), file.out = paste0(results_path, "seqs_ord_for_tree.fasta"))


######################################## Phangorn ########################################


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


#chytrid_tree <- readRDS(paste0(dada2_data_path, "Tree/chytrid_tree.rds"))


############################ create phyloseq object ##############################


# add to phyloseq object
phy_tree(chytrids) <- fitGTR$tree

#taxa_names(chytrids) <- refseq(chytrids)

# save phyloseq object to be imported into analysis scripts
saveRDS(dada2,paste0(results_path,"physeqob_DADA2_chytrids_complete.rds"))


## end of script


######################## Extracting chytrids ##################################


## end of script

