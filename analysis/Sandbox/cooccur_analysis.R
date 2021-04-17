##############################################################################
################### DADA2 and Esto Cooccurence Analysis ######################
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

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) # setup to accept arguments from command line
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied in the form of an R-script containing the following arguments: 
       1) dada2_data_path - path to directory containing input data (phyloseq object)
       2) results_path - path to results directory")
}
# source analysis arguments
source(args)
# print arguments as check
print(paste0("Data path: ", dada2_data_path))
print(paste0("Results path: ", results_path))

# load phyloseq object
print("Loading data...")
cooccur_full <- readRDS(paste0(results_path, "Cooccur/cooccur_full_TW16.rds"))
dada2 <- readRDS(paste0(dada2_data_path, "Cooccur/physeqob_DADA2_TW16.rds"))


##############################################################################
################################# Bd Analysis ################################


# look at cooccurence patterns just for Bd
cooccur_bd <- pair(cooccur_full, spp="bd")

# view all significant results for bd
print(cooccur_bd)

# add columns for +ve/-ve coocurence and save taxonomic assignment to separate df
tax <- as.data.frame(matrix(0, nrow(cooccur_bd), 7))
rownames(cooccur_bd) <- cooccur_bd$sp2
for (i in 1:nrow(cooccur_bd)) {
  if (cooccur_bd$p_gt[i] < 0.05) {
    cooccur_bd$cooccurance[i] <- "Positive"
  }
  if (cooccur_bd$p_lt[i] < 0.05) {
    cooccur_bd$cooccurance[i] <- "Negative"
  }
  tax[i,] <- tax_table(dada2)[rownames(tax_table(dada2)) == cooccur_bd$sp2[i]]
}
# set columns names
names(tax) <- colnames(tax_table(dada2))
# set rownames to match ASVs in cooccur object
rownames(tax) <- rownames(cooccur_bd)
# merge dataframes to include taxonomic assignment
cooccur_bd_tax <- merge(cooccur_bd, tax, by="row.names")
# remove unneccesary columns
cooccur_bd_tax_fil <- cooccur_bd_tax[,-c(2:6,10)]
# set first column name to ASV
names(cooccur_bd_tax_fil)[1] <- "ASV"

# extract Bd column from sample data
samps <- as(sample_data(dada2), "data.frame")
for (sam in 1:nrow(samps)) {
  if (!is.na(samps$Bd[sam]) & samps$Bd[sam] == "S") {samps$Bd[sam] <- 0}
}
bd <- samps$Bd

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

# plot heatmap for sginificant ASVs by sample
# subset to only include those significant with regards to Bd
c <- subset(cooccur_object_t, rownames(cooccur_object_t) %in% cooccur_bd_tax_fil$ASV)
# add column for ASVs
ASV <- row.names(c)
df <- cbind(ASV, c)
library(tidyr)
# gather into format required for heat map
melted_df <- gather(df, key="Sample", value = "Abundance", names(df)[2]:names(df)[ncol(df)])
# view data
head(melted_df)
library(ggplot2)
# plot heatmap for Bd cooccurance significant taxa
ggplot(data = melted_df, aes(x=Sample, y=ASV, fill=Abundance)) + 
  geom_tile()

# random package functions
summary(cooccur_bd)
plot(cooccur_bd)
pa <- pair.attributes(cooccur_full)
subset(pa, sppname == "bd")
pp <- pair.profile(cooccur_full)

# subset by positive and negative cooccurence
co_tax_pos <- subset(cooccur_bd_tax_fil, cooccurance == "Positive")
co_tax_neg <- subset(cooccur_bd_tax_fil, cooccurance == "Negative")

# negatively and positively associated taxa ordered by p-value
co_tax_neg_ord <- co_tax_neg[order(co_tax_neg$p_lt),]
co_tax_pos_ord <- co_tax_pos[order(co_tax_pos$p_gt),]

# negatively and positively associated taxa grouped by taxa
co_tax_neg_group <- co_tax_neg[order(co_tax_neg$Phylum, co_tax_neg$Class, co_tax_neg$Order, co_tax_neg$Family, co_tax_neg$Genus, co_tax_neg$Species),]
co_tax_pos_group <- co_tax_pos[order(co_tax_pos$Phylum, co_tax_pos$Class, co_tax_pos$Order, co_tax_pos$Family, co_tax_pos$Genus, co_tax_pos$Species),]

# plot full cooccurance result as a heatmap
plot(cooccur_full)


## end of script
