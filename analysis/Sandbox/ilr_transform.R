##############################################################################
############################## ilr transformation #############################
##############################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


##############################################################################
################################## Set up ####################################


# load packages
library("dplyr")
library("tidyr") # for separate function
library("phyloseq")
library("microbiome") # for abundances() function
library("compositions") # for ilr transform

# import arguments to run script on specific country data
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) # setup to accept arguments from command line
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied in the form of an R-script containing the following arguments: 
       1) dada2_data_path - path to directory containing input data (phyloseq object)
       2) results_path - path to results directory
       3) samp_vars - vector containing strings of variables to be analysed")
}
# load arguments into script
source(args)
# print arguments as check
print(paste0("Data path: ", dada2_data_path))
print(paste0("Results path: ", results_path))
print(paste0("Sample variables to analyse: ", samp_vars))

# load phyloseq object
print("Loading data...")
dada2 <- readRDS(paste0(dada2_data_path, "physeqob_DADA2.rds"))

# extract otu_table
abun <- abundances(dada2)
# extract sample data
sampdata <-as(sample_data(dada2),"data.frame")
sampdata <- sampdata %>%
  separate(Genus_Species, c("Genus", "Species"), " ")


###################################################################################
############################## Data transformations ###############################


print("Transforming data...")
#Subset phyloseq object to remove all samples with NA in Bd column
dada2_bd <- prune_samples(!is.na(sample_data(dada2)$Bd), dada2)
# remove any taxa with sum 0 abundance across all samples
dada2_bd_nonzero <- prune_taxa(taxa_sums(dada2_bd) > 0, dada2_bd)
# extract otu_table
abun_bd <- abundances(dada2_bd_nonzero)
# extract sample data
sampdata_bd <-as(sample_data(dada2_bd_nonzero),"data.frame")
# transform data to compositional
abun_bd_c <- transform(abun_bd, "compositional")
# account for zeros by performing replacement with small scaled numeric
detectlims_bd <- matrix(min(as.vector(abun_bd_c)[as.vector(abun_bd_c) > 0]/2), nrow(abun_bd_c), ncol(abun_bd_c))
abun0bd <- zeroreplace(abun_bd_c, d=detectlims_bd, a=0.65)
# transpose matrix
abun0t_bd <- t(abun0bd)
# transform data using isometric log-ratio transform
dada2_ilr_bd <- as.data.frame(clr(abun0t_bd))
# save transformed data for future analyses
saveRDS(dada2_ilr_bd, paste0(results_path, "dada2_ilr_bd.rds"))


## end script
