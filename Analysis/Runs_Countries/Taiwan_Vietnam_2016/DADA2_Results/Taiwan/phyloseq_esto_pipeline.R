##############################################################################
###################### Creation of Phyloseq object ###########################
##############################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


################################## Set up ####################################


# set wd if needed
setwd("~/Documents/MRes/MycobiomeProject/Analysis/Runs_Countries/Taiwan_Vietnam_2016/Taiwan")

# load packages
library(dplyr)
library(phyloseq)
library(ggplot2)

# load filtered Taiwan data (swabs only)
TWdata <- read.csv("Filtering/Taiwanswabsonly_plate2-4_filteredOTUtable.csv", 
                   header=T, 
                   stringsAsFactor=FALSE, 
                   check.names = F)

# load plate data
plates <- read.csv("Data/TW16_plate_data.csv", 
                   stringsAsFactor=F, 
                   header=F)
names(plates) <- c("barcode","sample","DNAqual","plate")

# load metadata
metadata <- read.csv("Data/TW16_metadata.csv",
                     stringsAsFactor=F, 
                     header=T)
rownames(metadata) <- metadata$miseq_id
metadata <- metadata[,-1]

# add column for Bd +ve/-ve to metadata
metadata$Bd <- 0
for (i in (1:nrow(metadata))){
  if (metadata$qPCR[i] > 0){
    metadata$Bd[i] <- 1
  }
}
# set as factor for plotting
metadata$Bd <- as.factor(metadata$Bd)


###################### Create phyloseq object of all data #####################


# subset all info and save as separate dataframe
info <- TWdata[,1:22] 

# subset samples only and save as matrix
samp_only <- as.matrix(TWdata[,23:ncol(TWdata)])
# add OTU ids as row names
rownames(samp_only) <- TWdata$otuid

# sum all read columns and save as vector
sampsum <- colSums(samp_only)
# filter by sums not equal to zero
samp_only <- samp_only[,sampsum > 10]

# subset taxa only and save as matrix
tax_data <- as.matrix(TWdata[,10:16])
# add OTU ids as row names
rownames(tax_data) <- TWdata$otuid  

# convert sample and taxa matrices into the tables needed
# to create a phyloseq object
OTU <- otu_table(samp_only, taxa_are_rows = TRUE)
TAX <- tax_table(tax_data)
Sample <- sample_data(metadata)
# combine tables into one phyloseq object
physeq <- phyloseq(OTU,TAX, Sample)


############## Create phyloseq object of phyla subsetted data ##########


# create another subset with certain phyla removed
# remove all OTUs of the followng phlya:
# - "" (blank)
# - "environmental samples"
phy_fil <- TWdata %>% 
  filter(phylum != "" 
         & phylum != "environmental samples")
# subset samples only and save as matrix
phy_fil_samp <- as.matrix(phy_fil[,23:ncol(phy_fil)])
# add OTU ids as row names
rownames(phy_fil_samp) <- phy_fil$otuid

# subset all info and save as separate dataframe
info_fil <- phy_fil[,1:22] 

# subset taxa only and save as matrix
phy_fil_tax <- as.matrix(phy_fil[,10:16])
# add OTU ids as row names
rownames(phy_fil_tax) <- phy_fil$otuid  

# check to see if there samples with no OTUs
# sum all read columns and save as vector
#sampsum <- colSums(phy_fil_samp)
# filter by sums not equal to zero
#phy_fil_samp1 <- phy_fil_samp[,sampsum != 0]

# filter metdata by IDs present in current MiSeq run
metadata_fil <- metadata[rownames(metadata) %in% colnames(phy_fil_samp),]

# convert sample and taxa matrices into the tables needed
# to create a phyloseq object
OTU_fil <- otu_table(phy_fil_samp, taxa_are_rows = TRUE)
TAX_fil <- tax_table(phy_fil_tax)
Sample_fil <- sample_data(metadata_fil)
# combine tables into one phyloseq object
physeq_fil <- phyloseq(OTU_fil,TAX_fil, Sample_fil)


############################# Filter data ###############################


#plot_bar(physeq_fil, fill = "phylum")

wh0 = genefilter_sample(physeq_fil, filterfun_sample(function(x) x > 1), A=0.05*nsamples(physeq_fil))
TW_physeq_filtered = prune_taxa(wh0, physeq_fil)

#phylum.sum = tapply(taxa_sums(TW_physeq_filtered), tax_table(TW_physeq_filtered)[, "Phylum"], sum, na.rm=TRUE)
#top5phyla = names(sort(phylum.sum, TRUE))[1:5]
#TW_physeq_filtered = prune_taxa((tax_table(TW_physeq_filtered)[, "Phylum"] %in% top5phyla), TW_physeq_filtered)

# save phyloseq object to be imported into analysis scripts
saveRDS(physeq, paste0("physeqob_esto.rds"))

## end of script
