##############################################################################
###################### Prelimary Taiwan Data Analysis ########################
##############################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


################################## Set up ####################################


# set wd if needed
setwd("~/Documents/CMEECourseWork/Project/Data_Scripts/Taiwan_Vietnam/")

# load packages
library(dplyr)
library(phyloseq)
library(ggplot2)

# load filtered Taiwan data (swbas only)
TWdata <- read.csv("Filtering/Taiwanswabsonly_plate2-4_filteredOTUtable.csv", 
                   header=T, 
                   stringsAsFactor=FALSE, 
                   check.names = F)

# load plate data
plates <- read.csv("TW16_plate_data.csv", 
                   stringsAsFactor=F, 
                   header=F)
names(plates) <- c("barcode","sample","DNAqual","plate")

# load metadata
metadata <- read.csv("TW16_metadata.csv",
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
# - "Fungi_phy_Incertae_sedis"
phy_fil <- TWdata %>% 
  filter(phylum != "" 
         & phylum != "environmental samples" 
         & phylum != "Fungi_phy_Incertae_sedis")
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


######################## Plot alpha-diversity ############################


# estimate value for alpha for each sample
est_alpha_shannon <- estimate_richness(physeq_fil, split = TRUE, measures = "Shannon")
metadata_fil$alpha <- est_alpha_shannon[,1]

########## Compare pair-wise and plot with ggplot

##### Lifestage

# subset Sample_fil twice to include only shannon diversity,
# split by lifestage
metadata_fil_alpha_life <- as.data.frame(metadata_fil %>% select(alpha, Lifestage))

# perform Mann-Whitney U test
wilcox.test(alpha ~ Lifestage, data = metadata_fil_alpha_life)

# plot box plot
p_lifestage_box <- ggplot(metadata_fil, 
       aes(x = Lifestage, y = alpha)) + 
  labs(y = "Shannon Alpha Diversity") +
  geom_boxplot() +
  scale_x_discrete(labels=c(paste("Adult (n=", table(metadata_fil_alpha_life$Lifestage)[1], ")"),
                            paste("Tadpole (n=", table(metadata_fil_alpha_life$Lifestage)[2], ")")))

pdf("Lifestage_box.pdf")
p_lifestage_box
dev.off()

# plot box plot with Bd comparison
ggplot(metadata_fil, aes(x = Lifestage, y = alpha, colour = Bd)) + 
  geom_boxplot()

# plot a violin plot with Bd comparison
ggplot(metadata_fil, aes(x = Lifestage, y = alpha, colour = Bd)) + 
  geom_violin()

##### Altitude_m

# subset Sample_fil twice to include only shannon diversity,
# split by altitude
metadata_fil_alpha_alt <- as.data.frame(metadata_fil %>% select(alpha, Altitude_m))

# perform Mann-Whitney U test
wilcox.test(alpha ~ Altitude_m, data = metadata_fil_alpha_alt)

# plot box plot
p_altitude_box <- ggplot(metadata_fil, 
                          aes(x = Altitude_m, y = alpha)) + 
  labs(y = "Shannon Alpha Diversity") +
  geom_boxplot() +
  #scale_x_discrete(labels=c(paste("Adult (n=", table(metadata_fil_alpha$Lifestage)[1], ")"),
  #                          paste("Tadpole (n=", table(metadata_fil_alpha$Lifestage)[2], ")")))

pdf("Altitude_box.pdf")
p_altitude_box
dev.off()

# plot box plot
ggplot(metadata_fil, aes(x = as.factor(Altitude_m), y = alpha, colour = Bd)) + 
  geom_boxplot()

# plot a violin plot
ggplot(metadata_fil, aes(x = as.factor(Altitude_m), y = alpha, colour = Bd)) + 
  geom_violin()

##### Latitude

# plot box plot
ggplot(metadata_fil, aes(x = as.factor(Latitude), y = alpha, colour = Bd)) + 
  geom_boxplot()

# plot a violin plot
ggplot(metadata_fil, aes(x = as.factor(Latitude), y = alpha, colour = Bd)) + 
  geom_violin()

##### Longitude (same as Latitude)

# plot box plot
ggplot(metadata_fil, aes(x = as.factor(Longitude), y = alpha, colour = Bd)) + 
  geom_boxplot()

# plot a violin plot
ggplot(metadata_fil, aes(x = as.factor(Longitude), y = alpha, colour = Bd)) + 
  geom_violin()

##### Species

# plot box plot
p_species_box <- ggplot(metadata_fil, aes(x = Species, y = alpha)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf("Species_box.pdf")
p_species_box
dev.off()

# plot a violin plot
ggplot(metadata_fil, aes(x = Species, y = alpha, colour = Bd)) + 
  geom_violin()

##### Bd +ve/-ve

# subset Sample_fil twice to include only shannon diversity,
# split by Bd +ve/-ve
metadata_fil_alpha_bd <- as.data.frame(metadata_fil %>% select(alpha, Bd))

# perform Mann-Whitney U test
wilcox.test(alpha ~ Bd, data = metadata_fil_alpha_bd)

# plot box plot
p_bd_box <- ggplot(metadata_fil, 
                         aes(x = Bd, y = alpha)) + 
  labs(y = "Shannon Alpha Diversity") +
  geom_boxplot()
  #scale_x_discrete(labels=c(paste("Adult (n=", table(metadata_fil_alpha$Lifestage)[1], ")"),
  #                          paste("Tadpole (n=", table(metadata_fil_alpha$Lifestage)[2], ")")))
  
pdf("Bd_box.pdf")
p_bd_box
dev.off()

# plot box plot
ggplot(metadata_fil, aes(x = Bd, y = alpha)) + 
  geom_boxplot()

# plot a violin plot
ggplot(metadata_fil, aes(x = Bd, y = alpha)) + 
  geom_violin()

########## using phyloseq

TW_physeq_fil_non0 <- prune_taxa(taxa_sums(physeq_fil) > 0, physeq_fil)

plot_richness(TW_physeq_fil_non0, x = "Bd", color = "Species", measures = "shannon")


######################### plot beta-diversity ###########################


physeq.ord <- ordinate(physeq, "NMDS", "bray")
plot_ordination(physeq, physeq.ord, type="samples", color="Bd", title="taxa")


# vegan trials
# TW_NMDS_k2 <- metaMDS(samp_only_trans, k=2) # ---> stress = 0.2304577
TW_NMDS_k3 <- metaMDS(samp_only_trans, k=3) # ---> stress = 0.1933014
# TW_NMDS <- metaMDS(samp_only_trans, k=4) ---> stress = 0.1595641
TW_NMDS <- metaMDS(samp_only_trans, k=5) # ---> stress = 0.1360477, 0.1359906, 0.1359468
# TW_NMDS <- metaMDS(samp_only_trans, k=5, trymax = 1000) # ---> stress = 0.1359034
TW_NMDS_k12 <- metaMDS(samp_only_trans, k=12) # ---> stress = 0.07332787
ordiplot(TW_NMDS_k12, type = "n")
orditorp(TW_NMDS_k12, display="species", col="red", air=0.01)
orditorp(TW_NMDS_k12, display="sites", cex=0.5, air=0.01)

# load metadata
plates <- read.csv("2colformatSampleSheet_additional.csv", stringsAsFactor=F, header=F)
names(plates) <- c("barcode","sample","DNAqual","plate")
environment <- read.csv("dataspreadsheet_Taiwan.csv", stringsAsFactor=F, header=F)

# split data by plate and do NMDS
plate2 <- TWdata[,names(TWdata) %in% c(plates$sample[plates$plate == "plate2"])]
samp_only_trans_plate2 <- t(plate2[,23:ncol(plate2)])
colnames(samp_only_trans_plate2) <- plate2$otuid

# remove OTUs without reads
# sum all read columns and save as vector
sampsum <- colSums(samp_only_trans_plate2)
# filter by sums not equal to zero
samp_only_trans_plate2 <- samp_only_trans_plate2[,sampsum != 0]

# NMDS
plate2_NMDS_k3 <- metaMDS(samp_only_trans, k=3) # ---> stress = 0.1933076
ordiplot(plate2_NMDS_k3, type = "n")
orditorp(plate2_NMDS_k3, display="species", col="red", air=0.01)
orditorp(plate2_NMDS_k3, display="sites", cex=0.5, air=0.01)


########## Kruskall-Wallis

# first add column for Bd +ve/-ve
metadata_fil$Bd <- 0
for (i in (1:nrow(metadata_fil))){
  if (metadata_fil$qPCR[i] > 0){
      metadata_fil$Bd[i] <- 1
  }
}

# load packages
require(TBrach)


########################### Kruskall-Wallis #############################


test_differential_abundance_Kruskal(physeq, group = "Bd", compare = NULL,
                                    block = NULL, p.adjust.method = "fdr")



        