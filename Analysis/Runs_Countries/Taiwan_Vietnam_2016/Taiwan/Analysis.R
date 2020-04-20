##############################################################################
###################### Prelimary Taiwan Data Analysis ########################
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
library(gridExtra)

# load filtered phyloseq object for Jen's data
esto <- readRDS("physeqob_jen.rds")

# create phyloseq object for DADA2
tax <- read.table("HPC_Results/Previous_Run/Taxa_Table.txt", 
                  stringsAsFactors = F)
seqtab <- read.table("HPC_Results/Previous_Run/Seq_Abun_Table.txt", 
                     stringsAsFactors = F, 
                     header = T)
rownames(seqtab) <- gsub("-", "", rownames(seqtab))
plates <- read.csv("Data/TW16_plate_data.csv", 
                   stringsAsFactor=F, 
                   header=F)
names(plates) <- c("barcode","sample","DNAqual","plate")
plates_fil <- plates[-c(3,4)]
trial <- merge(plates_fil, seqtab, by.x="barcode", by.y = "row.names")

rownames(trial) <- trial$sample
trial <- trial[-c(1,2)]

metadata <- read.csv("metadata_bd.csv",
                     stringsAsFactor=F, 
                     header=T)
rownames(metadata) <- metadata$X
metadata <- metadata[,-1]
metadata$Bd <- as.factor(metadata$Bd)

rownames(tax) <- colnames(seqtab)

tax_mat <- as.matrix(tax) # required to create phyloseq object
trial_mat <- as.matrix(trial)

dada2 <- phyloseq(tax_table(tax_mat), 
                  otu_table(trial_mat, taxa_are_rows = FALSE), 
                  sample_data(metadata))


######################## Plot alpha-diversity ############################


# estimate value for alpha for each sample
shannon_esto <- estimate_richness(esto, split = TRUE, measures = "Shannon")
sample_data(esto)$alpha <- as.numeric(shannon_esto[,1])

shannon_dada2 <- estimate_richness(dada2, split = TRUE, measures = "Shannon")
sample_data(dada2)$alpha <- as.numeric(shannon_dada2[,1])


########## Compare pair-wise and plot with ggplot

###################### Lifestage ##########################

# subset to include only shannon diversity, split by lifestage

alpha_life_esto <- as.data.frame(sample_data(esto) %>% select(alpha, Lifestage))

alpha_life_dada2 <- as.data.frame(sample_data(dada2) %>% select(alpha, Lifestage))

alpha_life_esto <- as.data.frame(as.matrix(alpha_life_esto))
alpha_life_esto$alpha <- sapply(alpha_life_esto[,1], as.double)
alpha_life_esto$Lifestage <- as.factor(alpha_life_esto$Lifestage)

alpha_life_dada2 <- as.data.frame(as.matrix(alpha_life_dada2))
alpha_life_dada2$alpha <- sapply(alpha_life_dada2[,1], as.double)
alpha_life_dada2$Lifestage <- as.factor(alpha_life_dada2$Lifestage)

# perform Mann-Whitney U test on esto
wilcox.test(alpha ~ Lifestage, data = alpha_life_esto)
require("effsize")
cohen.d(alpha ~ Lifestage, data = alpha_life_esto)

# perform Mann-Whitney U test on dada2
wilcox.test(alpha ~ Lifestage, data = alpha_life_dada2)
require("effsize")
cohen.d(alpha ~ Lifestage, data = alpha_life_dada2)

# plot box plot
p_lifestage_box_esto <- ggplot(sample_data(esto), 
       aes(x = Lifestage, y = alpha)) + 
  labs(y = "Shannon Alpha Diversity - Estonian Pipeline") +
  geom_boxplot() +
  scale_x_discrete(labels=c(paste("Adult (n=", table(alpha_life_esto$Lifestage)[1], ")"),
                            paste("Tadpole (n=", table(alpha_life_esto$Lifestage)[2], ")"))) +
  theme_bw()

p_lifestage_box_dada2 <- ggplot(sample_data(dada2), 
                               aes(x = Lifestage, y = alpha)) + 
  labs(y = "Shannon Alpha Diversity - DADA2 Pipeline") +
  geom_boxplot() +
  scale_x_discrete(labels=c(paste("Adult (n=", table(alpha_life_dada2$Lifestage)[1], ")"),
                            paste("Tadpole (n=", table(alpha_life_dada2$Lifestage)[2], ")"))) +
  theme_bw()

pdf("Lifestage_box.pdf")
grid.arrange(p_lifestage_box_esto,p_lifestage_box_dada2, ncol=2)
dev.off()

# plot box plot with Bd comparison
p_lifestage_box_bd_esto <- ggplot(sample_data(esto), 
                               aes(x = Lifestage, y = alpha, colour = Bd)) + 
  labs(y = "Shannon Alpha Diversity - Estonian Pipeline") +
  geom_boxplot() +
  scale_x_discrete(labels=c(paste("Adult (n=", table(alpha_life_esto$Lifestage)[1], ")"),
                            paste("Tadpole (n=", table(alpha_life_esto$Lifestage)[2], ")"))) +
  theme_bw()

p_lifestage_box_bd_dada2 <- ggplot(sample_data(dada2), 
                                aes(x = Lifestage, y = alpha, colour = Bd)) + 
  labs(y = "Shannon Alpha Diversity - DADA2 Pipeline") +
  geom_boxplot() +
  scale_x_discrete(labels=c(paste("Adult (n=", table(alpha_life_dada2$Lifestage)[1], ")"),
                            paste("Tadpole (n=", table(alpha_life_dada2$Lifestage)[2], ")"))) +
  theme_bw()

pdf("Lifestage_box_Bd.pdf")
grid.arrange(p_lifestage_box_bd_esto, p_lifestage_box_bd_dada2, ncol=2)
dev.off()

# plot a violin plot with Bd comparison
ggplot(metadata_fil, aes(x = Lifestage, y = alpha, colour = Bd)) + 
  geom_violin()

##### Altitude_m

# subset Sample_fil twice to include only shannon diversity,
# split by altitude
metadata_fil_alpha_alt <- as.data.frame(metadata_fil %>% select(alpha, Altitude_m))

# perform Mann-Whitney U test to look at statistical significance
wilcox.test(alpha ~ Altitude_m, data = metadata_fil_alpha_alt)
# calculate effect size
cohen.d(alpha ~ Altitude_m, data = metadata_fil_alpha_alt)

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


############################ Bd +ve/-ve #############################

# subset to include only shannon diversity, split by lifestage

alpha_bd_esto <- as.data.frame(sample_data(esto) %>% select(alpha, Bd))

alpha_bd_dada2 <- as.data.frame(sample_data(dada2) %>% select(alpha, Bd))

alpha_bd_esto <- as.data.frame(as.matrix(alpha_bd_esto))
alpha_bd_esto$alpha <- sapply(alpha_bd_esto[,1], as.double)
alpha_bd_esto$Bd <- as.factor(alpha_bd_esto$Bd)

alpha_bd_dada2 <- as.data.frame(as.matrix(alpha_bd_dada2))
alpha_bd_dada2$alpha <- sapply(alpha_bd_dada2[,1], as.double)
alpha_bd_dada2$Bd <- as.factor(alpha_bd_dada2$Bd)

# perform Mann-Whitney U test on esto
wilcox.test(alpha ~ Bd, data = alpha_bd_esto)
require("effsize")
cohen.d(alpha ~ Bd, data = alpha_bd_esto)

# perform Mann-Whitney U test on dada2
wilcox.test(alpha ~ Bd, data = alpha_bd_dada2)
require("effsize")
cohen.d(alpha ~ Bd, data = alpha_bd_dada2)

# plot box plot
p_bd_box_esto <- ggplot(sample_data(esto), 
                               aes(x = Bd, y = alpha)) + 
  labs(y = "Shannon Alpha Diversity - Estonian Pipeline") +
  geom_boxplot() +
  scale_x_discrete(labels=c(paste("Bd positive (n=", table(alpha_bd_esto$Bd)[1], ")"),
                            paste("Bd negative (n=", table(alpha_bd_esto$Bd)[2], ")"))) +
  theme_bw()

p_bd_box_dada2 <- ggplot(sample_data(dada2), 
                                aes(x = Bd, y = alpha)) + 
  labs(y = "Shannon Alpha Diversity - DADA2 Pipeline") +
  geom_boxplot() +
  scale_x_discrete(labels=c(paste("Bd positive (n=", table(alpha_bd_dada2$Bd)[1], ")"),
                            paste("Bd negative (n=", table(alpha_bd_dada2$Bd)[2], ")"))) +
  theme_bw()

pdf("Bd_box.pdf")
grid.arrange(p_bd_box_esto, p_bd_box_dada2, ncol=2)
dev.off()

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



        