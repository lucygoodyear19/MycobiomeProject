##############################################################################
###################### DADA2 and Esto Results Analysis #######################
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
library(dplyr)
library(phyloseq)
library(Biostrings)
library(ggtree)
library(scales)
#library(nlme)
library(lme4)
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


#################################################################################
############################## Co-occurence #####################################


library(cooccur)

# add Bd binary column to otu_table
cooccur_object <- cbind(as.data.frame(otu_table(dada2)), t(sample_data(dada2)$Bd))
# set Bd column to numeric
cooccur_object[,ncol(cooccur_object)] <- sapply(cooccur_object[,ncol(cooccur_object)], as.character)
cooccur_object[,ncol(cooccur_object)] <- sapply(cooccur_object[,ncol(cooccur_object)], as.numeric)
# set all abundances > 0 to 1 to get a presence absence table
cooccur_object[cooccur_object > 0] <- 1
# main coocur function to calculate pairwise cooccurence patterns
cooccur_full <- cooccur(cooccur_object, 
                        type="site_spp",
                        thresh=TRUE, 
                        spp_names=TRUE)

# look at cooccurence patterns just for Bd
#cooccur_bd <- pair(cooccur_full, "s__dendrobatidis")

# view all significant results for bd
#sig_results_bd <- print.cooccur(cooccur_bd)

# plot as a heatmap
#heatmap_bd <- plot.cooccur(cooccur_bd)


###################################################################################
################################# Beta-diversity ##################################


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


###############################################################################
############################## Kruskall-Wallis ################################


# load packages
require(microbiomeSeq)

# run Kruskal_Wallis test on Bd for dada2
(kw_dada2 <- kruskal_abundance(dada2, group = "Bd", pvalue.threshold = 0.05))
plot_signif(kw_dada2$plotdata)
# run Kruskal_Wallis test on Bd for dada2
(kw_esto <- kruskal_abundance(esto, group = "Bd", pvalue.threshold = 0.05))
plot_signif(kw_dada2$plotdata)

# create a list of the significant taxa
sig_asv <- kw_dada2$importance[,1]
sig_tax <- tax_table(dada2)[rownames(tax_table(dada2)) %in% sig_asv,]

# run Kruskal_Wallis using base R stats
kruskal.test(alpha ~ Bd, data = alpha_prep(dada2, "Bd"))
kruskal.test(alpha ~ Bd, data = alpha_prep(esto, "Bd"))

# detach packages
detach("package:microbiomeSeq", unload=TRUE)

        