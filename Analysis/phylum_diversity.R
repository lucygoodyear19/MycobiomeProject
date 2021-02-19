################################################################################################
####################################### Summary Plots ##########################################
################################################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


##############################################################################################
######################################### Set up #############################################


# load packages
library("phyloseq")
library("ggplot2")
theme_set(theme_bw())
library("microbiome") # for abundances()
library("dplyr")
library("tidyr")
library("vegan")
library("compositions")

# set wd
setwd("/Users/lucy/Documents/MRes/MycobiomeProject/Analysis/Results/Global_14/")

# load phyloseq object
dada2 <- readRDS(paste0(dada2_data_path, "physeqob_DADA2.rds"))

# extract sample data from phyloseq object
samps <- as(sample_data(dada2), "data.frame")
# extract taxa table as dataframe from phyloseq object
taxa <- as.data.frame(tax_table(dada2))
#taxa <- as.data.frame(lapply(taxa, function(x) gsub(".__", "", x)))
# extract ASV abundance table as dataframe from phyloseq object
seqs <- as.data.frame(t(abundances(dada2)))
# transpose abundance data frame
seqs_t <- as.data.frame(t(seqs))


########################################################################################
################################### Stats prep #########################################


#################################### Compositional matrix #####################################


# Add phylum to abundance data
perm_df <- seqs_t
perm_df$Phylum <- taxa$Phylum
# merge by samples by phylum
phylum_grid <- plyr::ddply(perm_df, "Phylum", plyr::numcolwise(sum))
# remove any NAs
phylum_grid <- phylum_grid[!is.na(phylum_grid$Phylum),]
# set rownames to phyla
rownames(phylum_grid) <- phylum_grid$Phylum
# remove phyla column
phylum_grid <- phylum_grid[,-1]
# sum each phylum
phylum_grid$sum <- rowSums(phylum_grid)
# remove any phyla with a sum equal to zero
phyla0 <- phylum_grid[phylum_grid$sum == 0,]
phylum_grid <- phylum_grid[!rownames(phylum_grid) %in% rownames(phyla0),]
# remove sum column
phylum_grid <- phylum_grid[,-ncol(phylum_grid)]
# set transposed abundances as data frame
phylum_grid_t <- as.data.frame(t(phylum_grid))
# sum each sample
phylum_grid_t$sum <- rowSums(phylum_grid_t)
# remove any samples with a sum equal to 0
samps0 <- phylum_grid_t[phylum_grid_t$sum == 0,]
phylum_grid_t <- phylum_grid_t[!rownames(phylum_grid_t) %in% rownames(samps0),]
# get compositional matrix by diviving by sums
phylum_grid_t <- phylum_grid_t/phylum_grid_t$sum
# remove sum column
phylum_grid_t <- phylum_grid_t[,-ncol(phylum_grid_t)]
# transform to matirx and transpose
phylum_mat <- as.matrix(phylum_grid_t)
phylum_mat_t <- t(phylum_mat)

# remove all NA entries for Bd and amphibian Family
samps_bd <- samps[!is.na(samps$Bd),]


##################################### ilr transform ######################################


# account for zeros by performing replacement with small scaled numeric
detectlims_phylum <- matrix(min(as.vector(phylum_mat)[as.vector(phylum_mat) > 0]/2), nrow(phylum_mat), ncol(phylum_mat))
phylum_mat0 <- zeroreplace(phylum_mat, detectlims_phylum, a=0.65) 
# transpose matrix
phylum_mat0t <- t(phylum_mat0)
# transform data using isometric log-ratio transform
#phylum_clr <- as.data.frame(clr(phylum_mat0))
phylum_ilr <- as.data.frame(ilr(phylum_mat0))
# subset to non-NA samples for Bd and amphibian Family
phylum_ilr_nona <- phylum_ilr[rownames(phylum_ilr) %in% rownames(samps_bd),]
samps_bd <- samps_bd[rownames(samps_bd) %in% rownames(phylum_ilr),]
# Just keep the necessary columns
samps_perm <- samps_bd[,c("Bd", "A_Family", "Country")]


########################################################################################
#################################### PERMANOVA #########################################


print("Performing PERMANOVA...")

# perform PERMANOVA on 3 variables
permanova <- adonis2(phylum_ilr_nona ~ get(samp_vars_1host[1])+get(samp_vars_1host[2])+get(samp_vars_1host[3]), 
                     by = "margin", data = data=samps_perm, permutations=999,method="euclidean")
# save to text file
cat("\n\nPERMANOVA: ", capture.output(permanova_out),
    file=paste0(path_out, "permanovas_phyla.txt"), sep = "\n", append=TRUE)

# calculate dispersion significance for the 3 variables
disp_ob <- cbind.data.frame(samps_perm, phylum_ilr_nona)

# create distance matrix
dist <- dist(disp_ob[,5:20], "euclidean")
# calculate dispersion for the 3 variables
for (var in 1:length(samp_vars_1host)){
  dispersion <- betadisper(dist, group=disp_ob[[samp_vars_1host[var]]],type="median",bias.adjust=TRUE)
  # extract ANOVA
  anova_disp <- anova(dispersion)
  # save dispersion ANOVAs to text file
  cat(paste0("\n\nANOVA for dispersion of ", samp_vars_1host[var], ": "), capture.output(anova_disp),
      file=paste0(path_out, "permanovas_phyla.txt"), sep = "\n", append=TRUE)
}


########################################################################################
################################## Inferences and numbers ##############################


# calculate percentage of reads of different phyla
print("Percentage of reads assigned to Ascomycota:")
sum(phylum_grid_t$p__Ascomycota)/nrow(phylum_grid_t)
print("Percentage of reads assigned to Basidiomycota:")
sum(phylum_grid_t$p__Basidiomycota)/nrow(phylum_grid_t)
print("Percentage of reads assigned to Chytriodiomycota:")
sum(phylum_grid_t$p__Chytridiomycota)/nrow(phylum_grid_t)
print("Percentage of reads assigned to Rozellomycota:")
sum(phylum_grid_t$p__Rozellomycota)/nrow(phylum_grid_t)

# view break down of chytrid genera
chyt <- taxa[!is.na(taxa$Phylum),]
chyt <- chyt[chyt$Phylum == "p__Chytridiomycota",]
chyt$Genus <- as.character(chyt$Genus)
chyt$Genus[is.na(chyt$Genus)] <- "NA"
print("Table of chytrid genera:")
table(chyt$Genus)

# find how many chytrid reads were found in Bd positive samples
bdpos <- samps[samps$Bd == 1,]
bdpos_phy <- phylum_grid_t[rownames(phylum_grid_t)%in% bdpos$MiSeqCode,]
print("Percentage of Bd positive sample reads assigned to Chytyridiomycota:")
sum(bdpos_phy$p__Chytridiomycota)/nrow(bdpos_phy)

# subset by Pyrenees
pyr <- samps[samps$Country == "Pyrenees",]
pyr_phy <- phylum_grid_t[rownames(phylum_grid_t)%in% pyr$MiSeqCode,]
# calculate percentage of reads of different phyla in Pyrenees
print("Percentage of reads in Pyrenees assigned to Chytriodiomycota:")
sum(pyr_phy$p__Chytridiomycota)/nrow(pyr_phy)
print("Percentage of reads in Pyrenees assigned to Rozellomycota:")
sum(pyr_phy$p__Rozellomycota)/nrow(pyr_phy)
# find break down of chytrid genera in Pyrenees
pyr_seq <- seqs[rownames(seqs)%in% pyr$MiSeqCode,]
pyr_seq_t <- as.data.frame(t(pyr_seq))
pyr_seq_t$sum <- rowSums(pyr_seq_t)
pyr_seq_t <- pyr_seq_t[pyr_seq_t$sum!=0,]
pyr_seq_t <- pyr_seq_t[,-ncol(pyr_seq_t)]
pyr_seq <- as.data.frame(t(pyr_seq_t))
chyt_pyr <- chyt[rownames(chyt) %in% colnames(pyr_seq),]
print("Table of chytrid genera in Pyrenees:")
table(chyt_pyr$Genus)


## end of script
