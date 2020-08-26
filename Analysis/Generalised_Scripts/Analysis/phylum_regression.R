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
library("viridis")
library("microbiome") # for abundances()
library("dplyr")
library("tidyr")

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


samps_bd_life <- samps[!is.na(samps$Lifestage),]

# subset samples by groups we want to do multivariate regression on
#samps_reg <- samps[,c("Country", "A_Family"), drop=F]
#samps_reg <- samps_reg[!rownames(samps_reg) %in% rownames(samps0), ,drop=F]

#samp_phylum <- cbind.data.frame(samps_reg, phylum_grid_t)

#samples_country <- samp_phylum %>% group_by(Country) %>% sample_n(size = 11)
#phyla_comp <- samples_country[,c(3:21)]
#phyla_comp <- as.matrix(phyla_comp)
#phyla_comp_t <- t(phyla_comp)
#countries <- samples_country[,1,drop=F]


##################################### ilr transform ######################################



library("compositions")
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
samps_perm <- samps_bd[,c("Bd", "A_Family", "Country", "A_Genus_Species")]



phylum_ilr_nona_life <- phylum_ilr[rownames(phylum_ilr) %in% rownames(samps_bd_life),]
samps_bd_life <- samps_bd_life[rownames(samps_bd_life) %in% rownames(phylum_ilr),]
samps_perm_life <- samps_bd_life[,c("Bd", "A_Family", "Country", "Lifestage")]

samps_perm_life_gs <- samps_bd_life[,c("Bd", "A_Family", "A_Genus_Species", "Country", "Lifestage")]



#samps_lm <- samps[,c("Country", "A_Family", "A_Order", "A_Genus_Species"), drop=F]
#samps_lm <- samps_lm[rownames(samps_lm) %in% rownames(phylum_ilr), ,drop=F]
#samps_bd <- samps[,c("Bd"), drop=F]
#samps_bd <- samps_bd[rownames(samps_bd) %in% rownames(phylum_ilr), ,drop=F]
# combine phylum_clr and samps_reg
#data <- cbind.data.frame(phylum_ilr, samps_lm, samps_bd)


#data_raw <- cbind.data.frame(samps_lm, samps_bd)
#data_raw <- data_raw[!is.na(data_raw$Bd),]
#phyla_raw <- phylum_grid[rownames(phylum_grid_t) %in% rownames(data_raw)]


################################################################################3#########
############################# Compositional alpha regression #############################


#require("Compositional")

#opt_alfa <- alfareg.tune(phyla_comp, countries, a = seq(0.1, 1, by = 0.1), nfolds = 10, folds = NULL, nc = 1, seed = TRUE, graph = FALSE)

#phylum_reg <- alfa.reg(phyla_comp, countries, a=1, xnew = NULL, yb = NULL, seb = FALSE)

# alfa not the issue
#alfa_trans <- alfa(phylum_mat, a=1, h=FALSE)


##################################### GLM trials ##########################################


#require(lme4)
#model_glm <- glmer(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17+V18~Country+A_Family+Bd+(1|A_Genus_Species),family=gaussian,data=data_bd)
#summary(model_glm)
#Anova(model_glm)

# try removing pyrenees
data_p <- data[!data$Country == "Pyrenees",]
#data_bd <- data[!is.na(data$Bd),]
#data_bd_fam <- data_bd[!is.na(data_bd$A_Family),]
#phylum_ilr_bd <- phylum_ilr[rownames(phylum_ilr) %in% rownames(data_bd), ,drop=F]
######## WORKS
# doesn't aacount for compositional nature of data (including depedencies between all Vs)
# or does it beause data has been ilr transformed?
#model_lmer <- lmer(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17+V18 ~ A_Family+Bd+Country+(1|A_Genus_Species)+(1|A_Order),data=data_bd)
#library(car)
#Anova(model_lmer)



#model2 <- lm(cbind(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16,V17,V18) ~ Country+A_Order+A_Family,data=data)
#Anova(model2)


######################################## Permanovas #####################################


library(vegan)
#permanova_calc_country_fam <- adonis(phylum_ilr ~ Country+A_Family, data = samps_bd_fam, permutations=999,method="euclidean")
#as.data.frame(permanova_calc_country_fam[["aov.tab"]])

#permanova_calc_fambdcountry <- adonis2(phylum_ilr_nona ~ Bd+Country+A_Family+A_Genus_Species,by="margin",data=samps_perm, permutations=999,method="euclidean")
#permanova_calc_fambdcountry_reord_nofam <- adonis2(phylum_ilr_nona ~ A_Genus_Species+Bd+Country,by="margin",data=samps_perm, permutations=999,method="euclidean")
permanova_calc_fambdcountry_nofam <- adonis2(phylum_ilr_nona ~ Bd+Country+A_Genus_Species,by="margin",data=samps_perm, permutations=999,method="euclidean")

disp_ob <- cbind.data.frame(samps_perm, phylum_ilr_nona)

# view dispersions
dist <- dist(disp_ob[,5:20], "euclidean")
disp_bd <- betadisper(dist, group=disp_ob$Bd,type="median",bias.adjust=TRUE)
disp_country <- betadisper(dist, group=disp_ob$Country,type="median",bias.adjust=TRUE)
disp_gs <- betadisper(dist, group=disp_ob$A_Genus_Species,type="median",bias.adjust=TRUE)
anova(disp_bd)
anova(disp_country)
anova(disp_gs)


# perform PCA on non-NA Bd ilr transformed data
ilr_bd_pca <- prcomp(phylum_ilr_nona)
# calculate values per sample for all axes
ilr_bd_pca_scores <- scores(ilr_bd_pca)
# extract principle axes
ilr_bd_pca_scores_sub <- ilr_bd_pca_scores[,1:2]
# extract variance explained by principle axes 
axes_bd <- summary(ilr_bd_pca)$importance[,1:2]
# add sample data 
ilr_bd_pca_final <- cbind(ilr_bd_pca_scores_sub, samps_perm)



print("Plotting...")
# plot for non-NA Bd PCA
pdf(paste0(results_path, "Summary/ilr_phylum_bd.pdf"))
ilr_bd_pca_plot <- ggplot(ilr_bd_pca_final,
                          aes(x = PC1, y= PC2)) + 
  stat_ellipse(type = "t", aes(color=Country), level = 0.95, alpha = 0.5) + 
  geom_point(aes(colour = Country), size=0.25) + 
  theme_bw() + 
  labs(fill = "Country", x = paste0("PC1 (", round(axes_bd[,1], 3), "%)"), y = paste0("PC2 (", round(axes_bd[,2], 3), "%)"))
print(ilr_bd_pca_plot)
dev.off()


# plot for all other variables
for (aspect in samp_vars) {
  pdf(paste0(results_path, "Summary/phylum_ilr_pca_", aspect, ".pdf"))
  ilr_pca_plot <- ggplot(ilr_bd_pca_final,
                         aes(x = PC1, y= PC2)) + 
    stat_ellipse(type = "t", aes(color=get(aspect)), level = 0.95, alpha = 0.5) + 
    geom_point(aes(colour = get(aspect)), size=0.25) + 
    theme_bw() + 
    labs(colour = aspect, x = paste0("PC1 (", round(axes_bd[,1], 2), "%)"), y = paste0("PC2 (", round(axes_bd[,2], 2), "%)")) 
  print(ilr_pca_plot)
  dev.off()
}

# plot for all other variables
for (aspect in samp_vars) {
  pdf(paste0(results_path, "Summary/phylum_ilr_pca_", aspect, ".pdf"))
  ilr_pca_plot <- ggplot(ilr_bd_pca_final,
                         aes(x = PC1, y= PC2)) + 
    stat_ellipse(type = "t", aes(color=get(aspect)), level = 0.95, alpha = 0.5) + 
    geom_point(aes(colour = get(aspect)), size=0.25) + 
    theme_bw() + 
    theme(legend.position = "none") +
    labs(colour = aspect, x = paste0("PC1 (", round(axes_bd[,1], 2), "%)"), y = paste0("PC2 (", round(axes_bd[,2], 2), "%)")) 
  print(ilr_pca_plot)
  dev.off()
}

pdf(paste0(results_path, "Summary/phylum_ilr_pca_genusspecies.pdf"))
ilr_pca_plot <- ggplot(ilr_bd_pca_final,
                       aes(x = PC1, y= PC2)) + 
  stat_ellipse(type = "t", aes(color=A_Genus_Species), level = 0.95, alpha = 0.5) + 
  geom_point(aes(colour = A_Genus_Species), size=0.25) + 
  theme_bw() + 
  theme(legend.position = "none") +
  labs(colour = aspect, x = paste0("PC1 (", round(axes_bd[,1], 2), "%)"), y = paste0("PC2 (", round(axes_bd[,2], 2), "%)")) 
print(ilr_pca_plot)
dev.off()



sum(phylum_grid_t$p__Ascomycota)/nrow(phylum_grid_t)
sum(phylum_grid_t$p__Basidiomycota)/nrow(phylum_grid_t)
pyr <- samps[samps$Country == "Pyrenees",]
pyr_phy <- phylum_grid_t[rownames(phylum_grid_t)%in% pyr$MiSeqCode,]
View(pyr_phy)
sum(pyr_phy$p__Chytridiomycota)/nrow(pyr_phy)
sum(phylum_grid_t$p__Chytridiomycota)/nrow(phylum_grid_t)

sum(pyr_phy$p__Rozellomycota)/nrow(pyr_phy)
sum(phylum_grid_t$p__Rozellomycota)/nrow(phylum_grid_t)


chyt <- taxa[!is.na(taxa$Phylum),]
chyt <- chyt[chyt$Phylum == "p__Chytridiomycota",]
chyt$Genus <- as.character(chyt$Genus)
chyt$Genus[is.na(chyt$Genus)] <- "NA"
table(chyt$Genus)

bdpos <- samps[samps$Bd == 1,]
bdpos_phy <- phylum_grid_t[rownames(phylum_grid_t)%in% bdpos$MiSeqCode,]
sum(bdpos_phy$p__Chytridiomycota)/nrow(bdpos_phy)
sum(phylum_grid_t$p__Chytridiomycota)/nrow(phylum_grid_t)


pyr_seq <- seqs[rownames(seqs)%in% pyr$MiSeqCode,]
pyr_seq_t <- as.data.frame(t(pyr_seq))
pyr_seq_t$sum <- rowSums(pyr_seq_t)
pyr_seq_t <- pyr_seq_t[pyr_seq_t$sum!=0,]
pyr_seq_t <- pyr_seq_t[,-ncol(pyr_seq_t)]
pyr_seq <- as.data.frame(t(pyr_seq_t))
chyt_pyr <- chyt[rownames(chyt) %in% colnames(pyr_seq),]
table(chyt_pyr$Genus)                                                   

