##############################################################################
################################ Beta Diversity ##############################
##############################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


##############################################################################
################################## Set up ####################################


# load packages
library("dplyr")
library("phyloseq")
library("compositions")
library("microbiome")
library("ggplot2")
library("vegan")
library("compositions")
library("tidyr") # for separate function

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
dada2 <- readRDS(paste0(dada2_data_path, "physeqob_DADA2.rds"))

# to account for Costa Rica singles
for (samp in 1:nrow(sample_data(dada2))) {
  if (!is.na(sample_data(dada2)$Bd[samp]) & sample_data(dada2)$Bd[samp] == "S") {
    sample_data(dada2)$Bd[samp] <- 0
  }
}
sample_data(dada2)$Bd <- as.factor(as.character(sample_data(dada2)$Bd))

# extract otu_table
abun <- abundances(dada2)
# extract sample data
sampdata <-as(sample_data(dada2),"data.frame")
sampdata <- sampdata %>%
  separate(Genus_Species, c("Genus", "Species"), " ")


###################################################################################
############################## Data transformations ###############################


# transform data to compositional
abun_c <- transform(abun, "compositional")
# account for zeros by performing replacement with small scaled numeric
detectlims <- matrix(min(as.vector(abun_c)[as.vector(abun_c) > 0]/2), nrow(abun_c), ncol(abun_c))
abun_0 <- zeroreplace(abun_c, d=detectlims, a=0.65)
# transpose matrix
abun0t <- t(abun_0)
#saveRDS(abun0t, paste0(results_path, "BetaDiversity/abun0t_africa.rds"))


####### Bd
#Subset ECU to Complete Cases 
dada2_bd <- prune_samples(!is.na(sample_data(dada2)$Bd), dada2)
#Can't Have Zero-Sum taxa for distances 
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

dada2_ilr_bd <- as.data.frame(clr(abun0t_bd))

dada2_ilr <- readRDS(paste0(dada2_data_path, "BetaDiversity/dada2_ilr_SouthAm.rds"))

# ilr transform
#dada2_ilr <- as.data.frame(ilr(abun0t))
# clr tranform using compositions package
#dada2_clr_comp <- as.data.frame(clr(abun0t))

#dada2_clr_comp_phylo <- phyloseq(otu_table(dada2_clr_comp, taxa_are_rows = F),
#                            tax_table(dada2),
#                            sample_data(dada2),
#                            refseq(dada2))

# clr transform using microbiome package (contains zero replacement function)
#dada2_clr_micro_phylo <- microbiome::transform(dada2, "clr")
#dada2_clr_micro <- as.data.frame(otu_table(dada2_clr_micro_phylo))


###################################################################################
############################## Using in-built phyloseq ############################


########## non-metric multidimensional scaling (NMDS) on raw data

#physeq_ord_nmds <- ordinate(dada2, "NMDS", "bray")
#pdf(paste0(results_path, "BetaDiversity/NMDS_phyloseq_country.pdf"))
#print(plot_ordination(dada2, physeq_ord_nmds, type="samples", color="Country", title=""))
#dev.off()

########## redundancy analysis (PCoA) on transformed data

# redundancy analysis on transformed data
#physeq_ord_rda_clr <- ordinate(dada2_clr_comp_phylo, "RDA")
#pdf(paste0(results_path, "BetaDiversity/phyloseq_comp_clr_rda_country.pdf"))
#plot_ordination(dada2, physeq_ord_rda_clr, type="samples", color="Country", title="") +
#  stat_ellipse(type = "t", aes(color=Country), level = 0.95, alpha = 0.5) + 
#  theme_bw()
#dev.off()

#physeq_ord_rda_micro_clr <- ordinate(dada2_clr_micro_phylo, "RDA")
#pdf(paste0(results_path, "BetaDiversity/phyloseq_micro_clr_rda_country.pdf"))
#plot_ordination(dada2, physeq_ord_rda_micro_clr, type="samples", color="Country", title="") +
#  stat_ellipse(type = "t", aes(color=Country), level = 0.95, alpha = 0.5) + 
#  theme_bw()
#dev.off()


#####################################################################################
############################### Using microbiome package ############################


# perform principle components analysis on transformed data
#dada2_clr_pca <- prcomp(dada2_clr_comp)
#dada2_clr_pca_scores <- scores(dada2_clr_pca)
#dada2_clr_pca_scores_sub <- dada2_clr_pca_scores[,1:2]
#Variance Explained by FIrst 2 axes  
#axes <- summary(dada2_clr_pca)$importance[,1:2] # (not a lot)
# Add Sample Data 
#final_clr_pca <- cbind(dada2_clr_pca_scores_sub, sampdata)

#pdf(paste0(results_path, "BetaDiversity/microcomp_clr_pca_country.pdf"))
#pca_clr_comp_plot <- ggplot(final_clr_pca,
#                            aes(x = PC1, y= PC2)) + 
#  stat_ellipse(type = "t", aes(color=Country), level = 0.95, alpha = 0.5) + 
#  geom_point(aes(colour = Country), size=0.25) + 
#  theme_bw() + 
#  labs(fill = "Country", x = paste0("PC1 (", round(axes[,1], 3), "%)"), y = paste0("PC2 (", round(axes[,2], 3), "%)")) #+ 
#theme(legend.position = "top", 
#axis.text = element_text(size=18),
#axis.title = element_text(size=20),
#legend.text = element_text(size=12),
#legend.title = element_text(size=18))
#print(pca_clr_comp_plot)
#dev.off()



dada2_ilr_pca <- prcomp(dada2_ilr)
dada2_ilr_pca_scores <- scores(dada2_ilr_pca)
dada2_ilr_pca_scores_sub <- dada2_ilr_pca_scores[,1:2]
#Variance Explained by FIrst 2 axes  
axes <- summary(dada2_ilr_pca)$importance[,1:2] # (not a lot)
# Add Sample Data  
final_ilr_pca <- cbind(dada2_ilr_pca_scores_sub, sampdata)

pdf(paste0(results_path, "BetaDiversity/microcomp_ilr_pca_country.pdf"))
pca_ilr_plot <- ggplot(final_ilr_pca,
                       aes(x = PC1, y= PC2)) + 
  stat_ellipse(type = "t", aes(color=Country), level = 0.95, alpha = 0.5) + 
  geom_point(aes(colour = Country), size=0.25) + 
  theme_bw() + 
  labs(fill = "Country", x = paste0("PC1 (", round(axes[,1], 3), "%)"), y = paste0("PC2 (", round(axes[,2], 3), "%)")) #+ 
#theme(legend.position = "top", 
#axis.text = element_text(size=18),
#axis.title = element_text(size=20),
#legend.text = element_text(size=12),
#legend.title = element_text(size=18))
print(pca_ilr_plot)
dev.off()


############################ Bd
final_ilr_pca_bd <- final_ilr_pca[!is.na(final_ilr_pca$Bd),]
pdf(paste0(results_path, "BetaDiversity/ilr_pca_bd_noNA.pdf"))
pca_ilr_plot <- ggplot(final_ilr_pca_bd,
                       aes(x = PC1, y= PC2)) + 
  stat_ellipse(type = "t", aes(color=Bd), level = 0.95, alpha = 0.5) + 
  geom_point(aes(colour = Bd), size=0.25) + 
  theme_bw() + 
  labs(fill = "Bd Status", x = paste0("PC1 (", round(axes[,1], 3), "%)"), y = paste0("PC2 (", round(axes[,2], 3), "%)"))
print(pca_ilr_plot)
dev.off()

# Then transform to compositional data
#abun <- abundances(dada2)
#abun_c <- transform(abun, "compositional")
#colnames(abun_c) <- colnames(abun)

#if (any(abun_c == 0)) {
#  v <- as.vector(abun_c)
#  minval <- min(v[v > 0])/2
#  abun_c <- abun_c + minval
#}

#no_zero <- multLN(X=abun_c, label=min(as.vector(abun_c)[as.vector(abun_c) > 0]), dl=p)
#non_zero <- apply(abun_c, 2, function (x) zeroreplace(x, d=min(as.vector(x)[as.vector(x) > 0])/2, a=0.65))
#abun_0 <- zeroreplace(abun_c, d=p, a=0.65)
#p <- matrix(min(as.vector(abun_c)[as.vector(abun_c) > 0]/2), 9292, 335)
# pca
#dada2_clr_micro_pca <- prcomp(dada2_clr_micro)

#Cool Biplot Showing How Diff ASVs affect the primary axes of the ordinatiton
#biplot(dada2_clr_pca, type = "points")
#Scree plot of relative variance explained by sequential axes
#plot(dada2_clr_pca)
#Variance Explained by FIrst 2 axes  
#axes <- summary(dada2_clr_micro_pca)$importance[,1:2] # (not a lot)

#dada2_clr_micro_pca_scores <- scores(dada2_clr_micro_pca)
#dada2_clr_micro_pca_scores_sub <- dada2_clr_micro_pca_scores[,1:2]
# Add Sample Data  
#final_clr_micro_pca <- cbind(dada2_clr_micro_pca_scores_sub, sampdata)

# plot
# plot
#pdf(paste0(results_path, "BetaDiversity/micro_clr_pca_country.pdf"))
#pca_plot <- ggplot(final_clr_micro_pca,
#                   aes(x = PC1, y= PC2)) + 
#  stat_ellipse(type = "t", aes(color=Country), level = 0.95, alpha = 0.5) + 
#  geom_point(aes(colour = Country), size=0.25) + 
#  theme_bw() + 
#  labs(fill = "Country", x = paste0("PC1 (", round(axes[,1], 3), "%)"), y = paste0("PC2 (", round(axes[,2], 3), "%)")) #+ 
#theme(legend.position = "top", 
#axis.text = element_text(size=18),
#axis.title = element_text(size=20),
#legend.text = element_text(size=12),
#legend.title = element_text(size=18))
#print(pca_plot)
#dev.off()

for (aspect in samp_vars) {
  pdf(paste0(results_path, "BetaDiversity/ilr_pca_", aspect, ".pdf"))
  pca_plot <- ggplot(final_ilr_pca,
                   aes(x = PC1, y= PC2)) + 
  stat_ellipse(type = "t", aes(color=get(aspect)), level = 0.95, alpha = 0.5) + 
  geom_point(aes(colour = get(aspect)), size=0.25) + 
  theme_bw() + 
  labs(colour = aspect, x = paste0("PC1 (", round(axes[,1], 3), "%)"), y = paste0("PC2 (", round(axes[,2], 3), "%)")) 
  print(pca_plot)
  dev.off()
}


######################################################################################
#################################### PERMANOVA #######################################


permanova_results <- list()
for (aspect in samp_vars) {
  if (aspect != "Bd") {
    permanova_calc <- adonis(dada2_ilr ~ get(aspect), data = sampdata, permutations=999,method="euclidean")
    permanova_results <- append(permanova_results, permanova_calc)
    permanova_out <- capture.output(permanova_calc)
    cat(paste0("\nPERMANOVA for ", aspect, ":\n"), permanova_out,
        file=paste0(results_path, "BetaDiversity/permanovas.txt"), sep = "\n", append=TRUE)
  }
}
# add Bd permanova (NA issue)
ilr_perm_bd <- adonis(dada2_ilr_bd ~ Bd, data = sampdata, permutations=999,method="euclidean")
permanova_results <- append(permanova_results, ilr_perm_bd)
cat(paste0("\nPERMANOVA for ", Bd, ":\n"), ilr_perm_bd,
    file=paste0(results_path, "BetaDiversity/permanovas.txt"), sep = "\n", append=TRUE)

# calculate q-values with Holm method
p_vals <- c()
Rs <- c()
for (p_stat in 1:length(permanova_results)) {
  var_stats <- as.data.frame(permanova_results[p_stat][["aov.tab"]])
  Rs <- c(Rs, var_stats$R2[1])
  p_vals <- c(p_vals, var_stats$`Pr(>F)`[1])
}
q_vals <- p.adjust(p_vals, method ="holm")
for (no in 1:length(samp_vars)) {
cat(paste0("\n", samp_vars[no], "\nQval: ", q_vals[no], "\nR^2: "), Rs[no], 
    file=paste0(results_path, "BetaDiversity/permanovas.txt"), sep = "", append=TRUE)
}

#adonis(dada2_ilr ~ Country+Amphibian_Type+Genus_Species, data = sampdata, permutations=999,method="euclidean")


## end of script