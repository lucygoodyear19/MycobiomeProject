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
library("tidyr") # for separate function
library("phyloseq")
library("microbiome") # for abundances() function
library("vegan") # for adonis() function
library("compositions") # for ilr transform
library("ggplot2")

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

# check if results directory exists and if not, create it
ifelse(!dir.exists(file.path(paste0(results_path, "/BetaDiversity/"))), dir.create(file.path(paste0(results_path, "/BetaDiversity/"))), FALSE)
# set directory for results to be sent to
path_out <- paste0(results_path, "/BetaDiversity/")


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
sampdata_bd <- sampdata_bd %>%
  separate(A_Genus_Species, c("A_Genus", "A_Species"), " ", remove=FALSE)
# transform data to compositional
abun_bd_c <- transform(abun_bd, "compositional")
# account for zeros by performing replacement with small scaled numeric
detectlims_bd <- matrix(min(as.vector(abun_bd_c)[as.vector(abun_bd_c) > 0]/2), nrow(abun_bd_c), ncol(abun_bd_c))
abun0bd <- zeroreplace(abun_bd_c, d=detectlims_bd, a=0.65) # calculates detection limits within function
# transpose matrix
abun0t_bd <- t(abun0bd)
# transform data using isometric log-ratio transform
dada2_ilr_bd <- as.data.frame(clr(abun0t_bd))
# save tranformed data for future analyses
saveRDS(dada2_ilr_bd, paste0(results_path, "ilr_abun_bd.rds"))


#####################################################################################
################################### Perform PCA #####################################


print("Performing PCA...")

# perform PCA on non-NA Bd ilr transformed data
ilr_bd_pca <- prcomp(dada2_ilr_bd)
# calculate values per sample for all axes
ilr_bd_pca_scores <- scores(ilr_bd_pca)
# extract principle axes
ilr_bd_pca_scores_sub <- ilr_bd_pca_scores[,1:2]
# extract variance explained by principle axes 
axes_bd <- summary(ilr_bd_pca)$importance[,1:2]
# add sample data 
ilr_bd_pca_final <- cbind(ilr_bd_pca_scores_sub, sampdata_bd)


#####################################################################################
################################### Plots ###########################################


# plot for all variables
for (aspect in samp_vars) {
  pdf(paste0(path_out, "ilr_pca_", aspect, ".pdf"))
  ilr_pca_plot <- ggplot(ilr_bd_pca_final,
                         aes(x = PC1, y= PC2)) + 
    stat_ellipse(type = "t", aes(color=get(aspect)), level = 0.95, alpha = 0.5) + # plot ellipses
    geom_point(aes(colour = get(aspect)), size=1) +
    labs(colour = aspect, x = paste0("PC1 (", round(axes_bd[,1], 2), "%)"), y = paste0("PC2 (", round(axes_bd[,2], 2), "%)")) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), # remove grid lines
          panel.grid.minor = element_blank(),
          axis.line.y = element_line(colour = "black"), # make axes darker
          axis.line.x = element_line(colour = "black"))
  print(ilr_pca_plot)
  dev.off()
}


######################################################################################
#################################### PERMANOVA #######################################


print("Performing PERMANOVA...")

# perform PERMANOVA on 3 variables
permanova <- adonis2(dada2_ilr_bd ~ get(samp_vars_1host[1])+get(samp_vars_1host[2])+get(samp_vars_1host[3]), by = "margin", data = sampdata_bd, permutations=999,method="euclidean")
# save to text file
cat("\n\nPERMANOVA: ", capture.output(permanova_out),
    file=paste0(path_out, "permanovas_beta.txt"), sep = "\n", append=TRUE)

# calculate dispersion significance for the 3 variables
# subset sample data by the 3 variables
sampsub <- sampdata_bd[,samp_vars_1host]
disp_ob <- cbind.data.frame(sampsub, dada2_ilr_bd)
# create distance matrix
dist <- dist(disp_ob[,4:ncol(disp_ob)], "euclidean")
# calculate dispersion for the 3 variables
for (var in 1:length(samp_vars_1host)){
  dispersion <- betadisper(dist, group=disp_ob[[samp_vars_1host[var]]],type="median",bias.adjust=TRUE)
  # extract ANOVA
  anova_disp <- anova(dispersion)
  # save dispersion ANOVAs to text file
  cat(paste0("\n\nANOVA for dispersion of ", samp_vars_1host[var], ": "), capture.output(anova_disp),
      file=paste0(path_out, "permanovas_beta.txt"), sep = "\n", append=TRUE)
}
  

## end of script
