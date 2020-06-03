##############################################################################
############# Abundance of fungal species compared to Bd +ve/-ve #############
##############################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


##############################################################################
################################## Set up ####################################


############################## Load packages ################################


library(dplyr) # for data wrangling
library(phyloseq) # to work with phyloseq objects
library(Biostrings) # TEMP to rename ASVs
library(microbiomeSeq) # for Kruskal-Wallis function
library(effsize) # for Cohen D effect size function
library(ggplot2) # for plotting


################################ User inputs ################################


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) # setup to accept arguments from command line
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied in the form of an R-script containing the following arguments: 
       1) esto_data_path: full path to esto data
       2) dada2_data_path: full path to dada2 data
       3) results_path: full path to general results directory
       4) datasets: vector of names of datasets to be analysed")
}

# source analysis arguments
source(args)

# print args as chack
print(esto_data_path)
print(dada2_data_path)
print(results_path)
print(datasets)

# set wd
#setwd(paste0(root, "Runs_Countries/", run))

# load filtered phyloseq object for both Estonian pipeline and DADA2 pipeline
dada2 <- readRDS(dada2_data_path)
esto <- readRDS(esto_data_path)

# check if results directory exists and if not, create it
ifelse(!dir.exists(file.path(paste0(results_path, "/Abundances_vs_Bd/"))), dir.create(file.path(paste0(results_path, "/Abundances_vs_Bd/"))), FALSE)
# set directory for results to be sent to
path_out <- paste0(results_path, "/Abundances_vs_Bd/")

##### Temporary set up
# rename ASVs for convenience (now included in data filtering script 
# so will be deleted for future analyses)
asv_seqs <- Biostrings::DNAStringSet(taxa_names(dada2))
names(asv_seqs) <- taxa_names(dada2)
dada2 <- merge_phyloseq(dada2, asv_seqs)
taxa_names(dada2) <- paste0("ASV", seq(ntaxa(dada2)))
dada2

print("Set up complete")


###############################################################################
##################### Kruskal-Wallis using microbiomeSeq ######################


# run Kruskal_Wallis test on Bd on abundances for all datasets

# note taxa_as_rows in phyloseq object must be false for function 
# to kruskal_abundance to work

# plot box plots using base R and scatter plot using phyloseq function
print("Running microbiomeSeq Kruskal-Wallis")
for (j in datasets) {
  print(paste0("Starting analysis for ", j))
  # perform Kruskal_Wallis
  (assign(paste0("kw_", j), kruskal_abundance(get(j), group = "Bd", pvalue.threshold = 0.05)))
  # save taxa details of significant ASVs
  sig_asv <- get(paste0("kw_", j))$importance[,1]
  sig_tax <- capture.output(tax_table(get(j))[rownames(tax_table(get(j))) %in% sig_asv,])
  to_save <- capture.output(get(paste0("kw_", j)))
  cat(paste0("Kruskal Wallis for ", j, "\n"), to_save, file=paste0(path_out, "mbseq_Abundance_associated_with_Bd_for_", j, ".txt"), sep = "\n", append=TRUE)
  cat("\nTaxa Details for significant ASVs\n", sig_tax, file=paste0(path_out, "mbseq_Abundance_associated_with_Bd_for_", j, ".txt"), sep = "\n", append=TRUE)
  # plot
  pdf(paste0(path_out, "mbseq_", j, "_abun_bd_plot.pdf"))
  print(plot_signif(get(paste0("kw_", j))$plotdata))
  data <- get(paste0("kw_", j))$plotdata
  data_ord <- data[order(data$Rank),]
  for (i in unique(data_ord$Rank)) {
    to_plot <- subset(data_ord, Rank==i)
    print(ggplot(subset(to_plot),aes(Groups,Value,colour=Groups))+
        geom_point(size = 3)+
        ggtitle(to_plot$Taxa[1]) +
        theme_bw())
  }
  dev.off()
  print(paste0("Completed analysis for ", j))
}

# detach package
detach("package:microbiomeSeq", unload=TRUE)

print("Completed microbiomeSeq Kruskal-Wallis")


###############################################################################
######################## Kruskal-Wallis using base R ##########################


# see: https://wiki.bits.vib.be/index.php/Ecology_analysis_using_vegan#Alternative_normalization_methods

# perform Kruskal-Wallis for all ASVs grouped by Bd +ve/-ve
print("Running base R Kruskal-Wallis")
for (j in datasets) {
  print(paste0("Starting analysis for ", j))
  ## get abundance data into right format
  # transform to data frame and then matrix
  samps <- as.data.frame(otu_table(get(j)))
  samps <- as.matrix(samps)
  # create matrix of proportions of ASVs per sample
  samp_prop <- samps/rowSums(samps)
  # add Bd status to samp_prop
  bd <- as.character(sample_data(get(j))$Bd)
  samp_prop_bd <- cbind(bd, samp_prop)
  # convert all data to numeric type and turn back to dataframe
  mode(samp_prop_bd) <- "numeric"
  samp_prop_bd <- as.data.frame(samp_prop_bd)
  # set Bd status as factor
  samp_prop_bd[,1] <- samp_prop_bd %>% select(bd) %>% mutate_at(vars(bd), funs(factor))
  # save p-values to empty vector
  pvals <- c()
  # perform Kruskal_Wallis test and save p-values
  for (i in 2:ncol(samp_prop_bd)) {
    kw <- kruskal.test(samp_prop_bd[,i] ~ bd, data = samp_prop_bd)
    pvals <- c(pvals, kw$p.value)
  }
  # name p-values by ASV/OTU
  names(pvals) <- colnames(samp_prop_bd)[2:ncol(samp_prop_bd)]
  # create list of adjusted p-values via two methods
  # Multiple testing correction of p-values using Benjamini-Hochberg method
  qvals_hoch = p.adjust(pvals, method ="BH")
  # sort by most significant
  qvals_hoch_ord <- sort(qvals_hoch)
  # only keep p-values below 0.05
  qvals_hoch_sig <- Filter(function(x) x < 0.05, qvals_hoch_ord)
  # Multiple testing correction of p-values using Holm's method
  qvals_holm = p.adjust(pvals, method ="holm")
  # sort by most significant
  qvals_holm_ord <- sort(qvals_holm)
  # only keep p-values below 0.05
  qvals_holm_sig <- Filter(function(x) x < 0.05, qvals_holm_ord)
  
  # calculate effect size for all significant ASVs using Cohen D
  qval_ls <- c("qvals_hoch_sig", "qvals_holm_sig")
  methods <- c("Benjamini-Hochberg", "Holm")
  # calculate Cohen's d for effect size and save to text file
  for (k in 1:length(qval_ls)) {
    data <- get(qval_ls[k])
    method <- methods[k]
    taxas <- c()
    qvals <- c()
    cohens <- c()
    # view taxonomic assignemnt for significant ASVs
    sig_taxa <- as.data.frame(tax_table(get(j))[names(data),])
    for (g in 1:length(data)) {
      cohen_out <- capture.output(cohen.d(samp_prop_bd[[names(data)[g]]] ~ samp_prop_bd$bd))
      taxa_out <- sig_taxa[g,]
      pval_out <- capture.output(pvals[[names(data)[g]]])
      qval_out <- capture.output(data[g])
      
      taxas <- rbind(taxas, taxa_out)
      qvals <- c(qvals, qval_out[2])
      cohens <- c(cohens, cohen_out[4])

      cat(paste0("Significant ASVs using ", method, " method\n"), capture.output(taxa_out), file=paste0(path_out, "base_Abundance_associated_with_Bd_for_", j, ".txt"), sep = "\n", append=TRUE)
      cat(paste0("\np-value\n"), pval_out, file=paste0(path_out, "base_Abundance_associated_with_Bd_for_", j, ".txt"), sep = "\n", append=TRUE)
      cat(paste0("\nAdjusted p-value\n"), qval_out, file=paste0(path_out, "base_Abundance_associated_with_Bd_for_", j, ".txt"), sep = "\n", append=TRUE)
      cat("\nCohen D\n", cohen_out, file=paste0(path_out, "base_Abundance_associated_with_Bd_for_", j, ".txt"), sep = "\n", append=TRUE)
    }
  
    # plots
    print("Plotting...")
    pdf(paste0(path_out, "base_", j, "_", method, "_abun_bd_plot.pdf"))
    for (i in 1:length(data)) {
      plot_kw <- cbind(samp_prop_bd$bd, samp_prop_bd[names(data)[i]])
      print(ggplot(plot_kw,aes(x=samp_prop_bd$bd, y=samp_prop_bd[,names(data)[i]])) + 
              geom_point(size = 3) + 
              labs(x = "Bd +ve/-ve", y = names(data)[i]) +
              annotate(geom="text", x=1, y=max(samp_prop_bd[,names(data)[i]]),label=paste0("Taxonomy:\n", taxas[i,2], "\n", taxas[i,3], "\n", taxas[i,4], "\n", taxas[i,5], "\n", taxas[i,6], "\n", taxas[i,7], "\nAdjusted p-value: ", qvals[i], "\nCohen ", cohens[i])) + 
              theme_bw())
    }
    dev.off()
    print("Completed plots")
  }
  print(paste0("Completed analysis for ", j))
}

print("All analysis completed")


## end of script
