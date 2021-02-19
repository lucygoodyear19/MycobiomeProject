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


library("dplyr") # for data wrangling
library("phyloseq") # to work with phyloseq objects
library("microbiomeSeq") # for Kruskal-Wallis function
library("effsize") # for Cohen D effect size function
library("ggplot2") # for plotting
library("microbiome") # for abundances()
library("compositions") #for zeroreplace() and clr()


################################ User inputs ################################


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) # setup to accept arguments from command line
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied in the form of an R-script containing the following arguments: 
       1) dada2_data_path - path to directory containing input data (phyloseq object)
       2) results_path - path to results directory")
}
# source analysis arguments
source(args)
# print arguments as check
print(paste0("Data path: ", dada2_data_path))
print(paste0("Results path: ", results_path))

# load phyloseq object
print("Loading data...")
dada2 <- readRDS(paste0(dada2_data_path, "physeqob_DADA2.rds"))

# set Bd as factor
sample_data(dada2)$Bd <- as.factor(as.character(sample_data(dada2)$Bd))

# check if results directory exists and if not, create it
ifelse(!dir.exists(file.path(paste0(results_path, "/Abundances_vs_Bd/"))), dir.create(file.path(paste0(results_path, "/Abundances_vs_Bd/"))), FALSE)
# set directory for results to be sent to
path_out <- paste0(results_path, "Abundances_vs_Bd/")

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
  #Subset phyloseq object to remove all samples with NA in Bd column
  bd <- prune_samples(!is.na(sample_data(get(j))$Bd), get(j))
  assign(j, bd)
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


print("Completed microbiomeSeq Kruskal-Wallis")


###############################################################################
######################## Kruskal-Wallis using base R ##########################


################################ Functions ####################################


# function for plotting pie charts
# inputs:
# 1) taxa_df - data frame containing taxonomy data
# 2) tax_rank - taxonomic rank to split pie chart by 
# Outputs:
# 1) pie chart saved to path_out directory
plot_bd_pies <- function(taxa_df, tax_rank) {
  # subset by taxa ta
  rank_bd <- as.data.frame(table(taxa_df[[tax_rank]]))
  names(rank_bd) <- c(tax_rank, "Frequency")
  rank_bd <- rank_bd[!rank_bd$Frequency == 0,] # remove any taxa with no samples
  # plot pie chart and save to pdf
  pdf(paste0(path_out, "base_pie_", tax_rank, ".pdf"))
  print(ggplot(rank_bd, aes(x= "", y= Frequency, fill= get(tax_rank))) +
          geom_bar(stat = "identity", width = 1) +
          coord_polar("y", start = 0) +
          guides(fill=guide_legend(title=tax_rank)) +
          theme_void())
  dev.off()
}


############################## Transform data #################################

for (j in (1:length(datasets))) {
  #Subset phyloseq object to remoave all samples with NA in Bd column
  bd <- prune_samples(!is.na(sample_data(get(datasets[j]))$Bd), get(datasets[j]))
  # remove any taxa with sum 0 abundance across all samples
  bd_nonzero <- prune_taxa(taxa_sums(bd) > 0, bd)
  # extract otu_table
  abun_bd <- abundances(bd_nonzero)
  # extract sample data
  sampdata_bd <-as(sample_data(bd_nonzero),"data.frame")
  # transform data to compositional
  abun_bd_c <- transform(abun_bd, "compositional")
  # account for zeros by performing replacement with small scaled numeric
  detectlims_bd <- matrix(min(as.vector(abun_bd_c)[as.vector(abun_bd_c) > 0]/2), nrow(abun_bd_c), ncol(abun_bd_c))
  abun0bd <- zeroreplace(abun_bd_c, d=detectlims_bd, a=0.65)
  # transpose matrix
  abun0t_bd <- t(abun0bd)
  # transform data using isometric log-ratio transform
  clr_bd <- as.data.frame(clr(abun0t_bd))
  # replace otu table with clr transformed data
  otu_table(bd) <- otu_table(clr_bd, taxa_are_rows = FALSE)
  assign(datasets[j], bd)
}


#################################### Kruskal-Wallis ###################################


# perform Kruskal-Wallis for all ASVs grouped by Bd +ve/-ve
print("Running base R Kruskal-Wallis")
for (j in datasets) {
  print(paste0("Starting analysis for ", j))
  # extract otu_table from transformed phyloseq object
  abun <- t(abundances(get(j)))
  # add Bd status to transformed abundances
  bd <- as.character(sample_data(get(j))$Bd)
  samp_bd <- cbind(bd, abun)
  # convert all data to numeric type and turn back to dataframe
  mode(samp_bd) <- "numeric"
  samp_bd <- as.data.frame(samp_bd)
  # set Bd status as factor
  samp_bd[,1] <- samp_bd %>% select(bd) %>% mutate_at(vars(bd), funs(factor))
  # save p-values to empty vector
  pvals <- c()
  # perform Kruskal_Wallis test and save p-values
  for (i in 2:ncol(samp_bd)) {
    kw <- kruskal.test(samp_bd[,i] ~ bd, data = samp_bd)
    pvals <- c(pvals, kw$p.value)
  }
  # name p-values by ASV/OTU
  names(pvals) <- colnames(samp_bd)[2:ncol(samp_bd)]
  # create list of adjusted p-values
  # Multiple testing correction of p-values using Benjamini-Holm method
  qvals_bh = p.adjust(pvals, method ="BH")
  # sort by most significant
  qvals_bh_ord <- sort(qvals_bh)
  # only keep p-values below 0.05
  qvals_ls <- Filter(function(x) x < 0.05, qvals_bh_ord)
  # calculate Cohen's d for effect size and save to text file
  taxas <- c()
  qvals <- c()
  cohens <- c()
  # view taxonomic assignment for significant ASVs
  sig_taxa <- as.data.frame(tax_table(get(j))[names(qvals_ls),])
  # save all pertinent information in text file
  for (g in 1:length(qvals_ls)) {
    cohen_out <- capture.output(cohen.d(samp_bd[[names(qvals_ls)[g]]] ~ samp_bd$bd))
    taxa_out <- sig_taxa[g,]
    pval_out <- capture.output(pvals[[names(qvals_ls)[g]]])
    qval_out <- capture.output(qvals_ls[g])
      
    taxas <- rbind(taxas, taxa_out)
    qvals <- c(qvals, qval_out[2])
    cohens <- c(cohens, cohen_out[4])

    cat(paste0("Significant ASVs using KW and BH correction method: \n"), capture.output(taxa_out), file=paste0(path_out, "base_Abundance_associated_with_Bd.txt"), sep = "\n", append=TRUE)
    cat(paste0("\np-value\n"), pval_out, file=paste0(path_out, "base_Abundance_associated_with_Bd.txt"), sep = "\n", append=TRUE)
    cat(paste0("\nAdjusted p-value\n"), qval_out, file=paste0(path_out, "base_Abundance_associated_with_Bd.txt"), sep = "\n", append=TRUE)
    cat("\nCohen D\n", cohen_out, file=paste0(path_out, "base_Abundance_associated_with_Bd.txt"), sep = "\n", append=TRUE)
  
    # plot abundance compared Bd status for each sample for each significant ASV
    pdf(paste0(path_out, "base_abun_bd_plot.pdf"))
    for (i in 1:length(qvals_ls)) {
      plot_kw <- cbind(samp_bd$bd, samp_bd[names(qvals_ls)[i]])
      print(ggplot(plot_kw,aes(x=samp_bd$bd, y=samp_bd[,names(qvals_ls)[i]])) + 
              geom_point(size = 3) + 
              labs(x = "Bd +ve/-ve", y = names(qvals_ls)[i]) +
              annotate(geom="text", x=1, y=max(samp_bd[,names(qvals_ls)[i]]),label=paste0("Taxonomy:\n", taxas[i,2], "\n", taxas[i,3], "\n", taxas[i,4], "\n", taxas[i,5], "\n", taxas[i,6], "\n", taxas[i,7], "\nAdjusted p-value: ", qvals[i], "\nCohen ", cohens[i])) + 
              theme_bw())
    }
  dev.off()
  }
  # function to plot pie chart showing taxa rank splits of sig taxa
  # remove prefixes in taxa
  taxas <- as.data.frame(lapply(taxas, function(x) gsub(".__", "", x)))
  # plot pie charts
  plot_bd_pies(taxas, "Phylum")
  plot_bd_pies(taxas, "Class")
  plot_bd_pies(taxas, "Order")
  print(paste0("Completed analysis for ", j))
}

print("All analysis completed")


## end of script
