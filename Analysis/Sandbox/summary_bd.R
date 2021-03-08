################################################################################################
####################################### Bd Summary Data ########################################
################################################################################################


# Author: Luke Goodyear (leg19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


##############################################################################################
######################################### Set up #############################################


# load packages
library("phyloseq")
library("ggplot2")
theme_set(theme_bw())
library("dplyr")
library("viridis")
library("microbiome") # for abundances()
library("tidyr")
library("VennDiagram")


# import arguments to run script on specific Latitude data
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
ifelse(!dir.exists(file.path(paste0(results_path, "/bd_summary/"))), dir.create(file.path(paste0(results_path, "/bd_summary/"))), FALSE)
# set directory for results to be sent to
path_out <- paste0(results_path, "/bd_summary/")

# extract sample data from phyloseq object
samps <- as(sample_data(dada2), "data.frame")
# extract taxa table as dataframe from phyloseq object
taxa <- as.data.frame(tax_table(dada2))
#taxa <- as.data.frame(lapply(taxa, function(x) gsub(".__", "", x)))
# extract ASV abundance table as dataframe from phyloseq object
seqs <- as.data.frame(t(abundances(dada2)))
# transpose abundance data frame
seqs_t <- as.data.frame(t(seqs))

# set to presence absence
seqs_pres_abs <- seqs
seqs_pres_abs[seqs_pres_abs > 0] <- 1
seqs_pres_abs_t <- as.data.frame(t(seqs_pres_abs))

# set up colour blind friendly palatte
cbpalette <- c("#4EB265", "#6195CF", "#CAE0AB", "#994F88", "#1965B0", 
               "#F7F056", "#90C987", "#F1932D", "#D9CCE3", "#437DBF", 
               "#882E72", "#72190E", "#CAACCB", "#7BAFDE", "#DC050C", 
               "#BA8DB4", "#A5170E", "#E8601C", "#F6C141", "#AA6F9E")


##################################################################################################
#################################### Pie charts of taxa ##########################################

taxa$Phylum <- as.character(taxa$Phylum)
for (i in 1:nrow(taxa)) {
  if (is.na(taxa$Phylum[i])) {taxa$Phylum[i] = "NA"}
}
taxa$Phylum <- as.factor(taxa$Phylum)
phyla <- as.data.frame(table(taxa$Phylum))
names(phyla) <- c("Phyla", "Frequency")

# plot pie chart
pdf(paste0(path_out, "phyla_pie.pdf"))
print(ggplot(phyla, aes(x= "", y= Frequency, fill= Phyla)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        theme_void())
dev.off()


#################################################################################################
############################# Taxa composition by presence/absence ##############################


# function for plot by presence absence by chosen factor
# inputs:
# 1) taxa_data - data frame containing taxonomy data
# 2) seq_data - data frame containing abundance data
# 3) samp_data - data frame containing sample data
# 4) tax_rank - taxonomic rank to plot, as string (matching column in taxa_data)
# 5) condition - name of variable chose, as string
# Outputs:
# 1) bar chart showing break down by tax_rank per condition, saved as pdf
find_taxa_breakdown <- function(taxa_data, seq_data, samp_data, taxa_rank, condition) {
  taxon_table <- as.data.frame(c())
  labs <- c() # to store plot labels
  factors <- c() # to ensure bars are in the right order
  # for loop over factor levels in condition
  for (var in unique(samp_data[[condition]])) {
    if (!is.na(var)) {
      samps_sub <- samp_data %>% filter(samp_data[[condition]] == var) # subset by factor level in condition
      labs <- c(labs, paste0(samps_sub[[condition]][1], " (n = ", nrow(samps_sub), ")")) 
      factors <- c(factors, var) 
      samp_seqs <- seq_data[colnames(seq_data) %in% rownames(samps_sub)] # subset seqs
      samp_seqs$sum <- rowSums(samp_seqs) # sum ASV abundance in factor
      samp_seqs <- samp_seqs[samp_seqs$sum > 0,] # remove any ASV not equal to 0
      # subset taxa by taxa_rank
      samp_taxa <- taxa_data[rownames(taxa_data) %in% rownames(samp_seqs),]
      taxon <- as.character(samp_taxa[[taxa_rank]])
      taxon[is.na(taxon)] <- "NA" # set NA to string to ensure inclusion in count
      count <- as.data.frame(table(taxon))
      count$taxon[count$taxon == "NA"] <- NA # set "NA" back to NA
      count$Freq <- count$Freq/sum(count$Freq) # calculate percentages
      count$taxon <- gsub(".__", "", count$taxon) # remove prefix in taxa
      count$Variable <- var
      taxon_table <- rbind(taxon_table, count)
    }
  }
  # ensure bars are in correct order
  taxon_table$Variable <- factor(taxon_table$Variable, ordered=TRUE, levels = factors)
  # plot stacked bar chart for all conditions
  pdf(paste0(results_path, "bd_summary/", condition, "_", taxa_rank, "_breakdown.pdf"))
  print(ggplot(taxon_table, aes(fill=taxon, y=Freq, x=Variable)) + 
          geom_bar(position="stack", stat="identity") +
          labs(fill = taxa_rank) +
          ylab("Proportion") +
          xlab(condition) +
          scale_x_discrete(labels = labs) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
}

find_taxa_breakdown(taxa, seqs_pres_abs_t, samps, "Phylum", "Bd")
find_taxa_breakdown(taxa, seqs_pres_abs_t, samps, "Class", "Bd")


#################################################################################################
################################ Taxa composition by abundance ##################################


# function for plot by abundance by chosen factor
# inputs:
# 1) taxa_data - data frame containing taxonomy data
# 2) seq_data - data frame containing abundance data
# 3) samp_data - data frame containing sample data
# 4) tax_rank - taxonomic rank to plot, as string (matching column in taxa_data)
# 5) condition - name of variable chose, as string
# Outputs:
# 1) bar chart showing break down by tax_rank per condition, saved as pdf
find_taxa_breakdown_abun  <- function(taxa_data, seq_data, samp_data, taxa_rank, condition) {
  taxon_table <- as.data.frame(c())
  labs <- c() # to store plot labels
  factors <- c() # to ensure bars are in the right order
  # for loop over factor levels in condition
  for (var in unique(samp_data[[condition]])) {
    if (!is.na(var)) {
      samps_sub <- samp_data %>% filter(samp_data[[condition]] == var) # subset by factor level in condition
      labs <- c(labs, paste0(samps_sub[[condition]][1], " (n = ", nrow(samps_sub), ")"))
      factors <- c(factors, var)
      samp_seqs <- seq_data[colnames(seq_data) %in% rownames(samps_sub)] # subset seqs
      samp_seqs$sum <- rowSums(samp_seqs) # sum ASV abundance in factor
      samp_seqs <- samp_seqs[samp_seqs$sum > 0,] # remove any ASV not equal to 0
      total <- sum(samp_seqs$sum)
      # subset taxa by taxa_rank
      samp_taxa <- taxa_data[rownames(taxa_data) %in% rownames(samp_seqs),]
      taxon <- samp_taxa[[taxa_rank]]
      pcts <- samp_seqs$sum/total # calculate percentages
      count <- cbind.data.frame(taxon, pcts)
      # merge similar taxa and sum percenatges
      require("dplyr")
      count <- count %>%
        group_by(taxon) %>%
        summarise(Percentage = sum(pcts))
      count$Variable <- var
      count$taxon <- gsub(".__", "", count$taxon) # remove prefix in taxa
      taxon_table <- rbind.data.frame(taxon_table, count)
    }
  }
  # set NA to string to ensure inclusion in count
  taxon_table$taxon <- as.character(taxon_table$taxon)
  taxon_table$taxon[is.na(taxon_table$taxon)] <- "NA"
  taxon_table$taxon <- as.factor(taxon_table$taxon)
  # set NA to the end of list of factor levels
  taxon_table$taxon <- fct_relevel(taxon_table$taxon, "NA", after = Inf)
  # ensure bars are in correct order
  taxon_table$Variable <- factor(taxon_table$Variable, ordered=TRUE, levels = factors)
  # write breakdown by taxa by condition to csv
  write.csv(taxon_table, paste0(path_out, condition, "_", taxa_rank, "_breakdown.csv"))
  # plot stacked bar chart for all conditions
  pdf(paste0(results_path, "bd_summary/", condition, "_", taxa_rank, "_abun_breakdown.pdf"))
  print(ggplot(taxon_table, aes(fill=taxon, y=Percentage, x=Variable)) + 
          geom_bar(position="stack", stat="identity") +
          coord_fixed(ratio=18) + # change ratio between x and y axes
          scale_fill_manual(values = cbpalette) +
          # remove white space above and below bars
          scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0,0))) +
          labs(fill = taxa_rank) +
          scale_x_discrete(labels = labs) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                #axis.text = element_text(size = 16),   
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line.y = element_line(colour = "black")))
  dev.off()
}

find_taxa_breakdown_abun(taxa, seqs_t, samps, "Phylum", "Bd")

# calculate total abundance breakdown by phyla
samp_seqs <- seqs_t
samp_seqs$sum <- rowSums(samp_seqs) # sum each ASV
total <- sum(samp_seqs$sum) # total read count
# subset by phylum
taxon <- taxa[["Phylum"]]
# find percentages per ASV
pcts <- samp_seqs$sum/total
count <- cbind.data.frame(taxon, pcts)
# generate percenatges per phyla
count <- count %>%
  group_by(taxon) %>%
  summarise(Percentage = sum(pcts))
# remove taxa prefix
count$taxon <- gsub(".__", "", count$taxon)
# save total phylum breakdwn to ASV
write.csv(taxon, paste0(path_out, "total_phylum_breakdown.csv"))


#################################################################################################
######################################### Bd status #############################################


# plot stacked bar chart for Bd status per Latitude, including NA
bd_table <- as.data.frame(c())
labs <- c() # for plot labels
factors <- c() # for correct label order
for (Latitude in unique(samps$Latitude)) {
  samps_sub <- samps %>% filter(samps$Latitude == Latitude)
  samps_sub$Latitude <- as.character(samps_sub$Latitude)
  labs <- c(labs, paste0(Latitude, " (n = ", nrow(samps_sub), ")"))
  factors <- c(factors, Latitude)
  subs <- subset(samps_sub, select = c("Latitude", "Bd"))
  count <- as.data.frame(table(subs, exclude = NULL))
  count$Freq <- count$Freq/sum(count$Freq)
  bd_table <- rbind(bd_table, count)
}
bd_table$Bd <- as.character(bd_table$Bd)
bd_table$Bd[is.na(bd_table$Bd)] <- "NA"
bd_table$Bd <- as.factor(bd_table$Bd)
bd_table$Latitude <- factor(bd_table$Latitude, ordered=TRUE, levels = factors)
# plot stacked bar chart for all countries
pdf(paste0(results_path, "bd_summary/bd_breakdown.pdf"))
ggplot(bd_table, aes(fill=Bd, y=Freq, x=Latitude)) + 
  geom_bar(position="stack", stat="identity") +
  coord_fixed(ratio=5) +
  scale_fill_manual(values = c("#90C987", "#6195CF", "#F6C141")) + # colour blind friendly colours
  scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0,0))) +
  labs(fill = "Bd Status") +
  ylab("Proportion") +
  scale_x_discrete(labels = labs) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(colour = "black"))
dev.off()


############################################################################################
################################ Pie charts for Bd status ##################################


# remove any samples where Bd is NA
noNAbd_phylo <- prune_samples(!is.na(sample_data(dada2)$Bd), dada2)

# separate out positive Bd samples and associated taxa
bd_phylo <- prune_samples(sample_data(noNAbd_phylo)$Bd == 1, noNAbd_phylo)
bd_phylo <- filter_taxa(bd_phylo, function(x) mean(x) > 0, TRUE)
bd_tax <- as.data.frame(tax_table(bd_phylo))
bd_tax <- as.data.frame(lapply(bd_tax, function(x) gsub(".__", "", x)))

# separate out negative Bd samples and associated taxa
nobd_phylo <- prune_samples(sample_data(noNAbd_phylo)$Bd == 0, noNAbd_phylo)
nobd_phylo <- filter_taxa(nobd_phylo, function(x) mean(x) > 0, TRUE)
nobd_tax <- as.data.frame(tax_table(nobd_phylo))
nobd_tax <- as.data.frame(lapply(nobd_tax, function(x) gsub(".__", "", x)))

# plot pies charts function
plot_bd_pies <- function(taxa_df, condition, tax_rank) {
  rank_bd <- as.data.frame(table(taxa_df[[tax_rank]]))
  names(rank_bd) <- c(tax_rank, "Frequency")
  rank_bd <- rank_bd[!rank_bd$Frequency == 0,]
  
  pdf(paste0(results_path, "bd_summary/", tax_rank, "_pie_bd", condition, ".pdf"))
  print(ggplot(rank_bd, aes(x= "", y= Frequency, fill= get(tax_rank))) +
          geom_bar(stat = "identity", width = 1) +
          coord_polar("y", start = 0) +
          theme(legend.title = tax_rank) +
          guides(fill=guide_legend(title=tax_rank)) +
          theme_void()) 
  dev.off()
}

# plot pie charts
plot_bd_pies(bd_tax, 1, "Phylum")
plot_bd_pies(nobd_tax, 0, "Phylum")


########################################################################################
############################## Summary pie charts ######################################


samps_bd_all <- samps
samps_bd_all$Bd <- as.character(samps_bd_all$Bd)
samps_bd_all$Bd[is.na(samps_bd_all$Bd)] <- "NA"
pie_bd <- as.data.frame(table(samps_bd_all$Bd))
names(pie_bd) <- c("BdStatus", "Frequency")

pdf(paste0(results_path, "bd_summary/pie_bd.pdf"))
print(ggplot(pie_bd, aes(x= "", y= Frequency, fill= BdStatus)) +
        geom_bar(stat = "identity", width = 0.7) +
        scale_fill_viridis_d() +
        coord_polar("y", start = 0) +
        labs(fill = "Bd Status") +
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              panel.background=element_blank(),
              panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank(),
              legend.title = element_text(size = 20),
              legend.text = element_text(size = 20),
              legend.box.spacing = unit(0.001, "cm"),
              legend.spacing = unit(0.8, "cm"))) 
dev.off()

