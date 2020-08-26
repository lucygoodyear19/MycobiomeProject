##########################################################################################
################################### Co-occurrence Analysis #################################
##########################################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


#########################################################################################
######################################## Set up #########################################


# load packages
library("phyloseq")
library("cooccur")
library("microbiome") # for abundances() function
library("dplyr")
library("tidyr")
library("ggplot2")
theme_set(theme_bw())

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

print("Loading data...")

# check if results directory exists and if not, create it
ifelse(!dir.exists(file.path(paste0(results_path, "Cooccur/"))), dir.create(file.path(paste0(results_path, "Cooccur/"))), FALSE)
# set path for outputs
path_out <- paste0(results_path, "Cooccur/")

# load data
cooccur_full <- readRDS(paste0(path_out, "cooccur_genus_dada2.rds"))
dada2 <- readRDS(paste0(dada2_data_path, "physeqob_DADA2.rds"))

# extract abundances and taxa
abun <- abundances(dada2)
abun_t <- as.data.frame(t(abun))
taxa <- as.data.frame(tax_table(dada2))
# rename all taxa to remove prefix
taxa <- as.data.frame(lapply(taxa, function(x) gsub(".__", "", x)))

# look at cooccurence patterns just for Bd
cooccur_bd <- pair(cooccur_full, spp="bd")

# add columns for +ve/-ve coocurence
rownames(cooccur_bd) <- cooccur_bd$sp2
for (i in 1:nrow(cooccur_bd)) {
  if (cooccur_bd$p_gt[i] < 0.05) {
    cooccur_bd$cooccurance[i] <- "Positive"
  }
  if (cooccur_bd$p_lt[i] < 0.05) {
    cooccur_bd$cooccurance[i] <- "Negative"
  }
}

# set first column name to ASV
names(cooccur_bd)[1] <- "Genus"
# save to csv
write.csv(cooccur_bd, paste0(path_out, "cooccur_bd_table.csv"))

# remove unnecessary columns
cooccur_bd_fil <- cooccur_bd[,-c(2:5)]


###########################################################################################
#################################### Plotting Prep ########################################


#################################### Presence/Absence ######################################


# add genus to abundance table
genus <- taxa$Genus
genus_abun <- cbind.data.frame(as.character(genus), abun)

# merge by samples by genus
genus_grid <- plyr::ddply(genus_abun, "genus", plyr::numcolwise(sum))
# remove NA row
genus_grid <- subset(genus_grid, !is.na(genus))
# rename rows to genus 
rownames(genus_grid) <- genus_grid$genus
#remove genus column
genus_grid <- genus_grid[,-1]
# transpose to add Bd column
genus_grid_t <- as.data.frame(t(genus_grid))

# extract Bd column from sample data
bd <- samps$Bd
bd <- as.numeric(as.character(bd))

# add Bd binary column to otu_table
cooccur_object <- cbind(genus_grid_t, bd)

# remove samples with NA Bd qPCR
cooccur_object <- cooccur_object[!is.na(cooccur_object$bd),]
bd_samps <- cooccur_object$bd

# set all abundances > 0 to 1 to get a presence absence table
cooccur_object[cooccur_object > 0] <- 1

# remove all taxa found only in one sample
cooccur_object_t <- as.data.frame(t(cooccur_object))
cooccur_object_t$sum <- rowSums(cooccur_object_t)
cooccur_object_t <- subset(cooccur_object_t, sum >= 2)
cooccur_object_t <- subset(cooccur_object_t, select=-c(sum))

# subset to only include significant genera with regards to Bd
cooccur_sig <- subset(cooccur_object_t, rownames(cooccur_object_t) %in% cooccur_bd_fil$Genus)
cooccur_sig_t <- as.data.frame(t(cooccur_sig))


########################################### Taxa #########################################


# remove kingdom and species columns
taxa_nospec_noking <- taxa[,-c(1,ncol(taxa))]
# subset to only include significant genuses with regards to Bd
taxa_sig <- subset(taxa_nospec_noking, taxa_nospec_noking$Genus %in% cooccur_bd_fil$Genus)

# merge row to only include uniqe taxa
taxa_sig_merged <- unique(taxa_sig)
# set genus as character instead of factor
taxa_sig_merged$Genus <- as.character(taxa_sig_merged$Genus)
# put genus in alphabetical order to match cooccurrence object
taxa_sig_merged <- taxa_sig_merged[order(taxa_sig_merged$Genus),]
# add pos/neg cooccurrence status to taxa
taxa_sig_merged$Cooccurrence <- cooccur_bd_fil$cooccurance


###########################################################################################
############################### Cooccurrence by country plot ##############################


# create vector of countries for Bd tested samples
countries <- sample_data(dada2)$Country[!is.na(sample_data(dada2)$Bd)]
# add countries to cooccur object
cooccur_sig_t$Country <- countries
# count number of samples per country
no_samples <- as.data.frame(table(cooccur_sig_t$Country))
names(no_samples) <- c("Country, Frequency")

# sum rows of the same country
country_plot_merged <- cooccur_sig_t %>%
                        group_by_at("Country") %>%
                        summarise_at(c(1:ncol(cooccur_sig_t)-1), sum)
country_plot_merged <- as.data.frame(country_plot_merged)
# set countries as rownames and remove country column
rownames(country_plot_merged) <- country_plot_merged$Country
country_plot_merged <- country_plot_merged[,-1]
# transpose
country_plot_merged_t <- as.data.frame(t(country_plot_merged))

# add co-occurrence pos/neg status to country data frame
Cooccurence <- cooccur_bd$cooccurance
to_plot <- cbind(Cooccurence, country_plot_merged_t)
# merge rows of the same cooccurrence status
to_plot_merged <- to_plot %>%
                  group_by_at("Cooccurence") %>%
                  summarise_at(c(1:ncol(to_plot)-1), sum)
# transpose and set names
to_plot_merged_t <- as.data.frame(t(to_plot_merged))
names(to_plot_merged_t) <- c("Negative", "Positive")
# remove first row ("Cooccurence" due to transposition)
to_plot_merged_t <- to_plot_merged_t[-1,]
# set columns to numeric
to_plot_merged_t$Positive <- as.numeric(as.character(to_plot_merged_t$Positive))
to_plot_merged_t$Negative <- as.numeric(as.character(to_plot_merged_t$Negative))
# sum rows to get total count per country
to_plot_merged_t$sum <- rowSums(to_plot_merged_t)
# get proportions of pos/neg for each country
to_plot_merged_t <- to_plot_merged_t/to_plot_merged_t$sum
# remove sum column
to_plot_merged_t <- to_plot_merged_t[,-3]
# create country vector
Country <- rownames(to_plot_merged_t)
# create new dataframe with country, postive and negative columns
for_melting <- cbind.data.frame(Country, to_plot_merged_t)
# transpose
for_melting_t <- as.data.frame(t(for_melting))
# # remove first row ("Country" due to transposition)
for_melting_t <- for_melting_t[-1,]
# add Cooccurrence column
for_melting_t$Cooccurence <- rownames(for_melting_t)

# gather into format required for plotting
melted_to_plot <- gather(for_melting_t, key="Country", value = "Abundance", names(for_melting_t)[1]:names(for_melting_t)[ncol(for_melting_t)-1])
# set country as factor
melted_to_plot$Country <- factor(melted_to_plot$Country, levels=unique(melted_to_plot$Country))
# set abundance to numeric
melted_to_plot$Abundance <- as.numeric(melted_to_plot$Abundance)

# set labels for plotting
labs <- c()
for (row in 1:nrow(no_samples)) {
  labs <- c(labs, paste0(no_samples[row,1], " (n = ", no_samples[row,2], ")"))
}

# plot bar chart showing proportion of pos/neg bd co-occurrence per country
pdf(paste0(path_out, "cooccur_by_country.pdf"))
ggplot(melted_to_plot, aes(fill=Cooccurence, y=Abundance, x=Country), cex = 6) + 
  geom_bar(position="stack", stat="identity") +
  scale_y_continuous(breaks=c(0.2,0.4,0.6,0.5,0.8,1.0)) +
  scale_x_discrete(labels = labs) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")
dev.off()


###########################################################################################
############################## Co-occurrence phylum bar Chart #############################


# subset by phylum
phylum <- as.data.frame(table(taxa_sig_merged$Phylum, taxa_sig_merged$Cooccurrence))
names(phylum) <- c("Phylum", "Cooccurrence", "Frequency")
# remove any phylum with a count equal to zero
phylum <- phylum[!phylum$Frequency == 0,]

# subset by positively and negatively co-occurring
neg <- phylum$Frequency[phylum$Cooccurrence == "Negative"]
pos <- phylum$Frequency[phylum$Cooccurrence == "Positive"]
negsum <- sum(neg)
possum <- sum(pos)
phylum$Frequency[phylum$Cooccurrence == "Negative"] <- neg/negsum
phylum$Frequency[phylum$Cooccurrence == "Positive"] <- pos/possum

# set up colour blind friendly palatte
cbpalette <- c("#6195CF", "#994F88", "#1965B0", "#90C987")

# plot bar chart comparing phylum composition of positive and negative cooccurrence taxa
pdf(paste0(path_out, "bar_phylum_Bd.pdf"))
print(ggplot(phylum, aes(x = Cooccurrence, y = Frequency, fill = Phylum)) +
        geom_bar(position="stack", stat="identity") +
        scale_fill_manual(values = cbpalette) +
        scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0,0))) +
        ylab("Proportion") +
        theme(text = element_text(size = 17),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.text = element_text(size = 16),   
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line.y = element_line(colour = "black")))  
dev.off()


###########################################################################################
################################# Co-occurrence Pie Charts #################################


# separate out postively and negatively cooccuring taxa
pos <- taxa_sig_merged[taxa_sig_merged$Cooccurrence == "Positive",]
pos <- pos[,-ncol(pos)]
neg <- taxa_sig_merged[taxa_sig_merged$Cooccurrence == "Negative",]
neg <- neg[,-ncol(neg)]
# save as list to use in for loop
taxa_sig_ls <- list(pos=pos, neg=neg)
taxa_ls <- c("Phylum", "Class", "Order")

# set list_index for pdf naming
list_index <- 1
for (item in taxa_sig_ls) {
  for (tax in taxa_ls) {
    list_index <- 1-list_index
    # subset df by taxa
    subs <- as.data.frame(table(item[[tax]]))
    names(subs) <- c(tax, "Frequency")
    # remove any taxa with 0 abundance
    subs <- subs[!subs$Frequency == 0,]
    # plot pie chart
    pdf(paste0(path_out, "pie_", tax, "_", if (list_index==0) "pos" else "neg", "withBd.pdf"))
    print(ggplot(subs, aes(x= "", y= Frequency, fill= get(tax))) +
          geom_bar(stat = "identity", width = 1) +
          coord_polar("y", start = 0) +
            labs(fill = tax) +
          theme_void())
    dev.off()
  }
}


###########################################################################################
################################### Cooccurring taxa ######################################


# combine taxa data with cooccurrence results
tax_p <- merge(taxa_sig_merged, cooccur_bd_fil, by.x = "Genus", by.y = "Genus")
tax_p <- tax_p[,-ncol(tax_p)]
# subset by positive and negative cooccurence
co_tax_pos <- subset(tax_p, tax_p$Cooccurrence == "Positive")
co_tax_neg <- subset(tax_p, tax_p$Cooccurrence == "Negative")
# negatively and positively associated taxa ordered by p-value
co_tax_neg_ord <- co_tax_neg[order(co_tax_neg$p_lt),]
co_tax_pos_ord <- co_tax_pos[order(co_tax_pos$p_gt),]
# negatively and positively associated taxa grouped by taxa
co_tax_neg_group <- co_tax_neg[order(co_tax_neg$Phylum, co_tax_neg$Class, co_tax_neg$Order, co_tax_neg$Family, co_tax_neg$Genus),]
co_tax_pos_group <- co_tax_pos[order(co_tax_pos$Phylum, co_tax_pos$Class, co_tax_pos$Order, co_tax_pos$Family, co_tax_pos$Genus),]
# save to csv
write.csv(co_tax_neg_ord, paste0(path_out, "negatively_cooccuring_genus.csv"))
write.csv(co_tax_pos_ord, paste0(path_out, "positively_cooccuring_genus.csv"))


## end of script
