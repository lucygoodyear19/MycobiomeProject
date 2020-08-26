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
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("dplyr")
library("viridis")
library("microbiome") # for abundances()
library("tidyr")


# set wd
setwd("/Users/lucy/Documents/MRes/MycobiomeProject/Analysis/Results/Global_14/")

# load phyloseq object
dada2 <- readRDS("/Users/lucy/Documents/MRes/MycobiomeProject/Analysis/Results/Global_14/physeqob_DADA2.rds")

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


################################################################################################
##################################### Map preparation ##########################################


# extract unique location for each country
sites <- unique(samps[,c("Latitude", "Longitude", "Country")])
# set country as factor
sites$Country <- as.factor(sites$Country)
# count number of sites per country
no_sites_per_country <- table(sites$Country)
# count number of samples per site
no_samps_per_site <- samps %>% group_by(Country, Latitude, Longitude) %>% summarize(count=n())
no_samps_per_site <- as.data.frame(no_samps_per_site)
# count number of samples per country 
# first extract first long/lat value per country
single_samps <- data.frame(c())
for (country in unique(samps$Country)) {
  samps_sub <- samps %>% filter(samps$Country == country)
  samp_row <- samps_sub[1,]
  samp_row <- samp_row %>% select(Longitude, Latitude, Country)
  single_samps <- rbind.data.frame(single_samps, samp_row)
}
# find number of samples per country
no_samps_per_country <- samps %>% group_by(Country) %>% summarize(count=n())
no_samps_per_country <- as.data.frame(no_samps_per_country)
# add longitude and latutide data for plotting
samps_per_country_plot <- merge(no_samps_per_country, single_samps, by.x = "Country", by.y = "Country")
samps_per_country_plot$Longitude <- as.numeric(as.character(samps_per_country_plot$Longitude))
samps_per_country_plot$Latitude <- as.numeric(as.character(samps_per_country_plot$Latitude))
samps_per_country_plot$count <- as.numeric(as.character(samps_per_country_plot$count))

# remove Taiwan since sample number is much larger so plot separately
tw_plot <- samps_per_country_plot[samps_per_country_plot$Country == "Taiwan",]
no_tw_plot <- samps_per_country_plot[!samps_per_country_plot$Country == "Taiwan",]

# load world map data
world <- ne_countries(scale = "medium", returnclass = "sf")
#pops <- ne_download(scale = 110, type = 'populated_places', category = 'cultural', returnclass='sf')


###############################################################################################
#################################### Plot map of sites ########################################


pdf(paste0(results_path, "Summary/map_samples.pdf"))
ggplot(data = world) +
  geom_sf(colour = "gray44", size = 0.1, fill = NA) +
  geom_point(data = no_tw_plot, 
             aes(x = Longitude, y = Latitude, colour = count),
             size = 3, 
             shape = 20) +
  scale_color_viridis_c(direction=-1) +
  geom_point(data = tw_plot,
             aes(x = Longitude, y = Latitude, fill = "590"),
             size = 3,
             shape = 20) +
  labs(colour = "No. of Samples") +
  labs(fill = element_blank()) +
  theme(legend.title = element_text(size=8),
        legend.text = element_text(size=7)) +
  guides(colour = guide_colourbar(order = 1), 
         fill = guide_legend(order = 2)) +
  #coord_sf(xlim = c(118, 123), ylim = c(21, 26), expand = FALSE) +
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
        plot.background=element_blank())
dev.off()
        

#################################################################################################
############################### Plot map of alpha diversity #####################################


alpha_table <- as.data.frame(c())
for (country in unique(samps$Country)) {
  samps_sub <- samps %>% filter(samps$Country == country)
  subs <- subset(samps_sub, select = c("Country", "Alpha_Shannon"))
  subs$Alpha_Shannon <- as.numeric(as.character(subs$Alpha_Shannon))
  count <- mean(subs$Alpha_Shannon)
  alpha_table <- rbind.data.frame(alpha_table, cbind(country, count))
}

# rename columns
names(alpha_table) <- c("Country", "Alpha")

# merge with long/lat data for plotting
alpha_plot <- merge(alpha_table, single_samps, by.x = "Country", by.y = "Country")
alpha_plot$Alpha <- as.numeric(as.character(alpha_plot$Alpha))
alpha_plot$Longitude <- as.numeric(as.character(alpha_plot$Longitude))
alpha_plot$Latitude <- as.numeric(as.character(alpha_plot$Latitude))

# plot map
pdf(paste0(results_path, "Summary/map_alpha.pdf"))
ggplot(data = world) +
  geom_sf(colour = "gray44", size = 0.1, fill = NA) +
  #geom_sf(data = pops) +
  geom_point(data = alpha_plot, 
             aes(x = Longitude, y = Latitude, colour = Alpha),
             size = 3, 
             shape = 20,
             show.legend = TRUE) +
  scale_colour_viridis() +
  labs(colour = "Shannon Alpha") +
  theme(legend.title = element_text(size=8),
        legend.text = element_text(size=7)) +
  #coord_sf(xlim = c(118, 123), ylim = c(21, 26), expand = FALSE) +
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
        plot.background=element_blank())
dev.off()


#################################################################################################
############################### Plot map of ASV richness ########################################


richness_table <- as.data.frame(c())
for (country in unique(samps$Country)) {
  samps_sub <- samps %>% filter(samps$Country == country)
  seqs_sub <- seqs[rownames(seqs) %in% rownames(samps_sub)]
  seqs_sub$sum <- rowSums(seqs_sub)
  #seqs_sub <- seqs_sub[!seqs_sub$sum < 100,]
  count <- mean(seqs_sub$sum)
  richness_table <- rbind.data.frame(richness_table, cbind(country, count))
}

# rename columns
names(richness_table) <- c("Country", "ASV_Richness")

# merge with long/lat data for plotting
richness_plot <- merge(richness_table, single_samps, by.x = "Country", by.y = "Country")
richness_plot$ASV_Richness <- as.numeric(as.character(richness_plot$ASV_Richness))
richness_plot$Longitude <- as.numeric(as.character(richness_plot$Longitude))
richness_plot$Latitude <- as.numeric(as.character(richness_plot$Latitude))

# remove Cameroon since sample number is much larger so plot separately
cam_plot <- richness_plot[samps_per_country_plot$Country == "Cameroon",]
no_cam_plot <- richness_plot[!samps_per_country_plot$Country == "Cameroon",]

pdf(paste0(results_path, "Summary/map_richness.pdf"))
ggplot(data = world) +
  geom_sf(colour = "gray44", size = 0.1, fill = NA) +
  geom_point(data = no_cam_plot, 
             aes(x = Longitude, y = Latitude, colour = ASV_Richness),
             size = 3, 
             shape = 20) +
  scale_color_viridis_c(direction=-1) +
  geom_point(data = cam_plot,
             aes(x = Longitude, y = Latitude, fill = "12,200"),
             size = 3,
             shape = 20) +
  labs(colour = "ASV Richness") +
  labs(fill = element_blank()) +
  theme(legend.title = element_text(size=8),
        legend.text = element_text(size=7)) +
  guides(colour = guide_colourbar(order = 1), 
         fill = guide_legend(order = 2)) +
  #coord_sf(xlim = c(118, 123), ylim = c(21, 26), expand = FALSE) +
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
        plot.background=element_blank())
dev.off()

##################################################################################################
#################################### Pie charts of taxa ##########################################


phyla <- as.data.frame(table(taxa$Phylum))
names(phyla) <- c("Phyla", "Frequency")

# ggplot pie chart
pdf("Summary/phyla_pie_ggplot.pdf")
print(ggplot(phyla, aes(x= "", y= Frequency, fill= Phyla)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void())
dev.off()

# base R pie chart
slices <- phyla$Frequency
lbls <- phyla$Phyla
pct <- round(slices/sum(slices)*100)
lbls <- paste0(lbls, ", ", pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pdf("Summary/phyla_pie_base.pdf")
pie(slices, labels = NA, col=rainbow(length(lbls)))
    #main="Pie Chart of Countries") 
legend("topright", lbls, fill=rainbow(length(slices)), cex = 0.5)
dev.off()


#################################################################################################
############################# Taxa composition by presence/absence ##############################


# plot by presence absence
find_taxa_breakdown <- function(taxa_data, seq_data, samp_data, taxa_rank, condition, NAs) {
  taxon_table <- as.data.frame(c())
  labs <- c()
  factors <- c()
  for (var in unique(samp_data[[condition]])) {
    if (!is.na(var)) {
      samps_sub <- samp_data %>% filter(samp_data[[condition]] == var)
      labs <- c(labs, paste0(samps_sub[[condition]][1], " (n = ", nrow(samps_sub), ")"))
      factors <- c(factors, var)
      samp_seqs <- seq_data[colnames(seq_data) %in% rownames(samps_sub)]
      samp_seqs$sum <- rowSums(samp_seqs)
      samp_seqs <- samp_seqs[samp_seqs$sum > 0,]
  
      samp_taxa <- taxa_data[rownames(taxa_data) %in% rownames(samp_seqs),]
      taxon <- as.character(samp_taxa[[taxa_rank]])
      if (NAs == "NAinc")
      {taxon[is.na(taxon)] <- "NA"}
      count <- as.data.frame(table(taxon))
      count$taxon[count$taxon == "NA"] <- NA
      count$Freq <- count$Freq/sum(count$Freq)
      count$taxon <- gsub(".__", "", count$taxon)
      count$Variable <- var
      taxon_table <- rbind(taxon_table, count)
    }
  }
  taxon_table$Variable <- factor(taxon_table$Variable, ordered=TRUE, levels = factors)
  #is.na(taxon_table) <- 0
  # plot stacked bar chart for all countries
  pdf(paste0(results_path, "Summary/", condition, "_", NAs, "_", taxa_rank, "_breakdown.pdf"))
  print(ggplot(taxon_table, aes(fill=taxon, y=Freq, x=Variable)) + 
    geom_bar(position="stack", stat="identity") +
    labs(fill = taxa_rank) +
    ylab("Proportion") +
    xlab(condition) +
    scale_x_discrete(labels = labs) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
}

find_taxa_breakdown(taxa, seqs_pres_abs_t, samps, "Phylum", "Country", "NAinc")
find_taxa_breakdown(taxa, seqs_pres_abs_t, samps, "Phylum", "Country", "noNA")
find_taxa_breakdown(taxa, seqs_pres_abs_t, samps, "Phylum", "Bd", "NAinc")
find_taxa_breakdown(taxa, seqs_pres_abs_t, samps, "Phylum", "Bd", "noNA")
find_taxa_breakdown(taxa, seqs_pres_abs_t, samps, "Phylum", "A_Order", "NAinc")
find_taxa_breakdown(taxa, seqs_pres_abs_t, samps, "Phylum", "A_Order", "noNA")
samps$Lifestage <- as.character(samps$Lifestage)
samps$Lifestage[samps$Lifestage == ""] <- NA
samps$Lifestage <- as.factor(samps$Lifestage)
find_taxa_breakdown(taxa, seqs_pres_abs_t, samps, "Phylum", "Lifestage", "NAinc")
find_taxa_breakdown(taxa, seqs_pres_abs_t, samps, "Phylum", "Lifestage", "noNA")

find_taxa_breakdown(taxa, seqs_pres_abs_t, samps, "Class", "Country", "NAinc")
find_taxa_breakdown(taxa, seqs_pres_abs_t, samps, "Class", "Country", "noNA")
find_taxa_breakdown(taxa, seqs_pres_abs_t, samps, "Class", "Bd", "NAinc")
find_taxa_breakdown(taxa, seqs_pres_abs_t, samps, "Class", "Bd", "noNA")
find_taxa_breakdown(taxa, seqs_pres_abs_t, samps, "Class", "A_Order", "NAinc")
find_taxa_breakdown(taxa, seqs_pres_abs_t, samps, "Class", "A_Order", "noNA")
samps$Lifestage <- as.character(samps$Lifestage)
samps$Lifestage[samps$Lifestage == ""] <- NA
samps$Lifestage <- as.factor(samps$Lifestage)
find_taxa_breakdown(taxa, seqs_pres_abs_t, samps, "Class", "Lifestage", "NAinc")
find_taxa_breakdown(taxa, seqs_pres_abs_t, samps, "Class", "Lifestage", "noNA")


#################################################################################################
################################ Taxa composition by abundance ##################################


seqs_t <- as.data.frame(t(seqs))

find_taxa_breakdown_abun  <- function(taxa_data, seq_data, samp_data, taxa_rank, condition) {
  taxon_table <- as.data.frame(c())
  labs <- c()
  factors <- c()
  for (var in unique(samp_data[[condition]])) {
    if (!is.na(var)) {
      samps_sub <- samp_data %>% filter(samp_data[[condition]] == var)
      labs <- c(labs, paste0(samps_sub[[condition]][1], " (n = ", nrow(samps_sub), ")"))
      factors <- c(factors, var)
      samp_seqs <- seq_data[colnames(seq_data) %in% rownames(samps_sub)]
      samp_seqs$sum <- rowSums(samp_seqs)
      samp_seqs <- samp_seqs[samp_seqs$sum > 0,]
      total <- sum(samp_seqs$sum)
    
      samp_taxa <- taxa_data[rownames(taxa_data) %in% rownames(samp_seqs),]
      taxon <- samp_taxa[[taxa_rank]]
      pcts <- samp_seqs$sum/total
      count <- cbind.data.frame(taxon, pcts)

      require(dplyr)
      count <- count %>%
               group_by(taxon) %>%
               summarise(Percentage = sum(pcts))
    
      count$Variable <- var
      count$taxon <- gsub(".__", "", count$taxon)
  
      taxon_table <- rbind.data.frame(taxon_table, count)
    }
  }
  taxon_table$taxon <- as.character(taxon_table$taxon)
  taxon_table$taxon[is.na(taxon_table$taxon)] <- "NA"
  taxon_table$taxon <- as.factor(taxon_table$taxon)
  taxon_table$taxon <- fct_relevel(taxon_table$taxon, "NA", after = Inf)
  taxon_table$Variable <- factor(taxon_table$Variable, ordered=TRUE, levels = factors)
  write.csv(taxon_table, paste0("Summary/", condition, "_", taxa_rank, "breakdown.csv"))
  #is.na(taxon_table) <- 0
  # plot stacked bar chart for all countries
  pdf(paste0(results_path, "Summary/", condition, "_", taxa_rank, "_abun_breakdown.pdf"))
  print(ggplot(taxon_table, aes(fill=taxon, y=Percentage, x=Variable)) + 
    geom_bar(position="stack", stat="identity") +
    coord_fixed(ratio=18) +
    scale_fill_manual(values = cbpalette) +
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

find_taxa_breakdown_abun(taxa, seqs_t, samps, "Phylum", "Country")
find_taxa_breakdown_abun(taxa, seqs_t, samps, "Phylum", "Bd")
find_taxa_breakdown_abun(taxa, seqs_t, samps, "Phylum", "A_Order")
samps$Lifestage <- as.character(samps$Lifestage)
samps$Lifestage[samps$Lifestage == ""] <- NA
samps$Lifestage <- as.factor(samps$Lifestage)
find_taxa_breakdown_abun(taxa, seqs_t, samps, "Phylum", "Lifestage")
find_taxa_breakdown_abun(taxa, seqs_t, samps, "Phylum", "A_Family")

# total abundance breakdown by phyla
samp_seqs <- seqs_t
samp_seqs$sum <- rowSums(samp_seqs)
total <- sum(samp_seqs$sum)

taxon <- taxa[["Phylum"]]

pcts <- samp_seqs$sum/total
count <- cbind.data.frame(taxon, pcts)

count <- count %>%
  group_by(taxon) %>%
  summarise(Percentage = sum(pcts))

count$taxon <- gsub(".__", "", count$taxon)

write.csv(taxon_table, paste0("Summary/total_phylum_breakdown.csv"))



#################################################################################################
######################################### Bd status #############################################


# extract samples with Bd not NA
samps_bd <- samps[!is.na(samps$Bd),]

# plot stacked bar chart for Bd status proportion per country
bd_table <- as.data.frame(c())
labs <- c() # for plot labels
factors <- c() # for correct label order
for (country in unique(samps_bd$Country)) {
  samps_sub <- samps_bd %>% filter(samps_bd$Country == country)
  samps_sub$Country <- as.character(samps_sub$Country)
  labs <- c(labs, paste0(country, " (n = ", nrow(samps_sub), ")"))
  factors <- c(factors, country)
  subs <- subset(samps_sub, select = c("Country", "Bd"))
  count <- as.data.frame(table(subs, exclude = NULL))
  count$Freq <- count$Freq/sum(count$Freq)
  bd_table <- rbind(bd_table, count)
}
bd_table$Country <- factor(bd_table$Country, ordered=TRUE, levels = factors)
# plot stacked bar chart for all countries
pdf(paste0(results_path, "Summary/bd_breakdown_noNA.pdf"))
ggplot(bd_table, aes(fill=Bd, y=Freq, x=Country)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = cbpalette) +
  scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0,0))) +
  labs(fill = "Bd Status") +
  scale_x_discrete(labels = labs) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text = element_text(size = 16),   
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(colour = "black"))
dev.off()


# plot stacked bar chart for Bd status per country, including NA
bd_table <- as.data.frame(c())
labs <- c() # for plot labels
factors <- c() # for correct label order
for (country in unique(samps$Country)) {
  samps_sub <- samps %>% filter(samps$Country == country)
  samps_sub$Country <- as.character(samps_sub$Country)
  labs <- c(labs, paste0(country, " (n = ", nrow(samps_sub), ")"))
  factors <- c(factors, country)
  subs <- subset(samps_sub, select = c("Country", "Bd"))
  count <- as.data.frame(table(subs, exclude = NULL))
  count$Freq <- count$Freq/sum(count$Freq)
  bd_table <- rbind(bd_table, count)
}
bd_table$Bd <- as.character(bd_table$Bd)
bd_table$Bd[is.na(bd_table$Bd)] <- "NA"
bd_table$Bd <- as.factor(bd_table$Bd)
bd_table$Country <- factor(bd_table$Country, ordered=TRUE, levels = factors)
# plot stacked bar chart for all countries
pdf(paste0(results_path, "Summary/bd_breakdown.pdf"))
ggplot(bd_table, aes(fill=Bd, y=Freq, x=Country)) + 
  geom_bar(position="stack", stat="identity") +
  coord_fixed(ratio=5) +
  scale_fill_manual(values = c("#90C987", "#6195CF", "#F6C141")) +
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
  
  pdf(paste0(results_path, "Summary/", tax_rank, "_pie_bd", condition, ".pdf"))
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

pdf(paste0(results_path, "Summary/pie_bd.pdf"))
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

samps_A_Order_all <- samps
samps_A_Order_all$A_Order <- as.character(samps_A_Order_all$A_Order)
pie_A_Order <- as.data.frame(table(samps_A_Order_all$A_Order))
names(pie_A_Order) <- c("A_OrderStatus", "Frequency")

pdf(paste0(results_path, "Summary/pie_A_Order.pdf"))
print(ggplot(pie_A_Order, aes(x= "", y= Frequency, fill= A_OrderStatus)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        scale_fill_viridis_d() +
        theme(legend.title = "Order") +
        #guides(fill=guide_legend(title=tax_rank)) +
        theme_void()) 
dev.off()

samps_Lifestage_all <- samps
samps_Lifestage_all$Lifestage <- as.character(samps_Lifestage_all$Lifestage)
samps_Lifestage_all$Lifestage[is.na(samps_Lifestage_all$Lifestage)] <- "NA"
pie_Lifestage <- as.data.frame(table(samps_Lifestage_all$Lifestage))
names(pie_Lifestage) <- c("LifestageStatus", "Frequency")

pdf(paste0(results_path, "Summary/pie_Lifestage.pdf"))
print(ggplot(pie_Lifestage, aes(x= "", y= Frequency, fill= LifestageStatus)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        scale_fill_viridis_d() +
        theme(legend.title = "Lifestage") +
        #guides(fill=guide_legend(title=tax_rank)) +
        theme_void()) 
dev.off()

samps_A_Family_all <- samps
samps_A_Family_all$A_Family <- as.character(samps_A_Family_all$A_Family)
samps_A_Family_all$A_Family[is.na(samps_A_Family_all$A_Family)] <- "NA"
pie_A_Family <- as.data.frame(table(samps_A_Family_all$A_Family))
names(pie_A_Family) <- c("A_FamilyStatus", "Frequency")

pdf(paste0(results_path, "Summary/pie_A_Family.pdf"))
print(ggplot(pie_A_Family, aes(x= "", y= Frequency, fill= A_FamilyStatus)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        scale_fill_viridis_d() +
        theme(legend.title = "Family") +
        #guides(fill=guide_legend(title=tax_rank)) +
        theme_void()) 
dev.off()


###################################################################################################
############################# Venn Diagram of shared ASVs/species #################################


library("VennDiagram")


######################################## functions ###############################################


# function to draw venn diagram
draw_venn <- function(sets, cat_names, file_prefix) {
  venn.diagram(x = sets,
               category.names = cat_names,
               filename = paste0(results_path, "Summary/", file_prefix, "_global_venn_diagram.png"), 
               output=TRUE)
}

# function to extract vectors of ASVs according to Bd status
find_asvs_for_venn <- function(samp_data, seq_data, condition, condition_result) {
  samps_sub <- samp_data %>% filter(samp_data[[condition]] == condition_result)
  cond_asv <- seq_data %>% select(colnames(seq_data[,colnames(seq_data) %in% rownames(samps_sub)]))
  cond_asv$Sum <- rowSums(cond_asv)
  cond_asv <- cond_asv %>% filter(cond_asv$Sum > 0)
  cond <- rownames(cond_asv)
}

# function to extract vectors of ASVs according to Bd status
# function preparation
taxa_species <- taxa %>% filter(!is.na(taxa$Species)) # subset taxa to species assigned ASVs only
# subset seqs to species assigned ASVs only
species_seqs <- seqs_t %>% filter(rownames(seqs_t) %in% rownames(taxa_species))
Species <- taxa_species$Species # save species names as vector
Genus <- taxa_species$Genus # save species names as vector
species_df <- cbind(Genus, Species, species_seqs) # create dataframe of seqs and species names

find_specs_for_venn <- function(samp_data, species_data, condition, condition_result) {
  samps_sub <- samp_data %>% filter(samp_data[[condition]] == condition_result)
  species_cond <- species_data %>% select(colnames(species_data[,colnames(species_data) %in% rownames(samps_sub)]))
  
  species_cond$Sum <- rowSums(species_cond[,-2])
  species_cond$Species <- paste0(Genus, " ", Species)
  species_cond <- species_cond %>% filter(species_cond$Sum > 0)
  specls <- c(as.character(unique(species_cond$Species)))
}


############################### extract results and plot venns ###################################


# continent
conts_samps <- samps
for (row in 1:nrow(conts_samps)) {
  if (conts_samps$Country[row] %in% c("China", "Kazakhstan", "Mongolia",
                                "North Korea", "South Korea", "Russia",
                                "Taiwan", "Vietnam")) {
    conts_samps$Continent[row] <- "Asia"
  } else {conts_samps$Continent[row] <- "non-Asia"}
}
table(conts_samps$Continent)

# random sample
set.seed(26)
asia_samp <- conts_samps[conts_samps$Continent == "Asia",]
noasia_samp <- conts_samps[conts_samps$Continent == "non-Asia",]

rand_samps_cont <- replicate(n = 100, sample(rownames(asia_samp), nrow(noasia_samp)), simplify = T)
asia_totals <- c()
noasia_totals <- c()
com_totals_cont <- c()
for (col in 1:ncol(rand_samps_cont)) {
  sampdata_rand_cont <- conts_samps[rownames(conts_samps) %in% rand_samps_cont[,col],]
  random_samples_all_cont <- rbind.data.frame(sampdata_rand_cont, noasia_samp)
  random_seqs_all_cont <- seqs_t[,colnames(seqs_t) %in% rownames(random_samples_all_cont)]
  # plot Venns for random sample
  asia_rand <- find_asvs_for_venn(random_samples_all_cont, random_seqs_all_cont, "Continent", "Asia")
  noasia_rand <- find_asvs_for_venn(random_samples_all_cont, random_seqs_all_cont, "Continent", "non-Asia")
  asia_tot <- length(asia_rand)
  noasia_tot <- length(noasia_rand)
  common_type_cont <- Reduce(intersect, list(asia_rand, noasia_rand))
  com_tot_cont <- length(common_type_cont)
  asia_totals <- c(asia_totals, asia_tot)
  noasia_totals <- c(noasia_totals, noasia_tot)
  com_totals_cont <- c(com_totals_cont, com_tot_cont)
}
com_av_cont <- round(mean(com_totals_cont), 0)
asia_av <- round(mean(asia_totals), 0)
noasia_av <- round(mean(noasia_totals), 0)




# Bd
bdpos <- find_asvs_for_venn(samps, seqs_t, "Bd", 1)
bdneg <- find_asvs_for_venn(samps, seqs_t, "Bd", 0)
draw_venn(list(bdpos, bdneg), c("Bd positive", "Bd negative"), "bd_asv")

# random sample
set.seed(26)
bd_sub <- samps[!is.na(samps$Bd),]
pos_samp <- bd_sub[bd_sub$Bd == 1,]
neg_samp <- bd_sub[bd_sub$Bd == 0,]

rand_samps_bd <- replicate(n = 100, sample(rownames(neg_samp), nrow(pos_samp)), simplify = T)
pos_totals <- c()
neg_totals <- c()
com_totals_bd <- c()
for (col in 1:ncol(rand_samps_bd)) {
  sampdata_rand_bd <- bd_sub[rownames(bd_sub) %in% rand_samps_bd[,col],]
  random_samples_all_bd <- rbind.data.frame(sampdata_rand_bd, pos_samp)
  random_seqs_all_bd <- seqs_t[,colnames(seqs_t) %in% rownames(random_samples_all_bd)]
  # plot Venns for random sample
  pos_rand <- find_asvs_for_venn(random_samples_all_bd, random_seqs_all_bd, "Bd", 1)
  neg_rand <- find_asvs_for_venn(random_samples_all_bd, random_seqs_all_bd, "Bd", 0)
  pos_tot <- length(pos_rand)
  neg_tot <- length(neg_rand)
  common_type_bd <- Reduce(intersect, list(pos_rand, neg_rand))
  com_tot_bd <- length(common_type_bd)
  pos_totals <- c(pos_totals, pos_tot)
  neg_totals <- c(neg_totals, neg_tot)
  com_totals_bd <- c(com_totals_bd, com_tot_bd)
}
com_av_bd <- round(mean(com_totals_bd), 0)
pos_av <- round(mean(pos_totals), 0)
neg_av <- round(mean(neg_totals), 0)


# Venn diagram for assigned species
#bd_pos_spec <- find_specs_for_venn(samps, species_df, "Bd", 1)
#bd_neg_spec <- find_specs_for_venn(samps, species_df, "Bd", 0)
#draw_venn(list(bd_pos_spec, bd_neg_spec), c("Bd positive", "Bd negative"), "bd_species")

# Lifestage
lifestage_adult <- find_asvs_for_venn(samps, seqs_t, "Lifestage", "Adult")
lifestage_tadpole <- find_asvs_for_venn(samps, seqs_t, "Lifestage", "Tadpole")
draw_venn(list(lifestage_adult, lifestage_tadpole), c("Adult", "Tadpole"), "lifestage_asv")
table(samps$Lifestage)

# random sample
set.seed(26)
lifestage_sub <- samps[!is.na(samps$Lifestage),]
adult_samp <- lifestage_sub[lifestage_sub$Lifestage == "Adult",]
tad_samp <- lifestage_sub[lifestage_sub$Lifestage == "Tadpole",]

rand_samps_life <- replicate(n = 100, sample(rownames(adult_samp), nrow(tad_samp)), simplify = T)
ad_totals <- c()
ta_totals <- c()
com_totals_life <- c()
for (col in 1:ncol(rand_samps_life)) {
  sampdata_rand <- lifestage_sub[rownames(lifestage_sub) %in% rand_samps_life[,col],]
  random_samples_all <- rbind.data.frame(sampdata_rand, tad_samp)
  random_seqs_all <- seqs_t[,colnames(seqs_t) %in% rownames(random_samples_all)]
  # plot Venns for random sample
  adult_rand <- find_asvs_for_venn(random_samples_all, random_seqs_all, "Lifestage", "Adult")
  tadpole_rand <- find_asvs_for_venn(random_samples_all, random_seqs_all, "Lifestage", "Tadpole")
  adult_tot <- length(adult_rand)
  tad_tot <- length(tadpole_rand)
  common_type_life <- Reduce(intersect, list(adult_rand, tadpole_rand))
  com_tot_life <- length(common_type_life)
  ad_totals <- c(ad_totals, adult_tot)
  ta_totals <- c(ta_totals, tad_tot)
  com_totals_life <- c(com_totals_life, com_tot_life)
}
com_av_life <- round(mean(com_totals_life), 0)
ad_av <- round(mean(ad_totals), 0)
ta_av <- round(mean(ta_totals), 0)


# Venn diagram for assigned species
#lifestage_tadpole_spec <- find_specs_for_venn(samps, species_df, "Lifestage", "Tadpole")
#lifestage_adult_spec <- find_specs_for_venn(samps, species_df, "Lifestage", "Adult")
#draw_venn(list(lifestage_adult_spec, lifestage_tadpole_spec), c("Adult", "Tadpole"), "lifestage_species")

# amphibian order
# subset by countries that have both caudata and anura samples
samps_caudata_countries <- unique(as.character(samps$Country[samps$A_Order == "Caudata"]))
samps_subs_order <- samps[samps$Country %in% samps_caudata_countries,]
seqs_subs_order <- seqs_t[,colnames(seqs_t) %in% rownames(samps_subs_order)]
# plot Venns
type_anura <- find_asvs_for_venn(samps_subs_order, seqs_subs_order, "A_Order", "Anura")
type_caudata <- find_asvs_for_venn(samps_subs_order, seqs_subs_order, "A_Order", "Caudata")
draw_venn(list(type_anura, type_caudata), c("Anura", "Caudata"), "type_asv")
table(samps_subs_order$A_Order)

# random sample
set.seed(26)
anurans_samp <- samps_subs_order[samps_subs_order$A_Order == "Anura",]
caudata_samp <- samps_subs_order[samps_subs_order$A_Order == "Caudata",]

rand_samps <- replicate(n = 100, sample(rownames(anurans_samp), nrow(caudata_samp)), simplify = T)
c_totals <- c()
a_totals <- c()
com_totals <- c()
for (col in 1:ncol(rand_samps)) {
  sampdata_rand <- samps_subs_order[rownames(samps_subs_order) %in% rand_samps[,col],]
  random_samples_all <- rbind.data.frame(sampdata_rand, caudata_samp)
  random_seqs_all <- seqs_subs_order[,colnames(seqs_subs_order) %in% rownames(random_samples_all)]
  # plot Venns for random sample
  type_anura_rand <- find_asvs_for_venn(random_samples_all, random_seqs_all, "A_Order", "Anura")
  type_caudata_rand <- find_asvs_for_venn(random_samples_all, random_seqs_all, "A_Order", "Caudata")
  caud_tot <- length(type_caudata_rand)
  anur_tot <- length(type_anura_rand)
  common_type <- Reduce(intersect, list(type_anura_rand, type_caudata_rand))
  com_tot <- length(common_type)
  c_totals <- c(c_totals, caud_tot)
  a_totals <- c(a_totals, anur_tot)
  com_totals <- c(com_totals, com_tot)
}
com_av <- round(mean(com_totals), 0)
c_av <- round(mean(c_totals), 0)
a_av <- round(mean(a_totals), 0)



venn.diagram(x = sets,
              category.names = cat_names,
              filename = paste0(results_path, "Summary/", file_prefix, "_global_venn_diagram.png"), 
              output=TRUE)
  
draw_venn(list(type_anura_rand, type_caudata_rand), c("Anura", "Caudata"), "randtype_asv")
table(random_samples_all$A_Order)


# Venn diagram for assigned species
#type_anura_spec <- find_specs_for_venn(samps, species_df, "A_Order", "Anura")
#type_caudata_spec <- find_specs_for_venn(samps, species_df, "A_Order", "Caudata")
#draw_venn(list(type_anura_spec, type_caudata_spec), c("Anura", "Caudata"), "type_species")


################# country 


############# ASVs


# create list of ASVs per country
country_prop_asv <- list()
for (country in samps$Country) {
  sub <- find_asvs_for_venn(samps, seqs_t, "Country", country)
  country_prop_asv[[country]] <- sub
}
# find intersecting ASVs for all countries
common_asvs <- Reduce(intersect, country_prop_asv)
# find taxa assignment for each common ASV
common_asvs_species <- taxa[rownames(taxa) %in% common_asvs,]
common_asvs_species_out <- capture.output(common_asvs_species)
# find proportion of ASVs common to all samples
no_common_asvs <- capture.output(nrow(common_asvs_species))
prop_asv_total <- capture.output(nrow(common_asvs_species)/nrow(taxa)) 

#taxa_shared_asv <- taxa[rownames(taxa) == common_asvs_species,]
seqs_shared_asv <- seqs[,colnames(seqs) %in% rownames(common_asvs_species)]
seqs_shared_asv$shared_sum <- rowSums(seqs_shared_asv)
per_total_asv <- capture.output(round(sum(seqs_shared_asv$shared_sum)/sum(seqs$total_sum), 4))

# find number of samples containing the common ASVs
props_asv <- list()
for (asv in (1:(ncol(seqs_shared_asv)-1))) {
  asv_sum <- sum(seqs_shared_asv[,asv])
  tot <- asv_sum/sum(seqs$total_sum)
  props_asv[[colnames(seqs_shared_asv)[asv]]] <- tot
}


############ Species


# create list of assigned species per country
country_prop_specs <- list()
for (country in samps$Country) {
  sub <- find_specs_for_venn(samps, species_df, "Country", country)
  country_prop_specs[[country]] <- sub
}
# find intersecting species for all countries
common_specs <- Reduce(intersect, country_prop_specs)
common_specs_out <- capture.output(common_specs)
# find proportion of assigned species common to all samples
taxa$full <- paste0(taxa$Genus, " ", taxa$Species)
common_species_asvs <- taxa[taxa$full %in% common_specs,]
no_common_specs <- capture.output(length(common_specs))
no_common_specs_asvs <- capture.output(nrow(common_species_asvs))
# percentage of common species out of all assigned species 
prop_spec_total <- capture.output(length(common_specs)/nrow(taxa_species))
# percentage of commmon species ASVs out of all ASVs
prop_spec_asv_total <- capture.output(nrow(common_species_asvs)/nrow(taxa))

# find average abundance of common species as a percentage of total number of reads
#seqs$total_sum <- rowSums(seqs)
#taxa_noNA <- taxa[!is.na(taxa$Species),]
taxa_shared <- taxa[taxa$full %in% common_specs,]
seqs_shared <- seqs[,colnames(seqs) %in% rownames(taxa_shared)]
seqs_shared$shared_sum <- rowSums(seqs_shared)
per_total <- capture.output(round(sum(seqs_shared$shared_sum)/sum(seqs$total_sum), 4))

# save output to text file
cat(paste0("Number of common ASVs: ", no_common_asvs, 
           "\nProportion of all ASVs that are common to all countries: ", prop_asv_total,
           "\nProportion of total reads for common ASVs: ", per_total_asv,
           "\nProportion of total reads for each common asv: ", props_asv,
           "\nTaxa assignment for common ASVs: "), common_asvs_species_out,
           file = paste0(results_path, "Summary/Common_taxa_summary_data.txt"),
                  sep = "\n", 
                  append=TRUE)
cat(paste0("\nNumber of common assigned species: ", no_common_specs,
           "\nProportion of all assigned species that are common to all countries: ", prop_spec_total,
           "\nProportion of ASVs that constitute common assigned species:", prop_spec_asv_total,
           "\nProportion of total reads for common species:", per_total,
           "\nCommon assigned species: "), common_specs_out,
           file = paste0(results_path, "Summary/Common_taxa_summary_data.txt"),
                  sep = "\n", 
                  append=TRUE)


##########################################################################################
################################### Most abundant ASV ####################################


# find the most abundant ASV in terms of presence/absence
seqs_0_df <- as.data.frame(seqs_pres_abs_t)
seqs_0_df$sum <- rowSums(seqs_0_df)
most_present_asvs <- rownames(seqs_0_df)[seqs_0_df$sum == max(seqs_0_df$sum)]
no_of_samples_max_present <- max(seqs_0_df$sum)
# find number of samples highest abundant taxa is present in
most_pres_all_samples <- seqs[,colnames(seqs) == most_present_asvs,drop=F]
most_pres_samps <- most_pres_all_samples[!most_pres_all_samples == 0, , drop=F]
# view breakdown of metadata for samples containing most abundant pres/abs ASV
samps_bd_most_pres <- samps[rownames(samps) %in% rownames(most_pres_samps),]
samps_bd_most_pres$Bd <- as.character(samps_bd_most_pres$Bd)
samps_bd_most_pres$Bd[is.na(samps_bd_most_pres$Bd)] <- "NA"
most_pres_bd <- capture.output(as.data.frame(table(samps_bd_most_pres$Bd)))
most_pres_country <- capture.output(as.data.frame(table(samps_bd_most_pres$Country)))
most_pres_order <- capture.output(as.data.frame(table(samps_bd_most_pres$A_Order)))
most_pres_life <- capture.output(as.data.frame(table(samps_bd_most_pres$Lifestage)))

cat(paste0("\nMost abundant ASVs (presence/absence): ", most_present_asvs,
           "\nNumber of samples most abundant ASV (presence/absence) is present in: ", nrow(most_pres_samps),
           "\nSpread of samples that contain most abundant ASV (presence/absence): ",
           most_pres_bd, "\n",
           most_pres_country, "\n"),
           most_pres_order,
    file = paste0(results_path, "Summary/Most_abundant_ASV_summary_data.txt"),
    sep = "\n", 
    append=TRUE)

# find the most abundant asv in terms of numbers of reads
seqs_df_t <- as.data.frame(seqs_t)
seqs_df_t$sum <- rowSums(seqs_df_t)
sum_reads <- sum(seqs_df_t$sum)
most_abun_asv <- rownames(seqs_df_t)[seqs_df_t$sum == max(seqs_df_t$sum)]
highest_abun <- max(seqs_df_t$sum)
highest_abun_per <- highest_abun/sum_reads
# find highest abundant taxa
highest_abun_taxa <- taxa[rownames(taxa) == most_abun_asv,]
# find number of samples highest abundant taxa is present in
most_abun_all_samples <- seqs[,colnames(seqs) == most_abun_asv,drop=F]
most_abun_samps <- most_abun_all_samples[!most_abun_all_samples == 0, , drop=F]
# view breakdown of metadata for samples containing most abundant ASV
samps_bd_most_abun <- samps[rownames(samps) %in% rownames(most_abun_samps),]
samps_bd_most_abun$Bd <- as.character(samps_bd_most_abun$Bd)
samps_bd_most_abun$Bd[is.na(samps_bd_most_abun$Bd)] <- "NA"
most_abun_bd <- capture.output(table(samps_bd_most_abun$Bd))
most_abun_country <- capture.output(table(samps_bd_most_abun$Country))
most_abun_order <- capture.output(table(samps_bd_most_abun$A_Order))

cat(paste0("\nMost abundant ASVs: ", most_abun_asv,
           "\nTaxa of most abundant ASV: ", highest_abun_taxa,
           "\nNumber of samples most abundant ASV is present in: ", nrow(most_abun_samps),
           "\nProportion of total reads: ", highest_abun_per,
           "\nSpread of samples that contain most abundant ASV: ",
           most_abun_bd, "\n",
           most_abun_country, "\n"),
    most_abun_order,
    file = paste0(results_path, "Summary/Most_abundant_ASV_summary_data.txt"),
    sep = "\n", 
    append=TRUE)


########################################################################################
################################# Summary stats ########################################


# summarise data collected
samps$Bd <- as.character(samps$Bd)
samps$Bd[is.na(samps$Bd)] <- "NA"
table(samps$Bd)

samps$Country <- as.character(samps$Country)
samps$Country[is.na(samps$Country)] <- "NA"
table(samps$Country)

samps$A_Order <- as.character(samps$A_Order)
samps$A_Order[is.na(samps$A_Order)] <- "NA"
table(samps$A_Order)

samps$Lifestage <- as.character(samps$Lifestage)
samps$Lifestage[samps$Lifestage == ""] <- NA
samps$Lifestage[is.na(samps$Lifestage)] <- "NA"
table(samps$Lifestage)

samps$Elevation_m <- as.character(samps$Elevation_m)
samps$Elevation_m[is.na(samps$Elevation_m)] <- "NA"
table(samps$Elevation_m)

samps$A_Family <- as.character(samps$A_Family)
samps$A_Family[is.na(samps$A_Family)] <- "NA"
gs_tab <- table(samps$A_Family)

samps$A_Genus_Species <- as.character(samps$A_Genus_Species)
samps$A_Genus_Species[is.na(samps$A_Genus_Species)] <- "NA"
gs_tab <- table(samps$A_Genus_Species)

# summarise phyla breakdown in table
phyla <- as.character(taxa$Phylum)
phyla_abun <- cbind.data.frame(phyla, seqs_pres_abs_t)
# merge by samples by genus
phyla_grid <- plyr::ddply(phyla_abun, "phyla", plyr::numcolwise(sum))

neo <- as.data.frame(phyla_grid[phyla_grid$Phylum == "Neocallimastigomycota",])

# remove NA row
phyla_grid <- subset(phyla_grid, !is.na(phyla))
# rename rows to genus 
rownames(phyla_grid) <- phyla_grid$phyla
#remove genus column
phyla_grid <- phyla_grid[,-1]
phyla_grid[phyla_grid > 0] <- 1
phyla_grid$sum <- rowSums(phyla_grid)
Sums <- phyla_grid$sum
Phyla <- rownames(phyla_grid)
phyla_sums <- cbind.data.frame(Phyla, Sums)

phylum_table <- as.data.frame(table(taxa$Phylum))
names(phylum_table) <- c("Phylum", "Number of ASVs")
write.csv(phylum_table, paste0(results_path, "Summary/phylum_breakdown.csv"))


## end of script
