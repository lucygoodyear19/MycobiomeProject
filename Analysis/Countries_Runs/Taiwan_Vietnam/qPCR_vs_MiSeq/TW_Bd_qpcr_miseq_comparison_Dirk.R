#######################################################################
###################### Taiwan Bd qPCR vs MiSeq ########################
#######################################################################

# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


############# initial set up, library and data loading ################


# set wd if needed
setwd("Documents/CMEECourseWork/Project/Data-scripts/Taiwan_Vietnam/qPCR_vs_MiSeq")

# load packages
library(dplyr)
library(ggplot2)
library(scales)
# load required functions
source("qPCR_vs_MiSeq_functions.R")

# load filtered Taiwan data (swabs only)
TWdata <- read.csv("../Filtering/Taiwanswabsonly_plate2-4_filteredOTUtable.csv", 
                   header=T, 
                   stringsAsFactor=FALSE, 
                   check.names = F)

# load plate data
plates <- read.csv("../TW16_plate_data.csv", 
                   stringsAsFactor=F, 
                   header=F)
names(plates) <- c("barcode","sample","DNAqual","plate")

# load qPCR data
qpcr <- read.csv("../TW_full_qpcr_data.csv", 
                 header=T, 
                 stringsAsFactor=FALSE, 
                 check.names = F)


###################### qpcr data filtering #######################


# remove underscores in qpcr sample names to match MiSeq data
qpcr$id <- gsub("_", "", qpcr$id)

# subset qpcr results to include only id and bd_results
qpcr_Bd <- qpcr %>% select(id, bd_results)


##################### miseq data filtering ######################


# filter OTU table by Bd only and remove meta data
miseq_Bd <- TWdata[TWdata$species == "Batrachochytrium_dendrobatidis",]
miseq_Bd <- miseq_Bd[,23:ncol(miseq_Bd)]


################# qpcr/miseq data preparation ###################


# filter Bd only OTU table by samples with qpcr measurement
samples_miseq_qpcr <- miseq_Bd[,names(miseq_Bd) %in% c(qpcr_Bd$id)]

# create new dataframe with overlapping ids as first column
samples_Bd <- data.frame(names(samples_miseq_qpcr))
# add columns for total Bd count per sample and qpcr Bd boolean (1-0)
samples_Bd <- cbind(samples_Bd, 
                    colSums(samples_miseq_qpcr), 
                    qpcr_Bd[qpcr_Bd$id %in% names(samples_miseq_qpcr),2])
# set column names
names(samples_Bd) <- c("id", "miseq", "qpcr")
samples_Bd$qpcr <- as.factor(samples_Bd$qpcr) # for plotting purposes


####################### Plotting data ########################


###### plotting overall data set

# remove all 0-0 matches for plotting
samples_Bd_plot <- samples_Bd[samples_Bd$qpcr != 0 | samples_Bd$miseq != 0,]
samples_Bd_plot$labels <- gsub("TW16", "", samples_Bd_plot$id)

# plot miseq vs qpcr
plot_all <- qplot(labels, miseq, data = samples_Bd_plot,
                  xlab = "Sample ID",
                  ylab = "MiSeq count",
                  main = paste("MiSeq and qPCR Bd reads by sample (n =", nrow(samples_Bd), ")"),
                  colour = qpcr) +
            theme_bw() +
            annotate(geom="text", x=55, y=1300, 
                     label=paste(percent(total_match(samples_Bd)), 
                                 "total matches\n",
                                 percent(total_match100(samples_Bd)), 
                                 "total matches if MiSeq reads < 100 are assumed 0\n",
                                 percent(qpcr_miseq_zero(samples_Bd)),
                                 "total qPCR-/MiSeq- matches (not displayed)\n",
                                 percent(qpcr_neg_miseq_pos(samples_Bd)/nrow(samples_Bd)), 
                                 "(",qpcr_neg_miseq_pos(samples_Bd), "samples)",
                                 "qPCR-/MiSeq+ for all reads\n",
                                 percent(qpcr_neg_miseq_pos100(samples_Bd)/nrow(samples_Bd)),
                                 "(",qpcr_neg_miseq_pos100(samples_Bd), "samples)",
                                 "qPCR-/MiSeq+ if MiSeq reads < 100 are assumed 0"),
                     color="black", cex = 4, hjust = 0) +
            scale_colour_discrete(name = "qPCR", labels = c("Bd-", "Bd+"))

###### plotting by plate

# subset data by plate
samples_Bd_plate2 <- samples_Bd[samples_Bd$id %in% c(plates$sample[plates$plate == "plate2"]),]
samples_Bd_plate3 <- samples_Bd[samples_Bd$id %in% c(plates$sample[plates$plate == "plate3"]),]
samples_Bd_plate4 <- samples_Bd[samples_Bd$id %in% c(plates$sample[plates$plate == "plate4"]),]

### plate 2

# remove all 0-0 matches for plotting
samples_Bd_plate2_plot <- samples_Bd_plate2[samples_Bd_plate2$qpcr != 0 | samples_Bd_plate2$miseq != 0,]
samples_Bd_plate2_plot$labels <- gsub("TW16", "", samples_Bd_plate2_plot$id)

# plot miseq vs qpcr by plate
plot_plate2 <- qplot(labels, miseq, data = samples_Bd_plate2_plot,
                     xlab = "Sample ID",
                     ylab = "MiSeq count",
                     main = paste("Plate 2: MiSeq and qPCR Bd reads by sample (n =", nrow(samples_Bd_plate2), ")"),
                     colour = qpcr) +
               theme_bw() +
               annotate(geom="text", x=5, y=1300, 
                        label=paste(percent(total_match(samples_Bd_plate2)), 
                                    "total matches\n",
                                    percent(total_match100(samples_Bd_plate2)), 
                                    "total matches if MiSeq reads < 100 are assumed 0\n",
                                    percent(qpcr_miseq_zero(samples_Bd_plate2)),
                                    "total qPCR-/MiSeq- matches (not displayed)\n",
                                    percent(qpcr_neg_miseq_pos(samples_Bd_plate2)/nrow(samples_Bd_plate2)), 
                                    "(",qpcr_neg_miseq_pos(samples_Bd_plate2), "samples)",
                                    "qPCR-/MiSeq+ for all reads\n",
                                    percent(qpcr_neg_miseq_pos100(samples_Bd_plate2)/nrow(samples_Bd_plate2)),
                                    "(",qpcr_neg_miseq_pos100(samples_Bd_plate2), "samples)",
                                    "qPCR-/MiSeq+ if MiSeq reads < 100 are assumed 0"),
                        color="black", cex = 4, hjust = 0) +
               scale_colour_discrete(name = "qPCR", labels = c("Bd-", "Bd+"))

### plate 3

# remove all 0-0 matches for plotting
samples_Bd_plate3_plot <- samples_Bd_plate3[samples_Bd_plate3$qpcr != 0 | samples_Bd_plate3$miseq != 0,]
samples_Bd_plate3_plot$labels <- gsub("TW16", "", samples_Bd_plate3_plot$id)

# plot miseq vs qpcr by plate
plot_plate3 <- qplot(labels, miseq, data = samples_Bd_plate3_plot,
                     xlab = "Sample ID",
                     ylab = "MiSeq count",
                     main = paste("Plate 3: MiSeq and qPCR Bd reads by sample (n =", nrow(samples_Bd_plate3), ")"),
                     colour = qpcr) +
               theme_bw() +
               annotate(geom="text", x=9, y=30, 
                        label=paste(percent(total_match(samples_Bd_plate3)), 
                                    "total matches\n",
                                    percent(total_match100(samples_Bd_plate3)), 
                                    "total matches if MiSeq reads < 100 are assumed 0 \n",
                                    percent(qpcr_miseq_zero(samples_Bd_plate3)),
                                    "total qPCR-/MiSeq- matches (not displayed)\n",
                                    percent(qpcr_neg_miseq_pos(samples_Bd_plate3)/nrow(samples_Bd_plate3)), 
                                    "(",qpcr_neg_miseq_pos(samples_Bd_plate3), "samples)",
                                    "qPCR-/MiSeq+ for all reads\n",
                                    percent(qpcr_neg_miseq_pos100(samples_Bd_plate3)/nrow(samples_Bd_plate3)),
                                    "(",qpcr_neg_miseq_pos100(samples_Bd_plate3), "samples)",
                                    "qPCR-/MiSeq+ if MiSeq reads < 100 are assumed 0"),
                        color="black", cex = 4, hjust = 0) +
               scale_colour_discrete(name = "qPCR", labels = c("Bd-", "Bd+"))

### plate 4

# remove all 0-0 matches for plotting
samples_Bd_plate4_plot <- samples_Bd_plate4[samples_Bd_plate4$qpcr != 0 | samples_Bd_plate4$miseq != 0,]
samples_Bd_plate4_plot$labels <- gsub("TW16", "", samples_Bd_plate4_plot$id)

# plot miseq vs qpcr by plate
plot_plate4 <- qplot(labels, miseq, data = samples_Bd_plate4_plot,
                     xlab = "Sample ID",
                     ylab = "MiSeq count",
                     main = paste("Plate 4: MiSeq and qPCR Bd reads by sample (n =", nrow(samples_Bd_plate4), ")"),
                     colour = qpcr) +
               theme_bw() +
               annotate(geom="text", x=5, y=450, 
                        label=paste(percent(total_match(samples_Bd_plate4)), 
                                    "total matches\n",
                                    percent(total_match100(samples_Bd_plate4)), 
                                    "total matches if MiSeq reads < 100 are assumed 0 \n",
                                    percent(qpcr_miseq_zero(samples_Bd_plate4)),
                                    "total qPCR-/MiSeq- matches (not displayed)\n",
                                    percent(qpcr_neg_miseq_pos(samples_Bd_plate4)/nrow(samples_Bd_plate4)),
                                    "(",qpcr_neg_miseq_pos(samples_Bd_plate4), "samples)",
                                    "qPCR-/MiSeq+ for all reads\n",
                                    percent(qpcr_neg_miseq_pos100(samples_Bd_plate4)/nrow(samples_Bd_plate4)),
                                    "(",qpcr_neg_miseq_pos100(samples_Bd_plate4), "samples)",
                                    "qPCR-/MiSeq+ if MiSeq reads < 100 are assumed 0"),
                        color="black", cex = 4, hjust = 0) +
               scale_colour_discrete(name = "qPCR", labels = c("Bd-", "Bd+"))


# save resulting plots as pdfs for review
pdf("qpcr_vs_miseq_all_Dirk.pdf", width = 14, height = 7, paper = "a4r")
plot_all
dev.off()

pdf("qpcr_vs_miseq_plate2_Dirk.pdf", width = 14, height = 7, paper = "a4r")
plot_plate2
dev.off()

pdf("qpcr_vs_miseq_plate3_Dirk.pdf", width = 14, height = 7, paper = "a4r")
plot_plate3
dev.off()

pdf("qpcr_vs_miseq_plate4_Dirk.pdf", width = 14, height = 7, paper = "a4r")
plot_plate4
dev.off()
