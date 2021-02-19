##############################################################################
###################### Prelimary Taiwan Data Analysis ########################
##############################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


################################## Set up ####################################


# set wd if needed
#setwd("~/Documents/MRes/MycobiomeProject/Analysis/Runs_Countries/Taiwan_Vietnam_2016/Taiwan")

# load packages
library(dplyr)
library(phyloseq)
#library(ggplot2)
#library(gridExtra)

# load filtered phyloseq object for both Estonian pipeline and DADA2 pipeline
esto <- readRDS("physeqob_jen.rds")
#dada2 <- readRDS("physeqob_dada2.rds")


######################## Plot alpha-diversity ############################


# estimate value for alpha for each sample

# esto
shannon_esto <- estimate_richness(esto, split = TRUE, measures = "Shannon")
sample_data(esto)$alpha <- as.numeric(shannon_esto[,1])

# dada2
#shannon_dada2 <- estimate_richness(dada2, split = TRUE, measures = "Shannon")
#sample_data(dada2)$alpha <- as.numeric(shannon_dada2[,1])


########## Compare pair-wise and plot with ggplot


# function subset to include only shannon diversity, split by desired factor
alpha_prep <- function(sample, var){
  df <- as.data.frame(sample_data(sample) %>% select(alpha,var))
  df <- as.data.frame(as.matrix(df))
  df$alpha <- sapply(df[,1], as.double)
  df$var <- as.factor(df$var)
  return(df)
}


###################### Lifestage ##########################


# subset by alpha and lifestage
alpha_prep(esto, Lifestage)

# subset by alpha and lifestage without function
alpha_life_esto <- as.data.frame(sample_data(esto) %>% select(alpha, Lifestage))
alpha_life_esto <- as.data.frame(as.matrix(alpha_life_esto))
alpha_life_esto$alpha <- sapply(alpha_life_esto[,1], as.double)
alpha_life_esto$Lifestage <- as.factor(alpha_life_esto$Lifestage)

#alpha_life_dada2 <- as.data.frame(sample_data(dada2) %>% select(alpha, Lifestage))
#alpha_life_dada2 <- as.data.frame(as.matrix(alpha_life_dada2))
#alpha_life_dada2$alpha <- sapply(alpha_life_dada2[,1], as.double)
#alpha_life_dada2$Lifestage <- as.factor(alpha_life_dada2$Lifestage)


# perform Mann-Whitney U test on esto
wilcox.test(alpha ~ Lifestage, data = alpha_life_esto)

# perform Mann-Whitney U test on dada2
#wilcox.test(alpha ~ Lifestage, data = alpha_life_dada2)

