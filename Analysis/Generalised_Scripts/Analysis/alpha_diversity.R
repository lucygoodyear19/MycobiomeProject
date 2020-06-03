##############################################################################
###################### DADA2 and Esto Alpha Diversity ########################
##############################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())



##############################################################################
################################## Set up ####################################


# load packages
library(dplyr)
library(phyloseq)
library(Biostrings)
#library(nlme)
library(lme4)
library(effsize) # load effect size package
library(ggplot2)
library(gridExtra)


# import arguments to run script on specific country data
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) # setup to accept arguments from command line
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied in the form of an R-script containing the following arguments: 
       1) path to results folder, results_path
       2) path to dada2 data folder, dada2_path 
       3) path to esto data folder, esto_path")
}

# load arguments into script
source(args)

######## TEMP
source("analysis_args.R")

# print arguments as check
print(results_path) # path to data folder
print(dada2_path) # path to folder that will contain filtered data files
print(esto_path) # path to folder that will contain filtered data files

# define path out
path_out <- paste0(results_path, "Alpha_Diversity/")

# load filtered phyloseq object for both DADA2 pipeline and Estonianpipeline
dada2 <- readRDS(paste0(dada2_path, "physeqob_DADA2_complete.rds"))
esto <- readRDS(paste0(esto_path, "physeqob_esto.rds"))

# define vector of variables of interest
myvars <- c("Site_code", "Species", "Age", "Lifestage", "Habitat_type", "Latitude", "Longitude", "Altitude_m", "Bd", "alpha")
# define vector of independent variables
ind_vars <- myvars[myvars != "alpha"]
# define vector of random effects variables
rand_effs <- c("Notes", "Sample_type", "Date", "ID")
# vecotr of datasets
datasets <- c("dada2", "esto")


# temp until new data filtering script is used
asv_seqs <- Biostrings::DNAStringSet(taxa_names(dada2))
names(asv_seqs) <- taxa_names(dada2)
dada2 <- merge_phyloseq(dada2, asv_seqs)
taxa_names(dada2) <- paste0("ASV", seq(ntaxa(dada2)))
dada2


#################################################################################
########################### Calculate alpha diversity ###########################


# estimate value for alpha for each sample

shannon_esto <- estimate_richness(esto, split = TRUE, measures = "Shannon")
sample_data(esto)$alpha <- as.numeric(shannon_esto[,1])

# dada2
shannon_dada2 <- estimate_richness(dada2, split = TRUE, measures = "Shannon")
sample_data(dada2)$alpha <- as.numeric(shannon_dada2[,1])


################################################################################
################################ Functions #####################################


# function subset to include only shannon diversity, split by desired factor
alpha_prep <- function(sample, var){
  # muliple conversions required due to format of sample_data(df)
  df <- as.data.frame(as.matrix(sample_data(sample) %>% select(alpha,var)))
  # convert to alpha character before double or output would be factor levels
  df$alpha <- sapply(df[,1], as.character) 
  df$alpha <- sapply(df[,1], as.double)
  # convert var to factor
  df[[var]] <- as.factor(df[[var]])
  return(df)
}


# function to plot labeled boxplot
box_plot <- function(sample, var, subset) {
  labs <- c()
  for (i in 1:length(unique(levels(subset[[var]])))) {
    labs <- c(labs, paste0(levels(subset[[var]])[i], " (n = ", table(sample_data(sample)[[var]])[i], ")"))
  }
  ggplot(sample_data(sample), 
         aes(x = factor(sample_data(sample)[[var]]), y = alpha)) + 
    labs(y = "Shannon Alpha Diversity", x = var) +
    geom_boxplot() +
    scale_x_discrete(labels = labs) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 270, hjust = 0)) 
}


################################################################################
########################### Calculate correlations #############################


# get data in the correct format
for (j in datasets) {
  corr_total <- as.data.frame(sample_data(get(j))) # convert to data frame
  corr_total <- corr_total %>% mutate_if(is.character,as.factor) # convert all chr columns to factors
  corr_total <- corr_total[,-c(1,2,5,8,10,11,12,13,16)] # remove irrelevant columns
  assign(paste0("corr_total_fac_", j), data.matrix(corr_total)) # convert all factors to underlying integer levels
  # run correlation function and capture output
  (assign(paste0("corr_out_", j), capture.output(cor(get(paste0("corr_total_fac_", j))))))
  # save output to a text file
  cat(paste0("Correlations for ", j), get(paste0("corr_out_", j)), file=paste0(path_out, "Correlations_", j, ".txt"), append=TRUE)
}

# print scatterplots for all pairs of varibles
for (j in datasets) { # loop over each dataset
  # get data in the right format
  df <- as.data.frame(sapply(sample_data(get(j)), as.factor))
  df$alpha <- as.numeric(df$alpha)
  df <- df[,-c(1,2,5,8,10,11,12,13,16)] # remove irrelevant columns
  # check all variables are the correct class
  str(df)
  pdf(paste0(path_out, j, "_pairs.pdf"))
  pairs(df)
  dev.off()
}

################################################################################
################################## Box plots ###################################


# for loop to create many two column dataframes with the required variables and alpha
for (j in datasets) {
  # create empty list to store dataframes
  lst <- c()
  for (i in colnames(get(paste0("corr_total_fac_", j))[,1:9])){
    # run alpha_prep on all varibales and save according to variable name
    assign(paste0("alpha_", i, "_", j), alpha_prep(get(j), i))
    # save all alpha_prep results in a list
    name <- paste0("alpha_", i, "_", j)
    df <- get(paste0("alpha_", i, "_", j))
    lst[[name]] <- df
  }
  assign(paste0("ls_alpha_", j), lst)
}


# plot box plots using base R and scatter plot using phyloseq function
for (j in datasets) {
  pdf(paste0(path_out, j, "_box.pdf"))
  for (i in 1:(ncol(get(paste0("corr_total_fac_", j)))-1)) {
    sub_df <- as.data.frame(get(paste0("ls_alpha_", j))[i])
    names(sub_df) <- c("alpha", colnames(get(paste0("corr_total_fac_", j)))[i])
    grid.arrange(box_plot(get(j), colnames(get(paste0("corr_total_fac_", j)))[i], sub_df),
             plot_richness(get(j), 
                           measures = "shannon", 
                           x = colnames(get(paste0("corr_total_fac_", j)))[i]),
             ncol = 2)
  }
  dev.off()
}


################################################################################
################################ Simple regressions ############################


# linear regression models for each variable singly
for (j in datasets) { # loop over each dataset
  # get data in the right format
  df <- as.data.frame(sapply(sample_data(get(j)), as.factor))
  df$alpha <- as.numeric(df$alpha)
  df <- df[,-c(1,2,5,8,10,11,12,13,16)] # remove irrelevant columns
  # check all variables are the correct class
  str(df)
  for (i in 1:(ncol(df)-1)) { # loop over each independent variable
    # calulate linear regression
    lmod <- lm(alpha ~ df[,i],
                  data = df)
    # capture summary output
    summ <- capture.output(summary(lmod))
    # save output to a text file
    cat(paste0("\n\nSimple regressions for ", j, " - ", colnames(df)[i]), summ, file=paste0(path_out, "Simple_regressions_for_", j, ".txt"), sep = "\n", append=TRUE)
    }
}


################################################################################
##################################### Stats ####################################


#################################### DADA2 #####################################


# loop to run stats tests on all variables one by one
for (i in 1:(ncol(corr_total_fac)-1)) {
  sub_df <- as.data.frame(ls_alpha[i])
  names(sub_df) <- c("alpha", colnames(corr_total_fac)[i])
  # perform Mann-Whitney U test and capture output
  (wilcox_out <- pairwise.wilcox.test(sub_df$alpha, sub_df[[colnames(corr_total_fac)[i]]]))
  (cohen_out <- cohen.d(sub_df$alpha ~ sub_df[[colnames(corr_total_fac)[i]]]))
  # save outputs to a text files
  cat(paste0("Wilcox for DADA2 - ", colnames(corr_total_fac)[i]), wilcox_out, file="Wilcox_DADA2.txt", append=TRUE)
  cat(paste0("Cohen for DADA2 - ", colnames(corr_total_fac)[i]), cohen_out, file="Cohen_DADA2.txt", append=TRUE)
}


#################################################################################
############################## using phyloseq ###################################


################################## DADA2 ########################################


plot_richness(dada2, measures = c("shannon", "simpson"), x = "Altitude_m", color = "Bd")
plot_richness(dada2, measures = c("shannon", "simpson"), x = "Latitude", color = "Bd")
plot_richness(dada2, measures = c("shannon", "simpson"), x = "Longitude", color = "Bd")
plot_richness(dada2, measures = c("shannon", "simpson"), x = "Lifestage", color = "Bd")
plot_richness(dada2, measures = c("shannon", "simpson"), x = "Habitat_type", color = "Bd")
plot_richness(dada2, measures = c("shannon", "simpson"), x = "Species", color = "Bd")
plot_richness(dada2, measures = c("shannon", "simpson"), x = "Site_code", color = "Bd")
plot_richness(dada2, measures = c("shannon", "simpson"), x = "Age", color = "Bd")
plot_richness(dada2, measures = c("shannon", "simpson"), x = "Notes", color = "Bd")
plot_richness(dada2, measures = c("shannon", "simpson"), x = "Sample_type", color = "Bd")

# random variance
plot_richness(dada2, measures = c("shannon", "simpson"), x = "Date", color = "Bd")
plot_richness(dada2, measures = c("shannon", "simpson"), x = "ID", color = "Bd")

# only one level
#plot_richness(dada2, measures = "shannon", x = "Sex", color = "Bd")
#plot_richness(dada2, measures = "shannon", x = "qPCR", color = "Bd")
#plot_richness(dada2, measures = "shannon", x = "wild_captive", color = "Bd")


#################################### Esto ######################################


plot_richness(esto, measures = c("shannon", "simpson"), x = "Bd")

plot_richness(esto, measures = c("shannon", "simpson"), x = "Altitude_m", color = "Bd")
plot_richness(esto, measures = c("shannon", "simpson"), x = "Latitude", color = "Bd")
plot_richness(esto, measures = c("shannon", "simpson"), x = "Longitude", color = "Bd")
plot_richness(esto, measures = c("shannon", "simpson"), x = "Lifestage", color = "Bd")
plot_richness(esto, measures = c("shannon", "simpson"), x = "Habitat_type", color = "Bd")
plot_richness(esto, measures = c("shannon", "simpson"), x = "Species", color = "Bd")

# redundant
plot_richness(esto, measures = c("shannon", "simpson"), x = "Age", color = "Bd") # lifestage
plot_richness(esto, measures = c("shannon", "simpson"), x = "Site_code", color = "Bd") # habitat type

# non-random variance to ignore (of no interest)
plot_richness(esto, measures = c("shannon", "simpson"), x = "Notes", color = "Bd")
plot_richness(esto, measures = c("shannon", "simpson"), x = "Sample_type", color = "Bd")

# random variance
plot_richness(esto, measures = c("shannon", "simpson"), x = "Date", color = "Bd")
plot_richness(esto, measures = c("shannon", "simpson"), x = "ID", color = "Bd")

# only one level
#plot_richness(esto, measures = "shannon", x = "Sex", color = "Bd")
#plot_richness(esto, measures = "shannon", x = "qPCR", color = "Bd")
#plot_richness(esto, measures = "shannon", x = "wild_captive", color = "Bd")


#################################################################################
############################ Mixed effect models ################################


dada2_sample_df <- as(sample_data(dada2),"data.frame")
dada2_sample_df <- dada2_sample_df[,colnames(dada2_sample_df) %in% myvars | colnames(dada2_sample_df) %in% rand_effs]

for (i in 1:length(ind_vars)){
  dada2_sample_df[[ind_vars[i]]] <- as.factor(dada2_sample_df[[ind_vars[i]]])
}

for (i in 1:length(rand_effs)){
  dada2_sample_df[[rand_effs[i]]] <- as.factor(dada2_sample_df[[rand_effs[i]]])
  dada2_sample_df[[rand_effs[i]]]<- as.numeric(dada2_sample_df[[rand_effs[i]]])
}
# check all variables are the correct class
str(dada2_sample_df)

# scale alpha to between 1 and 0
dada2_sample_df$alpha_scale <- scale(dada2_sample_df$alpha)[,1]

#model_nlme <- lme(alpha ~ Lifestage + Habitat_type + Lifestage*Habitat_type,
#         random = ~1|Species,
#         data = as(sample_data(dada2),"data.frame"))

model_lme4_1 <- lmer(alpha ~ Lifestage + Habitat_type + Lifestage*Habitat_type + (1|Date),
                     data = dada2_sample_df)
summary(model_lme4_1)

(model_lme4_2 <- glmer(alpha ~ 1 + Lifestage*Habitat_type*Latitude*Longitude*Altitude_m*Bd*Species + (1|Date) + (1|ID) + (1|Sample_type) + (1|Notes),
                       data = dada2_sample_df))
lmtry <- lm(alpha ~ Lifestage*Habitat_type*Latitude*Longitude*Altitude_m*Bd*Species,
            data = dada2_sample_df)

model <- glm(alpha_scale ~ 1 + 
                     Lifestage*Habitat_type*Latitude*Longitude*Altitude_m*Bd*Species + 
                     (1|dada2_sample_df[[rand_effs[1]]]) + 
                     (1|dada2_sample_df[[rand_effs[2]]]) + 
                     (1|dada2_sample_df[[rand_effs[3]]]) + 
                     (1|dada2_sample_df[[rand_effs[4]]]),
             family = binomial,
             data = dada2_sample_df)
require("olsrr")
ols_step_all_possible(model)


## end of script
