############################################################################################################
############################################### Associations with Bd #######################################
############################################################################################################


# Author: Luke Goodyear (leg19@ic.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())



##########################################################################################################
############################################### Set up ###################################################


# load packages
library(dplyr)
library(phyloseq)
library(Biostrings)
#library(nlme)
library(lme4)
library(effsize) # load effect size package
library(ggplot2)
library(gridExtra)

setwd("~/Dropbox/luke/documents/academia/mres/mycobiome_project/analysis/results/Cameroon/")

source("analysis_args.R")

# print arguments as check
print(results_path) # path to data folder

# check if results directory exists and if not, create it
ifelse(!dir.exists(file.path(paste0(results_path, "bd_associations/"))), dir.create(file.path(paste0(results_path, "bd_associations/"))), FALSE)
# define path out
path_out <- paste0(results_path, "/bd_associations/")

# load filtered phyloseq object for DADA2 pipeline
dada2 <- readRDS(paste0(dada2_data_path, "physeqob_DADA2.rds"))
mets <- as.data.frame(as.matrix(sample_data(dada2)))
mets <- as.data.frame(sapply(mets, as.factor))
mets$A_Genus_Species <- as.double(mets$A_Genus_Species)
mets$Elevation_m <- as.double(as.character(mets$Elevation_m))
mets$Latitude <- as.double(as.character(mets$Latitude))
mets$Longitude <- as.double(as.character(mets$Longitude))
#mets<-mets[!(mets$A_Order=="Gymnophiona"),]
#mets<-mets[!(mets$Lifestage=="Metamorph"),]
#mets<-mets[!(mets$Locality==""),]





# define vector of variables of interest
myvars <- c("A_Order", "A_Family", "A_Genus", "A_Genus_Species", "Latitude", "Longitude", "Elevation_m", "Lifestage")



##########################################################################################################
############################################### Functions ################################################


# function subset to include only bd, split by desired factor
bd_prep <- function(sample, var){
  # muliple conversions required due to format of sample_data(df)
  df <- as.data.frame(as.matrix(sample_data(sample))) %>% select(Bd_GE,var)
  # convert to alpha character before double or output would be factor levels
  df$Bd_GE <- sapply(df[,1], as.character) 
  df$Bd_GE <- sapply(df[,1], as.double)
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
         aes(x = factor(sample_data(sample)[[var]]), y = Bd_GE)) + 
    labs(y = "Bd GE", x = var) +
    geom_boxplot() +
    scale_x_discrete(labels = labs) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 270, hjust = 0)) 
}


#####################################################################################################
##################################### Calculate correlations ########################################


corr_total <- as(mets, "data.frame")
corr_total <- corr_total %>% mutate_if(is.character,as.factor) # convert all chr columns to factors
corr_total <- corr_total[,-c(1:4, 9, 15:20, 22, 25)] # remove irrelevant columns
corr_total_fac <- data.matrix(corr_total) # convert all factors to underlying integer levels
# run correlation function and capture output
corr_out <- capture.output(cor(corr_total_fac))
# save output to a text file
cat(corr_out, file=paste0(path_out, "Correlations_dada2.txt"), append=TRUE)

# get data in the right format
df <- as.data.frame(sapply(mets, as.factor))
df <- df[,-c(1:4, 9, 15:20, 22, 25)] # remove irrelevant columns
# check all variables are the correct class
str(df)
pdf(paste0(path_out, "pairs.pdf"))
pairs(df, col = df$Bd)
dev.off()


################################################################################################
########################################### Box plots ##########################################


# for loop to create many two column dataframes with the required variables and alpha
# create empty list to store dataframes
ls <- c()
for (i in colnames(corr_total_fac)[1:11]){
  # run alpha_prep on all variables and save according to variable name
  assign(paste0("bd_", i), bd_prep(dada2, i))
  # save all bd_prep results in a list
  name <- paste0("bd_", i)
  df <- get(paste0("bd_", i))
  ls[[name]] <- df
}

# plot box plots using base R
pdf(paste0(path_out, "box.pdf"))
for (i in 1:(ncol(corr_total_fac)-1)) {
    sub_df <- as.data.frame(ls[i])
    names(sub_df) <- c("bd", colnames(corr_total_fac)[i])
    print(box_plot(dada2, colnames(corr_total_fac)[i], sub_df))
}
dev.off()


###################################################################################################
######################################## Simple regressions #######################################


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


####################################################################################################
################################################# Stats ############################################


############################################# DADA2 ##############################################


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








########################################################################################################
############################################### Plots ##################################################


# plot stacked bar chart for Bd status per attribute
plot_bd_bar <- function(samps, attrib, cont) {
  bd_table <- as.data.frame(c())
  labs <- list() # for plot labels
  factors <- c() # for correct label order
  for (country in unique(samps[[attrib]])) {
    samps_sub <- samps %>% filter(samps[[attrib]] == country)
    samps_sub[[attrib]] <- as.character(samps_sub[[attrib]])
    labs[[as.character(country)]] <- nrow(samps_sub)
    factors <- c(factors, country)
    subs <- subset(samps_sub, select = c(attrib, "Bd"))
    count <- as.data.frame(table(subs, exclude = NULL))
    count$Freq <- count$Freq/sum(count$Freq)
    bd_table <- rbind(bd_table, count)
  }
  if (cont == F) {
    labs <- labs[order(names(labs), decreasing=F)]
  } else {
      nos <- as.double(names(labs))
      labs <- labs[order(nos)]
  }
  bd_table$Bd <- as.character(bd_table$Bd)
  bd_table$Bd[is.na(bd_table$Bd)] <- "NA"
  bd_table$Bd <- as.factor(bd_table$Bd)
  bd_table[[attrib]] <- factor(bd_table[[attrib]], levels = factors)
  names(bd_table)[1] <- "Variable"
  bd_table$Variable <- as.character(bd_table$Variable)
  bd_table <- bd_table[order(bd_table[,1], decreasing=F),]
  if (cont == T) {
    bd_table$Variable <- as.double(as.character(bd_table$Variable))
    bd_table <- bd_table[order(bd_table$Variable),]
    bd_table$Variable <- as.factor(bd_table$Variable)
  }
  # plot stacked bar chart for all countries
  lab_ls <- c()
  for (i in 1:length(names(labs))) {
    lab_ls <- c(lab_ls, paste0(names(labs)[i], " (n = ", labs[[names(labs)[i]]], ")"))
  }
  pdf(paste0(results_path, "bd_breakdown_", attrib, ".pdf"))
  ggplot(bd_table, aes(fill=Bd, y=Freq, x=Variable)) + 
    geom_bar(position="stack", stat="identity") +
    coord_fixed(ratio=5) +
    scale_fill_manual(values = c("#90C987", "#6195CF", "#F6C141")) + # colour blind friendly colours
    scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0,0))) +
    labs(fill = "Bd Status") +
    ylab("Proportion") +
    scale_x_discrete(labels = lab_ls) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.y = element_line(colour = "black"))
  dev.off()
}

plot_bd_bar(mets, "Mountain", F)
plot_bd_bar(mets, "Locality", F)
plot_bd_bar(mets, "Latitude", T)
plot_bd_bar(mets, "Elevation_m", T)
plot_bd_bar(mets, "A_Family", F)
plot_bd_bar(mets, "A_Genus_Species", F)


# plot chart for Bd GE per attribute
plot(mets$Latitude, mets$Bd_GE)
plot(mets$Elevation_m, mets$Bd_GE)
plot(mets$Mountain, mets$Bd_GE)
plot(mets$Locality, mets$Bd_GE)
plot(mets$A_Family, mets$Bd_GE)
plot(mets$Bd_GE, mets$Genus_Species)

# Bd positive samples only
mets_bdpos <- mets[!mets$Bd_GE==0.00,]
plot(mets_bdpos$Latitude, mets_bdpos$Bd_GE)
plot(mets_bdpos$Elevation_m, mets_bdpos$Bd_GE)
plot(mets_bdpos$Mountain, mets_bdpos$Bd_GE)
plot(mets_bdpos$Locality, mets_bdpos$Bd_GE)
plot(mets_bdpos$A_Family, mets_bdpos$Bd_GE)
plot(mets_bdpos$Bd_GE, mets_bdpos$Genus_Species)



######################################################################################################################
############################################### Regressions ##########################################################


# regression set up
mets <- as.data.frame(sapply(mets, as.factor))
mets$A_Genus_Species <- as.double(mets$A_Genus_Species)
mets$A_Family <- as.double(mets$A_Family)
mets$Elevation_m <- as.double(mets$Elevation_m)
mets$Latitude <- as.double(mets$Latitude)
mets$Longitude <- as.double(as.character(mets$Longitude))
mets$Bd_GE<- as.double(as.character(mets$Bd_GE))


# Logistics Regression
glm_fit <- glm(Bd ~ A_Genus_Species + Elevation_m + Latitude, data = mets, family = binomial)
summary(glm_fit)

res <- residuals(glm_fit, type="deviance")
plot(log(predict(glm_fit)), res)
abline(h=0, lty=2)
qqnorm(res)
qqplot(res)
qqline(res)



###################### NOT SIGNIFICANT
# Linear Regression on all samples
glm_fit_bd_ge <- glm(Bd_GE ~ A_Genus_Species + Elevation_m + Latitude, data = mets, family = gaussian)
summary(glm_fit_bd_ge)
res_bd_ge <- residuals(glm_fit_bd_ge, type="deviance")
plot(log(predict(glm_fit_bd_ge)), res_bd_ge)
abline(h=0, lty=2)
qqnorm(res_bd_ge)
qqplot(res_bd_ge)
qqline(res_bd_ge)

# Linear regression on positive samples only
glm_fit_bdpos <- glm(Bd_GE ~ A_Genus_Species  + Elevation_m + Latitude, data = mets_bdpos, family = gaussian)
summary(glm_fit_bdpos)
res <- residuals(glm_fit_bdpos, type="deviance")
plot(log(predict(glm_fit_bdpos)), res)
abline(h=0, lty=2)
qqnorm(res)
qqplot(res)
qqline(res)






