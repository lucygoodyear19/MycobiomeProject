# Mycobiome Project Generalised Analysis Scripts Repository

&nbsp;

*Author: Lucy Goodyear*  
*Created: 03/06/20*

This repository contains all scripts required for analysing mycobiome data. Each script has been generalised so can easily be used for different datasets.

&nbsp;

## Requirements

All code has been created for Mac so there may be a few differences in commands with respect to Linux.

&nbsp;

## List of Scripts

### 0) [analysis_args.R](#0.analysis_args.R)

An R-script containing the arguments required for all analysis scripts below so that the scripts can be run from the command line.

### 1) [abundances_vs_bd.R](#1.abundances_vs_bd.R)

Uses Kruskal-Wallis to compare abundances of the different fungal species with Bd postiive/negative status.

### 2) [alpha_diversity.R](#2.alpha_diversity.R)

Uses several methods to compare alpha diversity cross different factors.

### 3) [beta_diversity.R](#3.-beta_dviersity.R)

Uses PCA and NMDS to compare beta diversity across the different factors.

### 4) [cooccur.R](#4.cooccur.R)

Uses the cooccur package to to look at cooccurance between fungal species present in two or more smamples, focusing on Bd.

### 5) [tree_plotting.R](#4.tree_plotting.R)

Plots a phylogenetic tree from both subsetting the complete by order and using a non-taxonomic assignment method to locate possible unknown chytrids.

&nbsp;

## Details of Scripts

### 0) analysis_args.R

#### *Requirements*

#### Arguments: 

None.

#### Packages:

None.

#### *How to use*

This script is stored in the relevent country folder and contains 4 variables:  

root - full path to directory containing scripts, data and results
run - run directory name
country - name of country
year - year in name of run

&nbsp;

### 1) abundances_vs_bd.R

#### *Requirements*

#### Arguments: 

analysis_args.R

#### Packages:

R:
   
dplyr (v.)  
phyloseq (v.)    
microbiomeSeq (v.)
effsize
ggplot2

#### *How to use*

This script can be run locally and results in 3 different output sets, corresponding to 3 different methods.

Two different methods have been used to perform a Kruskal-Wallis test on taxa abundance data. Both methods use the base R function ```kruskal.test()``` but the microbiomeSeq pacakge uses different methods to adjust the p-values (to account for multiple testing) and to select the significant taxa.

The ```kruskal_abundance()``` function from the microbiomeSeq package is used to generate adjusted p-values and plotting data for the significant taxa. This output is saved to a .txt. The data is then plotted and saved to a pdf.

The base R method uses two different methods to adjust the p-values: Benjamini_Hochberg and Holm.

Outputs:

text file microbiomeSeq
text file base R
plot micorbiomeSeq
plot base R Benjamini-Hochberg
plot base R Holm

Each of these will be repated for each dataset.

An example command to run the script:

```Rscript /Users/lucy/Documents/MRes/MycobiomeProject/Analysis/Generalised_Scripts/Analysis/abundances_vs_bd.R analysis_args.R ```

&nbsp;
