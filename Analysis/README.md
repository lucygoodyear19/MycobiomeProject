# Mycobiome Project Generalised Analysis Scripts Repository

&nbsp;

*Author: Lucy Goodyear*  
*Created: 03/06/20*

This repository contains all scripts required for analysing mycobiome data. Scripts 1-5 have been generalised so can easily be used for different datasets.

&nbsp;

## Requirements

All code has been created for Mac so there may be a few differences in commands with respect to Linux.

All scripts are written in R and use version 3.6.3.

The following R packages are required to run all analysis scripts (there is also a breakdown per script):

dplyr (v.1.0.1)  
tidyr (v.1.1.1)   
ggplot2 (v.3.3.2)    
viridis (v.0.5.1)     
sf (v.0.9.5)  
rnaturalearth (v.0.1.0)  
rnaturalearthdata (v.0.1.0)  
effsize (v.0.8.0)   
VennDiagram (v.1.6.20)  
compositions (v.2.0.0)  
phyloseq (v.1.30.0)    
microbiomeSeq (v.0.1)   
microbiome (v.1.8.0)    
vegan (v.2.5.6)   
microbiomeSeq (v.0.1)  
cooccur (v.1.3)

&nbsp;

## List of Scripts

### 0) [analysis_args.R](#0.-analysis_args)

An R-script containing the arguments required for all analysis scripts below so that the scripts can be run from the command line. This is not in the repository.

### 1) [abundances_vs_bd.R](#1-abundances_vs_bd)

Uses Kruskal-Wallis to compare abundances of the different fungal species with Bd postiive/negative status.

### 2) [beta_diversity.R](#2.-beta_dviersity)

Transforms the data using an isometric-log-ratio (ilr) transform then performs a Principle Components Analysis (PCA) to compare beta diversity across factors.

### 3) [cooccur.R](#3.-cooccur)

Uses the cooccur package to to look at cooccurance between fungal ASVs present in two or more samples, focusing on Bd.

### 4) [cooccur-genus.R](#4.-cooccur-genus)

Idenitcal to cooccur.R, expects performs co-occurrence at the the genus level.

### 5) [cooccur-genus_analysis.R](#5.-cooccur-genus_analysis)

Analyses results from co-occur-genus script and produces plots and tables.

### 6) [phylum_diversity.R](#6.-phylum_diversity)

Merges similar phyla, transforms the data using an isometric-log-ratio (ilr) transform then performs a Principle Components Analysis (PCA) to compare phylum diversity across factors.

### 7) [summary.R](#7.-summary)

Summarises the data and produces a number of plots and tables.



&nbsp;



## Details of Scripts

### 0) analysis_args

#### *Requirements*

#### Arguments: 

None.

#### Packages:

None.

#### *How to use*

Note, this script is NOT in the repository. This script is stored in the relevent results folder and must contain the following 5 variables:  

dada2_data_path <- path to input data  
results_path <- path to output folder  
samp_vars <- vector containing indepedent variable names as strings 
samp_vars_1host <- same as samp_vars but with only one host variable, e.g. Order OR Family, with a strict number of 3 values.
datasets <- vector containing name of dataset as string, if analysing multiple datasets at once, can include those (note, multiple dataset option is currently only available in abundances_vs_bd.R)

This script is needed for all the below analysis scripts as the argument containing the necessary variables.



&nbsp;



### 1. abundances_vs_bd

#### *Requirements*

#### Arguments: 

analysis_args.R

#### Packages:

R:

dplyr (v.1.0.1)  
phyloseq (v.1.30.0)  
microbiomeSeq (v.0.1)  
effsize (v.0.8.0)  
ggplot2 (v.3.3.2)  

#### *How to use*

A directory called "Abundances_vs_Bd" is created if it does not already exist; all results are sent here. 

This script can be run locally and results in 2 different output sets, corresponding to 2 different methods.

Two different methods have been used to perform a Kruskal-Wallis test on taxa abundance data. Both methods use the base R function ```kruskal.test()``` but the microbiomeSeq pacakge uses different methods to adjust the p-values (to account for multiple testing) and to select the significant taxa. The ```kruskal_abundance()``` function from the microbiomeSeq package is used to generate adjusted p-values and plotting data for the significant taxa. This output is saved to a .txt. The data is then plotted and saved to a pdf. The base R method first transforms the data using a centre-log-ratio transform, performs Kruskal-Wallis tests and uses the Benjamini-Hochberg method to adjust the p-values.

Outputs:

1) mbseq_Abundance_associated_with_Bd_for_dada2.txt:  
    Contains all data for significant ASVs, as described by the microbiomeSeq package.
2) mbseq_dada2_abun_bd_plot.pdf:  
    Contains plots of all significant ASVs, including p-values.
3) base_Abundance_associated_with_Bd.txt:  
    Contains Kruskal-Wallis results for significant ASVs.
4) base_abun_bd_plot.pdf:  
    Contains plots of all significant ASVs, including p-values.
5) base_pie_Phylum.pdf:  
    Pie chart for breakdown of significant ASVs with respect to phylum
6) base_pie_Class.pdf:  
    Pie chart for breakdown of significant ASVs with respect to class
7) base_pie_Order.pdf:  
    Pie chart for breakdown of significant ASVs with respect to order

Each of these will be repated for each dataset specified in analysis_args.R.

An example command to run the script:

```Rscript /Users/lucy/Documents/MRes/MycobiomeProject/Analysis/Generalised_Scripts/Analysis/abundances_vs_bd.R analysis_args.R ```



&nbsp;



### 2) beta_diversity

#### *Requirements*

#### Arguments: 

analysis_args.R

#### Packages:

R:
   
dplyr (v.1.0.1)  
tidyr (v.1.1.1)   
phyloseq (v.1.30.0)  
microbiome (v.1.8.0)  
vegan (v.2.5.6)  
compositions (v.2.0.0)  
ggplot2 (v.3.3.2)

#### *How to use*

A directory called "beta_diversity" is created if it does not already exist; all results are sent here. 

This script performs an isometric-log-ratio-transform on the abundance data in the loaded phyloseq object. A PCA is then performed on the transformed data and plots are created for every variable provided and saved as pdf files.

A PERMANOVA is performed on the 3 chosen variables (specified in analysis_args.R) and the results are printed to a text file. A distance matrix is then created and dispersion analysis is performed on the 3 variables; the resulting ANOVAs for each variable are saved to the same text file as the PERMANOVA.

At this time, PERMANOVA and dispersion analysis can oly be perfomed on **3 variables**.

Outputs:

1) PERMANOVA.txt:  
    Contains PERMANOVA and ANOVA details.
2) ilr_pca_VARIABLE.pdf:  
    Contains PCA plot for each varibale and named accordingly.

An example command to run the script:

```Rscript /Users/lucy/Documents/MRes/MycobiomeProject/Analysis/Generalised_Scripts/Analysis/beta_diversity.R analysis_args.R ```



&nbsp;



### 3) cooccur

#### *Requirements*

#### Arguments: 

analysis_args.R

#### Packages:

R:
   
phyloseq (v.1.30.0)  
cooccur (v.1.3)

#### *How to use*

Wrangles abundance data into presence absence format and adds Bd status according qPCR results in metadata.

Runs cooccur function from cooccur package and saves output to find which ASVs significantly co-occur with Bd.

Outputs:

1) cooccur_full_dada2.rds:  
    Output from cooccur function saved as Rdata to be loaded into analysis script

An example command to run the script:

```Rscript /Users/lucy/Documents/MRes/MycobiomeProject/Analysis/Generalised_Scripts/Analysis/cooccur.R analysis_args.R ```



&nbsp;



### 4. cooccur-genus

#### *Requirements*

#### Arguments: 

analysis_args.R

#### Packages:

R:
   
phyloseq (v.1.30.0)  
cooccur (v.1.3)

#### *How to use*

Wrangles abundance data into presence absence format, grouped by genus, and adds Bd status according qPCR results in metadata.

Runs cooccur function from cooccur package and saves output to find which genera significantly co-occur with Bd.

Outputs:

1) cooccur_genus_dada2.rds:  
    Output from cooccur function saved as Rdata to be loaded into analysis script

An example command to run the script:

```Rscript /Users/lucy/Documents/MRes/MycobiomeProject/Analysis/Generalised_Scripts/Analysis/cooccur-genus.R analysis_args.R ```



&nbsp;



### 5. cooccur-genus_analysis

#### *Requirements*

#### Arguments: 

analysis_args.R

#### Packages:

R:
   
phyloseq (v.1.30.0)  
cooccur (v.1.3)
dplyr (v.1.0.1)  
tidyr (v.1.1.1)   
microbiome (v.1.8.0)  
ggplot2 (v.3.3.2)

#### *How to use*

A directory called "cooccur" is created if it does not already exist; all results are sent here. 

Extracts co-occurence with regard to Bd and saves subsetted cooccur object to a csv. Then performs data wrangling in order to plot a bar chart showing propotion of genera positively/negatively co-occurring with Bd by country.

Plots positively and negatively co-occurring genera broken down by phyla and saves to pdf file. Also plots pie charts showing positively and negatively co-occurring genera broken down by phylum, class and order. Finally saves full taxonomic details of significantly co-occurring genera in two csv files.

Outputs:

1) cooccur_bd_table.csv:  
    Coocur object with respect to Bd only.  
2) cooccur_by_country.pdf:   
    Bar chart showing propotion of ASVs positively/negatively co-occurring with Bd by country
3) bar_phylum_Bd.pdf:  
    Bar chart of positively and negatively co-occurring genera broken down by phyla.  
4) pie_TAX_POS/NEGwithBd.pdf:  
    6 pie charts corresponding to break down Phylum, Class and Order for both negatively co-occurring and positively co-occurring genera.   
5) negatively_cooccuring_genus.csv:  
    Contains taxonomoic breakd down for negatively co-occurring genera.
6) positively_cooccuring_genus.csv:  
    Contains taxonomoic breakd down for negatively co-occurring genera.


An example command to run the script:

```Rscript /Users/lucy/Documents/MRes/MycobiomeProject/Analysis/Generalised_Scripts/Analysis/cooccur-genus_analysis.R analysis_args.R ```



&nbsp;



### 6. phylum_diversity

#### *Requirements*

#### Arguments: 

analysis_args.R

#### Packages:

R:

dplyr (v.1.0.1)  
phyloseq (v.1.30.0)  
ggplot2 (v.3.3.2)  
tidyr (v.1.1.1)   
microbiome (v.1.8.0)  
vegan (v.2.5.6)  
compositions (v.2.0.0)

#### *How to use*

This script cannot currently be run from the command line. The analysis_args.R script needs to be run first in order for this script to access the arguments it needs. I plan to make this into a script that can be run on the command line but I haven't had time yet.

It creates a compositional matrix, grouped by phyla, and then performs an isometric-log-ratio transform on the data. A PERMANOVA is performed on the 3 chosen variables (specified in analysis_args.R) and the results are printed to a text file. A distance matrix is then created and dispersion analysis is performed on the 3 variables; the resulting ANOVAs for each variable are saved to the same text file as the PERMANOVA. 

Finally, some phyla percentages and breakdowns are printed to screen. these breakdowns have not been generalised and are highly project specific.

Outputs:

1) permanovas_phyla.txt
    Contains PERMANOVA and ANOVA details.



&nbsp;



### 7. summary

#### *Requirements*

#### Arguments: 

analysis_args.R

#### Packages:

R:

dplyr (v.1.0.1)  
phyloseq (v.1.30.0)  
microbiomeSeq (v.0.1) 
ggplot2 (v.3.3.2)  
tidyr (v.1.1.1)   
microbiome (v.1.8.0)  
vegan (v.2.5.6)  
sf (v.0.9.5)
rnaturalearth (v.0.1.0)
rnaturalearthdata (v.0.1.0)
viridis (v.0.5.1)
VennDiagram (v.1.6.20)

#### *How to use*

This script is my only non-generalised analysis script (apart from the last section in the phylum_diversity.R script below) but I hope to generalise it at a later date.

This script cannot currently be run from the command line. The analysis_args.R script needs to be run first in order for this script to access the arguments it needs. I plan to make this into a script that can be run on the command line but I haven't had time yet.

This script is the longest of my analysis scripts and is mostly output based. In addition to the outputs below, it also calculates ASV richness and shared ASVs, which are printed to screen. This based on random sampling approach where 100 random samples are drawn from the larger of the two groups so that the groups are of equal numbers. The script finally prints to screen some tables summarising the sample data.

Outputs:

1) map_samples.pff  
    A global map of sample locations with a coloured scale based on how many samples were collected from each country.  Plotting is split by Taiwan and non-Taiwan since Taiwan contains many more data points than anywhere else.
2) map_alpha.pdf  
    A global map of average Shannon alpha diversity per country.
3) map_richness.pdf  
    A global map of average ASV richness per country.
4) phyla_pie.pdf  
    Pie chart summarising overall phylum break down.
5) CONDITION_TAXA_breakdown.pdf  
    8 bar charts summarising taxa breakdown by presence/absence with respect to TAXA (phylum and class) for each of CONDITION (lifestage, Bd status, country and order of amphibian host).
6) CONDITION_TAXA_breakdown.csv  
    5 tables containing breakdowns per CONDITION per TAXA
7) CONDITION_TAXA_abun_breakdown.pdf   
    5 bar charts summarising taxa breakdown by abundance with respect to phylum (although other taxonomic ranks can easily be added) for each of CONDITION (lifestage, Bd status, country and order and family of amphibian host).
8)  total_phylum_breakdown.csv  
    Table containing total phylum breakdown.
9) bd_breakdown.pdf  
    Bar chart showing Bd status breakdown per country.
10) TAXA_pie_bd_CONDITION.pdf  
    2 pie charts showing taxa (phyla only) breakdown by CONDITION: (a) Bd positive samples and (b) Bd negative samples.
11) pie_CONDITION.pdf  
    4 pie charts showing breakdown of samples of Bd, amphibina host Order and Family, and life stage.
12) CONDITION_asv/species_global_venn_diagram.png   
    Generates 6 Venn diagrams for each of ASV and species breakdown per CONDITION (Bd, life stage and amphibian host order).
13) Common_taxa_summary_data.txt  
    Contains details on the ASVs and species shared by all countries, such as taxonomic information and percentage of reads.
14) Most_abundant_ASV_summary_data.txt  
    Contains details on the most abundant ASV (in terms of (a) presence/absence, and (b) read numbers), such as taxonomic information and percentage of reads.