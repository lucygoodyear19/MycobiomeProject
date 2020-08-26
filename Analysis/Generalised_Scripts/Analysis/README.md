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

### 1) [summary.R](#1.summary.R)

Uses Kruskal-Wallis to compare abundances of the different fungal species with Bd postiive/negative status.

### 2) [phylum_diversity.R](#2.phylum_diversity.R)

Uses Kruskal-Wallis to compare abundances of the different fungal species with Bd postiive/negative status.

### 3) [beta_diversity.R](#3.-beta_dviersity.R)

Transforms the data using an isometric-log-ratio (ilr) transform then performs a Principle Components Analysis (PCA) to compare beta diversity across factors.

### 4) [abundances_vs_bd.R](#4.abundances_vs_bd.R)

Uses Kruskal-Wallis to compare abundances of the different fungal species with Bd postiive/negative status.

### 5) [cooccur.R](#5.cooccur.R)

Uses the cooccur package to to look at cooccurance between fungal ASVs present in two or more samples, focusing on Bd.

### 6) [cooccur-genus.R](#6.cooccur-genus.R)

Idenitcal to cooccur.R, expects performs co-occurrence at the the genus level.

### 7) [cooccur-genus_analysis.R](#7.cooccur-genus_analysis.R)

Analyses results from co-occur-genus script and produces plots and tables.

&nbsp;

## Details of Scripts

### 0) analysis_args.R

#### *Requirements*

#### Arguments: 

None.

#### Packages:

None.

#### *How to use*

This script is stored in the relevent results folder and must contain the following 5 variables:  

dada2_data_path <- path to input data  
results_path <- path to output folder  
samp_vars <- vector containing indepedent variable names as strings 
samp_vars_1host <- same as samp_vars but with only one host variable, e.g. Order OR Family, with a strict number of 3 values.
datasets <- vector containing name of dataset as string, if analysing multiple datasets at once, can include those (note, multiple dataset option is currently only available in abundances_vs_bd.R)

This script is needed for all the below analysis scripts as the argument containing the necessary variables.

&nbsp;

### 1) summary.R

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



#### *How to use*

This script is my only non-generalised analysis script (apart from the last section in the phylum_diversity.R script below)but I hope to generalise it at a later date.

This script cannot currently be run from the command line. The analysis_args.R script needs to be run first in order for this script to access the arguments it needs. I plan to make this into a script that can be run on the command line but I haven't had time yet.

This script is the longest of my analysis scripts and can be broken down into 4 parts.

Part 1  

Plots a map of sample locations with a coloured scale based on how many samples were collected from each country. Plotting is split by Taiwan and non-Taiwan since Taiwan contains many more data points than anywhere else.
Plots a map of average Shannon alpha diversity per country.
Plots a map of average ASV richness per country.


Plots pie chart summarising overall phylum break down.
Plots bar chart summarising taxa breakdown by presence/absence. 8 charts are produced with phylum, class and order for each of lifestage, Bd status, country and order of amphibian host.
Plots bar chart summarising phylum breakdown by abundance. 5 for each of lifestage, Bd status, country, and order and family of amphibian host.


Outputs:

1) 

&nbsp;


### 2) phylum_diversity.R

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

### 3) abundances_vs_bd.R

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

### 4) beta_diversity.R

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

### 5) cooccur.R

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

### 6) cooccur-genus.R

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

### 7) cooccur-genus_analysis.R

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