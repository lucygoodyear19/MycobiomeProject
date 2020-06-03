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

### 2) run_DADA2_HPC_pipeline.sh

#### *How to use*

This script has been written to run on a PBS HPC so adjustments may need to be made to accommodate other machines.

The following command should be used to run the script, making sure to check your own releative paths:

```qsub -v "arg_path=run/country" ../../../../Generalised_Scripts/run_DADA2_HPC_pipeline.sh ```

Where ```run/country``` should be replaced with the relevent path. For example:

```qsub -v "arg_path=Taiwan_Vietnam_2016/Taiwan" ../../../../Generalised_Scripts/run_DADA2_HPC_pipeline.sh ```

&nbsp;

### 3) DADA2_HPC_pipeline.R

#### *Requirements*

#### Arguments: 

DADA2_args.R 

#### Packages:

R:

dada2 (v.1.14.1)   
ShortRead (v.1.44.3)  
Biostrings (v.2.54.0)  

#### *How to use*

This is a script following the DADA2 ITS2 tutorial with some modifications to use only the forward reads and to run on HPC. Thanks to Phil Jervis for giving me his script to edit and to the creators of the dada2 pipeline (https://benjjneb.github.io/dada2/ITS_workflow.html).

No modifications should be required in order for this script to run with any R-script argument containing the correct variables and with the sequences in the correct format.

&nbsp;

### 4) Data_filtering.R

#### *Requirements*

#### Arguments: 

DADA2_args.R

#### Packages:

R:
   
phyloseq (v.1.30.0)  

#### *How to use*

Filters and creates phyloseq object.

&nbsp;

### 5) phylogenetic_trees.R

#### *Requirements*

#### Arguments: 

DADA2_args.R

#### Packages:

R:
   
DECIPHER (v.2.14.0)  
phangorn (v.2.5.5)    
phyloseq (v.1.30.0)

#### *How to use*

Produces a phylogenetic tree, which is saved to the same phyloseq object as generated above.

&nbsp;

### 6) DADA2_args.R

#### *Requirements*

#### Arguments: 

None.

#### Packages:

None.

#### *How to use*

This script is stored in the relevent country folder and contains 12 variables:  

REV - reverse primer  
FWD - forward primer  
base_prefix - base name of sequences  
run  - run directory
country - country directory  
root_path - your path to the directory containing all runs  
path - path to directory containing sequences, using the root_path  
path2 - path to filtN directory, using the root_path  
path_out - path to directory where output will be sent to, using root_path  
metadata_path - path to metadata file, using root_path  
cutadapt - path to cutadapt directory, using root_path  
unite.ref - path to UNITE database, using root_path 