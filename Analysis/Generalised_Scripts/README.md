# Mycobiome Project Generalised Scripts Repository

&nbsp;

*Author: Lucy Goodyear*  
*Created: 04/02/20*

This repository contains all scripts required for prcoessing mycobiome data. Each script has been generalised so can easily be used for different datasets.

&nbsp;

## Requirements

All code has been created for Mac so there may be a few differences in commands with respect to Linux.

&nbsp;

## Workflow

### 1) [sample_seqs_prep.sh](#1.sample_seqs_prep.sh)

Given a tarball containing data from an Illumina MiSeq run, this script processes the sequences into a form ready to be run through the DADA2 pipeline.

### 2) [DADA2_HPC_pipeline.R](#2.-DADA2_HPC_pipeline.R)

Processes seqeunce data to produce useful outputs for analysis, including quality profiles and a phyloseq object.

### 3) [run_DADA2_HPC_pipleine.sh](#3.run_DADA2_HPC_pipleine.sh)

Runs the DADA2 R-script above on the HPC.

&nbsp;

## Scripts

### 1) sample_seqs_prep.sh

#### *Requirements*

#### Arguments: 

Mandatory arguments:

1) root name of tarball containing sequences
2) name of country 
3) number of sequences from run linked to that country  

Up to 4 countries in total can be listed as arguments but must be accompanied by the number of sequences for that country.

#### *How to use*

The script should be run from the directory containing the tar ball. I have saved this within the relevent **run** directory under a subdirectory **original_data**.

An example command to run the script:

```bash sample_seqs_prep.sh Taiwan_Vietnam_2016_IGFfolders Vietnam 52 Taiwan 332```

&nbsp;

### 2) DADA2_HPC_pipeline.R


#### *Requirements*

#### Arguments: 

An R-script "DADA2_args.R" stored in the relevent country folder and containing 11 variables:  

REV - reverse primer  
FWD - forward primer  
base_prefix - base name of sequences  
run_country - run directory/country directory  
root_path - your path to the directory containing all runs  
path - path to directory containing sequences, using the root_path  
path2 - path to filtN directory, using the root_path  
path_out - path to directory where output will be sent to, using root_path  
metadata_path - path to metadata file, using root_path  
cutadapt - path to cutadapt directory, using root_path  
unite.ref - path to UNITE database, using root_path  

#### Packages:

R:

dada2 (v.1.14.1)   
ShortRead (v.1.44.3)  
Biostrings (v.2.54.0)  
DECIPHER (v.2.14.0)  
phangorn (v.2.5.5)    
phyloseq (v.1.30.0)  

#### *How to use*

This is a script following the DADA2 ITS2 tutorial with some modifications to use only the forward reads and to run on HPC. Thanks to Phil Jervis for giving me his script to edit and to the creators of the dada2 pipeline (https://benjjneb.github.io/dada2/ITS_workflow.html).

No modifications should be required in order for this script to run with any R-script argument containing the correct variables and with the sequences in the correct format.

&nbsp;

### 3) run_DADA2_HPC_pipeline.sh

#### *How to use*

This script has been written to run on a PBS HPC so adjustments may need to be made to accommodate other machines.

The following command should be used to run the script, making sure to check your own releative paths:

```qsub -v "arg_path=run/country" ../../../../Generalised_Scripts/run_DADA2_HPC_pipeline.sh ```

Where ```run/country``` should be replaced with the relevent path. For example:

```qsub -v "arg_path=Taiwan_Vietnam_2016/Taiwan" ../../../../Generalised_Scripts/run_DADA2_HPC_pipeline.sh ```

