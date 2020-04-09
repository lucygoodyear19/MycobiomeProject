# Mycobiome Project Boilerplates Repository

*Author: Lucy Goodyear*  
*Created: 04/02/20*

This repository contains all relevant boilerplate scripts for my analysis.

## Requirements

All code has been created for Mac so there may be a few differences in commands with respect to Linux. MacTeX or equivalent is required for LaTeX related files.

## Scripts

All folders contain the following scripts:

dada2_hpc_pipeline.R

This is a script following the dada2 IS2 tutorial with some modifications to use only the forward reads and to run on HPC. Thanks to Phil Jervis for giving me his script to edit and to the creators of the dada2 pipeline (https://benjjneb.github.io/dada2/ITS_workflow.html).

First check that all files are in the corect format. 

Use 

'''
gunzip IGF0006739_*/*001.fastq.gz
'''
