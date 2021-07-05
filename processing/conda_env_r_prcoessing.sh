#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:mem=10gb

# Author: Luke Goodyear leg19@imperial.ac.uk
# Script: env_manager.sh

# Desc: 1) cleans conda environment 
#	2) creates new environment and installs packages
#	3) cleans conda environment
# Arguments: None
# Date: Jun 2021


# load anaconda
module load anaconda3/personal

# remove any unused packages and caches
conda clean -a -y

# create new R conda environment
conda create -n r_processing -c conda-forge r-base=4.1.0
# activate new environment
source activate r_processing

# install required packages
conda install -c bioconda bioconductor-dada2=1.20.0 -y
conda install -c bioconda bioconductor-shortread=1.50.0 -y
conda install -c bioconda bioconductor-biostrings=2.60.0 -y
conda install -c bioconda bioconductor-phyloseq=1.36.0 -y
conda install -c conda-forge r-dplyr=1.0.6 -y
conda install -c bioconda cutadapt=2.8 -y

# remove any unused packages and caches after installations
conda clean -a -y

# view packages in environment as check
conda list

# deactivate environment
conda deactivate


# end of script
