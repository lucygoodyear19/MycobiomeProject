#!/bin/bash
# Author: Lucy Goodyear lucy.goodyear19@imperial.ac.uk
# Script: plate_preparation.sh
# Desc: Unzips plate data and saves all fastq files to one directory
# Arguments: directory containing plate data
# Date: Mar 2020

# untar folder
tar -xf $1

# unzip all fastq files in folder
gunzip IGF*/*001.fastq.gz

# make a single directory to store all plate data
mkdir Plates

# move all fastq files to Plates directory
mv IGF*/*.fastq Plates/

# then separate out plates from different countries into different directories
# initially to compare Phil's plat data and plate data using this extraction method, below command was used:
# diff -r Plates/Ecuador_Plates/ Plates_Phil/

# end of script