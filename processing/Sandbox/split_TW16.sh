#!/bin/bash
# Author: Lucy Goodyear lucy.goodyear19@imperial.ac.uk
# Script: split_TW16.sh
# Desc: Moves samples into directories by plate number
# Arguments:
# 2) csv file containing Sample_ID (matching sample fastq names) and Plate (number)
# Date: Jul 2020


# activate conda environment containing csvcut
module load anaconda3/personal
source activate bash_env

# name argument variables
csv=$1

echo "Moving samples into subdirectories by plate number..."

first_line=true # to ignore reading headings in csv parse

# while loop to split and move samples into plate directories
while IFS=, read -r Sample_ID Plate Origin # read the three relevant columns with each value separated by ,
do
    if [ "$first_line" = true ] ; then
        first_line=false
        echo "Skipping first line..."
        continue
    fi
    
    mkdir -p ../../$Plate/Sample_Seqs/filtN/ # create plate directory if doesn't exist
    cp "${Sample_ID}"* ../../$Plate/Sample_Seqs/ # move sample to plate

done < <(csvcut -c Sample_ID,Plate,Origin ${csv}.csv) # only parse the relevant columns


## end of script
