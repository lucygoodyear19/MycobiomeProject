#!/bin/bash
# Author: Lucy Goodyear lucy.goodyear19@imperial.ac.uk
# Script: sample_seq_prep.sh
# Desc: Unzips sample sequence data and saves all fastq files to country directory
# Arguments:
# 1) root name of .tar.gz file containing data
# 2) csv file containing Sample_ID (matching directory IDs), Plate (number) 
#    and Origin (country matching country directory or Mock or PosC or NC) columns
# Date: Mar 2020

# activate conda environment containing csvcut
module load anaconda3/personal
source activate bash_env

# check for correct number of arguments
if [[ $# -lt 2 || $# -gt 2 ]] ; then
    to_print="
    Incorrect number of arguments supplied. 
    Requires 2 arguments:
    1) root name of .tar.gz file containing data
    2) csv file containing Sample_ID (matching directory IDs), Plate (number) 
       and Origin (country matching country directory or Mock or PosC or NC) columns
    "
    echo "$to_print"
    exit 1 # exit script
fi

# name argument variables
raw_data=$1
csv=$2

# untar folder
tar -xvf $raw_data.tar.gz

# cd to folder with sequence subdirectories
mv */*/*/*/*/*/*/ ./
ls
cd IGF*/

# delete extra folders 
#rm -rf *.*/

# declare arrays to store the relevant plate numbers per country
# and the Sample_IDs for the controls present per plate
declare -A plates_per_country
declare -A controls_per_plate

first_line=true # to ignore reading headings in csv parse

# while loop to split and move samples into country and to store plate numbers per country
# and control IDs per plate
while IFS=, read -r Sample_ID Plate Origin # read the three relevant columns with each value separated by ,
do
    echo "Moving country samples into country folders"
    if [ "$first_line" = true ] ; then
        first_line=false
        echo "Skipping first line..."
        continue
    fi 

    if [[ "$Origin" == "NC" || "$Origin" == "Mock" || "$Origin" == "PosC" ]] ; then
        echo "Origin is $Origin, skipping..."
        controls_per_plate[$Plate]="${controls_per_plate[$Plate]} 
        $Sample_ID" # store control ID and plate number (on a separate line)
        continue

    fi  

    mkdir -p $Origin/ # create country directory if doesn't exist
    mv $Sample_ID/ $Origin/ # move sample to country
    plates_per_country[$Origin]="${plates_per_country[$Origin]} 
    $Plate" # store plate number with country for that sample
    plates_per_country[$Origin]="$(echo "${plates_per_country[$Origin]}" | sort -un)" # delete if duplicate

done < <(csvcut -c Sample_ID,Plate,Origin ../${csv}.csv) # only parse the relevant columns

# print country and plate keys as a checks
echo "Country keys: ${!plates_per_country[@]}"
echo "Plate keys: ${!controls_per_plate[@]}"

# move controls to relevant countries corresponding to plate number
for country in ${!plates_per_country[@]}; do
	plates="${plates_per_country[$country]}" # store plate numbers for current value of iterator
	echo "Corresponding plates per country: $country -> ${plates}"
    for plate in $plates; do
        controls="${controls_per_plate[$plate]}" # store control IDs for current plate number
        echo "Controls present per plate: $plate -> ${controls}"
        cp -r $controls $country/ # copy paste each control to each relevant country corresponding to plate number
    done
done

# unzip all fastq files in folders
#gunzip */IGF*/*001.fastq.gz

# initially to compare Phil's data and my data using this extraction method, below command was used:
# diff -r ../$numsamples1/Sample_Seqs/ Sample_Seqs_Phil/

# move all unzipped fastq files per country to a new directory within each country directory
cd ../
for country in ${!plates_per_country[@]}; do
    echo "Moving all sample fastq files to resepctive run/country directories"
    mkdir -p ${country}/Sample_Seqs | mv IGF*/${country}/*/*.fastq ${country}/Sample_Seqs
    # create a directory to store filtered data files
    mkdir -p ${country}/Sample_Seqs/filtN # (must be created beforehand to allow path to be passed as external argument to DADA2 script)
done

# remove all unnessecary directories and files
rm -r IGF*/ work/

## end of script
