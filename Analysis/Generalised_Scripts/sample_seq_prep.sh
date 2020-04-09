#!/bin/bash
# Author: Lucy Goodyear lucy.goodyear19@imperial.ac.uk
# Script: sample_seq_prep.sh
# Desc: Unzips sample sequence data and saves all fastq files to country directory
# Arguments:
# 1) root name of .tar.gz file containing data
# 2),4),6),8) Directory names of relevant countries (in plate order, up to 4)
# 3),5),7),9) Number of samples per country (in plate order, up to 4)
# Date: Mar 2020

# name variables
raw_data=$1

# untar folder
tar -xvf $raw_data.tar.gz

# cd to folder with sequence subdirectories
mv */*/*/*/*/*/*/ ./
cd IGF*/

# delete extra folders 
rm -rf *.*/

# optional remove a whole plate
#subdir_prefix=$(ls -1d *_1 | cut -d'_' -f1)
#echo ${subdir_prefix}
#for i in {1..96}; do
#    echo ${subdir_prefix}_${i}
#    rm -r ${subdir_prefix}_${i}
#done

# see how many arguments there are
count=$((($#-1)/2))
echo $count
# create directories per country to store all country sample sequence data
declare -a countries
for (( j=1; j<=$count; j++ )); do
    countries=("${countries[@]}" "$(eval echo \$$((2*j)))")
done
mkdir ${countries[@]}
echo ${countries[@]}

numsamples_cuml=("0")
offset=0
for (( j=1; j<=$count; j++ )); do
    numsamples=$(eval echo \$$(((2*j)+1)))
    numsamples_cuml=("${numsamples_cuml[@]}" "$(($numsamples+${offset} ))")
    offset=${numsamples_cuml[j]}
done
echo ${numsamples_cuml[@]}

# separate the samples into their respective countries
echo ${subdir_prefix}
for (( k=0; k<=${#countries[@]}; k++ )); do
	start=${numsamples_cuml[k]}
	for i in $(eval echo "{$((start+1))..${numsamples_cuml[((k+1))]}}"); do
		echo ${countries[k]}
		echo mv ${subdir_prefix}_${i} ${countries[k]}
  		mv ${subdir_prefix}_${i} ${countries[k]}
    	done
done

# unzip all fastq files in folders
#gunzip */IGF*/*001.fastq.gz

# initially to compare Phil's data and my data using this extraction method, below command was used:
# diff -r ../$numsamples1/Sample_Seqs/ Sample_Seqs_Phil/

# move all unzipped fastq files to their respective country directory in the main run directory
for country in ${countries[@]}; do
    mkdir ../../${country}/Sample_Seqs | mv ${country}/*/*.fastq ../../${country}/Sample_Seqs
    # create a directory to store filtered data files
    mkdir ../../${country}/Sample_Seqs/filtN # (must be created beforehand to allow path to be passed as external argument to DADA2 script)
done

## end of script
