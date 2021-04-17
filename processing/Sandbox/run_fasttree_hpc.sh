#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=32:mem=124gb

# Author: Lucy Goodyear lucy.goodyear19@imperial.ac.uk
# Script: run_fasttree.sh
# Desc: Using fasttree to create a rough phylogenetic tree
# Arguments:
# 1) file containing sequences (fasta or interleaved phylip format)
# Date: Jul 2020

# get fasttree c code
#wget "http://meta.microbesonline.org/fasttree/FastTree.c"
wget "http://meta.microbesonline.org/fasttree/FastTree"

# compile fasttree
#gcc -DNO_SSE -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm

# make FastTree executable
chmod +x FastTree

pwd
ls -l

# run fasttree on input file
./FastTree -gtr -nt < "${arg}" > /rds/general/user/leg19/home/MRes/MycobiomeProject/Analysis/Results/Global_14/fasttree_tree


## end of script
