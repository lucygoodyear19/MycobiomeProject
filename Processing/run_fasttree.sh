#!/bin/bash
# Author: Lucy Goodyear lucy.goodyear19@imperial.ac.uk
# Script: run_fasttree.sh
# Desc: Using fasttree to create a rough phylogenetic tree
# Arguments:
# 1) file containing sequences (fasta or interleaved phylip format)
# Date: Jul 2020

# get fasttree c code
wget "http://meta.microbesonline.org/fasttree/FastTree.c"

# compile fasttree
gcc -DNO_SSE -O3 -funroll-loops -Wall -o FastTree FastTree.c -lm

# run fasttree on input file
./FastTree -gtr -nt < $1 > fasttree_tree

## end of script