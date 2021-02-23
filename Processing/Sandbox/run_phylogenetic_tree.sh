#!/bin/bash
#PBS -l walltime=100:00:00
#PBS -l select=1:ncpus=8:mem=96gb

# load required environment
module load anaconda3/personal
source activate R_tree_env

# checks and prep
echo ${arg_path}

echo "R is about to run phylogenetic tree script"
Rscript --vanilla $HOME/MRes/MycobiomeProject/Analysis/Generalised_Scripts/Processing/phylogenetic_trees.R "$HOME/MRes/MycobiomeProject/Analysis/Results/${arg_path}/analysis_args.R"
echo "R has finished running phylogenetic tree script"

## end of script
