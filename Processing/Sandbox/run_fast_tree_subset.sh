#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=32:mem=124gb

# load required environment
module load anaconda3/personal
source activate Renv

# checks and prep
echo ${arg_path}

echo "R is about to run phylogenetic tree script"
Rscript --vanilla $HOME/MRes/MycobiomeProject/Analysis/Generalised_Scripts/Processing/fast_tree_subset.R  "$HOME/MRes/MycobiomeProject/Analysis/Results/${arg_path}/analysis_args.R"
echo "R has finished running phylogenetic tree script"

## end of script
