#!/bin/bash
#PBS -l walltime=100:00:00
#PBS -l select=1:ncpus=8:mem=96gb
module load anaconda3/personal
source activate R_beta_env
echo "R is about to run beta diversity script"
Rscript --vanilla $HOME/MRes/MycobiomeProject/Analysis/Generalised_Scripts/Analysis/beta_diversity.R "$HOME/MRes/MycobiomeProject/Analysis/Results/${arg_path}/analysis_args.R"
echo "R has finished running beta diversity script"
# end of file
