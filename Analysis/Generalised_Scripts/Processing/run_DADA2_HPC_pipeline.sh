#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=120:mem=1440gb

# load required environment
module load anaconda3/personal
source activate Renv

# checks and prep
cutadapt --version
(cd "/rds/general/user/leg19/home/MRes/MycobiomeProject/Analysis/Runs_Countries/${arg_path}/Sample_Seqs/filtN/" && rm *)
echo ${arg_path}

echo "R is about to run DADA2 pipeline"
Rscript --vanilla $HOME/MRes/MycobiomeProject/Analysis/Generalised_Scripts/DADA2_HPC_pipeline.R "$HOME/MRes/MycobiomeProject/Analysis/Runs_Countries/${arg_path}/DADA2_args.R"
echo "R has finished running DADA2 pipeline"

echo "R is about to run filtering script on DADA2 results"
Rscript --vanilla $HOME/MRes/MycobiomeProject/Analysis/Generalised_Scripts/Data_Filtering_DADA2.R "$HOME/MRes/MycobiomeProject/Analysis/Runs_Countries/${arg_path}/DADA2_args.R"
echo "R has finished running filtering script on DADA2 results"

echo "R is about to run phylogenetic tree script"
Rscript --vanilla $HOME/MRes/MycobiomeProject/Analysis/Generalised_Scripts/phylogenetic_trees.R "$HOME/MRes/MycobiomeProject/Analysis/Runs_Countries/${arg_path}/DADA2_args.R"
echo "R has finished running phylogenetic tree script"

echo "All scripts have finsihed running. End."

# end of file
