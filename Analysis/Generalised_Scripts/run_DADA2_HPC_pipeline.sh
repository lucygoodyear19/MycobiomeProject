#!/bin/bash
#PBS -l walltime=16:00:00
#PBS -l select=1:ncpus=48:mem=124gb
module load anaconda3/personal
source activate Renv
cutadapt --version
(cd "/rds/general/user/leg19/home/MRes/MycobiomeProject/Analysis/Runs_Countries/${arg_path}/Sample_Seqs/filtN/" && rm *)
echo ${arg_path}
echo "R is about to run"
Rscript --vanilla $HOME/MRes/MycobiomeProject/Analysis/Generalised_Scripts/DADA2_HPC_pipeline.R "$HOME/MRes/MycobiomeProject/Analysis/Runs_Countries/${arg_path}/DADA2_args.R"
echo "R has finished running"
# end of file
