#!/bin/bash
#PBS -l walltime=16:00:00
#PBS -l select=1:ncpus=1:mem=16gb
module load anaconda3/personal
source activate Renv
echo "R is about to run"
R --vanilla < $HOME/Project/Ecuador_pipelineITS/Ecuador_pipeline.R
echo "R has finished running"
# end of file