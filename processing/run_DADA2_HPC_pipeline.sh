#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=120:mem=1440gb

# Author: Luke Goodyear leg19@imperial.ac.uk
# Script: run_DADA2_HPC_pipeline.sh
# Desc: runs DADA2 pipeline (for use when submitting batch job)
# Arguments:
# 1) name of run as string
# 2) name of country as string
# Date: Mar 2020


# load required environment
module load anaconda3/personal
source activate r_processing

# checks and prep
cutadapt --version
(cd "/rds/general/user/leg19/home/mres/mycobiome_project/analysis/runs_countries/${run}/original_data/${country}/sample_seqs/filtN/" && rm *)
echo ${run}
echo ${country}

echo "R is about to run DADA2 pipeline"
Rscript --vanilla $HOME/mres/mycobiome_project/analysis/scripts/processing/DADA2_pipeline.R "$HOME/mres/mycobiome_project/analysis/runs_countries/${run}/DADA2_results/${country}/DADA2_args.R"
echo "R has finished running DADA2 pipeline"

echo "R is about to run filtering script on DADA2 results"
Rscript --vanilla $HOME/mres/mycobiome_project/analysis/scripts/processing/DADA2_data_filtering.R "$HOME/mres/mycobiome_project/analysis/runs_countries/${run}/DADA2_results/${country}/filtering_args.R"
echo "R has finished running filtering script on DADA2 results"

echo "All scripts have finsihed running. End."

# end of file
