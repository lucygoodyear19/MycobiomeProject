######################################################################
#################### Arguments for DADA2 pipeline ####################
######################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1


### assign primers for import into R scripts

REV <- "CGCGCGCACGTTTCAHCGATGAAGAACGCAG"
FWD <- "GCATATHANTAAGSGGAGGTGACTGGCCGCCT"


### store base name of fastq files
base_prefix <- "160823_M03291_0044_000000000-ATG89_"


### assign path variables

# assign run/country to a variable
run <- "Taiwan_Vietnam_2016/"
country <- "Taiwan/"
# assign root paths to a variable
root_path <- "/rds/general/user/leg19/home/MRes/MycobiomeProject/Analysis/Runs_Countries/"

# assign full paths to arguments for imprt into R scripts
path <- paste0(root_path, run, country, "Sample_Seqs/")
path2 <- paste0(root_path, run, country, "Sample_Seqs/filtN/")
path_out <- paste0(root_path, run, country, "HPC_Results/")
metadata_path <- paste0(root_path, run, country, "metadata.csv")

cutadapt <- "/rds/general/user/leg19/home/anaconda3/envs/Renv/bin/cutadapt"
unite.ref <- "/rds/general/user/leg19/home/MRes/MycobiomeProject/Analysis/UNITE_database/sh_general_release_dynamic_s_02.02.2019.fasta"


## end of script
