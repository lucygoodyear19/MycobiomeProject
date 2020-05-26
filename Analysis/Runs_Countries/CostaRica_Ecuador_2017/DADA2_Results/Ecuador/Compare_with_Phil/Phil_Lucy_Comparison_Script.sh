#!/bin/bash
# Author: Lucy Goodyear lucy.goodyear19@imperial.ac.uk
# Script: Phil_Lucy_Comparison_Script.sh
# Desc: Compares all DADA2 pipleine results to make sure they are the same
# Arguments: None
# Date: Mar 2020

echo "comparing filering output tables"
diff ~/Documents/MRes/MycobiomeProject/References/FromPhil/Ecuador_pipelineITS/filtN_Phil/filtering_output8.csv ~/Documents/MRes/MycobiomeProject/Analysis/Countries_Runs/Ecuador_CR/HPC_Results_Ecuador/filtering_output_lg.csv | wc -l

echo "comparing abundance tables"
diff ~/Documents/MRes/MycobiomeProject/References/FromPhil/Ecuador_pipelineITS/filtN_Phil/Ecuador\ Abundance\ Table\ new.txt ~/Documents/MRes/MycobiomeProject/Analysis/Countries_Runs/Ecuador_CR/HPC_Results_Ecuador/Abundance_Table_lg.txt | wc -l

echo "comparing metadata"
diff ~/Documents/MRes/MycobiomeProject/References/FromPhil/Ecuador_pipelineITS/filtN_Phil/Ecuador\ Metadata\ new.txt ~/Documents/MRes/MycobiomeProject/Analysis/Countries_Runs/Ecuador_CR/HPC_Results_Ecuador/Metadata_lg.txt | wc -l

echo "comparing taxa tables"
diff ~/Documents/MRes/MycobiomeProject/References/FromPhil/Ecuador_pipelineITS/filtN_Phil/Ecuador\ Tax\ Table\ new.txt ~/Documents/MRes/MycobiomeProject/Analysis/Countries_Runs/Ecuador_CR/HPC_Results_Ecuador/Tax_Table_lg.txt | wc -l

# end of script