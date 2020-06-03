######################################################################
####################### Arguments for Analysis #######################
######################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1


######################### SET UP FOR EASE ############################


# only change these for ease
root <- "~/Documents/MRes/MycobiomeProject/Analysis/"
run <- "Taiwan_Vietnam_2016/"
country <- "Taiwan"
year <- "2016"


##################### INPUTS REQUIRED FOR SCRIPTS ####################

# paths
dada2_data_path <- paste0(root, 
                          "Runs_Countries/", 
                          run, 
                          "DADA2_Results/", 
                          country, 
                          "/physeqob_DADA2_complete.rds")
esto_data_path <- paste0(root, 
                    "Runs_Countries/", 
                    run, 
                    "Esto_Results/physeqob_esto.rds")
results_path <- paste0(root, 
                  "Results/", 
                  country, 
                  "_", 
                  year)

# vector of datasets
datasets <- c("dada2", "esto")

### NCBI database?
#ncbi.ref <- paste(root, "NCBI_database")


## end of script
