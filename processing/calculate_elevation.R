###########################################################################################################
################################# Calculate elevations using geonames #####################################
###########################################################################################################


# Author: Luke Goodyear
# Date 14/01/2021


########################################### Set up ########################################################


# packages
library("rgbif")
library("tidyr")

# set path if running directly from here
path <- "~/Dropbox/luke/documents/academia/mres/mycobiome_project/analysis/runs_countries/Taiwan_Vietnam_2016/DADA2_Results/Taiwan/"

# import data
metadata <- read.csv(paste0(path, "metadata.csv"))
# save original csv for possible future use
write.csv(metadata, paste0(path, "metadata_pre_elevs.csv"))

# remove any samples with no lat/long data
metadata <- metadata %>% drop_na("Latitude")
# import geonames username
user <- Sys.getenv("GEONAMES_USER")


################################### calculate elevation ####################################################


# calculate elevations using geonames
elevs <- elevation(latitude = metadata$Latitude,
                   longitude = metadata$Longitude,
                   elevation_model = "srtm1", # calculated over 90m x 90m area
                   username = user) # username to connect to geonames database

# add elevation data to metadata as additional column
metadata$Elevation_m <- elevs$elevation_geonames

# save to csv
write.csv(metadata, paste0(path, "metadata.csv"))


## end of script
