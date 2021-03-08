###########################################################################################################
################################# Calculate elevations using geonames #####################################
###########################################################################################################

# Author: Lucy Goodyear
# Date 14/01/2021


# packages
library("rgbif")
library("tidyr")

# set path if running directly from here
path <- "~/Dropbox/luke/documents/academia/mres/mycobiome_project/analysis/runs_countries/Cameroon_SouthAfrica_2018/DADA2_Results/Cameroon/"

# import data
metadata <- read.csv(paste0(path, "metadata.csv"))
write.csv(metadata, paste0(path, "metadata_pre_elevs.csv"))

metadata <- metadata %>% drop_na("Latitude")
# import geonames username
user <- Sys.getenv("GEONAMES_USER")

# calculate elevation
elevs <- elevation(latitude = metadata$Latitude,
                   longitude = metadata$Longitude,
                   elevation_model = "srtm1", # calculated over 90m x 90m area
                   username = user) # username to connect to geonames database

# add elevation data to metadata as additional column
metadata$Elevation_m <- elevs$elevation_geonames

# save to csv
write.csv(metadata, paste0(path, "metadata.csv"))


## end of script
