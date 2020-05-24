#################
### Filtering ###
#################

# MiSeq plate plan
setwd("~/Documents/MRes/MycobiomeProject/Analysis/Runs_Countries/Taiwan_Vietnam_2016/Taiwan/Filtering")
plates <- read.csv("../TW16_plate_data.csv", stringsAsFactor=FALSE, header=F)
names(plates) <- c("barcode","sample","DNAqual","plate")
# assign source to samples
plates$location                                 <- "Vietnam"
plates$location[grep("TW16", plates$sample)]    <- "Taiwan"
plates$location[grep("NTC", plates$sample)]     <- "NC_PCR"
plates$location[grep("posC", plates$sample)]    <- "posC"
plates$location[grep("mock", plates$sample)]    <- "mock"
table(plates$sample, plates$location);table(plates$location)

# read in completed information for OTUs combined with OTU table
otutabinfo <- read.csv("OTU_seq_species_taxonomy_OTUtable_samples.csv", header=T, stringsAsFactor=FALSE, check.names=FALSE)

#-------FILTERED FILES--------#

###Post-pipeline, per-filter OTU table has 49085 OTUs x 384 samples (plus 22 info colunns)

## (2.1) Remove non-fungal OTUs
#setwd("C:/Users/jenny/Dropbox/Jen's desktop/Results/sequencing/MiSeq/3-Vietnam_Taiwan/2x300bp/7-filtered_files/pear")
F_only <- read.csv("fungalOTUsamp.csv", header=T, stringsAsFactor=FALSE, check.names = F) #47396 OTUs x 384 samples

## (2.4a) Remove posC and mock community OTUs
#setwd("C:/Users/jenny/Dropbox/Jen's desktop/Results/sequencing/MiSeq/3-Vietnam_Taiwan/2x300bp/7-filtered_files/pear")
poscmock <- read.csv("posc_mock_samp.csv", header=T, stringsAsFactor=FALSE, check.names = F) #42154 OTUs x 384 samples

## (2.5a) Filter samples by removing remaining reads in 4 x control wells --> Remove OTUs and samples without reads
#setwd("C:/Users/jenny/Dropbox/Jen's desktop/Results/sequencing/MiSeq/3-Vietnam_Taiwan/2x300bp/7-filtered_files/pear")
ntc <- read.csv("A-sampfilt_subtractNTC.csv", header=T, stringsAsFactor=FALSE, check.names = F) #42154 OTUs x 367 samples (-control wells)
#setwd("C:/Users/jenny/Dropbox/Jen's desktop/Results/sequencing/MiSeq/3-Vietnam_Taiwan/2x300bp/7-filtered_files/pear")
twfinal234 <- read.csv("Taiwanonly_plate2-4_filteredOTUtable.csv", header=T, stringsAsFactor=FALSE, check.names = F)
swabfinal234 <- read.csv("Taiwanswabsonly_plate2-4_filteredOTUtable.csv", header=T, stringsAsFactor=FALSE, check.names = F)

###################################################################################################################

### (2.1) Remove non-fungal reads

# required: otutabinfo
# check which columns are info vs samples
names(otutabinfo)

# first filter step is to remove TW16032 and mock1 as these were mixed up during setup!
otutabinfo <- otutabinfo[!names(otutabinfo) %in% c("TW16032","mock1")]

# view the different taxa in the OTU table
table(otutabinfo$kingdom)

# remove non-fungal OTUs
# filter by kingdom equal to fungi
fungi <- otutabinfo[otutabinfo$kingdom == "Fungi",]
# separate out OTU info (first 22 columns) and sample info (all other columns)
info <- fungi[,1:22]
fungisamptot <- fungi[,23:ncol(fungi)]

# remove samples without reads in fungal OTUs
# sum all read columns and save as vector
sampsum <- colSums(fungisamptot)
# filter by sums not equal to zero
fungisamp <- fungisamptot[,sampsum != 0]

#combine with OTU info
F_only <- cbind(info, fungisamp)
#setwd("C:/Users/jenny/Dropbox/Jen's desktop/Results/sequencing/MiSeq/3-Vietnam_Taiwan/2x300bp/7-filtered_files/pear")
write.table(F_only, "fungalOTUsamp.csv", sep=",", row.names = F) #47396 OTUs x 382 samples

###################################################################################################################

### (2.4) Remove posC and mock OTUs

# required files: identity
# mocks & posC

# find out which OTUs for mock community need removing
# establish that mock community OTUs are appearing in samples in read numbers unlikely to be cross-contamination so will not be removed

# find out which OTUs for posC need removing
# create subset containing only posC and the OTU and species columns
posc <- cbind(F_only$otuid,F_only$species,F_only[grep("posC",names(F_only))])
names(posc) <- c("otu","species","posC1","posC2","posC3","posC4") # rename columns
# find out species names of the 4 posC (should all be same)
poscspec1 <- as.character(posc$species[posc$posC1 > 1000])
poscspec2 <- as.character(posc$species[posc$posC2 > 1000])
poscspec3 <- as.character(posc$species[posc$posC3 > 1000])
poscspec4 <- as.character(posc$species[posc$posC4 > 1000])

uniqposcspec <- unique(c(poscspec1, poscspec2, poscspec3, poscspec4)) # 1 species
# subset OTUs that are the same as the species of posC and count them
poscotu <- as.character(posc$otu[posc$species %in% unique(uniqposcspec)])
length(poscotu)

# remove posC OTUs
identityposc <- F_only[!F_only$otuid %in% poscotu,] 

#-----------

#Find out which OTUs for mock need removing
# create subset containing only posC and the OTU and species columns
mock <- cbind(identityposc$otuid,identityposc$species,identityposc[grep("mock",names(identityposc))]) 
names(mock) <- c("otu","species","mock2","mock3","mock4") # rename columns

# find out names of the species names in the mocks
mockspec2 <- as.character(mock$species[mock$mock2 > 1000])
mockspec3 <- as.character(mock$species[mock$mock3 > 1000])
mockspec4 <- as.character(mock$species[mock$mock4 > 1000])

uniqmockspec <- unique(c(mockspec2,mockspec3,mockspec4)) # 5 species
# subset OTUs that are the same as the species of mocks and count them
mockotu  <- as.character(mock$otu[mock$species %in% unique(uniqmockspec)]) #1020 OTUs
length(mockotu)

# remove mock OTUs from the dataframe we just removed posC OTUs from
identityposcmock <- identityposc[!identityposc$otuid %in% mockotu,]

# check which columns are info vs samples
names(identityposcmock)

# separate out OTU info (first 22 columns) and sample info (all other columns)
info <- identityposcmock[,1:22]
posmocksamptot <- identityposcmock[,23:ncol(identityposcmock)]
# remove samples without reads in remaining OTUs
# sum all read columns and save as vector
sampsum <- colSums(posmocksamptot)
# filter by sums not equal to zero
posmocksamp <- posmocksamptot[,sampsum != 0]

#combine with OTU info
poscmock <- cbind(info,posmocksamp)
#setwd("C:/Users/jenny/Dropbox/Jen's desktop/Results/sequencing/MiSeq/3-Vietnam_Taiwan/2x300bp/7-filtered_files/pear")
write.table(poscmock, "posc_mock_samp.csv", sep=",", row.names=F) #42154 OTUs x 382 samples

##!!Firework plot does not indicate high reads for P.citrinum in this run so all mock species should be excluded

###################################################################################################################

## (2.5a) Filter samples by removing remaining reads in 4 x control wells --> Remove OTUs and samples without reads

# required files: poscmock, plates
names(poscmock)

# separate out OTU info (first 22 columns) and plates
info <- poscmock[,1:22]

# split samples by plate number and create new column summing NTC and PosC reads per sample
plate1 <- poscmock[,names(poscmock) %in% c(plates$sample[plates$plate == "plate1"])]
plate1$sum1 <- plate1$NTC1 + plate1$posC1
head(plate1[c("NTC1", "posC1","sum1")]) # view data
length(names(plate1)) # count samples 
length(which(plate1$sum1 != 0)) # count OTUs with NTC + posC reads not equal to 0
sum(plate1$sum1) # sum overall NTC + posC reads

plate2 <- poscmock[,names(poscmock) %in% c(plates$sample[plates$plate == "plate2"])]
plate2$sum2 <- plate2$NTC2 + plate2$posC2 + plate2$mock2
head(plate2[c("NTC2", "posC2","mock2","sum2")])
length(names(plate2))
length(which(plate2$sum2 != 0))
sum(plate2$sum2)

plate3 <- poscmock[,names(poscmock) %in% c(plates$sample[plates$plate == "plate3"])]
plate3$sum3 <- plate3$NTC3 + plate3$posC3 + plate3$mock3
head(plate3[c("NTC3", "posC3","mock3","sum3")])
length(names(plate3))
length(which(plate3$sum3 != 0))
sum(plate3$sum3)

plate4 <- poscmock[,names(poscmock) %in% c(plates$sample[plates$plate == "plate4"])]
plate4$sum4 <- plate4$NTC4 + plate4$posC4 + plate4$mock4
head(plate4[c("NTC4", "posC4","mock4","sum4")])
length(names(plate4)); length(which(plate4$sum4 != 0))
sum(plate4$sum4)

# remove number of reads in NTC + posC OTUs from corresponding sample OTUs for each plate
plate1x <- as.data.frame(t(apply(plate1, 1, function(x) x - x["sum1"])))
plate1x[plate1x < 0] <- 0 # set any read value less than 0 to 0
sum(plate1x$sum1) # check: should be 0

plate2x <- as.data.frame(t(apply(plate2, 1, function(x) x - x["sum2"])))
plate2x[plate2x < 0] <- 0
sum(plate2x$sum2) # should be 0

plate3x <- as.data.frame(t(apply(plate3, 1, function(x) x - x["sum3"]))) 
plate3x[plate3x < 0] <- 0
sum(plate3x$sum3) # should be 0

plate4x <- as.data.frame(t(apply(plate4, 1, function(x) x - x["sum4"])))
plate4x[plate4x < 0] <- 0
sum(plate4x$sum4) # should be 0

# concatenate info and samples with NTC/posC reads deducted
all <- cbind(info, plate1x, plate2x, plate3x,plate4x)
all2 <- all[!names(all) %in% c("sum1","sum2","sum3","sum4")] # remove newly created sum columns

# check which columns are info vs samples
names(all2)

# separate out OTU info (first 22 columns) and sample info (all other columns)
info <- all2[,1:22]
ntcsamptot <- all2[,23:ncol(all2)]
# remove samples without reads in remaining OTUs
# sum all read columns and save as vector
sampsum <- colSums(ntcsamptot)
# filter by sums not equal to zero
ntcsamp <- ntcsamptot[,sampsum != 0]
names(ntcsamp) # check all mocks, NTCs and posCs have been removed

#combine with OTU info
ntcotutab <- cbind(info,ntcsamp)
ntc <- ntcotutab[!names(ntcotutab) %in% c("mock1","mock2","mock3","posC1","posC2","posC3","NC_swab1","NC_swab2","NC_swab3","NC_PCR1","NC_PCR2","NC_PCR3")] # remove these columns
#setwd("C:/Users/jenny/Dropbox/Jen's desktop/Results/sequencing/MiSeq/3-Vietnam_Taiwan/2x300bp/7-filtered_files/pear")
write.table(ntc, "A-sampfilt_subtractNTC.csv", sep=",", row.names=F) #42154 OTUs x 367 samples (-control wells)

#-----

## plates 2-4/Taiwan only

# required files: ntc & plates

# remove plate 1 reads from Tawian data because of contamination

# check which columns are info vs samples
names(ntc)

# separate out OTU info (first 22 columns) and sample info (all other columns)
info <- ntc[,1:22]
samp <- ntc[,23:ncol(ntc)]

# subset samples from plate 1
plate1samp <- plates$sample[plates$plate == "plate1"]
# remove all plate 1 samples from dataframe
taiwansamp <- samp[!names(samp) %in% plate1samp]

# remove OTUs without reads in Taiwan samples
twotusum <- rowSums(taiwansamp) # sum over each OTU
twall <- cbind(info,taiwansamp) # recombine dataframe with info
twfinal <- twall[twotusum != 0,] # remove all OTUs with 0 reads
#setwd("C:/Users/jenny/Dropbox/Jen's desktop/Results/sequencing/MiSeq/3-Vietnam_Taiwan/2x300bp/7-filtered_files/pear")
write.csv(twfinal, "Taiwanonly_plate2-4_filteredOTUtable.csv",row.names=FALSE)

#-----

## plates 2-4/Taiwan swabs only

# required files: twfinal

# check which columns are info vs samples
names(twfinal)

# separate out OTU info (first 22 columns) and sample info (all other columns)
info <- twfinal[,1:22]
samp <- twfinal[,23:ncol(twfinal)]

# search for all samples with "c" in the name to find swabs
swabsamp <- samp[,-grep("c",names(samp))]

# remove OTUs without reads in Taiwan swabs
swabotusum <- rowSums(swabsamp) # sum over each OTU
swaball <- cbind(info,swabsamp) # recombine dataframe with info
swabfinal <- swaball[swabotusum != 0,] # remove all OTUs with 0 reads
#setwd("C:/Users/js911/Dropbox/Jen's desktop/Results/sequencing/MiSeq/3-Vietnam_Taiwan/2x300bp/7-filtered_files/pear")
write.csv(swabfinal, "Taiwanswabsonly_plate2-4_filteredOTUtable.csv",row.names=FALSE)

