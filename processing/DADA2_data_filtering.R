##############################################################################
###################### Data filtering for DADA2 results ######################
##############################################################################


# Author: Luke Goodyear (leg19@ic.ac.uk)
# Version: 0.0.1

# clear worksapce
rm(list=ls())


################################## Set up ####################################


# on local computer: store all console output to an output file
#sink("DADA2_data_filtering_output.log", type=c("output", "message"))

# import arguments to run script on specific country data
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) # setup to accept arguments from command line
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied in the form of an R-script containing the following arguments: 
       1) run = name of run followed by /
       2) country = name of country
       3) root_path = path until directory above run folowed by /
       4) metadata_path = path to metadata csv
       5) path_in = path to input data
       6) path_out = path to output phyloseq object results to")
}
# load arguments into script
source(args)
# print arguments as check
print("Arguments:")
print(paste0("run: ", run))
print(paste0("country: ", country))
print(paste0("root_path: ", root_path))
print(paste0("metadata_path: ", metadata_path))
print(paste0("path_in: ", path_in))
print(paste0("path_out: ", path_out))

# load packages
print("Loading required packages...")
library("phyloseq")
library("dplyr")


################################# Load data #####################################


print("Loading data...")

# load taxa data
tax <- read.table(paste0(path_in, "Taxa_Table.txt"),
                  stringsAsFactors = F,
                  header = T)
# rename columns
colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# load abundance data
seqtab <- read.table(paste0(path_in, "Abun_Table_sample_row.txt"),#abun_table_samplesinrows.txt"), 
                     stringsAsFactors = F,
                     header = T)

# load plate data
plates <- read.csv(paste0(root_path, run, "Plate_Data.csv"), 
                   stringsAsFactor=F, 
                   header=T)

# load metadata
metadata <- read.csv(metadata_path,
                     stringsAsFactor=F, 
                     header=T)

# add column for Bd +ve/-ve to metadata
metadata$Bd <- NA
for (i in (1:nrow(metadata))){
  if (is.na(metadata$Bd_GE[i])){
    metadata$Bd[i] <- NA
  }
  else {
    if (metadata$Bd_GE[i] >= 0.1){
    metadata$Bd[i] <- 1
    } 
    if (metadata$Bd_GE[i] < 0.1){
    metadata$Bd[i] <- 0
    }
  }
}

# set as factor for plotting
metadata$Bd <- as.factor(metadata$Bd)


######################### remove PosC, Mocks and NC ###########################


#################### Set up


print("Performing set up...")

# change sample names in plate data to match processed data format
for (patt in 1:length(plate_patterns)) {
  plates$Sample_Name <- gsub(plate_patterns[patt], meta_patterns[patt], plates$Sample_Name)
}

# modification for Taiwan 2017 run to account for swapping of mock and posC labels
if (run == "Taiwan_2017/") {
  for (patt in 1:length(plate_patterns)) {
    colnames(tax) <- gsub(plate_patterns[patt], meta_patterns[patt], colnames(tax))
    rownames(seqtab) <- gsub(plate_patterns[patt], meta_patterns[patt], rownames(seqtab))
  }
}

# modification for Taiwan 2016 run to account for naming inconsistencies
if (run == "Taiwan_Vietnam_2016/" && country == "Taiwan") {
  # change sample names in plate data to match processed data format
  rownames(seqtab) <- gsub("-", "", rownames(seqtab))
  # replace sample indices with sample names in seqtab_t
  for (samp in 1:nrow(seqtab)) {
    rownames(seqtab)[samp] <- plates$Sample_Name[grep(rownames(seqtab)[samp], plates$indices)]
  }
}

# set rownames of taxonomy table to match column names of abundance table
rownames(tax) <- colnames(seqtab)

# transpose abundance table
seqtab_t <- as.data.frame(t(seqtab))

# subset taxonomy table
species_by_seq <- as.data.frame(cbind(rownames(tax), tax$Genus, tax$Species))
names(species_by_seq) <- c("Sequence", "Genus", "Species")

# merge subsetted taxa table and the abundance table
species_abun <- merge(species_by_seq, seqtab_t, by.x = "Sequence", by.y = "row.names")

# create vector of plate numbers that the country was present on
if(is.na(plate_nos)) {
  plate_nos <- c()
  for (samp in 1:nrow(plates)){
    if (plates$Origin[samp] == country) {
      plate_nos <- c(plate_nos, plates$Plate[samp])
    } 
    plate_nos <- unique(plate_nos)
  }
}


######################### PosC


print("Locating posC ASVs...")

# subset merged data frame to include only positive controls
posc <- cbind(species_abun$Sequence, paste(species_abun$Genus, species_abun$Species),species_abun[grep("posC",names(species_abun))])
# rename columns
names(posc)[1:2] <- c("Sequence", "Genus_Species")

# find posC species
posc_ls <- c()
for (plate in plate_nos) {
  assign(paste0("poscspec", plate), as.character(posc$Genus_Species[posc[[paste0("posC", plate)]] > 10]))
  posc_ls <- c(posc_ls, get(paste0("poscspec", plate)))
}

# find unique posC species
uniqposcspec <- unique(posc_ls)
print(paste("PosC species is", uniqposcspec))

# subset sequences that are the same as the species of posC and count them
posc_seqs <- as.character(posc$Sequence[posc$Genus_Species %in% unique(uniqposcspec)])
print(paste0("Number of ASVs that are same as posC species: ", length(posc_seqs)))

# remove PosC sequences
seqtab <- seqtab[,!names(seqtab) %in% posc_seqs]


###################### Mocks


print("Locating mock ASVs...")

# subset merged data frame to include only positive controls
mocks <- cbind(species_abun$Sequence, paste(species_abun$Genus, species_abun$Species),species_abun[grep("mock",names(species_abun))])
# rename columns
names(mocks)[1:2] <- c("Sequence", "Genus_Species")

# find mock species
mock_ls <- c()
for (plate in plate_nos) {
  assign(paste0("mockspec", plate), as.character(mocks$Genus_Species[mocks[[paste0("mock", plate)]] > 10]))
  mock_ls <- c(mock_ls, get(paste0("mockspec", plate)))
}

# find unique mock species
uniqmockspec <- unique(mock_ls)
print(paste("Mock species are", uniqmockspec))

# subset sequences that are the same as the species of the mocks and count them
mock_seqs <- as.character(mocks$Sequence[mocks$Genus_Species %in% unique(uniqmockspec)])
print(paste0("Number of ASVs that are mock species: ", length(mock_seqs)))

# remove mock sequences
seqtab <- seqtab[,!names(seqtab) %in% mock_seqs]
# remove mocks from samples
seqtab <- seqtab[!rownames(seqtab) %in% c("mock1", "mock2", "mock3", "mock4"),]


##################### Negative controls


print("Removing sequences found in negative controls...")

# create vector of negative control sample names
NCs <- c()
for (samp in 1:nrow(plates)){
  if (grepl("NC", plates$Origin[samp])) {NCs <- c(NCs, plates$Sample_Name[samp])} 
}


################### Removing control reads


# define function to remove posC and NC reads from samples by plate
# returns a dataframe containing samples from the specifed plate, in the same form as df_abun 
filter_controls <- function(df_abun, df_plates, plate_no) {
  
  # subset by plate number
  plate <- df_abun[,names(df_abun) %in% c(df_plates$Sample_Name[df_plates$Plate == plate_no])]
  # create new column summing NC and PosC reads per sample
  if (ncol(plate[,colnames(plate) %in% NCs, drop=FALSE]) == 0) { # check to see if there are any negative controls in the plate
    plate$sum <- plate[,(paste0("posC", plate_no))]
  } else {
    plate$sum <- rowSums(plate[,c(colnames(plate[,colnames(plate) %in% NCs, drop=FALSE]), paste0("posC", plate_no))])
  }
  
  # print summary statements
  print(paste0("Number of samples in plate ", plate_no, ": ", length(names(plate))-1)) # count samples - sum column
  print(paste0("Number of ASVs with NC/PosC reads not equal to 0: ", length(which(plate$sum != 0)))) # count ASVs with NC + posC reads not equal to 0
  print(paste0("Number of NC and posC reads overall before: ", sum(plate$sum))) # sum overall NC + posC reads
  
  # remove number of reads in NC + posC ASVs from corresponding ASVs for each sample
  platex <- as.data.frame(t(apply(plate, 1, function(x) x - x["sum"])))
  # set any negative values (create by the above) to 0
  platex[platex < 0] <- 0
  
  # print a check statement
  print(paste0("CHECK. Number of NC and posC reads overall after: ", sum(platex$sum), " (should be 0)")) # should be 0
  # remove sum column
  platex <- subset(platex, select=-sum)
  
  return(platex)
}

# transpose filtered abundance table
seqtab_t <- as.data.frame(t(seqtab))
# run the function on all plates
seqtab_filtered <- data.frame(matrix(ncol=ncol(seqtab),nrow=0))
colnames(seqtab_filtered) <- colnames(seqtab)
for (plate in plate_nos) {
  assign(paste0("plate", plate), t(filter_controls(seqtab_t, plates, plate)))
  # combine the plates into one dataframe
  seqtab_filtered <- rbind(seqtab_filtered, get(paste0("plate", plate)))
}

# remove empty all samples (including posC and NC) from data frame
# sum abundances by sample
seqtab_filtered <- cbind(seqtab_filtered, apply(seqtab_filtered, 1, sum))

# count how many samples have no reads
print(paste0("Number of sample with no reads after filtering: ", length(which(seqtab_filtered[ncol(seqtab_filtered),] == 0))))

print("Removing empty samples...")
# save empty samples as a vector
to_remove <- which(seqtab_filtered[,ncol(seqtab_filtered)] == 0)
# remove empty samples from dataframe
seqtab_filtered <- seqtab_filtered[!rownames(seqtab_filtered) %in% names(to_remove),]
# remove sum column
seqtab_filtered <- seqtab_filtered[,-ncol(seqtab_filtered)]


######################### Tidy up


# modification for one particular Taiwan run to remove all Bsal and cultures
if (run == "Taiwan_Vietnam_2016/" && country == "Taiwan") {
  require(tibble) # for `rownames_to_column` and `column_to_rownames`
  require(dplyr)
  
  # separate out Bsal from taxa table
  bsal <- tax %>%
    rownames_to_column('Seq') %>%
    filter(Species == "s__salamandrivorans") %>%
    column_to_rownames('Seq')

  # save Bsal sequences
  bsal_seqs <- rownames(bsal)

  # remove all Bsal sequences from abundance table
  seqtab_filtered <- as.data.frame(seqtab_filtered[,!colnames(seqtab_filtered) %in% bsal_seqs,])

  # remove cultures (denoted by "c" in sample name)
  seqtab_filtered <- seqtab_filtered[-c(grep("c", rownames(seqtab_filtered))), , drop=F]
}

# remove empty ASVs
print("Removing empty ASVs...")

# subset by non-empty ASVs
# sum abundances by ASV
seqtab_filtered <- rbind(seqtab_filtered, apply(seqtab_filtered, 2, sum))
# subset by non-empty ASVs
seqtab_filtered <- seqtab_filtered[,which(seqtab_filtered[nrow(seqtab_filtered),] != 0)]
# remove sum row
seqtab_filtered <- seqtab_filtered[c(1:(nrow(seqtab_filtered)-1)),]

# subset taxa table to include only the non-empty ASV assignments
tax <- tax[which(rownames(tax) %in% names(as.data.frame(seqtab_filtered))),]


################################ create phyloseq object ################################


print("Saving to phyloseq object...")

# transform to matrix
seqtab_mat <- as.matrix(seqtab_filtered)

tax_mat <- as.matrix(tax) # required to create phyloseq object

# create list of samples with no MiSeq code
rm_id <- c()
for (id in 1:nrow(metadata)){
  if (metadata$MiSeqCode[id] == "") {
    rm_id <- c(rm_id, id)
  }
}
# remove samples with no MiSeq code
if (!is.null(rm_id)) {
  metadata <- metadata[-c(rm_id),]
}

# set rownames of metadata to be MiSeq Code
rownames(metadata) <- metadata$MiSeqCode

# use the below command to solve the error:
# "Error in validObject(.Object) : invalid class “phyloseq” object: 
# Component sample names do not match.
# Try sample_names()"
# filling in the required substitution
#rownames(seqtab_mat) <- gsub("X", "", rownames(seqtab_mat))

dada2 <- phyloseq(tax_table(tax_mat), 
                  otu_table(seqtab_mat, taxa_are_rows = FALSE), 
                  sample_data(metadata))

# add shannon alpha diversity column to metadata
shannon_dada2 <- estimate_richness(dada2, split = TRUE, measures = "Shannon")
sample_data(dada2)$Alpha_Shannon <- as.numeric(shannon_dada2[,1])

# save with raw sequences for phylogenetic tree script
saveRDS(dada2,paste0(path_out,"physeqob_DADA2_rawseqs.rds"))

# save sequences to refseqs slot and rename ASVs for convenience
asv_seqs <- Biostrings::DNAStringSet(taxa_names(dada2))
names(asv_seqs) <- taxa_names(dada2)
dada2 <- merge_phyloseq(dada2, asv_seqs)
taxa_names(dada2) <- paste0("ASV", seq(ntaxa(dada2)))

# print out phyloseq object to screen
dada2

# save phyloseq object to be imported into analysis scripts
saveRDS(dada2,paste0(path_out,"physeqob_DADA2.rds"))

print("Script completed")


## end of script

