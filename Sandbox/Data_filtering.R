##############################################################################
###################### Data filtering for DADA2 results ######################
##############################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear worksapce
rm(list=ls())


################################## Set up ####################################


# import arguments to run script on specific country data
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) # setup to accept arguments from command line
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied in the form of an R-script containing the following arguments: 
       1) path containing results from DADA2 pipeline")
}
# load arguments into script
source(args)
# print arguments as check
print(path_out)



# assign run/country to a variable
run <- "Taiwan_Vietnam_2016/"
country <- "Taiwan" # no forward slash after country
# assign root paths to a variable
root_path <- "/Users/lucy/Documents//MRes/MycobiomeProject/Analysis/Runs_Countries/"

# assign metadata path to a variable
metadata_path <- paste0(root_path, run, "metadata.csv")

# assign full paths to arguments for import into R scripts
path_in <- paste0(root_path, run, "DADA2_Results/", country, "/phil_truncl/phil_truncl/")
path_out <- paste0(root_path, run, "DADA2_Results/", country, "/")


# load taxa data
tax <- read.table(paste0(path_in, "Taxa_Table.txt"), 
                  stringsAsFactors = F)
# rename columns
colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# load abundance data
seqtab <- read.table(paste0(path_in, "Abun_Table_sample_row.txt"), 
                     stringsAsFactors = F, 
                     header = T)

# load plate data
plates <- read.csv(paste0(root_path, run, "Plate_Data.csv"), 
                   stringsAsFactor=F, 
                   header=T)
# name plates columns
#names(plates) <- c("barcode","sample","DNAqual","plate")

# load metadata
metadata <- read.csv(paste0(root_path, run,"metadata.csv"),
                     stringsAsFactor=F, 
                     header=T)
# rename rownames to sample names
rownames(metadata) <- metadata$MiSeqCode
# remove defunct column
metadata <- metadata[,-1]
# add column for Bd +ve/-ve to metadata
metadata$Bd <- 0
for (i in (1:nrow(metadata))){
  if (metadata$Bd_GE[i] > 0){
    metadata$Bd[i] <- 1
  }
}
# set as factor for plotting
metadata$Bd <- as.factor(metadata$Bd)


##################### label abundance data with sample names ####################


# change row names to match plates format
rownames(seqtab) <- gsub("-", "", rownames(seqtab))

# merge plates and seqtab by sample name
seqtab_names <- merge(plates, seqtab, by.x="indices", by.y = "row.names")
# rename row names to sample names
rownames(seqtab_names) <- seqtab_names$Sample_Name
# remove defunct columns
seqtab_names <- seqtab_names[-c(1:9)]


######################### remove PosC, Mocks and NC ###########################


#################### Set up


# set rownames of taxonomy table to match column names of abundance table
rownames(tax) <- colnames(seqtab_names)

# transpose abundance table
seqtab_names_t <- as.data.frame(t(seqtab_names))
# subset taxonomy table
species_by_seq <- as.data.frame(cbind(rownames(tax), tax$Genus, tax$Species))
names(species_by_seq) <- c("Sequence", "Genus", "Species")

# merge subsetted taxa table and the abundance table
species_abun <- merge(species_by_seq, seqtab_names_t, by.x = "Sequence", by.y = "row.names")


######################### PosC


# subset merged data frame to include only positive controls
posc <- cbind(species_abun$Sequence, paste(species_abun$Genus, species_abun$Species),species_abun[grep("posC",names(species_abun))])
# rename columns
names(posc)[1:2] <- c("Sequence", "Genus_Species")

# find posC species
poscspec1 <- as.character(posc$Genus_Species[posc$posC1 > 1000])
poscspec2 <- as.character(posc$Genus_Species[posc$posC2 > 1000])
poscspec3 <- as.character(posc$Genus_Species[posc$posC3 > 1000])
poscspec4 <- as.character(posc$Genus_Species[posc$posC4 > 1000])

# find unique posC species
uniqposcspec <- unique(c(poscspec1, poscspec2, poscspec3, poscspec4))
print(paste("PosC species is", uniqposcspec))

# subset sequences that are the same as the species of posC and count them
posc_seqs <- as.character(posc$Sequence[posc$Genus_Species %in% unique(uniqposcspec)])
length(posc_seqs)

# remove PosC sequences
seqtab_names <- seqtab_names[,!names(seqtab_names) %in% posc_seqs]


###################### Mocks


# subset merged data frame to include only positive controls
mocks <- cbind(species_abun$Sequence, paste(species_abun$Genus, species_abun$Species),species_abun[grep("mock",names(species_abun))])
# rename columns
names(mocks)[1:2] <- c("Sequence", "Genus_Species")

# find mock species
mockspec1 <- as.character(mocks$Genus_Species[mocks$mock1 > 1000])
mockspec2 <- as.character(mocks$Genus_Species[mocks$mock2 > 1000])
mockspec3 <- as.character(mocks$Genus_Species[mocks$mock3 > 1000])
mockspec4 <- as.character(mocks$Genus_Species[mocks$mock4 > 1000])

# find unique mock species
uniqmockspec <- unique(c(mockspec1, mockspec2, mockspec3, mockspec4))
print(paste("Mock species are", uniqmockspec))

# subset sequences that are the same as the species of the mocks and count them
mock_seqs <- as.character(mocks$Sequence[mocks$Genus_Species %in% unique(uniqmockspec)])
length(mock_seqs)

# remove mock sequences
seqtab_names <- seqtab_names[,!names(seqtab_names) %in% mock_seqs]
# remove mocks from samples
seqtab_names <- seqtab_names[!rownames(seqtab_names) %in% c("mock1", "mock2", "mock3", "mock4"),]


##################### Controls


# define function to remove posC and NTC reads from samples by plate
# returns a dataframe containing samples from the speicifed plate, in the same form as df_abun 
filter_controls <- function(df_abun, df_plates, plate_no) {
  
  # subset by plate number
  plate <- df_abun[,names(df_abun) %in% c(df_plates$sample[df_plates$plate == paste0("plate", plate_no)])]
  # create new column summing NTC and PosC reads per sample
  plate$sum <- plate[[paste0("NTC", plate_no)]] + plate[[(paste0("posC", plate_no))]]
  
  # print summary statements
  print(paste0("Number of samples in plate: ", length(names(plate)))) # count samples
  print(paste0("Number of ASVs with NTC/PosC reads not equal to 0: ", length(which(plate$sum != 0)))) # count ASVs with NTC + posC reads not equal to 0
  print(paste0("Number of NTC and posC reads overall before: ", sum(plate$sum))) # sum overall NTC + posC reads
  
  # remove number of reads in NTC + posC ASVs from corresponding ASVs for each sample
  platex <- as.data.frame(t(apply(plate, 1, function(x) x - x["sum"])))
  # set any negative values (create by the above) to 0
  platex[platex < 0] <- 0
  
  # print a check statement
  print(paste0("CHECK. Number of NTC and posC reads overall after: ", sum(platex$sum), "(should be 0)")) # should be 0
  # remove sum column
  platex <- subset(platex, select=-sum)
  
  return(platex)
}

# transpose the filtered abundance table
seqtab_names_t <- as.data.frame(t(seqtab_names))
# run the function on all four plates
#plate1 <- filter_controls(seqtab_names_t, plates, "1")
plate2 <- filter_controls(seqtab_names_t, plates, "2")
plate3 <- filter_controls(seqtab_names_t, plates, "3")
plate4 <- filter_controls(seqtab_names_t, plates, "4")

# combine the plates above into one dataframe
seqtab_names_t_filtered <- cbind(plate2, plate3, plate4)

# remove empty all samples (including posC and NTC) from data frame
# sum abundances by sample
seqtab_names_t_filtered <- rbind(seqtab_names_t_filtered, apply(seqtab_names_t_filtered, 2, sum))
# transpose dataframe
seqtab_names_filtered <- t(seqtab_names_t_filtered)
# count how many samples have no reads
length(which(seqtab_names_filtered[,ncol(seqtab_names_filtered)] == 0))
# save empty samples as a vector
to_remove <- which(seqtab_names_filtered[,ncol(seqtab_names_filtered)] == 0)
# remove empty samples from dataframe
seqtab_names_filtered <- as.data.frame(seqtab_names_filtered[!rownames(seqtab_names_filtered) %in% names(to_remove),])
# remove sum column
seqtab_names_filtered <- seqtab_names_filtered[,-ncol(seqtab_names_filtered)]


######################### Tidy up


########### for Taiwan 2016 only, remove all Bsal and cultures


library(tibble) # for `rownames_to_column` and `column_to_rownames`
library(dplyr)

# separate out Bsal from taxa table
bsal <- tax %>%
  rownames_to_column('Seq') %>%
  filter(Species == "s__salamandrivorans") %>%
  column_to_rownames('Seq')

# save Bsal sequences
bsal_seqs <- rownames(bsal)

# remove all Bsal sequences from abundance table
seqtab_names_filtered <- as.data.frame(seqtab_names_filtered[,!colnames(seqtab_names_filtered) %in% bsal_seqs,])

# remove cultures (denoted by "c" in sample name)
seqtab_names_filtered <- seqtab_names_filtered[-c(grep("c", rownames(seqtab_names_filtered))), , drop=F]


################# For all


# remove empty ASVs

# subset by non-empty ASVs
# sum abundances by ASV
seqtab_names_filtered <- rbind(seqtab_names_filtered, apply(seqtab_names_filtered, 2, sum))
# subset by non-empty ASVs
seqtab_names_filtered <- seqtab_names_filtered[,which(seqtab_names_filtered[nrow(seqtab_names_filtered),] != 0)]
# remove sum row
seqtab_names_filtered <- seqtab_names_filtered[c(1:(nrow(seqtab_names_filtered)-1)),]

# subset taxa table to include only the non-empty ASV assignments
tax <- tax[which(rownames(tax) %in% names(seqtab_names_filtered)),]


################################ create phyloseq object ################################


print("Saving to phyloseq object...")

# transform to matrix
seqtab_mat <- as.matrix(seqtab_names_filtered)

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

require(phyloseq)
dada2 <- phyloseq(tax_table(tax_mat), 
                  otu_table(seqtab_mat, taxa_are_rows = FALSE), 
                  sample_data(metadata))

# save with raw sequences for phylogenetic tree script
saveRDS(dada2,paste0(path_out,"physeqob_DADA2_tree.rds"))

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
