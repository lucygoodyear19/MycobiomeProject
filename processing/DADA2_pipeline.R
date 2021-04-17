###########################################################################
################## DADA2 pipeline for processing sequences ################
###########################################################################


# Author: Luke Goodyear (leg19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


############################ Initial set up ##############################


# load packages
library("dada2")
library("ShortRead")
library("Biostrings")
# personal laptop
#library("cutadapt")

# import arguments to run script on specific country data
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) # setup to accept arguments from command line
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied in the form of an R-script containing the following arguments: 
       1) path = path to data folder 
       2) path2 = path to folder that will contain filtered data files 
       3) path_out = path to save outputs to 
       4) FWD = forward primer 
       5) RWD = reverse primer
       6) base_prefix = base name of fastq files 
       7) cutadapt = path to cutadapt installation
       8) unite.ref = path to local version of unite database")
}
# load arguments into script
source(args)
# print arguments as check
print("Arguments:")
print(paste0("path = ", path)) # path to data folder
print(paste0("path2 = ", path2)) # path to folder that will contain filtered data files
print(paste0("path_out = ", path_out)) # path to save outputs to
print(paste0("FWD = ", FWD)) # forward primer
print(paste0("REV = ", REV)) # reverse primer
print(paste0("base_prefix = ", base_prefix)) # sequences files base name
print(paste0("cutadapt = ", cutadapt)) # path to cutadapt installation
print(paste0("unite.ref = ", unite.ref)) # path to unite database

# check cutadapt is installed
print("Cutadapt version:")
system2(cutadapt, args = "--version")

# check folder has correct contents by viewing first item
print("Print first item file name in folder containing the data:")
list.files(path)[1]

# save all forward read path names as a list
fnFs <- sort(list.files(path, pattern = patFs, full.names = TRUE))
# save all reverse read path names as a list
fnRs <- sort(list.files(path, pattern = patRs, full.names = TRUE))

# check paths are correct by viewing first item
print("Print path to first item in folder containing the data:")
fnFs[1]


############################ Identify primers #################################


# define function to produce all orientations of the input sequence
allOrients <- function(primer) {
  dna <- DNAString(primer)  # Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString)) # convert back to character vector
}

# run function on primers
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# check forward primer orientations
print("Forward primer orientations:")
FWD.orients

## remove sequences with ambiguous bases (have "N"s present)

# empty filtN directory of any files from previsou runs
unlink(paste0(path2, "*"))
# put N-filtered files in filtN/ subdirectory for both forward and reverse reads
# and save in new lists
fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

# use dada2 filter and trim function on forward and reverse reads
print("Filtering to remove sequences with ambiguous bases")
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE)
# note maxN = 0 means any reads with one or more "N"s (representing an
# ambiguous base) will be removed

# reset working directory to the filtN subdirectory
setwd(path2)

# rename files to .gz for use in next steps
list2 <- list.files(path2)
if (grepl(".gz$", list2[1]) == FALSE) {
  file.rename(list2, paste0(list2,".gz"))
}

# check folder has correct contents (including file names)
print("Print name of first item in filtN subdirectory:")
list.files(path2)[1]

# save all renamed forward read path names as a new list
fnFs.filtN <- sort(list.files(path2, pattern = "R1_001.fastq.gz", full.names = TRUE))
# save all renamed forward read path names as a new list
fnRs.filtN <- sort(list.files(path2, pattern = "R2_001.fastq.gz", full.names = TRUE))

# define function to count number of reads in which the primer is found
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# check and adjust orientation of primers to ensure orientation matches sample reads
for (i in 1:length(fnFs.filtN)) {
  # view orientation of primers for sample i (wil be the same for all reads in all samples)
  primer_check <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[i]]), 
                        FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[i]]), 
                        REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[i]]), 
                        REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[i]]))
  print(paste0("Primer count for sample ", i))
  print(primer_check)
  # define break flag
  stopFWD = FALSE
  stopREV = FALSE
  # check orientations of forward primer and adjust so orientation in sample is set as FWD
  if (primer_check[,"Forward"][1] > 0) {
    stopFWD = TRUE
  } else {
    for (j in c("Complement", "Reverse", "RevComp")) {
      if (primer_check[,j][1] > 0) {
        print(paste0("Replacing FWD with ", j, " orientation"))
        FWD <- FWD.orients[[j]]
        stopFWD = TRUE
        break # break out of FWD for loop
      }
    }
  }
  if (stopFWD == FALSE & i < length(fnFs.filtN)) {
    print("All zeros in FWD.ForwardReads, use next sample to find orientation of forward primer")
  }
  
  # check orientations of reverse primer and adjust so orientation in sample is set as REV
  if (primer_check[,"Forward"][4] > 0) {
    stopREV = TRUE
  } else {
    for (j in c("Complement", "Reverse", "RevComp")) {
      if (primer_check[,j][4] > 0) {
        print(paste0("Replacing REV with ", j, " orientation"))
        REV <- REV.orients[[j]]
        stopREV = TRUE
        break # break out of REV loop
      }
    }
  }
  if (stopREV == FALSE & i < length(fnFs.filtN)) {
  print("All zeros in REV.ReverseReads, use next sample to find orientation of reverse primer")
  }
  # if sample orientations of both FWD and REV primers have been found, end loop
  # otherwise, continue to next sample
  if (stopFWD == TRUE & stopREV == TRUE) {
    break
  }
}
if (stopFWD == FALSE & stopREV == FALSE) {
  # error occurs if no primers are detected in FWD.ForwardReads or REV.ReverseReads in all samples
  stop("No primers detected in samples")
}


############################ Remove primers ################################


# create subdirectory called 'cutadapt' if it doesn't already exist
# and copy, paste and rename read files to the subdirectory
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# generate reverse complements of all forward and reverse reads
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# create flags to trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# create flags to trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# run Cutadapt
# note that the default in cutadapt is for only one primer sequence to
# be removed from each read, so to enable two (i.e. forward and reverse) 
# primer sequences to be removed, the "-n" flag (for number of times) 
# and "2" (for two times) need to be added as arguments
print("Begin cutadapt output:")
for (i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i], # input files
			     "--minimum-length 1"))
}
print("End of cutadapt output")

# test 1st sample to check all primers have been removed (output should all be 0 now)
print("Check all primers have been removed:")
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# save all primer-free forward read path and file names in a list
cutFs <- sort(list.files(path.cut, pattern = patFs, full.names = TRUE))
# add prefix to all file names
XcutFs <- sub(base_prefix,"",cutFs)

# save all primer-free reverse read path and file names in a list
cutRs <- sort(list.files(path.cut, pattern = patRs, full.names = TRUE))
# add prefix to all file names
XcutRs <- sub(base_prefix,"",cutRs)

# extract sample names
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(XcutFs, get.sample.name))
print("View first few sample names:")
head(sample.names)

print("From now on only forward reads will be processed because the reverse reads are always of lower quality") 
## (https://www.sciencedirect.com/science/article/pii/S1754504818302800)


########################### Inspect read quality profiles ############################


# plot quality profiles for forward reads
print("Plotting quality profiles for forward reads")
pdf(file = paste0(path_out,"Quality_Profiles.pdf"), paper = 'A4')
for (i in 1:length(cutFs)){
	print(plotQualityProfile(cutFs[i]))
}
dev.off()


####################### Filter and trim primer-free reads #######################


# assign filenames for the outputs of the filtered reads
filtFs <- file.path(path.cut, "filtered", basename(cutFs))

# filter reads by:
# - maximum number of "expected errors" allowed in a read 
# (rather than averaged quality scores)
# - enforce minimum length of 50bp to removespurious very low length sequences
# - truncate reads at the first instance of a quality score of less than 2 (default)
# - removes PhiX, which is a common control used by illumina (default)
# note does not compress the output files
print("Filtering reads")
out <- filterAndTrim(cutFs, filtFs, maxEE = 2, minLen = 50, compress = FALSE)

# view data
print("View number of reads before and after filtering/trimming per sample:")
out

########################### Learn the error rates ####################################


# learn the error rates (rate of error for each possible transition (A->C, A->G, …)
print("Learning error rates")
errF <- learnErrors(filtFs)
print("Visualising estimated error rates by plotting to pdf")
pdf(file = paste(path_out,"Error_Rates.pdf",sep = ''), paper = 'A4')
plotErrors(errF, nominalQ = TRUE) 
dev.off()


############################# Dereplicate identical reads ###########################


# dereplicate identical reads
print("Dereplicating identical reads")
derepFs <- derepFastq(filtFs)
# name the derep-class objects by the sample names
names(derepFs) <- sample.names


################################# Sample Inference ##################################


# apply the core sample inference algorithm to the dereplicated data
print("Running sample inference algorithm")
dadaFs <- dada(derepFs, err = errF) 


##################### Constructing ASV and taxonomy tables ########################


# construct an amplicon sequence variant table (ASV) table
# (this is a higher-resolution version of the OTU table produced by traditional methods.)
print("Constructing ASV table")
seqtab <- makeSequenceTable(dadaFs)
# view dimensions of the ASV table
print("Dimensions of ASV table:")
dim(seqtab)

# remove chimeras
print("Removing chimeras")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)

# inspect distribution of sequence lengths
print("Distribution of sequence lengths:")
table(nchar(getSequences(seqtab.nochim)))

# inspect number of reads that made it through each step of the pipeline to 
# verify everything worked as expected (i.e. not too many reads were lost after filtering)

# set a function to find the sum of uniques
getN <- function(x) sum(getUniques(x))

# create dataframe with total numbers of reads at each stage for each sample
track <- cbind(out, 
               sapply(dadaFs, getN), 
               rowSums(seqtab.nochim))
# note, if processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)

# name columns in dataframe
colnames(track) <- c("input", "filtered", "denoisedF", 
                     "nonchim")
# set sample names as row names in data frame
rownames(track) <- sample.names

# view dataframe
print("Number of reads that made it through each step of the pipeline:")
track
# write dataframe to a csv
write.csv(track,file=paste(path_out,"read_counts_during_pipeline_steps.csv",sep = ''))


####################### assign taxonomy ##########################


# run assignTaxonomy function
print("Assigning taxonomy")
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = FALSE, tryRC = TRUE)

# inspect taxonomic assignments
taxa.print <- taxa
# removing sequence rownames for display only
rownames(taxa.print) <- NULL
# view data
print("View first few taxonomic assignments:")
head(taxa.print)

# convert to a dataframe
taxa.print<-as.data.frame(taxa.print)

# print summary of dataframe to screen
print("Summary of taxonomic assignments:")
summary(taxa.print)

# remove names
taxa.slim <- unname(taxa)


####################### Rename any dulicates ########################


print("Renaming any duplicate sample names")
# extract sample names from final dataset and save as dataframe
samples.out <- data.frame(Sample=as.character(row.names(seqtab.nochim)))

# rename any duplicates
samples.out$Sample <- with(samples.out, make.unique(as.character(Sample)))

# reset row names in seqtab.nochim to remove duplicates for later imports
row.names(seqtab.nochim) <- samples.out$Sample


####################### View mock communities ########################


seqtab.slim <- t(seqtab.nochim) # transpose

print("View mock community assignments:")

# extract mock 1 from ASV table
if ("mock1" %in% colnames(seqtab.slim)){
        unqs.mock <- seqtab.nochim["mock1",]
        # drop ASVs absent in the Mock
        unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE)
        # print to screen the number of sample sequences present in the mock community
        cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock1 community.\n")
}

# extract mock 2 from ASV table
if ("mock2" %in% colnames(seqtab.slim)){
        unqs.mock <- seqtab.nochim["mock2",]
        # drop ASVs absent in the Mock
        unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE)
        # print to screen the number of sample sequences present in the mock community
        cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock2 community.\n")
}

# extract mock 3 from ASV table
if ("mock3" %in% colnames(seqtab.slim)){
        unqs.mock <- seqtab.nochim["mock3",]
        # drop ASVs absent in the Mock
        unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE)
        # print to screen the number of sample sequences present in the mock community
        cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock3 community.\n")
}

# extract mock 4 from ASV table, if it exists
if ("mock4" %in% colnames(seqtab.slim)){
        unqs.mock <- seqtab.nochim["mock4",]
        # drop ASVs absent in the Mock
        unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE)
        # print to screen the number of sample sequences present in the mock community
        cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock4 community.\n")
}


######################### Write files out ##########################


# save as .txt files
print("Saving final outputs")
write.table(taxa.slim,paste(path_out,"Taxa_Table.txt",sep = ''))
write.table(seqtab.nochim,paste(path_out,"Abun_Table_sample_row.txt",sep=''))

print("DADA2 pipeline complete")


## end of script
