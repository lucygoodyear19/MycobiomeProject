######################################################################
########## Alternative pipeline for processing sequences #############
######################################################################


# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# clear workspace
rm(list=ls())


############ initial set up, library and data loading ################


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
       1) path to data folder 
       2) path to folder that will contain filtered data files 
       3) path to save outputs to 
       4) forward primer 
       5) reverse primer
       6) base name of fastq files 
       7) path to cutadapt installation
       8) path to local version of unite database")
}
# load arguments into script
source(args)
# print arguments as check
print(path) # path to data folder
print(path2) # path to folder that will contain filtered data files
print(path_out) # path to save outputs to
print(FWD) # forward primer
print(REV) # reverse primer
print(base_prefix) # sequences files base name
print(cutadapt) # path to cutadapt installation
print(unite.ref) # path to unite database

# check cutadapt is installed
system2(cutadapt, args = "--version")

# check folder has correct contents
list.files(path)

# save all forward read path names as a list
fnFs <- sort(list.files(path, pattern = "R1_001.fastq$", full.names = TRUE))
# save all reverse read path names as a list
fnRs <- sort(list.files(path, pattern = "R2_001.fastq$", full.names = TRUE))

# check paths are correct
fnFs


######################### remove primers ###########################


# define function to produce all orientations of the input sequence
allOrients <- function(primer) {
  require(Biostrings) # in case function is called separately to main
  dna <- DNAString(primer)  # Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString)) # convert back to character vector
}

# run function on primers
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# check forward primer orientations
FWD.orients

# put N-filtered files in filtN_lg/ subdirectory for both forward and reverse reads
# and save in new lists
fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

# use dada2 filter and trim function on forward and reverse reads
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
# note maxN = 0 means any reads with one or more "N"s (representing an
# ambiguous base) will be removed

# reset working directory to the filtN_lg subdirectory
setwd(path2)

# rename files to .gz for use in next steps
list2 <- list.files(path2)
file.rename(list2, paste0(list2,".gz"))

# check folder has correct contents (including file names)
list.files(path2)

# save all renamed forward read path names as a new list
fnFs.filtN <- sort(list.files(path2, pattern = "R1_001.fastq.gz", full.names = TRUE))
# save all renamed forward read path names as a new list
fnRs.filtN <- sort(list.files(path2, pattern = "R2_001.fastq.gz", full.names = TRUE))

# define function to count number of reads in which the primer is found
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# only the first sample needs to be processed since 
# primers will be in the same place/orientation for all reads
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# flip reverse primer if it is in its reverse complement
REV <- REV.orients[["RevComp"]]

# create subdirectory called 'cutadapt' if it doesn't already exist
# and copy, paste and rename read files to the subdirectory
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# generate reverse complements of all forward and reverse reads
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# run Cutadapt
# note that the default in cutadapt is for only one primer sequence to
# be removed from each read, so to enable two (i.e. forward and reverse) 
# primer sequences to be removed, the "-n" flag (for number of times) 
# and "2" (for two times) need to be added as arguments
for (i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i], # input files
			     "--minimum-length 1"))
}

# test 1st sample to check it has worked (output should all be 0 now)
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# save all primer-free forward read path and file names in a list
cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq$", full.names = TRUE))
# add prefix to all file names
XcutFs <- sub(base_prefix,"",cutFs)

# save all primer-free reverse read path and file names in a list
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq$", full.names = TRUE))
# add prefix to all file names
XcutRs <- sub(base_prefix,"",cutRs)

# extract sample names
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(XcutFs, get.sample.name))
head(sample.names)

## From now on only forward reads will be processed because the reverse reads 
## are always of lower quality 
## (https://www.sciencedirect.com/science/article/pii/S1754504818302800)

# plot quality profiles for forward reads
pdf(file = paste(path_out,"Quality_Profiles.pdf",sep = ''), paper = 'A4')
plotQualityProfile(cutFs[1:2]) # HPC will save this in a pdf with the error plot
dev.off()


################# filter and trim primer-free reads ##################


# assign filenames for the outputs of the filtered reads
filtFs <- file.path(path.cut, "filtered", basename(cutFs))

# filter reads by:
# - maximum number of "expected errors" allowed in a read 
# (rather than averaged quality scores)
# - enforce minimum length of 50bp to removespurious very low length sequences
# - truncate lengths to 200bp
# - removes PhiX, which is a common control used by illumina
out <- filterAndTrim(cutFs, filtFs, maxN = 0, maxEE = (2), truncLen = (200),
                     minLen = 50, rm.phix = TRUE, compress = FALSE, multithread = TRUE)

# view data
head(out)
# save filtered data to a csv
write.csv(out,file=paste(path_out,"filtering_output.csv",sep = ''))
# view csv
out

# learn the error rates (rate of error for each possible transition (A→C, A→G, …)
errF <- learnErrors(filtFs, multithread = TRUE,verbose=TRUE)
# visualise estimated error rates
pdf(file = paste(path_out,"Error_Rates.pdf",sep = ''), paper = 'A4')
plotErrors(errF, nominalQ = TRUE) # HPC will save this in a pdf with the quality plot 
dev.off()

# dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
# name the derep-class objects by the sample names
names(derepFs) <- sample.names

# apply the core sample inference algorithm to the dereplicated data
dadaFs <- dada(derepFs, err = errF, multithread = TRUE) 

# construct an amplicon sequence variant table (ASV) table
# (this is a higher-resolution version of the OTU table produced by traditional methods.)
seqtab <- makeSequenceTable(dadaFs)
# view dimensions of the ASV table
dim(seqtab)

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# inspect distribution of sequence lengths
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
head(track)
track
# write dataframe to a csv
write.csv(track,file=paste(path_out,"read_counts_during_pipeline_steps.csv",sep = ''))


####################### assign taxonomy ##########################


# run assignTaxonomy function
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

# inspect taxonomic assignments
taxa.print <- taxa
# removing sequence rownames for display only
rownames(taxa.print) <- NULL
# view data
head(taxa.print)
taxa.print

# convert to a dataframe
taxa.print<-as.data.frame(taxa.print)

# print summary of dataframe to screen
summary(taxa.print)

# extract mock 1 from ASV table
if ("mock1" %in% colnames(seqtab.nochim)){
        unqs.mock <- seqtab.nochim["mock1",]
        # drop ASVs absent in the Mock
        unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE)
        # print to screen the number of sample sequences present in the mock community
        cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
}

# extract mock 2 from ASV table
if ("mock2" %in% colnames(seqtab.nochim)){
        unqs.mock <- seqtab.nochim["mock2",]
        # drop ASVs absent in the Mock
        unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE)
        # print to screen the number of sample sequences present in the mock community
        cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
}

# extract mock 3 from ASV table
if ("mock3" %in% colnames(seqtab.nochim)){
        unqs.mock <- seqtab.nochim["mock3",]
        # drop ASVs absent in the Mock
        unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE)
        # print to screen the number of sample sequences present in the mock community
        cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
}

# extract mock 4 from ASV table, if it exists
if ("mock4" %in% colnames(seqtab.nochim)){
	unqs.mock <- seqtab.nochim["mock4",]
	# drop ASVs absent in the Mock
	unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE)
	# print to screen the number of sample sequences present in the mock community
	cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
}


####################### rename any dulicates ########################


# extract sample names from final dataset and save as dataframe
samples.out <- data.frame(Sample=as.character(row.names(seqtab.nochim)))

# rename any duplicates
samples.out$Sample <- with(samples.out, make.unique(as.character(Sample)))

# reset row names in seqtab.nochim to remove duplicates for later imports
row.names(seqtab.nochim) <- samples.out$Sample


######################### write files out ##########################


# taxa
# remove names
taxa.slim <- unname(taxa)
# rename rows
rownames(taxa.slim) <- paste0("SV",seq(nrow(taxa)))  

# abundances
seqtab.slim <- t(seqtab.nochim) # transpose
# rename rows
rownames(seqtab.slim) <- paste0("SV",seq(nrow(taxa)))

# save as .txt files
write.table(taxa.slim,paste(path_out,"Taxa_Table.txt",sep = ''))
write.table(seqtab.slim,paste(path_out,"Abundance_Table.txt",sep = ''))
write.table(seqtab.nochim,paste(path_out,"Seq_Abun_Table.txt",sep=''))


## end of script
