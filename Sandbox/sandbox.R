
samples.out <- data.frame(Sample=as.character(c(167,168,168,169,170,"NC","NC","NC")), stringsAsFactors = FALSE)
# for loop to rename any duplicates
for (i in 1:nrow(samples.out)) {
  if (duplicated(samples.out)[i] == TRUE){
    samples.out[i,1] = paste0(samples.out[i,1], ".1")
  }
}

samples.out$Sample <- with(samples.out, make.unique(as.character(Sample)))
samples.out


# reset row names for comparison
rownames(samples.out) <- rownames(seqtab.nochim)

samples.out <- data.frame(Sample=as.character(c(167,168,"168c",169,"170c","NC","NC","NC")), stringsAsFactors = FALSE)
seqtab.nochim <- data.frame(Stuff=as.character(c("ACT","AGCT","GTC","A","T","G","C","ATG")), stringsAsFactors = FALSE)
row.names(seqtab.nochim) <- c(9,10,11,12,13,14,15,16)
samples.out
seqtab.nochim
# for loop to rename any duplicates
samples.out$Sample <- with(samples.out, make.unique(as.character(Sample)))
samples.out
# reset row names for comparison
row.names(samples.out) <- rownames(seqtab.nochim)
samples.out
row.names(seqtab.nochim) <- samples.out$Sample
seqtab.nochim

# remove cultures (denoted by "c" in sample name)
seqtab <- seqtab.nochim[-c(grep("c", rownames(seqtab.nochim))), , drop=F]
seqtab

