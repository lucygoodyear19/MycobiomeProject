#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

source(args)

print(path)


###############################################################

seqtab <- rbind(c(2,3,4,1,2,3,7), c(8,5,4,3,2,1,6))
seqtab <- rbind(seqtab, apply(seqtab, 2, sum))

