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

#############################################################

Forward <- c(4214, 0, 0, 4200)
Complement <- c(0, 0, 0, 0)
Reverse <- c(0, 0, 0, 0)
RevComp <- c(0, 3743, 3590, 0)

primer_check <- cbind(Forward, Complement, Reverse, RevComp)
primer_check

# flip reverse primer if it is in its reverse complement
REV <- c()
if (primer_check[,"RevComp"][2] > 0){
  REV <- REV.orients[["RevComp"]]
}
# flip forward primer if it is in its complement
if (primer_check[["Complement"]] > 0 {
  FWD <- FWD.orients[["Complement"]]
}