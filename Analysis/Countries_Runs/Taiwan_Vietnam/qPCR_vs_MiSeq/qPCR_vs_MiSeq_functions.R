############################################################################
################# Functions for qPCR vs MiSeq analysis #####################
############################################################################

# Author: Lucy Goodyear (lucy.goodyear19@imperial.ac.uk)
# Version: 0.0.1

# total exact matches between miseq and qpcr
total_match <- function(df){
  totalmatch <- 0
  for (i in 1:nrow(df)) {
    if (((df$miseq[i] != 0) & (df$qpcr[i] == 1)) | ((df$miseq[i] == 0) & (df$qpcr[i] == 0))) {
      totalmatch <- totalmatch + 1
    }
  }
  return(totalmatch/nrow(df))
}

# total matches when no. reads <= 100 is assumed 0
total_match100 <- function(df) {
  totalmatch100 <- 0
  for (i in 1:nrow(df)) {
    if (((df$miseq[i] >= 100) & (df$qpcr[i] == 1)) | ((df$miseq[i] < 100) & (df$qpcr[i] == 0))) {
      totalmatch100 <- totalmatch100 + 1
    }
  }
  return(totalmatch100/nrow(df))
}

# total 0-0 matches
qpcr_miseq_zero <- function(df) {
  match0 <- 0
  for (i in 1:nrow(df)) {
    if ((df$miseq[i] == 0) & (df$qpcr[i] == 0)) {
      match0 <- match0 + 1
    }
  }
  return(match0/nrow(df))
}

# proportion of samples where qpcr is negative but miseq is positive
qpcr_neg_miseq_pos <- function(df) {
  nomatch0 <- 0
  for (i in 1:nrow(df)) {
    if ((df$miseq[i] != 0) & (df$qpcr[i] == 0)) {
      nomatch0 <- nomatch0 + 1
    }
  }
  return(nomatch0)
}

# proportion of samples where qpcr is negative but miseq 
# is positive >= 100 reads
qpcr_neg_miseq_pos100 <- function(df) {
  nomatch100 <- 0
  for (i in 1:nrow(df)) {
    if ((df$miseq[i] >= 100) & (df$qpcr[i] == 0)) {
      nomatch100 <- nomatch100 + 1
    }
  }
  return(nomatch100)
}