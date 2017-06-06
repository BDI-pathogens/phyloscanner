# Blacklister for exact duplicates. The entries argument is a list, the entries of each being groups of tips whose reads are identical.

blacklist.exact.duplicates <- function(entries, raw.threshold, ratio.threshold, regexp, verbose = F){

  pairs.table	<- data.frame(V1=character(0), V2=character(0))
  
  for(entry in entries){
    # get rid of anything that doesn't match the regexp (outgroup etc)
    tmp <-  entry[!is.na(sapply(entry, read.count.from.label, regexp = regexp))]
    if(length(tmp)>1){
      tmp <- t(combn(tmp,2))
      colnames(tmp) <- c('V1','V2')
      pairs.table <- rbind(pairs.table, tmp)  
    }  
  }
  
  pairs.table$V1 <- as.character(pairs.table$V1)
  pairs.table$V2 <- as.character(pairs.table$V2)
  
  #	OR suggest: make get.count return integer(0) to avoid issues when V1 is character(0)
  tmp <- unlist(sapply(pairs.table$V1, function(x) read.count.from.label(x, regexp)))
  if(!is.null(tmp))
    pairs.table$reads.1 <- tmp
  if(is.null(tmp))
    pairs.table$reads.1 <- integer(0)
  
  tmp <- unlist(sapply(pairs.table$V2, function(x) read.count.from.label(x, regexp)))
  if(!is.null(tmp))
    pairs.table$reads.2 <- tmp
  if(is.null(tmp))
    pairs.table$reads.2 <- integer(0)
  
  # reverse the order so the read with the greater count is in the first column
  pairs.table[pairs.table$reads.2 > pairs.table$reads.1, c("V1", "V2", "reads.1", "reads.2")] <- pairs.table[pairs.table$reads.2 > pairs.table$reads.1, c("V2", "V1", "reads.2", "reads.1")] 
  
  # calculate the counts
  pairs.table$ratio <- pairs.table$reads.2 / pairs.table$reads.1
  if (verbose) cat("Making blacklist with a ratio threshold of ",ratio.threshold," and a raw threshold of ",raw.threshold,"\n",sep="")
  blacklisted <- pairs.table[which(pairs.table$ratio<ratio.threshold | (pairs.table$reads.2<raw.threshold)),2]
  
  return(blacklisted)
}