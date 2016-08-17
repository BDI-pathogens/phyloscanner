args = commandArgs(trailingOnly=TRUE)

get.count <- function(string){
  if(length(grep(regexp, string)>0)) {
    return(as.numeric(sub(regexp, "\\3", string)))
  } else {
    return(NA)
  }
}

# This script blacklists tree tips that are likely to be contaminants, based on being identical
# to tips from another patient. Four arguments:

# 1. A raw read threshold. The script will blacklist any tips that match those from another
# patient if the read count for that tip is below this value.
# 2. A ratio read threshold. The script will blacklist any tips that match those from another
# patient if the ratio of tips from that patient to the tip's patient is below this value
# 3. An input file name (A CSV file, output from XXXXXXXX)
# 4. An output file name (CSV format)

# Parse arguments

raw.threshold <- as.numeric(args[1])
ratio.threshold <- as.numeric(args[2])
regexp <- args[3]
input.name <- args[4]
output.name <- args[5]




# Read the lines of the input file as a list

cat("Opening file: ", input.name, "\n", sep="")

entries <- strsplit(readLines(input.name, warn=F),",")

first <- T

# The input file has variable row lengths - sometimes three or more patients have
# identical reads. We build a data frame for each pairwise combination of reads in each
# row

for(entry in entries){
  # get rid of anything that doesn't match the regexp (outgroup etc)
  entry <- entry[!is.na(get.count(entry))]

  if(first){
    first <- F
    raw.tab <- t(combn(entry,2))
  } else {
    raw.tab <- rbind(raw.tab, t(combn(entry,2)))
  }
}


pairs.table <- data.frame(raw.tab, stringsAsFactors = F)

pairs.table$reads.1 <- sapply(pairs.table[,1], get.count)
pairs.table$reads.2 <- sapply(pairs.table[,2], get.count)


# reverse the order so the read with the greater count is in the first column

pairs.table[pairs.table$reads.2 > pairs.table$reads.1, c("X1", "X2", "reads.1", "reads.2")] <- pairs.table[pairs.table$reads.2 > pairs.table$reads.1, c("X2", "X1", "reads.2", "reads.1")] 

# calculate the counts

pairs.table$sharedCount <- pairs.table$reads.2 + pairs.table$reads.1
pairs.table$seqOverShared <- pairs.table$reads.2/pairs.table$sharedCount

cat("Making blacklist with a ratio threshold of ",ratio.threshold," and a raw threshold of ",raw.threshold,"\n",sep="")

blacklisted <- pairs.table[which(pairs.table$seqOverShared<ratio.threshold | (pairs.table$reads.2<raw.threshold)),2]

cat("Blacklist length: ",length(blacklisted),"\n",sep="")

if(length(blacklisted)>0){
  write.table(blacklisted,output.name, sep=",", row.names=FALSE, col.names=FALSE)
}
