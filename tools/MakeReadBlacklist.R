list.of.packages <- "argparse"
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")

library(argparse)

get.count <- function(string){
  if(length(grep(regexp, string)>0)) {
    return(as.numeric(sub(regexp, "\\3", string)))
  } else {
    return(NA)
  }
}

arg_parser = ArgumentParser(description="Identify phylogeny tips for blacklisting as representing suspected contaminants, based on being exact duplicates of more numerous reads from other patients")

arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", 
                        help="Regular expression identifying tips from the dataset. Three groups: patient ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name.")
arg_parser$add_argument("rawThreshold", action="store", type="double", help="Raw threshold; tips with read counts less than this that are identical to a tip from another patient will be blacklisted, regardless of the count of the other read.")
arg_parser$add_argument("ratioThreshold", action="store", type="double", help="Ratio threshold; tips will be blacklisted if the ratio of their tip count to that of of another, identical tip from another patient is less than this value.")
arg_parser$add_argument("inputFileName", action="store", help="A CSV file outlining groups of tips that have identical sequences, each forming a single line.")
arg_parser$add_argument("outputFileName", action="store", help="The file to write the output to, a list of tips to be blacklisted.")

# Parse arguments

args <- arg_parser$parse_args()

raw.threshold <- args$rawThreshold
ratio.threshold <- args$ratioThreshold
regexp <- args$tipRegex
input.name <- args$inputFileName
output.name <- args$outputFileName

# Read the lines of the input file as a list

cat("MakeReadBlacklist.R run on: ", input.name, "\n", sep="")
entries <- strsplit(readLines(input.name, warn=F),",")
# The input file has variable row lengths - sometimes three or more patients have
# identical reads. We build a data frame for each pairwise combination of reads in each
# row
pairs.table	<- data.frame(V1=character(0), V2=character(0))
for(entry in entries)
{
  # get rid of anything that doesn't match the regexp (outgroup etc)
  tmp <- t(combn(entry[!is.na(get.count(entry))],2))
  colnames(tmp) <- c('V1','V2')
  pairs.table <- rbind(pairs.table, tmp)
}
pairs.table$V1 <- as.character(pairs.table$V1)
pairs.table$V2 <- as.character(pairs.table$V2)
#	OR suggest: make get.count return integer(0) to avoid issues when V1 is character(0)
tmp <- unlist(sapply(pairs.table$V1, get.count))
if(!is.null(tmp))
	pairs.table$reads.1 <- tmp
if(is.null(tmp))
	pairs.table$reads.1 <- integer(0)
tmp <- unlist(sapply(pairs.table$V2, get.count))
if(!is.null(tmp))
	pairs.table$reads.2 <- tmp
if(is.null(tmp))
	pairs.table$reads.2 <- integer(0)
# reverse the order so the read with the greater count is in the first column
pairs.table[pairs.table$reads.2 > pairs.table$reads.1, c("V1", "V2", "reads.1", "reads.2")] <- pairs.table[pairs.table$reads.2 > pairs.table$reads.1, c("V2", "V1", "reads.2", "reads.1")] 
# calculate the counts
pairs.table$sharedCount <- pairs.table$reads.2 + pairs.table$reads.1
pairs.table$seqOverShared <- pairs.table$reads.2/pairs.table$sharedCount

cat("Making blacklist with a ratio threshold of ",ratio.threshold," and a raw threshold of ",raw.threshold,"\n",sep="")

blacklisted <- pairs.table[which(pairs.table$seqOverShared<ratio.threshold | (pairs.table$reads.2<raw.threshold)),2]

cat("Blacklist length: ",length(blacklisted),"\n",sep="")


write.table(blacklisted,output.name, sep=",", row.names=FALSE, col.names=FALSE, quote=F)

