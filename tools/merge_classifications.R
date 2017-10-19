#!/usr/bin/env Rscript

suppressMessages(require(data.table, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(phyloscannerR, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(kimisc, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))

arg_parser = ArgumentParser()

arg_parser$add_argument("-rda", "--asRDA", action="store_true", help="Write output as an RDA file rather than CSV.")
arg_parser$add_argument("-s", "--summaryFile", action="store", help="The full output file from summary_statistics.R; tip and read count information will be read from this if it is given.")
arg_parser$add_argument("-cfe", "--csvFileExtension", action="store", default="csv", help="The file extension for table files (default .csv).")
arg_parser$add_argument("inputFiles", action="store", help="Either a list of all input files (output from classify_relationships.R), separated by colons, or a single string that begins every input file name.")
arg_parser$add_argument("outputFile", action="store", help="A .csv file to write the output to.")
arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what I'm doing.")

args <- arg_parser$parse_args()

summary.file             <- args$summaryFile
output.file              <- args$outputFile
verbose                  <- args$verbose
csv.fe                   <- args$csvFileExtension
input.file.name          <- args$inputFiles

input.files <- list.files.mod(dirname(input.file.name), pattern=paste(basename(input.file.name)), full.names=TRUE)

if(length(input.files)==0){
  stop("No input files found.")
}

all.tree.info <- list()

for(file in input.files){
  tree.info <- list()
  tree.info$classification.file.name <- file
  
  suffix <- get.suffix(file, input.file.name, csv.fe)
  tree.info$suffix <- suffix
  
  all.tree.info[[file]] <- tree.info
}

reads.table	<- NULL		# can replace prev give.denom by is.null(reads.table)
if(!is.null(summary.file)){
  if (verbose) cat("Getting window counts per patient from ",summary.file,"...\n", sep="")
  reads.table 	<- read.csv(summary.file, stringsAsFactors = F)
  required.columns <- c("file.suffix", "id", "reads", "tips")
  missing.columns <-
    required.columns[! required.columns %in% colnames(reads.table)]
  if (length(missing.columns) > 0) {
    stop(paste0("The following required columns are missing from the summary",
                " csv file: ", paste(missing.columns, collapse=', '), ". Quitting."))
  }
  reads.table 	<- reads.table[,c("file.suffix", "id", "reads", "tips")]
  reads.table 	<- as.data.table(reads.table)
  setnames(reads.table, c("file.suffix"), c("SUFFIX"))
  reads.table$SUFFIX <- as.character(reads.table$SUFFIX)
  reads.table[, present:= reads>0]  
} 

tt <- merge.classifications(all.tree.info)

if(!is.null(summary.file)){
  if(verbose) cat("Combining window data...\n")
  
  dp	<-  reads.table[,{
    #	combine by window only for present patients to avoid large mem
    dp			<- data.table(HOST.1=NA_character_, HOST.2=NA_character_)
    if(length(which(present))>1){
      dp			<- t( combn( id[present], 2 ) )
      colnames(dp)<- c('HOST.1','HOST.2')
      dp			<- as.data.table(dp)	
      # make sure combinations are in the direction as in 'tt'
      tmp			<- copy(dp)
      setnames(tmp, c('HOST.1','HOST.2'), c('HOST.2','HOST.1'))
      dp			<- rbind(dp, tmp)				
      dp			<- subset(dp, HOST.1<HOST.2)
    }
    dp						
  }, by=c('SUFFIX')]
  dp <- subset(dp, !is.na(HOST.1) & !is.na(HOST.2))
  tmp	<- subset(reads.table, present, c(SUFFIX, id, reads, tips))
  setnames(tmp, c("id","reads","tips"), c("HOST.1","HOST.1_reads","HOST.1_tips"))
  
  
  dp <- merge(dp, tmp, by=c('SUFFIX','HOST.1'))
  setnames(tmp, c("HOST.1","HOST.1_tips","HOST.1_reads"), c("HOST.2","HOST.2_tips","HOST.2_reads"))
  dp <- merge(dp, tmp, by=c('SUFFIX','HOST.2'))
  
  #
  #	merge reads/leaves with tt
  #
  tt			<- merge(tt, dp, by=c('HOST.1','HOST.2','SUFFIX'))
}

if(nrow(tt)==0){
  cat("Failed to merge tables; e.g. file suffix ",tt$SUFFIX[1]," not found in ",summary.file,"\n",sep="")
  quit(save="no", status=1)
}
if (verbose) cat('Writing summary to file',output.file,'\n')
if(args$asRDA){
  save(tt, file=output.file)
} else {
  write.csv(tt, file=output.file, row.names=FALSE, quote=FALSE)
}

