#	install missing packages

suppressMessages(require(tools, quietly=TRUE, warn.conflicts=FALSE))

command.line <- T

if(command.line){
  require(argparse, quietly=TRUE, warn.conflicts=FALSE)
  
  arg_parser = ArgumentParser(description="Look at identified multiple infections across all windows, and identify which are to be ignored entirely or reduced to largest subtree only.")

  arg_parser$add_argument("-b", "--existingBlacklistsPrefix", action="store", help="A file path and initial string identifying existing (.csv) blacklists whose output is to be appended to. Only such files whose suffix (the string after the prefix, without the file extension) matches a suffix from the dual reports will be appended to.")
  arg_parser$add_argument("-w", "--windowCount", action="store", help="The total number of windows in the analysis (will default to the total number of dual infection files if absent.)")
  arg_parser$add_argument("threshold", action="store", help="The proportion of windows that a dual infection needs to appear in to be considered genuine. Those patients whose dual infections are considered genuine are blacklisted entirely. Those who are not are reduced to their largest subtrees only.")
  arg_parser$add_argument("dualReportsPrefix", action="store", help="A file path and initial string identifying all dual infection files output from ParsimonyBasedBlacklister.R.")
  arg_parser$add_argument("newBlacklistsPrefix", action="store", help="A file path and initial string for all output blacklist files.")

  # Read in the arguments
  
  args <- arg_parser$parse_args()
  
  existing.bl.prefix <- args$existingBlacklistsPrefix
  duals.prefix <- args$dualReportsPrefix
  output.prefix <- args$newBlacklistsPrefix 
  threshold <- args$threshold
  total.windows <- as.numeric(args$windowCount)
  #print(duals.prefix)  
  dual.files <- list.files(dirname(duals.prefix), pattern=paste('^',basename(duals.prefix),sep=""), full.names=TRUE)
  dual.files.sans.ext <- file_path_sans_ext(dual.files)
  #print(dual.files)
  #print(dual.files.sans.ext)
 
  suffixes <- substr(dual.files.sans.ext, nchar(duals.prefix)+1, nchar(dual.files.sans.ext))
  #print(suffixes) 
  expected.blacklists <- paste(existing.bl.prefix, suffixes, ".csv", sep="")
  
  observed.bl.files <- list.files(dirname(existing.bl.prefix), pattern=paste('^',basename(existing.bl.prefix),sep=""), full.names=TRUE)
  expected.but.not.seen <- setdiff(expected.blacklists, observed.bl.files)
  seen.but.not.expected <- setdiff(observed.bl.files, expected.blacklists)
  
  if(length(expected.but.not.seen)>0){
    for(ebns in expected.but.not.seen){
      warning("Blacklist file ",ebns," not found. Will make a new blacklist file with this extension.")
    }
  }
  
  if(length(seen.but.not.expected)>0){
    for(sbne in expected.but.not.seen){
      warning("Blacklist file ",sbne," found but does not match a dual infections file; will be ignored.")
    }
  }
} else {

}

window.count.by.patient <- list()

for(suffix in suffixes){
  cat('Reading file',paste0(duals.prefix, suffix, ".csv"),'\n')
  dual.file <- read.csv(paste(duals.prefix, suffix, ".csv", sep=""), stringsAsFactors = F)
  for(patient in unique(dual.file$patient)){
    if(is.null(window.count.by.patient[[patient]])){
      window.count.by.patient[[patient]] <-  1
    } else {
      window.count.by.patient[[patient]] <- window.count.by.patient[[patient]]+1
    }
  }
  tmp	<- paste0(existing.bl.prefix, suffix, ".csv")
  if(file.exists(tmp) & file.size(tmp)>0){
	cat('Copy existing blacklist to',paste0(output.prefix, suffix, ".csv\n"))
    file.copy(tmp, paste(output.prefix, suffix, ".csv", sep=""))
  }
}

if(is.null(total.windows)){
  total.windows <- length(suffixes)
}
  
count.dual <- count.not <- 0
for(patient in labels(window.count.by.patient)[order(labels(window.count.by.patient))]){
  if(window.count.by.patient[[patient]]/total.windows > threshold){
    count.dual <- count.dual + 1
    cat("Patient ",patient," meets the threshold for dual infection (",window.count.by.patient[[patient]]," out of ",total.windows,") and is blacklisted entirely.\n", sep="")
    # this is a dual infection and we want it removed in its entirety
    for(suffix in suffixes){
      dual.file <- read.csv(paste(duals.prefix, suffix, ".csv", sep=""), stringsAsFactors = F)
      pat.tips <- dual.file$tip.name[which(dual.file$patient==patient)]
      if(length(pat.tips)>0){
        if(file.exists(paste(output.prefix, suffix, ".csv", sep=""))){
		  cat('Reading table to update',paste0(output.prefix, suffix, ".csv"))
          existing.bl <- read.table(paste(output.prefix, suffix, ".csv", sep=""), sep=",", stringsAsFactors = F, header=F)$V1
          new.bl <- unique(c(existing.bl, pat.tips))
  
        } else {
          new.bl <- pat.tips
        }
		cat('Writing table',paste0(output.prefix, suffix, ".csv\n"))
        write.table(new.bl, paste(output.prefix, suffix, ".csv", sep=""), sep=",", row.names=FALSE, col.names=FALSE, quote=F)
      }
    }
  } else {
    count.not <- count.not + 1
    cat("Patient ",patient," does not meet the threshold for dual infection (",window.count.by.patient[[patient]]," out of ",total.windows,") and smaller subtrees are being blacklisted.\n", sep="")
    # this is to have smaller subtrees blacklisted away on windows where they occur
    for(suffix in suffixes){
      
      dual.file <- read.csv(paste(duals.prefix, suffix, ".csv", sep=""), stringsAsFactors = F)
      pat.rows <- dual.file[which(dual.file$patient==patient),]
      if(nrow(pat.rows)>0){
        max.reads <- max(pat.rows$reads.in.subtree)
        smaller.tips <- pat.rows$tip.name[which(pat.rows$reads.in.subtree!=max.reads)]
      
        if(file.exists(paste(output.prefix, suffix, ".csv", sep=""))){
          existing.bl <- read.table(paste(output.prefix, suffix, ".csv", sep=""), sep=",", stringsAsFactors = F, header=F)$V1
          new.bl <- unique(c(existing.bl, smaller.tips))
          
        } else {
          new.bl <- smaller.tips
        }
        write.table(new.bl, paste(output.prefix, suffix, ".csv", sep=""), sep=",", row.names=FALSE, col.names=FALSE, quote=F)
      }
    }
  }
}
