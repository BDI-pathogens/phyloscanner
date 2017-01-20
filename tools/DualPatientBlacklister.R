#	install missing packages

suppressMessages(require(tools, quietly=TRUE, warn.conflicts=FALSE))

if (command.line) {
  require(argparse, quietly=TRUE, warn.conflicts=FALSE)
  
  arg_parser = ArgumentParser(description="Summarise patient statistics from each phylogeny across all windows.")
  
  arg_parser$add_argument("threshold", "What proportion of windows that a dual infection needs to appear in to be considered genuine.
                          Those patients whose dual infections are considered genuine are blacklisted entirely. Those who are not
                          are reduced to their largest subtrees only.")
  arg_parser$add_argument("dualReportsPrefix", help="A file path and initial string identifying all dual infection files output
                          from ParsimonyBasedBlacklister.R.")
  arg_parser$add_argument(-b, "existingBlacklistsPrefix", action="store", help="A file path and initial string identifying 
                          existing (.csv) blacklists whose output is to be appended to. Only such files whose suffix (the string after 
                          the prefix, without the file extension) matches a suffix from the dual reports will be appended to.")
  arg_parser$add_argument("newBlacklistsPrefix", help="A file path and initial string for all output blacklist files.")

  
  
  # Read in the arguments
  
  args <- arg_parser$parse_args()
  
  existing.bl.prefix <- args$existingBlacklistsPrefix
  duals.prefix <- args$dualReportsPrefix
  output.prefix <- args$newBlacklistsPrefix 
  threshold <- args$threshold
  
  dual.files <- list.files(dirname(duals.prefix), pattern=paste('^',basename(duals.prefix),sep=""), full.names=TRUE)
  dual.files.sans.ext <- file_path_sans_ext(dual.files)
  suffixes <- substr(dual.files.sans.ext, nchar(duals.prefix)+1, nchar(dual.files.sans.ext))
  
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
  dual.file <- read.csv(paste(duals.prefix, suffix, ".csv", sep=""), stringsAsFactors = F)
  for(patient in unique(dual.file$patient)){
    if(is.null(window.count.by.patient[[patient]])){
      window.count.by.patient[[patient]] <-  1
    } else {
      window.count.by.patient[[patient]] <- window.count.by.patient[[patient]]+1
    }
  }
  if(file.exists(paste(existing.bl.prefix, suffix, ".csv", sep=""))){
    file.copy(paste(existing.bl.prefix, suffix, ".csv", sep=""), paste(output.prefix, suffix, ".csv", sep=""))
  }
}

total.windows <- length(suffixes)

for(patient in labels(window.count.by.patient)){
  if(window.count.by.patient[[patient]]/total.windows > threshold){
    cat("Patient ",patient," meets the threshold for dual infection and is blacklisted entirely.\n", sep="")
    # this is a dual infection and we want it removed in its entirety
    for(suffix in suffixes){
      dual.file <- read.csv(paste(duals.prefix, suffix, ".csv", sep=""), stringsAsFactors = F)
      pat.tips <- dual.file$tip.name[which(dual.file$patient==patient)]
      if(length(pat.tips)>0){
        if(file.exists(paste(output.prefix, suffix, ".csv", sep=""))){
          existing.bl <- read.csv(paste(output.prefix, suffix, ".csv", sep=""), stringsAsFactors = F)$x
          new.bl <- unique(c(existing.bl, pat.tips))
  
        } else {
          new.bl <- pat.tips
        }
        write.table(new.bl, paste(output.prefix, suffix, ".csv", sep=""), sep=",", row.names=FALSE, col.names=FALSE, quote=F)
      }
    }
  } else {
    cat("Patient ",patient," does not meet the threshold for dual infection and smaller subtrees are being blacklisted.\n", sep="")
    # this is to have smaller subtrees blacklisted away on windows where they occur
    for(suffix in suffixes){
      
      dual.file <- read.csv(paste(duals.prefix, suffix, ".csv", sep=""), stringsAsFactors = F)
      pat.rows <- dual.file[which(dual.file$patient==patient),]
      if(nrow(pat.rows)>0){
        max.reads <- max(pat.rows$reads.in.subtree)
        smaller.tips <- pat.rows$tip.name[which(pat.rows$reads.in.subtree!=max.reads)]
      
        if(file.exists(paste(output.prefix, suffix, ".csv", sep=""))){
          existing.bl <- read.csv(paste(output.prefix, suffix, ".csv", sep=""), stringsAsFactors = F)$x
          new.bl <- unique(c(existing.bl, smaller.tips))
          
        } else {
          new.bl <- smaller.tips
        }
        write.table(new.bl, paste(output.prefix, suffix, ".csv", sep=""), sep=",", row.names=FALSE, col.names=FALSE, quote=F)
      }
    }
  }
}