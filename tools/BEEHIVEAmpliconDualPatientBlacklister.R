#	install missing packages

suppressMessages(require(tools, quietly=TRUE, warn.conflicts=FALSE))

command.line <- T

if(command.line){
  require(argparse, quietly=TRUE, warn.conflicts=FALSE)
  
  arg_parser = ArgumentParser(description="Look at identified multiple infections across all windows, and identify which are to be ignored entirely or reduced to largest subtree only.")

  arg_parser$add_argument("-b", "--existingBlacklistsPrefix", action="store", help="A file path and initial string identifying existing (.csv) blacklists whose output is to be appended to. Only such files whose suffix (the string after the prefix, without the file extension) matches a suffix from the dual reports will be appended to.")
  arg_parser$add_argument("threshold", action="store", help="The proportion of windows in a single amplicon that a dual infection needs to appear in to be considered genuine. Those patients whose dual infections are considered genuine are blacklisted entirely. Those who are not are reduced to their largest subtrees only.")
  arg_parser$add_argument("dualReportsPrefix", action="store", help="A file path and initial string identifying all dual infection files output from ParsimonyBasedBlacklister.R.")
  arg_parser$add_argument("newBlacklistsPrefix", action="store", help="A file path and initial string for all output blacklist files.")

  
  
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
  # BEEHIVE
  
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20161013/Blacklists/")
  
  existing.bl.prefix <- "RogueBlacklist.InWindow"
  
  duals.prefix <- "MultipleInfections.InWindow" 
  output.prefix <- "DualsBlacklist_new.InWindow"
  
  threshold <- 0.75
  
  dual.files <- list.files(dirname(duals.prefix), pattern=paste('^',basename(duals.prefix),sep=""), full.names=TRUE)
  dual.files.sans.ext <- file_path_sans_ext(dual.files)
  suffixes <- substr(dual.files.sans.ext, nchar(duals.prefix)+3, nchar(dual.files.sans.ext))
  
  expected.blacklists <- paste(existing.bl.prefix, suffixes, ".csv", sep="")
  
  observed.bl.files <- list.files(dirname(existing.bl.prefix), pattern=paste('^',basename(existing.bl.prefix),sep=""), full.names=TRUE)
  expected.but.not.seen <- setdiff(expected.blacklists, observed.bl.files)
  seen.but.not.expected <- setdiff(observed.bl.files, expected.blacklists)
  
  
}

window.count.by.patient <- list()

amps <- data.frame(starts = c(800, 1505, 4805, 5985), ends = c(2385, 5037, 7831, 9400))

total.windows.per.amp <- c(0,0,0,0)

for(suffix in suffixes){
  start <- as.numeric(unlist(strsplit(suffix, '_'))[2])
  end <- as.numeric(unlist(strsplit(suffix, '_'))[4])
  in.amp <- sapply(seq(1, 4), function (x) start >= amps$starts[x] & end <= amps$ends[x]  )
  total.windows.per.amp <- total.windows.per.amp + as.numeric(in.amp)
}

for(suffix in suffixes){
  start <- as.numeric(unlist(strsplit(suffix, '_'))[2])
  end <- as.numeric(unlist(strsplit(suffix, '_'))[4])
  
  in.amp <- sapply(seq(1, 4), function (x) start >= amps$starts[x] & end <= amps$ends[x]  )
  
  dual.file <- read.csv(paste(duals.prefix, suffix, ".csv", sep=""), stringsAsFactors = F)
  for(patient in unique(dual.file$patient)){
    if(is.null(window.count.by.patient[[patient]])){
      window.count.by.patient[[patient]] <-  as.numeric(in.amp)
    } else {
      window.count.by.patient[[patient]] <- window.count.by.patient[[patient]] + as.numeric(in.amp)
    }
  }
  if(file.exists(paste(existing.bl.prefix, suffix, ".csv", sep=""))){
    file.copy(paste(existing.bl.prefix, suffix, ".csv", sep=""), paste(output.prefix, suffix, ".csv", sep=""))
  }
}

first <- T

for(patient in labels(window.count.by.patient)[order(labels(window.count.by.patient))]){
  thresholds <- total.windows.per.amp*threshold
  
  dual.amplicons <- which(window.count.by.patient[[patient]] > thresholds)
  
  if(first){
    first <- F
    
    out.table <- c(patient, window.count.by.patient[[patient]], window.count.by.patient[[patient]]/total.windows.per.amp, length(dual.amplicons))
  } else {
    out.table <- rbind(out.table, c(patient, window.count.by.patient[[patient]], window.count.by.patient[[patient]]/total.windows.per.amp, length(dual.amplicons)))
  }
  
  if(length(dual.amplicons)>1){
    
    dual.patients <- c(dual.patients, patient)
    
    cat("Patient",patient,"meets the threshold for dual infection on amplicons",dual.amplicons,"(",window.count.by.patient[[patient]]/total.windows.per.amp,") and is blacklisted entirely.\n", sep=" ")
    # this is a dual infection and we want it removed in its entirety
    for(suffix in suffixes){
      dual.file <- read.csv(paste(duals.prefix, suffix, ".csv", sep=""), stringsAsFactors = F)
      pat.tips <- dual.file$tip.name[which(dual.file$patient==patient)]
      if(length(pat.tips)>0){
        if(file.exists(paste(output.prefix, suffix, ".csv", sep=""))){
          existing.bl <- read.table(paste(output.prefix, suffix, ".csv", sep=""), sep=",", stringsAsFactors = F, header=F)$V1
          new.bl <- unique(c(existing.bl, pat.tips))
  
        } else {
          new.bl <- pat.tips
        }
        write.table(new.bl, paste(output.prefix, suffix, ".csv", sep=""), sep=",", row.names=FALSE, col.names=FALSE, quote=F)
      }
    }
  } else {
    cat("Patient",patient,"does not meet the threshold for dual infection on any amplicon (",window.count.by.patient[[patient]]/total.windows.per.amp,") and smaller subtrees are being blacklisted.\n", sep=" ")
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