#	install missing packages

suppressMessages(require(tools, quietly=TRUE, warn.conflicts=FALSE))

command.line <- T

if(command.line){
  require(argparse, quietly=TRUE, warn.conflicts=FALSE)
  
  arg_parser = ArgumentParser(description="Look at identified multiple infections across all windows, and identify which are to be ignored entirely or reduced to largest subtree only.")
  
  arg_parser$add_argument("-s", "--summaryFile", action="store", help="A file to write a summary of the results to, detailing window counts for all patients")
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
  summary.file <- args$summaryFile
  
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
  
  existing.bl.prefix <- "RogueBlacklist_alt.InWindow"
  
  duals.prefix <- "FullInfoRB_alt.InWindow" 
  output.prefix <- "DualsBlacklist_alt.InWindow"
  
  threshold <- 0.75
  
  dual.files <- list.files(dirname(duals.prefix), pattern=paste('^',basename(duals.prefix),sep=""), full.names=TRUE)
  dual.files.sans.ext <- file_path_sans_ext(dual.files)
  suffixes <- substr(dual.files.sans.ext, nchar(duals.prefix)+3, nchar(dual.files.sans.ext))
  
  expected.blacklists <- paste(existing.bl.prefix, suffixes, ".csv", sep="")
  
  observed.bl.files <- list.files(dirname(existing.bl.prefix), pattern=paste('^',basename(existing.bl.prefix),sep=""), full.names=TRUE)
  expected.but.not.seen <- setdiff(expected.blacklists, observed.bl.files)
  seen.but.not.expected <- setdiff(observed.bl.files, expected.blacklists)
  
  summary.file <- "duals_summary_alt.csv"
  
}

numerators.by.patient <- list()
denominators.by.patient <- list()

for(patient.no in 1:length(dual.file$patient)){
  patient <- dual.file$patient[patient.no]
  numerators.by.patient[[patient]] <- rep(0,8)
  denominators.by.patient[[patient]] <- rep(0,8)
}

amps <- data.frame(starts = c(480, 1485, 4783, 5967), ends = c(2407, 5058, 7848, 9517))

for(suffix in suffixes){
  start <- as.numeric(unlist(strsplit(suffix, '_'))[2])
  end <- as.numeric(unlist(strsplit(suffix, '_'))[4])
  
  in.amp <- vector()
  
  in.amp[1] <- start < amps$starts[2]
  in.amp[2] <- start >= amps$starts[2] & end <= amps$ends[1]
  in.amp[3] <- end > amps$ends[1] & start < amps$starts[3]
  in.amp[4] <- start >= amps$starts[3] & end <= amps$ends[2]
  in.amp[5] <- end > amps$ends[2] & start < amps$starts[4]
  in.amp[6] <- start >= amps$starts[4] & end <= amps$ends[3]
  in.amp[7] <- end > amps$ends[3]
  in.amp[8] <- T
  
 if(sum(in.amp)!=2){
   stop(suffix)
 }

  dual.file <- read.csv(paste(duals.prefix, suffix, ".csv", sep=""), stringsAsFactors = F)
  for(patient.no in 1:length(dual.file$patient)){
    patient <- dual.file$patient[patient.no]
    if(dual.file$status[patient.no] == "dual"){
      numerators.by.patient[[patient]] <- numerators.by.patient[[patient]] + as.numeric(in.amp)
    }
    if(dual.file$status[patient.no] %in% c("dual", "kept")){
      denominators.by.patient[[patient]] <- denominators.by.patient[[patient]] + as.numeric(in.amp)
    }
  }
  if(file.exists(paste(existing.bl.prefix, suffix, ".csv", sep=""))){
    file.copy(paste(existing.bl.prefix, suffix, ".csv", sep=""), paste(output.prefix, suffix, ".csv", sep=""))
  }
}

first <- T

dual.patients <- vector()

for(patient in dual.file$patient){
  thresholds <- total.windows.per.amp*threshold
  
  dual.amplicons <- which(window.count.by.patient[[patient]] > thresholds)
  
  if(first){
    first <- F
    out.table <- c(numerators.by.patient[[patient]], denominators.by.patient[[patient]], length(dual.amplicons))
  } else {
    out.table <- rbind(out.table, c(numerators.by.patient[[patient]], denominators.by.patient[[patient]], length(dual.amplicons)))
  }
  
  # if(length(dual.amplicons)>1){
  #   
  #   dual.patients <- c(dual.patients, patient)
  #   
  #   cat("Patient",patient,"meets the threshold for dual infection on amplicons",dual.amplicons,"(",window.count.by.patient[[patient]]/total.windows.per.amp,") and is blacklisted entirely.\n", sep=" ")
  #   # this is a dual infection and we want it removed in its entirety
  #   for(suffix in suffixes){
  #     dual.file <- read.csv(paste(duals.prefix, suffix, ".csv", sep=""), stringsAsFactors = F)
  #     pat.tips <- dual.file$tip.name[which(dual.file$patient==patient)]
  #     if(length(pat.tips)>0){
  #       if(file.exists(paste(output.prefix, suffix, ".csv", sep=""))){
  #         existing.bl <- read.table(paste(output.prefix, suffix, ".csv", sep=""), sep=",", stringsAsFactors = F, header=F)$V1
  #         new.bl <- unique(c(existing.bl, pat.tips))
  #         
  #       } else {
  #         new.bl <- pat.tips
  #       }
  #       write.table(new.bl, paste(output.prefix, suffix, ".csv", sep=""), sep=",", row.names=FALSE, col.names=FALSE, quote=F)
  #     }
  #   }
  # } else {
  #   cat("Patient",patient,"does not meet the threshold for dual infection on any amplicon (",window.count.by.patient[[patient]]/total.windows.per.amp,") and smaller subtrees are being blacklisted.\n", sep=" ")
  #   # this is to have smaller subtrees blacklisted away on windows where they occur
  #   for(suffix in suffixes){
  #     
  #     dual.file <- read.csv(paste(duals.prefix, suffix, ".csv", sep=""), stringsAsFactors = F)
  #     pat.rows <- dual.file[which(dual.file$patient==patient),]
  #     if(nrow(pat.rows)>0){
  #       max.reads <- max(pat.rows$reads.in.subtree)
  #       smaller.tips <- pat.rows$tip.name[which(pat.rows$reads.in.subtree!=max.reads)]
  #       
  #       if(file.exists(paste(output.prefix, suffix, ".csv", sep=""))){
  #         existing.bl <- read.table(paste(output.prefix, suffix, ".csv", sep=""), sep=",", stringsAsFactors = F, header=F)$V1
  #         new.bl <- unique(c(existing.bl, smaller.tips))
  #         
  #       } else {
  #         new.bl <- smaller.tips
  #       }
  #       write.table(new.bl, paste(output.prefix, suffix, ".csv", sep=""), sep=",", row.names=FALSE, col.names=FALSE, quote=F)
  #     }
  #   }
  # }
}

if(!is.null(summary.file)){
  row.names(out.table) <- 1:nrow(out.table)
  out.df <- data.frame(out.table[,1:16], stringsAsFactors = F)
  out.df$patient <- dual.file$patient
  out.df <- out.df[,c(17, 1:16)]
  colnames(out.df) <- c("patient", "amp.1u.num", "amp.12o.num", "amp.2u.num", "amp.23o.num", "amp.3u.num", "amp.34o.num", "amp.4u.num", "total.num", 
                        "amp.1u.denom", "amp.12o.denom", "amp.2u.denom", "amp.23o.denom", "amp.3u.denom", "amp.34o.denom", "amp.4u.denom", "total.denom")
  write.csv(out.df, summary.file, row.names = F, quote = F)
}