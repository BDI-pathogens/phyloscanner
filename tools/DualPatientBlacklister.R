list.of.packages <- c("argparse", "tools", "ape", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat("Please run PackageInstall.R to continue")
  quit(save='n')
}

command.line <- T

if(command.line){
  require(argparse, quietly=TRUE, warn.conflicts=FALSE)
  
  arg_parser = ArgumentParser(description="Look at identified multiple infections across all windows, and identify which are to be ignored entirely or reduced to largest subtree only.")
  
  arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what I'm doing.")
  arg_parser$add_argument("-b", "--existingBlacklistsPrefix", action="store", help="A file path and initial string identifying existing (.csv) blacklists whose output is to be appended to. Only such files whose suffix (the string after the prefix, without the file extension) matches a suffix from the dual reports will be appended to.")
  arg_parser$add_argument("-t", "--treePrefix", action="store", help="A file path and initial string identifying read trees of all analyzed windows. From these, we determine the number of windows in which each individual is present.")
  arg_parser$add_argument("-w", "--windowCount", action="store", help="The total number of windows in the analysis (will default to the total number of dual infection files if absent.)")
  arg_parser$add_argument("-s", "--summaryFile", action="store", help="A file to write a summary of the results to, detailing window counts for all patients")
  arg_parser$add_argument("threshold", action="store", help="The proportion of windows that a dual infection needs to appear in to be considered genuine. Those patients whose dual infections are considered genuine are blacklisted entirely. Those who are not are reduced to their largest subtrees only.")
  arg_parser$add_argument("dualReportsPrefix", action="store", help="A file path and initial string identifying all dual infection files output from ParsimonyBasedBlacklister.R.")
  arg_parser$add_argument("newBlacklistsPrefix", action="store", help="A file path and initial string for all output blacklist files.")
  
  # Read in the arguments
  
  args <- arg_parser$parse_args()
  tree.prefix <- args$treePrefix
  existing.bl.prefix <- args$existingBlacklistsPrefix
  duals.prefix <- args$dualReportsPrefix
  output.prefix <- args$newBlacklistsPrefix 
  threshold <- args$threshold
  summary.file <- args$summaryFile
  verbose <- args$verbose
  total.windows	<- NULL

  
  dual.files <- list.files(dirname(duals.prefix), pattern=paste('^',basename(duals.prefix),sep=""), full.names=TRUE)
  
  if(!is.null(args$windowCount)) {
    total.windows	<- as.numeric(args$windowCount)
  } else {
    total.windows <- length(dual.files)
  }
  
  tree.files	<- data.table(F=rep(NA_character_,0))	
  if(!is.null(tree.prefix)){
    tree.files	<- data.table(F=list.files(dirname(tree.prefix), pattern=paste0('^',basename(tree.prefix)), full.names=TRUE))
    cat('Found tree.files to determine total.windows per patient, n=', nrow(tree.files))
  }  
  
  suffixes	<- substr(dual.files, nchar(duals.prefix)+1, nchar(dual.files))
  expected.blacklists <- paste(existing.bl.prefix, suffixes, sep="")
  observed.bl.files <- list.files(dirname(existing.bl.prefix), pattern=paste('^',basename(existing.bl.prefix),sep=""), full.names=TRUE)
  expected.but.not.seen <- setdiff(expected.blacklists, observed.bl.files)
  seen.but.not.expected <- setdiff(observed.bl.files, expected.blacklists)
  if(length(expected.but.not.seen)>0){
    for(ebns in expected.but.not.seen){
      warning("Blacklist file ",ebns," not found. Will make a new blacklist file with this extension.")
    }
  }
} else {
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20161013/")
  
  tree.prefix <- "AllTrees/RAxML_bestTree.InWindow_"
  existing.bl.prefix <- "Blacklists/Rogue blacklists/RogueBlacklist.InWindow_"
  duals.prefix <- "Blacklists/Multiple infection files/MultipleInfections.InWindow_"
  output.prefix <- "Blacklists/Dual blacklists/DualBlacklist.InWindow_"
  threshold <- 1
  summary.file <- "Blacklists/Dual blacklists/DualsSummary.csv"
  verbose <- TRUE
  total.windows	<- NULL
  
  dual.files <- list.files(dirname(duals.prefix), pattern=paste('^',basename(duals.prefix),sep=""), full.names=TRUE)
  

  tree.files	<- data.table(F=rep(NA_character_,0))	
  if(!is.null(tree.prefix)){
    tree.files	<- data.table(F=list.files(dirname(tree.prefix), pattern=paste0('^',basename(tree.prefix)), full.names=TRUE))
    cat('Found tree.files to determine total.windows per patient, n=', nrow(tree.files))
  }  
  
  suffixes	<- substr(dual.files, nchar(duals.prefix)+1, nchar(dual.files))
  expected.blacklists <- paste(existing.bl.prefix, suffixes, sep="")
  observed.bl.files <- list.files(dirname(existing.bl.prefix), pattern=paste('^',basename(existing.bl.prefix),sep=""), full.names=TRUE)
  expected.but.not.seen <- setdiff(expected.blacklists, observed.bl.files)
  seen.but.not.expected <- setdiff(observed.bl.files, expected.blacklists)
  if(length(expected.but.not.seen)>0){
    for(ebns in expected.but.not.seen){
      warning("Blacklist file ",ebns," not found. Will make a new blacklist file with this extension.")
    }
  }
}

suppressMessages(library(data.table, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(ape, quietly=TRUE, warn.conflicts=FALSE))

#	read dual files output by ParsimonyBasedBlacklister

dd	<- data.table(PATIENT=rep(NA_character_,0), READS_IN_SUBTREE=rep(NA_integer_,0), TIPS_IN_SUBTREE=rep(NA_integer_,0), W_INFO=rep(NA_character_,0),  W_POTENTIAL_DUAL=rep(NA_integer_,0))
if(length(suffixes)>0){
  dd	<- lapply(suffixes, function(suffix){
    file.name <- paste0(duals.prefix, suffix)
    
    if(file.size(file.name) == 0){
      cat('Skipping file ',file.name,' as it is empty\n')
    } else {
    
      cat('Reading file',file.name,'\n')
      dual.file <- as.data.table(read.csv(paste(duals.prefix, suffix, sep=""), stringsAsFactors = FALSE))	
      dual.file <- unique(dual.file, by=c('patient','reads.in.subtree','tips.in.subtree'))
      set(dual.file, NULL, 'tip.name', NULL)
      setnames(dual.file, colnames(dual.file), gsub('\\.','_',toupper(colnames(dual.file))))
      set(dual.file, NULL, 'W_INFO', suffix)
      set(dual.file, NULL, 'W_POTENTIAL_DUAL', 1L)
    }
  })
  dd	<- do.call('rbind',dd)
}

#	write detailed info to file

if(!is.null(summary.file)){	
  tmp	<- gsub('csv$','rda',summary.file)
  cat("\nWrite detailed dual summary file to", tmp)
  save(dd, file=tmp)
}

#	Count number of potential dual windows by patient

tmp								<- dd[, list(N= length(unique(W_INFO))), by='PATIENT']
window.count.by.patient 		<- as.list(tmp$N)
names(window.count.by.patient)	<- tmp$PATIENT

#	Copy existing blacklists

for(suffix in suffixes){
  tmp	<- paste0(existing.bl.prefix, suffix, ".csv")
  if(file.exists(tmp) & file.size(tmp)>0)
  {
    cat('Copying existing blacklist to',paste0(output.prefix, suffix, ".csv\n"))
    file.copy(tmp, paste(output.prefix, suffix, ".csv", sep=""))
  }
}

#	If window total in input args, convert window total to named numeric 

if(!is.null(total.windows))
{
  total.windows			<- rep(total.windows, length(unique(dual.file$patient)))
  names(total.windows)	<- unique(dual.file$patient)
}
#
#	determine window total if no input arg and tree files: 
#		total.windows= #windows an individual is present
#
if(nrow(tree.files)>0 & is.null(total.windows)){
  #	count number of unique reads in each tree file for every dual candidate patient
  tmp						<- tree.files[, {
    tmp		<- data.table(TAXA=read.tree(F)$tip.label)
    tmp2	<- labels(window.count.by.patient)
    tmp		<- sapply(tmp2, function(patient)	nrow(subset(tmp, regexpr(patient, TAXA)>0))	)
    list(POT_DUAL_ID= tmp2, UNIQUE_READS=tmp)
  }, by='F']						
  #	count number of windows in which dual candidate patient is present
  tmp						<- tmp[, list( WIN_N= length(which(UNIQUE_READS>0)) ), by='POT_DUAL_ID']
  #	convert to named numeric
  total.windows			<- tmp$WIN_N
  names(total.windows)	<- tmp$POT_DUAL_ID
}

#	Determine window total if no input arg and no tree files: 
#	by default length of dual files

if(nrow(tree.files)==0 & is.null(total.windows)){
  total.windows			<- rep(length(suffixes), length(unique(dual.file$patient)))
  names(total.windows)	<- unique(dual.file$patient)
}

# Output function

blacklist.reads.for.patient <- function(patient){
  if(window.count.by.patient[[patient]]/total.windows[patient] > threshold){
    if(verbose){
      cat("Patient ",patient," meets the threshold for dual infection (",window.count.by.patient[[patient]]," out of ",total.windows[patient],") and is blacklisted entirely.\n", sep="")
    }	
    # this is a dual infection and we want it removed in its entirety
    for(suffix in suffixes){
      dual.file <- read.csv(paste(duals.prefix, suffix, sep=""), stringsAsFactors=FALSE)
      pat.tips <- dual.file$tip.name[which(dual.file$patient==patient)]
      if(length(pat.tips)>0){
        if(file.exists(paste(output.prefix, suffix, sep=""))){
          cat('Reading table to update',paste0(output.prefix, suffix),'\n')
          existing.bl <- read.table(paste(output.prefix, suffix, sep=""), sep=",", stringsAsFactors=FALSE, header=FALSE)$V1
          new.bl <- unique(c(existing.bl, pat.tips))
        } else {
          new.bl <- pat.tips
        }
        if(verbose){
          cat('Writing table',paste0(output.prefix, suffix, "\n"))
        }
        write.table(new.bl, paste(output.prefix, suffix, sep=""), sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)
      }
    }
  } else {
    if(verbose){
      cat("Patient ",patient," does not meet the threshold for dual infection (",window.count.by.patient[[patient]]," out of ",total.windows[patient],") and smaller subtrees are being blacklisted.\n", sep="")
    }
    # this is to have smaller subtrees blacklisted away on windows where they occur
    for(suffix in suffixes){
      
      dual.file <- read.csv(paste(duals.prefix, suffix, sep=""), stringsAsFactors=FALSE)
      pat.rows <- dual.file[which(dual.file$patient==patient),]
      if(nrow(pat.rows)>0){
        max.reads <- max(pat.rows$reads.in.subtree)
        smaller.tips <- pat.rows$tip.name[which(pat.rows$reads.in.subtree!=max.reads)]
        
        if(file.exists(paste(output.prefix, suffix, sep=""))){
          existing.bl <- read.table(paste(output.prefix, suffix, sep=""), sep=",", stringsAsFactors=FALSE, header=FALSE)$V1
          new.bl <- unique(c(existing.bl, smaller.tips))
          
        } else {
          new.bl <- smaller.tips
        }
        if(verbose){
          cat('Writing table',paste0(output.prefix, suffix, "\n"))
        }
        write.table(new.bl, paste(output.prefix, suffix, sep=""), sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)
      }
    }
  }
  c(window.count.by.patient[[patient]]/total.windows[patient], window.count.by.patient[[patient]])
}

pat.vector <- labels(window.count.by.patient)[order(labels(window.count.by.patient))]

output <- lapply(pat.vector, blacklist.reads.for.patient)

proportions <- unlist(lapply(output, "[[", 1))
counts <- unlist(lapply(output, "[[", 2))

if(!is.null(summary.file)){
  out.df <- data.frame(patient = labels(window.count.by.patient)[order(labels(window.count.by.patient))], count = counts, proportion = proportions, stringsAsFactors = F)
  cat("\nWrite dual summary file to", summary.file)
  write.csv(out.df, file=summary.file, row.names=FALSE, quote=FALSE)
}