list.of.packages <- c("argparse", "ape", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat("Please run PackageInstall.R to continue\n")
  quit(save='no')
}

command.line <- T

suppressMessages(library(data.table, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(ape, quietly=TRUE, warn.conflicts=FALSE))

if(command.line){
  suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))
  
  arg_parser = ArgumentParser(description="Look at identified multiple infections across all windows, and identify which are to be ignored entirely or reduced to largest subtree only. Requires the output of ParsimonyBasedBlacklister.R.")
  
  arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", help="Regular expression identifying tips from the dataset. Three groups, in order: patient ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name.")
  arg_parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Talk about what I'm doing.")
  arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.")
  arg_parser$add_argument("-b", "--existingBlacklistsPrefix", action="store", help="A file path and initial string identifying existing (.csv) blacklists whose output is to be appended to. Only such files whose suffix (the string after the prefix, without the file extension) matches a suffix from the dual reports will be appended to.")
  arg_parser$add_argument("-s", "--summaryFile", action="store", help="A file to write a summary of the results to, detailing window counts for all patients")
  arg_parser$add_argument("threshold", action="store", help="The proportion of windows that a dual infection needs to appear in to be considered genuine. Those patients whose dual infections are considered genuine are blacklisted entirely. Those who are not are reduced to their largest subtrees only.")
  arg_parser$add_argument("treePrefix", action="store", help="A file path and initial string identifying read trees of all analyzed windows.")
  arg_parser$add_argument("dualReportsPrefix", action="store", help="A file path and initial string identifying all dual infection files output from ParsimonyBasedBlacklister.R.")
  arg_parser$add_argument("newBlacklistsPrefix", action="store", help="A file path and initial string for all output blacklist files.")
  
  # Read in the arguments
  
  args <- arg_parser$parse_args()
  
  script.dir <- args$scriptdir
  
  source(file.path(script.dir, "TreeUtilityFunctions.R"))
  
  tree.prefix <- args$treePrefix
  existing.bl.prefix <- args$existingBlacklistsPrefix
  duals.prefix <- args$dualReportsPrefix
  output.prefix <- args$newBlacklistsPrefix 
  threshold <- args$threshold
  summary.file <- args$summaryFile
  verbose <- args$verbose
  tip.regex <- args$tipRegex

  
  dual.files <- list.files.mod(dirname(duals.prefix), pattern=paste('^',basename(duals.prefix),sep=""), full.names=TRUE)
  
  if(!is.null(tree.prefix)){
    tree.files	<- list.files.mod(dirname(tree.prefix), pattern=paste0('^',basename(tree.prefix)), full.names=TRUE)
    cat('Found tree files to determine total windows present per patient, n=', nrow(tree.files), "\n", sep="")
  }  
  
  suffixes	<- substr(dual.files, nchar(duals.prefix)+1, nchar(dual.files))
  
  suffixes <- suffixes[order(suffixes)]
  tree.files <- tree.files[order(suffixes)]
  dual.files <- dual.files[order(suffixes)]
  
  output.bls <- sapply(suffixes, function(x) paste0(output.prefix, x, sep=""))
  
  if(!is.null(existing.bl.prefix)){
    expected.blacklists <- paste(existing.bl.prefix, suffixes, sep="")
    observed.bl.files <- list.files.mod(dirname(existing.bl.prefix), pattern=paste('^',basename(existing.bl.prefix),sep=""), full.names=TRUE)
    expected.but.not.seen <- setdiff(expected.blacklists, observed.bl.files)
    seen.but.not.expected <- setdiff(observed.bl.files, expected.blacklists)
    if(length(expected.but.not.seen)>0){
      for(ebns in expected.but.not.seen){
        cat("Blacklist file ",ebns," not found. Will make a new blacklist file with this extension.\n",sep="")
      }
    }
  }
  
  fn.df <- data.frame(row.names = suffixes, tree.input = tree.files, dual.input = dual.files, blacklist.output = output.bls, stringsAsFactors = F)

  if(!is.null(existing.bl.prefix)){
    fn.df$blacklist.input <- paste(existing.bl.prefix, suffixes, sep="")
  }
  file.details <- split(fn.df, rownames(fn.df))
  
} else {
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20161013/")
  script.dir <- "/Users/twoseventwo/Documents/phylotypes/tools"
  tree.prefix <- "AllTrees/RAxML_bestTree.InWindow_"
  duals.prefix <- "Blacklists/Multiple infection files/MultipleInfections.InWindow_"
  existing.bl.prefix <- "Blacklists/Rogue blacklists/RogueBlacklist.InWindow_"
  tip.regex <- "^(.*)-[0-9].*_read_([0-9]+)_count_([0-9]+)$"
  summary.file <- "summarytest.csv"
  output.prefix <- "blacklist_test2_"
  threshold <- 0.25
}

#	read dual files output by ParsimonyBasedBlacklister

dd	<- data.table(PATIENT=rep(NA_character_,0), READS_IN_SUBTREE=rep(NA_integer_,0), TIPS_IN_SUBTREE=rep(NA_integer_,0), W_INFO=rep(NA_character_,0),  W_POTENTIAL_DUAL=rep(NA_integer_,0))
if(length(suffixes)>0){
  dd	<- lapply(suffixes, function(suffix){
    file.name <- file.details[[suffix]]$dual.input
    
    if(file.size(file.name) == 0){
      cat('Skipping file ',file.name,' as it is empty\n', sep="")
    } else {
      if(verbose) cat('Reading file',file.name,'\n')
      dual.file <- as.data.table(read.csv(file.name, stringsAsFactors = FALSE))	
      dual.file <- unique(dual.file, by=c('patient','reads.in.subtree','tips.in.subtree'))
      set(dual.file, NULL, 'tip.name', NULL)
      setnames(dual.file, colnames(dual.file), gsub('\\.','_',toupper(colnames(dual.file))))
      set(dual.file, NULL, 'W_INFO', suffix)
      set(dual.file, NULL, 'W_POTENTIAL_DUAL', 1L)
    }
  })
  dd	<- do.call('rbind',dd)
}

# #	write detailed info to file
# 
# if(!is.null(summary.file)){	
#   tmp	<- gsub('csv$','rda',summary.file)
#   cat("\nWrite detailed dual summary file to", tmp)
#   save(dd, file=tmp)
# }

#	Count number of potential dual windows by patient

tmp								<- dd[, list(N = length(unique(W_INFO))), by = 'PATIENT']
window.count.by.patient 		<- as.list(tmp$N)
names(window.count.by.patient)	<- tmp$PATIENT

#	Copy existing blacklists

if(!is.null(existing.bl.prefix)){
  for(suffix in suffixes){
    tmp	<- file.details[[suffix]]$blacklist.input
    if(file.exists(tmp) & file.size(tmp)>0){
      cat('Copying existing blacklist to ',file.details[[suffix]]$blacklist.output,"\n",sep="")
      file.copy(tmp, file.details[[suffix]]$blacklist.output )
    }
  }
}

#	
#		total.windows = #windows an individual is present
#

#	count number of unique reads in each tree file for every dual candidate patient

pat.vector <- labels(window.count.by.patient)[order(labels(window.count.by.patient))]

tree.tips <- lapply(file.details, function(x){
  read.tree(x$tree.input)$tip.label
})


total.windows <- lapply(pat.vector, function(x){
  windows.present <- sapply(tree.tips, function(y){
    any(patient.from.label(y, tip.regex)==x)
  })
  sum(as.numeric(windows.present))
})
names(total.windows) <- pat.vector


# Output function

blacklist.reads.for.patient <- function(patient){
  if(window.count.by.patient[[patient]]/total.windows[[patient]] > threshold){
    if(verbose){
      cat("Patient ",patient," meets the threshold for dual infection (",window.count.by.patient[[patient]]," out of ",total.windows[[patient]],") and will be blacklisted entirely.\n", sep="")
    }	
    # this is a dual infection and we want it removed in its entirety
    return(TRUE)
  } else {
    if(verbose){
      cat("Patient ",patient," does not meet the threshold for dual infection (",window.count.by.patient[[patient]]," out of ",total.windows[[patient]],") and smaller subtrees are being blacklisted.\n", sep="")
    }
    # this is to have smaller subtrees blacklisted away on windows where they occur
    for(suffix in suffixes){
      if(!file.exists(file.details[[suffix]]$blacklist.output)){
        file.create(file.details[[suffix]]$blacklist.output)
      }
      
      if(file.size(file.details[[suffix]]$dual.input) > 0){
        duals.table <- read.csv(file.details[[suffix]]$dual.input, stringsAsFactors=FALSE)
        pat.rows <- duals.table[which(duals.table$patient==patient),]
        if(nrow(pat.rows)>0){
          max.reads <- max(pat.rows$reads.in.subtree)
          smaller.tips <- pat.rows$tip.name[which(pat.rows$reads.in.subtree!=max.reads)]
          
          if(file.size(file.details[[suffix]]$blacklist.output)>0){
            existing.bl <- read.table(file.details[[suffix]]$blacklist.output, sep=",", stringsAsFactors=FALSE, header=FALSE)$V1
            new.bl <- unique(c(existing.bl, smaller.tips))
            
          } else {
            new.bl <- smaller.tips
          }
          if(verbose){
            cat('Writing table ',file.details[[suffix]]$blacklist.output,"\n",sep="")
          }
          write.table(new.bl, file.details[[suffix]]$blacklist.output, sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)
        }
      }
    }
    return(FALSE)
  }
}

output <- sapply(pat.vector, blacklist.reads.for.patient)

fully.blacklisted <- pat.vector[which(output)]

if(length(fully.blacklisted)>0){
  cat("Blacklisting excluded patients from all tree files...\n")
  for(suffix in suffixes){
    tip.names <- tree.tips[[suffix]]
    tip.patients <- sapply(tip.names, function(x) patient.from.label(x, tip.regex))
    if(length(intersect(tip.patients, fully.blacklisted))>0){
      if(!file.exists(file.details[[suffix]]$blacklist.output)){
        file.create(file.details[[suffix]]$blacklist.output)
      }
      tips.to.go <- tip.names[which(tip.patients %in% fully.blacklisted)]
      if(file.size(file.details[[suffix]]$blacklist.output)>0){
        existing.bl <- read.table(file.details[[suffix]]$blacklist.output, sep=",", stringsAsFactors=FALSE, header=FALSE)$V1
        new.bl <- unique(c(existing.bl, tips.to.go))
      } else {
        new.bl <- tips.to.go
      }
      if(verbose){
        cat('Writing table ',file.details[[suffix]]$blacklist.output,"\n",sep="")
      }
      write.table(new.bl, file.details[[suffix]]$blacklist.output, sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)
    }
  }
}



if(!is.null(summary.file)) {
  
  output <- lapply(pat.vector, function(x) c(window.count.by.patient[[x]]/total.windows[[x]], window.count.by.patient[[x]]))
  
  proportions <- unlist(lapply(output, "[[", 1))
  counts <- unlist(lapply(output, "[[", 2))
  
  out.df <- data.frame(patient = labels(window.count.by.patient)[order(labels(window.count.by.patient))], count = counts, proportion = proportions, stringsAsFactors = F)
  if(verbose) cat("\nWriting dual summary file to ", summary.file, "\n", sep="")
  write.csv(out.df, file=summary.file, row.names=FALSE, quote=FALSE)
}
