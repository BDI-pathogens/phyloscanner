#!/usr/bin/env Rscript

source("/Users/twoseventwo/Documents/phylotypes/TransmissionUtilityFunctions.R")
source("/Users/twoseventwo/Documents/phylotypes/SummariseTrees_funcs.R")



command.line <- FALSE
if (command.line) {
  require(optparse)
  
  option_list = list(
    make_option(c("-i", "--idFile"), type="character", default=NULL, 
                help="The path of a file containing a list of all the patient IDs to be summarised"),
    make_option(c("-t", "--treeFiles"), type="character", default=NULL, 
                help="A list of phylogenies in Newick format, separated by colons"),
    make_option(c("-b", "--blacklistFiles"), type="character", default=NULL, 
                help="A list of CSV files outlining tips to ignore in each phylogeny. Multiple files 
                per phylogeny are allowed. Semicolons separate blacklists for the same phylogeny,
                colons separate the sets of blacklists for each phylogeny."),
    make_option(c("-s", "--splitsFiles"), type="character", default=NULL, 
                help="A list of CSV files, output from SplitPatientsToSubtrees.R, detailing the groups 
                of tips into which each patients' reads have been split in order to enable consistent 
                partitioning tree nodes."),
    make_option(c("-w", "--windows"), type="character", default=NULL, 
                help="The window in the genome which each tree file, in turn, is constructed from. 
                In each window this is given as n-m where n is the start and m the end; coordinates for 
                each window are separated by a colon"),
    make_option(c("-x", "--tipRegex"), type="character", default=NULL, 
                help="Regular expression identifying tips from the dataset. Three groups: patient ID,
                read ID, and read count."),
    make_option(c("-r", "--refSeqName"), type="character", default=NULL, 
                help="Reference sequence label (if unspecified, tree will be assumed to be already rooted)"),
  )
  
  # Read in the arguments

  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)

  id.file <- opt$idFile
  root.name <- opt$refSeqName
  tip.regex <- opt$tipRegex
  
  tree.files <- unlist(strsplit(opt$treeFiles, ":"))
  splits.files <- unlist(strsplit(opt$splitsFiles, ":"))
  
  blacklist.files <- lapply(strsplit(opt$blacklistFiles, ":"), function(x) strsplit(x, ";"))
  
  windows <- unlist(strsplit(opt$windows, ":"))
  
  if(length(tree.files)!=length(splits.files)){
    stop("Number of tree files and number of splits files differ")
  }
  if(length(tree.files)!=length(blacklist.files)){
    stop("Number of tree files and number of sets of blacklists differ")
  }
  if(length(tree.files)!=length(windows)){
    stop("Number of tree files and number of windows differ")
  }
  
  
} else {
  # Some options to run within R
  #setwd("/Users/cfraser/Dropbox (Infectious Disease)/PROJECTS/BEEHIVE/phylotypes")
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/")
  #id.file <- "inputs/ListOfIDs_IDUs_IncNotOK.txt"
  #id.file <- "inputs/ListOfIDs_SubtypeCclade.txt"
  #id.file <- "inputs/ListOfIDs_CFRandom.txt"
  id.file <- "ListOfIDs_run20160517.txt"

  seed <- 1
  id.delimiter <- "-"
  window.width <- 350
  #output.dir <- "CF_TREES/IDUs"
  #output.dir <- "CF_TREES/SubtypeCclade"
  #output.dir <- "CF_TREES/CF_RandomSelection"
  output.dir <- "/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160517/R_output"
  #input.dir<-"PhylotypeOutput/SubtypeCclade"
  #input.dir<-"/Users/cfraser/Dropbox (Infectious Disease)/PROJECTS/BEEHIVE/Christophe_notebook/phylotypes/run20151011/"
  input.dir<-"/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160517/"
  root.name<- "C.BW.00.00BW07621.AF443088"
  #tree.files <- list.files("PhylotypeOutput/IDUs", pattern="*.tree")
  tree.files <- sort(list.files(input.dir, pattern="RAxML_bestTree.*\\.tree"))
  blacklist.files <-  sort(list.files(input.dir, pattern="Blacklist.*\\.csv")) 
  splits.files.cons <-  sort(list.files(input.dir, pattern="*\\_consSubtrees.csv"))
  splits.files.RS <-  sort(list.files(input.dir, pattern="*\\_rsSubtrees.csv"))
  # has to include \\. otherwise ignores the dot
  tip.regex <- "^(BEE[0-9][0-9][0-9][0-9])-[0-9].*_read_([0-9]+)_count_([0-9]+)$"
}

AlignPlots <- function(...) {
  LegendWidth <- function(x) x$grobs[[8]]$grobs[[1]]$widths[[4]]
  
  plots.grobs <- lapply(list(...), ggplotGrob)
  
  max.widths <- do.call(unit.pmax, lapply(plots.grobs, "[[", "widths"))
  plots.grobs.eq.widths <- lapply(plots.grobs, function(x) {
    x$widths <- max.widths
    x
  })
  
  legends.widths <- lapply(plots.grobs, LegendWidth)
  max.legends.width <- do.call(max, legends.widths)
  plots.grobs.eq.widths.aligned <- lapply(plots.grobs.eq.widths, function(x) {
    if (is.gtable(x$grobs[[8]])) {
      x$grobs[[8]] <- gtable_add_cols(x$grobs[[8]],
                                      unit(abs(diff(c(LegendWidth(x),
                                                      max.legends.width))),
                                           "mm"))
    }
    x
  })
  
  plots.grobs.eq.widths.aligned
}

require(phytools)
require(dplyr)
require(ggplot2)
require(reshape)
require(gtable)
require(grid)
require(gridExtra)
require(RColorBrewer)

# Read in the IDs. Remove duplicates. Shuffle their order if desired.
ids <- scan(id.file, what="", sep="\n", quiet=TRUE, skip = 1)


num.ids <- length(ids)
if (num.ids == 0) {
  cat(paste("No IDs found in ", id.file, ". Quitting.\n", sep=""))
  quit("no", 1)
}
ids <- unique(ids)

ids.col <- vector()
windows.col <- vector()
leaves.col <- vector()
reads.col <- vector()
subtrees.counts.col.cons <- vector()
subtrees.counts.col.RS <- vector()
clades.counts.col <- vector()
overall.rtt.col <- vector()
largest.rtt.col.RS <- vector()
largest.rtt.col.cons <- vector()
largest.pat.dist.col.RS <- vector()
largest.pat.dist.col.cons <- vector()
patristic.variance.col <- vector()

read.proportions.cons <- list()
read.proportions.RS <- list()

max.splits.cons <- 0
max.splits.RS <- 0

for(window.no in seq(1, length(tree.files))){

  tree.file.name <- tree.files[window.no]
  blacklist.file.names <- blacklist.files[[window.no]]
  splits.file.name <- splits.files[window.no]
  
  window <- windows[window.no]
  start <- unlist(strsplit(window, "-"))[1]
  end <- unlist(strsplit(window, "-"))[2]
  
  # Read the tree
  
  tree <- read.tree(file=tree.file.name)
  
  # Root the tree
  
  tree<-root(phy = tree,outgroup = root.name)
  num.tips <- length(tree$tip.label)
  
  blacklist <- vector()
  
  for(blacklist.file.name in blacklist.file.names){
    if(file.exists(blacklist.file.name)){
      blacklisted.tips <- read.table(blacklist.file.name, sep=",", header=F, stringsAsFactors = F)
      if(length(blacklisted.tips)>0){
        blacklist <- c(blacklist, sapply(blacklisted.tips, get.tip.no, tree=tree))
      }
    } else {
      warning(paste("File ",blacklist.file.name," does not exist; skipping.",paste=""))
    }
  }
  
  splits.table <- read.table(splits.file.name, sep=",", stringsAsFactors = F, header=T)

  
  splits.table$reads <- as.numeric(unlist(strsplit(splits.table$tip.names,"count_"))[seq(2,2*nrow(splits.table),2)])

  dummy <- resolveTreeIntoPatientClades(tree, ids, tip.regex, blacklist)
  clade.mrcas.by.patient <- dummy$clade.mrcas.by.patient
  all.clades.by.patient <- dummy$clades.by.patient
  
  # Associate each tip to its patient.
  patient.tips <- list()
  for (id in ids) patient.tips[[id]] <- vector()
  for (tip.no in seq(1, length(tree$tip.label))) {
    tip.label <- tree$tip.label[tip.no]
    id <- patient.from.label(tip.label, tip.regex)
    if (id %in% ids & !(tip.no %in% blacklist)) {
      patient.tips[[id]] <- c(patient.tips[[id]], tip.label)
    }
  }
  
  # This is exceptionally bad code - how else to pull out window coordinate though?
  window.start <- as.numeric(strsplit(tree.file.name, "_")[[1]][3])
  window.end <- window.start + window.width
  
  ids.col <- c(ids.col, ids)
  windows.col <- c(windows.col, rep(window.start, num.ids))
  
  read.proportions.cons.this.window <- list()
  read.proportions.RS.this.window <- list()
  
  read.proportions.cons.this.window$window.start <- window.start
  read.proportions.RS.this.window$window.start <- window.start
  
  for (i in 1:num.ids) {
    cat(window.start, " - ", ids[i], "\n", sep="")
    id <- ids[i]

    num.leaves <- length(patient.tips[[id]])
    
    leaves.col <- c(leaves.col, num.leaves)
    
    reads.per.tip <- vector()
    
    for (tip in patient.tips[[id]]) reads.per.tip <- c(reads.per.tip, as.numeric(unlist(strsplit(tip,"count_"))[2]))
    
    num.reads <- sum(reads.per.tip)
    
    reads.col <- c(reads.col, num.reads)
    
    num.subtrees.cons <- length(unique(splits.table.cons[which(splits.table.cons$orig.patients==id),]$patient.splits))
    num.subtrees.RS <- length(unique(splits.table.RS[which(splits.table.RS$orig.patients==id),]$patient.splits))
    
    subtrees.counts.col.cons <- c(subtrees.counts.col.cons, num.subtrees.cons)
    subtrees.counts.col.RS <- c(subtrees.counts.col.RS, num.subtrees.RS)
    
    num.clades <- length(all.clades.by.patient[[id]])
    clades.counts.col <- c(clades.counts.col, num.clades)
    
    all.tips <- patient.tips[[id]]

    if(length(all.tips)>1){
      
      subtree.all <- drop.tip(tree, tip=tree$tip.label[!(tree$tip.label %in% all.tips)])
      
      #just in case pruning changes the tip order
      
      reads.per.tip <- vector()
      
      for (tip in subtree.all$tip.label) reads.per.tip <- c(reads.per.tip, as.numeric(unlist(strsplit(tip,"count_"))[2]))
      
      names(reads.per.tip) <- subtree.all$tip.label
      
      overall.rtt <- calcMeanRootToTip(subtree.all, reads.per.tip)
      
      pat.distances <- cophenetic(subtree.all)
      patristic.variance <- var(pat.distances[upper.tri(pat.distances)])/((mean(pat.distances[upper.tri(pat.distances)]))^2)
    } else {
      overall.rtt <- 0
      patristic.variance <- NA
    }
    
    overall.rtt.col <- c(overall.rtt.col, overall.rtt)
    patristic.variance.col <- c(patristic.variance.col, patristic.variance)
    
    if(num.subtrees.RS <= 1){
      largest.rtt.RS <- overall.rtt
      if(num.subtrees.RS==0){
        mean.pat.RS <- NA
      } else if(num.subtrees.RS==1) {
        if(length(all.tips)==1){
          mean.pat.RS <- NA
        } else {
          pat.distances <- cophenetic(subtree.all)
          mean.pat.RS <- mean(pat.distances[upper.tri(pat.distances)])
        }
      }
    } else {
      relevant.reads <- splits.table.RS[which(splits.table.RS$orig.patients==id),]
      splits <- unique(relevant.reads$patient.splits)
      reads.per.split <- sapply(splits, function(x) sum(splits.table.RS[which(splits.table.RS$patient.splits==x),]$reads ) )
      winner <- splits[which(reads.per.split==max(reads.per.split))]
      if(length(winner)>1){
        cat("Patient ",id," has a joint winner in window ", window.start, " (RS)\n", sep="" )
        #need to talk about this...
        winner <- winner[1]
      }
      winner.tips <- relevant.reads[which(relevant.reads$patient.splits==winner),]$tip.names
      if(length(winner.tips)==1){
        largest.rtt.RS <- 0
        mean.pat.RS <- NA
      } else {
        subtree.RS <- drop.tip(tree, tip=tree$tip.label[!(tree$tip.label %in% winner.tips)])
        
        #just in case pruning changes the tip order
        
        reads.per.tip <- vector()
        
        for (tip in subtree.RS$tip.label) reads.per.tip <- c(reads.per.tip, as.numeric(unlist(strsplit(tip,"count_"))[2]))
        
        names(reads.per.tip) <- subtree.RS$tip.label
        
        largest.rtt.RS <- calcMeanRootToTip(subtree.RS, reads.per.tip)
        pat.distances <- cophenetic(subtree.RS)
        mean.pat.RS <- mean(pat.distances[upper.tri(pat.distances)])
      }
    }
    
    largest.rtt.col.RS <- c(largest.rtt.col.RS, largest.rtt.RS)
    largest.pat.dist.col.RS <- c(largest.pat.dist.col.RS, mean.pat.RS)
    
    if(num.subtrees.cons <= 1){
      largest.rtt.cons <- overall.rtt
      if(num.subtrees.cons==0){
        mean.pat.cons <- NA
      } else if(num.subtrees.cons==1) {
        if(length(all.tips)==1){
          mean.pat.cons <- NA
        } else {
          pat.distances <- cophenetic(subtree.all)
          mean.pat.cons <- mean(pat.distances[upper.tri(pat.distances)])
        }
      }
    } else {
      relevant.reads <- splits.table.cons[which(splits.table.cons$orig.patients==id),]
      splits <- unique(relevant.reads$patient.splits)
      reads.per.split <- sapply(splits, function(x) sum(splits.table.cons[which(splits.table.cons$patient.splits==x),]$reads ) )
      winner <- splits[which(reads.per.split==max(reads.per.split))]
      if(length(winner)>1){
        cat("Patient ",id," has a joint winner in window ", window.start, " (conservative)\n", sep="" )
        #need to talk about this...
        winner <- winner[1]
      }
      winner.tips <- relevant.reads[which(relevant.reads$patient.splits==winner),]$tip.names
      if(length(winner.tips)==1){
        largest.rtt.cons <- 0
        mean.pat.cons <- NA
      } else {
        subtree.cons <- drop.tip(tree, tip=tree$tip.label[!(tree$tip.label %in% winner.tips)])
        
        #just in case pruning changes the tip order
        
        reads.per.tip <- vector()
        
        for (tip in subtree.cons$tip.label) reads.per.tip <- c(reads.per.tip, as.numeric(unlist(strsplit(tip,"count_"))[2]))
        
        names(reads.per.tip) <- subtree.cons$tip.label
        
        largest.rtt.cons <- calcMeanRootToTip(subtree.cons, reads.per.tip)
        pat.distances <- cophenetic(subtree.cons)
        mean.pat.cons <- mean(pat.distances[upper.tri(pat.distances)])
      }
      
      

    }

    largest.rtt.col.cons <- c(largest.rtt.col.cons, largest.rtt.cons)
    largest.pat.dist.col.cons <- c(largest.pat.dist.col.cons, mean.pat.cons)
    
    this.pat.splits.cons <- splits.table.cons[which(splits.table.cons$orig.patients==id),]
    
    if(nrow(this.pat.splits.cons)>0){
      this.pat.reads.by.split.cons <- aggregate(this.pat.splits.cons$reads, by=list(Category=this.pat.splits.cons$patient.splits), sum)
      
      if(nrow(this.pat.reads.by.split.cons) > max.splits.cons){
        max.splits.cons <- nrow(this.pat.reads.by.split.cons>max.splits.cons)
      }
      
      this.pat.reads.by.split.cons <- this.pat.reads.by.split.cons[order(this.pat.reads.by.split.cons$x, decreasing = T),]
      this.pat.reads.by.split.cons$proportion <- this.pat.reads.by.split.cons$x/sum(this.pat.reads.by.split.cons$x)
      read.proportions.cons.this.window[[id]] <- this.pat.reads.by.split.cons$proportion
    } else {
      read.proportions.cons.this.window[[id]] <- NA
    }
    
    this.pat.splits.RS <- splits.table.RS[which(splits.table.RS$orig.patients==id),]
    
    if(nrow(this.pat.splits.RS)>0){
      this.pat.reads.by.split.RS <- aggregate(this.pat.splits.RS$reads, by=list(Category=this.pat.splits.RS$patient.splits), sum)
      
      if(nrow(this.pat.reads.by.split.RS) > max.splits.RS){
        max.splits.RS <- nrow(this.pat.reads.by.split.RS>max.splits.RS)
      }
      
      this.pat.reads.by.split.RS <- this.pat.reads.by.split.RS[order(this.pat.reads.by.split.RS$x, decreasing = T),]
      this.pat.reads.by.split.RS$proportion <- this.pat.reads.by.split.RS$x/sum(this.pat.reads.by.split.RS$x)
      read.proportions.RS.this.window[[id]] <- this.pat.reads.by.split.RS$proportion
    } else {
      read.proportions.RS.this.window[[id]] <- NA
    }
  
  }
  read.proportions.cons[[window.no]] <- read.proportions.cons.this.window
  read.proportions.RS[[window.no]] <- read.proportions.RS.this.window

}

pat.stats <- data.frame(window = windows.col, patient = ids.col, leaves = leaves.col, 
                        reads = reads.col, subtrees.cons = subtrees.counts.col.cons, 
                        subtrees.RS = subtrees.counts.col.RS, clades = clades.counts.col,
                        overall.rtt = overall.rtt.col, largest.rtt.cons = largest.rtt.col.cons,
                        largest.rtt.RS = largest.rtt.col.RS, largest.pat.dist.cons = largest.pat.dist.col.cons,
                        largest.pat.dist.RS = largest.pat.dist.col.RS, patristic.variance = patristic.variance.col)

first <- T
for(split in seq(1, max.splits.cons)){
  split.col <- vector()
  for(window.no in seq(1, length(tree.files))){
    for(id in ids){
      split.col <- c(split.col, read.proportions.cons[[window.no]][[id]][split])
    }
  }
  if(first){
    first <- F
    cons.splits.props <- split.col
  } else {
    cons.splits.props <- cbind(cons.splits.props, split.col)
  }
}

cons.splits.props <- data.frame(cons.splits.props)
cons.splits.props$window <- windows.col
cons.splits.props$patient <- ids.col
cons.splits.props <- cons.splits.props[,c(ncol(cons.splits.props)-1, ncol(cons.splits.props), seq(1, ncol(cons.splits.props)-2))]

colnames(cons.splits.props) <- c("window", "patient", paste("prop.gp.cons.",seq(1,max.splits.cons),sep=""))

first <- T
for(split in seq(1, max.splits.RS)){
  split.col <- vector()
  for(window.no in seq(1, length(tree.files))){
    for(id in ids){
      split.col <- c(split.col, read.proportions.RS[[window.no]][[id]][split])
    }
  }
  if(first){
    first <- F
    RS.splits.props <- split.col
  } else {
    RS.splits.props <- cbind(RS.splits.props, split.col)
  }
}

RS.splits.props <- data.frame(RS.splits.props)
RS.splits.props$window <- windows.col
RS.splits.props$patient <- ids.col
RS.splits.props <- RS.splits.props[,c(ncol(RS.splits.props)-1, ncol(RS.splits.props), seq(1, ncol(RS.splits.props)-2))]

colnames(RS.splits.props) <- c("window", "patient", paste("prop.gp.RS.",seq(1,max.splits.RS),sep=""))

for (i in seq(1, length(ids))) {
  patient <- sort(ids)[i]
  
  this.pat.stats <- pat.stats[which(pat.stats$patient==patient),]
  
  this.pat.stats <- this.pat.stats[order(this.pat.stats$window),]
  this.pat.stats.temp <- this.pat.stats[which(this.pat.stats$reads>0),]
  if(nrow(this.pat.stats.temp)>0){
    
    #Get the zero runs for the grey rectangles
    
    last.window.zero <- F
    previous.zero <- NA
    starts <- vector()
    ends <- vector()
    
    for(item in seq(1, nrow(this.pat.stats))){
      if(!last.window.zero & this.pat.stats$reads[item] == 0){
        starts <- c(starts, this.pat.stats$window[item])
        
        last.window.zero <- T
        
      } else if(last.window.zero & this.pat.stats$reads[item] > 0){
        ends <- c(ends, this.pat.stats$window[item-1])
        
        last.window.zero <- F
      }
    }
    if(last.window.zero){
      ends <- c(ends, this.pat.stats$window[nrow(this.pat.stats)])
    }
    
    rectangles <- data.frame(start = starts, end = ends)
    
    #graph 1: read count and tip count
    
    this.pat.stats.1col <- melt(this.pat.stats.temp[1:4], id=c("patient", "window"))
    
    graph.1 <- ggplot(this.pat.stats.1col, aes(x=window, y=value, col=variable))
    
    graph.1 <- graph.1 + geom_point() +
      theme_bw() + 
      scale_y_log10(breaks=c(1,10,100,1000,10000)) +
      ylab("Count") +
      xlab("Window start") +
      scale_x_continuous(limits=c(675, 9175)) +
      scale_color_discrete(name="Variable", labels=c("Leaves", "Reads")) + 
      theme(text = element_text(size=8))
    
    y.limits <- ggplot_build(graph.1)$panel$ranges[[1]]$y.range
    
    if(nrow(rectangles)>0){
      for(rect.no in seq(1, nrow(rectangles))){
        graph.1 <- graph.1 + annotate("rect", xmin=(rectangles$start[rect.no] - (250/2)), xmax = (rectangles$end[rect.no] + (250/2)),
                                      ymin=10^(y.limits[1]), ymax=10^(y.limits[2]), fill="grey", alpha=0.5)
      }
    }
    
    
    #graph 2: subtree counts (Conservative, RS)
    
    this.pat.stats.1col <- melt(this.pat.stats.temp[c(1,2,5,6,7)], id=c("patient", "window"))
    
    graph.2a <- ggplot(this.pat.stats.1col[which(this.pat.stats.1col$variable!="subtrees.RS"),], aes(x=window, y=value))
    
    graph.2a <- graph.2a +
      geom_point(aes(shape=variable, size=variable)) +
      aes(col = variable) +
      theme_bw() + 
      ylab("Count") +
      xlab("Window start") +
      scale_x_continuous(limits=c(675, 9175)) +
      scale_shape_manual(values=c(1,19), name="Variable", labels=c("Subtrees", "Clades")) +  
      scale_size_manual(values=c(2,1), name="Variable", labels=c("Subtrees", "Clades")) +
      expand_limits(y=0) + 
      scale_color_discrete(name="Variable", labels=c("Subtrees", "Clades")) + 
      theme(text = element_text(size=8))
    
    y.limits <- ggplot_build(graph.2a)$panel$ranges[[1]]$y.range
    
    if(nrow(rectangles)>0){
      for(rect.no in seq(1, nrow(rectangles))){
        graph.2a <- graph.2a + annotate("rect", xmin=(rectangles$start[rect.no] - (250/2)), xmax = (rectangles$end[rect.no] + (250/2)),
                                      ymin=y.limits[1], ymax=y.limits[2], fill="grey", alpha=0.5)
      }
    }
    
    graph.2b <- ggplot(this.pat.stats.1col[which(this.pat.stats.1col$variable!="subtrees.cons"),], aes(x=window, y=value))
    
    graph.2b <- graph.2b +
      geom_point(aes(shape=variable, size=variable)) +
      aes(col = variable) +
      theme_bw() + 
      ylab("Count") +
      xlab("Window start") +
      scale_x_continuous(limits=c(675, 9175)) +
      expand_limits(y=0) + 
      scale_color_discrete(name="Variable", labels=c("Subtrees", "Clades")) + 
      scale_shape_manual(values=c(1,19), name="Variable", labels=c("Subtrees", "Clades")) +  
      scale_size_manual(values=c(2,1), name="Variable", labels=c("Subtrees", "Clades")) +
      theme(text = element_text(size=8))
    
    y.limits <- ggplot_build(graph.2b)$panel$ranges[[1]]$y.range
    
    if(nrow(rectangles)>0){
      for(rect.no in seq(1, nrow(rectangles))){
        graph.2b <- graph.2b + annotate("rect", xmin=(rectangles$start[rect.no] - (250/2)), xmax = (rectangles$end[rect.no] + (250/2)),
                                        ymin=y.limits[1], ymax=y.limits[2], fill="grey", alpha=0.5)
      }
    }
    
    
    this.pat.stats.1col <- melt(this.pat.stats.temp[c(1,2,8,9,10)], id=c("patient", "window"))
  
    graph.3a <- ggplot(this.pat.stats.1col[which(this.pat.stats.1col$variable!="largest.rtt.RS"),], aes(x=window, y=value))
    
    graph.3a <- graph.3a +
      geom_point(aes(shape=variable, size=variable)) +
      aes(col = variable) +
      theme_bw() + 
      ylab("Root-to-tip-distance") +
      xlab("Window start") +
      scale_x_continuous(limits=c(675, 9175)) +
      expand_limits(y=0) + 
      scale_color_discrete(name="Tip set", labels=c("All", "Largest subtree")) + 
      scale_shape_manual(values=c(1,19), name="Tip set", labels=c("All", "Largest subtree")) +
      scale_size_manual(values=c(2,1), name="Tip set", labels=c("All", "Largest subtree")) +
      theme(text = element_text(size=8))
    
    y.limits <- ggplot_build(graph.3a)$panel$ranges[[1]]$y.range
    
    if(nrow(rectangles)>0){
      for(rect.no in seq(1, nrow(rectangles))){
        graph.3a <- graph.3a + annotate("rect", xmin=(rectangles$start[rect.no] - (250/2)), xmax = (rectangles$end[rect.no] + (250/2)),
                                      ymin=y.limits[1], ymax=y.limits[2], fill="grey", alpha=0.5)
      }
    }
    
    graph.3b <- ggplot(this.pat.stats.1col[which(this.pat.stats.1col$variable!="largest.rtt.cons"),], aes(x=window, y=value))
    
    graph.3b <- graph.3b +
      geom_point(aes(shape=variable, size=variable)) +
      aes(col = variable) +
      theme_bw() + 
      ylab("Root-to-tip-distance") +
      xlab("Window start") +
      scale_x_continuous(limits=c(675, 9175)) +
      expand_limits(y=0) + 
      scale_color_discrete(name="Tip set", labels=c("All", "Largest subtree")) + 
      scale_shape_manual(values=c(1,19), name="Tip set", labels=c("All", "Largest subtree")) +
      scale_size_manual(values=c(2,1), name="Tip set", labels=c("All", "Largest subtree")) +
      theme(text = element_text(size=8))
    
    y.limits <- ggplot_build(graph.3b)$panel$ranges[[1]]$y.range
    
    if(nrow(rectangles)>0){
      for(rect.no in seq(1, nrow(rectangles))){
        graph.3b <- graph.3b + annotate("rect", xmin=(rectangles$start[rect.no] - (250/2)), xmax = (rectangles$end[rect.no] + (250/2)),
                                        ymin=y.limits[1], ymax=y.limits[2], fill="grey", alpha=0.5)
      }
    }
    
    
    graph.4a <- ggplot(this.pat.stats.temp, aes(x=window, y=largest.pat.dist.cons))
    
    graph.4a <- graph.4a +
      geom_point(alpha=0.5) +
      theme_bw() + 
      ylab("Mean patristic distance in\nlargest connected subtree") +
      xlab("Window start") +
      scale_x_continuous(limits=c(675, 9175)) +
      expand_limits(y=0) +
      theme(text = element_text(size=8))
    
    y.limits <- ggplot_build(graph.4a)$panel$ranges[[1]]$y.range
    
    if(nrow(rectangles)>0){
      for(rect.no in seq(1, nrow(rectangles))){
        graph.4a <- graph.4a + annotate("rect", xmin=(rectangles$start[rect.no] - (250/2)), xmax = (rectangles$end[rect.no] + (250/2)),
                                      ymin=y.limits[1], ymax=y.limits[2], fill="grey", alpha=0.5)
      }
    }
    
    graph.4b <- ggplot(this.pat.stats.temp, aes(x=window, y=largest.pat.dist.RS))
    
    graph.4b <- graph.4b +
      geom_point(alpha=0.5) +
      theme_bw() + 
      ylab("Mean patristic distance in\nlargest connected subtree") +
      xlab("Window start") +
      scale_x_continuous(limits=c(675, 9175)) +
      expand_limits(y=0) +
      theme(text = element_text(size=8))
    
    y.limits <- ggplot_build(graph.4b)$panel$ranges[[1]]$y.range
    
    if(nrow(rectangles)>0){
      for(rect.no in seq(1, nrow(rectangles))){
        graph.4b <- graph.4b + annotate("rect", xmin=(rectangles$start[rect.no] - (250/2)), xmax = (rectangles$end[rect.no] + (250/2)),
                                        ymin=y.limits[1], ymax=y.limits[2], fill="grey", alpha=0.5)
      }
    }
    
    cons.splits.props.1col <- melt(cons.splits.props[which(cons.splits.props$patient==patient),], id=c("patient", "window"))
    cons.splits.props.1col <- cons.splits.props.1col[!is.na(cons.splits.props.1col$value),]
    
    cons.splits.props.1col$ngroup <- sapply(cons.splits.props.1col$variable, function(x) which(levels(cons.splits.props.1col$variable)==x))
    
    colourCount = length(unique(cons.splits.props.1col$ngroup))
    getPalette = colorRampPalette(brewer.pal(9, "Greens"))
    
    graph.5a <- ggplot(cons.splits.props.1col, aes(x=window, weight=value, fill=as.factor(ngroup)))
    
    graph.5a <- graph.5a +
      geom_bar(width=200, colour="black", size=0.25) +
      theme_bw() + 
      ylab("Proportion of reads\nin discrete subtrees") +
      xlab("Window start") +
      scale_fill_manual(values = rev(getPalette(colourCount))) +
      theme(text = element_text(size=8)) + 
      guides(fill = guide_legend(title = "Subtree rank\n(by tip count)", keywidth = 1, keyheight = 0.4))
    
    y.limits <- ggplot_build(graph.5a)$panel$ranges[[1]]$y.range
    
    if(nrow(rectangles)>0){
      for(rect.no in seq(1, nrow(rectangles))){
        graph.5a <- graph.5a + annotate("rect", xmin=(rectangles$start[rect.no] - (250/2)), xmax = (rectangles$end[rect.no] + (250/2)),
                                        ymin=y.limits[1], ymax=y.limits[2], fill="grey", alpha=0.5)
      }
    }
    
    
    
    RS.splits.props.1col <- melt(RS.splits.props[which(RS.splits.props$patient==patient),], id=c("patient", "window"))
    RS.splits.props.1col <- RS.splits.props.1col[!is.na(RS.splits.props.1col$value),]
    
    RS.splits.props.1col$ngroup <- sapply(RS.splits.props.1col$variable, function(x) which(levels(RS.splits.props.1col$variable)==x))
    
    
    colourCount = length(unique(RS.splits.props.1col$ngroup))
    getPalette = colorRampPalette(brewer.pal(9, "Greens"))
    
    graph.5b <- ggplot(RS.splits.props.1col, aes(x=window, weight=value, fill=as.factor(ngroup)))
    
    graph.5b <- graph.5b +
      geom_bar(width=200, colour="black", size=0.25) +
      theme_bw() + 
      ylab("Proportion of reads\nin discrete subtrees") +
      xlab("Window start") +
      scale_fill_manual(values = rev(getPalette(colourCount))) +
      theme(text = element_text(size=8)) + 
      guides(fill = guide_legend(title = "Subtree rank\n(by tip count)", keywidth = 1, keyheight = 0.4))
    
    y.limits <- ggplot_build(graph.5b)$panel$ranges[[1]]$y.range
    
    if(nrow(rectangles)>0){
      for(rect.no in seq(1, nrow(rectangles))){
        graph.5b <- graph.5b + annotate("rect", xmin=(rectangles$start[rect.no] - (250/2)), xmax = (rectangles$end[rect.no] + (250/2)),
                                      ymin=y.limits[1], ymax=y.limits[2], fill="grey", alpha=0.5)
      }
    }
    
    
    graph.6 <- ggplot(this.pat.stats.temp, aes(x=window, y=patristic.variance))
    
    graph.6 <- graph.6 +
      geom_point(alpha=0.5) +
      theme_bw() + 
      ylab("Normalised variance of\npatristic distances between tips") +
      xlab("Window start") +
      scale_x_continuous(limits=c(675, 9175)) +
      expand_limits(y=0) +
      theme(text = element_text(size=8))
    
    y.limits <- ggplot_build(graph.6)$panel$ranges[[1]]$y.range
    
    if(nrow(rectangles)>0){
      for(rect.no in seq(1, nrow(rectangles))){
        graph.6 <- graph.6 + annotate("rect", xmin=(rectangles$start[rect.no] - (250/2)), xmax = (rectangles$end[rect.no] + (250/2)),
                                      ymin=y.limits[1], ymax=y.limits[2], fill="grey", alpha=0.5)
      }
    }
    
    plots1 <- list()
    plots1$main <- textGrob(patient,gp=gpar(fontsize=20))
    plots1 <- c(plots1, AlignPlots(graph.1, graph.2a, graph.3a, graph.4a, graph.5a, graph.6))
    plots1$ncol <- 1
    plots1$heights <- unit(c(0.25, rep(1,6)), "null")
    
    hi <- do.call(grid.arrange, plots1)
    
    ggsave(paste(output.dir, "/ConsClass/sumStats_",patient,".pdf", sep=""), hi, width=210, height=297, units="mm")
    
    plots2 <- list()
    plots2$main <- textGrob(patient,gp=gpar(fontsize=20))
    plots2 <- c(plots2, AlignPlots(graph.1, graph.2b, graph.3b, graph.4b, graph.5b, graph.6))
    plots2$ncol <- 1
    plots2$heights <- unit(c(0.25, rep(1,6)), "null")
    
    hi <- do.call(grid.arrange, plots2)
    
    ggsave(paste(output.dir, "/RSClass/sumStats_",patient,".pdf", sep=""), hi, width=210, height=297, units="mm")
    
  }
}
