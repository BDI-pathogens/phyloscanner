#!/usr/bin/env Rscript

script.dir <- "/Users/twoseventwo/Documents/phylotypes/"
source(file.path(script.dir, "TransmissionUtilityFunctions.R"))
source(file.path(script.dir, "SummariseTrees_funcs.R"))


command.line <- T
if (command.line) {
  require(optparse)
  
  option_list = list(
    make_option(c("-i", "--idFile"), type="character", default=NULL, 
                help="The path of a file containing a list of all the patient IDs to be summarised"),
    make_option(c("-o", "--outputRoot"), type="character", default=NULL, 
                help="Root string for output file names."),
    make_option("--treeFileRoot", type="character", default=NULL, 
                help="A string which every tree file begins with."),
    make_option("--treeFileList", type="character", default=NULL, 
                help="A list of tree files (single string, with file names separated by with colons). Ignored if
                --treeFileRoot is present"),
    make_option("--blacklistFileRoots", type="character", default=NULL, 
                help="One or more strings which files listing blacklisted tips begin with, separated by 
                colons."),
    make_option("--blacklistFileLists", type="character", default=NULL, 
                help="A list of lists of blacklist CSV files; each top level list corresponds to a tree file. Top 
                level lists are separated by colons, bottom level by semicolons. For example, in the string 
                'A.csv;B.csv:C.csv;D.csv', A.csv and B.csv are blacklist files for the first tree file, and 
                C.csv and D.csv are for the second. Ignored if -br is present."),
    make_option("--splitsFileRoot", type="character", default=NULL, 
                help="A string which every splits file begins with."),
    make_option("--splitsFileList", type="character", default=NULL, 
                help="A list of splits files (single string, with file names separated by with colons). Ignored if 
                -sr is present"),
    make_option(c("-w", "--windows"), type="character", default=NULL, 
                help="The window in the genome which each tree file, in alphabetical order by file name, 
                is constructed from.In each window this is given as n-m where n is the start and m the end; 
                coordinates for each window are separated by a colon. If not given, the script will attempt 
                to obtain these from the file names."),
    make_option(c("-x", "--tipRegex"), type="character", default=NULL, 
                help="Regular expression identifying tips from the dataset. Three groups: patient ID,
                read ID, and read count."),
    make_option(c("-r", "--refSeqName"), type="character", default=NULL, 
                help="Reference sequence label (if unspecified, tree will be assumed to be already rooted)")
  )
  
  # Read in the arguments

  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)

  id.file <- opt$idFile
  root.name <- opt$refSeqName
  tip.regex <- opt$tipRegex
  tree.file.root <- opt$treeFileRoot
  splits.file.root <- opt$splitsFileRoot
  blacklist.file.roots.raw <- opt$blacklistFileRoots
  output.root <- opt$outputRoot
  
  if(!is.null(tree.file.root)){
    tree.files <- sort(list.files(getwd(), pattern=paste(tree.file.root,".*\\.tree",sep="")))
  } else {
    tree.files <- unlist(strsplit(opt$treeFileList, ":"))
  }
  if(!is.null(splits.file.root)){
    splits.files <- sort(list.files(getwd(), pattern=paste(splits.file.root,".*\\.csv",sep="")))
  } else {
    splits.files <- unlist(strsplit(opt$splitsFileList, ":"))
  }
  
  if(length(tree.files)!=length(splits.files)){
    stop("Number of tree files and number of splits files differ")
  }
  
  blacklist.roots <- unlist(strsplit(blacklist.file.roots.raw, ":"))
  
  blacklist.files <- list()
  
  if(!is.null(blacklist.file.roots.raw)){
    for(root in blacklist.roots){
      blacklist.files[[root]] <- sort(list.files(getwd(), pattern=paste(root,".*\\.csv",sep="")))
      if(length(tree.files)!=length(blacklist.files[[root]])){
        stop("Number of tree files and number of blacklist files differ")
      }
    }
  } else {
    blacklist.file.lists.raw <- opt$blacklistFileLists
    groups <- unlist(strsplit(blacklist.file.lists.raw, ":"))
    if(length(groups)!=length(tree.files)){
      stop("Number of tree files and number of blacklist file groups differ")
    }
    groups.split <- strsplit(groups, ";")
    if(length(unique(lengths(groups.split)))!=1){
      stop("Variable numbers of blacklist files for different windows (please use blank files if necessary)")
    }
    for(i in seq(1, unique(lengths(groups.split)))){
      blacklist.files[[i]] <- unlist(lapply(groups.split, `[[`, i))
    }
  }
  
  
  if(!is.null(opt$windows)){
    windows <- unlist(strsplit(opt$windows, ":"))
    if(length(tree.files)!=length(windows)){
      stop("Number of tree files and number of windows differ")
    }
  } else {
    windows <- NULL
  }
} else {
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160517_clean/")
  id.file <- "patientIDList.txt"
  root.name<- "C.BW.00.00BW07621.AF443088"
  tip.regex <- "^(BEE[0-9][0-9][0-9][0-9])-[0-9].*_read_([0-9]+)_count_([0-9]+)$"
  tree.files <- sort(list.files(getwd(), pattern="RAxML_bestTree.InWindow_.*\\.tree"))
  splits.files <- sort(list.files(getwd(), pattern="Subtrees_c_run20160517_.*\\.csv"))
  blacklist.files <- list()
  blacklist.files$ReadBlacklist_ <- sort(list.files(getwd(), pattern="ReadBlacklist_.*\\.csv"))
  blacklist.files$PatientBlacklist_ <- sort(list.files(getwd(), pattern="PatientBlacklist_.*\\.csv"))
  windows <- NULL
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

add.no.data.rectangles <- function(graph, rectangle.coords, log = F){
  y.limits <- ggplot_build(graph)$panel$ranges[[1]]$y.range
  
  if(nrow(rectangle.coords)>0){
    for(rect.no in seq(1, nrow(rectangles))){
      if(log){
        graph <- graph + annotate("rect", xmin=rectangle.coords$start[rect.no], xmax = rectangle.coords$end[rect.no], 
                                ymin=10^(y.limits[1]), ymax=10^(y.limits[2]), fill="grey", alpha=0.5)
      } else {
        graph <- graph + annotate("rect", xmin=rectangle.coords$start[rect.no], xmax = rectangle.coords$end[rect.no], 
                                  ymin=y.limits[1], ymax=y.limits[2], fill="grey", alpha=0.5)
      }
    }
  }
  
  return(graph)
  
}

unfactorDataFrame <- function(x) {
  x <- data.frame(lapply(x, as.character), stringsAsFactors = F)
  x <- data.frame(lapply(x, function(y) type.convert(y, as.is=T)),
                  stringsAsFactors = F)
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
window.start.col <- vector()
window.end.col <- vector()
window.middle.col <- vector()
leaves.col <- vector()
reads.col <- vector()
subtrees.counts.col <- vector()
clades.counts.col <- vector()
overall.rtt.col <- vector()
largest.rtt.col <- vector()
largest.pat.dist.col <- vector()
patristic.variance.col <- vector()

read.proportions <- list()

max.splits <- 0

for(window.no in seq(1, length(tree.files))){

  tree.file.name <- tree.files[window.no]
  blacklist.file.names <- unlist(lapply(blacklist.files, `[[`, window.no))
  splits.file.name <- splits.files[window.no]
  
  if(!is.null(windows)){
    window <- windows[window.no]
    start <- as.numeric(unlist(strsplit(window, "-"))[1])
    end <- as.numeric(unlist(strsplit(window, "-"))[2])
    middle <- (start+end)/2
  } else {
    start <- as.numeric(strsplit(tree.file.name, "[_\\.]")[[1]][4])
    end <- as.numeric(strsplit(tree.file.name, "[_\\.]")[[1]][6])
    middle <- (start+end)/2
  }
  
  # Read the tree
  
  tree <- read.tree(file=tree.file.name)
  
  # Root the tree
  if(!is.null(root.name)){
    tree<-root(phy = tree,outgroup = root.name)
  }
  num.tips <- length(tree$tip.label)
  
  # Load the blacklists
  
  blacklist <- vector()
  
  for(blacklist.file.name in blacklist.file.names){
    if(file.exists(blacklist.file.name)){
      blacklisted.tips <- read.table(blacklist.file.name, sep=",", stringsAsFactors = F, col.names="read")
      if(length(blacklisted.tips)>0){
        blacklist <- c(blacklist, sapply(blacklisted.tips, get.tip.no, tree=tree))
      }
    } else {
      warning(paste("File ",blacklist.file.name," does not exist; skipping.",paste=""))
    }
  }
  
  # Load the splits
  
  splits.table <- read.table(splits.file.name, sep=",", stringsAsFactors = F, header=T)

  splits.table$reads <- as.numeric(unlist(strsplit(splits.table$tip.names,"count_"))[seq(2,2*nrow(splits.table),2)])

  # Find clades
  
  dummy <- resolveTreeIntoPatientClades(tree, ids, tip.regex, blacklist)
  clade.mrcas.by.patient <- dummy$clade.mrcas.by.patient
  all.clades.by.patient <- dummy$clades.by.patient
  
  # Associate each tip with its patient.
  patient.tips <- list()
  for (id in ids) patient.tips[[id]] <- vector()
  for (tip.no in seq(1, length(tree$tip.label))) {
    tip.label <- tree$tip.label[tip.no]
    id <- patient.from.label(tip.label, tip.regex)
    if (id %in% ids & !(tip.no %in% blacklist)) {
      patient.tips[[id]] <- c(patient.tips[[id]], tip.label)
    }
  }
  
  ids.col <- c(ids.col, ids)
  window.start.col <- c(window.start.col, rep(start, num.ids))
  window.end.col <- c(window.end.col, rep(end, num.ids))
  window.middle.col <- c(window.middle.col, rep(middle, num.ids))
  
  read.proportions.this.window <- list()
  
  read.proportions.this.window$window.middle <- middle
  
  # Get each statistic for each patient
  
  for (i in 1:num.ids) {
#    cat("Calculating statistics in window",start, " for patient ", ids[i], "\n", sep="")
    id <- ids[i]

    num.leaves <- length(patient.tips[[id]])
    
    leaves.col <- c(leaves.col, num.leaves)
    
    reads.per.tip <- vector()
    
    for (tip.label in patient.tips[[id]]) reads.per.tip <- c(reads.per.tip, as.numeric(sub(tip.regex, "\\3", tip.label)))
    
    num.reads <- sum(reads.per.tip)
    
    reads.col <- c(reads.col, num.reads)
    
    num.subtrees <- length(unique(splits.table[which(splits.table$orig.patients==id),]$patient.splits))
   
    subtrees.counts.col <- c(subtrees.counts.col, num.subtrees)
    
    num.clades <- length(all.clades.by.patient[[id]])
    clades.counts.col <- c(clades.counts.col, num.clades)
    
    all.tips <- patient.tips[[id]]

    if(length(all.tips)>1){
      
      subtree.all <- drop.tip(tree, tip=tree$tip.label[!(tree$tip.label %in% all.tips)])
      
      # just in case pruning changes the tip order
      
      reads.per.tip <- vector()
      
      for (tip.label in subtree.all$tip.label) reads.per.tip <- c(reads.per.tip, as.numeric(sub(tip.regex, "\\3", tip.label)))
      
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
    
    if(num.subtrees <= 1){
      largest.rtt <- overall.rtt
      if(num.subtrees==0){
        mean.pat <- NA
      } else if(num.subtrees==1) {
        if(length(all.tips)==1){
          mean.pat <- NA
        } else {
          pat.distances <- cophenetic(subtree.all)
          mean.pat <- mean(pat.distances[upper.tri(pat.distances)])
        }
      }
    } else {
      relevant.reads <- splits.table[which(splits.table$orig.patients==id),]
      splits <- unique(relevant.reads$patient.splits)
      reads.per.split <- sapply(splits, function(x) sum(splits.table[which(splits.table$patient.splits==x),]$reads ) )
      winner <- splits[which(reads.per.split==max(reads.per.split))]
      if(length(winner)>1){
        cat("Patient ",id," has a joint winner in window ", middle, "\n", sep="" )
        #need to talk about this...
        winner <- winner[1]
      }
      winner.tips <- relevant.reads[which(relevant.reads$patient.splits==winner),]$tip.names
      if(length(winner.tips)==1){
        largest.rtt <- 0
        mean.pat <- NA
      } else {
        subtree <- drop.tip(tree, tip=tree$tip.label[!(tree$tip.label %in% winner.tips)])
        
        # just in case pruning changes the tip order
        
        reads.per.tip <- vector()
        
        for (tip in subtree$tip.label) reads.per.tip <- c(reads.per.tip, as.numeric(unlist(strsplit(tip,"count_"))[2]))
        
        names(reads.per.tip) <- subtree$tip.label
        
        largest.rtt <- calcMeanRootToTip(subtree, reads.per.tip)
        pat.distances <- cophenetic(subtree)
        mean.pat <- mean(pat.distances[upper.tri(pat.distances)])
      }
    }
    
    largest.rtt.col <- c(largest.rtt.col, largest.rtt)
    largest.pat.dist.col <- c(largest.pat.dist.col, mean.pat)
    
    this.pat.splits <- splits.table[which(splits.table$orig.patients==id),]
    
    if(nrow(this.pat.splits)>0){
      this.pat.reads.by.split <- aggregate(this.pat.splits$reads, by=list(Category=this.pat.splits$patient.splits), sum)
      
      if(nrow(this.pat.reads.by.split) > max.splits){
        max.splits <- nrow(this.pat.reads.by.split>max.splits)
      }
      
      this.pat.reads.by.split <- this.pat.reads.by.split[order(this.pat.reads.by.split$x, decreasing = T),]
      this.pat.reads.by.split$proportion <- this.pat.reads.by.split$x/sum(this.pat.reads.by.split$x)
      read.proportions.this.window[[id]] <- this.pat.reads.by.split$proportion
    } else {
      read.proportions.this.window[[id]] <- NA
    }
  }
  read.proportions[[window.no]] <- read.proportions.this.window
}

pat.stats <- data.frame(window.start = window.start.col, window.middle = window.middle.col,
                        window.end = window.end.col, patient = ids.col, leaves = leaves.col, 
                        reads = reads.col, subtrees = subtrees.counts.col, clades = clades.counts.col,
                        overall.rtt = overall.rtt.col, largest.rtt = largest.rtt.col,
                        largest.pat.dist = largest.pat.dist.col, patristic.variance = patristic.variance.col)

first <- T
for(split in seq(1, max.splits)){
  split.col <- vector()
  for(window.no in seq(1, length(tree.files))){
    for(id in ids){
      split.col <- c(split.col, read.proportions[[window.no]][[id]][split])
    }
  }
  if(first){
    first <- F
    splits.props <- split.col
  } else {
    splits.props <- cbind(splits.props, split.col)
  }
}

splits.props <- data.frame(splits.props)
splits.props$window.middle <- window.middle.col
splits.props$patient <- ids.col
splits.props <- splits.props[,c(ncol(splits.props)-1, ncol(splits.props), seq(1, ncol(splits.props)-2))]

colnames(splits.props) <- c("window.middle", "patient", paste("prop.gp.",seq(1,max.splits),sep=""))

ews <- min(pat.stats$window.start)
lwe <- max(pat.stats$window.end)

pat.stats$prop.reads.largest.subtree <- splits.props$prop.gp.1

# todo - the subtree read proportions aren't here, but how relevant are the means of those anyway?

mean.na.rm <- function(x) mean(x, na.rm = T)

pat.stats <- unfactorDataFrame(pat.stats)
by.patient <- pat.stats %>% group_by(patient)
pat.stats.summary <-
  as.data.frame(by.patient %>% summarise_each(funs(mean.na.rm)))

write.csv(pat.stats.summary, file.path(paste(output.root,"_patStatsSummary.csv",sep="")))

pdf(paste(output.root,"_plots.pdf",sep=""), width=8.26772, height=11.6929)

for (i in seq(1, length(ids))) {
  
  patient <- sort(ids)[i]
#  cat("Drawing graphs for patient ",patient,"\n", sep="")
  
  this.pat.stats <- pat.stats[which(pat.stats$patient==patient),]
  
  this.pat.stats <- this.pat.stats[order(this.pat.stats$window.middle),]
  this.pat.stats.temp <- this.pat.stats[which(this.pat.stats$reads>0),]
  if(nrow(this.pat.stats.temp)>0){
    
    #Get the zero runs for the grey rectangles
    
    last.window.zero <- F
    previous.zero <- NA
    starts <- vector()
    ends <- vector()
    
    display.starts <- this.pat.stats$window.start + (0.3*(this.pat.stats$window.middle - this.pat.stats$window.start))
    display.ends <- this.pat.stats$window.end - (0.3*(this.pat.stats$window.end - this.pat.stats$window.middle))
    
    for(item in seq(1, nrow(this.pat.stats))){
      if(!last.window.zero & this.pat.stats$reads[item] == 0){
        starts <- c(starts, display.starts[item])
        
        last.window.zero <- T
        
      } else if(last.window.zero & this.pat.stats$reads[item] > 0){
        ends <- c(ends, display.ends[item-1])
        
        last.window.zero <- F
      }
    }
    if(last.window.zero){
      ends <- c(ends, display.ends[nrow(this.pat.stats)])
    }
    
    rectangles <- data.frame(start = starts, end = ends)
    
    #graph 1: read count and tip count
    
    this.pat.stats.1col <- melt(this.pat.stats.temp[c("window.middle","patient","leaves","reads")], id=c("patient", "window.middle"))
    
    graph.1 <- ggplot(this.pat.stats.1col, aes(x=window.middle, y=value, col=variable))
    
    graph.1 <- graph.1 + geom_point() +
      theme_bw() + 
      scale_y_log10(breaks=c(1,10,100,1000,10000)) +
      ylab("Count") +
      xlab("Window centre") +
      scale_x_continuous(limits=c(ews, lwe)) +
      scale_color_discrete(name="Variable", labels=c("Leaves", "Reads")) + 
      theme(text = element_text(size=8))
    
    graph.1 <- add.no.data.rectangles(graph.1, rectangles, TRUE)
    
    #graph 2: subtree counts (Conservative, RS)
    
    this.pat.stats.1col <- melt(this.pat.stats.temp[c("window.middle", "patient", "subtrees", "clades")], id=c("patient", "window.middle"))
    
    graph.2 <- ggplot(this.pat.stats.1col, aes(x=window.middle, y=value))
    
    graph.2 <- graph.2 +
      geom_point(aes(shape=variable, size=variable)) +
      aes(col = variable) +
      theme_bw() + 
      ylab("Count") +
      xlab("Window centre") +
      scale_x_continuous(limits=c(ews, lwe)) +
      scale_shape_manual(values=c(1,19), name="Variable", labels=c("Subtrees", "Clades")) +  
      scale_size_manual(values=c(2,1), name="Variable", labels=c("Subtrees", "Clades")) +
      expand_limits(y=0) + 
      scale_color_discrete(name="Variable", labels=c("Subtrees", "Clades")) + 
      theme(text = element_text(size=8))
    
    graph.2 <- add.no.data.rectangles(graph.2, rectangles)
    
    this.pat.stats.1col <- melt(this.pat.stats.temp[c("window.middle","patient","overall.rtt","largest.rtt")], id=c("patient", "window.middle"))
  
    graph.3 <- ggplot(this.pat.stats.1col, aes(x=window.middle, y=value))
    
    graph.3 <- graph.3 +
      geom_point(aes(shape=variable, size=variable)) +
      aes(col = variable) +
      theme_bw() + 
      ylab("Root-to-tip-distance") +
      xlab("Window centre") +
      scale_x_continuous(limits=c(ews, lwe)) +
      expand_limits(y=0) + 
      scale_color_discrete(name="Tip set", labels=c("All", "Largest subtree")) + 
      scale_shape_manual(values=c(1,19), name="Tip set", labels=c("All", "Largest subtree")) +
      scale_size_manual(values=c(2,1), name="Tip set", labels=c("All", "Largest subtree")) +
      theme(text = element_text(size=8))
    
    graph.3 <- add.no.data.rectangles(graph.3, rectangles)
    
    graph.4 <- ggplot(this.pat.stats.temp, aes(x=window.middle, y=largest.pat.dist))
    
    graph.4 <- graph.4 +
      geom_point(alpha=0.5) +
      theme_bw() + 
      ylab("Mean patristic distance in\nlargest connected subtree") +
      xlab("Window centre") +
      scale_x_continuous(limits=c(ews, lwe)) +
      expand_limits(y=0) +
      theme(text = element_text(size=8))
    
    graph.4 <- add.no.data.rectangles(graph.4, rectangles)
    
    splits.props.1col <- melt(splits.props[which(splits.props$patient==patient),], id=c("patient", "window.middle"))
    splits.props.1col <- splits.props.1col[!is.na(splits.props.1col$value),]
    
    splits.props.1col$ngroup <- sapply(splits.props.1col$variable, function(x) which(levels(splits.props.1col$variable)==x))
    
    colourCount = length(unique(splits.props.1col$ngroup))
    getPalette = colorRampPalette(brewer.pal(9, "Greens"))
    
    graph.5 <- ggplot(splits.props.1col, aes(x=window.middle, weight=value, fill=as.factor(ngroup)))
    
    graph.5 <- graph.5 +
      geom_bar(width=200, colour="black", size=0.25) +
      theme_bw() + 
      ylab("Proportion of reads\nin discrete subtrees") +
      xlab("Window centre") +
      scale_x_continuous(limits=c(ews, lwe)) +
      scale_fill_manual(values = rev(getPalette(colourCount))) +
      theme(text = element_text(size=8)) + 
      guides(fill = guide_legend(title = "Subtree rank\n(by tip count)", keywidth = 1, keyheight = 0.4))
    
    graph.5 <- add.no.data.rectangles(graph.5, rectangles)
    
    
    graph.6 <- ggplot(this.pat.stats.temp, aes(x=window.middle, y=patristic.variance))
    
    graph.6 <- graph.6 +
      geom_point(alpha=0.5) +
      theme_bw() + 
      ylab("Normalised variance of\npatristic distances between tips") +
      xlab("Window centre") +
      scale_x_continuous(limits=c(ews, lwe)) +
      expand_limits(y=0) +
      theme(text = element_text(size=8))
    
    graph.6 <- add.no.data.rectangles(graph.6, rectangles)
    
    plots1 <- list()
    plots1$main <- textGrob(patient,gp=gpar(fontsize=20))
    plots1 <- c(plots1, AlignPlots(graph.1, graph.2, graph.3, graph.4, graph.5, graph.6))
    plots1$ncol <- 1
    plots1$heights <- unit(c(0.25, rep(1,6)), "null")
    
    do.call(grid.arrange, plots1)
    
#    ggsave(paste("sumStats_",patient,".pdf", sep=""), hi, width=210, height=297, units="mm")
    
  }
}

dev.off()