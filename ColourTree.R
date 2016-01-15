#!/usr/bin/env Rscript

command.line <- FALSE
if (command.line) {
  # Read in the arguments
  args <- commandArgs(TRUE)
  if (length(args) < 7) {
    cat(paste("At least 7 arguments must be specified:\n* a file containing",
    "the IDs to be coloured, one per line;\n* the string/character that",
    "follows the ID part of the tip name (i.e. those tip names to be coloured",
    "should equal one of the IDs followed by this string/character then",
    "anything at all);\n* the font size for the tip labels;\n* the line width", 
    "for the tree;\n* a seed for randomising the order of colour allocation to", 
    "the IDs file (use -1 to turn off shuffle;\n* the directory where output",
    "tree pdfs will be produced;\n* finally, any number of tree",
    "files.\nQuitting.\n"))
    quit("no", 1)
  }
  id.file <- args[1]
  id.delimiter <- args[2]
  font.size <- as.numeric(args[3])
  line.width <- as.numeric(args[4])
  seed <- as.integer(args[5])
  output.dir <- args[6]
  tree.files <- args[-(1:6)]
  for (tree.file in tree.files) {
    if (! file.exists(tree.file)) {
      print(paste(tree.file, "does not exist. Quitting."))
      quit("no", 1)
    }
  }
} else {
  # Some options to run within R
  setwd("/Users/cfraser/Dropbox (Infectious Disease)/PROJECTS/BEEHIVE/phylotypes")
  rm(list=ls())
  #id.file <- "inputs/ListOfIDs_IDUs_IncNotOK.txt"
  id.file <- "inputs/ListOfIDs_SubtypeCclade.txt"
  id.delimiter <- "-1"
  font.size <- 0.3
  line.width <- 1
  seed <- 1
  #output.dir <- "CF_TREES/IDUs"
  output.dir <- "CF_TREES/SubtypeCclade"
  input.dir<-"PhylotypeOutput/SubtypeCclade"
  #tree.files <- list.files("PhylotypeOutput/IDUs", pattern="*.tree")
  tree.files <- list.files(input.dir, pattern="*.tree")
}

require(phytools)
require(grid)

dir.create(output.dir, showWarnings=F)
#StatsList <- list()

# Read in the IDs. Remove duplicates. Shuffle their order if desired.
ids <- scan(id.file, what="", sep="\n", quiet=TRUE)
num.ids <- length(ids)
if (num.ids == 0) {
  cat(paste("No IDs found in ", id.file, ". Quitting.\n", sep=""))
  quit("no", 1)
}
ids <- unique(ids)
if (seed != -1) {
  set.seed(seed)
  ids <- sample(ids)
}


# Define a colour for each unique ID.
# Needs to be written twice to work properly for some reason!
id.colours <- setNames(palette(rainbow(num.ids))[1:num.ids], ids)
id.colours <- setNames(palette(rainbow(num.ids))[1:num.ids], ids)

#dummy <- rep(0, num.ids)
#pat.stats <- data.frame(patient.id=ids, num.leaves=dummy, monophyletic=dummy, mean.size=dummy, coeff.of.var.size=dummy)
#pat.stats <- data.frame(patient.id=character(), window=integer(), num.leaves=integer(), monophyletic=integer(), mean.size=double(), coeff.of.var.size=double())
pat.stats <- data.frame(patient.id="test", window=0, num.leaves=1, monophyletic=1,
mean.size=1, coeff.of.var.size=1)
pat.stats$patient.id <- as.character(pat.stats$patient.id)

for (tree.file in tree.files) {
  
  tree.file.name<-paste(input.dir,tree.file,sep="/")

  tree <- read.tree(file=tree.file.name)
  num.tips <- length(tree$tip.label)
  
  # Assign each tip its colour
  this.tree.colours <- vector()
  for (i in 1:num.tips) {
    tip.label <- tree$tip.label[i]
    id <- unlist(strsplit(tip.label, id.delimiter))[1]
    if (id %in% ids) {
      this.tree.colours[i] <- id.colours[which(ids==id)]
    } else {
      this.tree.colours[i] <- "black"
    }
  }

  # Rotate internal nodes to have right daughter clades more species rich than
  # left daughter clades. 
  tree <- ladderize(tree)  

  # Plot the tree.
  # First make the automatic tip labels transparent, then replot with colour.
  out.tree.basename <- paste(basename(tree.file), ".pdf", sep="")
  out.tree <- file.path(output.dir, out.tree.basename)
  pdf(out.tree)
  opar <- par()
  par(fg="transparent")
  plotTree(tree, color="black", fsize=font.size, lwd=line.width, ylim=c(-1, num.tips))
  lastPP <- get("last_plot.phylo", env=.PlotPhyloEnv)
  par(fg="black")
  text(lastPP$xx[1:num.tips], lastPP$yy[1:num.tips], tree$tip.label,
       cex=font.size, col=this.tree.colours, pos=4, offset=0.1, font=0.1)
  dev.off()
  par <- opar

  # Associate each tip to its patient.
  patient.tips <- list()
  for (id in ids) patient.tips[[id]] <- vector()
  for (tip.label in tree$tip.label) {
    id <- unlist(strsplit(tip.label, id.delimiter))[1]
    if (id %in% ids) {
      patient.tips[[id]] <- c(patient.tips[[id]], tip.label)
    }
  }

  # This is exceptionally bad code - how else to pull out window coordinate though?
  window <- as.numeric(strsplit(tree.file, "_")[[1]][3])

  # For each patient: record whether his/her tips are monophyletic, find the
  # pairwise patristic distances between the tips - the 'cophenetic distances' -
  # and characterise those distances. 
  for (i in 1:num.ids) {
    id <- ids[i]
    num.leaves <- length(patient.tips[[id]])
    monophyletic <- as.numeric(is.monophyletic(tree, patient.tips[[id]]))
    if (num.leaves > 1) {
      subtree <- drop.tip(phy=tree,
      tip=tree$tip.label[!(tree$tip.label %in% patient.tips[[id]])])
      subtree.dist.matrix <- cophenetic(subtree)
      subtree.dist.cw <- subtree.dist.matrix[lower.tri(subtree.dist.matrix)]
      subtree.dist.cf <- as.vector(subtree.dist.matrix)
      # Distances are not weighted for the number of reads for each tip. 
      # Corrections are included because mean and variance of distance matrices
      # include zeros on diagonal, but shouldn't.
      mean.size.cf <- mean(subtree.dist.cf)/(1-1/num.leaves)
      variance.cf <- var(subtree.dist.cf)*(1+1/num.leaves)-mean.size.cf^2/num.leaves
      mean.size.cw <- mean(subtree.dist.cw)
      variance.cw <- var(subtree.dist.cw)
      #cat("CF mean & var: ", mean.size.cf, variance.cf, '\n')
      #cat("CW mean & var: ", mean.size.cw, variance.cw, '\n')
      #cat('\n')
      mean.size <- mean.size.cw
      variance <- variance.cw
      coeff.of.var.size <- ifelse(num.leaves > 2, sqrt(variance)/mean.size, 0)
    } else {
      mean.size <- 0
      coeff.of.var.size <- 0
    }
    pat.stats <- rbind(pat.stats, c(id, window, num.leaves, monophyletic,
    mean.size, coeff.of.var.size))
  }
}

pat.stats$mean.size <- as.numeric(pat.stats$mean.size)
pat.stats$coeff.of.var.size <- as.numeric(pat.stats$coeff.of.var.size)
pat.stats <- pat.stats[!pat.stats$patient.id=="test", ]
pat.stats <- pat.stats[order(pat.stats$patient.id), ]

pdf(file.path(output.dir, "Patient.Statistics.pdf"))
for (id in ids) {
  par(mfrow=c(2, 3))
  pat <- pat.stats[pat.stats$patient.id==id, ]
  plot(pat$window, pat$num.leaves, main=pat$patient.id[1],
  xlab="Genome location", ylab="Number of tips")
  plot(pat$window, pat$monophyletic, main=pat$patient.id[1],
  xlab="Genome location", ylab="Monophyletic (1=yes, 0=no)", ylim=c(0, 1))
  plot(pat$window, pat$mean.size, main=pat$patient.id[1],
  xlab="Genome location", ylab="Mean cophenetic distance",
  ylim=c(0, max(pat.stats$mean.size)))
  plot(pat$window, pat$coeff.of.var.size, main=pat$patient.id[1],
  xlab="Genome location",
  ylab="Coefficient of variation of cophenetic distances",
  ylim=c(0, max(pat.stats$coeff.of.var.size)))
}
dev.off()
