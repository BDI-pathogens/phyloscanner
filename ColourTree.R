#!/usr/bin/env Rscript

options(warn=1, error=traceback)
warnings()

# Some style points from https://google.github.io/styleguide/Rguide.xml enforced
# by Chris:
# Exactly one space either side of +, -, <-, <, >, ==, etc.
# Always a space after a comma, never before.
# Opening brackets should have a space before but not after, closing brackets
# should have a space after but not before.
# Line length <= 80 characters.

# TODO:
# 2. Update plot of patients statistics & cross window mean of patient
# statistics
# 3. Code is now slow:
# Main problem is that repeatedly generating and searching all possible
# subclades is just plain stupid
# Solution: generate subclades only once, document their size, and then
# thin them as the loop goes

# TODO: Rich Tips:
# Print things consistently - message? cat? print? stdout/stderr

# TODO: find a system-independent way to make R source a file in the same dir
# or a sub dir. (Wrap up into a package?) Answers at
# http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
# do not work.
# The directory in which this code (and therefore SummariseTrees_funcs.R) lives:
path.here <- '~/Dropbox (Infectious Disease)/PhylotypesCode/'
source(file.path(path.here, 'SummariseTrees_funcs.R'))

# We flip the value of the bool command.line to determine whether arguments are
# read in from the command line (TRUE) or set internally in the code (FALSE).
command.line <- TRUE
if (command.line) {
	args <- commandArgs(TRUE)
  checkArgs(args)
	id.file <- args[1]
	font.size <- as.numeric(args[2])
	line.width <- as.numeric(args[3])
	seed <- as.integer(args[4])
	output.dir <- args[5]
  #whose.code <- args[6]
	tree.files <- args[-(seq_len(5))]
  tree.files.basenames <- lapply(tree.files, basename)
  how.to.read.ids <- "List" # other options are deprecated to be removed later
  
} else {
  # In this scope, we set by hand the values of the variables read in from the
  # command line in the scope above. Other variables should be set below this
  # scope.
  rm(list=ls())
  # Folder locations
  root.folder <- "/Users/Christophe/Dropbox (Infectious Disease)/PROJECTS"
  nameFolder <- function(name) paste(root.folder, name, sep = "")
  setwd(nameFolder("/HIV/phylotypes"))
  input.dir <- nameFolder("/BEEHIVE/phylotypes/run20160209/RAxML_bestTree")
  output.dir <- nameFolder("/BEEHIVE/phylotypes/run20160209/R_Output")
  
  # Find the tree files. NB \\. is needed in the pattern or it ignores the dot.
  tree.files.basenames <- list.files(input.dir, pattern="*\\.tree")
  tree.files <- lapply(list.files(input.dir, pattern="*\\.tree"),
                       function(x) file.path(input.dir, x))
	font.size <- 0.15 # for the tip labels in the tree
	line.width <- 1 # in the tree
	seed <- -1 # seed for randomising the colours associated with IDs
	# no randomisation if seed <- -1

  # Deprecated, to be removed later:
  # How do we read in all patient IDs?
  # Options are: "List", "Alignment" and "Trees"
  # "List" reads in all the patient IDs from a list in file 'id.file'
  # "Alignment" reads in the patient IDs from the master alignment
  # 'alignment.file'
  # "Trees" reads IDs from the tip labels of the trees
  how.to.read.ids <- "Trees"
  #id.file <- "ListOfIDs_run20160111.txt"
  alignment.file <- nameFolder("/BEEHIVE/phylotypes/run20160209/RefsAln.fasta")
  alignment.file.tag <- "-1_ref.fasta"

}

################################################################################
# Some more user input variables:

# The tip name associated with the root of the phylogeny:
#root.name <- "Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
#root.name <- "B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
root.name <- "C.BW.00.00BW07621.AF443088"
print.trees <- T # whether or not to print separate PDF files for each tree
exit.after.colouring.trees <- F
colour.edges <- TRUE
print.per.patient.pdfs <- F
test.clade.ordering.only <- F
num.clades.for.output <- 5
# Regex for how input tree files should be named:
tree.file.regex <- "RAxML_bestTree\\.InWindow_([0-9]+)_to_([0-9]+).*"
# Regex for how tip names from patients should be named:
tip.regex <- "^(.*)_read_([0-9]+)_count_([0-9]+)$"


# Some deprecated variables, for backwards compatibility, to be removed later:
whose.code <- "Rich" # "Rich", "Christophe", "Chris" or "Chris_slow"
id.delimiter <- '_' 
################################################################################

require(phytools)
require(dplyr)
require(ggplot2)
require(plotrix)
require(RColorBrewer)
require(phangorn)

# Utlities ##################
unfactorDataFrame <- function(x) {
  x <- data.frame(lapply(x, as.character), stringsAsFactors = F)
  x <- data.frame(lapply(x, function(y) type.convert(y, as.is=T)),
                  stringsAsFactors = F)
}

#############################

if (!dir.exists(output.dir)) dir.create(output.dir)

checkTreeFiles(tree.files)
checkTreeFileNames(tree.files.basenames, tree.file.regex)

# Read in the trees.
# TODO: once the deprecated how.to.read.ids == "Trees" option is removed, move
# this block inside the iteration over trees.
num.trees <- length(tree.files)
trees <- list(length=num.trees)
for (tree.number in seq_len(num.trees)) {
  tree.file <- tree.files[[tree.number]]
  tree <- read.tree(file=tree.file)
  trees[[tree.file]] <- tree
}
 
# Deprecated, for backwards compatibility, to be removed later:
ids <- vector()
if (how.to.read.ids == "Trees") {
  for (tree.file in tree.files) {
    id.list <- sapply(trees[[tree.file]]$tip.label, function(x) strsplit(x,
                                                                         split = "_")[[1]][1])
    id.list <- sapply(id.list, function(x) strsplit(x,
                                                    split = id.delimiter)[[1]][1])
    id.list <- unique(id.list)
    ids <- c(ids, id.list)
  }
  ids <- unique(ids)
  ids <- ids[order(ids)]
  ids <- ids[grep("BEE", ids)]
  
} else if (how.to.read.ids == "Alignment") {
  # !!!!!This option seems messed up, and likely won't work!!!!
  # reads in the IDs from the list of files in the folder - unusually named
  # files might get a bit garbled.
  ids <- list.files(id.folder, pattern = "*\\.fasta")
  ids <- sapply(ids, function (x) unlist(strsplit(x,
                                                  split = id.folder.tag)[[1]][1]))
} else {
	# alternatively reads in the IDs from a specific file
	ids <- readLines(id.file)
}

# Remove duplicate IDs. Shuffle the ID order if desired (for colouring).
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

if (print.trees) {
  # Define a colour for each unique ID.
  # Needs to be written twice to work properly for some reason!
  id.colours <- setNames(palette(rainbow(num.ids))[seq_len(num.ids)], ids)
  id.colours <- setNames(palette(rainbow(num.ids))[seq_len(num.ids)], ids)
}

# Define the data frame to which summary statistics will be written.
# (The first row is later deleted.)
pat.stats <- data.frame(patient.id = "test",
		window = 0,
		num.reads.total = 1,
		num.leaves = 1,
		num.clades = 1,
		overall.root.to.tip = 0,
    # TODO: how can the ten lines below be replaced by something compatible with
    # an arbitrary num.clades.for.output (currently set to 5)?
		prop.reads.clade.1 = 0,
		prop.reads.clade.2 = 0,
		prop.reads.clade.3 = 0,
		prop.reads.clade.4 = 0,
		prop.reads.clade.5 = 0,
		root.to.tip.clade.1 = 0,
		root.to.tip.clade.2 = 0,
		root.to.tip.clade.3 = 0,
		root.to.tip.clade.4 = 0,
		root.to.tip.clade.5 = 0,
		stringsAsFactors = F)

# Iterate through the trees
for (tree.number in seq_len(num.trees)) {
  
  tree.file <- tree.files[[tree.number]]
  tree.file.basename <- tree.files.basenames[[tree.number]]
  tree <- trees[[tree.file]]
	num.tips <- length(tree$tip.label)
  message(paste('Now analysing tree', tree.file.basename, '\n'))

  # Root the tree
  if (! root.name %in% tree$tip.label) {
    cat(paste('Error: the specified root, ', root.name, ', is not in ',
    tree.file.basename, '. Continuing to the next tree.\n', sep=''))
    next
  }
	tree <- root(phy = tree, outgroup = root.name)

  # Resolve the tree into clades each consisting of only one patient:
  if (whose.code == "Rich") {
    dummy <- resolveTreeIntoPatientClades(tree, ids, tip.regex)
    clade.mrcas.by.patient <- dummy$clade.mrcas.by.patient
    all.clades.by.patient <- dummy$clades.by.patient
    #print(clade.mrcas.by.patient)
  } else if (whose.code == "Chris_slow") {
    all.clades.by.patient <- resolveTreeIntoPatientCladesChris(tree)
  }

  # Print the tree if desired
  if (print.trees) {
    colourTree(tree, tip.regex, id.colours, tree.file.basename, output.dir,
    font.size, line.width, ids, colour.edges, clade.mrcas.by.patient)
  }
  if (exit.after.colouring.trees) next

	# Associate each tip to its patient.
	patient.tips <- list()
	for (id in ids) patient.tips[[id]] <- vector()
	for (tip.label in tree$tip.label) {
		id <- sub(tip.regex, "\\1", tip.label)
		if (id %in% ids) {
			patient.tips[[id]] <- c(patient.tips[[id]], tip.label)
		}
	}
	
	# Pull out the window coordinates
	left.window.edge  <- as.integer(sub(tree.file.regex, "\\1",
  tree.file.basename))
	right.window.edge <- as.integer(sub(tree.file.regex, "\\2",
  tree.file.basename))
  window <- (left.window.edge + right.window.edge) / 2

  # Iterate through all patients
  for (id in ids) {

    if (whose.code == "Christophe") {
      summary <- analyseCladesChristophe(id, patient.tips, tree,
      test.clade.ordering.only)
    } else {

      # Find this patient's clades
      if (whose.code == "Chris") {
        clades <- findOnePatientsClades(tree, id)
      } else {
        clades <- all.clades.by.patient[[id]]
      }

      # For testing clade finding algorithms do the same thing:
      if (test.clade.ordering.only && length(clades) > 0) {
          printSortedClades(clades, tip.regex, num.clades.for.output)
      }

      summary <- summariseClades(tree, clades, num.clades.for.output, tip.regex,
      whose.code)
    }
    pat.stats <- rbind(pat.stats, c(id, window, summary$num.reads.total, 
    summary$num.tips, summary$num.clades, summary$overall.root.to.tip, 
    summary$prop.reads.per.clade, summary$root.to.tip.per.clade))
  }    

}

if (exit.after.colouring.trees) {
  print("All trees coloured; exiting successfully.")
  quit("no", 0)
}

if (test.clade.ordering.only) {
  print("Tested clade ordering; exiting successfully.")
  quit("no", 0)
}

pat.stats <- pat.stats[!pat.stats$patient.id=="test",]
# Delete the first dummy line
pat.stats <- pat.stats[order(pat.stats$patient.id),]
# TODO: replace pat.stats$num.leaves = 0 by NA

if (print.per.patient.pdfs) {
  pdf(file.path(output.dir, "Patient.Statistics.pdf"), paper = "a4")
  for (id in sort(ids)) {
    print(id)
    par(mfrow=c(3, 2))
    # putting this here also sets new page, in case there aren't six plots
    pat <- pat.stats[pat.stats$patient.id==id,]
    pat <- unfactorDataFrame(pat)
    pat <- pat[order(pat$window),]
    print(pat$num.reads.total)
    plot(pat$window, pat$num.reads.total, main=pat$patient.id[1],
         xlab="Genome location", ylab="Number of reads",
         log = "y", ylim = c(1, 1e+5))
    if (sum(pat$num.reads.total) > 0) {
      plot(pat$window, pat$num.leaves, main=pat$patient.id[1],
           xlab="Genome location", ylab="Number of tips")
      plot(pat$window, pat$num.clades, main = pat$patient.id[1],
           ylim = c(0, 10),
           xlab="Genome location", ylab="Number of clades")
      
      opar <- par(mar = c(2.1, 4.1, 2.1, 5.1) + 0.3)
      plot(pat$window, pat$prop.reads.clade.1, main = pat$patient.id[1],
           xlab="Genome location", ylab="Proportion of reads in largest clade",
           type = "b", pch = 16, cex = 0.6, ylim = c(0, 1))
      par(new = T)
      plot(pat$window, pat$overall.root.to.tip,
           type = "b", pch = 1, col = "red", ylim = c(0, 0.4), cex = 0.6,
           axes = F, xlab = NA, ylab = NA)
      axis(4, ylim = c(0, 0.4), col = "red")
      mtext(side = 4, line = 3, "Root to tip distance", cex=0.7)
      points(pat$window, pat$root.to.tip.clade.1,
             type = "l", col = "darkgreen", lty = 2)
      par(opar)
      
      plot(pat$window, pat$root.to.tip.clade.1, main = pat$patient.id[1],
           type = "b", pch = 16, cex = 0.6,
           xlab="Genome location", ylab="Root to tip distance in largest clade")
      
      ydat <- matrix(
        c(pat$prop.reads.clade.1,
          pat$prop.reads.clade.2,
          pat$prop.reads.clade.3,
          pat$prop.reads.clade.4,
          pat$prop.reads.clade.5
        ),
        nrow = length(pat$prop.reads.clade.1),
        ncol = 5
      )
      ydat <- t(ydat)
      
      barplot(ydat,
              xlab = "Genome location",
              ylab = "Proportion of reads in each of 5 largest clades",
              col = brewer.pal(5, "Blues")[5:1],
              names.arg = pat$window)
      
    }
  }
  dev.off()
}

# create the ggplot2 data, using Year and thousands (of people) and filled by
# age group.
#    clade.1 <- pat[,c("patient.id", "window", "prop.reads.clade.1")]
#    clade.2 <- pat[,c("patient.id", "window", "prop.reads.clade.2")]
#    clade.3 <- pat[,c("patient.id", "window", "prop.reads.clade.3")]
#    clade.4 <- pat[,c("patient.id", "window", "prop.reads.clade.4")]
#    clade.5 <- pat[,c("patient.id", "window", "prop.reads.clade.5")]
#    colnames(clade.1)[3] <- "prop.reads"
#    colnames(clade.2)[3] <- "prop.reads"
#    colnames(clade.3)[3] <- "prop.reads"
#    colnames(clade.4)[3] <- "prop.reads"
#    colnames(clade.5)[3] <- "prop.reads"
#    clade.1$clade <- 1
#    clade.2$clade <- 2
#    clade.3$clade <- 1
#    clade.4$clade <- 2
#    clade.5$clade <- 1
#    plot.data <- rbind(clade.1, clade.2, clade.3, clade.4, clade.5)
#    plot.data$window <- as.numeric(plot.data$window)
#    plot.data$prop.reads <- as.numeric(plot.data$prop.reads)
#    plot.data$clade <- as.factor(plot.data$clade)

#   #plot(pat$window, pat$mean.size, main=pat$patient.id[1],
#   #  xlab="Genome location", ylab="Mean cophenetic distance",
#   #  ylim=c(0, max(pat.stats$mean.size, na.rm=T)))
#   boxplot(mean.size~window, data = pat.stats[pat.stats$monophyletic==1,],
#           main=pat$patient.id[1],
#           xlab="Genome location", ylab="Mean cophenetic distance",
#           ylim=c(0, max(pat.stats$mean.size, na.rm=T)))
#   points(1:length(pat$mean.size), pat$mean.size, col="red", pch=19, cex=1)
#   significant.stars <- ifelse(pat$size.significantly.big == 1, 1, NA)
#   points(1:length(pat$mean.size), -0.01*significant.stars, col="red", pch=8,
#   cex=1)
#
#   plot(pat$window, pat$coeff.of.var.size, main=pat$patient.id[1],
#        xlab="Genome location",
#        ylab="Coefficient of variation of cophenetic distances",
#        ylim=c(0, max(pat.stats$coeff.of.var.size, na.rm=T)))
#   plot(pat$window, pat$root.to.tip, main=pat$patient.id[1],
#        xlab="Genome location",
#        ylab="Mean root to tip patristic distance",
#        ylim=c(0, max(pat.stats$root.to.tip, na.rm=T)))
#
# write.csv(pat.stats, file.path(output.dir, "Patient.Statistics.csv"))
#
mean.na.rm <- function(x) mean(x, na.rm = T)
#
# TODO: Update this list for turning pat.stats columns into numeric
# pat.stats$window.has.data <- as.numeric(pat.stats$num.leaves > 0)
# pat.stats$num.leaves <- as.numeric(as.character(pat.stats$num.leaves))
# pat.stats$num.reads <- as.numeric(as.character(pat.stats$num.reads))
# pat.stats$monophyletic <- as.numeric(as.character(pat.stats$monophyletic))
#
pat.stats <- unfactorDataFrame(pat.stats)
by.patient <- pat.stats %>% group_by(patient.id)
pat.stats.summary <-
  as.data.frame(by.patient %>% summarise_each(funs(mean.na.rm)))

write.csv(pat.stats.summary, file.path(output.dir,
                                       "Patient.Statistics.Summary.csv"))
