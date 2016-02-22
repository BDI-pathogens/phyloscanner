#!/usr/bin/env Rscript

# Some style points from https://google.github.io/styleguide/Rguide.xml enforced
# by Chris:
# Exactly one space either side of +, -, <-, <, >, ==, etc.
# Always a space after a comma, never before.
# Opening brackets should have a space before but not after, closing brackets
# should have a space after but not before.
# Line length <= 80 characters.

# TODO:
# 1. Code reads in lengthy patient IDs, that don't correspond to reads
# 2. Update plot of patients statistics & cross window mean of patient
# statistics
# 3. Code is now slow:
# Main problem is that repeatedly generating and searching all possible
# subclades is just plain stupid
# Solution: generate subclades only once, document their size, and then
# thin them as the loop goes
# 4. Plots should be by midpoint of the window, not start coordinate


command.line <- FALSE

if (command.line) {
	# Read in the arguments
	args <- commandArgs(TRUE)
	if (length(args) < 7) {
		cat(paste("At least 7 arguments must be specified:\n* a file containing",
						"the IDs to be coloured, one per line;\n* the string/character",
						"that follows the ID part of the tip name (i.e. those tip names to",
						"be coloured should equal one of the IDs followed by this",
						"string/character then anything at all);\n* the font size for the",
						"tip labels;\n* the line width for the tree;\n* a seed for",
						"randomising the order of colour allocation to the IDs file (use",
						"-1 to turn off shuffle;\n* the directory where output tree pdfs",
						"will be produced;\n* finally, any number of tree",
            "files.\nQuitting.\n"))
		quit("no", 1)
	}
	if(verbose) cat("Reading command line arguments\n")
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
  tree.files.basenames <- lapply(tree.files, basename)
	how.to.read.ids <- "List"
} else {
	# Some options to run within R
	rm(list=ls()) 

  root.folder <- "/Users/twoseventwo/Dropbox (Infectious Disease)"
  nameFolder <- function(name) paste(root.folder, name, sep = "")
	
	# Folder locations
  setwd(nameFolder("/BEEHIVE"))
	input.dir <- nameFolder("/BEEHIVE/phylotypes/run20160211/")
	tree.files.basenames <- list.files(input.dir, pattern="RAxML_bestTree.*\\.tree$")
	# list of all the trees to analyse
	# has to include \\. otherwise ignores the dot
  tree.files <- lapply(list.files(input.dir, pattern="RAxML_bestTree.*\\.tree$"),
                       function(x) file.path(input.dir, x))
	
	output.dir <- nameFolder("/BEEHIVE/phylotypes/run20160211/")
	
	how.to.read.ids <- "Trees"
	# Options are: "List", "Alignment" and "Trees"
	# "List" reads in all the patient IDs from a list in file 'id.file'
	# "Alignment" reads in the patient IDs from the master alignment
  # 'alignment.file'
	# !!! Option "Alignment" is somehow messed up
	# "Trees" reads IDs from the tip labels of the trees

	#id.file <- "ListOfIDs_run20160111.txt"
	alignment.file <- nameFolder("/BEEHIVE/phylotypes/run20160211/RefsAln.fasta")
	alignment.file.tag <- "-1_ref.fasta"
	# only reads IDs from file names that contain this string
	id.delimiter <- "_" # string that precedes this is the patient ID
	
	font.size <- 0.15 # for the tip labels in the tree
	line.width <- 1 # in the tree
	seed <- -1 # seed for randomising the colours associated with IDs
	# no randomisation if seed <- -1
}
verbose <- TRUE
# The tip name associated with the root of the phylogeny:
#root.name <- "Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
root.name <- "B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
#root.name <- "B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
print.trees <- T # whether or not to print separate PDF files for each tree
exit.after.colouring.trees <- FALSE

if(verbose) cat("Loading packages\n")
require(phytools)
require(dplyr)
require(ggplot2)
require(plotrix)
require(RColorBrewer)

if (!dir.exists(output.dir)) dir.create(output.dir)

# Utlities ##################
unfactorDataFrame <- function(x) {
	x <- data.frame(lapply(x, as.character), stringsAsFactors = F)
	x <- data.frame(lapply(x, function(y) type.convert(y, as.is=T)),
  stringsAsFactors = F)
}

#############################

# Read in the trees.

trees <- list()
num.trees <- length(tree.files)

if(verbose) cat("Loading trees\n")
for (tree.number in 1:num.trees) {
  if(verbose) cat("Opening file: ",tree.files.basenames[[tree.number]],"\n",sep = "")
  
  tree.file <- tree.files[[tree.number]]
  tree.file.basename <- tree.files.basenames[[tree.number]]
  tree <- read.tree(file=tree.file)
  trees[[tree.file]] <- tree

}

# Read in the IDs. Remove duplicates. Shuffle their order if desired.

ids <- vector()

if (how.to.read.ids == "Trees") {
  if(verbose) cat("Reading IDs from trees\n")
  for (tree.file in tree.files) {
    if(verbose) cat("Reading from tree: ",tree.file,"\n",sep = "")
    
    id.list <- sapply(trees[[tree.file]]$tip.label, function(x) strsplit(x,
    split = "_")[[1]][1])
#    id.list <- sapply(id.list, function(x) strsplit(x,
#    split = id.delimiter)[[1]][1])
    id.list <- unique(id.list)
    ids <- c(ids, id.list)
  }
  ids <- unique(ids)
  ids <- ids[ order(ids) ]
  ids <- ids[ grep("BEE", ids) ]

} else if (how.to.read.ids == "Alignment") {
  if(verbose) cat("Reading IDs from alignment\n")
  # !!!!!This option seems messed up, and likely won't work!!!!
	# reads in the IDs from the list of files in the folder - unusually named
  # files might get a bit garbled.
	ids <- list.files(id.folder, pattern = "*\\.fasta")
	ids <- sapply(ids, function (x) unlist(strsplit(x,
  split = id.folder.tag)[[1]][1]))
} else {
	# alternatively reads in the IDs from a specific file
  if(verbose) cat("Reading IDs from file\n")
	ids <- scan(id.file, what="", sep="\n", quiet=TRUE)
}

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
if(verbose) cat("Generating colours\n")
id.colours <- setNames(palette(rainbow(num.ids))[1:num.ids], ids)
id.colours <- setNames(palette(rainbow(num.ids))[1:num.ids], ids)

# Define the data frame to which summary statistics will be written
# First row is later deleted
pat.stats <- data.frame(patient.id = "test",
		window = 0,
		num.reads.total = 1,
		num.leaves = 1,
		num.clades = 1,
		overall.root.to.tip = 0,
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

if(verbose) cat("Processing trees\n")

for (tree.number in 1:num.trees) {
  
  if(verbose) cat("Tree: ",tree.number," of ",num.trees,"\n",sep = "")
	
  tree.file <- tree.files[[tree.number]]
  tree.file.basename <- tree.files.basenames[[tree.number]]
  tree <- trees[[tree.file]]
	
	tree <- root(phy = tree, outgroup = root.name)
	num.tips <- length(tree$tip.label)
	
	if (print.trees) {
		
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
		out.tree.basename <- paste(tree.file.basename, ".pdf", sep="")
		out.tree <- file.path(output.dir, out.tree.basename)
		pdf(out.tree, height = num.tips/10, width = 5)
		opar <- par()
		par(fg="transparent")
		plotTree(tree, color="black", fsize=font.size, lwd=line.width,
    ylim=c(-1, num.tips))
		lastPP <- get("last_plot.phylo", env=.PlotPhyloEnv)
		par(fg="black")
		text(lastPP$xx[1:num.tips], lastPP$yy[1:num.tips], tree$tip.label,
				cex=font.size, col=this.tree.colours, pos=4, offset=0.1, font=0.1)
		dev.off()
		par <- opar
	}

  if (exit.after.colouring.trees) {
    next
  }
	
	
	# Associate each tip to its patient.
	patient.tips <- list()
	for (id in ids) patient.tips[[id]] <- vector()
	for (tip.label in tree$tip.label) {
		id <- unlist(strsplit(tip.label, id.delimiter))[1]
		if (id %in% ids) {
			patient.tips[[id]] <- c(patient.tips[[id]], tip.label)
		}
	}
	
	# Is there a more robust way to pull out the window coordinate?
	window <- as.numeric(strsplit(tree.file, "_")[[1]][3])
	
	# For each patient:
	# Determine all the distinct monophyletic clades that the patient's tips make
	# N.b. a single isolated read is considered a type of monophyletic clade
	# Count the number of reads associated with each clade
	# ordered.clades is a list of these clades, ordered by decreasing number of
	# reads Compute the overall root-to-tip distance for all reads associated with
	# a patient and for each of the monophyletic subclades
	# Root to tip distance is considered zero for single isolated read
	# Record details for the five largest subclades in the summary statistic data
  # frame
	
	makeSingleTipTree <- function(tip.label, num.reads) {
		# creates a clade object based on a single tip label
		# num.reads is a named vector of all number of reads for all labels
		return(list(tree = tip.label, is.a.tip = T,
    num.reads = num.reads[tip.label]))
	}
	
	getLargestClade <- function(patient.subtree, tree, num.reads) {
		# from a 'patient.subtree' of 'tree', this function returns the largest
		# monophyletic subtree, as a clade object.
		# Clade objects can be a single tip.
		# This function only currently works if the patient.subtree is a tree (not
		# a clade object). This could be modified, which might simplify the main
		# code for finding all monophyletic subclades
		subtrees <- subtrees(patient.subtree )
		subtrees[[length(subtrees) + 1]] <- patient.subtree
		subtrees.monophyletic <- sapply(subtrees, function(x) is.monophyletic(tree,
    tips = x$tip.label))
		monophyletic.subtrees <- subtrees[ subtrees.monophyletic == T ]
		
		num.reads.largest.clade <- 0
		if (length(monophyletic.subtrees) > 0) {
			# if there is a monphyletic subtree
			for (subtree in monophyletic.subtrees) {
				num.reads.sub <- 0
				for (tip in subtree$tip.label) {
					num.reads.sub <- num.reads.sub + num.reads[tip]
				}
				if (num.reads.sub > num.reads.largest.clade) {
					num.reads.largest.clade <- num.reads.sub
					largest.mononophyletic.clade <- subtree
					largest.mononophyletic.clade.is.a.tip <- F
				}
			}
		}
		for (tip in patient.subtree$tip.label) {
			# this will either find the most common tip if there isn't a monophyletic
			# subtree, or find a tip which isn't in the monophyletic subtree which is
			# more common than that the sum of all the reads in that tree.
			# For that second option to work, this loop must come after the one above.
			num.reads.sub <- num.reads[tip]
			if (num.reads.sub > num.reads.largest.clade) {
				num.reads.largest.clade <- num.reads.sub
				largest.mononophyletic.clade <- tip
				largest.mononophyletic.clade.is.a.tip <- T
			}
		}
		return(list(tree = largest.mononophyletic.clade,
    is.a.tip = largest.mononophyletic.clade.is.a.tip,
    num.reads = unname(num.reads.largest.clade)))
	}
	
	calcMeanRootToTip <- function(clade, num.reads) {
		# Mean root-to-tip distance, weighted by the number of reads associated with
		# each tip. Returns 0 if the clade is a single tip.
		# Input is a clade object (a named list), with three items:
		# $is.a.tip is a logical, indicating whether the tree is a single tip of a
		# 'proper' tree
		# $tree is either the tree, or the label of the tip
		# $num.reads is the number of reads associated with the whole clade
		# The second input is num.reads, a vector with named entries, which
		# returns the number of reads associated with each tip label in the tree
		root.to.tip <- 0
		if (clade$is.a.tip == F) {
			for (i in 1:length(clade$tree$tip.label)) {
				root.to.tip <- root.to.tip +
						nodeheight(clade$tree, i) * num.reads[ clade$tree$tip.label[i] ]
			}
			root.to.tip <- root.to.tip/clade$num.reads
		}
		return(unname(root.to.tip))
	}
	
	
	for (i in 1:num.ids) {
		# The main looop for finding all sub-clades of the patient-specific subtree
		id <- ids[i]
		num.leaves <- length(patient.tips[[id]])
		num.reads <- list()
		for (tip in patient.tips[[id]]) {
      num.reads[[tip]] <- as.numeric(unlist(strsplit(tip, "count_"))[2])
    }
		num.reads <- unlist(num.reads)
		num.reads.total <- sum(num.reads)
		ordered.clades <- list()
		current.clade <- 1
		
		if (num.leaves > 0) {
			# otherwise, there were no reads
			if (num.leaves == 1) {
				# only one tip for that patient
				ordered.clades[[current.clade]] <- makeSingleTipTree(tip, num.reads)
				overall.root.to.tip <- 0
			} else if (num.leaves > 1) {
				# define the subtree of all tips associated with the patient
				patient.subtree <- drop.tip(phy = tree,
						tip = tree$tip.label[!(tree$tip.label %in% patient.tips[[id]])])
				patient.clade <- list(tree = patient.subtree, is.a.tip = F,
        num.reads = num.reads.total)
				overall.root.to.tip <- calcMeanRootToTip(patient.clade, num.reads)
				ordered.clades[[current.clade]] <- getLargestClade(patient.subtree,
        tree, num.reads)
				# first find the largest sub clade
				
				done <- F
				while (done == F) {
					# Then find all the other subclades, in order
					# It would probably be much more efficient to find all of the
					# subclades in one go, and then order them.
					if (ordered.clades[[current.clade]]$is.a.tip) {
						test.length <- 1
						tips <- ordered.clades[[current.clade]]$tree
					} else {
						test.length <-
            length(ordered.clades[[current.clade]]$tree$tip.label)
						tips <- ordered.clades[[current.clade]]$tree$tip.label
					}
					
					if (test.length < (length(patient.subtree$tip.label) - 1)) {
						patient.subtree <- drop.tip(patient.subtree,
            tip = patient.subtree$tip.label[(
            patient.subtree$tip.label %in% tips)])
						current.clade <- current.clade + 1
						next.clade <- getLargestClade(patient.subtree, tree, num.reads)
						ordered.clades <- c(ordered.clades, list(next.clade))
					} else if (test.length == (length(patient.subtree$tip.label) - 1)) {
						tip.label <-patient.subtree$tip.label[which(!(
            patient.subtree$tip.label %in% tips))]
						current.clade <- current.clade + 1
						ordered.clades <- c(ordered.clades,
								list(makeSingleTipTree(tip.label, num.reads)))
						done <- T
					} else {
						done <- T
					}
				}
			}
			num.clades <- length(ordered.clades)
			# The patient subtree has now been split into num.clades clades
			
			prop.reads.per.clade <- root.to.tip.per.clade <- vector(length = 5)
			# We only compute data on the five largest clades
			root.to.tip.per.clade
			for (i in 1:5) {
				if (i > num.clades) {
					prop.reads.per.clade[i] <- 0
					root.to.tip.per.clade[i] <- NA
					# zeros are plotted, NAs are ommitted
				} else {
					prop.reads.per.clade[i] <-
							ordered.clades[[i]]$num.reads/num.reads.total
					root.to.tip.per.clade[i] <-
							calcMeanRootToTip(ordered.clades[[i]], num.reads)
				}
			}
			print(paste(window, id, num.clades,  num.reads.total),
					round(prop.reads.per.clade[1], 3))
		} else { # if there are no reads for this patient in this window
			num.clades <- NA
			ordered.clades <- NA
			overall.root.to.tip <- NA
			prop.reads.per.clade <- NA
			root.to.tip.per.clade <- NA
			print(paste(i, id, "no reads"))
		}
		pat.stats <- rbind(pat.stats, c(id,
						window,
						num.reads.total,
						num.leaves,
						num.clades,
						overall.root.to.tip,
						prop.reads.per.clade[1],
						prop.reads.per.clade[2],
						prop.reads.per.clade[3],
						prop.reads.per.clade[4],
						prop.reads.per.clade[5],
						root.to.tip.per.clade[1],
						root.to.tip.per.clade[2],
						root.to.tip.per.clade[3],
						root.to.tip.per.clade[4],
						root.to.tip.per.clade[5]))
	}
	
}

if (exit.after.colouring.trees) {
  print("All trees coloured; exiting successfully.")
  quit("no", 0)
}

pat.stats <- pat.stats[!pat.stats$patient.id=="test", ]
# Delete the first dummy line
pat.stats <- pat.stats[order(pat.stats$patient.id), ]
# TODO: replace pat.stats$num.leaves = 0 by NA

pdf(file.path(output.dir, "Patient.Statistics.pdf"), paper = "a4")
for (id in sort(ids)) {
	par(mfrow=c(3, 2))
	# putting this here also sets new page, in case there aren't six plots
	pat <- pat.stats[pat.stats$patient.id==id, ]
	pat <- unfactorDataFrame(pat)
	pat <- pat[order(pat$window), ]
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


# create the ggplot2 data, using Year and thousands (of people) and filled by
# age group.
#    clade.1 <- pat[ ,c("patient.id", "window", "prop.reads.clade.1")]
#    clade.2 <- pat[ ,c("patient.id", "window", "prop.reads.clade.2")]
#    clade.3 <- pat[ ,c("patient.id", "window", "prop.reads.clade.3")]
#    clade.4 <- pat[ ,c("patient.id", "window", "prop.reads.clade.4")]
#    clade.5 <- pat[ ,c("patient.id", "window", "prop.reads.clade.5")]
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
#   boxplot(mean.size~window, data = pat.stats[pat.stats$monophyletic==1, ],
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
