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

# We flip the value of the bool command.line to determine whether arguments are
# read in from the command line (TRUE) or set internally in the code (FALSE).
command.line <- TRUE
if (command.line) {
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
	# only reads IDs from file names that contain this string
	id.delimiter <- "-1" # string that precedes this is the patient ID
	
	font.size <- 0.15 # for the tip labels in the tree
	line.width <- 1 # in the tree
	seed <- -1 # seed for randomising the colours associated with IDs
	# no randomisation if seed <- -1
}

# The tip name associated with the root of the phylogeny:
#root.name <- "Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
#root.name <- "B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
root.name <- "C.BW.00.00BW07621.AF443088"
print.trees <- T # whether or not to print separate PDF files for each tree
exit.after.colouring.trees <- F
print.per.patient.pdfs <- F
num.clades.for.output <- 5
use.christophes.code <- T

require(phytools)
require(dplyr)
require(ggplot2)
require(plotrix)
require(RColorBrewer)
require(phangorn)

if (!dir.exists(output.dir)) dir.create(output.dir)

# Utlities ##################
unfactorDataFrame <- function(x) {
	x <- data.frame(lapply(x, as.character), stringsAsFactors = F)
	x <- data.frame(lapply(x, function(y) type.convert(y, as.is=T)),
  stringsAsFactors = F)
}

#############################

# Check that tree files are named as expected, so that we know how to extract
# the window coordinates from them.
for (tree.file.basename in tree.files.basenames) {
 if (! grepl("^RAxML_bestTree\\.InWindow_[0-9]*_to_[0-9]*\\.tree",
  tree.file.basename)) {
    cat(paste("Unexpected name format for tree file ", tree.file.basename,
    '; expected\nRAxML_bestTree.InWindow_[0-9]*_to_[0-9]*.tree\nQuitting.\n',
    sep=''))
    quit("no", 1)
  }
}

# Read in the trees.
trees <- list()
num.trees <- length(tree.files)
for (tree.number in 1:num.trees) {
  tree.file <- tree.files[[tree.number]]
  tree <- read.tree(file=tree.file)
  trees[[tree.file]] <- tree
}

  

# Read in the IDs.
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
	ids <- scan(id.file, what="", sep="\n", quiet=TRUE)
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

# Define a colour for each unique ID.
# Needs to be written twice to work properly for some reason!
id.colours <- setNames(palette(rainbow(num.ids))[1:num.ids], ids)
id.colours <- setNames(palette(rainbow(num.ids))[1:num.ids], ids)

# Define the data frame to which summary statistics will be written.
# (The first row is later deleted.)
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


makeSingleTipTree <- function(tip.label, read.counts.per.tip) {
	# creates a clade object based on a single tip label
	# read.counts.per.tip is a named vector of all number of reads for all labels
	return(list(tree = tip.label, is.a.tip = T,
   num.reads = read.counts.per.tip[tip.label]))
}
	
getLargestClade <- function(patient.subtree, tree, read.counts.per.tip) {
	# from a 'patient.subtree' of 'tree', this function returns the largest
	# monophyletic subtree, as a clade object.
	# Clade objects can be a single tip.
  # This function only currently works if the patient.subtree is a tree (not
	# a clade object). This could be modified, which might simplify the main
	# code for finding all monophyletic subclades
	subtrees <- subtrees(patient.subtree)
	subtrees[[length(subtrees) + 1]] <- patient.subtree
	subtrees.monophyletic <- sapply(subtrees, function(x) is.monophyletic(tree,
  tips = x$tip.label))
	monophyletic.subtrees <- subtrees[subtrees.monophyletic == T]
		
	num.reads.largest.clade <- 0
	if (length(monophyletic.subtrees) > 0) {
		# if there is a monphyletic subtree
		for (subtree in monophyletic.subtrees) {
			num.reads.sub <- 0
			for (tip in subtree$tip.label) {
				num.reads.sub <- num.reads.sub + read.counts.per.tip[tip]
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
		num.reads.sub <- read.counts.per.tip[tip]
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
	
calcMeanRootToTip <- function(clade, read.counts.per.tip) {
	# Mean root-to-tip distance, weighted by the number of reads associated with
	# each tip. Returns 0 if the clade is a single tip.
	# Input is a clade object (a named list), with three items:
	# $is.a.tip is a logical, indicating whether the tree is a single tip of a
	# 'proper' tree
	# $tree is either the tree, or the label of the tip
	# $num.reads is the number of reads associated with the whole clade
	# The second input is read.counts.per.tip, a vector with named entries, which
	# returns the number of reads associated with each tip label in the tree
	root.to.tip <- 0
	if (clade$is.a.tip == F) {
		for (i in 1:length(clade$tree$tip.label)) {
			root.to.tip <- root.to.tip + nodeheight(clade$tree, i) * 
			read.counts.per.tip[clade$tree$tip.label[i]]
		}
		root.to.tip <- root.to.tip/clade$num.reads
	}
	return(unname(root.to.tip))
}

calcMeanRootToTipChris <- function(tree, read.counts.per.tip) {
	# Mean root-to-tip distance, weighted by the number of reads associated with
	# each tip. Returns 0 if the clade is a single tip.
	# The second input is read.counts.per.tip, a vector with named entries, which
	# returns the number of reads associated with each tip label in the tree.
	root.to.tip <- 0
  num.tips <- length(tree$tip.label)
	if (num.tips > 1) {
    num.reads.total <- 0
		for (tip.number in 1:length(tree$tip.label)) {
      tip.count <- read.counts.per.tip[tree$tip.label[tip.number]]
			root.to.tip <- root.to.tip + nodeheight(tree, tip.number) * tip.count
      num.reads.total <- num.reads.total + tip.count
		}
    stopifnot(num.reads.total > 0)
		root.to.tip <- root.to.tip / num.reads.total
	}
	return(root.to.tip)
}

getIdFromTip <- function(tip.label) {
  # Gets the part of a tip name that will be the patient id for patient tips.
  return(unlist(strsplit(tip.label, id.delimiter))[1])
}

resolveTreeIntoPatientClades <-
function(tree, first.call=TRUE, num.tips=NULL, node=NULL) {
  # Resolves a tree into a list (indexed by patient ID) of lists; each entry in
  # in one of the latter lists is a monophyletic set of tips for that patient, 
  # collected into a vector.

  # This function calls itself iteratively. If this is the first call, set the 
  # node to be the tree root and find the number of tips (used repeatedly).
  if (first.call) {
    node <- findMRCA(tree, tree$tip.label, type="node")
    num.tips <- length(tree$tip.label)
  }

  # This function shouldn't be called on a 'node' that's actually a tip.
  stopifnot(node > num.tips)

  # Create a list, indexed by patient ids, of empty lists.
	clades.at.this.level <- list()
	for (id in ids) clades.at.this.level[[id]] <- list()

  # If all tips descending from this node are from one patient, add their names
  # to this patient's list.
  descendant.tips <- Descendants(tree, node=node, type="tips")
  descendant.tips <- rapply(descendant.tips,c)
  first.id <- NULL
  all.one.patient <- TRUE
  for (tip.number in descendant.tips) {
    tip.name <- tree$tip.label[tip.number]
    id <- getIdFromTip(tip.name)
    if (! id %in% ids) {
      all.one.patient <- FALSE
      break
    } else if (is.null(first.id)) {
      first.id <- id
    } else if (id != first.id) {
      all.one.patient <- FALSE
      break
    }
  }
  if (all.one.patient) {
    clades.at.this.level[[first.id]][[1]] <- tree$tip.label[descendant.tips]
  } else {
    # Iterate through the immediate children of this node.
    children <- Descendants(tree, node=node, type="children")
    for (child in children) {
      if (child <= num.tips) {
        # The child is a tip.
        # If it's a patient tip, add it to that patient's vector.
        tip.name <- tree$tip.label[child]
        id <- getIdFromTip(tip.name)
        if (id %in% ids) {
          clades.at.this.level[[id]][length(clades.at.this.level[[id]])+1] <-
          c(tip.name)
        }
      } else {
        # The child is a node. Merge clades.at.this.level with the result of
        # applying this function to the child.
        clades.at.this.level <- Map(c, clades.at.this.level,
        resolveTreeIntoPatientClades(tree, FALSE, num.tips, child))
      }
    }
  }
  return(clades.at.this.level)
}

# Iterate through the trees
for (tree.number in 1:num.trees) {
	
  tree.file <- tree.files[[tree.number]]
  tree.file.basename <- tree.files.basenames[[tree.number]]
  tree <- trees[[tree.file]]
	num.tips <- length(tree$tip.label)

  # Root the tree
	tree <- root(phy = tree, outgroup = root.name)

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
    #nodelabels(bg="white", cex=font.size)
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
	
	# Pull out the window coordinates
  filename.trimmed <- gsub("\\.tree", "", tree.file.basename)
  filename.components <- strsplit(filename.trimmed, "_")[[1]]
	left.window.edge  <- as.numeric(filename.components[3])
	right.window.edge <- as.numeric(filename.components[5])
  print(paste('Now analysing window', left.window.edge, "to",
  right.window.edge))
  window <- (left.window.edge + right.window.edge) / 2

  if (! use.christophes.code) {
  print('calling function now...')
  all.clades.by.patient <- resolveTreeIntoPatientClades(tree)
  print('...done')

  # Iterate through all patients
  for (id in ids) {

    tips <- patient.tips[[id]]
    num.tips <- length(tips)
    if (num.tips == 0) {
      num.reads.total <- 0
      num.clades <- NA
			overall.root.to.tip <- NA
			prop.reads.per.clade  <- rep(NA, num.clades.for.output)
			root.to.tip.per.clade <- rep(NA, num.clades.for.output)
    } else {
    
      # Find this patient's clades
      clades = all.clades.by.patient[[id]]
      num.clades <- length(clades)
      read.counts.per.clade <- vector(mode="integer", length=num.clades)
      read.counts.per.tip <- list()
      current.clade <- 0
      for (clade in clades) {
        current.clade <- current.clade + 1
        num.reads.this.clade <- 0
        for (tip in clade) {
          num.reads.this.tip <- as.numeric(unlist(strsplit(tip, "count_"))[2])
          num.reads.this.clade <- num.reads.this.clade + num.reads.this.tip
          read.counts.per.tip[[tip]] <- num.reads.this.tip
        }
        read.counts.per.clade[current.clade] <- num.reads.this.clade
      }
      read.counts.per.tip <- unlist(read.counts.per.tip)
      num.reads.total <- sum(read.counts.per.clade)

      # Define the subtree of all tips associated with the patient. Find the root-
      # to-tip distance.
      if (num.tips == 1) {
        overall.root.to.tip <- 0
			  prop.reads.per.clade  <- c(1, rep(NA, num.clades.for.output -1))
			  root.to.tip.per.clade <- c(0, rep(NA, num.clades.for.output -1))
      } else {
		    patient.subtree <- drop.tip(phy=tree,
        tip=tree$tip.label[!(tree$tip.label %in% patient.tips[[id]])])
		    overall.root.to.tip <- calcMeanRootToTipChris(patient.subtree, 
        read.counts.per.tip)

        # Order clades by size (how many reads they contain). Find the root-to-tip
        # distance and propotion of reads in the largest few.
        clade.ordering <- order(read.counts.per.clade, decreasing=TRUE)
        ordered.clade.counts <- read.counts.per.clade[clade.ordering]
        ordered.clades <- clades[clade.ordering]
        prop.reads.per.clade <- root.to.tip.per.clade <- vector(length = 5)
        for (i in 1:num.clades.for.output) {
			    if (i > num.clades) {
				    prop.reads  <- 0
				    root.to.tip <- NA
			    } else {
            prop.reads <- ordered.clade.counts[i] / num.reads.total
            clade <- ordered.clades[i]
            if (length(clade) == 1) {
              root.to.tip <- 0
            } else {
              print('finding clade subtree')
              print('clade tips:')
              print(clade)
              print(tree$tip.label[!(tree$tip.label %in% clade)])
		          clade.subtree <- drop.tip(phy=tree,
              tip=tree$tip.label[!(tree$tip.label %in% clade)])
              root.to.tip <- calcMeanRootToTipChris(clade, read.counts.per.tip)
            }
          }
    		  prop.reads.per.clade[i]  <- prop.reads
			    root.to.tip.per.clade[i] <- root.to.tip
        }

      }
    }
		pat.stats <- rbind(pat.stats, c(id, window, num.reads.total, num.tips,
    num.clades, overall.root.to.tip, prop.reads.per.clade,
    root.to.tip.per.clade))
  }
      
	} else {
	# For each patient:
	# Determine all the distinct monophyletic clades that the patient's tips make.
	# (NB a single isolated read is considered a type of monophyletic clade.)
	# Count the number of reads associated with each clade. ordered.clades is a
	# list of these clades, ordered by decreasing number of reads.
  # Compute the overall root-to-tip distance for all reads associated with
	# a patient and for each of the monophyletic subclades.
	# (Root to tip distance is considered zero for single isolated read.)
	# Record details for the five largest subclades in the summary statistic data
  # frame.
  # Iterate through all ids (i.e. all patients):
	for (i in 1:num.ids) {

    # Find how many reads are associated with each tip, and the total number of
    # reads
		id <- ids[i]
		num.leaves <- length(patient.tips[[id]])
		read.counts.per.tip <- list()
		for (tip in patient.tips[[id]]) {
      read.counts.per.tip[[tip]] <- as.numeric(unlist(strsplit(tip,
      "count_"))[2])
    }
		read.counts.per.tip <- unlist(read.counts.per.tip)
		num.reads.total <- sum(read.counts.per.tip)

		ordered.clades <- list()
		current.clade <- 1
		if (num.leaves > 0) {
			if (num.leaves == 1) {
				ordered.clades[[current.clade]] <- makeSingleTipTree(tip,
        read.counts.per.tip)
				overall.root.to.tip <- 0
			} else {
				# define the subtree of all tips associated with the patient
				patient.subtree <- drop.tip(phy = tree,
						tip = tree$tip.label[!(tree$tip.label %in% patient.tips[[id]])])
				patient.clade <- list(tree = patient.subtree, is.a.tip = F,
        num.reads = num.reads.total)
				overall.root.to.tip <- calcMeanRootToTip(patient.clade,
        read.counts.per.tip)
				ordered.clades[[current.clade]] <- getLargestClade(patient.subtree,
        tree, read.counts.per.tip)
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
						next.clade <- getLargestClade(patient.subtree, tree,
            read.counts.per.tip)
						ordered.clades <- c(ordered.clades, list(next.clade))
					} else if (test.length == (length(patient.subtree$tip.label) - 1)) {
						tip.label <-patient.subtree$tip.label[which(!(
            patient.subtree$tip.label %in% tips))]
						current.clade <- current.clade + 1
						ordered.clades <- c(ordered.clades,
								list(makeSingleTipTree(tip.label, read.counts.per.tip)))
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
							calcMeanRootToTip(ordered.clades[[i]], read.counts.per.tip)
				}
			}
			#print(paste(window, id, num.clades,  num.reads.total),
			#		round(prop.reads.per.clade[1], 3))
		} else { # if there are no reads for this patient in this window
			num.clades <- NA
			ordered.clades <- NA
			overall.root.to.tip <- NA
			prop.reads.per.clade <- rep(NA, 5)
			root.to.tip.per.clade <- rep(NA, 5)
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
	
}

if (exit.after.colouring.trees) {
  print("All trees coloured; exiting successfully.")
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
