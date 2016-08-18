checkArgs <- function(args) {
  # Quits with an error message if the wrong number of args is given.
	if (length(args) < 6) {
		cat(paste("At least 6 arguments must be specified:\n* a file containing",
						"the IDs to be coloured, one per line;\n* the font size for the",
						"tip labels;\n* the line width for the tree;\n* a seed for",
						"randomising the order of colour allocation to the IDs file (use",
						"-1 to turn off shuffle;\n* the directory where output tree pdfs",
						"will be produced;\n* finally, any number of tree",
            "files.\nQuitting.\n"))
		quit("no", 1)
	}
}

checkArgs2 <- function(args) {
  # Quits with an error message if the wrong number of args is given.
	if (length(args) < 2) {
		cat(paste("At least 2 arguments must be specified:\n* a file containing",
						"the IDs to be coloured, one per line; and any number of tree",
            "files.\nQuitting.\n"))
		quit("no", 1)
	}
}

checkTreeFiles <- function(tree.files) {
  # Check files exist
	for (tree.file in tree.files) {
		if (! file.exists(tree.file)) {
			print(paste(tree.file, "does not exist. Quitting."))
			quit("no", 1)
		}
	}
}

checkTreeFileNames <- function(tree.files.basenames, tree.file.regex) {
  # Check that tree files are named as expected, so that we know how to extract
  # the window coordinates from them.
  for (tree.file.basename in tree.files.basenames) {
   if (! grepl(tree.file.regex, tree.file.basename)) {
      cat(paste("Unexpected name format for tree file ", tree.file.basename,
      '. Quitting.\n',
      sep=''))
      quit("no", 1)
    }
  }
}

colourTree <- function(tree, tip.regex, id.colours, tree.file.basename,
output.dir, font.size, line.width, ids, colour.edges, clade.mrcas.by.patient,
ladderise.trees) {
  # Prints a tree coloured by patient ID

  num.tips <- length(tree$tip.label)

  # Assign each tip its colour
	this.tree.colours <- vector(length=num.tips)
	for (i in seq_len(num.tips)) {
		tip.label <- tree$tip.label[i]
		id <- sub(tip.regex, "\\1", tip.label)
		if (id %in% ids) {
			this.tree.colours[i] <- id.colours[which(ids==id)]
		} else {
			this.tree.colours[i] <- "black"
		}
	}
		
	# Rotate internal nodes to have right daughter clades more species rich than
	# left daughter clades.
  if (ladderise.trees) {
	  tree <- ladderize(tree)
  }

  if (colour.edges) {
    # Generate a string guaranteed not to clash with any id: the character "a"
    # repeated more times than the length of any id. We use this to label the
    # ancestral state for edge colouring.
    longest.id.length <- max(nchar(ids))
    anc.state <- paste(rep("a", longest.id.length+1), collapse = "")
    edge.colours <- id.colours
    edge.colours[[anc.state]] <- "black"

    # Colour the edges associated with each patient clade
    for (id in ids) {
      for (clade.mrca in clade.mrcas.by.patient[[id]]) {
        tree <- paintSubTree(tree, node=clade.mrca, stem=1, state=id,
        anc.state=anc.state)
      }
    }
  }
		
	# Plot the tree.
	# First make the automatic tip labels transparent, then replot with colour.
	out.tree.basename <- paste(tree.file.basename, ".pdf", sep="")
	out.tree <- file.path(output.dir, out.tree.basename)
	pdf(out.tree, height = num.tips/10, width = 5)
	opar <- par()
	par(fg="transparent")
  if (colour.edges) {
	  plotSimmap(tree, edge.colours, fsize=font.size, lwd=line.width,
    ylim=c(-1, num.tips))
  } else {
    plotTree(tree, color="black", fsize=font.size, lwd=line.width,
    ylim=c(-1, num.tips))
    #,node.numbers=T)
  }
  lastPP <- get("last_plot.phylo", env=.PlotPhyloEnv)
	par(fg="black")
  #nodelabels(bg="white", cex=font.size)
	text(lastPP$xx[seq_len(num.tips)], lastPP$yy[seq_len(num.tips)], tree$tip.label,
			cex=font.size, col=this.tree.colours, pos=4, offset=0.1, font=0.1)

  #node.labels <- sapply(as.numeric(as.character(tree$node.label)), 
  #               function(bs){
  #                 if(!is.na(bs) & bs>30){
  #                   return(bs)
  #                 } else {
  #                   return("")
  #                 }
  #               })
  #nodelabels(node.labels, (length(tree$tip.label) + 1):(length(tree$tip.label) +
  #tree$Nnode), adj = c(0.5,0.05), frame = "none", pch = NULL, thermo = NULL,
  #pie = NULL, piecol = NULL, col = "black", cex = 0.6)
  #add.scale.bar(0,-1,length=0.05)
	dev.off()
	par <- opar
}


resolveTreeIntoPatientClades <- function(tree, ids, tip.regex, blacklisted.tips = vector()) {

  num.tips <- length(tree$tip.label)
  num.patients <- length(ids)

  # Reorder the nodes in such a manner that we always encounter a descendant 
  # node before its ancestor (not a unique process but that doesn't matter).
  tree <- reorder(tree, "postorder")

  # Make a vector of bools describing which tips are patients. Remove blacklisted tips from this vector.
  is.patient <- grep(tip.regex, tree$tip.label)
  is.patient <- is.patient[which(!(is.patient %in% blacklist))]
  
  tip.ids <- sub(tip.regex, "\\1", tree$tip.label[is.patient])
  is.patient <- is.patient[tip.ids %in% ids]
  tip.ids <- tip.ids[tip.ids %in% ids]

  # Make a matrix that has one row for every tip & one row for every node, and 
  # one column for each patient. All values are false, for columns where that 
  # patient and tip/node are associated (i.e. for tips, the tip IS that patient; 
  # for nodes, that patient has a tip descendending from this node). The tip 
  # rows we initialise immediately; the node rows we will update as we parse 
  # through the nodes.
  num.tips.and.nodes <- num.tips + tree$Nnode
  m <- matrix(FALSE, num.tips.and.nodes, num.patients)
  m[cbind(is.patient, match(tip.ids, ids))] <- TRUE

  # Our main object of interest: each patient will have a list, with one item in 
  # that list being a list of monophyletic tips for that patient.
  groups <- vector("list", num.patients)
  clade.mrcas.by.patient <- vector("list", num.patients)

  # Find all nodes in the tree, and the children of each node
  nodes <- unique(tree$edge[,1])
  children <- split(tree$edge[,2], factor(tree$edge[,1], nodes))

  # Make a vector of lists, one for each tip & one for each node. For tips, the
  # value in the vector is just a one-component list containing the ID of the 
  # tip itself; for nodes, the will contain the IDs of all descendants from that
  # node.
  descendants <- vector("list", num.tips.and.nodes)
  descendants[seq_len(num.tips)] <- seq_len(num.tips)

  # Make a vector of bools "keep.going", one bool for each tip & one for each 
  # node. The tip values are true for patient tips and false for external ref 
  # tips. The node values are initialised to TRUE; we will update them as we 
  # parse through the nodes, and they only remain TRUE if all their descendants 
  # are from one patient. To decide this, we look at the children of a node 
  # (which can be nodes or tips). If all children are TRUE and they all 
  # correspond to the same patient (which is implied by the rows in m associated 
  # with these children all being the same and containing exactly one TRUE 
  # value), we 'keep going': this node gets to keep its TRUE status, we update 
  # the row of m for this node to match its children, and we set the descendants 
  # for this node to be the descendants of its children. If we don't keep going,
  # then for any child which is 'keep going' we append its descendants to that
  # patient's list in groups.
  keep.going <- rep(TRUE, num.tips.and.nodes)
  keep.going[!(seq_len(num.tips) %in% is.patient)] <- FALSE

  # Iterate through the nodes (we always meet ancestors nodes after their
  # descendants):
  for (i in seq_along(nodes)) {

    node <- nodes[[i]]
    these.children <- children[[i]]
    keep.going[node] <- (all(keep.going[these.children]) && 
                       sum(colSums(m[these.children, ]) > 0) == 1)

    if (keep.going[node]) {
      m[node, ] <- m[these.children[1], ]
      descendants[[node]] <- unlist(descendants[these.children])
    } else {
      for (j in these.children[keep.going[these.children]]) {
        k <- which(m[j,])
        groups[[k]] <- c(groups[[k]], descendants[j])
        clade.mrcas.by.patient[[k]] <- c(clade.mrcas.by.patient[[k]], j)
      }
    }
  }

  # Label tips and ids by strings instead of ints
  names(groups) <- ids
  names(clade.mrcas.by.patient) <- ids
  clades.by.patient <- lapply(groups, lapply, function(el) tree$tip.label[el])

  #read.counts.per.tip <- rep(NA_integer_, length(num.tips))
  #read.counts.per.tip[is.patient] <- as.integer(sub(re, "\\3",
  #tree$tip.label[is.patient]))
  #group.counts <- lapply(groups, vapply,
  #function(el) sum(read.counts.per.tip[el]), integer(1))

  return(list(clades.by.patient=clades.by.patient,
  clade.mrcas.by.patient=clade.mrcas.by.patient))
}


summariseClades <- function(tree, clades, num.clades.for.output, tip.regex,
whose.code) {
  # TODO: auto indent this function!

    num.clades <- length(clades)
    num.tips <- length(unlist(clades))

    # No tips: data is missing.
    if (num.tips == 0) {
      num.reads.total <- 0
      num.clades <- NA
			overall.root.to.tip <- NA
			prop.reads.per.clade  <- rep(NA, num.clades.for.output)
			root.to.tip.per.clade <- rep(NA, num.clades.for.output)

    # Only one tip; proportions and root-to-tip distances are trivial.
    } else if (num.tips == 1) {
      num.clades <- 1
  	  prop.reads.per.clade  <- c(1, rep(0, num.clades.for.output -1))
		  root.to.tip.per.clade <- c(0, rep(NA, num.clades.for.output -1))
      overall.root.to.tip <- 0
      tip.name <- unlist(clades)
      num.reads.total <- as.integer(sub(tip.regex, "\\2", tip.name))

    # At least two tips:
    } else {
      read.counts.per.clade <- vector(mode="integer", length=num.clades)
      read.counts.per.tip <- list()
      current.clade <- 0
      for (clade in clades) {
        current.clade <- current.clade + 1
        num.reads.this.clade <- 0
        for (tip in clade) {
          num.reads.this.tip <- as.integer(sub(tip.regex, "\\2", tip))
          num.reads.this.clade <- num.reads.this.clade + num.reads.this.tip
          read.counts.per.tip[[tip]] <- num.reads.this.tip
        }
        read.counts.per.clade[current.clade] <- num.reads.this.clade
      }
      read.counts.per.tip <- unlist(read.counts.per.tip)
      num.reads.total <- sum(read.counts.per.clade)

      # Define the subtree of all tips associated with the patient, and find
      # the root-to-tip distance.
		  patient.subtree <- drop.tip(phy=tree,
      tip=tree$tip.label[!(tree$tip.label %in% patient.tips[[id]])])
		  overall.root.to.tip <- calcMeanRootToTip(patient.subtree, 
      read.counts.per.tip)

      # If there is only one clade, proportions and root-to-tip distances are
      # trivial:
      if (num.clades == 1) {
			  prop.reads.per.clade  <- c(1, rep(0, num.clades.for.output -1))
			  root.to.tip.per.clade <- c(overall.root.to.tip,
        rep(NA, num.clades.for.output -1))
      } else {

        # There are multiple clades. Order them by size (how many reads they
        # contain). Find the root-to-tip distance and propotion of reads in
        # the largest few.
        clade.ordering <- order(read.counts.per.clade, decreasing=TRUE)
        ordered.clade.counts <- read.counts.per.clade[clade.ordering]
        ordered.clades <- clades[clade.ordering]
        prop.reads.per.clade <- root.to.tip.per.clade <- vector(length = 5)
        for (i in seq_len(num.clades.for.output)) {
			    if (i > num.clades) {
				    prop.reads  <- 0
				    root.to.tip <- NA
			    } else {
            prop.reads <- ordered.clade.counts[i] / num.reads.total
            clade <- unlist(ordered.clades[i])
            if (length(clade) == 1) {
              root.to.tip <- 0
            } else {
		          clade.subtree <- drop.tip(phy=tree,
              tip=tree$tip.label[!(tree$tip.label %in% clade)])
              root.to.tip <- calcMeanRootToTip(clade.subtree,
              read.counts.per.tip)
            }
          }
      		prop.reads.per.clade[i]  <- prop.reads
			    root.to.tip.per.clade[i] <- root.to.tip
        }

      }
    }
    return(list(num.reads.total=num.reads.total, num.tips=num.tips,
    num.clades=num.clades, overall.root.to.tip=overall.root.to.tip,
    prop.reads.per.clade=prop.reads.per.clade,
    root.to.tip.per.clade=root.to.tip.per.clade))
}

printSortedClades <- function(clades, tip.regex, num.clades.for.output) {
  # Order each clade by the read number of the tips, and the set of clades by
  # their smallest read number. This is unique (unlike read count sorting) so
  # helps test clade finding algorithms.
  # TODO: we're only printing the first (num.clades.for.output) clades, since
  # Christophe's algorithm will only calculate that many clades, but we print
  # the first (num.clades.for.output) clades after sorting in the manner done
  # here (essentially alphabetically), whereas Christophe's algorithm will find 
  # the first (num.clades.for.output) clades after sorting by read count then
  # stop, even if smaller clades would be amongst the first alphabetically.
  min.read.number.from.each.clade <- vector(length=length(clades))
  ordered.clades <- list("list")
  current.clade <- 0
  for (clade in clades) {
    current.clade <- current.clade + 1
    read.numbers <- as.integer(sub(tip.regex, "\\2", clade))
    min.read.number.from.each.clade[current.clade] <- min(read.numbers)
    ordered.clades[[current.clade]] <- clade[order(read.numbers)]
  }
  ordered.clades <- ordered.clades[order(min.read.number.from.each.clade)]
  ordered.clades <- ordered.clades[seq_len(min(length(ordered.clades),
  num.clades.for.output))]
  print(ordered.clades)
  cat('\n\n')
}

calcMeanRootToTip <- function(tree, read.counts.per.tip) {
  if(is.null(read.counts.per.tip)){
    read.counts.per.tip <- rep(1, length(tree$tip.label))
  }
	# Mean root-to-tip distance, optionally weighted by the number of reads associated with
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



################################################################################
# All functions below deprecated, for backward compatibility...
################################################################################

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
	
calcMeanRootToTipChristophe <- function(clade, read.counts.per.tip) {
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

getIdFromTip <- function(tip.label, id.delimiter) {
  # Gets the part of a tip name that will be the patient id for patient tips.
  return(unlist(strsplit(tip.label, id.delimiter))[1])
}
# calculate this once for all tips...

resolveTreeIntoPatientCladesChris <-
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
    id <- getIdFromTip(tip.name, id.delimter)
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
        id <- getIdFromTip(tip.name, id.delimter)
        if (id %in% ids) {
          clades.at.this.level[[id]][length(clades.at.this.level[[id]])+1] <-
          c(tip.name)
        }
      } else {
        # The child is a node. Merge clades.at.this.level with the result of
        # applying this function to the child.
        clades.at.this.level <- Map(c, clades.at.this.level,
        resolveTreeIntoPatientCladesChris(tree, FALSE, num.tips, child))
      }
    }
  }
  return(clades.at.this.level)
}



findOnePatientsClades <-
function(tree, id, first.call=TRUE, num.tips=NULL, node=NULL) {
  # From a tree, returns a list in which each entry is a monophyletic set of
  # tips for the named patient, collected into a vector.

  # This function calls itself iteratively. If this is the first call, set the 
  # node to be the patient's mrca and find the number of tips (used repeatedly).
  if (first.call) {
    their.tips <- patient.tips[[id]]
    node <- findMRCA(tree, their.tips, type="node")
    num.tips <- length(tree$tip.label)
  }

	clades.at.this.level <- list()

  # Are all the tips that descend from this node from this patient?
  descendant.tips <- Descendants(tree, node=node, type="tips")[[1]]
  all.this.patient <- TRUE
  for (tip.number in descendant.tips) {
    if (getIdFromTip(tree$tip.label[tip.number]) != id) {
      all.this.patient <- FALSE
      break
    }
  }
  if (all.this.patient) {
    clades.at.this.level[[1]] <- tree$tip.label[descendant.tips]
  } else {
    # Iterate through the immediate children of this node.
    children <- Descendants(tree, node=node, type="children")
    for (child in children) {
      if (child <= num.tips) {
        # The child is a tip.
        # If it's a tip from our patient, add it to our list.
        tip.name <- tree$tip.label[child]
        if (getIdFromTip(tip.name, id.delimter) == id) {
          clades.at.this.level[length(clades.at.this.level)+1] <- c(tip.name)
        }
      } else {
        # The child is a node. Merge clades.at.this.level with the result of
        # applying this function to the child.
        clades.at.this.level <- c(clades.at.this.level,
        findOnePatientsClades(tree, id, FALSE, num.tips, child))
      }
    }
  }
  return(clades.at.this.level)
}




analyseCladesChristophe <- function(id, patient.tips, tree,
test.clade.ordering.only) {
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

    # TODO: auto-indent this function
		num.tips <- length(patient.tips[[id]])
		read.counts.per.tip <- list()
		for (tip in patient.tips[[id]]) {
      read.counts.per.tip[[tip]] <- as.numeric(unlist(strsplit(tip,
      "count_"))[2])
    }
		read.counts.per.tip <- unlist(read.counts.per.tip)
		num.reads.total <- sum(read.counts.per.tip)

		ordered.clades <- list()
		current.clade <- 1
		if (num.tips > 0) {
			if (num.tips == 1) {
				ordered.clades[[current.clade]] <- makeSingleTipTree(tip,
        read.counts.per.tip)
				overall.root.to.tip <- 0
			} else {
				# define the subtree of all tips associated with the patient
				patient.subtree <- drop.tip(phy = tree,
						tip = tree$tip.label[!(tree$tip.label %in% patient.tips[[id]])])
				patient.clade <- list(tree = patient.subtree, is.a.tip = F,
        num.reads = num.reads.total)
				overall.root.to.tip <- calcMeanRootToTipChristophe(patient.clade,
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
							calcMeanRootToTipChristophe(ordered.clades[[i]], read.counts.per.tip)
				}
			}

    # Order each clade by the read number of the tips, and the set of clades by
    # their smallest read number. This is unique (unlike read count sorting) so
    # helps test clade finding algorithms.
    if (test.clade.ordering.only) {
      clades.as.list <- list()
      current.clade <- 0
      for (clade in ordered.clades) {
        current.clade <- current.clade + 1
        if (clade$is.a.tip) {
          clades.as.list[[current.clade]] <- clade$tree
        } else {
          clades.as.list[[current.clade]] <- clade$tree$tip.label
        }
      }
      sorted.clades <- printSortedClades(clades.as.list, tip.regex,
      num.clades.for.output)
      print(id)
      print(sorted.clades)
      cat('\n\n')
    }

		} else { # if there are no reads for this patient in this window
			num.clades <- NA
			ordered.clades <- NA
			overall.root.to.tip <- NA
			prop.reads.per.clade <- rep(NA, 5)
			root.to.tip.per.clade <- rep(NA, 5)
		}

  return(list(num.reads.total=num.reads.total, num.tips=num.tips,
  num.clades=num.clades, overall.root.to.tip=overall.root.to.tip,
  prop.reads.per.clade=prop.reads.per.clade,
  root.to.tip.per.clade=root.to.tip.per.clade))
}

