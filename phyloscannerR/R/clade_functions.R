#' @keywords internal
#' @export resolve.tree.into.host.clades

resolve.tree.into.host.clades <- function(tree, ids, tip.regex, blacklisted.tips = vector(), no.read.counts = T) {
  
  num.tips <- length(tree$tip.label)
  num.patients <- length(ids)

  # Reorder the nodes in such a manner that we always encounter a descendant 
  # node before its ancestor (not a unique process but that doesn't matter).
  tree <- reorder(tree, "postorder")

  # Make a vector of bools describing which tips are patients. Remove blacklisted tips from this vector.
  is.patient <- grep(tip.regex, tree$tip.label)
  is.patient <- is.patient[which(!(is.patient %in% blacklisted.tips))]
  
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
                       sum(colSums(m[these.children, , drop=F]) > 0) == 1)

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

  return(list(clades.by.patient=clades.by.patient, clade.mrcas.by.patient=clade.mrcas.by.patient))
}

#' @keywords internal
#' @export calc.mean.root.to.tip

calc.mean.root.to.tip <- function(tree, read.counts.per.tip) {
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