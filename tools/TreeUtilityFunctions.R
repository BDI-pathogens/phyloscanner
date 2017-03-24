# list.files seems to behave differently on different systems

list.files.mod <- function(path = ".", pattern = NULL, all.files = FALSE,
                          full.names = FALSE, recursive = FALSE,
                          ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE){
  
  original <- list.files(path, pattern, all.files, full.names, recursive, ignore.case, include.dirs, no..)
  return(sapply(original, function(x) if(substr(x, 1, 2)=="./") {substr(x, 3, nchar(x))} else {x}))
}


# Get the tip number of this label

get.tip.no <- function(tree, name) {
  return(match(name, tree$tip.label))
}

# Get all tips associated with a patient

get.tips.for.patient <-
  function(tree, patient.string, patient.ids, blacklist) {
    node.numbers <- which(grepl(paste0('^',patient.string),patient.ids))
    
    node.numbers <-
      node.numbers[which(!(node.numbers %in% blacklist))]
    return(node.numbers)
  }

get.tips.for.sample <-
  function(tree, sample.string){
    return(which(grepl(paste0('^',sample.string),tree$tip.label)))
  }

# Edge length by node number (probably in some library somewhere)

get.edge.length <- function(tree, node) {
  index <- which(tree$edge[,2] == node)
  return(tree$edge.length[index])
}

# Node is a tip (likewise)

is.tip <- function(tree, node) {
  return(node <= length(tree$tip.label))
}

# Node is root (likewise)

is.root <- function(tree, node) {
  return(Ancestors(tree, node, type = "parent") == 0)
}

# Get the sequence of ancestors from one node to the other

get.ancestral.sequence <- function(tree, desc, anc) {
  if (!(anc %in% Ancestors(tree, desc, type = "all"))) {
    stop("anc is not an ancestor of desc at all")
  }
  current <- desc
  out <- desc
  
  while (current != anc) {
    current <- Ancestors(tree, current, type = "parent")
    
    out <- c(out, current)
    
  }
  return(out)
}

# Output the patient that this tip belongs to. The vector patient.ids 
# records patients in the same order as the tips in the tree

get.patient.from.tip <- function(tree, node, patient.ids) {
  # is it a tip?
  if (!(is.tip(tree, node))) {
    stop("Not a tip")
  }
  patient.id <- patient.ids[node]
  return(patient.id)
  
}

# Output the set of patients whose tips are descended from a node 
# (or the host of the node itself if a tip)

get.patients.from.this.clade <- function(tree, node, patient.ids) {
  if (is.tip(tree,node)) {
    return(get.patient.from.tip(tree, node, patient.ids))
  }
  
  out <- vector()
  tips <- Descendants(tree, node, type = "tips")[[1]]
  for (tip in tips) {
    if (!(get.patient.from.tip(tree,tip) %in% out)) {
      out <- c(out, get.patient.from.tip(tree, tip, patient.ids))
    }
  }
  return(out)
}

# Is desc _unambiguously_ a descendant of anc? I.e. is the MRCA node of desc 
# a descendant of the MRCA node of anc, but no tips of anc are descended from
# the MRCA node of desc?
# mrca.list is a list of MRCA nodes indexed by patient
# tip.list is a list of vectors of tips indexed by patient

is.descendant.of <- function(tree, desc, anc, mrca.list, tip.list) {
  mrca.desc <- mrca.list[[desc]]
  mrca.anc <- mrca.list[[anc]]
  
  if (!(mrca.anc %in% Ancestors(tree, mrca.desc, type = "all"))) {
    return(FALSE)
  } else {
    for (tip in tip.list[[anc]]) {
      if (mrca.desc %in% Ancestors(tree, tip, type = "all")) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

# Is the intersection between vectors x and y empty?

nonempty.intersection <- function(x,y) {
  return(length(intersect(x,y)) > 0)
}

# Returns whether the tips corresponding to patients pat.1 and pat.2 are intermingled
# in tree, i.e. if the MRCA of one set is both descended from (or equal to) the MRCA of 
# the other and an ancestor of at least one tip from the other
# mrca.list is a list of MRCA nodes indexed by patient
# tip.list is a list of vectors of tips indexed by patient

are.intermingled <-
  function(tree, pat.1, pat.2, mrca.list, tip.list) {
    mrca.1 <- mrca.list[[pat.1]]
    mrca.2 <- mrca.list[[pat.2]]
    
    if (mrca.1 == mrca.2) {
      return(TRUE)
    }
    
    if ((mrca.1 %in% Ancestors(tree, mrca.2, type = "all"))) {
      for (tip in tip.list[[pat.1]]) {
        if (mrca.2 %in% Ancestors(tree, tip, type = "all")) {
          return(TRUE)
        }
      }
    }
    
    if ((mrca.2 %in% Ancestors(tree, mrca.1, type = "all"))) {
      for (tip in tip.list[[pat.2]]) {
        if (mrca.1 %in% Ancestors(tree, tip, type = "all")) {
          return(TRUE)
        }
      }
    }
    
    return(FALSE)
  }

# Test if two patients form sibling clades (MM in Romero-Severson terminology). This also
# requires that there is no third patient whose position in the tree suggests that it was
# either a common source of both infections, or that the transmission chain passed through
# it to get from the common source to either patient.

are.siblings <-
  function(tree, pat.1, pat.2, mrca.list, tip.list, zero.length.tips.count =
             TRUE) {
    mrca.1 <- mrca.list[[pat.1]]
    mrca.2 <- mrca.list[[pat.2]]
    
    if (mrca.1 == mrca.2) {
      return(F)
    }
    if (mrca.1 %in% Ancestors(tree, mrca.2, type = "all")) {
      return(F)
    }
    if (mrca.2 %in% Ancestors(tree, mrca.1, type = "all")) {
      return(F)
    }
    
    mrca.mrca <- mrca.phylo(tree, c(mrca.1, mrca.2))
    
    anc.sequence.1 <- get.ancestral.sequence(tree, mrca.1, mrca.mrca)
    anc.sequence.2 <- get.ancestral.sequence(tree, mrca.2, mrca.mrca)
    
    return(
      test.direct.ancestry(tree, anc.sequence.1, pat.1, zero.length.tips.count) &
        test.direct.ancestry(tree, anc.sequence.2, pat.2, zero.length.tips.count)
    )
    
  }

# This returns the distance along the branches from the TMRCA of one patient to the TMRCA
# of the other, if the patients are siblings

sibling.distance <- function(tree, pat.1, pat.2, mrca.list, tip.list, zero.length.tips.count){
  
  mrca.1 <- mrca.list[[pat.1]]
  mrca.2 <- mrca.list[[pat.2]]
  
  if (!are.siblings(tree, pat.1, pat.2, mrca.list, tip.list, zero.length.tips.count)){
    stop("Not siblings")
  }
  
  mrca.mrca = mrca.phylo(tree, c(mrca.1, mrca.2))
  
  dist.1 <- 0
  current.node <- mrca.1
  while(current.node != mrca.mrca){
    dist.1 <- dist.1 + get.edge.length(tree, current.node)
    current.node <- Ancestors(tree, current.node, type="parent")
  }
  
  dist.2 <- 0
  current.node <- mrca.2
  while(current.node != mrca.mrca){
    dist.2 <- dist.2 + get.edge.length(tree, current.node)
    current.node <- Ancestors(tree, current.node, type="parent")
  }
  
  return(dist.1+dist.2)
  
}

# Output the set of patients whose mrca is this node

get.patients.with.these.mrcas <-
  function(tree, node, patient.mrcas) {
    out <- vector()
    mrca.vec <- sapply(patient.mrcas, "[[", 1)
    
    numbers <- which(mrca.vec == node)
    
    return(names(numbers))
  }

# Get the distance between the MRCA nodes of these two patients

get.mrca.distance <- function(tree, pat.1, pat.2, mrca.list) {
  mrca.1 <- mrca.list[[pat.1]]
  mrca.2 <- mrca.list[[pat.2]]
  
  return(total.length.between(tree, mrca.1, mrca.2))
  
}

# Get the distance between these two nodes

total.length.between <- function(tree, node.1, node.2) {
  if (node.1 == node.2) {
    return(0)
  }
  total.length <- 0
  if (node.1 %in% Descendants(tree, node.2, type = "all")) {
    current.node <- node.1
    while (current.node != node.2) {
      total.length <- total.length + get.edge.length(tree, current.node)
      current.node <- Ancestors(tree, current.node, type = "parent")
    }
    return(total.length)
    
  } else if (node.2 %in% Descendants(tree, node.1, type = "all")) {
    current.node <- node.2
    while (current.node != node.1) {
      total.length <- total.length + get.edge.length(tree, current.node)
      current.node <- Ancestors(tree, current.node, type = "parent")
    }
    return(total.length)
    
  } else {
    joint.mrca <- mrca.phylo(tree, c(node.1, node.2))
    return(
      total.length.between(tree, node.1, joint.mrca) + total.length.between(tree, node.2, joint.mrca)
    )
  }
  
}

# mrca.phylo applied to a tip only will not return that tip. This function will.

mrca.phylo.or.unique.tip <-
  function(tree, node, zero.length.tips.count = FALSE) {
    if (length(node) == 1) {
      if (!zero.length.tips.count | !is.tip(tree,node)) {
        return(node)
      } else {
        length <- get.edge.length(tree, node)
        while (length < 1E-5) {
          node <- Ancestors(tree, node, type = "parent")
          if (is.root(tree, node)) {
            break
          }
          length <- get.edge.length(tree, node)
          
        }
        return(node)
      }
    } else {
      mrca <- mrca.phylo(tree, node)
      
      return(mrca)
      
    }
  }

# This function returns TRUE wherever elements are the same, including NA's, and false 
# everywhere else.

compareNA <- function(v1,v2) {
  same <- (v1 == v2)  |  (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}

# Last element in a vector

get.last <- function(vec)
  return(vec[length(vec)])

# Get patient id from tip label

patient.from.label <- function(label, regexp){
  if(length(grep(regexp, label)>0)) {
    return(sub(regexp, "\\1", label))
  } else {
    return(NA)
  }
}

# Get read count from label

read.count.from.label <- function(label, regexp){
  if(length(grep(regexp, label)>0)) {
    return(as.numeric(sub(regexp, "\\3", label)))
  } else {
    return(NA)
  }
}

# Starting at node 1, determine whether the path to node 2 is blocked by node 3
# Should work on both rooted and unrooted trees I _hope_.

path.exists <- function(tree, node.1, node.2, node.3, last.node = -1){
  #  cat(node.1," ")
  if(node.1 == node.2){
    return(TRUE)
  }
  if(node.1 == node.3){
    return(FALSE)
  }
  neighbours = vector()
  if(!is.tip(tree, node.1)){
    neighbours <- Children(tree, node.1)
  }
  if(length(Ancestors(tree, node.1, type="parent"))!=0){
    if(Ancestors(tree, node.1, type="parent")!=0){
      neighbours <- c(neighbours, Ancestors(tree, node.1, type="parent"))
    }
  }
  for(neighbour in neighbours){
    if(neighbour != last.node){
      if(path.exists(tree, neighbour, node.2, node.3, node.1)){
        return(TRUE)
        break
      }
    }
  }
  #you've been everywhere and you can't find a way through
  return(FALSE)
}

tips.reachable <- function(tree, node, blocking.edges, last.node = -1){
  out <- vector()
  if(is.tip(tree, node)){
    out <- c(out, node)
  }
  neighbours = neighbouring.nodes(tree, node)
  for(neighbour in neighbours){
    if(neighbour != last.node){
      edge.vec <- c(node, neighbour)
      if(is.na(row.match(edge.vec, t(as.matrix(blocking.edges)))) &  is.na(row.match(rev(edge.vec), t(as.matrix(blocking.edges))))){
        out <- c(out, tips.reachable(tree, neighbour, blocking.edges, node))
      }
    }
  }
  return(out)
}

neighbouring.nodes <- function(tree, node){
  part.1 <- tree$edge[which(tree$edge[,2]==node), 1]
  part.2 <- tree$edge[which(tree$edge[,1]==node), 2]
  return(c(part.1, part.2))
}


# Drop a set of tips and return a vector which maps nodes from the full tree to the subtree

drop.tip.get.map <- function(phy, tip){
  if(length(unique(tree$tip.label))!=length(tree$tip.label)){
    stop("This won't work if there are duplicate tip names")
  }
  reference <- vector()
  
  phy.2 <- drop.tip(phy, tip)
  
  for(new.label in seq(1, length(phy.2$tip.label))){
    reference[new.label] = which(phy$tip.label==phy.2$tip.label[new.label])
  }
  
  mrcas.1 <- mrca(phy)
  mrcas.2 <- mrca(phy.2)
  
  for(tip.1 in seq(1, length(phy.2$tip.label))){
    for(tip.2 in seq(1, length(phy.2$tip.label))){
      if(tip.1 < tip.2){
        new.mrca <- mrcas.2[tip.1,tip.2]
        old.mrca <- mrcas.1[reference[tip.1], reference[tip.2]]
        reference[new.mrca] = old.mrca
      }
    }
  }
  return(list(tree = phy.2, reference=reference))
}

extract.subtrees.for.patients <- function(tree, patients, splits){
  labels.to.keep <- splits$tip.names[which(splits$orig.patients %in% patients)]
  
  pruned.tree <- drop.tip(tree, which(!(tree$tip.label %in% labels.to.keep)))
  
  return(pruned.tree)
}

prop.internal.longer.than.root <- function(tree, split, splits){
  tips <- splits$tip.names[which(splits$patient.splits==split)]
  
  if(length(tips)==1){
    return(0)
  }
  
  split.root <- mrca.phylo.or.unique.tip(tree, which(tree$tip.label %in% tips))
  dist.to.root <- 0 
  
  current.node <- split.root
  while(!is.root(tree, current.node)){
    dist.to.root <- dist.to.root + get.edge.length(tree, split.root)
    current.node <- Ancestors(tree, current.node, type="parent")
  }
  
  just.clade <- drop.tip(tree, which(!(tree$tip.label %in% tips)))
  
  edges <- just.clade$edge.length
  
  return(sum(edges > dist.to.root)/length(edges))
}