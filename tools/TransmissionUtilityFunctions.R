suppressMessages(require(phangorn, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(gdata, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(prodlim, quietly=TRUE, warn.conflicts=FALSE))

# Get the tip number of this label

get.tip.no <- function(tree, name) {
  return(match(name, tree$tip.label))
}

# Get all tips associated with a patient

get.tips.for.patient <-
  function(tree, patient.string, patient.ids, blacklist) {
    node.numbers <- which(startsWith(patient.ids, patient.string))
    
    node.numbers <-
      node.numbers[which(!(node.numbers %in% blacklist))]
    return(node.numbers)
  }

get.tips.for.sample <-
  function(tree, sample.string){
    return(which(startsWith(tree$tip.label, sample.string)))
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

# This returns the type of intermingling, according to a modified Romero-Severson
# classification

intermingled.class <-
  function(tree, pat.1, pat.2, mrca.list, tip.list) {
    
    if(!are.intermingled(tree, pat.1, pat.2, mrca.list, tip.list)){
      return(NULL)
    }
    
    keep.tips <- c(tip.list[[pat.1]], tip.list[[pat.2]])
    
    keep.labels.1 <- tree$tip.label[tip.list[[pat.1]]]
    keep.labels.2 <- tree$tip.label[tip.list[[pat.2]]]
    
    small.tree <- drop.tip(phy = tree, tip = tree$tip.label[-keep.tips])
    
    root <- getRoot(small.tree)
    
    return(classify(small.tree, root, pat.1, keep.labels.1, pat.2, keep.labels.2 ))
    
  }

# Classify this node from two patients with tips descending from it according to the 
# R-S rules (ignoring any other patients). Will return "*" if equivocal, the patient
# otherwise

classify <- function(tree, node, pat.1, labels.1, pat.2, labels.2){
  
  if(length(tree$tip.label)!=length(labels.1)+length(labels.2)){
    stop("This subtree hasn't been correctly pruned")
  }
  
  if(is.tip(tree, node)){
    if(tree$tip.label[node] %in% labels.1){
      return(pat.1)
    } else {
      return(pat.2)
    }
  } else {
    child.classes <- vector()
    
    for(child in Descendants(tree, node, type="children")){
      child.classes <- c(child.classes, classify(tree, child, pat.1, labels.1, pat.2, labels.2))
    }
    no.pat.1 <- sum(child.classes==pat.1)
    no.pat.2 <- sum(child.classes==pat.2)
    no.unclear <- sum(child.classes=="*")
    if(no.pat.1 == no.pat.2){
      return("*")
    } else if(no.pat.1 > no.pat.2) {
      return(pat.1)
    } else {
      return(pat.2)
    }
  }
  
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

# This checks if the MRCA of the tips from patient desc is a descendant of the 
# MRCA of the tips from patient anc, and that there is no node on the route from one to
# the other that is an ancestor of two tips from a third patient (which would suggest)
# that that third patient was a transmission intermediary. If zero.length.tips.count
# is true, then a node with a child separated from itself by a zero length branch
# (or one of length 1E-5, which is how RAxML represents branches indicating no
# mutation) is also taken to indicate a sampled intermediary.


is.direct.descendant.of <-
  function(tree, desc, anc, mrca.list, tip.list, zero.length.tips.count = TRUE) {
    if (verbose) {
      cat("Testing if any ancestry exists\n")
    }
    
    if (!is.descendant.of(tree, desc, anc, mrca.list, tip.list)) {
      return(FALSE)
    }
    if (verbose) {
      cat("Getting MRCA nodes involved\n")
    }
    
    mrca.desc <- mrca.list[[desc]]
    mrca.anc <- mrca.list[[anc]]
    
    
    # Check that there are no nodes with two children with descendent tips belonging to 
    # the same host (other than anc and desc) in the sequence of ancestors from desc to 
    # anc (taking polytomies into account) - this means that either an intervening sampled 
    # host cannot be ruled out, or one or the other MRCA is shared with another patient
    # and hence the situation is too ambiguous
    if (verbose) {
      cat("Getting ancestral sequence\n")
    }
    
    anc.sequence <-
      get.ancestral.sequence(tree, mrca.desc, mrca.anc)
    
    if (verbose) {
      cat("Checking for intervening hosts\n")
    }
    
    return(test.direct.ancestry(tree, anc.sequence, c(anc, desc), zero.length.tips.count))
  }

# Sub-method of the above; takes a sequence of nodes that represents a chain
# of ancestors and tests if it looks like there were no intervening patients

test.direct.ancestry <-
  function(tree, sequence, exceptions, zero.length.tips.count){
    for (node in sequence) {
      if (!is.tip(tree, node)) {
        children <- Children(tree, node)
        
        all.descendant.hosts <-
          unlist(sapply(children, get.patients.from.this.clade, tree = tree))
        
        # all.descendant.hosts will contain duplicates if any patient is attached 
        # to a descendant of more than one element of children
        
        duplicates <-
          unique(all.descendant.hosts[duplicated(all.descendant.hosts)])
        
        # exceptions are a list of patients whose presence can be discounted; anc
        # and desc should normally be exceptions 
        
        duplicates <-
          duplicates[!duplicates %in% exceptions]
        
        if (length(duplicates) > 0) {
          # there is at least one intervening node with two descendants from the same patient
          return(FALSE)
        }
        
        if (zero.length.tips.count) {
          # if zero length tips count, then a tip attached to an internal node on the ancestral sequence by a branch
          # of length zero counts as an ancestor
          for (child in Children(tree, node)) {
            if (is.tip(tree, child) & get.edge.length(tree, child) < 1E-5) {
              if (!(get.patient.from.tip(tree, child) %in% exceptions)) {
                return(FALSE)
                if (verbose) {
                  cat("Returning FALSE due to zero-length interruption")
                }
              }
            }
          }
        }
      }
    }
    return(TRUE)
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
  function(tree, node, zero.length.tips.count = TRUE) {
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

# Splits the patient given by id up until the subtree derived from its tips only. Currently redundant.
# has no branches longer than threshold. 
# The following are modified by this procedure:
# patients - vector of patient names
# patient.tips - list of vectors of tips indexed by patient
# patient.mrcas - list of nodes index by patient
# patient.ids - the patient that each tip is assigned to

split.patient <-
  function(tree, id, split.threshold, patients, patient.tips, patient.mrcas, patient.ids, 
           split.counts, original.id = id, already.split = FALSE) {
    
    if (is.null(original.id)) {
      original.id = id
    }
    
    if (!already.split) {
      #this patient has not previously been split
      split.counts[[id]] <- 0
    }
    
    
    if (length(patient.tips[[id]]) <= 1) {
      if (verbose) {
        cat(id,": single tip, no split possible\n", sep = "")
      }
    } else {
      patient.subtree <-
        drop.tip(phy = tree, tip = tree$tip.label[-patient.tips[[id]]])
      if (verbose) {
        cat(length(patient.subtree$tip.label)," tips\n",sep = "")
      }
      
      #can't unroot a tree with less than three edges (a bit weirdly)
      
      if (length(patient.subtree$edge.length) >= 3) {
        unrooted.patient.subtree <- unroot(patient.subtree)
        max.branch.length <- max(unrooted.patient.subtree$edge.length)
      } else {
        max.branch.length <- sum(patient.subtree$edge.length)
      }
      
      if (max.branch.length > split.threshold) {
        
        if (length(patient.subtree$edge.length) >= 3) {
          edge.index <-
            which(unrooted.patient.subtree$edge.length == max.branch.length)
          
          start.node.index <-
            unrooted.patient.subtree$edge[edge.index,2]
          end.node.index <-
            unrooted.patient.subtree$edge[edge.index,1]
          
          split.set <- list()
          split.set[[1]] <- vector()
          split.set[[2]] <- vector()
          
          for (tip in seq(1,length(unrooted.patient.subtree$tip.label))) {
            current.node <- tip
            
            while (!(current.node %in% c(start.node.index, end.node.index))) {
              current.node <- next.node(unrooted.patient.subtree, current.node)
            }
            if (current.node == start.node.index) {
              split.set[[1]] <-
                c(split.set[[1]], unrooted.patient.subtree$tip.label[tip])
            } else {
              split.set[[2]] <-
                c(split.set[[2]], unrooted.patient.subtree$tip.label[tip])
            }
            
          }
        } else {
          split.set <- list()
          split.set[[1]] <- patient.subtree$tip.label[1]
          split.set[[2]] <- patient.subtree$tip.label[2]
        }
        
        patients <- patients[patients != id]
        patient.tips[[id]] <- NULL
        patient.mrcas[[id]] <- NULL
        
        if (already.split) {
          #this gets the split count right, but is a bit of a hack
          split.counts[[original.id]] <- split.counts[[original.id]] - 1
        }
        
        for (tip.group in split.set) {
          split.counts[[original.id]] <- split.counts[[original.id]] + 1
          
          new.name <-
            paste(original.id,"_split",split.counts[[original.id]], sep = "")
          if (verbose) {
            cat("Child ",split.counts[[original.id]],": ", new.name,"\n", sep = "")
          }
          
          patients <- c(patients, new.name)
          
          tip.group.in.full.tree <-
            which(tree$tip.label %in% tip.group)
          
          patient.tips[[new.name]] <- tip.group.in.full.tree
          patient.mrcas[[new.name]] <-
            mrca.phylo.or.unique.tip(tree, tip.group.in.full.tree, zero.length.tips.count)
          
          for (tip.no in patient.tips[[new.name]]) {
            patient.ids[tip.no] <- new.name
          }
          
          split.patient(
            tree, new.name, split.threshold, patients, patient.tips, patient.mrcas, patient.ids,
            split.counts, original.id, already.split = TRUE
          )
        }
      }
      return(list(patient.names = patients, patient.tips = patient.tips, patient.mrcas = patient.mrcas,
                  patients.for.tips <- patient.ids))
    }
    
    # A function is required to traverse the tree. This is surprisingly difficult for unrooted 
    # trees; as a hack, pick an adjacent node randomly; the full tree will be explored in the end.
    
    # todo make this better
    
    next.node <- function(tree, node) {
      connections.1 <- tree$edge[which(tree$edge[,2] == node),1]
      connections.2 <- tree$edge[which(tree$edge[,1] == node),2]
      
      all.connections <-
        c(connections.1, connections.2)[!(is.na(c(connections.1, connections.2)))]
      
      return(all.connections[sample.int(length(all.connections), 1)])
    }
    
    # Does the full R-S classification from this node downwards
    
    # todo may be better off doing this for LikelyTransmissions
    
    classify.down <- function(tree, node, blacklist){
      print(node)
      if(node %in% blacklist){
        assignments[[node]] <<- "blacklisted"
        print(paste(node,": node is blacklisted",sep=""))
        return("blacklisted")
      }
      if(is.tip(tree, node)){
        assignments[[node]] <<- get.patient.from.tip(tree, node)
        print(paste(node,": node is a tip",sep=""))
        return(get.patient.from.tip(tree, node))
      }
      print(paste(node,": node is internal", sep=""))
      child.classifications <- sapply(Children(tree, node), function(node) classify.down(tree, node, blacklist))
      child.classifications <- child.classifications[which(child.classifications!="blacklisted")]
      if(length(child.classifications)==0){
        assignments[[node]] <<- "blacklisted"
        print(paste(node,": all children blacklisted", sep=""))
        return("blacklisted")
      }
      not.equivocal <- child.classifications[which(lengths(child.classifications)==1)]
      if(length(not.equivocal)==0){
        first <- TRUE
        for(pat.options in child.classifications){
          if(first){
            pat.intersection <- pat.options
            pat.union <- pat.options
            first <- FALSE
          } else {
            pat.intersection <- intersect(pat.intersection, pat.options)
            pat.union <- union(pat.union, pat.options)
          }
        }
        if(length(pat.intersection)==1){
          print(paste(node,": all children equivocal but one association with every child", sep=""))
          assignments[[node]] <<- pat.intersection[1]
          return(intersection[1])
        }
        print(paste(node,": all children equivocal and not unambiguous", sep=""))
        assignments[[node]] <<- unique(pat.union)
        return(unique(pat.union))
      }
      child.classifications <- not.equivocal
      counts <- as.data.frame(table(unlist(child.classifications)))
      max.counts <- max(counts$Freq)
      if(length(counts$Freq[which(counts$Freq==max.counts)])>1){
        assignments[node] <<- as.character(counts$Var1[which(counts$Freq==max.counts)])
        print(paste(node,": ",length(counts$Freq[which(counts$Freq==max.counts)])," patients appear the same number of times in children",sep=""))
        return(as.character(counts$Var1[which(counts$Freq==max.counts)]))
      } else {
        print(paste(node,": winner is ",as.character(counts$Var1[which(counts$Freq==max.counts)[1]]),sep=""))
        assignments[node] <<- as.character(counts$Var1[which(counts$Freq==max.counts)[1]])
        return(as.character(counts$Var1[which(counts$Freq==max.counts)[1]]))
      }
    }
  }

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
    return(sub(regexp, "\\3", label))
  } else {
    return(NA)
  }
}

# Transmission tree navigation

get.tt.parent <- function(tt, label){
  if(length(which(tt$unique.splits==label))>0){
    return(tt$parent.splits[which(tt$unique.splits==label)])
  } else {
    return(NULL)
  }
}

get.tt.ancestors <- function(tt, label){
  out <- vector()
  current.node <- get.tt.parent(tt, label)
  while(current.node != "root"){
    out <- c(out, current.node)
    current.node <- get.tt.parent(tt, current.node)
  }
  return(out)
}

get.tt.children <- function(tt, label){
  return(tt$unique.splits[which(tt$parent.splits==label)])
}

get.tt.parent.patient <- function(tt, label){
  if(length(which(tt$unique.splits==label))>0){
    return(tt$patients[which(tt$unique.splits==label)])
  } else {
    return(NULL)
  }
}

get.tt.mrca <- function(tt, label1, label2){
  # sanity check
  ca.vec <- intersect(c(label1, get.tt.ancestors(tt, label1)), c(label2, get.tt.ancestors(tt, label2)))
  if(length(ca.vec)>1){
    for(i in seq(1, length(ca.vec)-1)){
      if(get.tt.parent(tt, ca.vec[i]) != ca.vec[i+1]){
        stop("get.tt.mrca not working as intended")
      }
    }
  }
  
  return(ca.vec[1])
}

get.tt.mrca.patient <- function(tt, label1, label2){
   node <- intersect(c(label1,get.tt.ancestors(tt, label1)), c(label2, get.tt.ancestors(tt, label2)))[1]
   
   return(tt$patients[which(tt$unique.splits==node)])
}

get.tt.path <- function(tt, label1, label2){
  mrca <- get.tt.mrca(tt, label1, label2)
  first.half <- vector()
  current.node <- label1
  while(current.node != mrca){
    first.half <- c(first.half, current.node)
    current.node <- get.tt.parent(tt, current.node)
  }
  
  second.half <- vector()
  current.node <- label2
  while(current.node != mrca){
    second.half <- c(second.half, current.node)
    current.node <- get.tt.parent(tt, current.node)
  }
  
  return(c(first.half, mrca, rev(second.half)))
}

# starting at node 1, determine whether the path to node 2 is blocked by node 3
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
  
  return(sum(which(edges > dist.to.root))/length(edges))
}

