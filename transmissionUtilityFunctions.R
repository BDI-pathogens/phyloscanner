require(phangorn)

get.tip.no <- function(tree, name) {
  return(match(name, tree$tip.label))
}

# Get all tips associated with a patient

get.tips.for.patient <- function(tree, patient.string, patient.ids, blacklist) {
  node.numbers <- which(patient.ids == patient.string)
  
  node.numbers <- node.numbers[which(!(node.numbers %in% blacklist))]
  return(node.numbers)
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
  return(Ancestors(tree, node, type="parent")==0)
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

get.patient.from.tip <- function(tree, node) {
  # is it a tip?
  if (!(is.tip(tree, node))) {
    stop("Not a tip")
  }
  patient.id <- patient.ids[node]
  return(patient.id)
  
}

# Output the set of patients whose tips are descended from node (or the host of the node itself if a tip)

get.patients.from.this.clade <- function(tree, node, patient.tips) {
  if (is.tip(tree,node)) {
    return(get.patient.from.tip(tree, node))
  }
  
  out <- vector()
  tips <- Descendants(tree, node, type = "tips")[[1]]
  for (tip in tips) {
    if (!(get.patient.from.tip(tree,tip) %in% out)) {
      out <- c(out, get.patient.from.tip(tree, tip))
    }
  }
  return(out)
}

# Is desc _unambiguously_ a descendant of anc? I.e. is the MRCA node of desc a descendant of the MRCA
# node of anc, but no tips of anc are descended from the MRCA node of desc?

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

nonempty.intersection <- function(x,y) {
  return(length(intersect(x,y)) > 0)
}

is.direct.descendant.of <-
  function(tree, desc, anc, mrca.list, tip.list, zero.length.tips.count=TRUE) {
    if (verbose) {
      cat("Testing if any ancestry exists\n")
    }
    
    if (!is.descendant.of(tree,desc,anc,mrca.list,tip.list)) {
      return(FALSE)
    }
    if (verbose) {
      cat("Getting MRCA nodes involved\n")
    }
    
    mrca.desc <- mrca.list[[desc]]
    mrca.anc <- mrca.list[[anc]]
    
    
    #check that there are no nodes with two children with descendent tips belonging to the same host (other than
    #anc and desc) in the sequence of ancestors from desc to anc (taking polytomies into account) - this means that
    #either an intervening sampled host cannot be ruled out, or one or the other MRCA is shared with another patient
    #and hence the situation is too ambiguous
    if (verbose) {
      cat("Getting ancestral sequence\n")
    }
    
    anc.sequence <- get.ancestral.sequence(tree, mrca.desc, mrca.anc)
    
    if (verbose) {
      cat("Checking for intervening hosts\n")
    }
    
    for (node in anc.sequence) {
      if (!is.tip(tree, node)) {
        children <- Children(tree, node)
        
        all.descendant.hosts <-
          unlist(sapply(children, get.patients.from.this.clade, tree = tree))
        
        # all.descendant.hosts will contain duplicates if any patient is attached to a descendant of more than
        # one element of children
        
        duplicates <-
          unique(all.descendant.hosts[duplicated(all.descendant.hosts)])
        
        # anc or desc being in the set is OK
        
        duplicates <-
          duplicates[!duplicates %in% c(patients[anc], patients[desc])]
        
        if (length(duplicates) > 0) {
          if(verbose){
            cat("Returning FALSE due to problems on the sequence (",node,")\n")
          }
          return(FALSE)
          
          
        }
        
        if (zero.length.tips.count) {
          # if zero length tips count, then a tip attached to an internal node on the ancestral sequence by a branch
          # of length zero counts as an ancestor
          for (child in Children(tree, node)) {
            if (is.tip(tree, child) & get.edge.length(tree, child) < 1E-5) {
              if (!(get.patient.from.tip(tree, child) %in% c(patients[anc], patients[desc]))) {
                return(FALSE)
                if(verbose){
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

get.mrca.distance <- function(tree, pat.1, pat.2, mrca.list){
  mrca.1 <- mrca.list[[pat.1]]
  mrca.2 <- mrca.list[[pat.2]]
  
  return(total.length.between(tree, mrca.1, mrca.2))
  
}

total.length.between <- function(tree, node.1, node.2){
  if(verbose) cat("Calculating length between nodes ",node.1," and ",node.2,"\n",sep="")
  if(node.1==node.2){
    return(0)
  }
  total.length <- 0
  if(node.1 %in% Descendants(tree, node.2, type="all")){
    if(verbose) cat(node.1," descends from ",node.2,"\n",sep="")
    current.node <- node.1 
    while(current.node!=node.2){
      total.length <- total.length + get.edge.length(tree, current.node)
      current.node <- Ancestors(tree, current.node, type="parent")
    }
    return(total.length)
    
  } else if(node.2 %in% Descendants(tree, node.1, type="all")){
    if(verbose) cat(node.2," descends from ",node.1,"\n",sep="")
    current.node <- node.2 
    while(current.node!=node.1){
      total.length <- total.length + get.edge.length(tree, current.node)
      current.node <- Ancestors(tree, current.node, type="parent")
    }
    return(total.length)
    
  } else {
    if(verbose) cat(node.2," and ",node.1," are not directly related\n",sep="")
    joint.mrca <- mrca.phylo(tree, c(node.1, node.2))
    return(total.length.between(tree, node.1, joint.mrca) + total.length.between(tree, node.2, joint.mrca))
  }
  
}

# Need a MRCA function which if given a single tip, returns that tip rather than NA.

mrca.phylo.or.unique.tip <- function(tree, node, zero.length.tips.count=TRUE) {
  if (length(node) == 1) {
    if(verbose){
      cat("Node ",node," is unique for this patient.\n", sep="")
    }
    if(!zero.length.tips.count | !is.tip(tree,node)){
      return(node)
    } else {
      length <- get.edge.length(tree, node)
      while(length<1E-5){
        
        node <- Ancestors(tree, node, type="parent")
        if(is.root(tree, node)){
          break
        }
        if(verbose){
          cat("Length short, moving up to node ",node,"\n",sep="")
        }
        length <- get.edge.length(tree, node)
        
      }
      return(node)
    }
  } else {
    
    mrca <- mrca.phylo(tree, node)
    if(verbose){
      cat("Node ",mrca," is the MRCA of this set.\n", sep="")
    }
    
    return(mrca)
    
  }
}

#Turns the tree into a multifurcating tree; all internal nodes with zero edge lengths are removed and their children
#directly attached to their earliest ancestor with nonzero edge length



compareNA <- function(v1,v2) {
  # This function returns TRUE wherever elements are the same, including NA's,
  # and false everywhere else.
  same <- (v1 == v2)  |  (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}

get.last <- function(vec)
  return(vec[length(vec)])

split.patient <- function(tree, id, split.threshold, patients, patient.tips, patient.mrcas, patient.ids, split.counts, 
                          original.id=id, already.split=FALSE){
  
  patients <<- patients
  patient.tips <<- patient.tips
  patient.mrcas <<- patient.mrcas
  patient.ids <<- patient.ids
  split.counts <<- split.counts
  
  if(verbose){
    cat("Testing if ",id," needs to be split\n",sep="")
  }
  
  if(is.null(original.id)){
    original.id = id
  }
  
  if(!already.split){
    #this patient has not previously been split
    split.counts[[id]] <- 0
  }
  
  if(verbose){
    cat("Extracting subtree...")
  }
  
  if(length(patient.tips[[id]])<=1){
    if(verbose){
      cat(id,": single tip, no split possible\n", sep="")
    }
  } else {
    patient.subtree <- drop.tip(phy = tree, tip = tree$tip.label[-patient.tips[[id]]])
    if(verbose){
      cat(length(patient.subtree$tip.label)," tips\n",sep="")
    }
    
    # temporary multifurcation check
    #     
    #     if(length(Children(patient.subtree, getRoot(patient.subtree)))>2){
    #       stop(paste("Found a multifurcating root ",id, sep=""))
    #     }
    #     
    #can't unroot a tree with less than three edges (a bit weirdly)
    
    if(length(patient.subtree$edge.length)>=3){
      unrooted.patient.subtree <- unroot(patient.subtree)
      max.branch.length <- max(unrooted.patient.subtree$edge.length)
    } else {
      max.branch.length <- sum(patient.subtree$edge.length)
    }
    
    if(max.branch.length>split.threshold){
      if(verbose){
        cat("Splitting ",id," due to exceeding longest branch length threshold (length ",max.branch.length,")\n", sep="")
      }
      
      if(length(patient.subtree$edge.length)>=3){
        
        if(verbose){
          cat("More than two edges\n")
        }
        edge.index <- which(unrooted.patient.subtree$edge.length==max.branch.length)
        
        start.node.index <- unrooted.patient.subtree$edge[edge.index,2]
        end.node.index <- unrooted.patient.subtree$edge[edge.index,1]
        
        split.set <- list()
        split.set[[1]] <- vector()
        split.set[[2]] <- vector()
        
        for(tip in seq(1,length(unrooted.patient.subtree$tip.label))){
          current.node <- tip
          if(verbose){
            cat("Dealing with tip ",unrooted.patient.subtree$tip.label[tip],"\n",sep="" )
          }

          while(!(current.node %in% c(start.node.index, end.node.index))){
            current.node <- next.node(unrooted.patient.subtree, current.node)
          }
          if(current.node==start.node.index){
            split.set[[1]] <- c(split.set[[1]], unrooted.patient.subtree$tip.label[tip] )
          } else {
            split.set[[2]] <- c(split.set[[2]], unrooted.patient.subtree$tip.label[tip] )
          }
          
        }
      } else {
        split.set <- list()
        split.set[[1]] <- patient.subtree$tip.label[1]
        split.set[[2]] <- patient.subtree$tip.label[2]
      }
      
      patients <- patients[patients!=id]
      patient.tips[[id]] <- NULL
      patient.mrcas[[id]] <- NULL
      
      if(already.split){
        #this gets the split count right, but is a bit of a hack
        split.counts[[original.id]] <- split.counts[[original.id]]-1
      }
      
      for(tip.group in split.set){
        split.counts[[original.id]] <- split.counts[[original.id]]+1
        
        new.name <- paste(original.id,"_split",split.counts[[original.id]], sep="")
        if(verbose){
          cat("Child ",split.counts[[original.id]],": ", new.name,"\n", sep="")
        } 
        
        patients <- c(patients, new.name)
        
        tip.group.in.full.tree <- which(tree$tip.label %in% tip.group)
        
        patient.tips[[new.name]] <- tip.group.in.full.tree
        patient.mrcas[[new.name]] <- mrca.phylo.or.unique.tip(tree, tip.group.in.full.tree, zero.length.tips.count)
        
        for(tip.no in patient.tips[[new.name]]){
          patient.ids[tip.no] <- new.name
        }
        
        split.patient(tree, new.name, split.threshold, patients, patient.tips, patient.mrcas, patient.ids,
                      split.counts, original.id, already.split = TRUE)
      }
    } else if(verbose){
      cat("No split needed (length ",max.branch.length,")\n",sep="")
    }
  }
}

#This is surprisingly difficult on unrooted trees, since Ancestors still gives a root even in an unrooted tree. 
#Need a procedure to explore nodes in a subtree until it reaches one end of the broken branch. As a hack, pick 
#the next one randomly; it will get there in the end.

next.node <- function(tree, node){
  
  connections.1 <- tree$edge[which(tree$edge[,2]==node),1]
  connections.2 <- tree$edge[which(tree$edge[,1]==node),2]
  
  all.connections <- c(connections.1, connections.2)[!(is.na(c(connections.1, connections.2)))]
  
  return(all.connections[sample.int(length(all.connections), 1)]) 
}
