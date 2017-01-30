split.and.annotate <- function(tree, patients, patient.tips, patient.mrcas, blacklist, tip.regex, method="r", k=NA, break.ties.unsampled = FALSE, sankhoff.mode = "total.lengths", verbose=F){
  
  if(method == "c"){
    cat("Finding nodes that would have to be associated with more than one patient with no splits...\n")
    
    node.assocs <- annotate.internal(tree, patients, patient.tips, patient.mrcas)
    
    patients.with.conflicts <- vector()
    
    for(node.item in node.assocs$details){
      if(length(node.item) > 1){
        for(problem.patient in node.item){
          patients.with.conflicts <- unique(c(patients.with.conflicts, problem.patient))
        }
      }
    }
    
    cat("Resolving patients with conflicts into separate groups...\n")
    
    # Copy the patients vector and patient tip list for splitting
    
    patients.copy <- patients
    patient.tips.copy <- patient.tips
    
    pat.length <- length(patients.copy)
    
    # Split each patient
    
    for(patient in patients){
      sp.res <- split.patient(tree, patient, patients, patient.tips, node.assocs$details, tip.regex)
      split.pats <- sp.res$patients
      if(length(split.pats)>0){
        patients.copy <- patients.copy[which(patients.copy!=patient)]
        patients.copy <- c(patients.copy, split.pats)
        split.tips <- sp.res$patient.tips
        patient.tips.copy[[patient]] <- NULL
        patient.tips.copy <- c(patient.tips.copy, split.tips)
      }
    }
    
    # Need to update the MRCA list and node association list to reflect the splitting
    
    cat("Recalculating MRCAs for groups...\n")
    
    patient.mrcas.copy <- lapply(patient.tips.copy, function(node) mrca.phylo.or.unique.tip(tree, node, zero.length.tips.count))
    
    cat("Recalculating node associations for groups...\n")
    
    node.assocs <- annotate.internal(tree, patients.copy, patient.tips.copy, patient.mrcas.copy)
    
    return(list(assocs = node.assocs$details, split.patients = patients.copy, split.tips = patient.tips.copy, 
                first.nodes = patient.mrcas.copy))
    
  } else if (method == "r") {
    
    cat("Applying the Romero-Severson parsimony classification to internal nodes...\n")
    
    tip.assocs <- annotate.tips(tree, patients, patient.tips)
    
    for(tip in seq(1, length(tree$tip.label))){
      if(tip %in% blacklist | is.na(patient.from.label(tree$tip.label[tip], tip.regex))){
        tip.assocs[[tip]] <- "*"
      }
    }
    
    new.assocs <- classify.down(getRoot(tree), tree, tip.assocs, list(), patient.mrcas)
    
    cat("Filling in stars...\n")
    
    star.runs <- get.star.runs(tree, new.assocs)
    
    star.bottoms <- which(lengths(star.runs)>0)
    
    resolved.assocs <- new.assocs
    
    for(bottom in star.bottoms){
      child.assocs <- unlist(new.assocs[Children(tree, bottom)])
      c.a.df <- as.data.frame(table(child.assocs), stringsAsFactors=F)
      c.a.df <- c.a.df[which(c.a.df[,1]!="*"),]
      winners <- c.a.df[which(c.a.df[,2]==max(c.a.df[,2])),]
      possibilities <- winners[,1]
      
      last.in.chain <- star.runs[[bottom]][length(star.runs[[bottom]])]
      parent.lic <- Ancestors(tree, last.in.chain, type="parent")
      if(parent.lic!=0){
        parental.assoc <- new.assocs[[parent.lic]]
        if(parental.assoc %in% possibilities){
          
          resolved.assocs[star.runs[[bottom]]] <- parental.assoc
        } 
      }
    }
    
    # Now the splits. A new split is where you encounter a node with a different association to its parent
    
    cat("Identifying split patients...\n")
    
    splits.count <- rep(0, length(patients))
    first.nodes <- list()
    
    splits <- count.splits(tree, getRoot(tree), resolved.assocs, patients, splits.count, first.nodes)
    
    counts.by.patient <- splits$counts
    first.nodes.by.patients <- splits$first.nodes
    
    patients.copy <- patients
    patient.tips.copy <- patient.tips
    
    split.assocs <- resolved.assocs
    
    #find the splits and reannotate the tree
    
    for(pat.no in seq(1, length(patients))){
      patient <- patients[pat.no]
      no.splits <- counts.by.patient[pat.no]
      
      if(no.splits>1){
        pat.tips <- patient.tips[[patient]]
        patient.tips.copy[[patient]] <- NULL
        patients.copy <- patients.copy[which(patients.copy!=patient)]
        patients.copy <- c(patients.copy, paste(patient,"-S",seq(1, no.splits),sep=""))
        
        for(tip in pat.tips){
          # go up until you find one of the first nodes
          current.node <- tip
          subtree.roots <- first.nodes.by.patients[[patient]]
          while(!(current.node %in% subtree.roots)){
            if(current.node == 0){
              stop("Reached the root?!")
            }
            current.node <- Ancestors(tree, current.node, type="parent")
          }
          split.index <- which(subtree.roots == current.node)
          new.name <- paste(patient,"-S", split.index, sep="")
          
          current.node <- tip
          split.assocs[[subtree.roots[split.index]]] <- new.name
          
          while(current.node != subtree.roots[split.index]){
            if(current.node == 0){
              stop("Reached the root?!")
            }
            split.assocs[[current.node]] <- new.name
            current.node <- Ancestors(tree, current.node, type="parent")
          }
          
          
          patient.tips.copy[[new.name]] <- c(patient.tips.copy[[new.name]], tip)
        }
      }
    }
    
    # No need for the stars anymore
    
    return(list(assocs = split.assocs, split.patients = patients.copy, split.tips = patient.tips.copy, 
                first.nodes = first.nodes.by.patients))
    
    
  } else if (method=="s") {
    
    if(is.na(k)){
      stop("k must be specified for Sankhoff reconstruction")
    }
    
    cat("Reconstructing internal node hosts with the Sankhoff algorithm...\n")
 
    non.patient.tips <- which(is.na(sapply(tree$tip.label, function(name) patient.from.label(name, tip.regex))))
    patient.ids <- sapply(tree$tip.label, function(x) patient.from.label(x, tip.regex))
    patient.ids[c(non.patient.tips, blacklist)] <- "unsampled"
    
    patients <- unique(patient.ids)
    
    patient.tips <- lapply(patients, function(x)  which(patient.ids==x))
    names(patient.tips) <- patients
    
    if(!("unsampled" %in% patients)){
      patients <- c(patients, "unsampled")
    }
    
    tip.assocs <- annotate.tips(tree, patients, patient.tips)
    
    if(verbose){
      cat("Calculating all costs...\n")
    }
    
    # This matrix is the cost of an infection for a given patient along the branch ENDING in a given
    # node. This is the distance from the START of the branch (could make it halfway at a later point) 
    # to the first tip from that patient
    
    individual.costs <- matrix(Inf, ncol=length(patients), nrow=length(tree$tip.label) + tree$Nnode) 
    
    if(sankhoff.mode == "close.tip"){
      for(tip in seq(1, length(tree$tip.label))){
        pat <- tip.assocs[[tip]]
        individual.costs[tip, which(patients==pat)] <- 0
        current.distance <- 0
        current.node <- tip
        repeat{
          # The cost if other children of the current node are ignored...
          
          current.distance <- current.distance + get.edge.length(tree, current.node)
          
          current.node <- Ancestors(tree, current.node, type="parent")
          
          if(individual.costs[current.node, which(patients==pat)] <= current.distance){
            # already looked at a closer tip
            break
          }
          individual.costs[current.node, which(patients==pat)] <- current.distance
          if(is.root(tree, current.node)){
            break
          }
        }
      }
    } else if(sankhoff.mode=="total.lengths"){
      
      # First, go up the tree and flag which nodes have finite costs for which patients (i.e. those
      # that are ancestral to a tip _from_ that patient)
      
      finite.cost <- matrix(FALSE, ncol=length(patients), nrow=length(tree$tip.label) + tree$Nnode) 
      
      for(tip in seq(1, length(tree$tip.label))){
        pat <- tip.assocs[[tip]]
        
        finite.cost[tip, which(patients==pat)] <- TRUE
        current.node <- tip
        repeat{
          # The cost if other children of the current node are ignored...
          
          current.node <- Ancestors(tree, current.node, type="parent")
          
          if(finite.cost[current.node, which(patients==pat)]){
            # already been here
            break
          }
          finite.cost[current.node, which(patients==pat)] <- TRUE
          if(is.root(tree, current.node)){
            break
          }
        }
      }
      
      # Then traverse
      
      for(pat.no in 1:length(patients)){
        start.col <- sapply(finite.cost[,pat.no], function(x) if(x) 0 else Inf)
        
        individual.costs[,pat.no] <- cost.of.subtree(tree, getRoot(tree), patients[pat.no], tip.assocs, finite.cost[,pat.no], start.col)
      }
      
    } else {
      stop("Unsupported Sankhoff cost function")
    }
    
    # The rows of the cost matrix are nodes. The columns are patients; the last column is the unsampled
    # state
    
    cost.matrix <- matrix(NA, ncol=length(patients), nrow=length(tree$tip.label) + tree$Nnode)
    
    cost.matrix <- make.cost.matrix(getRoot(tree), tree, patients, tip.assocs, individual.costs, cost.matrix, k, verbose)
    
    if(verbose){
      cat("Reconstructing...\n")
    }
    
    full.assocs <- reconstruct(tree, getRoot(tree), "unsampled", list(), tip.assocs, patients, cost.matrix, individual.costs, k, break.ties.unsampled, verbose)
    
    temp.ca <- rep(NA, length(tree$tip.label) + tree$Nnode)
    
    for(item in seq(1, length(full.assocs))){
      if(!is.null(full.assocs[[item]])){
        if(full.assocs[[item]] != "unsampled"){
          temp.ca[item] <- full.assocs[[item]]
        }
      }
    }
    
    # Now the splits. A new split is where you encounter a node with a different association to its parent
    
    actual.patients <- patients[which(patients!="unsampled")]
    
    if(verbose){
      cat("Identifying split patients...\n")
    }
    
    splits.count <- rep(0, length(actual.patients))
    first.nodes <- list()
    
    splits <- count.splits(tree, getRoot(tree), full.assocs, actual.patients, splits.count, first.nodes)
    
    counts.by.patient <- splits$counts
    first.nodes.by.patients <- splits$first.nodes
    
    patients.copy <- actual.patients
    patient.tips.copy <- patient.tips
    
    split.assocs <- full.assocs
    
    #find the splits and reannotate the tree
    
    for(pat.no in seq(1, length(actual.patients))){
      patient <- actual.patients[pat.no]
      no.splits <- counts.by.patient[pat.no]
      
      if(no.splits>1){
        
        patient.tips.copy[[patient]] <- NULL
        patients.copy <- patients.copy[which(patients.copy!=patient)]
        patients.copy <- c(patients.copy, paste(patient,"-S",seq(1, no.splits),sep=""))
        
        subtree.roots <- first.nodes.by.patients[[patient]]
        
        for(split.no in 1:no.splits){
          split.root <- subtree.roots[split.no]
          new.name <- paste(patient,"-S", split.no, sep="")
          split.assocs[[split.root]] <- new.name
          split.assocs <- assign.splits.down(split.root, tree, full.assocs, split.assocs)
          
          for(tip in 1:length(tree$tip.label)){
            if(split.assocs[[tip]] == new.name){
              patient.tips.copy[[new.name]] <- c(patient.tips.copy[[new.name]], tip)
            }
          }
        }
      }
    }
    
    return(list(assocs = split.assocs, split.patients = patients.copy, split.tips = patient.tips.copy, 
                first.nodes = first.nodes.by.patients))
    
    
    
  } else {
    stop("Unsupported splitting method")
  }
}

# CONSERVATIVE METHODS

# Work out which nodes must be assigned to each patient for full connectivity (ignoring conflicts)

annotate.internal <- function(tree, patients, patient.tips, patient.mrcas){
  # The patients associated with each node
  node.assocs <- list()
  # The number of patients associated with each node
  node.assoc.counts <- rep(0, length(tree$tip.label) + tree$Nnode)
  for(patient in patients){
    for(tip in patient.tips[[patient]]){
      keep.going <- T
      current.node <- tip
      while(keep.going){
        # if you're at this node, it should be given this association
        
        if(is.null(node.assocs[current.node][[1]])){
          # no existing associations
          node.assocs[current.node][1] <- patient
          node.assoc.counts[current.node] <- 1
        } else {
          if(patient %in% node.assocs[current.node][[1]]){
            # already associated, no need to go further up the tree
            break
          } else {
            # add it to the vector and increase the count
            node.assocs[current.node][[1]] <- c(node.assocs[current.node][[1]], patient)
            node.assoc.counts[current.node] <- node.assoc.counts[current.node] + 1
          }
        }
        # If the patient MRCA has been reached, can stop
        if(current.node == patient.mrcas[[patient]]){
          keep.going <- F
        }
        # Otherwise, go up
        current.node <- Ancestors(tree, current.node, type="parent")
      }
    }
  }
  return(list(details = node.assocs, counts = node.assoc.counts))
}

# Find all the splits for a given patient and update the patients vector and patient tips list

split.patient <- function(tree, patient, patients, patient.tips, associations, regexp){
  #  cat("Splitting ",patient,"...\n",sep="")
  list.of.splits <- find.splits(tree, patient, patients, patient.tips, associations, regexp)
  
  new.patients <- vector()
  new.patient.tips <- list()
  
  if(length(list.of.splits) > 1){
    for(item.no in seq(1, length(list.of.splits))){
      new.name <- paste(patient,"-S",item.no, sep="")
      new.patients <- c(new.patients, new.name)
      new.patient.tips[[new.name]] <- which(tree$tip.label %in% list.of.splits[[item.no]])
    }
  }
  return(list(patients = new.patients, patient.tips = new.patient.tips))
}

find.splits <- function(tree, patient, patients, patient.tips, associations, regexp){
  #  cat("Finding splits for ",patient,"...\n",sep="")
  
  tip.pool <- tree$tip.label[patient.tips[[patient]]]
  
  out <- list()
  
  last.ancestors <- sapply(tip.pool, function(x) last.nonconflicted.ancestor(tree, x, associations, regexp))
  
  for(unique.last.anc in unique(last.ancestors)){
    out[[paste("node_", unique.last.anc,sep="")]] <- tip.pool[which(last.ancestors==unique.last.anc)]
  }
  
  return(out)
}

# Tips in the same connected subtree will have the same earliest ancestor that is associated with
# that patient but is not conflicted

last.nonconflicted.ancestor <- function(tree, tip.label, associations, regexp){
  #  cat("Finding last nonconflicted ancestor for ", tip.label, "...\n", sep="")
  tip <- which(tree$tip.label == tip.label)
  patient <- patient.from.label(tip.label, regexp)
  
  current.node <- tip
  
  repeat{
    next.node <- Ancestors(tree, current.node, type="parent")
    if(next.node == 0){
      return(current.node)
    }
    
    # associations will only contain keys up to the highest number node that has any associations
    
    if(next.node > length(associations)){
      return(current.node)
    }
    if((length(associations[[next.node]]) > 1) | !(patient %in% associations[[next.node]])){
      return(current.node)
    }
    current.node <- next.node
  }
}

# ROMERO-SEVERSON METHODS

# Annotate just the tips with their patients

annotate.tips <- function(tree, patients, patient.tips){
  node.assocs <- list()
  for(patient in patients){
    for(tip in patient.tips[[patient]]){
      node.assocs[[tip]] <- patient
    }
  }
  return(node.assocs)
}

# Used to find every connected subtree in the whole tree by finding the earliest node in them

count.splits <- function(tree, node, assocs, patients, counts.vec, first.nodes.list){
  
  for(child in Children(tree, node)){
    child.results <- count.splits(tree, child, assocs, patients, counts.vec, first.nodes.list) 
    
    counts.vec <- child.results$counts
    first.nodes.list <- child.results$first.nodes
  }
  
  parent <- Ancestors(tree, node, type="parent")
  change <- F
  if(parent == 0){
    change <- T
  } else if (assocs[[parent]]!=assocs[[node]]){
    change <- T
  }
  
  if(!(assocs[[node]] %in% c("*", "unsampled"))){
    if(change){
      patient <- assocs[[node]]
      pat.index <- which(patients == patient)
      counts.vec[pat.index] <- counts.vec[pat.index] + 1
      first.nodes.list[[patient]] <- c(first.nodes.list[[patient]], node)
    }
  }
  return(list(counts = counts.vec, first.nodes = first.nodes.list))
}

# paints new associations

assign.splits.down <- function(node, tree, unsplit.assocs, split.assocs){
  for(child in Children(tree, node)){
    if(unsplit.assocs[[child]] == unsplit.assocs[[node]]){
      split.assocs[[child]] <- split.assocs[[node]]
      split.assocs <- assign.splits.down(child, tree, unsplit.assocs, split.assocs)
    }
  }
  return(split.assocs)
}



# Does the RS classification (todo move to TUF?)

classify.down <- function(node, tree, tip.assocs, temp.assocs, patient.mrcas){
  
  current.assocs <- temp.assocs
  
  if(is.tip(tree, node)){
    result <- tip.assocs[[node]]
  } else {
    child.assocs.vector <- rep(NA, length(Children(tree, node)))
    
    for(i in seq(1,length(Children(tree, node)))){
      child <- Children(tree, node)[i]
      new.assocs <- classify.down(child, tree, tip.assocs, current.assocs, patient.mrcas)
      child.assocs.vector[i] <- new.assocs[[child]]
      current.assocs <- new.assocs
    }
    
    c.a.df <- as.data.frame(table(child.assocs.vector), stringsAsFactors=F )
    
    if(nrow(c.a.df)==1){
      result <- c.a.df[1,1]
    } else {
      c.a.df <- c.a.df[which(c.a.df[,1]!="*"),]
      
      winners <- c.a.df[which(c.a.df[,2]==max(c.a.df[,2])),]
      if(nrow(winners)==1){
        result <- winners[1,1]
      } else {
        result <- "*"
      }
    }
  }
  
  #check that this isn't above the MRCA for the patient, which leads to weird results
  
  if(result != "*"){
    if(node %in% Ancestors(tree, patient.mrcas[[result]], type="all")){
      result <- "*"
    }
  }
  
  current.assocs[[node]] <- result
  return(current.assocs)
}

# "Star runs" are sequences of ancestors that are assigned "*" by the RS classification. Once the
# classification is done, they can be converted to having assigned nodes if the ancestor of the
# earliest node is a possible host for all nodes in the run

get.star.runs <- function(tree, assocs){
  out <- list()
  for(int.node in seq(length(tree$tip.label) + 1,  length(tree$tip.label) + tree$Nnode)){
    if(assocs[[int.node]]=="*"){
      #need to check that it got the star because of an ambiguity, not because there are only stars around there
      child.assocs <- unlist(assocs[Children(tree, int.node)])
      c.a.df <- as.data.frame(table(child.assocs), stringsAsFactors=F)
      c.a.df <- c.a.df[which(c.a.df[,1]!="*"),]
      if(nrow(c.a.df)>1){
        run <- int.node
        current.node <- Ancestors(tree, int.node, type="parent")
        if(current.node != 0){
          while(assocs[[current.node]]=="*"){
            run <- c(run, current.node)
            current.node <- Ancestors(tree, current.node, type="parent")
            if(current.node == 0){
              break
            }
          }
        }
        out[[int.node]] <- run
      }
    }
  }
  return(out)
}

# SANKHOFF METHODS

# 0 is completely ignore (no important descendants); 1 is ignore when passing; 2 is don't ignore

find.ignored.descendants <- function(node, tree, temp.ignore.list, temp.surrogate.list, tip.assocs){
  current.ignore.list <- temp.ignore.list
  current.surrogate.list <- temp.surrogate.list
  
  if(is.tip(tree, node)){
    if(tip.assocs[[node]] == "*"){
      result <- 0
      current.surrogate.list[node] <- NA
    } else {
      result <- 2
      current.surrogate.list[node] <- node
    }
  } else {
    counted.children <- 0
    last.counted.child <- NA
    for(child in Children(tree, node)){
      child.results <- find.ignored.descendants(child, tree, current.ignore.list, current.surrogate.list, tip.assocs)
      
      new.ignore.list <- child.results$ignore
      new.surrogate.list <- child.results$surrogates
      if(new.ignore.list[[child]]!=0){
        counted.children <- counted.children + 1
        last.counted.child <- child
      }
      current.ignore.list <- new.ignore.list
      current.surrogate.list <- new.surrogate.list
    }
    if(counted.children == 0){
      result <- 0
      current.surrogate.list[node] <- NA
    } else if(counted.children == 1){
      result <- 1
      current.surrogate.list[node] <- last.counted.child
    } else {
      result <- 2 
      current.surrogate.list[node] <- node
    }
  }
  
  current.ignore.list[[node]] <- result
  return(list(ignore = current.ignore.list, surrogates = current.surrogate.list))
}

make.cost.matrix <- function(node, tree, patients, tip.assocs, individual.costs, current.matrix, k, verbose = T){
  if(verbose){
    cat("Node number ",node,":", sep="")
  }
  if(is.tip(tree, node)){
    if(verbose){
      cat(" is a tip (host = ",tip.assocs[[node]],")\n", sep="")
    }
    infinity.vector <- rep(Inf, length(patients))
    infinity.vector[which(patients == tip.assocs[[node]])] <- 0
    current.matrix[node,] <- infinity.vector
  } else {
    if(verbose){
      cat(" looking at children (")
      cat(Children(tree, node),sep=" ")
      cat(")\n")
    }
    this.row <- vector()
    child.nos <- Children(tree, node)
    nodes.for.calcs <- vector()
    for(child in child.nos){
      current.matrix <- make.cost.matrix(child, tree, patients, tip.assocs, individual.costs, current.matrix, k, verbose)
    }
    if(length(child.nos)==0){
      stop("Reached an internal node with no children(?)")
    }
    if(verbose){
      cat("Done; back to node ",node,"\n", sep="")
    }

    row <- vapply(seq(1, length(patients)), function(x) node.cost(x, patients, current.matrix, individual.costs, k, child.nos), 0)
    current.matrix[node,] <- row
  }
  if(verbose){
    cat("Lowest cost (",min(current.matrix[node,]),") for node ",node," goes to ",patients[which(current.matrix[node,] == min(current.matrix[node,]))],"\n",sep="") 
    cat("\n")
  }
  return(current.matrix)
  
}

node.cost <- function(patient.index, patients, current.matrix, individual.costs, k, child.nos){
  return(sum(vapply(child.nos, function (x) child.min.cost(x, patients, patient.index, current.matrix, individual.costs, k), 0)))
}

child.min.cost <- function(child.index, patients, top.patient.no, current.matrix, individual.costs, k){
  scores <- vapply(seq(1, length(patients)), function(x) child.cost(child.index, top.patient.no, x, current.matrix, individual.costs, k), 0)   
  return(min(scores))
}

child.cost <- function(child.index, top.patient.no, bottom.patient.no, current.matrix, individual.costs, k){
  bottom.patient <- patients[bottom.patient.no]
  if(top.patient.no == bottom.patient.no) {
    out <- current.matrix[child.index, bottom.patient.no]
  } else if(patients[bottom.patient.no] != "unsampled") {
    out <- current.matrix[child.index, bottom.patient.no] + 1 + k*individual.costs[child.index, bottom.patient.no]
    if(is.nan(out)) {
      out <- Inf
    }
  } else {
    out <- current.matrix[child.index, bottom.patient.no]
  }
  return(out)
}

reconstruct <- function(tree, node, node.state, node.assocs, tip.assocs, patients, full.cost.matrix, node.cost.matrix, k, break.ties.unsampled, verbose=F){
  node.assocs[[node]] <- node.state
  if(verbose){
    cat("Node ",node," reconstructed as ",node.state,"\n", sep="")
  }
  if(is.tip(tree, node)){
    if(node.state != tip.assocs[[node]]){
      stop("Something is badly wrong")
    }
    decision <- node.state
  } else {
    for(child in Children(tree, node)){
      costs <- vector()
      
      for(patient.no in seq(1, length(patients))){
        patient <- patients[patient.no]
        if(patient == node.state){
          # No cost for the transition because it doesn't happen here
          costs[patient.no] <- full.cost.matrix[child, patient.no]
        } else if(patient == "unsampled"){
          # "unsampled" - no cost for the infection
          costs[patient.no] <- full.cost.matrix[child, patient.no]
        } else {
          # the cost of "patient" being at the child plus the cost of "patient" being infected along this branch 
          costs[patient.no] <- full.cost.matrix[child, patient.no] + 1 + k*node.cost.matrix[child, patient.no]
          if(is.nan(costs[patient.no])){
            costs[patient.no] <- Inf
          }
        }
      }
      
      min.cost <- min(costs)
      if(length(which(costs == min.cost))==1){
        
        decision <- patients[which(costs == min.cost)]
        if(verbose){
          cat("Single minimum cost belongs to ", decision, "\n", sep="")
        }
      } else {
        if(verbose){
          cat("Tie at node",node,"between",patients[which(costs == min.cost)],"...", sep=" ")
        }
        if(break.ties.unsampled){
          if("unsampled" %in% patients[which(costs == min.cost)]){
            if(verbose){
              cat("broken in favour of unsampled\n")
            }
            decision <- "unsampled"
          } else if (node.state %in% patients[which(costs == min.cost)]) {
            if(verbose){
              cat("broken in favour of",node.state,"\n")
            }
            decision <- node.state
          } else {
            if(verbose){
              cat(patients[which(costs == min.cost)], '\n')
            }
            paste("WARNING: some ties broken at random\n", sep="")
            decision <- patients[sample(which(costs == min.cost), 1)]
          }
        } else {
          if(node.state %in% patients[which(costs == min.cost)]){
            if(verbose){
              cat("broken in favour of",node.state,"\n")
            }
            decision <- node.state
          } else if("unsampled" %in% patients[which(costs == min.cost)]) {
            if(verbose){
              cat("broken in favour of unsampled\n")
            }
            decision <- "unsampled"
          } else {
            if(verbose){
              cat(patients[which(costs == min.cost)], '\n')
            }
            paste("WARNING: some ties broken at random\n", sep="")
            decision <- patients[sample(which(costs == min.cost), 1)]
          }
        }
      }
      node.assocs <- reconstruct(tree, child, decision, node.assocs, tip.assocs, patients, full.cost.matrix, node.cost.matrix, k, break.ties.unsampled)
    }
  }
  return(node.assocs)
}

cost.of.subtree <- function(tree, node, patient, tip.assocs, finite.cost.col, results){
  if(is.tip(tree, node)){
    if(patient == tip.assocs[[node]]){
      results[node] <- 0
    } else {
      stop("This should not happen")
    }
  } else {
    for(child in Children(tree, node)){
      if(finite.cost.col[child]){
        results <- cost.of.subtree(tree, child, patient, tip.assocs, finite.cost.col, results)
        results[node] <- results[node] + results[child] + get.edge.length(tree, child)
      }
    }
  }
  return(results)
}

# collapsed tree methods

output.trans.tree <- function(tree, assocs, file.name = NULL){
  # find the association of each node
  
  assocs.vec <- vector()
  
  for(node.no in seq(1, tree$Nnode + length(tree$tip.label))){
    if(node.no > length(assocs)){
      # annoying and hacky
      assocs.vec[node.no] <- "none"
    } else if (is.null(assocs[[node.no]])){
      assocs.vec[node.no] <- "none"
    } else if (is.na(assocs[[node.no]])){
      assocs.vec[node.no] <- "none"
    } else if (assocs[[node.no]] %in% c("*", "unsampled")) {
      assocs.vec[node.no] <- "none"
    } else {
      if(length(assocs[[node.no]])>1){
        stop("This method is intended for use on trees with conflicts resolved")
      }
      assocs.vec[node.no] <- assocs[[node.no]]
    }
  }
  
  first.of.split <- vector()
  
  for(node.no in seq(1, tree$Nnode + length(tree$tip.label))){
    parent.no <- Ancestors(tree, node.no, type="parent")
    if(parent.no == 0){
      first.of.split[node.no] <- T
    } else {
      first.of.split[node.no] <- assocs.vec[parent.no] != assocs.vec[node.no]
    }
  }
  
  splits.vec <- assocs.vec
  
  splits.vec[which(splits.vec=="none" & first.of.split)] <- 
    paste("none", which(splits.vec=="none" & first.of.split), sep="-S")
  
  for(node.no in seq(1, tree$Nnode + length(tree$tip.label))){
    if(assocs.vec[node.no]=="none" & !first.of.split[node.no]){
      current.node.no <- node.no
      while(!first.of.split[current.node.no]){
        current.node.no <- Ancestors(tree, current.node.no, type="parent")
      }
      if(!startsWith(splits.vec[current.node.no], "none-")){
        print(node.no)
        stop("Uh-oh")
      }
      splits.vec[node.no] <- splits.vec[current.node.no]
    }
  }
  
  unique.splits <- unique(splits.vec)
  parent.splits <- vector()
  lengths <- vector()
  root.nos <- vector()
  
  for(unique.split in unique.splits){
    root <- which(splits.vec == unique.split & first.of.split)
    parent.node <- Ancestors(tree, root, type="parent")
    if(parent.node != 0){
      parent.splits <- c(parent.splits, splits.vec[parent.node])
      
      root.edge.no <- which(tree$edge[,2]==root)
      lengths <- c(lengths, tree$edge.length[root.edge.no])
    } else {
      parent.splits <- c(parent.splits, "root")
      lengths <- c(lengths, 0)
    }
    root.nos <- c(root.nos, root)
  }
  
  patients <-  unlist(lapply(strsplit(unique.splits, "-"), `[[`, 1)) 
  parent.patients <-  unlist(lapply(strsplit(parent.splits, "-"), `[[`, 1)) 
  
  cytoscape.input <- data.frame(unique.splits, parent.splits, patients, parent.patients, lengths, root.nos, stringsAsFactors = F)
  
  if(!is.null(file.name)){
    write.csv(cytoscape.input[,1:5], file.name, row.names = F, quote=F)
  }
  
  return(cytoscape.input)
}

check.contiguous <- function(tt, patients, splits.for.patients, patients.for.splits){
  if(length(patients)!=2){
    stop("Not implemented")
  }
  pat.1.id <- patients[1]
  pat.2.id <- patients[2]
  
  OK <- TRUE
  all.nodes <-  c(splits.for.patients[[pat.1.id]], splits.for.patients[[pat.2.id]])
  
  for(node.1 in seq(1, length(all.nodes))){
    for(node.2 in seq(1, length(all.nodes))){
      if(node.1 < node.2){
        node.1.id <- all.nodes[node.1]
        node.2.id <- all.nodes[node.2]
        path <- get.tt.path(tt, node.1.id, node.2.id)
        for(node in path){
          if(!startsWith(node, "none")){
            if(!(patients.for.splits[[node]] %in% c(pat.1.id, pat.2.id))){
              OK <- FALSE
              break
            }
          }
        }
      }
      if(!OK){
        break
      }
    }
    if(!OK){
      break
    }
  }
  
  return(OK)
}

extract.tt.subtree <- function(tt, patients, splits.for.patients, patients.for.splits){
  # for now, at least
  
  if(length(patients)!=2){
    stop("Not implemented")
  }
  if(!check.contiguous(tt, patients, splits.for.patients, patients.for.splits)){
    stop("Not contiguous")
  }
  
  pat.1.id <- patients[1]
  pat.2.id <- patients[2]
  
  pat.1.splts <- splits.for.patients[[pat.1.id]]
  pat.2.splts <- splits.for.patients[[pat.2.id]]
  
  sub.tt <- tt[which(tt$unique.splits %in% c(pat.1.splts, pat.2.splts)),]
  unsampled.below <- tt[which(tt$patients == "none" & (tt$parent.splits %in% c(pat.1.splts, pat.2.splts))),]
  unsampled.above <- tt[which(tt$patients == "none" & (tt$unique.splits %in% sub.tt$parent.splits)),]
  
  none.but.maybe.relevant <- c(unsampled.above$unique.splits, unsampled.below$unique.splits)
  
  adjacent.relevance.count <- sapply(none.but.maybe.relevant, function(x) length(intersect(c(pat.1.splts, pat.2.splts), get.tt.adjacent(tt, x) ) ))
  
  return(c(sub.tt$unique.splits, unique(none.but.maybe.relevant[which(adjacent.relevance.count > 1)])))
}

# every distance between subtrees

all.subtree.distances <- function(tree, tt, splits, assocs){
  tree.dist <- dist.nodes(tree)
  
  temp <- matrix(ncol = length(splits), nrow=length(splits))
  for(spt.1.no in 1:length(splits)){
    for(spt.2.no in 1:length(splits)){
      if(spt.1.no==spt.2.no){
        temp[spt.1.no, spt.2.no] <- 0
      } else if(spt.1.no<spt.2.no){
        spt.1 <- splits[spt.1.no]
        spt.2 <- splits[spt.2.no]
        
        chain.1 <- get.tt.ancestors(tt, spt.1)
        chain.2 <- get.tt.ancestors(tt, spt.2)
        
        if(spt.1 %in% chain.2){
          mrca.2 <- tt$root.nos[which(tt$unique.splits==spt.2)]
          current.node <- mrca.2
          length <- 0
          while(assocs[[current.node]]!=spt.1){
            length <- length + get.edge.length(tree, current.node)
            current.node <- Ancestors(tree, current.node, type="parent")
            if(is.root(tree, current.node)){
              stop("Big problem here (1)")
            }
          }
          temp[spt.1.no, spt.2.no] <- length
          
        } else if(spt.2 %in% chain.1){
          mrca.1 <- tt$root.nos[which(tt$unique.splits==spt.1)]
          current.node <- mrca.1
          length <- 0
          while(assocs[[current.node]]!=spt.2){
            length <- length + get.edge.length(tree, current.node)
            current.node <- Ancestors(tree, current.node, type="parent")
            if(is.root(tree, current.node)){
              stop("Big problem here (1)")
            }
          }
          temp[spt.1.no, spt.2.no] <- length
          
        } else {
          
          mrca.1 <- tt$root.nos[which(tt$unique.splits==spt.1)]
          mrca.2 <- tt$root.nos[which(tt$unique.splits==spt.2)]
          temp[spt.1.no, spt.2.no] <- tree.dist[mrca.1, mrca.2]
        }
        temp[spt.2.no, spt.1.no] <- temp[spt.1.no, spt.2.no]
      }
    }
  }
  
  
  colnames(temp) <- splits
  rownames(temp) <- splits
  return(temp)
}

# are each pair of subtrees adjacent? If unsampled does not matter then two nodes separated only by "none"
# are still adjacent

subtrees.adjacent <- function(tt, splits, none.matters = F){
  out <- matrix(ncol = length(splits), nrow=length(splits))
  for(spt.1.no in 1:length(splits)){
    for(spt.2.no in 1:length(splits)){
      if(spt.1.no==spt.2.no){
        out[spt.1.no, spt.2.no] <- NA
      } else if(spt.1.no<spt.2.no){
        spt.1 <- splits[spt.1.no]
        spt.2 <- splits[spt.2.no]
        
        if(spt.1 %in%  get.tt.adjacent(tt, spt.2)){
          out[spt.1.no, spt.2.no] <- T
          out[spt.2.no, spt.1.no] <- T
        } else if(!none.matters) {
          path <- get.tt.path(tt, spt.1, spt.2)
          internal.path <- path[2:(length(path)-1)]
          adj <- length(internal.path)==1 & startsWith(internal.path[1], "none")
          out[spt.1.no, spt.2.no] <- adj
          out[spt.2.no, spt.1.no] <- adj
        } else {
          out[spt.1.no, spt.2.no] <- F
          out[spt.2.no, spt.1.no] <- F
        }
      }
    }
  }
  colnames(out) <- splits
  rownames(out) <- splits
  
  return(out)
}

# are pairs of subtrees from two patients not separated by any other subtrees from either of those patients?

subtrees.unblocked <- function(tt, splits){
  out <- matrix(ncol = length(splits), nrow=length(splits))
  for(spt.1.no in 1:length(splits)){
    for(spt.2.no in 1:length(splits)){
      if(spt.1.no==spt.2.no){
        out[spt.1.no, spt.2.no] <- NA
      } else if(spt.1.no<spt.2.no){
        spt.1 <- splits[spt.1.no]
        spt.2 <- splits[spt.2.no]
        pat.1 <- strsplit(spt.1, "-S")[[1]][1]
        pat.2 <- strsplit(spt.2, "-S")[[1]][1]
        
        if(spt.1 %in%  get.tt.adjacent(tt, spt.2)){
          out[spt.1.no, spt.2.no] <- T
          out[spt.2.no, spt.1.no] <- T
        } else {
          path <- get.tt.path(tt, spt.1, spt.2)
          internal.path <- path[2:(length(path)-1)]
          blockers <- which(startsWith(internal.path, pat.1) | startsWith(internal.path, pat.2))
          
          adj <- length(blockers) == 0
          out[spt.1.no, spt.2.no] <- adj
          out[spt.2.no, spt.1.no] <- adj
        } 
      }
    }
  }
  colnames(out) <- splits
  rownames(out) <- splits
  
  return(out)
}