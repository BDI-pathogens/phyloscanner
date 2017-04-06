split.and.annotate <- function(tree, patients, patient.tips, patient.mrcas, blacklist, tip.regex, method="r", k=NA, nonancestry.penalty = 0, ties.rule = "c", useff=F, verbose=F){

  if (method == "r") {
    
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
    
    sankhoff.mode <- "total.lengths"
    
    if(is.na(k)){
      stop("k must be specified for Sankhoff reconstruction")
    }
    
    cat("Reconstructing internal node hosts with the Sankhoff algorithm...\n")
    
    non.patient.tips <- which(is.na(sapply(tree$tip.label, function(name) patient.from.label(name, tip.regex))))
    patient.ids <- sapply(tree$tip.label, function(x) patient.from.label(x, tip.regex))
    patient.ids[c(non.patient.tips, blacklist)] <- "unsampled"
    
    patients <- unique(patient.ids)
    # patients <- patients[order(patients)]
    
    patient.tips <- lapply(patients, function(x)  which(patient.ids==x))
    names(patient.tips) <- patients
    
    if(!("unsampled" %in% patients)){
      patients <- c(patients, "unsampled")
    }
    
    tip.assocs <- annotate.tips(tree, patients, patient.tips)
    
    cat("Calculating node costs...\n")
    
    # This matrix is the cost of an infection for a given patient along the branch ENDING in a given
    # node. This is the distance from the START of the branch (could make it halfway at a later point) 
    # to the first tip from that patient
    
    if(useff){
      individual.costs <- ff(Inf, dim=c(length(tree$tip.label) + tree$Nnode,length(patients))) 
    } else {
      individual.costs <- matrix(Inf, ncol=length(patients), nrow=length(tree$tip.label) + tree$Nnode) 
    }
    
    
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
      
      # Handily those which have an infinite costs as a recipient are those which have a nonancestry penalty
      # as a donor
      
      if(useff){
        finite.cost <- ff(FALSE, dim=c(length(tree$tip.label) + tree$Nnode, length(patients)))
      } else {
        finite.cost <- matrix(FALSE, ncol=length(patients), nrow=length(tree$tip.label) + tree$Nnode) 
      }
      
      
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
      
      penalties <- as.numeric(!finite.cost)*nonancestry.penalty
      penalties[,which(patients=="unsampled")] <- nonancestry.penalty
      
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
    
    cat("Building full cost matrix...\n")
    
    if(useff){
      cost.matrix <- ff(NA, dim=c(length(tree$tip.label) + tree$Nnode, length(patients)), vmode = "single")
    } else {
      cost.matrix <- matrix(NA, nrow=length(tree$tip.label) + tree$Nnode, ncol=length(patients))
    }
    

    
    cost.matrix <- make.cost.matrix(getRoot(tree), tree, patients, tip.assocs, individual.costs, penalties, cost.matrix, k)
    
    cat("Reconstructing...\n")
    
    full.assocs <- reconstruct(tree, getRoot(tree), "unsampled", list(), tip.assocs, patients, cost.matrix, individual.costs, k, ties.rule, verbose)
    
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
    
    cat("Identifying split patients...\n")
    
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

# Takes splits and annotates the tree with them

assign.splits.down <- function(node, tree, unsplit.assocs, split.assocs){
  for(child in Children(tree, node)){
    if(unsplit.assocs[[child]] == unsplit.assocs[[node]]){
      split.assocs[[child]] <- split.assocs[[node]]
      split.assocs <- assign.splits.down(child, tree, unsplit.assocs, split.assocs)
    }
  }
  return(split.assocs)
}

# Does the RS classification

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

# Make the full Sankhoff cost matrix

make.cost.matrix <- function(node, tree, patients, tip.assocs, individual.costs, penalties, current.matrix, k, verbose = F){
  if(verbose){
    cat("Node number ",node,":", sep="")
  }
  if(is.tip(tree, node)){
    if(verbose){
      cat(" is a tip (host = ",tip.assocs[[node]],")\n", sep="")
    }
    current.matrix[node,] <- rep(Inf, length(patients))
    current.matrix[node, which(patients == tip.assocs[[node]])] <- 0
    
  } else {
    if(verbose){
      cat(" looking at children (")
      cat(Children(tree, node),sep=" ")
      cat(")\n")
    }
    this.row <- vector()
    child.nos <- Children(tree, node)
    for(child in child.nos){
      current.matrix <- make.cost.matrix(child, tree, patients, tip.assocs, individual.costs, penalties, current.matrix, k, verbose)
    }
    if(length(child.nos)==0){
      stop("Reached an internal node with no children(?)")
    }
    if(verbose){
      cat("Done; back to node ",node,"\n", sep="")
    }
    penalties <- penalties[node,]
    current.matrix[node,] <- vapply(seq(1, length(patients)), function(x) node.cost(x, patients, current.matrix, individual.costs, penalties[x], k, child.nos), 0)
    
    if(verbose){
      cat("Row for node ",node," calculated\n", sep="")
    }
    
  }
  if(verbose){
    cat("Lowest cost (",min(current.matrix[node,]),") for node ",node," goes to ",patients[which(current.matrix[node,] == min(current.matrix[node,]))],"\n",sep="") 
    cat("\n")  
  }
  if(length(which(!is.na(current.matrix[,1]))) %% 100 == 0){
    cat(length(which(!is.na(current.matrix[,1]))), " of ", nrow(current.matrix), " matrix rows calculated.\n", sep="")
  }
  
  return(current.matrix)
}

node.cost <- function(patient.index, patients, current.matrix, individual.costs, penalty, k, child.nos){
  return(sum(vapply(child.nos, function (x) child.min.cost(x, patients, patient.index, current.matrix, individual.costs, penalty, k), 0)))
}

child.min.cost <- function(child.index, patients, top.patient.no, current.matrix, individual.costs, penalty, k){
  # if(impossible[Ancestors(tree, child.index, type="parent"), top.patient.no]){
  #   return(Inf)
  # }
  
  scores <- rep(Inf, length(patients))
  
  finite.scores <- c(top.patient.no, which(patients=="unsampled"), which(is.finite(individual.costs[child.index,])))
  finite.scores <- finite.scores[order(finite.scores)]
  
  scores[finite.scores] <- vapply(finite.scores, function(x) child.cost(child.index, patients, top.patient.no, x, current.matrix, individual.costs, penalty, k), 0)
  
  # scores <- vapply(seq(1, length(patients)), function(x) child.cost(child.index, patients, top.patient.no, x, current.matrix, individual.costs, k), 0)   
  
  return(min(scores))
}

# Submethods of the above

child.cost <- function(child.index, patients, top.patient.no, bottom.patient.no, current.matrix, individual.costs, penalty, k){
  if(top.patient.no == bottom.patient.no) {
    out <- current.matrix[child.index, bottom.patient.no]
  } else if (patients[bottom.patient.no] != "unsampled") {
    out <- current.matrix[child.index, bottom.patient.no] + 1 + k*individual.costs[child.index, bottom.patient.no]
    if(is.nan(out)) {
      out <- Inf
    }
    out <- out + penalty
  } else {
    # No cost or penalty for transition to unsampled
    out <- current.matrix[child.index, bottom.patient.no]
  }
  return(out)
}

# Reconstruct node states based on the cost matrix. 

reconstruct <- function(tree, node, node.state, node.assocs, tip.assocs, patients, full.cost.matrix, node.cost.matrix, k, ties.rule, verbose=F){
  node.assocs[[node]] <- node.state
  if(verbose){
    cat("Node ",node," reconstructed as ",node.state,"\n", sep="")
  }
  if(is.tip(tree, node)){
    if(node.state != tip.assocs[[node]]){
      stop("Attempting to reconstruct the wrong state onto a tip")
    }
    decision <- node.state
  } else {
    for(child in Children(tree, node)){
      costs <- vapply(seq(1, length(patients)), function(x) calc.costs(x, patients, node.state, child, node.cost.matrix, full.cost.matrix, k), 0)
      
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
        if(ties.rule == "u"){
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
        } else if(ties.rule == "c")  {
          if(node.state %in% patients[which(costs == min.cost)]){
            cat("HOI\n")
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
        } else if(ties.rule=="b"){
          # if this branch is longer than twice its parent branch, give it "unsampled". Otherwise, give it its parent's state
          
          if(node.state %in% patients[which(costs == min.cost)] & "unsampled" %in% patients[which(costs == min.cost)]){
            child.branch.length <- get.edge.length(tree, child)
            parent.branch.length <- get.edge.length(tree, node)
            if(child.branch.length > 2*parent.branch.length){
              decision <- "unsampled"
            } else {
              decision <- node.state
            }
            
          } else {
            if(verbose){
              cat(patients[which(costs == min.cost)], '\n')
            }
            paste("WARNING: some ties broken at random\n", sep="")
            decision <- patients[sample(which(costs == min.cost), 1)]
          }
          
        } else {
          stop("Unknown tie-breaking method.")
        }
      }
      node.assocs <- reconstruct(tree, child, decision, node.assocs, tip.assocs, patients, full.cost.matrix, node.cost.matrix, k, ties.rule, verbose)
    }
  }
  return(node.assocs)
}

# Submethods of the above

calc.costs <- function(patient.no, patients, node.state, child.node, node.cost.matrix, full.cost.matrix, k){
  patient <- patients[patient.no]
  if(patient == node.state){
    # No cost for the transition because it doesn't happen here
    out <- full.cost.matrix[child.node, patient.no]
  } else if(patient == "unsampled"){
    # "unsampled" - no cost for the infection
    out <- full.cost.matrix[child.node, patient.no]
  } else {
    # the cost of "patient" being at the child plus the cost of "patient" being infected along this branch 
    out <- full.cost.matrix[child.node, patient.no] + 1 + k*node.cost.matrix[child.node, patient.no]
    if(is.nan(out)){
      out <- Inf
    }
  }
  return(out)
}

cost.of.subtree <- function(tree, node, patient, tip.assocs, finite.cost.col, results){
  if(is.tip(tree, node)){
    if(patient == tip.assocs[[node]]){
      results[node] <- 0
    } else {
      stop("A tip with the wrong association has been reached")
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
