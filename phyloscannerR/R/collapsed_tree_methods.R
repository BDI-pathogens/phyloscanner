# Collapsed tree navigation

# The parent of a node in the collapsed tree

#' @keywords internal
#' @export get.tt.parent

get.tt.parent <- function(tt, label){
  if(length(which(tt$unique.splits==label))>0){
    return(tt$parent.splits[which(tt$unique.splits==label)])
  } else {
    return(NULL)
  }
}

# All ancestors of a node in the collapsed tree.

#' @keywords internal
#' @export get.tt.ancestors

get.tt.ancestors <- function(tt, label){
  out <- vector()
  current.node <- get.tt.parent(tt, label)
  while(current.node != "root"){
    out <- c(out, current.node)
    current.node <- get.tt.parent(tt, current.node)
  }
  return(out)
}

# All children of a node in the collapsed tree

#' @keywords internal
#' @export get.tt.children

get.tt.children <- function(tt, label){
  return(tt$unique.splits[which(tt$parent.splits==label)])
}

# All adjacent nodes of a node in the collapsed tree (parent and children)

#' @keywords internal
#' @export get.tt.adjacent

get.tt.adjacent <- function(tt, label){
  return(c(get.tt.children(tt, label), get.tt.parent(tt, label)))
}

# The grandparent of a node in the collapsed tree

#' @keywords internal
#' @export get.tt.parent.host

get.tt.parent.host <- function(tt, label){
  if(length(which(tt$unique.splits==label))>0){
    return(tt$hosts[which(tt$unique.splits==label)])
  } else {
    return(NULL)
  }
}

# MRCA of a pair of nodes in the collapsed tree

#' @keywords internal
#' @export get.tt.mrca

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

# The host corresponding to the MRCA of two nodes

#' @keywords internal
#' @export get.tt.mrca.host

get.tt.mrca.host <- function(tt, label1, label2){
  node <- intersect(c(label1,get.tt.ancestors(tt, label1)), c(label2, get.tt.ancestors(tt, label2)))[1]
  
  return(tt$hosts[which(tt$unique.splits==node)])
}

# The path from one node in the collapsed tree to another

#' @keywords internal
#' @export get.tt.path

get.tt.path <- function(tt, label1, label2){
  mrca <- get.tt.mrca(tt, label1, label2)
  
  if(mrca == label1){
    result <- vector()
    current.node <- label2
    while(current.node != label1){
      result <- c(result, current.node)
      current.node <- get.tt.parent(tt, current.node)
    }
    result <- c(result, label1)
    return(result)
  }
  
  if(mrca == label2){
    result <- vector()
    current.node <- label1
    while(current.node != label2){
      result <- c(result, current.node)
      current.node <- get.tt.parent(tt, current.node)
    }
    result <- c(result, label2)
    return(result)
  }
  
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

# Output the collapsed tree from a phylogeny and its node associations

#' @keywords internal
#' @export output.trans.tree

output.trans.tree <- function(tree, assocs){
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
    } else if (assocs[[node.no]] %in% c("*", "unassigned")) {
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
  
  unassigned.roots <- which(splits.vec=="none" & first.of.split)
  
  splits.vec[unassigned.roots] <- 
    paste("unassigned_region-SPLIT", 1:length(unassigned.roots), sep="")
  
  for(node.no in seq(1, tree$Nnode + length(tree$tip.label))){
    if(assocs.vec[node.no]=="none" & !first.of.split[node.no]){
      current.node.no <- node.no
      while(!first.of.split[current.node.no]){
        current.node.no <- Ancestors(tree, current.node.no, type="parent")
      }
      if(!grepl("^unassigned_region", splits.vec[current.node.no])){      
        stop("Parent of unassigned node is not one of the unassigned subgraph roots")
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
  
  hosts <-  unlist(lapply(strsplit(unique.splits, "-SPLIT"), `[[`, 1)) 
  parent.hosts <-  unlist(lapply(strsplit(parent.splits, "-SPLIT"), `[[`, 1)) 
  
  tt.table <- data.frame(unique.splits, parent.splits, hosts, parent.hosts, lengths, root.nos, stringsAsFactors = F)
  
  return(tt.table)
}

#' @keywords internal
#' @export prune.unassigned.tips

prune.unassigned.tips <- function(tt.table){
  
  for.output <- tt.table[,1:6]
  
  unassigned.tips <- which(grepl("unassigned",for.output$unique.splits) &
                            !(for.output$unique.splits %in% for.output$parent.splits))
  
  if(length(unassigned.tips) > 0){
    for.output <- for.output[-unassigned.tips,]
  }
  #renumber
  unassigned.rows <- which(grepl("^unassigned_region",for.output$unique.splits))
  
  unassigned.labels <- for.output$unique.splits[which(grepl("^unassigned_region",for.output$unique.splits))]
  
  for(x in 1:length(unassigned.rows)) {
    old.label <- unassigned.labels[x]
    new.label <- paste("UnassignedRegion-SPLIT",x,sep="")
    for.output$unique.splits[unassigned.rows[x]] <- new.label
    for.output$parent.splits[which(for.output$parent.splits==old.label)] <- new.label
  } 
  
  
  for.output$hosts[which(grepl("^unassigned_region",for.output$hosts))] <- "UnassignedRegion"
  
  for.output$parent.hosts[which(grepl("^unassigned_region",for.output$parent.hosts))] <- "UnassignedRegion"
  
  return(for.output)
}

#' @keywords internal
#' @export check.contiguous

check.contiguous <- function(tt, hosts, splits.for.hosts, hosts.for.splits){
  if(length(hosts)!=2){
    stop("Not implemented")
  }
  pat.1.id <- hosts[1]
  pat.2.id <- hosts[2]
  
  OK <- TRUE
  all.nodes <-  c(splits.for.hosts[[pat.1.id]], splits.for.hosts[[pat.2.id]])
  
  for(node.1 in seq(1, length(all.nodes))){
    for(node.2 in seq(1, length(all.nodes))){
      if(node.1 < node.2){
        node.1.id <- all.nodes[node.1]
        node.2.id <- all.nodes[node.2]
        path <- get.tt.path(tt, node.1.id, node.2.id)
        for(node in path){
          if(!grepl("^unassigned_region",node)){
            if(!(hosts.for.splits[[node]] %in% c(pat.1.id, pat.2.id))){
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

#' @keywords internal
#' @export check.uninterrupted

check.uninterrupted <- function(tt, hosts, splits.for.hosts, hosts.for.splits){
  # this could certainly be faster
  if(length(hosts)!=2){
    stop("Not implemented")
  }
  pat.1.id <- hosts[1]
  pat.2.id <- hosts[2]
  
  nodes.1 <- splits.for.hosts[[pat.1.id]]
  nodes.2 <- splits.for.hosts[[pat.2.id]]
  
  any.contiguity <- F
  any.interruption <- F
  
  for(node.1 in nodes.1){
    current.node <- get.tt.parent(tt, node.1)
    while(current.node!="root" & (grepl("^unassigned_region",current.node) | (if(is.null(hosts.for.splits[[current.node]])) {T} else {hosts.for.splits[[current.node]]==pat.1.id}))){
      
      current.node <- get.tt.parent(tt, current.node)
    }
    if(current.node != "root"){
      if(hosts.for.splits[[current.node]]==pat.2.id){
        any.contiguity <- T
      } else {
        # this is a blocking node
        chain <- get.tt.ancestors(tt, current.node)
        if(length(intersect(c(pat.1.id, pat.2.id), unlist(hosts.for.splits[chain])))>0){
          any.interruption <- T
          break
        }
      }
    }
  }
  
  if(!any.interruption){
    for(node.2 in nodes.2){
      current.node <- get.tt.parent(tt, node.2)
      while(current.node!="root" & (grepl("unassigned_region",current.node) | (if(is.null(hosts.for.splits[[current.node]])) {T} else {hosts.for.splits[[current.node]]==pat.2.id}))){
        current.node <- get.tt.parent(tt, current.node)
      }
      if(current.node != "root"){
        if(hosts.for.splits[[current.node]]==pat.1.id){
          any.contiguity <- T
        } else {
          # this is a blocking node
          chain <- get.tt.ancestors(tt, current.node)
          if(length(intersect(c(pat.1.id, pat.2.id), unlist(hosts.for.splits[chain])))>0){
            any.interruption <- T
            break
          }
        }
      }
      
    }
  }
  return(any.contiguity & !any.interruption)
}

#' @keywords internal
#' @export extract.tt.subgraph

extract.tt.subgraph <- function(tt, hosts, splits.for.hosts, hosts.for.splits){
  # for now, at least
  
  if(length(hosts)!=2){
    stop("Not implemented")
  }
  if(!check.contiguous(tt, hosts, splits.for.hosts, hosts.for.splits)){
    stop("Not contiguous")
  }
  
  pat.1.id <- hosts[1]
  pat.2.id <- hosts[2]
  
  pat.1.splts <- splits.for.hosts[[pat.1.id]]
  pat.2.splts <- splits.for.hosts[[pat.2.id]]
  
  sub.tt <- tt[which(tt$unique.splits %in% c(pat.1.splts, pat.2.splts)),]
  unassigned.below <- tt[which(tt$hosts == "unassigned_region" & (tt$parent.splits %in% c(pat.1.splts, pat.2.splts))),]
  unassigned.above <- tt[which(tt$hosts == "unassigned_region" & (tt$unique.splits %in% sub.tt$parent.splits)),]
  
  none.but.maybe.relevant <- c(unassigned.above$unique.splits, unassigned.below$unique.splits)
  
  adjacent.relevance.count <- sapply(none.but.maybe.relevant, function(x) length(intersect(c(pat.1.splts, pat.2.splts), get.tt.adjacent(tt, x) ) ))
  
  return(c(sub.tt$unique.splits, unique(none.but.maybe.relevant[which(adjacent.relevance.count > 1)])))
}

# every distance between subgraphs

#' @keywords internal
#' @export all.subgraph.distances
#' @importFrom ape dist.nodes node.depth.edgelength

all.subgraph.distances <- function(tree, tt, splits, assocs, slow=F, total.pairs, verbose = F, no.progress.bars = F){
  
  if (verbose & !no.progress.bars) progress.bar <- txtProgressBar(width=50, style=3)
  
  if(!slow){
    tree.dist <- dist.nodes(tree)
  } else {
    depths <- node.depth.edgelength(tree)
  }
  
  count <- 0
  
  temp <- matrix(ncol = length(splits), nrow=length(splits))
  
  if(total.pairs==0 & verbose & !no.progress.bars){
    setTxtProgressBar(progress.bar, 1)
  } else {
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
                stop("Reached the root?")
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
                stop("Reached the root?")
              }
            }
            temp[spt.1.no, spt.2.no] <- length
            
          } else {
            
            mrca.1 <- tt$root.nos[which(tt$unique.splits==spt.1)]
            mrca.2 <- tt$root.nos[which(tt$unique.splits==spt.2)]
            if(!slow){
              temp[spt.1.no, spt.2.no] <- tree.dist[mrca.1, mrca.2]
            } else {
              temp[spt.1.no, spt.2.no] <- pat.dist(tree, depths, mrca.1, mrca.2)
            }
          }
          
          count <- count + 1
          
          if(verbose & !no.progress.bars){
            setTxtProgressBar(progress.bar, count/total.pairs)
          }
          temp[spt.2.no, spt.1.no] <- temp[spt.1.no, spt.2.no]
        }
      }
    }
  }
  
  if(verbose & !no.progress.bars) close(progress.bar)
  colnames(temp) <- splits
  rownames(temp) <- splits
  return(temp)
}

#' @keywords internal
#' @export pat.dist

pat.dist <- function(tree, depths, node.1, node.2){
  node.1.chain <- Ancestors(tree, node.1, type="all")
  node.2.chain <- Ancestors(tree, node.2, type="all")
  if(node.1 %in% node.2.chain){
    mrca <- node.1
    dist <- depths[node.2] - depths[node.1]
  } else if(node.2 %in% node.1.chain){
    mrca <- node.2
    dist <- depths[node.1] - depths[node.2]
  } else {
    mrca <- node.1.chain[which(node.1.chain %in% node.2.chain)][1]
    dist <- (depths[node.1] - depths[mrca]) + (depths[node.2] - depths[mrca])
  }
  return(dist)
}

# are each pair of subgraphs adjacent? If !none.matters then two nodes separated only by "none" are still adjacent

#' @keywords internal
#' @export subgraphs.adjacent

subgraphs.adjacent <- function(tt, splits, none.matters = F){
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
          adj <- length(internal.path)==1 & grepl("unassigned_region",internal.path[1])
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

# are pairs of subgraphs from two hosts not separated by any other subgraphs from either of those hosts?

#' @keywords internal
#' @export subgraphs.unblocked

subgraphs.unblocked <- function(tt, splits, total.pairs, verbose = F, no.progress.bars = F){

  if (verbose & !no.progress.bars) progress.bar <- txtProgressBar(width=50, style=3)
  count <- 0
  
  out <- matrix(ncol = length(splits), nrow=length(splits))
  if(total.pairs==0 & verbose & !no.progress.bars){
    setTxtProgressBar(progress.bar, 1)
  } else {
    
    for(spt.1.no in 1:length(splits)){
      for(spt.2.no in 1:length(splits)){
        if(spt.1.no==spt.2.no){
          out[spt.1.no, spt.2.no] <- NA
        } else if(spt.1.no<spt.2.no){
          spt.1 <- splits[spt.1.no]
          spt.2 <- splits[spt.2.no]
          
          pat.1 <- strsplit(spt.1, "-SPLIT")[[1]][1]
          pat.2 <- strsplit(spt.2, "-SPLIT")[[1]][1]
          
          if(spt.1 %in%  get.tt.adjacent(tt, spt.2)){
            out[spt.1.no, spt.2.no] <- T
            out[spt.2.no, spt.1.no] <- T
          } else {
            path <- get.tt.path(tt, spt.1, spt.2)
            internal.path <- path[2:(length(path)-1)]
            blockers <- which(grepl(paste0('^',pat.1),internal.path) | grepl(paste0('^',pat.2),internal.path))
            
            adj <- length(blockers) == 0
            out[spt.1.no, spt.2.no] <- adj
            out[spt.2.no, spt.1.no] <- adj
          } 
          count <- count + 1
          if (verbose & !no.progress.bars) {
            setTxtProgressBar(progress.bar, count/total.pairs)
          }
        }
      }
    }
  }
  if (verbose & !no.progress.bars) close(progress.bar)
  
  colnames(out) <- splits
  rownames(out) <- splits
  
  return(out)
}


#' @keywords internal
#' @export check.adjacency

check.adjacency <- function(tt, hosts, splits.for.hosts){
  # this could certainly be faster
  if(length(hosts)!=2){
    stop("Not implemented")
  }
  pat.1.id <- hosts[1]
  pat.2.id <- hosts[2]
  
  nodes.1 <- splits.for.hosts[[pat.1.id]]
  nodes.2 <- splits.for.hosts[[pat.2.id]]
  
  for(node.1 in nodes.1){
    for(node.2 in nodes.2){
      if(check.tt.node.adjacency(tt, node.1, node.2, T)){
        return(T)
      }
    }
  } 
  
  return(F)
}

#' @keywords internal
#' @export check.tt.node.adjacency

check.tt.node.adjacency <- function(tt, label1, label2, allow.unassigned = F){
  path <- get.tt.path(tt, label1, label2)
  
  if(!allow.unassigned){
    return(length(path)==2)
  }
  
  if(length(path)==2){
    # they must be next to each other
    return(T)
  }
  
  if(length(path)>3){
    # they can't be - a path length greater than 2 can only go through an unassigned region, and adjacent unassigned regions are not allowed
    return(F)
  }

  # #START TEMPORARY BIT - if the middle node is the region around the root this adjacency is not interesting
  # 
  # if(substr(path[2], 1, 17) == "unassigned_region" & get.tt.parent(tt, path[2])=="root"){
  #   return(F)
  # }
  # 
  # #END TEMPORARY BIT
  
  return(grepl("^unassigned_region",path[2]))
  
}

#' @keywords internal
#' @export classify

classify <- function(tree.info, verbose = F, no.progress.bars = F) {	
  
  if(is.null(tree.info[["tree"]])){
    
    if (verbose) cat("Reading tree file ", tree.info$tree.file.name, "...\n", sep = "")
    
    pseudo.beast.import <- read.beast(tree.info$tree.file.name)
    tree <- attr(pseudo.beast.import, "phylo")
    
    if (verbose) cat("Reading annotations...\n")
    
    annotations <- attr(pseudo.beast.import, "stats")
    
    annotations$INDIVIDUAL <- as.character(annotations$INDIVIDUAL)
    annotations$SPLIT <- as.character(annotations$SPLIT)
    
    # to deal with a problem with older versions of ggtree
    
    strip.quotes <- function(string) if(substr(string, 1, 1)=="\"") substr(string, 2, nchar(string)-1) else string
    
    annotations$INDIVIDUAL <- sapply(annotations$INDIVIDUAL, strip.quotes)
    annotations$SPLIT <- sapply(annotations$SPLIT, strip.quotes)
    
  } else {
    
    tree <- tree.info$tree
    
    annotations <- data.frame(node = seq(1, length(tree$tip.label) + tree$Nnode), 
                              INDIVIDUAL = as.character(attr(tree, "INDIVIDUAL")), 
                              SPLIT = as.character(attr(tree, "SPLIT")), 
                              stringsAsFactors = F)
  }
  
  
  if(is.null(tree.info$splits.table)){
    if (verbose) cat("Reading splits file", tree.info$splits.file.name, "...\n")
    
    splits <- read.csv(tree.info$splits.file.name, stringsAsFactors = F)
  } else {
    splits <- tree.info$splits.table
  }
  
  if (verbose) cat("Collecting tips for each host...\n")
  
  hosts <- unique(splits$host)
  
  hosts <- hosts[hosts!="unassigned"]
  
  all.splits <- unique(splits$subgraph)
  all.splits <- all.splits[all.splits!="unassigned"]
  
  in.order <- match(seq(1, length(tree$tip.label) + tree$Nnode), annotations$node)
  
  assocs <- annotations$SPLIT[in.order]
  assocs <- lapply(assocs, function(x) replace(x, is.na(x), "none"))
  assocs <- lapply(assocs, function(x) replace(x, x=="unassigned", "none"))
  
  splits.for.hosts <- lapply(hosts, function(x) unique(splits$subgraph[which(splits$host==x)] ))
  names(splits.for.hosts) <- hosts
  
  hosts.for.splits <- lapply(all.splits, function(x) unique(splits$host[which(splits$subgraph==x)] ))
  names(hosts.for.splits) <- all.splits
  
  hosts.included <- hosts
  
  total.split.pairs <- (length(all.splits) ^ 2 - length(all.splits))/2
  total.host.pairs <- (length(hosts)^2 - length(hosts))/2
  
  if (verbose) cat("Collapsing subgraphs...\n")
  
  tt <- output.trans.tree(tree, assocs)
  
  if (verbose) cat("Identifying pairs of unblocked splits...\n")
  
  collapsed.adjacent <- subgraphs.unblocked(tt, all.splits, total.split.pairs, verbose, no.progress.bars)
  
  if (verbose) cat("Calculating pairwise distances between splits...\n")
  
  split.distances <- tryCatch(
    all.subgraph.distances(tree, tt, all.splits, assocs, FALSE, total.split.pairs, verbose, no.progress.bars), warning=function(w){return(NULL)}, error=function(e){return(NULL)})
  
  if(is.null(split.distances)){
    split.distances <- all.subgraph.distances(tree, tt, all.splits, assocs, TRUE, total.split.pairs, verbose, no.progress.bars)
  }
  
  if (verbose) cat("Testing pairs...\n")
  
  progress.bar <- NULL
  
  if (verbose & !no.progress.bars) progress.bar <- txtProgressBar(width=50, style=3)
  
  count <- 0
  adjacency.matrix <- matrix(NA, length(hosts.included), length(hosts.included))
  contiguity.matrix <- matrix(NA, length(hosts.included), length(hosts.included))
  path.matrix <- matrix(NA, length(hosts.included), length(hosts.included))
  nodes.1.matrix <- matrix(NA, length(hosts.included), length(hosts.included))
  nodes.2.matrix <- matrix(NA, length(hosts.included), length(hosts.included))
  dir.12.matrix <- matrix(NA, length(hosts.included), length(hosts.included))
  dir.21.matrix <- matrix(NA, length(hosts.included), length(hosts.included))
  min.distance.matrix <- matrix(NA, length(hosts.included), length(hosts.included))
  
  if(total.host.pairs==0 & verbose & !no.progress.bars){
    setTxtProgressBar(progress.bar, 1)
  } else {
    
    for(pat.1 in seq(1, length(hosts.included))){
      for(pat.2 in  seq(1, length(hosts.included))){
        if (pat.1 < pat.2) {
          
          count <- count + 1
          
          pat.1.id <- hosts.included[pat.1]
          pat.2.id <- hosts.included[pat.2]
          
          nodes.1 <- splits.for.hosts[[pat.1.id]]
          nodes.2 <- splits.for.hosts[[pat.2.id]]
          
          all.nodes <-  c(nodes.1, nodes.2)
          
          adjacency.matrix[pat.1, pat.2] <- check.adjacency(tt, c(pat.1.id, pat.2.id), splits.for.hosts)		
          contiguity.matrix[pat.1, pat.2] <- check.contiguous(tt, c(pat.1.id, pat.2.id), splits.for.hosts, hosts.for.splits)		
          
          
          count.12 <- 0
          count.21 <- 0
          
          for(node.2 in nodes.2){
            ancestors <- get.tt.ancestors(tt, node.2)
            if(length(intersect(ancestors, nodes.1)) > 0){
              count.12 <- count.12 + 1
            }
          }
          
          for(node.1 in nodes.1){
            ancestors <- get.tt.ancestors(tt, node.1)
            if(length(intersect(ancestors, nodes.2)) > 0){
              count.21 <- count.21 + 1
            }
          }
          
          prop.12 <- count.12/length(nodes.2)
          prop.21 <- count.21/length(nodes.1)
          
          nodes.1.matrix[pat.1, pat.2] <- length(nodes.1)
          nodes.2.matrix[pat.1, pat.2] <- length(nodes.2)
          
          dir.12.matrix[pat.1, pat.2] <- count.12
          dir.21.matrix[pat.1, pat.2] <- count.21
          
          if(count.12 == 0 & count.21 == 0){
            path.matrix[pat.1, pat.2] <- "none"
          } else if(count.12 != 0 & count.21 == 0 & prop.12 == 1) {
            if(count.12 == 1){
              path.matrix[pat.1, pat.2] <- "anc"
            } else {
              path.matrix[pat.1, pat.2] <- "multiAnc"
            }
          } else if(count.21 != 0 & count.12 == 0 & prop.21 == 1) {
            if(count.21 == 1){
              path.matrix[pat.1, pat.2] <- "desc"
            } else {
              path.matrix[pat.1, pat.2] <- "multiDesc"
            }
          } else {
            path.matrix[pat.1, pat.2] <- "complex"
          }
          
          pairwise.distances <- vector()
          
          for(node.1 in nodes.1){
            for(node.2 in nodes.2){
              if(collapsed.adjacent[node.1, node.2]){
                pairwise.distances <- c(pairwise.distances, split.distances[node.1, node.2])
              }
            }
          }
          
          min.distance.matrix[pat.1, pat.2] <- min(pairwise.distances)
          
          if (verbose & !no.progress.bars) {
            setTxtProgressBar(progress.bar, count/total.host.pairs)
          }
        }
      }
    }
  }
  
  if (verbose & !no.progress.bars) close(progress.bar)
  
  normalisation.constant <- tree.info$normalisation.constant
  
  normalised.distance.matrix<- min.distance.matrix/normalisation.constant
  
  adjacency.df <- convert.to.columns(adjacency.matrix, hosts.included)
  contiguity.df <- convert.to.columns(contiguity.matrix, hosts.included)
  dir.12.df <- convert.to.columns(dir.12.matrix, hosts.included)
  dir.21.df <- convert.to.columns(dir.21.matrix, hosts.included)
  nodes.1.df <- convert.to.columns(nodes.1.matrix, hosts.included)
  nodes.2.df <- convert.to.columns(nodes.2.matrix, hosts.included)
  path.df <- convert.to.columns(path.matrix, hosts.included)
  min.distance.df <- convert.to.columns(min.distance.matrix, hosts.included)
  if(normalisation.constant!=1){
    normalised.distance.df <- convert.to.columns(normalised.distance.matrix, hosts.included)
  }
  
  classification <- cbind(adjacency.df, contiguity.df[,3], dir.12.df[,3], dir.21.df[,3], nodes.1.df[,3], nodes.2.df[,3], path.df[,3], min.distance.df[,3])
  
  column.names <- c("host.1", "host.2", "adjacent", "contiguous", "paths12", "paths21", "nodes1", "nodes2", "path.classification", "min.distance.between.subgraphs")
  
  if(normalisation.constant!=1){
    classification <- cbind(classification, normalised.distance.df[,3])
    column.names <- c(column.names, "normalised.min.distance.between.subgraphs")
  }
  
  colnames(classification) <- column.names
  return(list(classification = classification, collapsed=tt))
}

convert.to.columns <- function(matrix, names){
  a.table <- as.table(matrix)
  colnames(a.table) <- names
  rownames(a.table) <- names
  a.df <- as.data.frame(a.table, stringsAsFactors=F)
  keep <- complete.cases(a.df)
  a.df <- a.df[keep,]
  return(a.df)
}

#' @keywords internal
#' @export merge.classifications



merge.classifications <- function(all.tree.info, allow.mt = T, verbose = F){
  tt	<-	lapply(all.tree.info, function(tree.info){
    if(is.null(tree.info$classification.results$classification) & is.null(tree.info$classification.file.name)){
      NULL
    }
    
    if(!is.null(tree.info$classification.results$classification)){
      tt <- as.data.table(tree.info$classification.results$classification)
    } else {
      if (verbose) cat("Reading window input file ",tree.info$classification.file.name,"\n", sep="")
      tt <- as.data.table(read.table(tree.info$classification.file.name, sep=",", header=TRUE, stringsAsFactors=FALSE))
    }
    tt[, ID:=tree.info$id]			
    tt
  })	
  
  if(length(tt)==0){
    stop("No classification results present in any window; cannot continue.\n")
  }
  
  # need to worry about transmissions listed in the wrong direction in the input file
  # to avoid very large memory, consolidate by FILE
  
  if(verbose) cat("Rearranging host pairs...\n")
  
  tt	<- lapply(tt, function(x){
    #x	<- tt[[1]]
    tmp	<- copy(x)
    setnames(tmp, c('host.1','host.2','paths12','paths21'), c('host.2','host.1','paths21','paths12'))
    set(tmp, tmp[, which(path.classification=="anc")], 'path.classification', 'TMP')
    set(tmp, tmp[, which(path.classification=="desc")], 'path.classification', 'anc')
    set(tmp, tmp[, which(path.classification=="TMP")], 'path.classification', 'desc')
    set(tmp, tmp[, which(path.classification=="multiAnc")], 'path.classification', 'TMP')
    set(tmp, tmp[, which(path.classification=="multiDesc")], 'path.classification', 'multiAnc')
    set(tmp, tmp[, which(path.classification=="TMP")], 'path.classification', 'multiDesc')
    x	<- rbind(x, tmp)
    setkey(x, host.1, host.2)
    subset(x, host.1 < host.2)
  })
  #
  # rbind consolidated files
  #
  
  if(verbose) cat("Consolidating file contents...\n")
  
  tt	<- do.call('rbind',tt)
  
  if(verbose) cat("Finding patristic distance columns...\n")
  
  # reset names depending on which Classify script was used
  if(any('normalised.min.distance.between.subgraphs'==colnames(tt))){
    setnames(tt, 'normalised.min.distance.between.subgraphs', 'PATRISTIC_DISTANCE')
    set(tt, NULL, 'min.distance.between.subgraphs', NULL)
  } else if(any('min.distance.between.subgraphs'==colnames(tt))){
    setnames(tt, 'min.distance.between.subgraphs', 'PATRISTIC_DISTANCE')
  }
  
  setnames(tt, c('host.1','host.2','path.classification','paths21','paths12','adjacent','contiguous'), c('HOST.1','HOST.2','TYPE','PATHS.21','PATHS.12','ADJACENT','CONTIGUOUS'))
  
  # change type name depending on allow.mt
  if(!allow.mt){
    if(verbose) cat("Allowing only single lineage transmission to be used to infer directionality\n")
    set(tt, tt[, which(TYPE%in%c("multiAnc", "multiDesc"))], 'TYPE', 'complex')
  }
  
  #	check we have patristic distances, paths
  
  stopifnot( !nrow(subset(tt, is.na(PATRISTIC_DISTANCE))) )
  stopifnot( !nrow(subset(tt, is.na(PATHS.12))) )
  stopifnot( !nrow(subset(tt, is.na(PATHS.21))) )
  
  #	set to numeric
  set(tt, NULL, 'PATRISTIC_DISTANCE', tt[, as.numeric(PATRISTIC_DISTANCE)])
  # 	add window coordinates
  
  #	rename TYPE
  set(tt, tt[,which(TYPE=='anc')], 'TYPE', 'anc_12')
  set(tt, tt[,which(TYPE=='desc')], 'TYPE', 'anc_21')
  set(tt, tt[,which(TYPE=='multiAnc')], 'TYPE', 'multi_anc_12')
  set(tt, tt[,which(TYPE=='multiDesc')], 'TYPE', 'multi_anc_21')
  
  if(verbose) cat("Reordering...\n")
  
  #	reorder
  setkey(tt, ID, HOST.1, HOST.2)
  
  return(tt)
}

#' @keywords internal
#' @export summarise.classifications

summarise.classifications <- function(all.tree.info, min.threshold, dist.threshold, allow.mt = T, close.sib.only = F, verbose = F, contiguous = F){
  
  tt <- merge.classifications(all.tree.info, allow.mt, verbose)
  
  if (verbose) cat("Making summary output table...\n")
  
  set(tt, NULL, c('PATHS.12','PATHS.21'),NULL)
  
  existence.counts <- tt[, list(both.exist=length(ID)), by=c('HOST.1','HOST.2')]
  
  tt <- merge(tt, existence.counts, by=c('HOST.1', 'HOST.2'))
  
  if(!close.sib.only){
    if(!contiguous){
      tt.close <- tt[which(tt$ADJACENT & tt$PATRISTIC_DISTANCE < dist.threshold ),]
    } else {
      tt.close <- tt[which(tt$CONTIGUOUS & tt$PATRISTIC_DISTANCE < dist.threshold ),]
    }
  } else {
    if(!contiguous){
      tt.close <- tt[which(tt$ADJACENT & (tt$TYPE != "none" | tt$PATRISTIC_DISTANCE < dist.threshold)),]
    } else {
      tt.close <- tt[which(tt$CONTIGUOUS & (tt$TYPE != "none" | tt$PATRISTIC_DISTANCE < dist.threshold)),]
    }
  }
  
  # How many windows have this relationship, ADJACENT and PATRISTIC_DISTANCE below the threshold?
  type.counts	<- tt.close[, list(trees.with.this.relationship=length(ID)), by=c('HOST.1','HOST.2','TYPE')]
  # How many windows have ADJACENT and PATRISTIC_DISTANCE below the threshold?
  any.counts  <- tt.close[, list(trees.with.any.relationship=length(ID)), by=c('HOST.1','HOST.2')]

  tt.close		<- merge(tt.close, type.counts, by=c('HOST.1','HOST.2','TYPE'))
  tt.close		<- merge(tt.close, any.counts, by=c('HOST.1','HOST.2'))

  tt.close[, fraction:=paste(trees.with.this.relationship,'/',both.exist,sep='')]
  
  #	convert "anc_12" and "ans_21" to "anc" depending on direction
  tt.close[, DUMMY:=NA_character_]
  tmp			<- tt.close[, which(TYPE=="anc_12")]
  set(tt.close, tmp, 'TYPE', "trans")
  tmp			<- tt.close[, which(TYPE=="anc_21")]
  set(tt.close, tmp, 'DUMMY', tt.close[tmp, HOST.1])
  set(tt.close, tmp, 'HOST.1', tt.close[tmp, HOST.2])
  set(tt.close, tmp, 'HOST.2', tt.close[tmp, DUMMY])
  set(tt.close, tmp, 'TYPE', "trans")
  
  tmp			<- tt.close[, which(TYPE=="multi_anc_12")]
  set(tt.close, tmp, 'TYPE', "multi_trans")
  tmp			<- tt.close[, which(TYPE=="multi_anc_21")]
  set(tt.close, tmp, 'DUMMY', tt.close[tmp, HOST.1])
  set(tt.close, tmp, 'HOST.1', tt.close[tmp, HOST.2])
  set(tt.close, tmp, 'HOST.2', tt.close[tmp, DUMMY])
  set(tt.close, tmp, 'TYPE', "multi_trans")
  
  tt.close[, DUMMY:=NULL]
  
  set(tt.close, NULL, c('ID', 'ADJACENT','PATRISTIC_DISTANCE', "CONTIGUOUS", "nodes1", "nodes2"), NULL)
  tt.close <- tt.close[!duplicated(tt.close),]
  
  setkey(tt.close, HOST.1, HOST.2, TYPE)
  
  setnames(tt.close, c('HOST.1','HOST.2','TYPE', 'fraction', 'trees.with.this.relationship', 'trees.with.any.relationship'), 
           c('host.1', 'host.2', 'ancestry', 'fraction', 'ancestry.tree.count', 'any.relationship.tree.count'))
  
  setcolorder(tt.close, c('host.1', 'host.2', 'ancestry', 'ancestry.tree.count', 'both.exist', 'fraction', 'any.relationship.tree.count'))
  
  return(subset(tt.close, any.relationship.tree.count>min.threshold))
}

#' @keywords internal
#' @importFrom network as.network.matrix
#' @export simplify.summary

simplify.summary <- function(summary, arrow.threshold, total.trees, plot = F){
  
  done <- rep(FALSE, nrow(summary))
  
  for(line in 1:nrow(summary)){
    if(!done[line]){
      forwards.rows <- which(summary$host.1 == summary$host.1[line] & summary$host.2 == summary$host.2[line])
      backwards.rows <- which(summary$host.1 == summary$host.2[line] & summary$host.2 == summary$host.1[line])
      
      done[forwards.rows] <- T
      done[backwards.rows] <- T
      
      summary$ancestry[intersect(forwards.rows, which(summary$ancestry=="trans"))] <- "trans12"
      summary$ancestry[intersect(forwards.rows, which(summary$ancestry=="multi_trans"))] <- "multi_trans12"
      
      summary$ancestry[intersect(backwards.rows, which(summary$ancestry=="trans"))] <- "trans21"
      summary$ancestry[intersect(backwards.rows, which(summary$ancestry=="multi_trans"))] <- "multi_trans21"
      
      summary[backwards.rows,c(1,2)] <- summary[backwards.rows,c(2,1)]
    }
  }
  
  summary.wide <- reshape(summary, direction="w", idvar=c("host.1", "host.2"), timevar = "ancestry", v.names = "ancestry.tree.count",  drop=c("fraction"))
  
  summary.wide[is.na(summary.wide)] <- 0

  summary.wide$total.equiv <- 0
  
  if("ancestry.tree.count.none" %in% colnames(summary.wide)){
    summary.wide$total.equiv <- summary.wide$total.equiv + summary.wide$ancestry.tree.count.none
  }
  if("ancestry.tree.count.complex" %in% colnames(summary.wide)){
    summary.wide$total.equiv <- summary.wide$total.equiv + summary.wide$ancestry.tree.count.complex
  }
  
  if("ancestry.tree.count.trans12" %in% colnames(summary.wide)){
    summary.wide$total.12 <- summary.wide$ancestry.tree.count.trans12
  } else {
    summary.wide$total.12 <- rep(0, nrow(summary.wide))
  }
  
  if(!is.null(summary.wide$ancestry.tree.count.multi_trans12)){
    summary.wide$total.12 <- summary.wide$total.12 + summary.wide$ancestry.tree.count.multi_trans12
  }
  
  
  if("ancestry.tree.count.trans21" %in% colnames(summary.wide)){
    summary.wide$total.21 <- summary.wide$ancestry.tree.count.trans21
  } else {
    summary.wide$total.21 <- rep(0, nrow(summary.wide))
  }
  
  if(!is.null(summary.wide$ancestry.tree.count.multi_trans21)){
    summary.wide$total.21 <- summary.wide$total.21 + summary.wide$ancestry.tree.count.multi_trans21
  }
  
  summary.wide$total <- summary.wide$total.21 + summary.wide$total.12 + summary.wide$total.equiv

  dir <- summary.wide$total.12 >= arrow.threshold*total.trees | summary.wide$total.21 >= arrow.threshold*total.trees
  
  summary.wide$arrow[!dir] <- "none"
  
  if(length(which(dir)>0)){
    summary.wide$arrow[dir] <- sapply(which(dir), function(x)  if(summary.wide$total.12[x]>summary.wide$total.21[x]) "forwards" else "backwards")
    summary.wide$label[dir] <- paste0(round(pmax(summary.wide$total.12[dir],summary.wide$total.21[dir])/total.trees, 2),"/", round(summary.wide$total[dir]/total.trees, 2))
  }

  summary.wide$label[!dir] <- as.character(round(summary.wide$total[!dir]/total.trees, 2))
  
  out.table <- summary.wide[,c("host.1","host.2","arrow","label")]

  out.table[which(out.table$arrow=="backwards"),c(1,2)] <- out.table[which(out.table$arrow=="backwards"),c(2,1)] 
  
  out.table$arrow <- out.table$arrow!="none"
  
  out <- list(simp.table = out.table)
  
  if(plot){
    network.obj <- as.network.matrix(out.table[,c(1,2)], matrix.type = "edgelist")
    
    arrangement <- ggnet2(network.obj)$data[,c("label", "x", "y")]

    out.table$x.start <- sapply(out.table$host.1, function(x) arrangement$x[match(x, arrangement$label)]) 
    out.table$y.start <- sapply(out.table$host.1, function(x) arrangement$y[match(x, arrangement$label)]) 
    out.table$x.end <- sapply(out.table$host.2, function(x) arrangement$x[match(x, arrangement$label)]) 
    out.table$y.end <- sapply(out.table$host.2, function(x) arrangement$y[match(x, arrangement$label)]) 
    out.table$x.midpoint <- (out.table$x.end + out.table$x.start)/2
    out.table$y.midpoint <- (out.table$y.end + out.table$y.start)/2
  
    out.diagram <- ggplot() + 
      geom_segment(data=out.table[which(out.table$arrow),], aes(x=x.start, xend = x.end, y=y.start, yend = y.end), arrow = arrow(length = unit(0.01, "npc"), type="closed"), col="steelblue3", size=1.5, lineend="round") +
      geom_segment(data=out.table[which(!out.table$arrow),], aes(x=x.start, xend = x.end, y=y.start, yend = y.end), col="chartreuse3", size=1.5, lineend="round") +
      geom_label(aes(x=arrangement$x, y=arrangement$y, label=arrangement$label), alpha=0.25, fill="darkgoldenrod3") + 
      geom_text(data=out.table, aes(x=x.midpoint, y=y.midpoint, label=label)) + 
      theme_void()
    
    out$simp.diagram <- out.diagram
  }
  
  return(out)
}
