library(phangorn)

# Get the tip number of the named taxon

get.tip.no <- function(tree, name){
  return(match(name, tree$tip.label))
}

# Get all tips associated with a patient

get.tips.for.patient <- function(patient.string, tree, patient.ids){
  node.numbers <- which(patient.ids==patient.string)
  return(node.numbers)
}

# Edge length by node number (probably in some library somewhere)

get.edge.length <- function(node, tree){
  index <- which(tree$edge[,2]==node)
  return(tree$edge.length[index])
}

# Node is a tip (likewise)

is.tip <- function(node, tree){
  return(node <= length(tree$tip.label))
}

# Get the sequence of ancestors from one node to the other

get.ancestral.sequence <- function(desc, anc, tree){
  if(!(anc %in% Ancestors(tree, desc, type="all"))){
    stop("anc is not an ancestor of desc at all")
  }
  current <- desc
  out <- desc
  
  while(current!=anc){
    
    current <- Ancestors(tree, current, type="parent")
    
    out <- c(out, current)
    
  }
  return(out)
}

get.patient.from.tip <- function(node, tree){
  # is it a tip?
  if(!(is.tip(node,tree))){
    stop("Not a tip")
  }
  patient.id <- strsplit(tree$tip.label[node], "_")[[1]][1]
  return(patient.id)
  
}

# Output the set of patients whose tips are descended from node (or the host of the node itself if a tip)

get.patients.from.this.clade <- function(node, tree, patient.tips){
  if(is.tip(node,tree)){
    return(get.patient.from.tip(node, tree))
  }
  
  out <- vector()
  tips <- Descendants(tree, node, type="tips")[[1]]
  for(tip in tips){
    if(!(get.patient.from.tip(tip,tree) %in% out)){
      out <- c(out, get.patient.from.tip(tip, tree))
    }
  }
  return(out)
}

# Is desc _unambiguously_ a descendant of anc? I.e. is the MRCA node of desc a descendant of the MRCA
# node of anc, but no tips of anc are descended from the MRCA node of desc?

is.descendant.of <- function(desc, anc, tree, mrca.list, tip.list){
  mrca.desc <- mrca.list[[desc]]
  mrca.anc <- mrca.list[[anc]]
  
  if(!(mrca.anc %in% Ancestors(tree, mrca.desc, type="all"))){
    return(FALSE)
  } else {
    for(tip in tip.list[[anc]]){
      if(mrca.desc %in% Ancestors(tree, tip, type="all")){
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

is.direct.descendant.of <- function(desc, anc, tree, mrca.list, tip.list, all.mrcas){
  
  if(!is.descendant.of(desc,anc,tree,mrca.list,tip.list)){
    return(FALSE)
  }
  
  mrca.desc <- mrca.list[[desc]]
  mrca.anc <- mrca.list[[anc]]
  
  
  #check that there are no other MRCAs in the sequence of ancestors from desc to anc - this means that either
  #there is an intervening sampled host, or one or the other MRCA is shared with another patient and hence
  #situation is too ambiguous
  
  anc.sequence <- get.ancestral.sequence(mrca.desc, mrca.anc, tree)
  
  anc.nodes <- unique(unlist(lapply(anc.sequence, all.nodes.in.polytomy, tree=tree, include.tips=TRUE)))
  
  
  
  for(node in anc.nodes){
    children <- Children(tree, node)
    if(length(children)>0){ 
      patients.1 <- get.patients.from.this.clade(children[1], tree, tip.list)
      patients.2 <- get.patients.from.this.clade(children[2], tree, tip.list)
      
      patients.1 <- patients.1[! patients.1 %in% c(patients[anc], patients[desc])]
      patients.2 <- patients.2[! patients.2 %in% c(patients[anc], patients[desc])]
      
      if(length(intersect(patients.1,patients.2))>0){
        return(FALSE)
      }
    } else {
      # if node is a tip, then need to check it isn't associated with another host
      patient <- get.patient.from.tip(node, tree)
      if(patient!=patients[anc] & patient!=patients[desc]){
        return(FALSE)
      }
    }
  }
  
  return(TRUE)
}

# Output the set of patients whose mrca is this node

get.patients.with.these.mrcas <- function(node, tree, patient.mrcas){
  out <- vector()
  mrca.vec <- sapply(patient.mrcas, "[[", 1)
  
  numbers <- which(mrca.vec==node)
  
  return(names(numbers))
}

# Need a MRCA function which if given a single tip, returns that tip rather than NA. Also, since the tree is basically
# multifurcating but cheats with very short branch lengths, the MRCA node should the last one that can be reached by moving
# up 1E-6 from the MRCA than phangorn gives.

mrca.phylo.or.unique.tip <- function(x, node){
  #  cat("Node set: ",node,"\n")
  if(length(node)==1){
    #    cat("Node ",node," is unique for this patient.\n")
    mrca = node
  } else {
    mrca = mrca.phylo(x, node)
    #    cat("Node ",mrca," is the MRCA of this set.\n")
  }
  too.long <- FALSE
  while(!too.long){
    if(mrca==getRoot(tree)){
      too.long <- TRUE
    }else if(get.edge.length(mrca, tree)>1E-5){
      too.long <- TRUE
    } else {
      mrca <- Ancestors(tree, mrca, type="parent")
    }
  }
  return(mrca)
}

#return all of the bifuracting nodes in what looks like a polytomy. If include.tips = TRUE, also include any tips
#attached to the "polytomy" by "zero" length branches.

all.nodes.in.polytomy <- function(tree, node, include.tips=FALSE, upwards=TRUE){
  
  out <- vector()
  
  if(upwards){
    too.long <- FALSE
    while(!too.long){
      if(node==getRoot(tree)){
        too.long <- TRUE
        mrca <- node 
      } else if(get.edge.length(node, tree)>1E-5) {
        too.long <- TRUE
        mrca <- node 
      } else {
        node <- Ancestors(tree, node, type="parent")
      }
    }
  } else {
    mrca <- node
  }
  
  out <- mrca
  
  children <- Children(tree, mrca)
  
  for(child in children){
    if(get.edge.length(child, tree)<1E-5){
      if(include.tips | !is.tip(child, tree)){
        out <- c(out, all.nodes.in.polytomy(tree, child, include.tips, FALSE))
      }
    }
  }
  return(out)
}

setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160111/")

for(start in seq(501, 9301, by=100)){
  
  end <- start+349 
  
  if(file.exists(paste("RAxML_bestTree/RAxML_bestTree.InWindow_",start,"_to_",end,".tree",sep=""))){
    
    cat("Opening file: RAxML_bestTree.InWindow_",start,"_to_",end,".tree\n", sep="")
    
    tree <- read.tree(paste("RAxML_bestTree/RAxML_bestTree.InWindow_",start,"_to_",end,".tree",sep=""))
    
    root.name <- "Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
    
    tree <- root(tree, root.name)
    
    tip.labels <- tree$tip.label
    
    outgroup.no <- which(substr(tip.labels, 1, 3)=="Ref")
    
    tip.labels[outgroup.no] <- "Ref_read_NA_count_NA"
    
    split.tip.labels <- strsplit(tip.labels, "_")
    
    patient.ids <- sapply(split.tip.labels, "[[", 1)
    
    get.last <- function(vec) return(vec[length(vec)])
    
    read.counts <- as.numeric(sapply(split.tip.labels, get.last))
    
    depths <- node.depth.edgelength(tree)
    
    patients <- unique(patient.ids)
    patients <- patients[which(patients!="Ref")]
    
    patient.tips <- sapply(patients, get.tips.for.patient, tree=tree, patient.ids = patient.ids)
    
    patient.mrcas <- lapply(patient.tips, function(node) mrca.phylo.or.unique.tip(tree, node))
    
    all.mrca.nodes <- unique(unlist(patient.mrcas))
    
    total.pairs <- length(patients)^2 - length(patients)
    
    count <- 0
    direct.descendant.matrix <- matrix(NA, length(patients), length(patients))
    for(desc in seq(1, length(patients))){
      for(anc in seq(1, length(patients))){
        if(desc!=anc){
          #     cat("Testing if ",patients[desc]," is directly descended from ",patients[anc],".\n")
          
          count <- count + 1
    
          direct.descendant.matrix[desc, anc] <- is.direct.descendant.of(desc, anc, tree, patient.mrcas, patient.tips, all.mrca.nodes)
          if(direct.descendant.matrix[desc, anc]){
            cat(patients[desc],"is a descendant of",patients[anc],"\n")
          }
          if(count %% 100==0){
            cat(paste("Done ",count," of ",total.pairs," pairwise calculations\n", sep=""))
          }
        }
      }
    }
    
    direct.descendant.table <- as.table(direct.descendant.matrix)
    
    colnames(direct.descendant.table) <- patients
    rownames(direct.descendant.table) <- patients
    
    
    dddf = as.data.frame(direct.descendant.table)
    dddf = dddf[complete.cases(dddf),]
    
    colnames(dddf) <- c("Descendant", "Ancestor", "Present")
    
    dddf <- dddf[dddf$Present,]
    
    write.table(dddf, file=paste("transmissions/transmissions.InWindow_",start,"_to_",end,".csv", sep=""), sep=",", row.names=FALSE, col.names=TRUE)
  }
}