library(phangorn)

verbose <- TRUE
# Get the tip number of the named taxon

zero.length.tips.count <- TRUE

get.tip.no <- function(tree, name) {
  return(match(name, tree$tip.label))
}

# Get all tips associated with a patient

get.tips.for.patient <- function(tree, patient.string, patient.ids) {
  node.numbers <- which(patient.ids == patient.string)
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
  patient.id <- strsplit(tree$tip.label[node], "_")[[1]][1]
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
  function(tree, desc, anc, mrca.list, tip.list) {
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
          cat("Duplicate found")
          return(FALSE)
        }
        
        if (zero.length.tips.count) {
          # if zero length tips count, then a tip attached to an internal node on the ancestral sequence by a branch
          # of length zero counts as an ancestor
          for (child in Children(tree, node)) {
            if (is.tip(tree, child) & get.edge.length(tree, child) < 1E-5) {
              if (!(get.patient.from.tip(tree, child) %in% c(patients[anc], patients[desc]))) {
                return(FALSE)
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

# Need a MRCA function which if given a single tip, returns that tip rather than NA.

mrca.phylo.or.unique.tip <- function(tree, node) {
  if (length(node) == 1) {
    #    cat("Node ",node," is unique for this patient.\n")
    if(!zero.length.tips.count){
      return(node)
    } else {
      length <- get.edge.length(tree, node)
      while(length<1E-5){
        node <- Ancestors(tree, node, type="parent")
        length <- get.edge.length(tree, node)
      }
      return(node)
    }
  } else {
    return(mrca.phylo(tree, node))
    #    cat("Node ",mrca," is the MRCA of this set.\n")
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

setwd(
  "/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160111/"
)

for (start in seq(501, 9301, by = 100)) {
  end <- start + 349
  
  if (file.exists(
    paste(
      "RAxML_bestTree/RAxML_bestTree.InWindow_",start,"_to_",end,".tree",sep =
      ""
    )
  )) {
    cat("Opening file: RAxML_bestTree.InWindow_",start,"_to_",end,".tree\n", sep =
          "")
    
    tree <-
      read.tree(
        paste(
          "RAxML_bestTree/RAxML_bestTree.InWindow_",start,"_to_",end,".tree",sep =
            ""
        )
      )
    
    root.name <- "Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
    
    tree <- root(tree, root.name)
    
    tree <- di2multi(tree, tol = 1E-5)
    
    tip.labels <- tree$tip.label
    
    outgroup.no <- which(substr(tip.labels, 1, 3) == "Ref")
    
    tip.labels[outgroup.no] <- "Ref_read_NA_count_NA"
    
    split.tip.labels <- strsplit(tip.labels, "_")
    
    patient.ids <- sapply(split.tip.labels, "[[", 1)
    
    get.last <- function(vec)
      return(vec[length(vec)])
    
    read.counts <- as.numeric(sapply(split.tip.labels, get.last))
    
    depths <- node.depth.edgelength(tree)
    
    patients <- unique(patient.ids)
    patients <- patients[which(patients != "Ref")]
    
    patient.tips <-
      sapply(patients, get.tips.for.patient, tree = tree, patient.ids = patient.ids)
    
    patient.mrcas <-
      lapply(patient.tips, function(node)
        mrca.phylo.or.unique.tip(tree, node))
    
  
    total.pairs <- length(patients) ^ 2 - length(patients)
  
    
    count <- 0
    direct.descendant.matrix <- matrix(NA, length(patients), length(patients))
    for (desc in seq(1, length(patients))) {
      for (anc in seq(1, length(patients))) {
        if (desc != anc) {
          if (verbose) {
            cat(
              "Testing if ",patients[desc]," is directly descended from ",patients[anc],".\n", sep =
                ""
            )
          }
          
          count <- count + 1
          
          
          if (!compareNA(direct.descendant.matrix[anc, desc],TRUE)) {
            direct.descendant.matrix[desc, anc] <-
              is.direct.descendant.of(tree, desc, anc, patient.mrcas, patient.tips)
            if (direct.descendant.matrix[desc, anc]) {
              cat(patients[desc],"is a descendant of",patients[anc],"\n")
            }
          } else {
            direct.descendant.matrix[desc, anc] <- FALSE
          }
          
          if (count %% 100 == 0) {
            cat(
              paste(
                "Done ",count," of ",total.pairs," pairwise calculations\n", sep = ""
              )
            )
          }
        }
      }
    }
    
    direct.descendant.table <- as.table(direct.descendant.matrix)
    
    colnames(direct.descendant.table) <- patients
    rownames(direct.descendant.table) <- patients
    
    
    dddf <- as.data.frame(direct.descendant.table)
    dddf <- dddf[complete.cases(dddf),]
    
    colnames(dddf) <- c("Descendant", "Ancestor", "Present")
    
    dddf <- dddf[dddf$Present,]
    
    write.table(
      dddf, file = paste(
        "transmissions/transmissions.InWindow_",start,"_to_",end,"_altTest.csv", sep = ""
      ), sep = ",", row.names = FALSE, col.names = TRUE
    )
  }
}