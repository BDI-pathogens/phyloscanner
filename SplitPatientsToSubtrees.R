command.line <- T

list.of.packages <- c("phangorn", "optparse", "phytools", "ggplot2", "gdata", "mvtnorm", "expm")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")

if(!("ggtree" %in% installed.packages()[,"Package"])){
  source("https://bioconductor.org/biocLite.R")
  biocLite("ggtree")
}

library(phangorn)
library(optparse)
library(phytools)
library(ggplot2)
library(ggtree)

source("/Users/twoseventwo/Documents/phylotypes/TransmissionUtilityFunctions.R")

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
  cat("Splitting ",patient,"...\n",sep="")
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
  cat("Finding splits for ",patient,"...\n",sep="")
  
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
  cat("Finding last nonconflicted ancestor for ", tip.label, "...\n", sep="")
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
  if(assocs[[node]]!="*"){
    if(parent==0 | assocs[[parent]]!=assocs[[node]]){
      patient <- assocs[[node]]
      pat.index <- which(patients == patient)
      counts.vec[pat.index] <- counts.vec[pat.index] + 1
      first.nodes.list[[patient]] <- c(first.nodes.list[[patient]], node)
    }
  }
  return(list(counts = counts.vec, first.nodes = first.nodes.list))
}

# Does the RS classification (todo move to TUF?)

classify.down <- function(node, tree, tip.assocs, patient.mrcas, blacklist){
  
  if(is.tip(tree, node)){
    if(node %in% blacklist | !(startsWith(tree$tip.label[node], "BEE"))){
      result <- "*"
    } else {
      result <- tip.assocs[[node]]
    }
  } else {
    
    child.assocs.vector <- sapply(Children(tree, node), classify.down, tree=tree, tip.assocs=tip.assocs, patient.mrcas = patient.mrcas, blacklist=blacklist)
    
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
  
  new.assocs[[node]] <<- result
  return(result)
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
        while(assocs[[current.node]]=="*"){
          run <- c(run, current.node)
          current.node <- Ancestors(tree, current.node, type="parent")
          if(current.node == 0){
            break
          }
        }
        out[[int.node]] <- run
      }
    }
  }
  return(out)
}


if(command.line){
  option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL, 
                help="Input tree file name", metavar="inputTreeFileName"),
    make_option(c("-o", "--outRoot"), type="character", default="out", 
                help="Output file name [default= %default]", metavar="outputFileName"),
    make_option(c("-r", "--refSeqName"), type="character", default=NULL, 
                help="Reference sequence label (if unspecified, tree will be assumed to be already rooted)"),
    make_option(c("-b", "--blacklists"), type="character", default=NULL, 
                help="Path to one or more .csv files listing tips to ignore; separated with colons"),
    make_option(c("-v", "--verbose"), type="logical", default=FALSE, 
                help="Talk about what I'm doing"),
    make_option(c("-x", "--tipRegex"), type="character", default=NULL, 
                help="Regular expression identifying tips from the dataset"),
    make_option(c("-z", "--zeroLengthTipsCount"), type="logical", default=FALSE, 
                help="If TRUE, a zero length terminal branch associates the parent at the tip with its parent
                node, interrupting any inferred transmission pathway between another pair of hosts that goes
                through the node."),
    make_option(c("-s", "--splitsRule"), type="character", default="c", help="The rules by which the sets
                of patients are split into groups in order to ensure that all groups can be members of 
                connected subtrees without causing conflicts. Currently available: c=conservative, 
                r=Romero-Severson")
    
  )
  
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  
  verbose <- opt$verbose
  zero.length.tips.count <- opt$zeroLengthTipsCount
  
  file.name <- opt$file
  out.root <- opt$outRoot
  blacklists.string <- opt$blacklist
  blacklist.files <- unlist(strsplit(blacklists.string, ":")) 
  root.name <- opt$refSeqName
  
  tip.regex <- opt$tipRegex
  
  mode <- opt$splitsRule
  if(!(mode %in% c("c", "r"))){
    stop("Unknown split classifier")
  }

} else {
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160517_clean/")
  file.name <- "RAxML_bestTree.InWindow_800_to_1150.tree"
  blacklist.files <-c("PatientBlacklist_InWindow_800_to_1150.csv", "ReadBlacklist_InWindow_800_to_1150.csv")
  out.root <- "BEEHIVE180716.InWindow_800_to_1150"
  root.name <- "C.BW.00.00BW07621.AF443088"
  tip.regex <- "^(.*)-[0-9].*_read_([0-9]+)_count_([0-9]+)$"
}

tree <- read.tree(file.name)
tree <- unroot(tree)
tree <- di2multi(tree, tol = 1E-5)
tree <- root(tree, outgroup = which(tree$tip.label==root.name), resolve.root = T)

tip.labels <- tree$tip.label

blacklist <- vector()

for(blacklist.file.name in blacklist.files){
  if(file.exists(blacklist.file.name)){
    blacklisted.tips <- read.table(blacklist.file.name, sep=",", header=F, stringsAsFactors = F)
    if(length(blacklisted.tips)>0){
      blacklist <- c(blacklist, sapply(blacklisted.tips, get.tip.no, tree=tree))
    }
  } else {
    warning(paste("File ",blacklist.file.name," does not exist; skipping.",paste=""))
  }
}

outgroup.no <- which(tip.labels == root.name)

patient.ids <- sapply(tip.labels, function(x) patient.from.label(x, tip.regex))

patients <- unique(patient.ids[-blacklist])

patients <- patients[!is.na(patients)]

patient.tips <-
  sapply(patients, function(x)  setdiff(which(patient.ids==x), blacklist))

patient.mrcas <- lapply(patient.tips, function(node) mrca.phylo.or.unique.tip(tree, node, zero.length.tips.count))

if(mode=="c"){

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
  
  cat("Drawing tree with conflicts displayed...\n")
  
  # list of tip statuses
  
  conflict.tips <- rep("normal", length(tree$tip.label))
  conflict.tips[which(patient.ids %in% patients.with.conflicts)] <- "conflict"
  conflict.tips[blacklist] <- "blacklist"
  
  # Need to make this the same length as the set of nodes for ggtree
  
  conflict.tips <- c(conflict.tips, rep(NA, tree$Nnode))
  
  # Only put geom_points on nodes with a conflict
  
  node.assoc.counts <- node.assocs$counts
  node.assoc.counts[which(node.assoc.counts<=1)] <- NA
  
  tree.display <- ggtree(tree) +
    geom_point2(shape = 21,aes(fill = node.assoc.counts, size=!is.na(node.assoc.counts))) +
    theme(legend.position="right") +
    scale_fill_gradient(low = "darkgreen", high="red") +
    theme(legend.title = element_text(size=40)) +
    theme(legend.text = element_text(size=40)) +
    theme(legend.key.size = unit(10, "lines")) +
    labs(size="Conflict exists", fill="Number of\nconflicting patients", col="Problem\npatient") +
    geom_tiplab(aes(col = conflict.tips)) + 
    scale_color_manual(values=c("blue", "red", "black"))
  tree.display
  
  ggsave(file=paste("tree_conflicts_",out.root,".pdf",sep=""), 
         height = 400, width = 100, limitsize = F)
  
  cat("Resolving patients with conflicts into separate groups...\n")
  
  # Copy the patients vector and patient tip list for splitting
  
  patients.copy <- patients
  patient.tips.copy <- patient.tips
  
  pat.length <- length(patients.copy)
  
  # Split each patient
  
  for(patient in patients){
    sp.res <- split.patient(tree, patient, patients,patient.tips, node.assocs$details, tip.regex)
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
  
  cat("Drawing tree...\n")
  
  temp.ca <- rep(NA, length(tree$tip.label) + tree$Nnode)
  
  for(item in seq(1, length(node.assocs$details))){
    if(!is.null(node.assocs$details[[item]])){
      temp.ca[item] <- node.assocs$details[[item]]
    }
  }
  
  temp.ca <- factor(temp.ca, levels = sample(levels(as.factor(temp.ca))))
  
  tree.display <- ggtree(tree) +
    geom_point(shape = 21, aes(fill = temp.ca)) +
    scale_fill_hue(na.value = "black") +
    theme(legend.position="none") +
    geom_tiplab(aes(col = conflict.tips)) + 
    scale_color_manual(values=c("blue", "red", "black"))
  
  tree.display
  
  ggsave(file=paste("tree_c_",out.root,".pdf",sep=""), 
         height = 400, width = 100, limitsize = F)
  
  
} else if(mode=="r"){
  
  cat("Applying the Romero-Severson classification to internal nodes...\n")
  
  tip.assocs <- annotate.tips(tree, patients, patient.tips)
  
  new.assocs <- list()
  
  classify.down(getRoot(tree), tree, tip.assocs, patient.mrcas, blacklist)
  
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
  
  #Now the splits. A new split is where you encounter a node with a different association to its parent
  
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
      patients.copy <- c(patients.copy, paste(patient,"_S",seq(1, no.splits),sep=""))
      
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
        new.name <- paste(patient,"_S", split.index, sep="")
        
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
  
  cat("Drawing tree...\n")
  
  temp.ca <- rep(NA, length(tree$tip.label) + tree$Nnode)
  
  for(item in seq(1, length(resolved.assocs))){
    if(split.assocs[[item]] != "*"){
      temp.ca[item] <- split.assocs[[item]]
    }
  }
  
  temp.ca <- factor(temp.ca, levels = sample(levels(as.factor(temp.ca))))
  
  tree.display <- ggtree(tree) +
    geom_point(shape = 21, aes(fill = temp.ca)) +
    scale_fill_hue(na.value = "black") +
    geom_tiplab() + 
    theme(legend.position="none")
  
  tree.display
  
  ggsave(file=paste("tree_r_",out.root,".pdf",sep=""), 
         height = 400, width = 100, limitsize = F)
  

}
  
cat("Writing output...\n")

orig.patients <- vector()
patient.splits <- vector()
tip.names <- vector()

for(patient.plus in patients.copy){
  tips <- patient.tips.copy[[patient.plus]]
  orig.patients <- c(orig.patients, rep(unlist(strsplit(patient.plus, "-S"))[1], length(tips)))
  patient.splits <- c(patient.splits, rep(patient.plus, length(tips)))
  tip.names <- c(tip.names, tree$tip.label[tips])
}

rs.subtrees <- data.frame(orig.patients, patient.splits, tip.names)

rs.file.name <- paste("Subtrees_",out.root,"_",mode,".csv", sep="")
write.csv(rs.subtrees, rs.file.name)


