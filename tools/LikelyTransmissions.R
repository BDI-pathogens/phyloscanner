command.line <- T

list.of.packages <- c("phangorn", "argparse", "phytools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")

if(command.line){
  library(argparse)
  
  arg_parser = ArgumentParser(description="Identify the likely direction of transmission between all patients present in a phylogeny.")
  arg_parser$add_argument("-x", "--tipRegex", action="store", default="^(.*)_read_([0-9]+)_count_([0-9]+)$", 
                          help="Regular expression identifying tips from the dataset. Three groups: patient ID, read ID, and read count. If absent, input will be assumed to be from the phyloscanner pipeline, and the patient ID will be the BAM file name. Necessary only if -s is specified.")
  arg_parser$add_argument("-r", "--outgroupName", action="store", help="Label of tip to be used as outgroup (if unspecified, tree will be assumed to be already rooted).")
  arg_parser$add_argument("-z", "--zeroLengthTipsCount", action="store_true", default=FALSE, help="If present, a zero length terminal branch associates the patient at the tip with its parent node, interrupting any inferred transmission pathway between another pair of hosts that goes through the node.")
  arg_parser$add_argument("-t", "--dualInfectionThreshold", type="double", help="Length threshold at which a branch of the subtree constructed just from reads from a single patient is considered long enough to indicate a dual infection (such a patient will be ignored, at present). If absent, keep all patients.")
  arg_parser$add_argument("-s", "--romeroSeverson", default=FALSE, action="store_true", help="If present, the Romero-Severson classification will be repeated on tree to annotate internal nodes beyond the MRCAs of each split. Recommended if R-S was used to determine splits.")
  arg_parser$add_argument("treeFileName", action="store", help="Tree file name")
  arg_parser$add_argument("splitsFileName", action="store", help="Splits file name")
  arg_parser$add_argument("outputFileName", action="store", help="Output file name (.csv format)")
  arg_parser$add_argument("-D", "--scriptdir", action="store", help="Full path of the script directory.", default="/Users/twoseventwo/Documents/phylotypes/")
  
  args <- arg_parser$parse_args()
  
  tree.file.name <- args$treeFileName
  splits.file.name <- args$splitsFileName
  script.dir <- args$scriptdir
  output.name <- args$outputFileName
  root.name <- args$outgroupName
  tip.regex <- args$tipRegex
  split.threshold <- args$dualInfectionThreshold
  if(is.null(split.threshold)){
    split.threshold <- Inf
  }
  zero.length.tips.count <- args$zeroLengthTipsCount
  romero.severson <- args$romeroSeverson
  
} else {
  setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160517_clean/")
  
  tree.file.name <- "RAxML_bestTree.InWindow_2550_to_2900.tree"
  splits.file.name <- "Subtrees_r_run20160517_inWindow_2550_to_2900.csv"
  output.name <- "LikelyTransmissions.InWindow_2550_to_2900.csv"
  root.name <- "C.BW.00.00BW07621.AF443088"
  split.threshold <- 0.08
  zero.length.tips.count <- F
  romero.severson <- T
  tip.regex <- "^(.*)-[0-9].*_read_([0-9]+)_count_([0-9]+)$"
}

library(phytools)
library(phangorn)

source(file.path(script.dir, "TransmissionUtilityFunctions.R"))
source(file.path(script.dir, "SubtreeMethods.R"))

cat("Opening file: ", tree.file.name, "...\n", sep = "")

tree <-read.tree(tree.file.name)
tree <- root(tree, root.name)

cat("Making tree multiphyletic...\n")

tree <- di2multi(tree, tol = 1E-5)

cat("Reading splits file...\n")

splits <- read.csv(splits.file.name, stringsAsFactors = F)

cat("Collecting tips for each patient...\n")

patients <- unique(splits$orig.patients)
patient.tips <- lapply(patients, function(x) which(tree$tip.label %in% splits[which(splits$orig.patients==x), "tip.names"]))
names(patient.tips) <- patients

if(romero.severson){
  
  # Really the only reason to load the splits in is to avoid needing the blacklist yet again. Maybe rethink this.
  
  patient.mrcas <- lapply(patient.tips, function(node) mrca.phylo.or.unique.tip(tree, node, zero.length.tips.count))
  
  tip.assocs.pats <- lapply(tree$tip.label, function(x) if(x %in% splits$tip.names) return(splits$orig.patients[which(splits$tip.names==x)]) else return("*"))
  
  new.assocs <- list()
  
  classify.down(getRoot(tree), tree, tip.assocs.pats, patient.mrcas, vector(), tip.regex)
  
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
  
  splits.for.patients <- list()
  patients.for.splits <- list()
  
  #find the splits and reannotate the tree
  
  for(pat.no in seq(1, length(patients))){
    patient <- patients[pat.no]
    no.splits <- counts.by.patient[pat.no]
    if(no.splits > 1){
      splits.for.patients[[patient]] <- paste(patient,"-S",seq(1, no.splits),sep="")
    } else {
      splits.for.patients[[patient]] <- patient
    }
    for(split.id in splits.for.patients[[patient]]){
      patients.for.splits[[split.id]] <- patient
    }
    
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
  
  split.ids <- patients.copy
  
  was.split <- unique(patients[!(patients %in% patients.copy)])
  
  split.tips <- patient.tips.copy
  assocs <- split.assocs
  
} else {
  cat("Collecting tips for each split...\n")
  
  all.splits <- unique(splits$patient.splits)
  split.tips <- lapply(all.splits, function(x) which(tree$tip.label %in% splits[which(splits$patient.splits==x), "tip.names"]))
  names(split.tips) <- all.splits
  
  tip.assocs.splits <- lapply(tree$tip.label, function(x) if(x %in% splits$tip.names) return(splits$patient.splits[which(splits$tip.names==x)]) else return("*"))
  
  split.mrcas <- lapply(split.tips, function(node) mrca.phylo.or.unique.tip(tree, node, zero.length.tips.count))
  
  node.assocs <- annotate.internal(tree, all.splits, split.tips, split.mrcas)
  
  splits.for.patients <- lapply(patients, function(x) unique(splits$patient.splits[which(splits$orig.patients==x)] ))
  names(splits.for.patients) <- patients
  
  patients.for.splits <- lapply(all.splits, function(x) unique(splits$orig.patients[which(splits$patient.splits==x)] ))
  names(patients.for.splits) <- all.splits
  
  was.split <- unique(splits$orig.patients[which(splits$orig.patients!=splits$patient.splits)])
  
  assocs <- node.assocs$details
  
}

# get rid of what look like dual infections

ignore.because.split <- vector()

for(split.patient in was.split){
  tips <- patient.tips[[split.patient]]
  subtree <- drop.tip(phy = tree, tip = tree$tip.label[-patient.tips[[split.patient]]])
  max.length <- max(subtree$edge.length)
  if(max.length > split.threshold){
    ignore.because.split <- c(ignore.because.split, split.patient)
  }
}

patients.included <- setdiff(patients, ignore.because.split)

total.pairs <- (length(patients.included) ^ 2 - length(patients.included))/2

cat("Collapsing subtrees...\n")

tt <- output.trans.tree(tree, assocs)

cat("Testing pairs\n")

count <- 0
direct.descendant.matrix <- matrix(NA, length(patients.included), length(patients.included))
for(pat.1 in seq(1, length(patients.included))){
  for(pat.2 in  seq(1, length(patients.included))){
    if (pat.1 < pat.2) {
      count <- count + 1
      
      pat.1.id <- patients.included[pat.1]
      pat.2.id <- patients.included[pat.2]
      
      all.nodes <-  c(splits.for.patients[[pat.1.id]], splits.for.patients[[pat.2.id]])
      
      #we want that they form a perfect blob with no intervening nodes for any other patients
      OK <- TRUE
      
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
      
      # then we want to see if any node from one patient interrupts a path between two nodes from the other
      if(OK){
        entangled <- F
        
        if(length(splits.for.patients[[pat.1.id]])>1){
          combns.1 <- t(combn(splits.for.patients[[pat.1.id]], 2))
          for(comb in seq(1, nrow(combns.1))){
            path <- get.tt.path(tt, combns.1[comb, 1], combns.1[comb, 2])
            for(node in path){
              if(!startsWith(node, "none")){
                if(patients.for.splits[[node]]==pat.2.id){
                  entangled <- T
                  break
                }
              }
            }
            if(entangled){
              break
            }
          }
        }
        if(!entangled & length(splits.for.patients[[pat.2.id]])>1){
          combns.2 <- t(combn(splits.for.patients[[pat.2.id]], 2))
          for(comb in seq(1, nrow(combns.2))){
            path <- get.tt.path(tt, combns.2[comb, 1], combns.2[comb, 2])
            for(node in path){
              if(!startsWith(node, "none")){
                if(patients.for.splits[[node]]==pat.1.id){
                  entangled <- T
                  break
                }
              }
            }
            if(entangled){
              break
            }
          }
        }
        rel.determined <- F
        
        if(!entangled){
          for(pat.1.splt in splits.for.patients[[pat.1.id]]){
            ancestors <- get.tt.ancestors(tt, pat.1.splt)
            if(length(intersect(ancestors, splits.for.patients[[pat.2.id]]))>0){
              direct.descendant.matrix[pat.1, pat.2] <- "desc"
              rel.determined <- T
              break
            }
          }
          if(!rel.determined){
            for(pat.2.splt in splits.for.patients[[pat.2.id]]){
              ancestors <- get.tt.ancestors(tt, pat.2.splt)
              if(length(intersect(ancestors, splits.for.patients[[pat.1.id]]))>0){
                direct.descendant.matrix[pat.1, pat.2] <- "anc"
                rel.determined <- T
                break
              }
            }
          }
          if(!rel.determined){
            direct.descendant.matrix[pat.1, pat.2] <- "sib"
          }
        } else {
          current.node <- all.nodes[1]
          for(node in all.nodes[2:length(all.nodes)]){
            current.node <- get.tt.mrca(tt, current.node, all.nodes[2])
          }
          if(startsWith(current.node, "none")){
            direct.descendant.matrix[pat.1, pat.2] <- "trueInt"
          } else if(patients.for.splits[[current.node]]==pat.1.id){
            direct.descendant.matrix[pat.1, pat.2] <- "intAnc"
          } else if(patients.for.splits[[current.node]]==pat.2.id){
            direct.descendant.matrix[pat.1, pat.2] <- "intDesc"
          } else {
            stop("Huh?")
          }
        }
      }
      
      if (count %% 100 == 0) {
        cat(
          paste(
            "Done ",format(count, scientific = F)," of ",format(total.pairs, scientific = F)," pairwise calculations\n", sep = ""
          )
        )
      }
    }
  }
}

direct.descendant.table <- as.table(direct.descendant.matrix)

colnames(direct.descendant.table) <- patients.included
rownames(direct.descendant.table) <- patients.included


dddf <- as.data.frame(direct.descendant.table)
dddf <- dddf[complete.cases(dddf),]

colnames(dddf) <- c("Patient_1", "Patient_2", "Relationship")

write.table(dddf, file = output.name, sep = ",", row.names = FALSE, col.names = TRUE, quote=F)