library(phangorn)
library(optparse)
library(phytools)

get.last <- function(vec)
  return(vec[length(vec)])

split.patient <- function(tree, id, split.threshold, patients, patient.tips, patient.mrcas, split.counts, 
                          original.id=NULL, already.split=FALSE){
  
  patients <<- patients
  patient.tips <<- patient.tips
  patient.mrcas <<- patient.mrcas
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
        cat("Splitting ",id," due to exceeding longest branch length threshold\n", sep="")
      }

      if(length(patient.subtree$edge.length)>=3){
        
        if(verbose){
          cat("More than two edges")
        }
        edge.index <- which(unrooted.patient.subtree$edge.length==max.branch.length)
        
        start.node.index <- unrooted.patient.subtree$edge[edge.index,2]
        end.node.index <- unrooted.patient.subtree$edge[edge.index,1]
        
        split.set <- list()
        split.set[[1]] <- vector()
        split.set[[2]] <- vector()
        
        for(tip in seq(1,length(unrooted.patient.subtree$tip.label))){
          current.node <- tip
          while(!(current.node %in% c(start.node.index, end.node.index))){
            current.node <- Ancestors(unrooted.patient.subtree, current.node, type="parent")
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
        #this gets to split count right, but is a bit of a hack
        split.counts[[original.id]] <- split.counts[[original.id]]-1
      }
      
      for(tip.group in split.set){
        split.counts[[original.id]] <- split.counts[[original.id]]+1
        
        new.name <- paste(original.id,"_split",split.counts[[original.id]], sep="")
        if(verbose){
          cat("Child ",counter,": ", new.name,"\n", sep="")
        } 
        
        patients <- c(patients, new.name)
     
        tip.group.in.full.tree <- which(tree$tip.label %in% tip.group)
        
        patient.tips[[new.name]] <- tip.group.in.full.tree
        patient.mrcas[[new.name]] <- mrca.phylo.or.unique.tip(tree, tip.group.in.full.tree, zero.length.tips.count)
        
        split.patient(tree, new.name, split.threshold, patients, patient.tips, patient.mrcas, split.counts, original.id, already.split = TRUE)
      }
    } else if(verbose){
      cat("No split needed\n")
    }
  }
}

source("TransmissionUtilityFunctions.R")


option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Input tree file name", metavar="inputTreeFileName"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="Output file name [default= %default]", metavar="outputFileName"),
  make_option(c("-r", "--refSeqName"), type="character", default=NULL, 
              help="Reference sequence label (if unspecified, tree will be assumed to be already rooted)"),
  make_option(c("-p", "--patientIDPosition"), type="integer", default=1, 
              help="Position in the tip label of the patient ID (default=1)"),
  make_option(c("-s", "--labelSeparator"), type="character", default="_", 
              help="Field separator for tip labels (default underscore)"),
  make_option(c("-b", "--blacklist"), type="character", default=NULL, 
              help="Path to a .csv file listing tips to ignore"),
  make_option(c("-t", "--splitThreshold"), type="double", default=Inf, 
              help="Length threshold at which the root branch will be cut to make multiple infections"),
  make_option(c("-v", "--verbose"), type="logical", default=FALSE, 
              help="Talk about what I'm doing"),
  make_option(c("-z", "--zeroLengthTipsCount"), type="logical", default=FALSE, 
              help="If TRUE, a zero length terminal branch associates the parent at the tip with its parent
              node, interrupting any inferred transmission pathway between another pair of hosts that goes
              through the node.")
  
)



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


verbose <- opt$verbose
zero.length.tips.count <- opt$zeroLengthTipsCount

file.name <- opt$file
blacklist.file <- opt$blacklist
root.name <- opt$refSeqName
label.separator <- opt$labelSeparator
patient.id.position <- opt$patientIDPosition
split.threshold <- opt$fullSplitThreshold

cat("Opening file: ", file.name, sep = "")

tree <-read.tree(file.name)

#blacklist is a list of the tip numbers to avoid

if(!is.null(blacklist.file)){
  if(verbose){
    cat("Processing blacklist")
  }
  blacklist.table <- read.table(blacklist.file, sep=",", stringsAsFactors = FALSE)
  if(ncol(blacklist)!=1){
    stop("Blacklist should be a single column CSV file")
  }
  blacklist <- sapply(blacklist.table[,1], get.tip.no, tree=tree)
} else {
  blacklist <- NULL
}

tree <- root(tree, root.name)

if(verbose){
  cat("Making tree multiphyletic\n")
}

tree <- di2multi(tree, tol = 1E-5)

if(verbose){
  cat("Getting patient IDs\n")
}


tip.labels <- tree$tip.label

outgroup.no <- which(tip.labels == root.name)

split.tip.labels <- strsplit(tip.labels, label.separator)

patient.ids <- sapply(split.tip.labels, "[[", patient.id.position)

patient.ids[outgroup.no] <- "Ref"

#if all tips from a given patient are on the blacklist, safe to ignore him/her

patients <- unique(patient.ids[-blacklist])
patients <- patients[which(substr(patients, 1,3) == "BEE")]

patient.tips <-
  sapply(patients, get.tips.for.patient, tree = tree, patient.ids = patient.ids, blacklist=blacklist)

patient.mrcas <- lapply(patient.tips, function(node) mrca.phylo.or.unique.tip(tree, node, zero.length.tips.count))

split.counts <- list()

#splitting
if(split.threshold<Inf){
  for(id in patients){
    split.patient(tree, id, split.threshold, patients, patient.tips, patient.mrcas, split.counts)
  }
}

total.pairs <- length(patients) ^ 2 - length(patients)

if(verbose){
  cat("Testing pairs\n")
}

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
          is.direct.descendant.of(tree, desc, anc, patient.mrcas, patient.tips, zero.length.tips.count)
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

write.table(dddf, file = opt$out, sep = ",", row.names = FALSE, col.names = TRUE)
