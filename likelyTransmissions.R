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

# Output the set of patients whose mrca is this node

get.patients.with.these.mrcas <- function(node, tree, patient.mrcas){
  out <- vector()
  mrca.vec <- sapply(patient.mrcas, "[[", 1)
  
  numbers <- which(mrca.vec==node)
  
  return(names(numbers))
}


tree <- read.tree("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160111/RAxML_bestTree/RAxML_bestTree.InWindow_4101_to_4450.tree")

root.name<-"Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455"

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

# Need a MRCA function which if given a single tip, returns that tip rather than NA

mrca.phylo.or.unique.tip <- function(x, node){
  if(length(node)==1){
    return(node)
  } else {
    return(mrca.phylo(x, node))
  }
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
  #situtation is too ambiguous
  
  anc.sequence <- get.ancestral.sequence(mrca.desc, mrca.anc, tree)
  
  patients.involved <- unlist(lapply(anc.sequence, get.patients.with.these.mrcas, patient.mrcas=patient.mrcas))
  
  if(length(patients.involved)<2){
    #something is seriously wrong
    stop("Sequence of ancestors contains less than two MRCAs")
  }
  
  if(length(patients.involved)>2){
    return(FALSE)
  }
  
  return(TRUE)
}

patient.tips <- sapply(patients, get.tips.for.patient, tree=tree, patient.ids = patient.ids)

patient.mrcas <- lapply(patient.tips, function(node) mrca.phylo.or.unique.tip(tree, node))

all.mrca.nodes <- unique(unlist(patient.mrcas))

total.pairs <- length(patients)^2 - length(patients)



count <- 0
descendant.matrix <- matrix(NA, length(patients), length(patients))
direct.descendant.matrix <- matrix(NA, length(patients), length(patients))
for(desc in seq(1, length(patients))){
  for(anc in seq(1, length(patients))){
    if(desc!=anc){
      count <- count + 1
      
      desc.patient <- patients[desc]
      anc.patient <- patients[anc]
      
      descendant.matrix[desc, anc] <- is.descendant.of(desc.patient, anc.patient, tree, patient.mrcas, patient.tips)
      direct.descendant.matrix[desc, anc] <- is.direct.descendant.of(desc.patient, anc.patient, tree, patient.mrcas, patient.tips, all.mrca.nodes)
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

dddf <- dddf[dddf$Present,]

write.table(dddf, file="/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160111/transmissions.InWindow_4101_to_4450.csv", sep=",", row.names=FALSE, col.names=TRUE)