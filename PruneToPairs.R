# For identified transmission pairs, prunes the relevant phylogeny down to just the tips involving the hosts in question
# (and the outgroup)

library(phangorn)

setwd("/Users/twoseventwo/Documents/phylotypes/")
source("transmissionUtilityFunctions.R")

setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160111/")

for(start in seq(501, 9301, by=100)){
  end <- start+349 
  if(file.exists(paste("RAxML_bestTree/RAxML_bestTree.InWindow_",start,"_to_",end,".tree",sep =""))){
    
    cat("Opening file: RAxML_bestTree.InWindow_",start,"_to_",end,".tree\n", sep = "")
    
    tree <- read.tree(paste("RAxML_bestTree/RAxML_bestTree.InWindow_",start,"_to_",end,".tree", sep =""))
    
    root.name <- "Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
    
    tree <- root(tree, root.name)
    
    tree <- di2multi(tree, tol = 1E-5)
    
    tip.labels <- tree$tip.label
    
    outgroup.no <- which(substr(tip.labels, 1, 3) == "Ref")
    
    tip.labels[outgroup.no] <- "Ref_read_NA_count_NA"
    
    split.tip.labels <- strsplit(tip.labels, "_")
    
    patient.ids <- sapply(split.tip.labels, "[[", 1)
    
    transmissions.table <- read.table(paste("transmissions/transmissions.InWindow_",start,"_to_",end,".csv",sep=""), 
                                      sep=",", 
                                      header=TRUE,
                                      stringsAsFactors=FALSE)
    
    all.identified <- unique(c(transmissions.table$Descendant, transmissions.table$Ancestor))
    
    tips <- sapply(all.identified, get.tips.for.patient, tree=tree, patient.ids=patient.ids)
    
    to.keep <- unique(unlist(tips))
    
    #keep the outgroup
    
    to.keep <- c(to.keep, outgroup.no)
  
    pruned.tree <- drop.tip(tree, tree$tip.label[-to.keep])
    
    pruned.tree <- di2multi(pruned.tree, tol = 1E-5)
    
    write.nexus(pruned.tree, file=paste("transmissionPairSubtrees/prunedTree.InWindow_",start,"_to_",end,".allIdentified.nex",sep=""))
    
    for(line in seq(1:nrow(transmissions.table))){
      tips <- sapply(transmissions.table[line,1:2], get.tips.for.patient, tree=tree, patient.ids=patient.ids)
      
      to.keep <- unique(unlist(tips))
      
      #keep the outgroup
      
      to.keep <- c(to.keep, outgroup.no)
      
      pruned.tree <- drop.tip(tree, tree$tip.label[-to.keep])
      
      pruned.tree <- di2multi(pruned.tree, tol = 1E-5)
      
      parent <- transmissions.table[line,2]
      child <- transmissions.table[line,1]
      
      write.nexus(pruned.tree, file=paste("transmissionPairSubtrees/prunedTree.InWindow_",start,"_to_",end,".",parent,"_infects_",child,".nex",sep=""))
      
    }
    
    
    
  }
}
