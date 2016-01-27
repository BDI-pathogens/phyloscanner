library(phangorn)

verbose <- FALSE
# Get the tip number of the named taxon



setwd("/Users/twoseventwo/Documents/phylotypes/")
source("transmissionUtilityFunctions.R")

setwd(
  "/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160111/"
)

for (start in seq(6701, 9301, by = 100)) {
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
    
    if(verbose){
      cat("Making tree multiphyletic\n")
    }
    
    tree <- di2multi(tree, tol = 1E-5)
    
    if(verbose){
      cat("Getting patient IDs\n")
    }
    
    
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
        "transmissions/transmissions.InWindow_",start,"_to_",end,".csv", sep = ""
      ), sep = ",", row.names = FALSE, col.names = TRUE
    )
  }
}