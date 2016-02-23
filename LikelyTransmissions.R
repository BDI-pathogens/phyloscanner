library(phangorn)
library(optparse)
library(phytools)

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
output.name <- opt$out
blacklist.file <- opt$blacklist
root.name <- opt$refSeqName
label.separator <- opt$labelSeparator
patient.id.position <- opt$patientIDPosition
split.threshold <- opt$splitThreshold

cat("Opening file: ", file.name, "\n", sep = "")

tree <-read.tree(file.name)

#blacklist is a list of the tip numbers to avoid

if(!is.null(blacklist.file)){
  if(verbose){
    cat("Processing blacklist")
  }
  blacklist.table <- read.table(blacklist.file, sep=",", stringsAsFactors = FALSE)
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

if(!is.null(blacklist)){
  patients <- unique(patient.ids[-blacklist])
} else {
  patients <- unique(patient.ids)
}

patients <- patients[which(substr(patients, 1,3) == "BEE")]

patient.tips <-
  sapply(patients, get.tips.for.patient, tree = tree, patient.ids = patient.ids, blacklist=blacklist)

patient.mrcas <- lapply(patient.tips, function(node) mrca.phylo.or.unique.tip(tree, node, zero.length.tips.count))

split.counts <- list()

#splitting
if(split.threshold<Inf){
  for(id in patients){
    split.patient(tree, id, split.threshold, patients, patient.tips, patient.mrcas, patient.ids, split.counts)
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

write.table(dddf, file = output.name, sep = ",", row.names = FALSE, col.names = TRUE)
