library(ape)

args = commandArgs(trailingOnly=TRUE)

patient.from.label <- function(label){
  return(unlist(strsplit(label, "_read_"))[1])
}

# This script blacklists all reads from all patients in a given window. It is currently very 
# BEEHIVE specific (unlike MakeReadBlacklist.R) and it would be suggested that other users
# simply provide custom read blacklists.

start <- as.numeric(args[1])
end <- as.numeric(args[2])
in.file <- args[3]
tree.root <- args[4]
out.root <- args[5]

boundaries <- c(800, 1505, 2385, 4805, 5037, 5985, 7831, 9400)

amps <- read.table(in.file, sep=",", header=T, stringsAsFactors = F)
tree <- read.tree(paste(tree.root,"_",start,"_to_",end,".tree", sep=""))
out.file <- paste(out.root,"_",start,"_to_",end,".csv", sep="")

cat("MakePatientBlacklist.R run on: ", paste(tree.root,"_",start,"_to_",end,".tree", sep=""), "\n", sep="")

bad.patients <- vector()

for(row in seq(1,nrow(amps))){
  if(amps[row,2]){
    bad.patients <- c(bad.patients, amps[row,1])
  } else {
    for(amp.id in seq(1,7)){
      start.amp.section <- boundaries[amp.id]
      end.amp.section <- boundaries[amp.id+1]
      if(!amps[row,2+amp.id] & ((start>start.amp.section & start<end.amp.section) | (end>start.amp.section & end < end.amp.section))){
        bad.patients <- c(bad.patients, amps[row,1])
      }
    }
  }
}

tip.is.bad <- sapply(tree$tip.label, function(x)  patient.from.label(x) %in% bad.patients)

bad.tips <- tree$tip.label[which(tip.is.bad)]


cat("Blacklist length: ",length(bad.tips),"\n",sep="")

if(length(bad.tips)>0){
  write.table(bad.tips,out.file, sep=",", row.names=FALSE, col.names=FALSE, quote=F)
}