args = commandArgs(trailingOnly=TRUE)

# This script blacklists all reads from all patients in a given window. It is currently very 
# BEEHIVE specific (unlike MakeReadBlacklist.R) and it would be suggested that other uses
# simply provide custom read blacklists

start <- as.numeric(args[1])
end <- as.numeric(args[2])
in.root <- args[3]
out.root <- args[4]

boundaries <- c(800, 1505, 2385, 4805, 5037, 5985, 7831, 9400)

amps <- read.table(in.root, sep=",", header=T, stringsAsFactors = F)

for(row in seq(1,nrow(amps))){
  if(amps[row,2]){
    out <- rbind(amps[row,1])
  } else {
    for(amp.id in seq(1,7)){
      start.amp.section <- boundaries[amp.id]
      end.amp.section <- boundaries[amp.id+1]
      if(!amps[row,2+amp.id] & ((start>start.amp.section & start<end.amp.section) | (end>start.amp.section & end < end.amp.section))){
        out <- rbind(out, amps[row,1])
      }
    }
  }
}

cat("Blacklist length: ",length(out),"\n",sep="")

if(length(out)>0){
  write.table(out,output.name, sep=",", row.names=FALSE, col.names=FALSE)
}
}