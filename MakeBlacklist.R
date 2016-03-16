args = commandArgs(trailingOnly=TRUE)

file.name <- args[1]
output.name <- args[2]
raw.threshold <- as.numeric(args[3])
ratio.threshold <- as.numeric(args[4])

table <- read.csv(file.name, stringsAsFactors = F)

cat("Making blacklist with a ratio threshold of ",args[4]," and a raw threshold of ",args[3],".\n",sep="")

new.col <- vector()

split.taxon.names <- strsplit(table[,1], "_")

count <- 0
for(item in split.taxon.names){
  count <- count+1
  
  new.col <- c(new.col, item[length(item)])
  
}

table$SeqCount <- as.numeric(new.col)

out <- table[which(table$SeqCount.SharedCount<ratio.threshold | (table$SeqCount.SharedCount<1 & table$SeqCount<raw.threshold)),1]
cat("Blacklist length: ",length(out),"\n",sep="")
if(length(out)>0){

  write.table(out,output.name, sep=",", row.names=FALSE, col.names=FALSE)
}