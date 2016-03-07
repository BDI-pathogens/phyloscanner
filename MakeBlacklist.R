args = commandArgs(trailingOnly=TRUE)

file.name <- args[1]
output.name <- args[2]
threshold <- args[3]

table <- read.csv(file.name, stringsAsFactors = F)

out <- table[which(table[,3]<threshold),1]
if(length(out)>0){

  write.table(table[which(table[,3]<threshold),1],output.name, sep=",", row.names=FALSE, col.names=FALSE)
}