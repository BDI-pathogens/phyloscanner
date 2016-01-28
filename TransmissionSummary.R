# summarises transmission information

setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160111/")

patient.seq <- seq(1,907)
patient.ids <- paste("BEE", formatC(patient.seq, width = 4, format = "d", flag = "0"), "-1", sep="")

parent.list <- vector("list", length(patient.ids))

names(parent.list) <- patient.ids

window.count <- 0

for(start in seq(501, 6501, by=100)){
  end <- start+349 
  
  if(file.exists(paste("transmissions/transmissions.InWindow_",start,"_to_",end,".csv",sep=""))){
    
    window.count <- window.count + 1
    
    transmissions.table <- read.table(paste("transmissions/transmissions.InWindow_",start,"_to_",end,".csv",sep=""), 
                                      sep=",", 
                                      header=TRUE,
                                      stringsAsFactors=FALSE)
    
    
    for(row in seq(1,nrow(transmissions.table))){
      parent.list[[transmissions.table[row,1]]] <- c(parent.list[[transmissions.table[row,1]]], transmissions.table[row,2])
    }
    
  }
  
}

parents <- vector()
children <- vector()
presence <- vector()
presence.where.recorded <- vector()

parents.2 <- vector()
children.2 <- vector()
presence.2 <- vector()
presence.where.recorded.2 <- vector()


for(child in patient.ids){
  for(parent in unique(parent.list[[child]])){
    count <- sum(parent.list[[child]]==parent)
    
    parents <- c(parents, parent)
    children <- c(children, child)
    presence <- c(presence, count)
    presence.where.recorded <- c(presence.where.recorded, count/length(parent.list[[child]]))
    
    if(count>=2){
      
      parents.2 <- c(parents.2, parent)
      children.2 <- c(children.2, child)
      presence.2 <- c(presence.2, count)
      presence.where.recorded.2 <- c(presence.where.recorded.2, count/length(parent.list[[child]]))
      
    }
    
  }
}


out <- data.frame(parents, children, presence, presence.where.recorded)
out.2 <- data.frame(parents.2, children.2, presence.2, presence.where.recorded.2)

write.table(out, file="transmissionSummary.csv", sep=",", col.names = TRUE, row.names = FALSE)
write.table(out.2, file="transmissionSummary_noSingletons.csv", sep=",", col.names = TRUE, row.names = FALSE)