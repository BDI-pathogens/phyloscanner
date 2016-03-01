# summarises transmission information

setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/JuicyCladeOrig_NoMerge/")

patient.seq <- seq(1,907)
patient.ids <- paste("BEE", formatC(patient.seq, width = 4, format = "d", flag = "0"), "-1", sep="")

parent.list <- vector("list", length(patient.ids))

names(parent.list) <- patient.ids

window.count <- 0

for(start in seq(800, 9200, by=200)){
  print(start)
  end <- start+300 
  
  if(file.exists(paste("LikelyTransmissions_",start,"_to_",end,".csv",sep=""))){
    
    window.count <- window.count + 1
    
    transmissions.table <- read.table(paste("LikelyTransmissions_",start,"_to_",end,".csv",sep=""), 
                                      sep=",", 
                                      header=TRUE,
                                      stringsAsFactors=FALSE)
    
    if(nrow(transmissions.table)>0){
      for(row in seq(1,nrow(transmissions.table))){
        parent <- transmissions.table[row,1]
        child <- transmissions.table[row,2]
        
        parent <- substr(parent, 1, 9)
        child <- substr(child, 1, 9)
        
        parent.list[[parent]] <- c(parent.list[[parent]], child)
      }
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

parents.3 <- vector()
children.3 <- vector()
presence.3 <- vector()
presence.where.recorded.3 <- vector()

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
    
    if(count>=3){
      
      parents.3 <- c(parents.3, parent)
      children.3 <- c(children.3, child)
      presence.3 <- c(presence.3, count)
      presence.where.recorded.3 <- c(presence.where.recorded.3, count/length(parent.list[[child]]))
      
    }
    
  }
}


out <- data.frame(parents, children, presence, presence.where.recorded)
out.2 <- data.frame(parents.2, children.2, presence.2, presence.where.recorded.2)
out.3 <- data.frame(parents.3, children.3, presence.3, presence.where.recorded.3)

colnames(out.2) <- colnames(out)
colnames(out.3) <- colnames(out)

write.table(out, file="transmissionSummary.csv", sep=",", col.names = TRUE, row.names = FALSE)
write.table(out.2, file="transmissionSummary_noSingletons.csv", sep=",", col.names = TRUE, row.names = FALSE)
write.table(out.3, file="transmissionSummary_min3.csv", sep=",", col.names = TRUE, row.names = FALSE)