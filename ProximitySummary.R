# summarises transmission information

setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/Matthew_notebook/")


get.neighbours <- function(table, patient, threshold){
  part.1 <- table[which(table$pat.1==patient & table$mrca.distance<=threshold ),2]
  part.2 <- table[which(table$pat.2==patient & table$mrca.distance<=threshold ),1]
  return(c(part.1, part.2))
}

get.nearest.neighbours <- function(table, patient, number){
  part.1 <- table[which(table$pat.1==patient),c(2,3)]
  part.2 <- table[which(table$pat.2==patient),c(1,3)]
  colnames(part.1) <- c("Neighbour", "Distance")
  colnames(part.2) <- c("Neighbour", "Distance")
  full <- rbind(part.1, part.2)
  full <- full[order(full$Distance),]
  return(full[1:number,])
}

contains.split <- function(string){
  return(grepl("split", string))
}

splits.list <- list()

for(start in seq(800, 9200, by=200)){
  print(start)
  end <- start+300 
  
  if(file.exists(paste("ProximityCheck_",start,"_to_",end,".csv",sep=""))){
    
    
    proximity.table <- read.table(paste("ProximityCheck_",start,"_to_",end,".csv",sep=""), 
                                  sep=",", 
                                  header=TRUE,
                                  stringsAsFactors=FALSE)
    
    proximity.table <- proximity.table[which(proximity.table$pat.1!=proximity.table$pat.2),]
    
    
    split.patients <- unique(c(proximity.table[which(contains.split(proximity.table$pat.1)),1],
                               proximity.table[which(contains.split(proximity.table$pat.2)),2]))
    
    split.originals <- substring(split.patients, 1, 9) 
    
    for(patient in unique(split.originals)){
      if(!is.null(splits.list[[patient]])){
        cat(patient," already exists in the list and has value ",splits.list[[patient]],".\n",sep="")
        splits.list[patient] <- splits.list[[patient]] + 1
      } else {
        cat(patient," is new to the list.\n", sep="")
        splits.list[patient] <- 1
      }
    }
  }
}

out <- vector()

for(start in seq(800, 9200, by=200)){
  
  end <- start+300 
  
  cat("Window: ",start," to ",end,".\n", sep="")
  
  if(file.exists(paste("ProximityCheck_",start,"_to_",end,".csv",sep=""))){
    
    
    proximity.table <- read.table(paste("ProximityCheck_",start,"_to_",end,".csv",sep=""), 
                                  sep=",", 
                                  header=TRUE,
                                  stringsAsFactors=FALSE)
    
    proximity.table <- proximity.table[which(proximity.table$pat.1!=proximity.table$pat.2),]
    
    for(patient in names(splits.list)){
      cat("Examining patient ", patient, ".\n",sep="")
      
      if(nrow(proximity.table[which(proximity.table$pat.1==patient),])>0 |
         nrow(proximity.table[which(proximity.table$pat.2==patient),])>0){
        cat("No split in this window.\n")
        
        neighbours <- get.nearest.neighbours(proximity.table, patient, 5)
        
        row <- c(start, end, patient, patient, neighbours[,1])
        out <- rbind(out, row)
      } else {
        cat("Split in this window:")
        splits.exist.here <- TRUE
        counter <- 1
        while(splits.exist.here){
          split.name <- paste(patient, "_split", counter, sep="")
          
          if(nrow(proximity.table[which(proximity.table$pat.1==split.name),])>0 |
             nrow(proximity.table[which(proximity.table$pat.2==split.name),])>0){ 
            cat(split.name,"...")
            
            neighbours <- get.nearest.neighbours(proximity.table, split.name, 5)
            
            row <- c(start, end, patient, split.name, neighbours[,1])
            out <- rbind(out, row)
            counter <- counter + 1
            
          } else {
            splits.exist.here <- FALSE
            cat("done.\n")
          }
        }
      }
      
    }
  }
}

out.df <- data.frame(out)

colnames(out.df) <- c("start", "end", "patient", "split.patient", "neighbour.1", "neighbour.2", 
                      "neighbour.3", "neighbour.4", "neighbour.5")

write.table(out.df, file="neighbours.csv", sep="", row.names=NA)