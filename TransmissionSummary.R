# summarises transmission information

library(prodlim)

setwd("/Users/twoseventwo/Dropbox (Infectious Disease)/BEEHIVE/phylotypes/run20160307/")

patient.seq <- seq(1,907)
patient.ids <- paste("BEE", formatC(patient.seq, width = 4, format = "d", flag = "0"), "-1", sep="")

parent.list <- vector("list", length(patient.ids))

names(parent.list) <- patient.ids

for(patient.id in patient.ids){
  parent.list[[patient.id]] <- rep(NA, 31)
}

window.count <- 0

transmissions.ever.suggested <- vector()

for(start in seq(800, 9050, by=275)){
  window.count <- window.count + 1
  
  end <- start+375 
  
  cat("Window ",start, " to ",end,"\n",sep="")

  
  if(file.exists(paste("LikelyTransmissions_",start,"_to_",end,".csv",sep=""))){
    
    transmissions.table <- read.table(paste("LikelyTransmissions_",start,"_to_",end,".csv",sep=""), 
                                      sep=",", 
                                      header=TRUE,
                                      stringsAsFactors=FALSE)
    
    if(nrow(transmissions.table)>0){
      for(row in seq(1,nrow(transmissions.table))){
        if(transmissions.table$Relationship[row]=="anc" | transmissions.table$Relationship[row]=="desc"){
          
          pat.1 <- transmissions.table[row,1]
          pat.2 <- transmissions.table[row,2]
          
          pat.1 <- substr(pat.1, 1, 9)
          pat.2 <- substr(pat.2, 1, 9)
          
          transmissions.ever.suggested <- rbind(transmissions.ever.suggested, c(pat.1, pat.2))

        } 
      }
    }
  }
}

transmissions.ever.suggested <- as.data.frame(unique(transmissions.ever.suggested), stringsAsFactors=FALSE)

reversed.t.e.s <- as.data.frame(cbind(transmissions.ever.suggested[,2],transmissions.ever.suggested[,1]), stringsAsFactors = FALSE)

matches <- row.match(transmissions.ever.suggested, reversed.t.e.s)

to.go <- rep(FALSE, nrow(transmissions.ever.suggested))

for(comp in seq(1, length(matches))){
  if(!(is.na(matches[comp])) & (comp < matches[comp])){
    to.go[comp] <- TRUE
  }
}

transmissions.ever.suggested <- transmissions.ever.suggested[which(!to.go), ]
window.count <- 0

first <-TRUE

for(start in seq(800, 9050, by=275)){
  window.vector <- rep(NA, nrow(transmissions.ever.suggested))
  
  window.count <- window.count + 1
  
  print(start)
  end <- start+375 
  
  if(file.exists(paste("LikelyTransmissions_",start,"_to_",end,".csv",sep=""))){
    
    
    
    transmissions.table <- read.table(paste("LikelyTransmissions_",start,"_to_",end,".csv",sep=""), 
                                      sep=",", 
                                      header=TRUE,
                                      stringsAsFactors=FALSE)
    
    transmissions.table[,1] <- substr(transmissions.table[,1],1,9)
    transmissions.table[,2] <- substr(transmissions.table[,2],1,9)
    
    if(nrow(transmissions.table)>0){
      for(row in seq(1,nrow(transmissions.table))){
        # if the row exists
        if(!is.na(row.match(transmissions.table[row,1:2], transmissions.ever.suggested))){
          #whatever the third column says, things are the right way round
          
          window.vector[row.match(transmissions.table[row,1:2], transmissions.ever.suggested)] <- transmissions.table[row,3] 
        } else if(!is.na(row.match(rev(transmissions.table[row,1:2]), transmissions.ever.suggested))){
          #transmissions will be the wrong way round
          
          if(transmissions.table[row,3] == "anc"){
            window.vector[row.match(rev(transmissions.table[row,1:2]), transmissions.ever.suggested)] <- "desc" 
          } else if (transmissions.table[row,3] == "desc"){
            window.vector[row.match(rev(transmissions.table[row,1:2]), transmissions.ever.suggested)] <- "anc" 
          } else {
            window.vector[row.match(rev(transmissions.table[row,1:2]), transmissions.ever.suggested)] <- transmissions.table[row,3]
          }
        }
      }
    }
  }
  if(first){
    first <- FALSE
    out <- window.vector
  } else {
    out <- cbind(out, window.vector) 
  }
}

out <- cbind(transmissions.ever.suggested, out)

colnames(out) <- c("pat.1", "pat.2", paste("window.start.",seq(800,9050,by=275), sep="")  ) 
  
write.table(out, file="quickout.csv", sep=",", col.names=NA)

first <- TRUE

for(row in seq(1, nrow(out))){
  row.list <- list()
  row.list["anc"] <- 0
  row.list["desc"] <- 0
  row.list["int"] <- 0
  row.list["sib"] <- 0
  
  for(col in seq(3, ncol(out))){
    if(!is.na(out[row,col])){
      row.list[[out[row,col]]] <- row.list[[out[row,col]]] + 1
    }
  }
  
  anc.row <- c(out[row,1], out[row,2], row.list[["anc"]], "trans", row.list[["anc"]]+row.list[["desc"]])
  desc.row <- c(out[row,2], out[row,1], row.list[["desc"]], "trans", row.list[["anc"]]+row.list[["desc"]])
  int.row <- c(out[row,1], out[row,2], row.list[["int"]], "int", row.list[["anc"]]+row.list[["desc"]])
  sib.row <- c(out[row,1], out[row,2], row.list[["sib"]], "sib", row.list[["anc"]]+row.list[["desc"]])
  
  new.rows <- data.frame(rbind(anc.row, desc.row, int.row, sib.row), stringsAsFactors = F)
  
  colnames(new.rows) <- c("pat.1", "pat.2", "windows", "type", "total.trans")
  
  if(first){
    first <- FALSE
    new.out <- new.rows[which(new.rows$windows>0),]
  } else {
    new.out <- rbind(new.out, new.rows[which(new.rows$windows>0),])
  }
  
}

new.out$windows <- as.numeric(new.out$windows)
new.out$total.trans <- as.numeric(new.out$total.trans)

write.table(new.out, file="transSummary.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out[which(new.out$total.trans>=2),], file="transSummary_2.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out[which(new.out$total.trans>=3),], file="transSummary_3.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out[which(new.out$total.trans>=4),], file="transSummary_4.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out[which(new.out$total.trans>=5),], file="transSummary_5.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out[which(new.out$total.trans>=6),], file="transSummary_6.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out[which(new.out$total.trans>=10),], file="transSummary_10.csv", row.names = F, sep=",", quote=FALSE)
write.table(new.out[which(new.out$total.trans>=16),], file="transSummary_16.csv", row.names = F, sep=",", quote=FALSE)
