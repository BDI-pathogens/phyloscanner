#!/usr/bin/env Rscript

command.line<-FALSE
if(command.line) {
  # Read in the arguments
  args <- commandArgs(TRUE)
  if (length(args) != 7) {
    cat(paste("Exactly 7 arguments must be specified:\na tree,\nthe name of",
              "the output pdf file that will be produced,\nthe string/character with which",
              "all tip-names-to-be-coloured begin (other tip names are left black),\nthe",
              "string/character that marks the end of the ID part of the tip name (i.e. if",
              "tip names have a common name up until that point, they get coloured the",
              "same),\nthe font size for the tip labels,\nthe line width for the tree,\n",
              "and a seed for randomising the colour shuffling between different",
              "patients.\nQuitting.\n"))
    quit("no", 1)
  }
  TreeFile   <- args[1]
  OutputTree <- args[2]
  TipNamePrefix <- args[3]
  IDdelimiter <- args[4]
  FontSize <- as.numeric(args[5])
  LineWidth <- as.numeric(args[6])
  seed <- as.integer(args[7])
  if (! file.exists(TreeFile)) {
    print(paste(TreeFile, "does not exist. Quitting."))
    quit("no", 1)
  }
} else{
  # Some options to run within R
  setwd("/Users/cfraser/Dropbox (Infectious Disease)/PROJECTS/BEEHIVE/phylotypes")
  rm(list=ls())
  TreeFiles   <- "TreeFiles.txt"
  InputTreeFolder <- "PhylotypeOutput/"
  OutputTreeFolder <- "CF_TREES/"
  TipNamePrefix <- "BEE"
  IDdelimiter <- "-1"
  FontSize <- 0.3
  LineWidth <- 1
  seed <- 12345
}

require(phytools)
require(grid)

TreeFiles<-scan(TreeFiles,what="",sep="\n")
dir.create(OutputTreeFolder,showWarnings=F)
#StatsList<-list()

#dummy<-rep(0,length(UniqueIDs))
#PatStats<-data.frame(PatientID=UniqueIDs,N.Leaves=dummy,Monophyletic=dummy,Mean.Size=dummy,CoV.Size=dummy)
#PatStats<-data.frame(PatientID=character(),Window=integer(),N.Leaves=integer(),Monophyletic=integer(),Mean.Size=double(),CoV.Size=double())
PatStats<-data.frame(PatientID="test",Window=0,N.Leaves=1,Monophyletic=1,Mean.Size=1,CoV.Size=1)
PatStats$PatientID<-as.character(PatStats$PatientID)

for(treefile in TreeFiles){
  #  treefile<-TreeFiles[[1]]
  OutputTreeFile<-paste(OutputTreeFolder,treefile,".pdf",sep="")
  Window<-as.numeric(strsplit(treefile,"_")[[1]][3]) #This is exceptionally bad code - how else to pull out window coordinate though?
  
  
  # Plot the tree.
  # Make the automatic tip labels transparent, so we can colour later.
  pdf(OutputTreeFile)
  opar<-par()
  MyTree <- read.tree(file=paste(InputTreeFolder,treefile,sep=""))
  NumTips <- length(MyTree$tip.label)
  par(fg="transparent")
  plotTree(MyTree, fsize=FontSize, lwd=LineWidth, ylim=c(-1, NumTips))
  
  # Find the longest tip label, and append "_" to it to create a string guaranteed
  # not to clash with any of the tip labels.
  LongestTipLabel <- ""
  for (TipLabel in MyTree$tip.label) {
    if (nchar(TipLabel) > nchar(LongestTipLabel)) {
      LongestTipLabel <- TipLabel
    }
  }
  NewTipLabel <- paste(LongestTipLabel, "_", sep="")
  
  # Find all unique IDs (the bit of the tip label coming before the delimiter, for
  # those tip labels beginning with the desired prefix).
  TipNameRegex <- paste("^", TipNamePrefix, sep="")
  UniqueIDs <- vector()
  NumUniqueIDs <- 0
  for (TipLabel in MyTree$tip.label) {
    if (grepl(TipNameRegex, TipLabel)) {
      ID <- unlist(strsplit(TipLabel, IDdelimiter))[1]
      if (! ID %in% UniqueIDs) {
        NumUniqueIDs <- NumUniqueIDs + 1
        UniqueIDs[NumUniqueIDs] <- ID
      }
    }
  }
  
  # Shuffle into a random order so that adjecent IDs don't have similar colours.
  set.seed(seed)
  UniqueIDs <- sample(UniqueIDs)
  
  # Define a colour for each unique ID.
  # Needs to be written twice to work properly for some reason!
  colors<-setNames(palette(rainbow(length(UniqueIDs)))[1:length(UniqueIDs)],
                   UniqueIDs)
  colors<-setNames(palette(rainbow(length(UniqueIDs)))[1:length(UniqueIDs)],
                   UniqueIDs)
  
  # Assign each tip its colour
  colors2 <- vector()
  for (i in 1:NumTips) {
    TipLabel <- MyTree$tip.label[i]
    if (grepl(TipNameRegex, TipLabel)) {
      ID <- unlist(strsplit(TipLabel, IDdelimiter))[1]
      colors2[i] <- colors[which(UniqueIDs==ID)]
    }
    else {
      colors2[i] <- "black"
    }
  }
  
  # Add the coloured tip labels
  lastPP <- get("last_plot.phylo", env=.PlotPhyloEnv)
  par(fg="black")
  text(lastPP$xx[1:NumTips], lastPP$yy[1:NumTips], MyTree$tip.label, cex=FontSize,
       col=colors2, pos=4, offset=0.1, font=0.1)
  
  dev.off()
  par<-opar
  
  
  patient.tips<-list()
  for(id in UniqueIDs) patient.tips[[id]]<-vector()
  
  for(TipLabel in MyTree$tip.label) {
    if (grepl(TipNameRegex, TipLabel)) {
      IDtip <- unlist(strsplit(TipLabel, IDdelimiter))[1]
      patient.tips[[IDtip]]<-c(patient.tips[[IDtip]],TipLabel)
    }
  }
  
  for(i in 1:length(UniqueIDs)) {
    ID<-UniqueIDs[i]
    N.Leaves<-length(patient.tips[[ID]])
    Monophyletic<-as.numeric(is.monophyletic(MyTree,patient.tips[[ID]]))
    if(N.Leaves>1) {
      SubTree<-drop.tip(phy = MyTree,tip = MyTree$tip.label[!(MyTree$tip.label %in% patient.tips[[ID]])])
      subtree.dist<-as.vector(cophenetic(SubTree))
      # Mean and variance of pairwise patristic distances between tips(a.k.a. cophenetic distances) 
      # These are not weighted for the number of reads for each tip 
      # Corrections are included because mean and variance of distance matrices include zeros on diagonal, but shouldn't
      Mean.Size<-mean(subtree.dist)/(1-1/N.Leaves)
      Variance<-var(subtree.dist)*(1+1/N.Leaves)-Mean.Size^2/N.Leaves
      CoV.Size<- ifelse(N.Leaves>2,sqrt(Variance)/Mean.Size,0)
    } else {
      Mean.Size<-0
      CoV.Size<-0
    }
    PatStats<-rbind(PatStats,c(ID,Window,N.Leaves,Monophyletic,Mean.Size,CoV.Size))
  }
}

PatStats$Mean.Size<-as.numeric(PatStats$Mean.Size)
PatStats$CoV.Size<-as.numeric(PatStats$CoV.Size)
PatStats<-PatStats[!PatStats$PatientID=="test",]
PatStats<-PatStats[order(PatStats$PatientID),]

pdf("Patient.Statistics.pdf")
for (ID in UniqueIDs){
  par(mfrow=c(2,3))
  pat<-PatStats[PatStats$PatientID==ID,]
  plot(pat$Window,pat$N.Leaves,main = pat$PatientID[1],xlab = "Genome location",ylab="Number of tips")
  plot(pat$Window,pat$Monophyletic,main = pat$PatientID[1],xlab = "Genome location",ylab="Monophyletic",ylim = c(0,1))
  plot(pat$Window,pat$Mean.Size,main = pat$PatientID[1],xlab = "Genome location",ylab="Mean Cophenetic Distance",ylim=c(0,max(PatStats$Mean.Size)))
  plot(pat$Window,pat$CoV.Size,main = pat$PatientID[1],xlab = "Genome location",ylab="CoV of Cophenetic Distance",ylim=c(0,max(PatStats$CoV.Size)))
}
dev.off()
