#!/usr/bin/env Rscript

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

require(phytools)

# Plot the tree.
# Make the automatic tip labels transparent, so we can colour later.
pdf(OutputTree)
MyTree <- read.tree(file=TreeFile)
NumTips <- length(MyTree$tip.label)
par(fg="transparent")
plotTree(MyTree, fsize=FontSize, lwd=LineWidth, ylim=c(-1, NumTips))

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

# Shuffle into a random order so that adjacent IDs don't have similar colours.
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

