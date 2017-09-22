#!/usr/bin/env Rscript

# Load required packages
list.of.packages <- c("ape", "argparse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat("Missing dependencies; replacing the path below appropriately, run\n[path to your phyloscanner code]/tools/package_install.R\nthen try again.\n")
  quit(save="no", status=1)
}
suppressMessages(library(ape, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(argparse, quietly=TRUE, warn.conflicts=FALSE))

# Set up the arguments
arg.parser = ArgumentParser(description=paste("Calculate the median of the",
"patristic distances between all possible pairs of tips in a tree."))
arg.parser$add_argument("tree file")
args <- arg.parser$parse_args()

tree <- read.tree(file=args[['tree file']])

pdm <- cophenetic(tree) # symmetric matrix of patristic distances

pat.dists.list <- pdm[row(pdm) > col(pdm)] # take one set of off-diagonals

cat(median(pat.dists.list))
