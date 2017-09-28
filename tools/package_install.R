#!/usr/bin/Rscript

every.package <- c("argparse", 
                   "ape", 
                   "data.table",
                   "dplyr",
                   "dtplyr",
                   "ff",
                   "GGally",
                   "ggplot2", 
                   "grid", 
                   "gridExtra", 
                   "gtable", 
                   "kimisc",
                   "pegas",
                   "phangorn", 
                   "phytools", 
                   "prodlim",
                   "RColorBrewer", 
                   "reshape", 
                   "reshape2",
                   "scales")

new.packages <- every.package[!(every.package %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")

tryCatch({
  if(!("ggtree" %in% installed.packages()[,"Package"])){
    source("https://bioconductor.org/biocLite.R")
    biocLite("ggtree")
  }
}, error=function(err){
  message("Failed to install ggtree. Check your Bioconductor installation, in particular BiocInstaller, is up-to-date.")
})
