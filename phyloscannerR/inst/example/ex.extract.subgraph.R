#
# Example on data from Rakai Community Cohort Study
#
require(phyloscannerR)
ph.phyloscanner <- readRDS(file=system.file(file.path('extdata','ptyr292_1600_to_1849.rds'),package='phyloscannerR'))

# add SIMMAP elements to tree
# because internally the extract.subgraph function uses functions written for
# trees in SIMMAP format
ph.phsc.plus.simmap <- phyloscanner.to.simmap(ph.phyloscanner, delete.phyloscanner.structures=FALSE)

# extract subgraph MRCAS and choose all those corresponding to a particular host
host <- 'RkA07714M'
mrcas <- which( attr(ph.phsc.plus.simmap, 'SUBGRAPH_MRCA') )
mrcas <- mrcas[ attr(ph.phsc.plus.simmap, 'INDIVIDUAL')[mrcas]==host ]	

# extract subgraphs
subgraphs <- lapply(mrcas, function(mrca) extract.subgraph(ph.phsc.plus.simmap, mrca))
# the first subgraph consists of just 1 taxon and is therefore hard to represent as an ape object
# note the dummy structures
# the second subgraph consists of 9 taxa and is in standard ape format
# note the extra elements "subgraph.name", "subgraph.root.edge", "subgraph.parent.state"
str( subgraphs[[2]] )