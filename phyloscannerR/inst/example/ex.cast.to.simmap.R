#
# Example on data from Rakai Community Cohort Study
#
require(phyloscannerR)
ph.phyloscanner <- readRDS(file=system.file(file.path('extdata','ptyr292_1600_to_1849.rds'),package='phyloscannerR'))

#	cast from phyloscanner format to SIMMAP format
#	so that the SIMMAP tree still has the phyloscanner attributes
ph.simmap <- phyloscanner.to.simmap(ph.phyloscanner, delete.phyloscanner.structures=FALSE)
attr(ph, "SPLIT")
#	to delete the phyloscanner attributes, use
ph.simmap <- phyloscanner.to.simmap(ph.phyloscanner, delete.phyloscanner.structures=TRUE)
#
#	cast from SIMMAP format to phyloscanner format
#	similarly, use delete.simmap.structures to keep/delete SIMMAP elements from the tree in phyloscanner format
ph.phyloscanner2 <- simmap.to.phyloscanner(ph.simmap, delete.simmap.structures=TRUE)
#
#	test that ph.phyloscanner2 is the same as ph.phyloscanner:
stopifnot(all( ph.phyloscanner$tip.label == ph.phyloscanner2$tip.label ))
stopifnot(all( ph.phyloscanner$edge == ph.phyloscanner2$edge ))
stopifnot(all( ph.phyloscanner$edge.length == ph.phyloscanner2$edge.length ))
stopifnot(all( ph.phyloscanner$Nnode == ph.phyloscanner2$Nnode ))
tmp <- as.character(attr(ph.phyloscanner,'SPLIT'))
tmp2 <- as.character(attr(ph.phyloscanner2,'SPLIT'))
tmp[is.na(tmp)] <- 'Unknown'
tmp2[is.na(tmp2)] <- 'Unknown'
stopifnot( all(tmp==tmp2) )
tmp <- as.character(attr(ph.phyloscanner,'INDIVIDUAL'))
tmp2 <- as.character(attr(ph.phyloscanner2,'INDIVIDUAL'))
tmp[is.na(tmp)] <- 'Unknown'
tmp2[is.na(tmp2)] <- 'Unknown'
stopifnot( all(tmp==tmp2) )
tmp <- as.character(attr(ph.phyloscanner,'BRANCH_COLOURS'))
tmp2 <- as.character(attr(ph.phyloscanner2,'BRANCH_COLOURS'))
tmp[is.na(tmp)] <- 'Unknown'
tmp2[is.na(tmp2)] <- 'Unknown'
stopifnot( all(tmp==tmp2) )
tmp <- attr(ph.phyloscanner,'SUBGRAPH_MRCA')
tmp2 <- attr(ph.phyloscanner2,'SUBGRAPH_MRCA')
stopifnot( all(tmp==tmp2) )
	