\dontrun{
require(data.table)
require(tidyverse)
require(phyloscannerR)

#
# reconstruct transmission networks and most likely transmission chain
#
infile <- system.file(file.path('extdata','Rakai_phscnetworks_allpairs_190706.RData'),package='phyloscannerR')
load(infile) 	# loads phyloscanner relationship counts 'dc'
tmp <- find.networks(dc, neff.cut=3, verbose=TRUE)
dnet <- copy(tmp$transmission.networks)
dchain <- copy(tmp$most.likely.transmission.chains)

# get meta data
meta.file <- file.path(workdir,"RakaiPopSample_data","Dataset_S2.csv")
tmp <- "https://datadryad.org/bitstream/handle/10255/dryad.208474/Dataset_S2.csv?sequence=2"
download.file(tmp, destfile=meta.file, method="curl")
dmeta <- as_tibble(read.csv(meta.file, stringsAsFactors=FALSE))

#
# plot with host names, and male/female in blue/pink nodes, 
# 	and highlight pairings with >60% support for phylogenetic linkage
#
idclus <- sort(unique(dnet$IDCLU))
di <- copy(dmeta)
setnames(di, 'ID', 'H')
df <- dnet %>% 
		filter(IDCLU == 34) %>%
		select(-c(H1_SEX,H2_SEX))		
control<- list()
control$point.size = 10
control$edge.gap = 0.04
control$edge.size = 2
control$curvature = -0.2
control$arrow = arrow(length = unit(0.04, "npc"), type = "open")
control$curv.shift = 0.06
control$label.size = 3
control$node.label = "H" 
control$node.fill = "SEX"
control$node.shape = NA_character_
control$node.shape.values = c(`NA` = 16)
control$node.fill.values = c(F = "hotpink2", M = "steelblue2")
control$threshold.linked = 0.6
p <- plot.network(df, di, control)	

}