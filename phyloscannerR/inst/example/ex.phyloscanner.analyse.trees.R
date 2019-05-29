#
# Example on data from Rakai Community Cohort Study
#
\dontrun{

require(phyloscannerR)

#	extract RCCS example data
tree.file.zip <- system.file(file.path('extdata','Rakai_run192_trees.zip'),package='phyloscannerR')
tree.file.directory <- tempdir()	
unzip(tree.file.zip, exdir=tree.file.directory, junkpaths=TRUE)
	
#	arguments used for RCCS analysis
file.name.regex <- "^\\D*([0-9]+)_to_([0-9]+)\\D*$"
max.reads.per.host <- 50
multifurcation.threshold <- 1e-5
norm.ref.file.name <- system.file('HIV_DistanceNormalisationOverGenome.csv',package='phyloscannerR')	
outgroup.name <- "REF_CPX_AF460972"
raw.blacklist.threshold <- 20
sankoff.k <- 20
sankoff.unassigned.switch.threshold <- 0
seed <- 42
splits.rule <- 's'
relaxed.ancestry <- TRUE
allow.mt <- TRUE
tip.regex <- "^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$"
tree.file.regex <- "^ptyr192_InWindow_([0-9]+_to_[0-9]+)\\.tree$"
verbosity <- 1

#	analyse deep sequence trees
phsc <- phyloscanner.analyse.trees(tree.file.directory,
             allow.mt=allow.mt,
             alignment.file.directory = NULL, 
             alignment.file.regex = NULL,
             blacklist.underrepresented = FALSE,
             count.reads.in.parsimony = TRUE,
             do.dual.blacklisting = FALSE,
             duplicate.file.directory = NULL,
             duplicate.file.regex = NULL,
             file.name.regex = file.name.regex,
             guess.multifurcation.threshold = FALSE,
             max.reads.per.host = max.reads.per.host,
             multifurcation.threshold = multifurcation.threshold,
             norm.constants = NULL,
             norm.ref.file.name = NULL,
             norm.standardise.gag.pol = TRUE,
             no.progress.bars = FALSE,
             outgroup.name = outgroup.name,
             parsimony.blacklist.k = sankoff.k,
             prune.blacklist = FALSE,
             ratio.blacklist.threshold = 0, 
             raw.blacklist.threshold = raw.blacklist.threshold,			
             recombination.file.directory = NULL,
             recombination.file.regex = NULL,
             relaxed.ancestry = relaxed.ancestry,
             sankoff.k = sankoff.k,
             sankoff.unassigned.switch.threshold = sankoff.unassigned.switch.threshold,
             seed = seed,
             splits.rule = splits.rule,
             tip.regex = tip.regex,
             tree.file.regex = tree.file.regex,
             use.ff = FALSE,
             user.blacklist.directory = NULL, 
             user.blacklist.file.regex = NULL,
             verbosity = verbosity 
             )
}