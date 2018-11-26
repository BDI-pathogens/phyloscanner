#' @title Select for further analysis relationship classifications by read and tip counts 
#' @export
#' @import tidyverse 
#' @author Oliver Ratmann, Matthew Hall
#' @param ptrees A list of class \code{phyloscanner.trees} produced by \code{phyloscanner.analyse.trees}.
#' @param dwin A data frame produced by \code{classify.pairwise.relationships}.
#' @param tip.regex The regular expression used to identify host IDs in tip names
#' @param min.reads The minimum number of reads from a host in a window needed in order for that window to count in determining relationships involving that patient
#' @param min.tips The minimum number of tips from a host in a window needed in order for that window to count in determining relationships involving that patient
#' @param verbose Verbose output
#' @return A data frame with viral phylogenetic classifications of pairwise host relationships in each deep sequence phylogeny
#' @examples
#' \dontrun{
#' require(phyloscannerR)
#' #
#' #	continue Rakai example,
#' #	load phyloscanner output from 'phyloscanner.analyse.trees'
#' #	
#' file	<- system.file(file.path('extdata','ptyr192_phsc_analyse_trees_output.R'),package='phyloscannerR')
#' load(file)	#loads 'phsc', output from 'phyloscanner.analyse.trees'
#' #	use distance thresholds found in analysis of Rakai couples
#' close.threshold <- 0.025
#' distant.threshold <- 0.05	
#' #	use relationship types based on adjacency
#' #	this also considers linkage etc between individuals who have dual infections, recombinants etc 
#' #	..and thus may not have *all* their subgraphs adjacent to each other
#' relationship.types <- c('proximity.3.way',					
#' 		'close.and.adjacent',					
#' 		'close.and.adjacent.and.directed',					
#' 		'close.and.adjacent.and.ancestry.cat')	
#' dwin <- classify.pairwise.relationships(phsc, allow.mt=TRUE, close.threshold=close.threshold, distant.threshold=distant.threshold,relationship.types=relationship.types, verbose=TRUE)
#' tip.regex <- "^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$"
#' min.reads <- 30
#' min.tips <- 1
#' dwin <- select.windows.by.read.and.tip.count(ptrees, dwin, tip.regex, min.reads, min.tips)
#' #
#' # 	end of Rakai example
#' #
#' }	 	
select.windows.by.read.and.tip.count	<- function(ptrees, dwin, tip.regex, min.reads, min.tips, verbose=TRUE)			
{
	host.tips.and.reads <- map(ptrees, function(x) phyloscannerR:::get.tip.and.read.counts(x, all.hosts.from.trees(ptrees), tip.regex, attr(ptrees, 'has.read.counts'), verbose))
	host.tips.and.reads <- bind_rows(host.tips.and.reads)
	
	if(verbose) cat('Merging tip and read counts...\n')
	
	dwin <- dwin %>% 
			inner_join(host.tips.and.reads, by=c("host.1"="host.id", "tree.id")) %>% 
			rename(tips.1 = tips, reads.1=reads)
	dwin <- dwin %>% 
			inner_join(host.tips.and.reads, by=c("host.2"="host.id", "tree.id")) %>% 
			rename(tips.2 = tips, reads.2=reads)	
	
	if(verbose) cat('Reducing transmission window stats to windows with at least',min.reads,'reads and at least',min.tips,'tips...\n')
	
	dwin	<- dwin %>% 
			filter(reads.1 >= min.reads & reads.2 >= min.reads & tips.1 >= min.tips & tips.2 >= min.tips)
	
	if(verbose) cat('Total number of windows with transmission assignments is ',nrow(dwin),'.\n', sep="")		
	
	dwin
}