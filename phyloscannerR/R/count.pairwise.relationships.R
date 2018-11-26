#' @title Count pairwise relationships across deep-sequence trees 
#' @export 
#' @import tidyverse 
#' @author Oliver Ratmann, Matthew Hall
#' @param dwin A data frame produced by \code{classify.pairwise.relationships}.
#' @param w.slide Increment between genomic windows. Default: NA.
#' @param verbose Verbose output. Default: TRUE.
#' @return A data frame with counts of viral phylogenetic classifications between pairs of individuals.
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
#' dwin <- select.windows.by.read.and.tip.count(phsc, dwin, tip.regex, min.reads, min.tips)
#' # count phylogenetic relationships across deep-sequence trees
#' dc <- count.pairwise.relationships(dwin)
#' #
#' # 	end of Rakai example
#' #
#' }	
count.pairwise.relationships<- function(dwin, w.slide=NA, verbose=TRUE)	
{
	stopifnot(c('host.1', 'host.2', 'window.start', 'window.end', 'basic.classification','categorical.distance') %in% colnames(dwin))
	
	if(verbose) cat('Finding window coordinates...\n')
	
	dwin <- dwin %>% 
			mutate(window.start = map_int(tree.id, function(x) as.integer(gsub('[^0-9]*([0-9]+)_to_([0-9]+).*','\\1', x))))
	dwin <- dwin %>% 
			mutate(window.end = map_int(tree.id, function(x) as.integer(gsub('[^0-9]*([0-9]+)_to_([0-9]+).*','\\2', x))))
	
	#	identify chunks of contiguous windows
	dwin <- dwin %>% arrange(host.1, host.2, window.start)
	if(is.na(w.slide)) {
		# There are warnings here but it's OK unless _every_ pair in the dataset occurs in only one window
		
		w.slide <- dwin %>% group_by(host.1, host.2) %>% summarise(slide.width = min(diff(window.start))) %>% ungroup()
		
		w.slide	<- w.slide %>% filter(!is.na(slide.width))
		w.slide	<- ifelse(nrow(w.slide), min(w.slide$slide.width), 1L)		
		
		if(w.slide == Inf){
			stop("Cannot calculate sliding window size from this data")
		}
	}
	
	#	define chunks	
	dwin <- dwin %>% 
			group_by(host.1, host.2) %>% 
			mutate(new.block = (function(x){
							if_else(is.na(lag(x,1)), TRUE, lag(x, 1) != x-w.slide) 
						})(window.start)) %>% 
			mutate(chunk.no = cumsum(new.block)) %>% 
			select(-new.block) %>% 
			ungroup()
	
	#	define chunk length in terms of non-overlapping windows	& number of windows in chunk	
	dwin <- dwin %>% 
			group_by(host.1, host.2, chunk.no) %>% 
			mutate(effective.n.windows = (max(window.end) + 1 - min(window.start))/(window.end + 1 - window.start), n.windows = n()) %>% 
			ungroup()
	
	#	for each chunk, count: windows by type and effective length of chunk
	#	then sum across chunks
	#	then convert to long format
	dc	<- dwin %>% 
			select(host.1, host.2, tree.id, chunk.no, window.start, basic.classification, categorical.distance, n.windows, effective.n.windows)	%>% 
			arrange(host.1, host.2, window.start) %>%
			group_by(host.1, host.2, chunk.no, basic.classification, categorical.distance) %>%
			summarise(k = n(), k.eff = n()/n.windows[1] * effective.n.windows[1]) %>%
			group_by(host.1, host.2, basic.classification, categorical.distance) %>%
			summarise(k = sum(k), k.eff=sum(k.eff)) %>%
			gather(statistic, value, k:k.eff) %>%
			ungroup()
	
	#	now determine effective counts of all other relationship types
	#	first determine mapping of basic.relationship with other relationship types	
	relationship.types	<- colnames(dwin)[ grepl('*cat$',colnames(dwin)) ]
	relationship.types	<- c("basic.classification","categorical.distance",relationship.types)
	tmp <- dwin %>% 
			select(relationship.types) %>% 
			distinct() %>%
			gather(cat, type, -basic.classification, -categorical.distance)
	#	second, merge mapping
	#	then count k and keff for each relationship type
	tmp	<- dc %>% 
			inner_join(tmp, by=c("basic.classification","categorical.distance")) %>%
			filter(!is.na(type)) %>%
			group_by(host.1, host.2, cat, type, statistic) %>%
			summarise(value = sum(value)) %>% 
			spread(statistic, value)
	#	third, cast pairwise.relationship and row bind with other results
	dc	<- dc %>%
			rename(type=basic.classification) %>%
			mutate(cat:='basic.classification') %>%
			spread(statistic, value) %>%
			bind_rows(tmp)
	
	# 	expand to include zeros	
	tmp <- dc %>% 
			select(cat, categorical.distance, type) %>% 
			distinct() %>%
			mutate(dummy:=1L)
	dc	<- dc %>% 
			select(host.1,host.2) %>% 
			distinct() %>%
			mutate(dummy:=1L) %>%
			group_by(dummy) %>%
			inner_join(tmp, by='dummy') %>%
			ungroup() %>%
			select(-dummy) %>%
			group_by(host.1,host.2,cat,categorical.distance,type) %>%
			left_join(dc, by=c('host.1','host.2','cat','categorical.distance','type')) %>%
			replace_na(list(k=0, k.eff=0)) %>%
			ungroup()
	
	# add totals and sort
	tmp	<- dc %>%
			group_by(host.1,host.2,cat) %>%
			summarise(n=sum(k), n.eff=sum(k.eff)) %>%
			ungroup() 
	dc	<- dc %>% 
			group_by(host.1,host.2,cat) %>%
			inner_join(tmp, by=c('host.1','host.2','cat')) %>%
			ungroup() %>%
			rename(categorisation=cat) %>%
			arrange(host.1,host.2,cat,categorical.distance,type) 			
	dc
}