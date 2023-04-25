#' @title Count pairwise relationships across deep-sequence trees 
#' @export 
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
#' dwin <- classify.pairwise.relationships(phsc, close.threshold=close.threshold, distant.threshold=distant.threshold,relationship.types=relationship.types, verbose=TRUE)
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
  stopifnot(c('host.1', 'host.2', 'basic.classification','categorical.distance') %in% colnames(dwin))
  
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
    arrange(host.1,host.2,categorisation,categorical.distance,type) %>%
    mutate(score = k.eff/n.eff)
  dc
}


#' @title Classify pairwise host relationships in deep sequence phylogenies 
#' @export
#' @author Oliver Ratmann, Matthew Hall
#' @param ptrees A list of class \code{phyloscanner.trees} produced by \code{phyloscanner.analyse.trees}.
#' @param close.threshold The (potentially normalised) patristic threshold used to determine if two patients' subgraphs are "close".
#' @param distant.threshold If present, a second distance threshold determines hosts that are "distant" from each other, with those lying between \code{close.threshold} and \code{dist.threshold} classed as "intermediate". The default is the same as \code{close.threshold}, so the intermediate class does not exist.
#' @param relationship.types Classification types.
#' \itemize{
#'  \item{"proximity.3.way"}{Classify individuals by phylogenetic distance between subgraphs. Suggested use: to exclude phylogenetic linkage based on distance alone.}
#'  \item{"close.and.contiguous"}{Classify individuals by phylogenetic distance and contiguity of subgraphs. Suggested use: to identify phylogenetically linked pairs.}
#'  \item{"close.and.adjacent"}{Classify individuals by phylogenetic distance and adjacency of subgraphs. Suggested use: to identify phylogenetically linked pairs.}
#'  \item{"close.and.contiguous.and.directed"}{Classify ancestry among contiguous subgraphs. Suggested use: to identify direction of transmission based on contiguous subgraphs.}
#'  \item{"close.and.adjacent.and.directed"}{Classify ancestry among adjacent subgraphs. Suggested use: to identify direction of transmission based on adjacent subgraphs.} 
#'  \item{"close.and.contiguous.and.ancestry.cat"}{Classify contiguity and ancestry between individuals. Suggested use: to determine probabilities for transmission networks.}
#'  \item{"close.and.adjacent.and.ancestry.cat"}{Classify adjacency and ancestry between individuals. Suggested use: to determine probabilities for transmission networks.}
#' }
#' @param verbose Verbose output
#' @return A data frame with viral phylogenetic classifications of pairwise host relationships in each deep sequence phylogeny
#' @examples
#' \dontrun{
#' require(phyloscannerR)
#' #
#' # Example on data from Rakai Community Cohort Study
#' # load phyloscanner output from 'phyloscanner.analyse.trees'
#' #	
#' file	<- system.file(file.path('extdata','ptyr192_phsc_analyse_trees_output.RData'),package='phyloscannerR')
#' load(file)	#loads 'phsc', output from 'phyloscanner.analyse.trees'
#' #	use distance thresholds found in analysis of Rakai couples
#' close.threshold <- 0.025
#' distant.threshold <- 0.05	
#' #	use relationship types based on adjacency
#' #	this also considers linkage etc between individuals who have dual infections, recombinants etc 
#' #	..and thus may not have *all* their subgraphs adjacent to each other
#' relationship.types <- c('close.and.adjacent',					
#' 		'close.and.adjacent.and.directed',					
#' 		'close.and.adjacent.and.ancestry.cat')	
#' dwin <- classify.pairwise.relationships(phsc, 
#'            close.threshold=close.threshold, 
#'            distant.threshold=distant.threshold,
#'            relationship.types=relationship.types, 
#'            verbose=TRUE)
#' }	
classify.pairwise.relationships<- function(ptrees, 
                                           close.threshold=0.025, 
                                           distant.threshold=0.05,	
                                           relationship.types=c('proximity.3.way',
                                                                'any.ancestry',
                                                                'close.x.contiguous',
                                                                'close.and.contiguous',
                                                                'close.and.adjacent',
                                                                'close.and.contiguous.and.directed',
                                                                'close.and.adjacent.and.directed',
                                                                'close.and.contiguous.and.ancestry.cat',
                                                                'close.and.adjacent.and.ancestry.cat'),	
                                           verbose = FALSE
)
{		
  
  dwin <- phyloscannerR:::merge.classifications(ptrees, verbose)	
  
  if(!all(dwin$host.1 < dwin$host.2)){
    stop("Failed to rearrange pair orders alphabetically.")
  }

  if(verbose) cat('Assigning discrete proximity categories...\n')
  
  dwin <- dwin %>% 
    mutate(categorical.distance = map_chr(patristic.distance, function(x){
      if(x < close.threshold){
        "close"
      } else if(x >= distant.threshold){
        "distant"
      } else {
        "intermediate"
      }
    }))
  
  if(verbose) cat('Calculating basic pairwise relationships for windows (n=',nrow(dwin),')...\n', sep="")
  
  dwin	<- phyloscannerR:::get.pairwise.relationships.basic.by.ancestry(dwin)
  
  if(verbose) cat('Calculating derived pairwise relationships for windows (n=',nrow(dwin),')...\n', sep="")
  
  if('proximity.3.way' %in% relationship.types)
  {
    dwin <- dwin %>% 
      select(host.1,host.2,tree.id,categorical.distance) %>%
      phyloscannerR:::get.pairwise.relationships.proximity.3.way.cat() %>%
      group_by(host.1,host.2,tree.id) %>%
      inner_join(dwin, by=c('host.1','host.2','tree.id')) %>%
      ungroup() %>%
      select(-proximity.3.way.cat,proximity.3.way.cat)				
  }
  if('any.ancestry' %in% relationship.types)
  {
    dwin <- dwin %>% 
      select(host.1,host.2,tree.id,contiguous,basic.classification) %>%
      phyloscannerR:::get.pairwise.relationships.any.ancestry.cat() %>%
      group_by(host.1,host.2,tree.id) %>%
      inner_join(dwin, by=c('host.1','host.2','tree.id')) %>%
      ungroup() %>%
      select(-any.ancestry.cat,any.ancestry.cat)
  }
  if('close.x.contiguous' %in% relationship.types)
  {
    dwin <- dwin %>% 
      select(host.1,host.2,tree.id,contiguous,categorical.distance) %>%
      phyloscannerR:::get.pairwise.relationships.close.x.contiguous.cat() %>%
      group_by(host.1,host.2,tree.id) %>%
      inner_join(dwin, by=c('host.1','host.2','tree.id')) %>%
      ungroup() %>%
      select(-close.x.contiguous.cat,close.x.contiguous.cat)
  }	
  if('close.and.contiguous' %in% relationship.types)
  {
    dwin <- dwin %>% 
      select(host.1,host.2,tree.id,contiguous,categorical.distance) %>%
      phyloscannerR:::get.pairwise.relationships.close.and.contiguous.cat() %>%
      group_by(host.1,host.2,tree.id) %>%
      inner_join(dwin, by=c('host.1','host.2','tree.id')) %>%
      ungroup() %>%
      select(-close.and.contiguous.cat,close.and.contiguous.cat)
  }
  if('close.and.adjacent' %in% relationship.types)
  {
    dwin <- dwin %>% 
      select(host.1,host.2,tree.id,adjacent,categorical.distance) %>%
      phyloscannerR:::get.pairwise.relationships.close.and.adjacent.cat() %>%
      group_by(host.1,host.2,tree.id) %>%
      inner_join(dwin, by=c('host.1','host.2','tree.id')) %>%
      ungroup() %>%
      select(-close.and.adjacent.cat,close.and.adjacent.cat)
  }
  if('close.and.contiguous.and.directed' %in% relationship.types)
  {
    dwin <- dwin %>% 
      select(host.1,host.2,tree.id,contiguous,categorical.distance,basic.classification) %>%
      phyloscannerR:::get.pairwise.relationships.close.and.contiguous.and.directed.cat() %>%
      group_by(host.1,host.2,tree.id) %>%
      inner_join(dwin, by=c('host.1','host.2','tree.id')) %>%
      ungroup() %>%
      select(-close.and.contiguous.and.directed.cat,close.and.contiguous.and.directed.cat)
  }
  if('close.and.adjacent.and.directed' %in% relationship.types)
  {
    dwin <- dwin %>% 
      select(host.1,host.2,tree.id,adjacent,categorical.distance,basic.classification) %>%
      phyloscannerR:::get.pairwise.relationships.close.and.adjacent.and.directed.cat() %>%
      group_by(host.1,host.2,tree.id) %>%
      inner_join(dwin, by=c('host.1','host.2','tree.id')) %>%
      ungroup() %>%
      select(-close.and.adjacent.and.directed.cat,close.and.adjacent.and.directed.cat)
  }
  if('close.and.contiguous.and.ancestry' %in% relationship.types)
  {
    dwin <- dwin %>% 
      select(host.1,host.2,tree.id,contiguous,categorical.distance,basic.classification) %>%
      phyloscannerR:::get.pairwise.relationships.close.and.contiguous.and.ancestry.cat() %>%
      group_by(host.1,host.2,tree.id) %>%
      inner_join(dwin, by=c('host.1','host.2','tree.id')) %>%
      ungroup() %>%
      select(-close.and.contiguous.and.ancestry.cat,close.and.contiguous.and.ancestry.cat)
  }
  if('close.and.adjacent.and.ancestry' %in% relationship.types)
  {
    dwin <- dwin %>% 
      select(host.1,host.2,tree.id,adjacent,categorical.distance,basic.classification) %>%
      phyloscannerR:::get.pairwise.relationships.close.and.adjacent.and.ancestry.cat() %>%
      group_by(host.1,host.2,tree.id) %>%
      inner_join(dwin, by=c('host.1','host.2','tree.id')) %>%
      ungroup() %>%
      select(-close.and.adjacent.and.ancestry.cat,close.and.adjacent.and.ancestry.cat)
  }
  dwin
}	


#' @keywords internal
#' @title Define basic topological relationships between two individuals by ancestry
#' @author Oliver Ratmann
get.pairwise.relationships.basic.by.ancestry <- function(all.classifications)
{
  stopifnot( all(c('contiguous','adjacent','paths12','paths21') %in% colnames(all.classifications)) )
  # define "XXX_contiguous"
  
  tmp			<- all.classifications %>% filter(contiguous & adjacent)
  class 		<- list()
  class[[1]]	<- tmp %>% filter(ancestry == "anc" | ancestry == "multiAnc") %>% mutate(basic.classification:= "anc_contiguous")
  class[[2]]	<- tmp %>% filter(ancestry == "desc" | ancestry == "multiDesc") %>% mutate(basic.classification:= "desc_contiguous")
  class[[3]]	<- tmp %>% filter(ancestry=='complex') %>% mutate(basic.classification:= "intermingled_contiguous")
  class[[4]]	<- tmp %>% filter(ancestry=='noAncestry') %>% mutate(basic.classification:= "sibling_contiguous")
  
  # define "XXX_noncontiguous"
  tmp			<- all.classifications %>% filter(!contiguous & adjacent)
  class[[5]]	<- tmp %>% filter(ancestry == "anc" | ancestry == "multiAnc") %>% mutate(basic.classification:= "anc_noncontiguous")
  class[[6]]	<- tmp %>% filter(ancestry == "desc" | ancestry == "multiDesc") %>% mutate(basic.classification:= "desc_noncontiguous")
  class[[7]]	<- tmp %>% filter(ancestry=='complex') %>% mutate(basic.classification:= "intermingled_noncontiguous")
  class[[8]]	<- tmp %>% filter(ancestry=='noAncestry') %>% mutate(basic.classification:= "sibling_noncontiguous")
  # define "other"
  class[[9]]	<- all.classifications %>% filter(!adjacent) %>% mutate(basic.classification:= "other")
  
  class		<- bind_rows(class)
  
  stopifnot(nrow(all.classifications)==nrow(class))
  class
}

#' @keywords internal
#' @title Group to determine 2x2 table of distance and topology
#' @author Oliver Ratmann, Matthew Hall
get.pairwise.relationships.any.ancestry.cat <- function(all.classifications)
{
  stopifnot( all(c('host.1','host.2','tree.id','basic.classification','contiguous')%in%colnames(all.classifications)) )	
  class 		<- list()
  class[[1]]	<- all.classifications %>% 
    filter(contiguous & grepl('anc|desc|intermingled',basic.classification)) %>% 
    mutate(any.ancestry.cat:= 'anc.or.complex')		
  class[[2]]	<- all.classifications %>% 
    filter(!(contiguous & grepl('anc|desc|intermingled',basic.classification))) %>% 
    mutate(any.ancestry.cat:= 'other')	
  class		<- bind_rows(class)
  stopifnot(nrow(class)==nrow(all.classifications))
  class 		<- class %>% select(host.1,host.2,tree.id,any.ancestry.cat)
  class	
}


#' @keywords internal
#' @title Group to determine 2x2 table of distance and topology
#' @author Oliver Ratmann, Matthew Hall
get.pairwise.relationships.close.x.contiguous.cat <- function(all.classifications)
{
  stopifnot( all(c('host.1','host.2','tree.id','categorical.distance','contiguous')%in%colnames(all.classifications)) )	
  class 		<- list()
  class[[1]]	<- all.classifications %>% filter(categorical.distance!='close' & !contiguous) %>% mutate(close.x.contiguous.cat:= 'not.close.other')	
  class[[2]]	<- all.classifications %>% filter(categorical.distance=='close' & !contiguous) %>% mutate(close.x.contiguous.cat:= 'noncontiguous_close')
  class[[3]]	<- all.classifications %>% filter(categorical.distance!='close' & contiguous) %>% mutate(close.x.contiguous.cat:= 'contiguous_not.close')
  class[[4]]	<- all.classifications %>% filter(categorical.distance=='close' & contiguous) %>% mutate(close.x.contiguous.cat:= 'contiguous_close')	
  class		<- bind_rows(class)
  stopifnot(nrow(class)==nrow(all.classifications))
  class 		<- class %>% select(host.1,host.2,tree.id,close.x.contiguous.cat)
  class	
}


#' @keywords internal
#' @title Group to determine phylogenetic linkage based on distance only
#' @author Oliver Ratmann, Matthew Hall
get.pairwise.relationships.proximity.3.way.cat <- function(all.classifications)
{
  stopifnot( all(c('host.1','host.2','tree.id','categorical.distance')%in%colnames(all.classifications)) )	
  class <- all.classifications %>% mutate(proximity.3.way.cat:= categorical.distance)
  class <- class %>% select(host.1,host.2,tree.id,proximity.3.way.cat)
  class
}

#' @keywords internal
#' @title Group to determine phylogenetic linkage based on adjacent subgraphs
#' @author Oliver Ratmann, Matthew Hall
get.pairwise.relationships.close.and.adjacent.cat <- function(all.classifications)
{
  stopifnot( all(c('host.1','host.2','tree.id','adjacent','categorical.distance')%in%colnames(all.classifications)) )
  class 		<- list()
  class[[1]]	<- all.classifications %>% filter(categorical.distance=='close' & adjacent) %>% mutate(close.and.adjacent.cat:= 'close.and.adjacent')
  class[[2]]	<- all.classifications %>% filter(!(categorical.distance=='close' & adjacent)) %>% mutate(close.and.adjacent.cat:= 'not.close.or.nonadjacent')
  class		<- bind_rows(class)
  stopifnot(nrow(class)==nrow(all.classifications))
  class 		<- class %>% select(host.1,host.2,tree.id,close.and.adjacent.cat)
  class
}

#' @keywords internal
#' @title Group to determine phylogenetic linkage based on contiguous subgraphs
#' @author Oliver Ratmann, Matthew Hall
get.pairwise.relationships.close.and.contiguous.cat <- function(all.classifications)
{
  stopifnot( all(c('host.1','host.2','tree.id','contiguous','categorical.distance')%in%colnames(all.classifications)) )
  class 		<- list()
  class[[1]]	<- all.classifications %>% filter(categorical.distance=='close' & contiguous) %>% mutate(close.and.contiguous.cat:= 'close.and.contiguous')
  class[[2]]	<- all.classifications %>% filter(!(categorical.distance=='close' & contiguous)) %>% mutate(close.and.contiguous.cat:= 'not.close.or.noncontiguous')
  class		<- bind_rows(class)
  stopifnot(nrow(class)==nrow(all.classifications))
  class 		<- class %>% select(host.1,host.2,tree.id,close.and.contiguous.cat)
  class
}

#' @keywords internal
#' @title Group to determine probabilities for transmission networks based on adjacent subgraphs
#' @description Define basic topological relationships 'close and 12', 'close and 21', 'close and ambiguous direction', 'not close or not adjacent'
#' @author Oliver Ratmann, Matthew Hall
get.pairwise.relationships.close.and.adjacent.and.ancestry.cat <- function(all.classifications)
{
  stopifnot( all(c('host.1','host.2','tree.id','adjacent','categorical.distance','basic.classification')%in%colnames(all.classifications)) )
  tmp			<- all.classifications %>% filter(categorical.distance=='close' & adjacent)
  class 		<- list()
  class[[1]] 	<- tmp %>% filter(grepl('anc',basic.classification)) %>% mutate(close.and.adjacent.and.ancestry.cat:= '12')
  class[[2]] 	<- tmp %>% filter(grepl('desc',basic.classification)) %>% mutate(close.and.adjacent.and.ancestry.cat:= '21')
  class[[3]] 	<- tmp %>% filter(grepl('sibling|intermingled',basic.classification)) %>% mutate(close.and.adjacent.and.ancestry.cat:= 'complex.or.no.ancestry')
  class[[4]] 	<- all.classifications %>% filter(!(categorical.distance=='close' & adjacent)) %>% mutate(close.and.adjacent.and.ancestry.cat:= 'not.close.or.nonadjacent')
  class		<- bind_rows(class)
  stopifnot(nrow(class)==nrow(all.classifications))
  class 		<- class %>% select(host.1,host.2,tree.id,close.and.adjacent.and.ancestry.cat)
  class			
}

#' @keywords internal
#' @title Group to determine direction of transmission in likely pairs based on adjacent subgraphs
#' @author Oliver Ratmann, Matthew Hall 
get.pairwise.relationships.close.and.adjacent.and.directed.cat <- function(all.classifications)
{
  stopifnot( all(c('host.1','host.2','tree.id','adjacent','categorical.distance','basic.classification')%in%colnames(all.classifications)) )
  tmp			<- all.classifications %>% filter(categorical.distance=='close' & adjacent)
  class 		<- list()
  class[[1]] 	<- tmp %>% filter(grepl('anc',basic.classification)) %>% mutate(close.and.adjacent.and.directed.cat:= '12')
  class[[2]] 	<- tmp %>% filter(grepl('desc',basic.classification)) %>% mutate(close.and.adjacent.and.directed.cat:= '21')
  class[[3]] 	<- tmp %>% filter(grepl('sibling|intermingled',basic.classification)) %>% mutate(close.and.adjacent.and.directed.cat:= NA_character_)
  class[[4]] 	<- all.classifications %>% filter(!(categorical.distance=='close' & adjacent)) %>% mutate(close.and.adjacent.and.directed.cat:= NA_character_)
  class		<- bind_rows(class)
  stopifnot(nrow(class)==nrow(all.classifications))
  class 		<- class %>% select(host.1,host.2,tree.id,close.and.adjacent.and.directed.cat)
  class			
}

#' @keywords internal
#' @title Group to determine direction of transmission in likely pairs based on contiguous subgraphs
#' @author Oliver Ratmann, Matthew Hall
get.pairwise.relationships.close.and.contiguous.and.directed.cat <- function(all.classifications)
{
  stopifnot( all(c('host.1','host.2','tree.id','contiguous','categorical.distance','basic.classification')%in%colnames(all.classifications)) )
  tmp			<- all.classifications %>% filter(categorical.distance=='close' & contiguous)
  class 		<- list()
  class[[1]] 	<- tmp %>% filter(grepl('anc',basic.classification)) %>% mutate(close.and.contiguous.and.directed.cat:= '12')
  class[[2]] 	<- tmp %>% filter(grepl('desc',basic.classification)) %>% mutate(close.and.contiguous.and.directed.cat:= '21')
  class[[3]] 	<- tmp %>% filter(grepl('sibling|intermingled',basic.classification)) %>% mutate(close.and.contiguous.and.directed.cat:= NA_character_)
  class[[4]] 	<- all.classifications %>% filter(!(categorical.distance=='close' & contiguous)) %>% mutate(close.and.contiguous.and.directed.cat:= NA_character_)
  class		<- bind_rows(class)
  stopifnot(nrow(class)==nrow(all.classifications))
  class 		<- class %>% select(host.1,host.2,tree.id,close.and.contiguous.and.directed.cat)
  class			
}

#' @keywords internal
#' @title Group to determine probabilities for transmission networks based on contiguous subgraphs
#' @description Define basic topological relationships 'close and 12', 'close and 21', 'close and ambiguous direction', 'not close or not contiguous'
#' @author Oliver Ratmann, Matthew Hall
get.pairwise.relationships.close.and.contiguous.and.ancestry.cat <- function(all.classifications)
{
  stopifnot( all(c('host.1','host.2','tree.id','contiguous','categorical.distance','basic.classification')%in%colnames(all.classifications)) )
  tmp			<- all.classifications %>% filter(categorical.distance=='close' & contiguous)
  class 		<- list()
  class[[1]] 	<- tmp %>% filter(grepl('anc',basic.classification)) %>% mutate(close.and.contiguous.and.ancestry.cat:= '12')
  class[[2]] 	<- tmp %>% filter(grepl('desc',basic.classification)) %>% mutate(close.and.contiguous.and.ancestry.cat:= '21')
  class[[3]] 	<- tmp %>% filter(grepl('sibling|intermingled',basic.classification)) %>% mutate(close.and.contiguous.and.ancestry.cat:= 'complex.or.no.ancestry')
  class[[4]] 	<- all.classifications %>% filter(!(categorical.distance=='close' & contiguous)) %>% mutate(close.and.contiguous.and.ancestry.cat:= 'not.close.or.noncontiguous')
  class		<- bind_rows(class)
  stopifnot(nrow(class)==nrow(all.classifications))
  class 		<- class %>% select(host.1,host.2,tree.id,close.and.contiguous.and.ancestry.cat)
  class			
}


#' @keywords internal
#' @export summarise.classifications

summarise.classifications <- function(ptrees, min.threshold, dist.threshold, tip.regex, min.reads = 1, min.tips = 1, close.sib.only = F, verbose = F, contiguous = F){

  tt <- merge.classifications(ptrees, verbose)

  if(min.reads > 1 | min.tips > 1){
    tt <- select.windows.by.read.and.tip.count(ptrees, tt, tip.regex, min.reads, min.tips, verbose)
    tt <- tt %>% select(-reads.1, -reads.2, -tips.1, -tips.2)
    
  }
  
  if (verbose) cat("Making summary output table...\n")
  
  tt <- tt %>% select(-c(paths12, paths21))
  
  existence.counts <- tt %>% count(host.1, host.2)
  
  tt <- tt %>% 
    group_by(host.1, host.2) %>% 
    mutate(both.exist = n()) %>%
    ungroup()
  
  if(!close.sib.only){
    if(!contiguous){
      tt.close <- tt %>% 
        filter(adjacent & patristic.distance < dist.threshold)
    } else {
      tt.close <- tt %>% 
        filter(contiguous & patristic.distance < dist.threshold)
    }
  } else {
    if(!contiguous){
      tt.close <- tt %>% 
        filter(adjacent & (patristic.distance < dist.threshold | ancestry != "noAncestry"))
    } else {
      tt.close <- tt %>% 
        filter(contiguous & (patristic.distance < dist.threshold | ancestry != "noAncestry"))
    }
  }
  
  # How many windows have this relationship, adjacent and patristic.distance below the threshold?
  tt.close	<- tt.close %>% 
    group_by(host.1, host.2, ancestry) %>% 
    mutate(ancestry.tree.count = n()) %>% 
    ungroup()
  
  # How many windows have adjacent and patristic.distance below the threshold?
  
  tt.close <- tt.close %>% group_by(host.1, host.2) %>% 
    mutate(any.relationship.tree.count = n()) %>% 
    ungroup()
  
  tt.close <- tt.close %>% 
    mutate(fraction = paste0(ancestry.tree.count,'/',both.exist))
  
  #	Flip directions - all ancestry is now "trans" or "multiTrans"
  
  tt.close <- tt.close %>% 
    mutate(new.host.1 = ifelse(ancestry %in% c("desc", "multiDesc"), host.2, host.1), 
           new.host.2 = ifelse(ancestry %in% c("desc", "multiDesc"), host.1, host.2), 
           host.1 = new.host.1,
           host.2 = new.host.2) %>% 
    select(-new.host.1, -new.host.2) %>% 
    mutate(ancestry = replace(ancestry, ancestry == "anc" | ancestry == "desc", "trans")) %>% 
    mutate(ancestry = replace(ancestry, ancestry == "multiAnc" | ancestry == "multiDesc", "multiTrans")) %>% 
    select(-tree.id, -adjacent, -contiguous, -nodes1, -nodes2, -patristic.distance) %>% 
    distinct() %>% 
    arrange(host.1, host.2, ancestry) %>% 
    select(host.1, host.2, ancestry, ancestry.tree.count, both.exist, fraction, any.relationship.tree.count) %>% 
    filter(any.relationship.tree.count > min.threshold)
  
  tt.close 
}

#' @keywords internal
#' @importFrom network as.network.matrix
#' @export simplify.summary


simplify.summary <- function(summary, arrow.threshold, total.trees, plot = F){
  
  done <- rep(FALSE, nrow(summary))
  
  for(line in 1:nrow(summary)){
    if(!done[line]){
      forwards.rows <- which(summary$host.1 == summary$host.1[line] & summary$host.2 == summary$host.2[line])
      backwards.rows <- which(summary$host.1 == summary$host.2[line] & summary$host.2 == summary$host.1[line])
      
      done[forwards.rows] <- T
      done[backwards.rows] <- T
      
      summary$ancestry[intersect(forwards.rows, which(summary$ancestry=="trans"))] <- "trans12"
      summary$ancestry[intersect(forwards.rows, which(summary$ancestry=="multiTrans"))] <- "multi_trans12"
      
      summary$ancestry[intersect(backwards.rows, which(summary$ancestry=="trans"))] <- "trans21"
      summary$ancestry[intersect(backwards.rows, which(summary$ancestry=="multiTrans"))] <- "multi_trans21"
      
      summary[backwards.rows,c(1,2)] <- summary[backwards.rows,c(2,1)]
    }
  }
  
  summary.wide <- summary %>% 
    select(host.1, host.2, both.exist, ancestry, ancestry.tree.count) %>% 
    spread(ancestry, ancestry.tree.count, fill = 0) %>% 
    mutate(total.equiv = 0, total.12 = 0, total.21 = 0)
  
  
  if("none" %in% names(summary.wide)){
    summary.wide <- summary.wide %>% mutate(total.equiv = total.equiv + none)
  }
  if("complex" %in% names(summary.wide)){
    summary.wide <- summary.wide %>% mutate(total.equiv = total.equiv + complex)
  }
  if("trans12" %in% names(summary.wide)){
    summary.wide <- summary.wide %>% mutate(total.12 = total.12 + trans12)
  } 
  if("multi_trans12" %in% names(summary.wide)){
    summary.wide <- summary.wide %>% mutate(total.12 = total.12 + multi_trans12)
  }
  if("trans21" %in% names(summary.wide)){
    summary.wide <- summary.wide %>% mutate(total.21 = total.21 + trans21)
  } 
  if("multi_trans21" %in% names(summary.wide)){
    summary.wide <- summary.wide %>% mutate(total.21 = total.21 + multi_trans21)
  }
  
  summary.wide <- summary.wide %>% 
    mutate(total = total.21 + total.12 + total.equiv,
           dir = total.12 >= arrow.threshold*total.trees | total.21 >= arrow.threshold*total.trees,
           arrow = pmap_chr(list(dir, total.12, total.21), function(x, y, z) {
             if(!x){
               return("none")
             }
             if(y > z){
               return("forwards")
             } else {
               return("backwards")
             }
           }),
           label = pmap_chr(list(dir, total.12, total.21, total), function(w, x, y, z) {
             if(!w){
               as.character(round(z/total.trees, 2))
             } else {
               paste0(round(max(x,y)/total.trees, 2), "/", round(z/total.trees, 2))
             }
             
           })
    ) %>%
    select(host.1, host.2, arrow, label)
  
  # this has not been tidied - but it's much more succinct in standard R language
  
  summary.wide[which(summary.wide$arrow=="backwards"),c(1,2)] <- summary.wide[which(summary.wide$arrow=="backwards"),c(2,1)] 
  
  summary.wide$arrow <- summary.wide$arrow!="none"
  
  out <- list(simp.table = summary.wide)
  
  if(plot){
    out$simp.diagram <- make.simplified.plot(summary.wide)
  }
  
  return(out)
}

#' @keywords internal
#' @export simplify.summary.multinomial

simplify.summary.multinomial <- function(summary, win.threshold, arrow.threshold, contiguous = F, plot = F){
  
  if(contiguous){
    relevant.lines <- summary %>%
      filter(categorisation == "close.and.contiguous.and.ancestry.cat")
  } else {
    relevant.lines <- summary %>%
      filter(categorisation == "close.and.adjacent.and.ancestry.cat")
  }
  
  relevant.lines <- relevant.lines %>% 
    filter(type != "not.close.or.nonadjacent" & type != "not.close.or.noncontiguous") %>%
    group_by(host.1, host.2) %>%
    mutate(total.score = sum(score)) %>% 
    filter(total.score >= win.threshold) 
  
  if(nrow(relevant.lines) == 0){
    return(NULL)
  }
  
  relevant.lines <- relevant.lines %>%
    select(host.1, host.2, type, score) %>%
    mutate(type = factor(type, levels = c("12", "21", "complex.or.no.ancestry"))) %>%
    spread(type, score, fill = 0, drop=F) %>%
    rename(total.12 = `12`, total.21 = `21`, total.equiv = complex.or.no.ancestry) %>% 
    mutate(total = total.21 + total.12 + total.equiv,
           dir = total.12 >= arrow.threshold | total.21 >= arrow.threshold,
           arrow = pmap_chr(list(dir, total.12, total.21), function(x, y, z) {
             if(!x){
               return("none")
             }
             if(y > z){
               return("forwards")
             } else {
               return("backwards")
             }
           }),
           label = pmap_chr(list(dir, total.12, total.21, total), function(w, x, y, z) {
             if(!w){
               as.character(round(z, 2))
             } else {
               paste0(round(max(x,y), 2), "/", round(z, 2))
             }
             
           })
    ) %>%
    filter(total > 0) %>%
    select(host.1, host.2, total.12, total.21, total, arrow, label)
  
  relevant.lines[which(relevant.lines$arrow=="backwards"),c(1,2)] <- relevant.lines[which(relevant.lines$arrow=="backwards"),c(2,1)] 
  
  relevant.lines$arrow <- relevant.lines$arrow!="none"
  
  out <- list(simp.table = relevant.lines)
  
  if(plot){
    out$simp.diagram <- make.simplified.plot(relevant.lines)
  }
  
  return(out)
}


#' @keywords internal
#' #' @importFrom network as.network.matrix
#' @export make.simplified.plot

make.simplified.plot <- function(summary){
  network.obj <- as.network.matrix(summary[,c(1,2)], matrix.type = "edgelist")
  
  arrangement <- ggnet2(network.obj)$data[,c("label", "x", "y")]
  
  summary$x.start <- map_dbl(summary$host.1, function(x) arrangement$x[match(x, arrangement$label)]) 
  summary$y.start <- map_dbl(summary$host.1, function(x) arrangement$y[match(x, arrangement$label)]) 
  summary$x.end <- map_dbl(summary$host.2, function(x) arrangement$x[match(x, arrangement$label)]) 
  summary$y.end <- map_dbl(summary$host.2, function(x) arrangement$y[match(x, arrangement$label)]) 
  summary$x.midpoint <- (summary$x.end + summary$x.start)/2
  summary$y.midpoint <- (summary$y.end + summary$y.start)/2
  
  out.diagram <- ggplot() + 
    geom_segment(data=summary[which(summary$arrow),], aes(x=x.start, xend = x.end, y=y.start, yend = y.end), arrow = arrow(length = unit(0.01, "npc"), type="closed"), col="steelblue3", size=1.5, lineend="round") +
    geom_segment(data=summary[which(!summary$arrow),], aes(x=x.start, xend = x.end, y=y.start, yend = y.end), col="chartreuse3", size=1.5, lineend="round") +
    geom_label(aes(x=arrangement$x, y=arrangement$y, label=arrangement$label), alpha=0.25, fill="darkgoldenrod3") + 
    geom_text(data=summary, aes(x=x.midpoint, y=y.midpoint, label=label)) + 
    theme_void()
  
  out.diagram
}

