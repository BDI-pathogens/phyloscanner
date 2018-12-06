#' @title Classify pairwise host relationships in deep sequence phylogenies Make viral phylogenetic classification 
#' @export
#' @import tidyverse 
#' @author Oliver Ratmann, Matthew Hall
#' @param ptrees A list of class \code{phyloscanner.trees} produced by \code{phyloscanner.analyse.trees}.
#' @param allow.mt If FALSE, directionality is only inferred between pairs of hosts where a single clade from one host is nested in one from the other; this is more conservative. Default: TRUE. 
#' @param close.threshold The (potentially normalised) patristic threshold used to determine if two patients' subgraphs are "close".
#' @param distant.threshold If present, a second distance threshold determines hosts that are "distant" from each other, with those lying between \code{close.threshold} and \code{dist.threshold} classed as "intermediate". The default is the same as \code{close.threshold}, so the intermediate class does not exist.
#' @param use.paths.to.define.relationships If true, the path variable only is used to define ancestry (paths12>0 & paths21==0), intermingled (paths12>0 & paths21>0), and sibling (paths12==0 & paths21==0) relationships. Default: FALSE. 
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
#' use.paths.to.define.relationships <- 1
#' relationship.types <- c('proximity.3.way',					
#' 		'close.and.adjacent',					
#' 		'close.and.adjacent.and.directed',					
#' 		'close.and.adjacent.and.ancestry.cat')	
#' dwin <- classify.pairwise.relationships(phsc, allow.mt=TRUE, close.threshold=close.threshold, distant.threshold=distant.threshold,use.paths.to.define.relationships=use.paths.to.define.relationships,relationship.types=relationship.types, verbose=TRUE)
#' #
#' # 	end of Rakai example
#' #
#' }	
classify.pairwise.relationships<- function(ptrees, 
		allow.mt=TRUE,
		close.threshold=0.025, 
		distant.threshold=0.05,	
		use.paths.to.define.relationships=FALSE,
		relationship.types=c('proximity.3.way',
				'any.ancestry',
				'close.x.contiguous',
				'close.and.contiguous',
				'close.and.adjacent',
				'close.and.contiguous.and.directed',
				'close.and.adjacent.and.directed',
				'close.and.contiguous.and.ancestry.cat',
				'close.and.adjacent.and.ancestry.cat'),	
		verbose= FALSE
)
{		
	dwin <- phyloscannerR:::merge.classifications(ptrees, allow.mt, verbose)	
	
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
	
	if(use.paths.to.define.relationships)
		dwin	<- phyloscannerR:::get.pairwise.relationships.basic.by.paths(dwin)
	if(!use.paths.to.define.relationships)
		dwin	<- phyloscannerR:::get.pairwise.relationships.basic.by.ancestry(dwin)
	
	if(verbose) cat('Calculating derived pairwise relationships for windows (n=',nrow(dwin),')...\n', sep="")
	
	if('proximity.3.way'%in%relationship.types)
	{
		dwin <- dwin %>% 
				select(host.1,host.2,tree.id,categorical.distance) %>%
				phyloscannerR:::get.pairwise.relationships.proximity.3.way.cat() %>%
				group_by(host.1,host.2,tree.id) %>%
				inner_join(dwin, by=c('host.1','host.2','tree.id')) %>%
				ungroup() %>%
				select(-proximity.3.way.cat,proximity.3.way.cat)				
	}
	if('any.ancestry'%in%relationship.types)
	{
		dwin <- dwin %>% 
				select(host.1,host.2,tree.id,contiguous,basic.classification) %>%
				phyloscannerR:::get.pairwise.relationships.any.ancestry.cat() %>%
				group_by(host.1,host.2,tree.id) %>%
				inner_join(dwin, by=c('host.1','host.2','tree.id')) %>%
				ungroup() %>%
				select(-any.ancestry.cat,any.ancestry.cat)
	}
	if('close.x.contiguous'%in%relationship.types)
	{
		dwin <- dwin %>% 
				select(host.1,host.2,tree.id,contiguous,categorical.distance) %>%
				phyloscannerR:::get.pairwise.relationships.close.x.contiguous.cat() %>%
				group_by(host.1,host.2,tree.id) %>%
				inner_join(dwin, by=c('host.1','host.2','tree.id')) %>%
				ungroup() %>%
				select(-close.x.contiguous.cat,close.x.contiguous.cat)
	}	
	if('close.and.contiguous'%in%relationship.types)
	{
		dwin <- dwin %>% 
				select(host.1,host.2,tree.id,contiguous,categorical.distance) %>%
				phyloscannerR:::get.pairwise.relationships.close.and.contiguous.cat() %>%
				group_by(host.1,host.2,tree.id) %>%
				inner_join(dwin, by=c('host.1','host.2','tree.id')) %>%
				ungroup() %>%
				select(-close.and.contiguous.cat,close.and.contiguous.cat)
	}
	if('close.and.adjacent'%in%relationship.types)
	{
		dwin <- dwin %>% 
				select(host.1,host.2,tree.id,adjacent,categorical.distance) %>%
				phyloscannerR:::get.pairwise.relationships.close.and.adjacent.cat() %>%
				group_by(host.1,host.2,tree.id) %>%
				inner_join(dwin, by=c('host.1','host.2','tree.id')) %>%
				ungroup() %>%
				select(-close.and.adjacent.cat,close.and.adjacent.cat)
	}
	if('close.and.contiguous.and.directed'%in%relationship.types)
	{
		dwin <- dwin %>% 
				select(host.1,host.2,tree.id,contiguous,categorical.distance,basic.classification) %>%
				phyloscannerR:::get.pairwise.relationships.close.and.contiguous.and.directed.cat() %>%
				group_by(host.1,host.2,tree.id) %>%
				inner_join(dwin, by=c('host.1','host.2','tree.id')) %>%
				ungroup() %>%
				select(-close.and.contiguous.and.directed.cat,close.and.contiguous.and.directed.cat)
	}
	if('close.and.adjacent.and.directed'%in%relationship.types)
	{
		dwin <- dwin %>% 
				select(host.1,host.2,tree.id,adjacent,categorical.distance,basic.classification) %>%
				phyloscannerR:::get.pairwise.relationships.close.and.adjacent.and.directed.cat() %>%
				group_by(host.1,host.2,tree.id) %>%
				inner_join(dwin, by=c('host.1','host.2','tree.id')) %>%
				ungroup() %>%
				select(-close.and.adjacent.and.directed.cat,close.and.adjacent.and.directed.cat)
	}
	if('close.and.contiguous.and.ancestry.cat'%in%relationship.types)
	{
		dwin <- dwin %>% 
				select(host.1,host.2,tree.id,contiguous,categorical.distance,basic.classification) %>%
				phyloscannerR:::get.pairwise.relationships.close.and.contiguous.and.ancestry.cat() %>%
				group_by(host.1,host.2,tree.id) %>%
				inner_join(dwin, by=c('host.1','host.2','tree.id')) %>%
				ungroup() %>%
				select(-close.and.contiguous.and.ancestry.cat,close.and.contiguous.and.ancestry.cat)
	}
	if('close.and.adjacent.and.ancestry.cat'%in%relationship.types)
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
	stopifnot( all(c('contiguous','adjacent','paths12','paths21')%in%colnames(all.classifications)) )
	# define "XXX_contiguous"
	tmp			<- all.classifications %>% filter(contiguous & adjacent)
	class 		<- list()
	class[[1]]	<- tmp %>% filter(grepl("anc|Anc",ancestry)) %>% mutate(basic.classification:= "anc_contiguous")
	class[[2]]	<- tmp %>% filter(grepl("desc|Desc",ancestry)) %>% mutate(basic.classification:= "desc_contiguous")
	class[[3]]	<- tmp %>% filter(ancestry=='complex') %>% mutate(basic.classification:= "intermingled_contiguous")
	class[[4]]	<- tmp %>% filter(ancestry=='none') %>% mutate(basic.classification:= "sibling_contiguous")
	# define "XXX_noncontiguous"
	tmp			<- all.classifications %>% filter(!contiguous & adjacent)
	class[[5]]	<- tmp %>% filter(grepl("anc|Anc",ancestry)) %>% mutate(basic.classification:= "anc_noncontiguous")
	class[[6]]	<- tmp %>% filter(grepl("desc|Desc",ancestry)) %>% mutate(basic.classification:= "desc_noncontiguous")
	class[[7]]	<- tmp %>% filter(ancestry=='complex') %>% mutate(basic.classification:= "intermingled_noncontiguous")
	class[[8]]	<- tmp %>% filter(ancestry=='none') %>% mutate(basic.classification:= "sibling_noncontiguous")
	# define "other"
	class[[9]]	<- all.classifications %>% filter(!adjacent) %>% mutate(basic.classification:= "other")
	class		<- bind_rows(class)
	stopifnot(nrow(all.classifications)==nrow(class))
	class
}

#' @keywords internal
#' @title Define basic topological relationships between two individuals by paths
#' @author Oliver Ratmann
get.pairwise.relationships.basic.by.paths <- function(all.classifications)
{
	stopifnot( all(c('contiguous','adjacent','paths12','paths21')%in%colnames(all.classifications)) )
	# define "XXX_contiguous"
	tmp			<- all.classifications %>% filter(contiguous & adjacent)
	class 		<- list()
	class[[1]]	<- tmp %>% filter(paths12>0 & paths21==0) %>% mutate(basic.classification:= "anc_contiguous")
	class[[2]]	<- tmp %>% filter(paths12==0 & paths21>0) %>% mutate(basic.classification:= "desc_contiguous")
	class[[3]]	<- tmp %>% filter(paths12>0 & paths21>0) %>% mutate(basic.classification:= "intermingled_contiguous")
	class[[4]]	<- tmp %>% filter(paths12==0 & paths21==0) %>% mutate(basic.classification:= "sibling_contiguous")
	# define "XXX_noncontiguous"
	tmp			<- all.classifications %>% filter(!contiguous & adjacent)
	class[[5]]	<- tmp %>% filter(paths12>0 & paths21==0) %>% mutate(basic.classification:= "anc_noncontiguous")
	class[[6]]	<- tmp %>% filter(paths12==0 & paths21>0) %>% mutate(basic.classification:= "desc_noncontiguous")
	class[[7]]	<- tmp %>% filter(paths12>0 & paths21>0) %>% mutate(basic.classification:= "intermingled_noncontiguous")
	class[[8]]	<- tmp %>% filter(paths12==0 & paths21==0) %>% mutate(basic.classification:= "sibling_noncontiguous")
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
