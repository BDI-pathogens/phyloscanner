#' @export
#' @author Oliver Ratmann
#' @import tidyverse 
#' @title Find pairs of individuals between whom linkage is not excluded phylogenetically
#' @param indir Full directory path to output of phyloscanner runs
#' @param batch.regex Regular expression that identifies the batch ID of multiple phyloscanner analyses. Default: '^ptyr([0-9]+)_.*'.
#' @param control List of control variables:
#' \itemize{
#' 		\item{\code{linked.group}} Phyloscanner classification used to identify pairs in networks. Default 'close.and.adjacent.cat'.
#' 		\item{\code{linked.no}} Phyloscanner classification type quantifying that pairs are not linked. Default 'not.close.or.nonadjacent'.
#' 		\item{\code{linked.yes}} Phyloscanner classification type quantifying that pairs are linked. Default 'close.and.adjacent'.
#' 		\item{\code{conf.cut}} Threshold on the proportion of deep-sequence phylogenies with distant/disconnected subgraphs above which pairs are considered phylogenetically unlinked. Default: 0.6
#' 		\item{\code{neff.cut}} Threshold on the minimum number of deep-sequence phylogenies with sufficient reads from two individuals to make any phylogenetic inferences. Default: 3.
#' }
#' @param verbose Flag to switch on/off verbose mode. Default: TRUE. 
#' @param dmeta Optional individual-level meta-data that is to be added to output. Can be NULL.
#' @return 
#' Three R objects are generated: 
#' \itemize{
#' 		\item{\code{network.pairs}} is a tibble that describes pairs of individuals between whom linkage is not excluded phylogenetically.
#' 		\item{\code{relationship.counts}} is a tibble that summarises the phylogenetic relationship counts for each pair.
#' 		\item{\code{windows}} is a tibble that describes the basic phyloscanner statistics (distance, adjacency, paths between subgraphs) in each deep-sequence phylogeny for each pair.
#' }
#' @description This function identifies pairs of individuals between whom linkage is not excluded phylogenetically in a large number of phyloscanner analyses, and provides detailed information on them.
#' @seealso \code{\link{phyloscanner.analyse.trees}}, \code{\link{cmd.phyloscanner.analyse.trees}}
#' @example inst/example/ex.transmission.networks.R
find.pairs.in.networks <- function(indir, batch.regex='^ptyr([0-9]+)_.*', control= list(linked.group='close.and.adjacent.cat',linked.no='not.close.or.nonadjacent',linked.yes='close.and.adjacent', conf.cut=0.6, neff.cut=3), dmeta=NULL, verbose=TRUE)
{	
	#	check if dmeta in right format
	if(!is.null(dmeta))
	{		
		if( !'tbl_df'%in%class(dmeta) )
			stop('"dmeta" must have class tbl_df (be a tibble)')
		if( !'ID'%in%colnames(dmeta) )
			stop('"dmeta" must have a column "ID"')
		if( class(dmeta$ID)!='character' )
			stop('"dmeta$ID" must be a character')		
	}	
	#	internal variables
	linked.group <- control$linked.group 
	linked.no <- control$linked.no 
	linked.yes <- control$linked.yes 
	conf.cut <- control$conf.cut
	neff.cut <- control$neff.cut	

	infiles	<- tibble(F=list.files(indir, pattern='workspace.rda', full.names=TRUE)) %>%
			mutate( PTY_RUN:= as.integer(gsub(batch.regex,'\\1',basename(F))) ) %>%
			arrange(PTY_RUN)	
	
	if(verbose) cat('\nFound phylogenetic relationship files, n=', nrow(infiles))
	if(verbose) cat('\nProcessing files...')
	
	dpls <- vector('list', nrow(infiles))
	dcs <- vector('list', nrow(infiles))
	dwins <- vector('list', nrow(infiles))
	for(i in seq_len(nrow(infiles)))
	{		
		if(verbose)	cat('\nReading ', infiles$F[i])
		tmp	<- load( infiles$F[i] )
		#	ensure we have multinomial output in workspace
		if(!'dc'%in%tmp)
			stop('Cannot find object "dc". Check that Analyze.trees was run with multinomial=TRUE.')
		if(!'dwin'%in%tmp)
			stop('Cannot find object "dwin". Check that Analyze.trees was run with multinomial=TRUE.')
		#	ensure IDs are characters
		if(!all(c( 	class( dc$host.1 )=='character',
						class( dc$host.2 )=='character',	
						class( dwin$host.1 )=='character',
						class( dwin$host.2 )=='character'
				)))
			stop('host.1 or host.2 not of character. This is unexpected, contact maintainer.')
		#	ensure host.1 < host.2		
		if( dc %>% filter( host.1 < host.2) %>% nrow() != dc %>% nrow() )
			stop('Not all host.1 < host.2 in "dc". This is unexpected, contact maintainer.')
		if( dwin %>% filter( host.1 < host.2) %>% nrow() != dwin %>% nrow() )
			stop('Not all host.1 < host.2 in "dwin". This is unexpected, contact maintainer.')
		
		#	get list of pairs that are not ph-unlinked
		dpl	<- dc %>% 
				filter(categorisation==linked.group & type==linked.no) %>% 
				filter(n.eff>=neff.cut & k.eff/n.eff < conf.cut) %>%
				select(host.1, host.2) %>%
				mutate(PTY_RUN:= infiles$PTY_RUN[i])
		dpl	<- dc %>% 
				filter(categorisation==linked.group & type==linked.yes) %>%
				right_join(dpl, by=c('host.1','host.2'))	
		#	gather all classifications counts for these pairs
		dwin	<- dpl %>%
				select(PTY_RUN, host.1, host.2) %>%
				left_join(dwin, by=c('host.1','host.2'))
		#	gather all classification counts for these pairs
		dc	<- dpl %>%
				select(PTY_RUN, host.1, host.2) %>%
				left_join(dc, by=c('host.1','host.2'))
		
		#	save
		dpls[[i]] <- dpl
		dwins[[i]] <- dwin
		dcs[[i]] <- dc		
	}		
	dpls <- do.call('rbind', dpls)
	dcs <- do.call('rbind', dcs)
	dwins <- do.call('rbind', dwins)
	#	the pairs may just be ambiguous, keep only those with some evidence for linkage	
	dpls <- dpls %>% filter(k.eff/n.eff > 0)	
	if(verbose) cat('\nFound (potentially duplicate) pairs between whom linkage is not excluded phylogenetically, n=', nrow(dpls))
	#	bring dwins into long format
	#	dwins <- dwins %>% gather('GROUP','TYPE', c('categorical.distance','basic.classification',ends_with('.cat')))
	
	#
	#	select analysis in which each pair has highest neff
	if(verbose) cat('\nIf pairs are in several batches, select batch with most deep-sequence phylogenies...')
	tmp <- dcs %>% 
			filter( categorisation==linked.group & type==linked.yes )	%>%
			group_by(host.1,host.2) %>%
			summarise(PTY_RUN= PTY_RUN[which.max(n.eff)]) %>%
			ungroup()
	dcs <- tmp %>% left_join(dcs, by=c('host.1','host.2','PTY_RUN'))
	dwins <- tmp %>% left_join(dwins, by=c('host.1','host.2','PTY_RUN'))
	dpls <- tmp %>% left_join(dpls, by=c('host.1','host.2','PTY_RUN'))	
	if(verbose) cat('\nLeft with pairs between whom linkage is not excluded phylogenetically, n=', nrow(dpls))
	
	#
	#	upper case col names
	setnames(dpls, colnames(dpls), toupper(gsub('\\.','_',gsub('host\\.','h',colnames(dpls)))))
	setnames(dcs, colnames(dcs), toupper(gsub('\\.','_',gsub('host\\.','h',colnames(dcs)))))
	setnames(dwins, colnames(dwins), toupper(gsub('\\.','_',gsub('host\\.','h',colnames(dwins)))))
	
	#
	#	add meta-data if provided
	if(!is.null(dmeta))
	{		
		if(verbose) cat('\nAdd meta-data...')
		tmp			<- unique(dmeta, by='ID')
		setnames(tmp, colnames(tmp), gsub('H1_ID','H1',paste0('H1_',colnames(tmp))))
		dpls <- dpls %>% left_join(tmp, by='H1')
		dcs <- dcs %>% left_join(tmp, by='H1')
		dwins <- dwins %>% left_join(tmp, by='H1')		
		setnames(tmp, colnames(tmp), gsub('H1','H2',colnames(tmp)))
		dpls <- dpls %>% left_join(tmp, by='H2')
		dcs <- dcs %>% left_join(tmp, by='H2')
		dwins <- dwins %>% left_join(tmp, by='H2')		
	}
	
	if(verbose) cat('\nDone. Found pairs, n=', nrow(dpls), '. Found relationship counts, n=', nrow(dcs), '. Found phyloscanner statistics, n=', nrow(dwins), '.')
	#	return
	list(network.pairs=dpls, relationship.counts=dcs, windows=dwins)
}


#' @export
#' @author Oliver Ratmann
#' @import tidyverse
#' @importFrom igraph graph.data.frame clusters
#' @title Find phylogenetic transmission networks and most likely transmission chain
#' @param dc Summary of phylogenetic relationship counts for each pair, stored as tibble.
#' @param control List of control variables:
#' \itemize{
#' 		\item{\code{linked.group}} Phyloscanner classification used to identify pairs in networks. Default 'close.and.adjacent.cat'.
#' 		\item{\code{linked.no}} Phyloscanner classification type quantifying that pairs are not linked. Default 'not.close.or.nonadjacent'.
#' 		\item{\code{linked.yes}} Phyloscanner classification type quantifying that pairs are linked. Default 'close.and.adjacent'.
#' 		\item{\code{neff.cut}} Threshold on the minimum number of deep-sequence phylogenies with sufficient reads from two individuals to make any phylogenetic inferences. Default: 3.
#' } 	
#' @param verbose Flag to switch on/off verbose mode. Default: TRUE. 
#' @return list of two R objects 
#' \itemize{
#' 		\item{\code{transmission.networks}} is a tibble that describes the edge list of pairs of individuals in a network, and corresponding phyloscanner scores
#' 		\item{\code{most.likely.transmission.chains}} is a tibble that describes the edge list of pairs of individuals in the most likely chain, and corresponding phyloscanner scores
#' }
#' See description.
#' @description 
#' This function computes transmission networks from phyloscanner output of a population-based deep sequence sample. 
#' A transmission network is defined as a set of individuals between whom phylogenetic linkage is not excluded.
#' Every individual in the network has at least one partner in the network between whom evidence for being phylogenetically unlinked is below a threshold.
#' These pairs of individuals are identified with a separate function, \code{\link{find.pairs.in.networks}}.  
#' Due to the nature of the deep-sequence
#' data, there are up to three edges between pairs of individuals, giving the strength of evidence
#' of spread in each direction (two possibilities) and the strength of evidence for phylogenetic linkage with 
#' the direction remaining unclear. Some of these pairs have limited evidence for phylogenetic linkage. The networks
#' are best interpreted as partially sampled transmission chain. 
#' 
#' This function also finds the most likely transmission 
#' chain among all chains spanning the nodes in a specified transmission network.
#' The transmission network consists of at most three edges between a set of
#' individuals (directed edge in either direction, and undirected edge). Chains 
#' are defined as spanning graphs through the set of nodes, without loops and with
#' in-degree equal to one for all nodes, except the start node. Each directed edge 
#' in a chain has a weight, which corresponds to the phylogenetic evidence of transmission
#' in this direction. It is set to the phyloscanner score for transmission in this 
#' direction plus half the phyloscanner score of the undirected edge, for transmission 
#' with direction unclear. The probability of the entire chain is given by the product
#' of the phyloscanner scores along each edge in the chain.
#'    
#' @seealso \code{\link{find.pairs.in.networks}}, \code{\link{plot.network}}, \code{\link{plot.chain}}
#' @example inst/example/ex.transmission.networks.R
find.networks<- function(dc, control= list(linked.group='close.and.adjacent.cat',linked.no='not.close.or.nonadjacent',linked.yes='close.and.adjacent', neff.cut=3), verbose=TRUE)
{	
	#	internal constants
	linked.group 	<- control$linked.group 
	linked.no 		<- control$linked.no 
	linked.yes 		<- control$linked.yes 	
	neff.cut 		<- control$neff.cut		
	scores.group	<- 'close.and.adjacent.and.ancestry.cat'
	scores.nolink	<- 'not.close.or.nonadjacent'
	scores.ambig	<- 'complex.or.no.ancestry' 
	dir.group		<- "close.and.adjacent.and.directed.cat"
	
	#
	#	construct tri-edge transmission network	
	dnet <- dc %>% 
			filter( CATEGORISATION==linked.group &  TYPE==linked.yes & N_EFF>neff.cut) %>%
			select(-c(CATEGORICAL_DISTANCE, TYPE, K, K_EFF, N, N_EFF, SCORE))
	#	define potential transmission network membership
	if(verbose) cat('\nReconstruct transmission networks among linked pairs, n=',nrow(dnet))
	tmp <- dnet %>% select(H1, H2)			
	tmp <- igraph:::graph.data.frame(tmp, directed=FALSE, vertices=NULL)
	rtc <- tibble(ID=V(tmp)$name, CLU=clusters(tmp, mode="weak")$membership)
	rtc <- rtc %>% 
			group_by(CLU) %>% 
			summarise( CLU_SIZE:=length(ID) ) %>% 
			arrange(desc(CLU_SIZE)) %>%
			ungroup() %>%
			mutate( IDCLU:=seq_along(CLU_SIZE) ) %>%
			inner_join(rtc, by='CLU') %>%
			select(-CLU)	
	#	add info on edges: network membership
	setnames(rtc, c('ID'), c('H1'))
	dnet <- dnet %>% inner_join(rtc, by='H1')
	rtc <- rtc %>% select(-CLU_SIZE)
	setnames(rtc, c('H1'), c('H2'))
	dnet <- dnet %>% inner_join(rtc, by=c('H2','IDCLU'))
	#	add info on edges: direction 12, direction 21, direction ambiguous, unlinked
	#	add this info in new rows
	tmp <- dc %>% 
			filter(CATEGORISATION==scores.group) %>%
			select(PTY_RUN, H1, H2, TYPE, SCORE)
	dnet <- dnet %>% inner_join(tmp, by=c('H1','H2','PTY_RUN'))		
	if(verbose) cat(	'\nFound transmission networks, n=', length(unique(dnet$IDCLU)),'. ', 
				'Number of links (either direction and ambiguous)=', dnet %>% filter(TYPE!=scores.nolink) %>% distinct(H1,H2) %>% nrow(), '.',
				'Number of individuals=', length(unique(c(dnet$H1, dnet$H2))),'.')
	#	generate most likely transmission chains
	if(verbose) cat('\nReconstruct most likely transmission chains...')	
	tmp <- dnet %>% 
			filter(TYPE!=scores.nolink) %>% 
			select(PTY_RUN,H1,H2,IDCLU,CLU_SIZE,TYPE,SCORE)
	dchain <- find.most.likely.chains.RBGLedmonds(tmp)
	
	#	merge pw linkage scores to the ML chain representation
	#	rationale: this describes prob of linkage. here, any 'inconsistent direction' is still considered as prob for linkage
	tmp <- dc %>% 
			filter( CATEGORISATION==linked.group & TYPE==linked.yes ) %>%
			select(H1,H2,PTY_RUN,SCORE) %>%
			rename(SCORE_LINKED=SCORE)
	dchain <- dchain %>% inner_join(tmp, by=c('H1','H2','PTY_RUN'))
	#	merge pw direction scores to the ML chain representation
	#	this is considering in denominator 12 + 21 before reducing probs to achieve self-consistency
	#	rationale: decide on evidence for direction based on comparing only the flows in either direction, 12 vs 21
	tmp <- dc %>% 
			filter( CATEGORISATION==dir.group ) %>%
			select(H1,H2,PTY_RUN,TYPE,SCORE) %>%
			mutate(TYPE:= paste0('SCORE_DIR_',TYPE)) %>%
			spread(TYPE,SCORE)
	dchain <- dchain %>% inner_join(tmp, by=c('H1','H2','PTY_RUN'))
	#	merge pw network scores to the ML chain representation
	#	this is considering in denominator 12 + 21 + unclear reducing probs to achieve self-consistency
	#	same as MX_PROB_12, MX_PROB_21, after the final step below that sets one of the two probs to zero
	tmp <- dc %>% 
			filter( CATEGORISATION==scores.group & TYPE!=scores.nolink) %>%
			select(H1,H2,PTY_RUN,TYPE,SCORE) %>%
			mutate(TYPE:= replace(TYPE, TYPE==scores.ambig, 'AMB')) %>%
			mutate(TYPE:= paste0('SCORE_NW_',TYPE)) %>%
			spread(TYPE,SCORE)
	dchain <- dchain %>% inner_join(tmp, by=c('H1','H2','PTY_RUN'))
	#	ensure scores are compatible with self-consistency in ML chain
	dchain <- dchain %>%
			filter(LINK_12==1 | LINK_21==1) %>%
			mutate(	SCORE_NW_12= replace(SCORE_NW_12, LINK_12==0 & LINK_21==1 & SCORE_NW_12>SCORE_NW_21, 0),
					SCORE_NW_21= replace(SCORE_NW_21, LINK_12==1 & LINK_21==0 & SCORE_NW_21>SCORE_NW_12, 0),
					SCORE_DIR_12= replace(SCORE_DIR_12, LINK_12==0 & LINK_21==1 & SCORE_DIR_12>SCORE_DIR_21, 0),
					SCORE_DIR_21= replace(SCORE_DIR_21, LINK_12==1 & LINK_21==0 & SCORE_DIR_21>SCORE_DIR_12, 0),
			)
	#	copy meta data to ML chain
	tmp <- dnet %>%  
			select(-c(CLU_SIZE, IDCLU, TYPE, SCORE)) %>% 
			distinct()
	dchain <- dchain %>% inner_join(tmp, by=c('H1','H2','PTY_RUN'))
	#	verbose
	if(verbose) cat('\nFound most likely transmission chains, n=', dchain %>% select(IDCLU) %>% distinct() %>% nrow(), '.',
					'Number of links=', nrow(dchain), '.',
					'Number of individuals=', length(unique(c(dchain$H1, dchain$H2))),'.')
	if(verbose) cat('\nDone.')
	# return
	list(transmission.networks=dnet, most.likely.transmission.chains=dchain)
}

#' @keywords internal
#' @author Oliver Ratmann
#' @importFrom igraph graph.data.frame igraph.to.graphNEL
#' @importFrom RBGL edmondsOptimumBranching
#' @title Find most likely transmission chains
#' @description This function finds the most likely transmission 
#' chain among all chains spanning the nodes in a specified transmission network.
#' The transmission network consists of at most three edges between a set of
#' individuals (directed edge in either direction, and undirected edge). Chains 
#' are defined as spanning graphs through the set of nodes, without loops and with
#' in-degree equal to one for all nodes, except the start node. Each directed edge 
#' in a chain has a weight, which corresponds to the phylogenetic evidence of transmission
#' in this direction. It is set to the phyloscanner score for transmission in this 
#' direction plus half the phyloscanner score of the undirected edge, for transmission 
#' with direction unclear. The probability of the entire chain is given by the product
#' of the phyloscanner scores along each edge in the chain.   
#' Edmonds algorithm is used to solve for the most probable chain.
#' @param rtnn tibble encoding the three edges of a phyloscanner transmission network. Must have columns 'H1','H2','IDCLU','TYPE','SCORE','K_EFF'.   
#' @return tibble encoding the most likely chain. Has columns 'H1','H2','IDCLU', 'LINK_12', 'LINK_21' (either 1 or 0 for a link in the corresponding direction), and 'MX_PROB_12', 'MX_PROB_21' (associated posterior probabilities)
#' @seealso \code{\link{find.networks}} 
find.most.likely.chains.RBGLedmonds<- function(rtnn, verbose=0)
{
	stopifnot(c('PTY_RUN','H1','H2','IDCLU','CLU_SIZE','TYPE','SCORE')%in%colnames(rtnn))
	if( length(setdiff(c('complex.or.no.ancestry','21','12'), rtnn %>% distinct(TYPE) %>% pull(TYPE)) ))
		stop('Unexpected classification types. Contact maintainer.')	
	#	define weights to convert tri-edges to bi-edges 
	#	and create bi-edges
	rtm <- rtnn %>% 
			mutate(	ID1_IN_WEIGHT= case_when(	TYPE=='12' ~ 0,
							TYPE=='complex.or.no.ancestry' ~ 0.5,
							TYPE=='21' ~ 1),
					ID2_IN_WEIGHT= case_when(	TYPE=='21' ~ 0,
							TYPE=='complex.or.no.ancestry' ~ 0.5,
							TYPE=='12' ~ 1)															
			) %>%
			group_by( PTY_RUN,IDCLU,CLU_SIZE,H1,H2 ) %>%
			summarise( 	PROB_21= sum(SCORE*ID1_IN_WEIGHT), 					
						PROB_12= sum(SCORE*ID2_IN_WEIGHT)
						) %>%
			ungroup()
	#	handle networks of size 2 - this is easy	
	ans <- rtm %>% 
			filter(CLU_SIZE==2) %>%
			mutate( LINK_12 = ifelse(PROB_12>PROB_21, 1L, 0L),
					LINK_21 = ifelse(PROB_21>PROB_12, 1L, 0L) )	
	#	handle networks of size >2 - use Edmonds algorithm 	
	rtm <- rtm %>% 
			filter(CLU_SIZE>2) %>%
			group_by(IDCLU) %>%
			group_split()
	for(i in seq_along(rtm))
	{
		#	
		#i	<- 13 
		if(verbose)
			cat('\nIDCLU ',i)
		#	convert to weighted edge list
		edgelist <- rtm[[i]] %>% 
				mutate(weight:=PROB_21) %>% 
				rename(H1:= H2, H2:= H1) %>%
				select(H1,H2,weight) %>%
				bind_rows(tmp) %>%
				filter(weight>0) 
		#	maximisation is over sum of weights, so transform to log prob
		#	since edmonds cannot handle negative values or zeros, add some constant
		edgelist <- edgelist %>% mutate(weight:= log(weight) - min(log(weight)) + 1 )
		#	run Edmonds 
		g <- graph.data.frame(edgelist)	
		g2 <- igraph.to.graphNEL(g)
		g3 <- edmondsOptimumBranching(g2)
		edgelist <- as_tibble(t(g3$edgeList)) %>%
				rename(H1:=from, H2:= to) %>%
				mutate(	LINK_12:= 1L,
						LINK_21:= 0L)
		edgelist <- edgelist %>% 
				rename(H1:= H2, H2:=H1, LINK_12:=LINK_21, LINK_21:=LINK_12) %>%
				bind_rows(edgelist) 
		rtm[[i]] <- rtm[[i]] %>% inner_join(edgelist, by=c('H1','H2'))						
	}
	rtm		<- do.call('rbind',rtm)
	ans		<- rbind(ans, rtm)
	ans		<- ans %>% 
			select(IDCLU, CLU_SIZE, PTY_RUN, H1, H2, LINK_12, LINK_21)
			#mutate(	PROB_21= replace(PROB_21, LINK_12==1, 0),
			#		PROB_12= replace(PROB_12, LINK_21==1, 0)
			#		) %>%
			#rename(MX_PROB_21=PROB_21, MX_PROB_12=PROB_12, MX_KEFF_21=KEFF_21, MX_KEFF_12=KEFF_12 )
	ans
}