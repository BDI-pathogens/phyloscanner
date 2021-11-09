#' @export
#' @author Oliver Ratmann
#' @import tidyverse
#' @title Find pairs of individuals between whom linkage is not excluded phylogenetically
#' @param dwin A \code{data.frame} describing pairwise relationships between the hosts in each tree; normally output of \code{classify.pairwise.relationships}
#' @param dc A \code{data.frame} summarising pairwise relationships between the hosts across all trees; normally output of \code{count.pairwise.relationships}
#' @param control List of control variables:
#' \itemize{
#'         \item{\code{linked.group}} Phyloscanner classification used to identify pairs in networks. Default 'close.and.adjacent.cat'.
#'         \item{\code{linked.no}} Phyloscanner classification type quantifying that pairs are not linked. Default 'not.close.or.nonadjacent'.
#'         \item{\code{linked.yes}} Phyloscanner classification type quantifying that pairs are linked. Default 'close.and.adjacent'.
#'         \item{\code{conf.cut}} Threshold on the proportion of deep-sequence phylogenies with distant/disconnected subgraphs above which pairs are considered phylogenetically unlinked. Default: 0.6
#'         \item{\code{neff.cut}} Threshold on the minimum number of deep-sequence phylogenies with sufficient reads from two individuals to make any phylogenetic inferences. Default: 3.
#' }
#' @param verbose Flag to switch on/off verbose mode. Default: TRUE.
#' @param dmeta Optional individual-level meta-data that is to be added to output. Can be NULL.
#' @return
#' Three R objects are generated:
#' \itemize{
#'         \item{\code{network.pairs}} is a tibble that describes pairs of individuals between whom linkage is not excluded phylogenetically.
#'         \item{\code{relationship.counts}} is a tibble that summarises the phylogenetic relationship counts for each pair.
#'         \item{\code{windows}} is a tibble that describes the basic phyloscanner statistics (distance, adjacency, paths between subgraphs) in each deep-sequence phylogeny for each pair.
#' }
#' @description This function identifies pairs of individuals between whom linkage is not excluded phylogenetically in a large number of phyloscanner analyses, and provides detailed information on them.
#' @seealso \code{\link{phyloscanner.analyse.trees}}, \code{\link{cmd.phyloscanner.analyse.trees}}
#' @example inst/example/ex.transmission.networks.R
find.pairs.in.networks <- function(dwin, dc, control= list(linked.group='close.and.adjacent.cat',linked.no='not.close.or.nonadjacent',linked.yes='close.and.adjacent', conf.cut=0.6, neff.cut=3), dmeta=NULL, verbose=TRUE)
{
    #    check if dmeta in right format
    if(!is.null(dmeta))
    {
        if( !'tbl_df'%in%class(dmeta) )
            stop('"dmeta" must have class tbl_df (be a tibble)')
        if( !'ID'%in%colnames(dmeta) )
            stop('"dmeta" must have a column "ID"')
        if( class(dmeta$ID)!='character' )
            stop('"dmeta$ID" must be a character')
    }
    #    internal variables
    linked.group <- control$linked.group
    linked.no <- control$linked.no
    linked.yes <- control$linked.yes
    conf.cut <- control$conf.cut
    neff.cut <- control$neff.cut
    
    #    ensure IDs are characters
    if(!all(c(     class( dc$host.1 )=='character',
                    class( dc$host.2 )=='character',
                    class( dwin$host.1 )=='character',
                    class( dwin$host.2 )=='character'
            )))
        stop('host.1 or host.2 not of character. This is unexpected, contact maintainer.')
    #    ensure host.1 < host.2
    if( dc %>% filter( host.1 < host.2) %>% nrow() != dc %>% nrow() )
        stop('Not all host.1 < host.2 in "dc". This is unexpected, contact maintainer.')
    if( dwin %>% filter( host.1 < host.2) %>% nrow() != dwin %>% nrow() )
        stop('Not all host.1 < host.2 in "dwin". This is unexpected, contact maintainer.')
    
    #    get list of pairs that are not ph-unlinked
    dpl    <- dc %>%
            filter(categorisation==linked.group & type==linked.no) %>%
            filter(n.eff>=neff.cut & k.eff/n.eff < conf.cut) %>%
            select(host.1, host.2)
    # %>%        mutate(PTY_RUN:= infiles$PTY_RUN[i])
    dpl    <- dc %>%
            filter(categorisation==linked.group & type==linked.yes) %>%
            right_join(dpl, by=c('host.1','host.2'))
    #    the pairs may just be ambiguous, keep only those with some evidence for linkage
    dpl <- dpl %>% filter(k.eff/n.eff > 0)
    #    gather all classifications counts for these pairs
    dwin    <- dpl %>%
           # select(PTY_RUN, host.1, host.2) %>%
            select(host.1, host.2) %>%
            left_join(dwin, by=c('host.1','host.2'))
    #    gather all classification counts for these pairs
    dc    <- dpl %>%
            # select(PTY_RUN, host.1, host.2) %>%
            select(host.1, host.2) %>%
            left_join(dc, by=c('host.1','host.2'))
    
    #
    #    upper case col names
    setnames(dpl, colnames(dpl), toupper(gsub('\\.','_',gsub('host\\.','h',colnames(dpl)))))
    setnames(dc, colnames(dc), toupper(gsub('\\.','_',gsub('host\\.','h',colnames(dc)))))
    setnames(dwin, colnames(dwin), toupper(gsub('\\.','_',gsub('host\\.','h',colnames(dwin)))))
    
    #
    #    add meta-data if provided
    if(!is.null(dmeta))
    {
        if(verbose) cat('\nAdd meta-data...')
        tmp            <- unique(dmeta, by='ID')
        setnames(tmp, colnames(tmp), gsub('H1_ID','H1',paste0('H1_',colnames(tmp))))
        dpl <- dpl %>% left_join(tmp, by='H1')
        dc <- dc %>% left_join(tmp, by='H1')
        dwin <- dwin %>% left_join(tmp, by='H1')
        setnames(tmp, colnames(tmp), gsub('H1','H2',colnames(tmp)))
        dpl <- dpl %>% left_join(tmp, by='H2')
        dc <- dc %>% left_join(tmp, by='H2')
        dwin <- dwin %>% left_join(tmp, by='H2')
    }
    
    if(verbose) cat('\nDone. Found pairs, n=', nrow(dpl), '. Found relationship counts, n=', nrow(dc), '. Found phyloscanner statistics, n=', nrow(dwin), '.')
    #    return
    list(network.pairs=dpl, relationship.counts=dc, windows=dwin)
}


#' @export
#' @author Oliver Ratmann
#' @import tidyverse glue
#' @importFrom igraph graph.data.frame clusters
#' @title Find phylogenetic transmission networks and most likely transmission chain
#' @param dc Summary of phylogenetic relationship counts for each pair, stored as tibble.
#' @param control List of control variables:
#' \itemize{
#'         \item{\code{linked.group}} Phyloscanner classification used to identify pairs in networks. Default 'close.and.adjacent.cat'.
#'         \item{\code{linked.no}} Phyloscanner classification type quantifying that pairs are not linked. Default 'not.close.or.nonadjacent'.
#'         \item{\code{linked.yes}} Phyloscanner classification type quantifying that pairs are linked. Default 'close.and.adjacent'.
#'         \item{\code{neff.cut}} Threshold on the minimum number of deep-sequence phylogenies with sufficient reads from two individuals to make any phylogenetic inferences. Default: 3.
#'         \item{\code{weight.complex.or.no.ancestry}} Weight given to score complex.or.no.ancestry. Default: 50%.
#' }
#' @param verbose Flag to switch on/off verbose mode. Default: TRUE.
#' @return list of two R objects
#' \itemize{
#'         \item{\code{transmission.networks}} is a tibble that describes the edge list of pairs of individuals in a network, and corresponding phyloscanner scores
#'         \item{\code{most.likely.transmission.chains}} is a tibble that describes the edge list of pairs of individuals in the most likely chain, and corresponding phyloscanner scores
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
find.networks<- function(dc, control= list(linked.group='close.and.adjacent.cat',linked.no='not.close.or.nonadjacent',linked.yes='close.and.adjacent', neff.cut=3, weight.complex.or.no.ancestry=0.5), verbose=TRUE)
{
    #    internal constants
    linked.group     <- control$linked.group
    linked.no         <- control$linked.no
    linked.yes         <- control$linked.yes
    neff.cut         <- control$neff.cut
    weight.complex.or.no.ancestry <- control$weight.complex.or.no.ancestry
    scores.group    <- 'close.and.adjacent.and.ancestry.cat'
    scores.nolink    <- 'not.close.or.nonadjacent'
    scores.ambig    <- 'complex.or.no.ancestry'
    dir.group        <- "close.and.adjacent.and.directed.cat"
    
    #
    #    construct tri-edge transmission network
    dnet <- dc %>%
            filter( CATEGORISATION==linked.group &  TYPE==linked.yes & N_EFF>neff.cut) %>%
            select(-c(CATEGORICAL_DISTANCE, TYPE, K, K_EFF, N, N_EFF, SCORE))
    #    define potential transmission network membership
    if(verbose) cat(glue('\nReconstructing {nrow(dnet)} transmission networks among linked pairs'))
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
    #    add info on edges: network membership
    rtc <- rtc %>% rename(H1:=ID)
    dnet <- dnet %>% inner_join(rtc, by='H1')
    rtc <- rtc %>% select(-CLU_SIZE)
    rtc <- rtc %>% rename(H2:=H1)
    dnet <- dnet %>% inner_join(rtc, by=c('H2','IDCLU'))
    #    add info on edges: direction 12, direction 21, direction ambiguous, unlinked
    #    add this info in new rows
    tmp <- dc %>%
            filter(CATEGORISATION==scores.group) %>%
            select(PTY_RUN, H1, H2, TYPE, SCORE)
    dnet <- dnet %>% inner_join(tmp, by=c('H1','H2','PTY_RUN'))
    if(verbose) cat(glue('\nFound {length(unique(dnet$IDCLU))} transmission networks with
                         {dnet %>% filter(TYPE!=scores.nolink) %>% distinct(H1,H2) %>% nrow()}
                         links (in either direction or undirected) and {length(unique(c(dnet$H1, dnet$H2)))}
                         individuals'))
    #    generate most likely transmission chains
    if(verbose) cat('\nReconstructing most likely transmission chains...')
    tmp <- dnet %>%
            filter(TYPE!=scores.nolink) %>%
            select(PTY_RUN,H1,H2,IDCLU,CLU_SIZE,TYPE,SCORE)
    dchain <- find.most.likely.chains.RBGLedmonds(tmp, weight.complex.or.no.ancestry=weight.complex.or.no.ancestry)
    
    #    merge pw linkage scores to the ML chain representation
    #    rationale: this describes prob of linkage. here, any 'inconsistent direction' is still considered as prob for linkage
    tmp <- dc %>%
            filter( CATEGORISATION==linked.group & TYPE==linked.yes ) %>%
            select(H1,H2,PTY_RUN,SCORE) %>%
            rename(SCORE_LINKED=SCORE)
    dchain <- dchain %>% inner_join(tmp, by=c('H1','H2','PTY_RUN'))
    #    merge pw direction scores to the ML chain representation
    #    this is considering in denominator 12 + 21 before reducing probs to achieve self-consistency
    #    rationale: decide on evidence for direction based on comparing only the flows in either direction, 12 vs 21
    tmp <- dc %>%
            filter( CATEGORISATION==dir.group ) %>%
            select(H1,H2,PTY_RUN,TYPE,SCORE) %>%
            mutate(TYPE:= paste0('SCORE_DIR_',TYPE)) %>%
            spread(TYPE,SCORE)
    dchain <- dchain %>% left_join(tmp, by=c('H1','H2','PTY_RUN'))
    #    merge pw network scores to the ML chain representation
    #    this is considering in denominator 12 + 21 + unclear reducing probs to achieve self-consistency
    #    same as MX_PROB_12, MX_PROB_21, after the final step below that sets one of the two probs to zero
    tmp <- dc %>%
            filter( CATEGORISATION==scores.group & TYPE!=scores.nolink) %>%
            select(H1,H2,PTY_RUN,TYPE,SCORE) %>%
            mutate(TYPE:= replace(TYPE, TYPE==scores.ambig, 'AMB')) %>%
            mutate(TYPE:= paste0('SCORE_NW_',TYPE)) %>%
            spread(TYPE,SCORE)
    dchain <- dchain %>% inner_join(tmp, by=c('H1','H2','PTY_RUN'))
    #    fill in NA
    dchain <- dchain %>%
            mutate_at(c('SCORE_DIR_12','SCORE_DIR_21','SCORE_NW_12','SCORE_NW_21'), function(x) replace(x, is.na(x), 0))
     #    ensure scores are compatible with self-consistency in ML chain
    dchain <- dchain %>%
            filter(LINK_12==1 | LINK_21==1) %>%
            mutate(    SCORE_NW_12= replace(SCORE_NW_12, LINK_12==0 & LINK_21==1 & SCORE_NW_12>SCORE_NW_21, 0),
                    SCORE_NW_21= replace(SCORE_NW_21, LINK_12==1 & LINK_21==0 & SCORE_NW_21>SCORE_NW_12, 0),
                    SCORE_DIR_12= replace(SCORE_DIR_12, LINK_12==0 & LINK_21==1 & SCORE_DIR_12>SCORE_DIR_21, 0),
                    SCORE_DIR_21= replace(SCORE_DIR_21, LINK_12==1 & LINK_21==0 & SCORE_DIR_21>SCORE_DIR_12, 0),
            )
    #    copy meta data to ML chain
    tmp <- dnet %>%
            select(-c(CLU_SIZE, IDCLU, TYPE, SCORE)) %>%
            distinct()
    dchain <- dchain %>% inner_join(tmp, by=c('H1','H2','PTY_RUN'))
    #    verbose
    if(verbose) cat(glue("\nIdentified {dchain %>% select(IDCLU) %>% distinct() %>% nrow()} most likely transmission chains
                          with {nrow(dchain)} links and {length(unique(c(dchain$H1, dchain$H2)))} individuals"))
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
#' @param weight.complex.or.no.ancestry weight given to score complex.or.no.ancestry, default: 50\%
#' @return tibble encoding the most likely chain. Has columns 'H1','H2','IDCLU', 'LINK_12', 'LINK_21' (either 1 or 0 for a link in the corresponding direction), and 'MX_PROB_12', 'MX_PROB_21' (associated posterior probabilities)
#' @seealso \code{\link{find.networks}}
find.most.likely.chains.RBGLedmonds<- function(rtnn, weight.complex.or.no.ancestry=0.5, verbose=0)
{
    stopifnot(c('PTY_RUN','H1','H2','IDCLU','CLU_SIZE','TYPE','SCORE')%in%colnames(rtnn))
    if( length(setdiff(c('complex.or.no.ancestry','21','12'), rtnn %>% distinct(TYPE) %>% pull(TYPE)) ))
        stop('Unexpected classification types. Contact maintainer.')
    #    define weights to convert tri-edges to bi-edges
    #    and create bi-edges
    rtm <- rtnn %>%
            mutate(    ID1_IN_WEIGHT= case_when(    TYPE=='12' ~ 0,
                            TYPE=='complex.or.no.ancestry' ~ weight.complex.or.no.ancestry,
                            TYPE=='21' ~ 1),
                    ID2_IN_WEIGHT= case_when(    TYPE=='21' ~ 0,
                            TYPE=='complex.or.no.ancestry' ~ weight.complex.or.no.ancestry,
                            TYPE=='12' ~ 1)
            ) %>%
            group_by( PTY_RUN,IDCLU,CLU_SIZE,H1,H2 ) %>%
            summarise(     PROB_21= sum(SCORE*ID1_IN_WEIGHT),
                        PROB_12= sum(SCORE*ID2_IN_WEIGHT)
                        ) %>%
            ungroup()
    #    handle networks of size 2 - this is easy
    ans <- rtm %>%
            filter(CLU_SIZE==2) %>%
            transform(TIES_RND_BRK= runif(length(CLU_SIZE))) %>%
            mutate( LINK_12 = ifelse((PROB_12>PROB_21) | (PROB_12==PROB_21 & TIES_RND_BRK>=0.5), 1L, 0L),
                    LINK_21 = ifelse((PROB_21>PROB_12) | (PROB_12==PROB_21 & TIES_RND_BRK<0.5), 1L, 0L) ) %>%
            select(-TIES_RND_BRK)
    #    handle networks of size >2 - use Edmonds algorithm
    rtm <- rtm %>%
            filter(CLU_SIZE>2) %>%
            group_by(IDCLU) %>%
            group_split()
    rtm2 <- list()
    for(i in seq_along(rtm))
    {
        #
        #i    <- 13
        if(verbose)
            cat('\nIDCLU ',i)
        #    convert to weighted edge list
        tmp <- rtm[[i]] %>%
                mutate(weight:=PROB_12) %>%
                select(H1,H2,weight)
        edgelist <- rtm[[i]] %>%
                mutate(weight:=PROB_21) %>%
                rename(H1:= H2, H2:= H1) %>%
                select(H1,H2,weight) %>%
                bind_rows(tmp) %>%
                filter(weight>0)
        #    maximisation is over sum of weights, so transform to log prob
        #    since edmonds cannot handle negative values or zeros, add some constant
        edgelist <- edgelist %>% mutate(weight:= log(weight) - min(log(weight)) + 1 )
        #    run Edmonds
        g <- graph.data.frame(edgelist)
        g2 <- igraph.to.graphNEL(g)
        g3 <- edmondsOptimumBranching(g2)
        edgelist <- as_tibble(t(g3$edgeList)) %>%
                rename(H1:=from, H2:= to) %>%
                mutate(    LINK_12:= 1L,
                        LINK_21:= 0L)
        edgelist <- edgelist %>%
                rename(H1:= H2, H2:=H1, LINK_12:=LINK_21, LINK_21:=LINK_12) %>%
                bind_rows(edgelist)
        # rtm[[i]] <- rtm[[i]] %>% inner_join(edgelist, by=c('H1','H2'))
        rtm2[[i]] <- rtm[[i]] %>% inner_join(edgelist, by=c('H1','H2'))
    }
    rtm2        <- do.call('rbind',rtm2)
    ans        <- rbind(ans, rtm2)
    ans        <- ans %>%
            select(IDCLU, CLU_SIZE, PTY_RUN, H1, H2, LINK_12, LINK_21)
            #mutate(    PROB_21= replace(PROB_21, LINK_12==1, 0),
            #        PROB_12= replace(PROB_12, LINK_21==1, 0)
            #        ) %>%
            #rename(MX_PROB_21=PROB_21, MX_PROB_12=PROB_12, MX_KEFF_21=KEFF_21, MX_KEFF_12=KEFF_12 )
    ans
}

