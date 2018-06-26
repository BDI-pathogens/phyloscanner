

#' @title Calculate parameters of the posterior density for pairwise host relationships
#' @export
#' @param phyloscanner.trees A list of class \code{phyloscanner.trees} produced by \code{phyloscanner.analyse.trees}.
#' @param close.threshold The (potentially normalised) patristic threshold used to determine if two patients' subgraphs are "close"
#' @param tip.regex The regular expression used to identify host IDs in tip names
#' @param allow.mt If FALSE, directionality is only inferred between pairs of hosts where a single clade from one host is nested in one from the other; this is more conservative.
#' @param min.reads The minimum number of reads from a host in a window needed in order for that window to count in determining relationships involving that patient
#' @param min.tips The minimum number of tips from a host in a window needed in order for that window to count in determining relationships involving that patient
#' @param distant.threshold If present, a second distance threshold determines hosts that are "distant" from each other, with those lying between \code{close.threshold} and \code{dist.threshold} classed as "intermediate". The default is the same as \code{close.threshold}, so the intermediate class does not exist.
#' @param verbose Verbose output
#' @return A list with two items: \code{dwin} giving information on the genome windows for each pair of hosts, and \code{rplkl} giving information on phylogenetic relationships between each pair of hosts.

multinomial.calculations <- function(ptrees, 
                                     close.threshold,
                                     prior.keff = 3,
                                     prior.neff = 4,
                                     prior.calibrated.prob = 0.66,
                                     tip.regex = "^(.*)_read_([0-9]+)_count_([0-9]+)$", 
                                     allow.mt = F, 
                                     min.reads = 0, 
                                     min.tips = 0, 
                                     distant.threshold = close.threshold,
                                     relationship.types	= c('proximity.3.way',
                                                            'any.ancestry',
                                                            'close.x.contiguous',
                                                            'close.and.contiguous',
                                                            'close.and.contiguous.and.directed',
                                                            'close.and.adjacent.and.directed',
                                                            'close.and.contiguous.and.ancestry.cat',
                                                            'close.and.adjacent.and.ancestry.cat',
                                                            'adjacent.and.proximity.cat'),
                                     verbose=F) {
  
  if(!attr(ptrees, "readable.coords")){
    stop("No window cooardinates detected in this tree set, cannot do multinomial calculations.")
  }
  
  all.classifications <- merge.classifications(ptrees, allow.mt, verbose)
  
  host.tips.and.reads <- map(ptrees, function(x) get.tip.and.read.counts(x, all.hosts.from.trees(ptrees), tip.regex, attr(ptrees, 'has.read.counts'), verbose))
  host.tips.and.reads <- bind_rows(host.tips.and.reads)
  
  all.classifications <- all.classifications %>% inner_join(host.tips.and.reads, by=c("host.1"="host.id", "tree.id")) %>% rename(tips.1 = tips, reads.1=reads)
  all.classifications <- all.classifications %>% inner_join(host.tips.and.reads, by=c("host.2"="host.id", "tree.id")) %>% rename(tips.2 = tips, reads.2=reads)

  all.classifications <- all.classifications %>% mutate(window.start = map_int(tree.id, function(x) as.integer(gsub('[^0-9]*([0-9]+)_to_([0-9]+).*','\\1', x))))
  all.classifications <- all.classifications %>% mutate(window.end = map_int(tree.id, function(x) as.integer(gsub('[^0-9]*([0-9]+)_to_([0-9]+).*','\\2', x))))
  
  all.classifications <- all.classifications %>% mutate(categorical.distance = map_chr(patristic.distance, function(x){
    if(x < close.threshold){
      "close"
    } else if(x >= distant.threshold){
      "distant"
    } else {
      "intermediate"
    }
  }))
  
  if(verbose) cat('\nReducing transmission window stats to windows with at least',min.reads,'reads and at least',min.tips,'tips ...')
  all.classifications	<- filter(all.classifications, reads.1>=min.reads & reads.2>=min.reads & tips.1>=min.tips & tips.2>=min.tips)
  if(verbose) cat('\nTotal number of windows with transmission assignments is ',nrow(all.classifications),'.', sep="")		
  
  if(verbose) cat('\nCalculating basic pairwise relationships for windows (n=',nrow(all.classifications),')...', sep="")
  all.classifications	<- all.classifications %>% categorise("basic.classification", "other",
                     list(contiguous=TRUE, adjacent=TRUE, ancestry=c("anc", "multiAnc"), label="anc_contiguous"),
                     list(contiguous=FALSE, adjacent=TRUE, ancestry=c("anc", "multiAnc"), label="anc_noncontiguous"),
                     list(contiguous=TRUE, adjacent=TRUE, ancestry=c("desc", "multiDesc"), label="desc_contiguous"),
                     list(contiguous=FALSE, adjacent=TRUE, ancestry=c("desc", "multiDesc"), label="desc_noncontiguous"),
                     list(contiguous=TRUE, adjacent=TRUE, ancestry="none", label="noancestry_contiguous"),
                     list(contiguous=FALSE, adjacent=TRUE, ancestry="none", label="noancestry_noncontiguous"),
                     list(contiguous=TRUE, adjacent=TRUE, ancestry="complex", label="complex_contiguous"),
                     list(contiguous=FALSE, adjacent=TRUE, ancestry="complex", label="complex_noncontiguous")
  )
  
  
  all.classifications	<- all.classifications %>% categorise("basic.classification", "other",
                     list(categorical.distance = "close", label="close"),
                     list(categorical.distance = "distant", label="distant"),
                     list(categorical.distance = "intermediate", label="intermediate"))
  
  if(verbose) cat('\nCalculating derived pairwise relationships for windows (n=',nrow(all.classifications),')...', sep="")
  all.classifications <- all.classifications %>% get.pairwise.relationships(get.groups=relationship.types)
  
  if(verbose) cat('\nCalculating KEFF and NEFF for windows (n=',nrow(all.classifications),')...', sep="")
  rplkl	<- dwin %>% get.keff.and.neff(relationship.types)
  
  if(verbose) cat('\nCalculating posterior state probabilities for pairs and relationship groups (n=',nrow(rplkl),')...', sep="")
  rplkl	<- phsc.get.pairwise.relationships.posterior(rplkl, n.type=prior.keff, n.obs=prior.neff, confidence.cut=prior.calibrated.prob)
  
  list(dwin=dwin, rplkl=rplkl)
}

#' @title Calculate pairwise relationships
#' @description This function calculates pairwise relationships between each pair of two individuals in any window. Several different relationship groups can be calculated, for example just using pairwise distance, or using both pairwise distance and topology to define likely pairs.
#' @export    
#' @param df data.table to which new columns of relationship groups will be added. Must contain a column with name TYPE_BASIC. This column contains fundamental relationship states for every window, from which other relationships are derived. 
#' @param get.groups names of relationship groups  
#' @keywords internal
#' @return input data.table with new columns. Each new column defines relationship states for a specific relationship group. 
 
get.pairwise.relationships <- function(df, get.groups=c('proximity.3.way',
                                                        'any.ancestry',
                                                        'close.x.contiguous',
                                                        'close.and.contiguous',
                                                        'close.and.contiguous.and.directed',
                                                        'close.and.adjacent.and.directed',
                                                        'close.and.contiguous.and.ancestry.cat',
                                                        'close.and.adjacent.and.ancestry.cat',
                                                        'adjacent.and.proximity.cat'
                                                        )){
  
  if('proximity.3.way' %in% get.groups) {
    df <- df %>% categorise("proximity.3.way", "intermediate",
                     list(categorical.distance = "close", label="close"),
                     list(categorical.distance = "distant", label="distant"))
  }	
  #	
  #	group to define likely pair just based on topology
  #	
  if('any.ancestry' %in% get.groups) {
    df <- df %>% categorise("any.ancestry", "other",
                     list(contiguous = TRUE, ancestry = c("anc_12", "anc_21", "multi_anc_12", "multi_anc_21", "complex"), label = "anc.or.complex"))
  }
  
  #
  #	group for full 2x2 table of distance and topology
  #
  if('close.x.contiguous' %in% get.groups) {		
    df <- df %>% categorise("close.x.contiguous", "not.close.other",
                     list(contiguous = TRUE, categorical.distance = "close", label = "contiguous_close"),
                     list(contiguous = TRUE, categorical.distance = c("intermediate", "distant"), label = "contiguous_not.close"),
                     list(contiguous = FALSE, categorical.distance = "close", label = "noncontiguous_close"))
  }
  #
  #	group to determine linkage - only 2 states, linked / unlinked
  #
  if('close.and.contiguous' %in% get.groups) {
    df <- df %>% categorise("close.and.contiguous", "not.close.or.noncontiguous",
                     list(contiguous = TRUE, categorical.distance = "close", label="close_contiguous"))
  }	
  #
  #	group to determine direction of transmission in likely pairs
  #	based on contiguous linkage
  #
  if('close.and.contiguous.and.directed' %in% get.groups) {
    df <- df %>% categorise("close.and.contiguous.and.directed", NA_character_,
                     list(contiguous = TRUE, ancestry=c("anc_12", "multi_anc_12"), categorical.distance = "close", label="12"),
                     list(contiguous = TRUE, ancestry=c("anc_21", "multi_anc_21"), categorical.distance = "close", label="21"))
  }
  #
  #	group to determine direction of transmission in likely pairs
  #	based on adjacent linkage
  #
  if('close.and.adjacent.and.directed' %in% get.groups) {
    df	<- df %>% categorise("close.and.adjacent.and.directed", NA_character_,
                     list(adjacent = TRUE, ancestry=c("anc_12", "multi_anc_12"), categorical.distance = "close", label="12"),
                     list(adjacent = TRUE, ancestry=c("anc_21", "multi_anc_21"), categorical.distance = "close", label="21"))
  }
  #
  #	group to determine probabilities for transmission networks
  #	based on contiguous linkage
  #
  if('close.and.contiguous.and.ancestry.cat' %in% get.groups) {				
    df	<- df %>% categorise("close.and.contiguous.and.ancestry.cat", "not.close.or.noncontiguous",
                     list(contiguous = TRUE, ancestry=c("anc_12", "multi_anc_12"), categorical.distance = "close", label = "12"),
                     list(contiguous = TRUE, ancestry=c("anc_21", "multi_anc_21"), categorical.distance = "close", label = "21"),
                     list(contiguous = TRUE, ancestry=c("none", "complex"), categorical.distance = "close", label = "complex.or.no.ancestry"))
  }
  #
  #	group to determine probabilities for transmission networks
  #	based on adjacent linkage
  #
  if('close.and.adjacent.and.ancestry.cat' %in% get.groups) {				
    df <- df %>% categorise("close.and.adjacent.and.ancestry.cat", "not.close.or.nonadjacent",
                     list(adjacent = TRUE, ancestry=c("anc_12", "multi_anc_12"), categorical.distance = "close", label = "12"),
                     list(adjacent = TRUE, ancestry=c("anc_21", "multi_anc_21"), categorical.distance = "close", label = "21"),
                     list(adjacent = TRUE, ancestry=c("none", "complex"), categorical.distance = "close", label = "complex.or.no.ancestry"))
  }		
  #
  #	group to transmission networks
  #	
  if('adjacent.and.proximity.cat' %in% get.groups) {
    df <- df %>% categorise("adjacent.and.proximity.cat", "distant",
                     list(categorical.distance = "close", adjacent=TRUE, label = "cluster"),
                     list(categorical.distance = "intermediate", adjacent=TRUE, label="ambiguous"))
  }
  
  df
}


#' @title Count observed relationship states
#' @description This function counts for each pair of individuals their relationship states across all windows (KEFF), and the total number of windows (NEFF). Since windows can be overlapping, the count is a real value.
#' @export    
#' @keywords internal
#' @param df data.table with basic relationship types for paired individuals across windows. Must contain columns 'ID1','ID2','W_FROM','W_TO','TYPE_BASIC'. 
#' @param get.groups names of relationship groups
#' @import data.table  
#' @return new data.table with columns host.1 host.2 GROUP TYPE K KEFF N NEFF.
#'  
get.keff.and.neff <- function(df, get.groups, w.slide=NA){
  
  stopifnot(c('host.1', 'host.2', 'window.start', 'window.end', 'basic.classification') %in% colnames(df))
  
  df <- df %>% arrange(host.1, host.2, window.start)
  
  #	identify chunks of contiguous windows
  
  if(is.na(w.slide)) {
    w.slide <- df %>% group_by(host.1, host.2) %>% summarise(slide.width = min(diff(window.start))) %>% ungroup()
    
    w.slide	<- w.slide %>% filter(!is.na(slide.width))
    w.slide	<- ifelse(nrow(w.slide), min(w.slide$slide.width), 1L)		
  }
  #	define chunks

  tmp		<- df[, {
    tmp<- as.integer( c(TRUE,(W_FROM[-length(W_FROM)]+w.slide)!=W_FROM[-1]) )
    list(W_FROM=W_FROM, W_TO=W_TO, CHUNK=cumsum(tmp))
  }, by=c('ID1','ID2')]
  df		<- merge(df,tmp,by=c('ID1','ID2','W_FROM','W_TO'))
  #	define chunk length in terms of non-overlapping windows	& number of windows in chunk
  tmp		<- df[, {
    list(W_FROM=W_FROM, W_TO=W_TO, CHUNK_L=(max(W_TO+1L)-min(W_FROM))/(W_TO[1]+1L-W_FROM[1]), CHUNK_N=length(W_FROM))
  }, by=c('ID1','ID2','CHUNK')]
  df		<- merge(df,tmp,by=c('ID1','ID2','CHUNK','W_FROM','W_TO'))
  
  categorisation <- df[,c('TYPE_BASIC', get.groups), with=F]
  setkey(categorisation)
  categorisation <- unique(categorisation)
  
  setkey(categorisation, "TYPE_BASIC")
  
  #	for each chunk, count: windows by type and effective length of chunk
  #	then sum chunks
  rplkl	<- df[, list(	K= length(W_FROM), KEFF= length(W_FROM)/CHUNK_N[1] * CHUNK_L[1]), by=c('ID1','ID2','CHUNK','TYPE_BASIC')]	
  rplkl	<- rplkl[, list(STAT=c('K','KEFF'), V=c(sum(K),sum(KEFF))), by=c('ID1','ID2','TYPE_BASIC')]
  #
  #	add relationship types
  #
  
  setkey(rplkl, "TYPE_BASIC")
  rplkl <- rplkl[categorisation]
  
  #	melt relationship groups
  rplkl	<- melt(rplkl, measure.vars=c(get.groups,'TYPE_BASIC'), variable.name='GROUP', value.name='TYPE')
  rplkl	<- subset(rplkl, !is.na(TYPE))
  #	sum K and KEFF of same relationship state
  rplkl	<- rplkl[, list(V=sum(V)), by=c('ID1','ID2','GROUP','TYPE','STAT')]
  #	add zero-count relationship states (change to wide table and set NA's to zero's)
  tmp		<- unique(subset(rplkl, select=c(GROUP,TYPE)))
  tmp2	<- unique(subset(rplkl, select=c(ID1,ID2, STAT)))
  tmp[, DUMMY:=1L]
  tmp2[, DUMMY:=1L]
  tmp		<- merge(tmp, tmp2, by='DUMMY',allow.cartesian=TRUE)
  set(tmp, NULL, 'DUMMY', NULL)
  rplkl	<- merge(tmp, rplkl, all.x=1, by=c('ID1','ID2','GROUP','TYPE','STAT'))
  set(rplkl, rplkl[,which(is.na(V))], 'V', 0)	
  #	expand KEFF and K columns now that everything is done
  rplkl	<- dcast.data.table(rplkl, ID1+ID2+GROUP+TYPE~STAT, value.var='V')	
  #	calculate N and NEFF
  tmp		<- rplkl[, list(N= sum(K), NEFF= sum(KEFF)), by=c('ID1','ID2','GROUP')]	
  rplkl	<- merge(rplkl, tmp, by=c('ID1','ID2','GROUP'))
  rplkl
}

#' @title Calculate marginal posterior probability for two individuals being in a particular relationship state
#' @description This function calculates the parameters that specify the marginal posterior probability for two individuals being in a particular relationship state. The marginal posterior is Beta distributed and this function calculates the ALPHA and BETA parameters.
#' @export  
#' @param df Input data.table   
#' @param n.type Calibration parameter for the prior: minimum number of windows of state to select a pair of individuals with confidence of at least at least confidence.cut, if the total number of windows is n.obs
#' @param n.obs Calibration parameter for the prior: total number of windows. 
#' @param confidence.cut Calibration parameter for the prior: confidence cut off.  
#' @return Input data.table with two additional columns POSTERIOR_ALPHA and POSTERIOR_BETA
#' @keywords internal
phsc.get.pairwise.relationships.posterior <- function(df, n.type=2, n.obs=3, n.type.dir=n.type, n.obs.dir=n.obs, confidence.cut=0.5) {
  
  stopifnot(c('GROUP')%in%colnames(df))
  tmp		<- phsc.get.pairwise.relationships.numbers()	
  tmp		<- tmp[, {
    #if(!grepl('_DIR',GROUP))
    #	z<- phsc.get.prior.parameter.n0(N_TYPE, keff=n.type, neff=n.obs, confidence.cut=confidence.cut)
    #if(grepl('_DIR',GROUP))
    #	z<- phsc.get.prior.parameter.n0(N_TYPE, keff=n.type.dir, neff=n.obs.dir, confidence.cut=confidence.cut)
    #list(PAR_PRIOR=z)
    list(PAR_PRIOR=N_TYPE)
  }, by=c('GROUP','N_TYPE')]
  df		<- merge(df, tmp, by=c('GROUP'))
  # df[, POSTERIOR_ALPHA:= PAR_PRIOR/N_TYPE+KEFF]
  df[, POSTERIOR_ALPHA:= PAR_PRIOR/N_TYPE+KEFF]
  df[, POSTERIOR_BETA:= PAR_PRIOR*(1-1/N_TYPE)+NEFF-KEFF]	
  df[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-N_TYPE)]
  df
}


#' @export
#' @keywords internal
#' @title Generate a new column, or append to an existing column, in a \code{data.frame} which determines whether each row belongs to one or more categories determiend by other columns in the table.
#' @param df A \code{data.frame}
#' @param col.name The name of the column to be added or appended to
#' @param no.match.result A string that will identify rows that do not match any of the categories.
#' @param ... One or more lists. Each entry in each list, except \code{label} should take the name of a column in \code{df} and a value or vector of values which are the acceptable values for that column in this category. The \code{label} item gives the name for this category, which will be automatically generated in a long-winded fashion if \code{label} is absent. Categories must be mutually exclusive. 
#' @return If \code{name} is not an existing column in \code{df}, then \code{df} with a new column \code{name} which contains the category names. If \code{name} does exist, then its entries are appended with the category names following an underscore.

categorise <- function(df, col.name, no.match.result = NA_character_, ...){
  
  # TODO this stuff is ugly and a complete overhaul at some point would be nice
  
  conditions <- list(...)
  adding.column <- !(col.name %in% colnames(df))
  
  if(adding.column){
    df <- df %>% mutate(!!col.name := NA_character_)
  }
  
  leftovers <- 1:nrow(df)
  
  for(condition in conditions){

    if(!is.null(condition$label)){
      value.string <- condition$label
      if(!adding.column){
        value.string <- paste0("_", value.string)
      }
    } else {
      if(adding.column){
        value.string <- ""
      } else {
        value.string <- "_"
      }
    }
    
    first <- T
    rows <- 1:nrow(df)
    
    for(column.no in 1:length(condition)){
      
      if(is.null(condition$label)){
        if(first){
          first <- F
        } else {
          value.string <- paste0(value.string, "_")
        }
      }
      queried.col.name <- labels(condition)[column.no]
      
      if(queried.col.name != "label"){
        col.value <- condition[[column.no]]
        
        if(is.null(condition$label)){
          if(is.logical(col.value)){
            if(col.value){
              value.string <- paste0(value.string, queried.col.name)
            } else {
              value.string <- paste0(value.string, "not.", queried.col.name)
            }
          } else {
            value.string <- paste0(value.string, queried.col.name, ".", paste(col.value, collapse="."))
          }
        }
        
        kept.rows <- which(pull(df, queried.col.name) %in% col.value)   
        rows <- intersect(rows, kept.rows)
      }
    }
    
    if(adding.column){
      df <- df %>% mutate(!!col.name := replace(get(col.name), rows, value.string))
    } else {
      df <- df %>% mutate(!!col.name := replace(get(col.name), rows, paste0(df %>% slice(rows) %>% pull(get(col.name)), value.string)))
    }
    leftovers <- setdiff(leftovers, rows)

  }
  
  if(length(leftovers)>0){
    if(adding.column){
      df <- df %>% mutate(!!col.name := replace(get(col.name), leftovers, no.match.result))
    } else {
      df <- df %>% mutate(!!col.name := replace(get(col.name), rows, paste0(df %>% slice(leftovers) %>% pull(get(col.name)), "_", no.match.result)))
    }
  }
  df
}


#' @title Calculate prior parameter n0
#' @description This function calculates the prior parameter n0 such that the minimum number of windows needed is at least n.obs in order to select a pair of individuals with confidence at least confidence.cut.    
#' @param n.states Number of relationship states.
#' @param n.obs Minimum number of windows to select a pair of individuals with confidence at least confidence.cut. 
#' @param confidence.cut Confidence cut off.  
#' @keywords internal
#' @return Prior parameter n0 
phsc.get.prior.parameter.n0 <- function(n.states, keff=2, neff=3, confidence.cut=0.66) {
  #phsc.find.n0.aux<- function(n0, n.states, keff, neff, confidence.cut)
  #{
  #	abs( (n0+n.states*(keff-1))/ (n.states*(neff+n0-2)) - confidence.cut )
  #}
  phsc.find.n0.aux<- function(n0, n.states, keff, neff, confidence.cut)
  {
    abs( (n0+n.states*keff)/ (n.states*(neff+n0)) - confidence.cut )
  }	
  
  ans	<- optimize(phsc.find.n0.aux, c(.001,1e2), n.states=n.states, keff=keff, neff=neff, confidence.cut=confidence.cut)
  ans	<- round(ans$minimum, d=4)
  ans
}

phsc.get.pairwise.relationships.numbers <- function() {
  tmp	<- matrix(c('TYPE_PAIR_DI2','3',
                  'TYPE_PAIR_TO','2',
                  'TYPE_PAIR_TODI2x2','4',					
                  'TYPE_PAIR_TODI2','2',
                  'TYPE_DIR_TODI2','2',															
                  'TYPE_CHAIN_TODI','3',
                  'TYPE_ADJ_DIR_TODI2','2',
                  'TYPE_NETWORK_SCORES','4',
                  'TYPE_ADJ_NETWORK_SCORES','4',
                  'TYPE_BASIC','24'), ncol=2,byrow=TRUE)
  colnames(tmp)	<- c('GROUP','N_TYPE')
  tmp				<- as.data.table(tmp)
  set(tmp, NULL, 'N_TYPE', tmp[, as.integer(N_TYPE)])
  tmp
}

