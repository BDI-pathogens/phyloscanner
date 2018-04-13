

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
multinomial.calculations <- function(phyloscanner.trees, 
                                     close.threshold,
                                     prior.keff = 3,
                                     prior.neff = 4,
                                     prior.calibrated.prob = 0.66,
                                     tip.regex = "^(.*)_read_([0-9]+)_count_([0-9]+)$", 
                                     allow.mt = F, 
                                     min.reads = 0, 
                                     min.tips = 0, 
                                     distant.threshold = close.threshold,
                                     relationship.types	= c('TYPE_PAIR_DI2','TYPE_PAIR_TO','TYPE_PAIR_TODI2x2','TYPE_PAIR_TODI2','TYPE_DIR_TODI2','TYPE_NETWORK_SCORES','TYPE_ADJ_NETWORK_SCORES','TYPE_CHAIN_TODI'),
                                     verbose=F){
  
  if(!attr(phyloscanner.trees, "readable.coords")){
    stop("No window cooardinates detected in this tree set, cannot do multinomial calculations.")
  }
  
  mc <- merge.classifications(phyloscanner.trees, allow.mt, verbose)
  
  ss <- lapply(phyloscanner.trees, function(x) get.tip.and.read.counts(x, all.hosts.from.trees(phyloscanner.trees), tip.regex, attr(phyloscanner.trees, 'has.read.counts'), verbose))
  ss <- rbindlist(ss)

  setkey(ss, id, file.id)
  setkey(mc, HOST.1, ID)
  merge.tab.1 <- mc[ss, nomatch=0]
  setnames(merge.tab.1, c('tips', 'reads'), c('tips.1', 'reads.1'))
  
  setkey(merge.tab.1, HOST.2, ID)
  merge.tab.2 <- merge.tab.1[ss, nomatch=0]
  setnames(merge.tab.2, c('tips', 'reads'), c('tips.2', 'reads.2'))
  
  dwin <- merge.tab.2
  
  setnames(dwin, c('HOST.1', 'HOST.2', 'tips.1', 'reads.1', 'tips.2', 'reads.2', 'PATHS.12', 'PATHS.21'), 
           c('ID1', 'ID2', 'ID1_L', 'ID1_R', 'ID2_L', 'ID2_R', 'PATHS_12', 'PATHS_21'))
  
  set(dwin, NULL, 'W_FROM', dwin[, as.integer(gsub('[^0-9]*([0-9]+)_to_([0-9]+).*','\\1', ID))])
  set(dwin, NULL, 'W_TO', dwin[, as.integer(gsub('[^0-9]*([0-9]+)_to_([0-9]+).*','\\2', ID))])
  
  dwin$CATDISTANCE <- sapply(dwin$PATRISTIC_DISTANCE, function(x){
    if(x < close.threshold){
      "close"
    } else if(x >= distant.threshold){
      "distant"
    } else {
      "intermediate"
    }
  })
  
  if(verbose) cat('\nReducing transmission window stats to windows with at least',min.reads,'reads and at least',min.tips,'tips ...')
  dwin	<- subset(dwin, ID1_R>=min.reads & ID2_R>=min.reads & ID1_L>=min.tips & ID2_L>=min.tips)
  if(verbose) cat('\nTotal number of windows with transmission assignments is',nrow(dwin),'...')		
  
  if(verbose) cat('\nCalculating basic pairwise relationships for windows n=',nrow(dwin),'...')
  dwin	<- categorise(dwin, "TYPE_BASIC", "other",
                      list(CONTIGUOUS=TRUE, ADJACENT=TRUE, TYPE=c("anc_12", "multi_anc_12"), label="trans12_contiguous"),
                      list(CONTIGUOUS=FALSE, ADJACENT=TRUE, TYPE=c("anc_12", "multi_anc_12"), label="trans12_noncontiguous"),
                      list(CONTIGUOUS=TRUE, ADJACENT=TRUE, TYPE=c("anc_21", "multi_anc_21"), label="trans21_contiguous"),
                      list(CONTIGUOUS=FALSE, ADJACENT=TRUE, TYPE=c("anc_21", "multi_anc_21"), label="trans21_noncontiguous"),
                      list(CONTIGUOUS=TRUE, ADJACENT=TRUE, TYPE="none", label="noancestry_contiguous"),
                      list(CONTIGUOUS=FALSE, ADJACENT=TRUE, TYPE="none", label="noancestry_noncontiguous"),
                      list(CONTIGUOUS=TRUE, ADJACENT=TRUE, TYPE="complex", label="complex_contiguous"),
                      list(CONTIGUOUS=FALSE, ADJACENT=TRUE, TYPE="complex", label="complex_noncontiguous")
  )
  

  dwin	<- categorise(dwin, "TYPE_BASIC", "OTHER",
                      list(CATDISTANCE = "close", label="close"),
                      list(CATDISTANCE = "distant", label="distant"),
                      list(CATDISTANCE = "intermediate", label="intermediate"))
  
  if(verbose) cat('\nCalculating derived pairwise relationships for windows n=',nrow(dwin),'...')
  dwin	<- phsc.get.pairwise.relationships(dwin, get.groups=relationship.types)

  if(verbose) cat('\nCalculating KEFF and NEFF for windows n=',nrow(dwin),'...')
  rplkl	<- phsc.get.pairwise.relationships.keff.and.neff(dwin, relationship.types)
  
  if(verbose) cat('\nCalculating posterior state probabilities for pairs and relationship groups n=',nrow(rplkl),'...')
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
phsc.get.pairwise.relationships<- function(df, get.groups=c('TYPE_PAIR_DI2','TYPE_PAIR_TO','TYPE_PAIR_TODI2x2','TYPE_PAIR_TODI2','TYPE_DIR_TODI2','TYPE_NETWORK_SCORES','TYPE_CHAIN_TODI','TYPE_ADJ_NETWORK_SCORES','TYPE_ADJ_DIR_TODI2'))
{

  if('TYPE_PAIR_DI2' %in% get.groups)
  {
    df	<- categorise(df, "TYPE_PAIR_DI2", "not.close",
                      list(CATDISTANCE = "close", label="close"))
  }	
  #	
  #	group to define likely pair just based on topology
  #	
  if('TYPE_PAIR_TO' %in% get.groups)
  {
    df	<- categorise(df, "TYPE_PAIR_TO", "other",
                      list(CONTIGUOUS = TRUE, TYPE = c("anc_12", "anc_21", "multi_anc_12", "multi_anc_21", "complex"), label = "anc.or.complex"))
  }

  #
  #	group for full 2x2 table of distance and topology
  #
  if('TYPE_PAIR_TODI2x2' %in% get.groups)
  {		
    df	<- categorise(df, "TYPE_PAIR_TODI2x2", "not.close.other",
                      list(CONTIGUOUS = TRUE, CATDISTANCE = "close", label = "contiguous_close"),
                      list(CONTIGUOUS = TRUE, CATDISTANCE = c("intermediate", "distant"), label = "contiguous_not.close"),
                      list(CONTIGUOUS = FALSE, CATDISTANCE = "close", label = "noncontiguous_close"))
  }
  #
  #	group to determine linkage - only 2 states, linked / unlinked
  #
  if('TYPE_PAIR_TODI2' %in% get.groups)
  {
    df	<- categorise(df, "TYPE_PAIR_TODI2", "not.close.or.noncontiguous",
                      list(CONTIGUOUS = TRUE, CATDISTANCE = "close", label="close_contiguous"))
  }	
  #
  #	group to determine direction of transmission in likely pairs
  #	based on contiguous linkage
  #
  if('TYPE_DIR_TODI2' %in% get.groups)
  {
    df	<- categorise(df, "TYPE_DIR_TODI2", NA_character_,
                      list(CONTIGUOUS = TRUE, TYPE=c("anc_12", "multi_anc_12"), CATDISTANCE = "close", label="12"),
                      list(CONTIGUOUS = TRUE, TYPE=c("anc_21", "multi_anc_21"), CATDISTANCE = "close", label="21"))
  }
  #
  #	group to determine direction of transmission in likely pairs
  #	based on adjacent linkage
  #
  if('TYPE_ADJ_DIR_TODI2' %in% get.groups)
  {
    df	<- categorise(df, "TYPE_ADJ_DIR_TODI2", NA_character_,
                      list(ADJACENT = TRUE, TYPE=c("anc_12", "multi_anc_12"), CATDISTANCE = "close", label="12"),
                      list(ADJACENT = TRUE, TYPE=c("anc_21", "multi_anc_21"), CATDISTANCE = "close", label="21"))
  }
  #
  #	group to determine probabilities for transmission networks
  #	based on contiguous linkage
  #
  if('TYPE_NETWORK_SCORES'%in%get.groups)
  {				
    df	<- categorise(df, "TYPE_NETWORK_SCORES", "not.close.or.noncontiguous",
                      list(CONTIGUOUS = TRUE, TYPE=c("anc_12", "multi_anc_12"), CATDISTANCE = "close", label = "12"),
                      list(CONTIGUOUS = TRUE, TYPE=c("anc_21", "multi_anc_21"), CATDISTANCE = "close", label = "21"),
                      list(CONTIGUOUS = TRUE, TYPE=c("none", "complex"), CATDISTANCE = "close", label = "complex.or.no.ancestry"))
  }
  #
  #	group to determine probabilities for transmission networks
  #	based on adjacent linkage
  #
  if('TYPE_ADJ_NETWORK_SCORES' %in% get.groups)
  {				
    df	<- categorise(df, "TYPE_ADJ_NETWORK_SCORES", "not.close.or.noncontiguous",
                      list(ADJACENT = TRUE, TYPE=c("anc_12", "multi_anc_12"), CATDISTANCE = "close", label = "12"),
                      list(ADJACENT = TRUE, TYPE=c("anc_21", "multi_anc_21"), CATDISTANCE = "close", label = "21"),
                      list(ADJACENT = TRUE, TYPE=c("none", "complex"), CATDISTANCE = "close", label = "complex.or.no.ancestry"))
  }		
  #
  #	group to transmission networks
  #	
  if('TYPE_CHAIN_TODI' %in% get.groups)
  {
    df	<- categorise(df, "TYPE_CHAIN_TODI", "distant",
                      list(CATDISTANCE = "close", label = "close"))
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
#' @return new data.table with columns ID1 ID2 GROUP TYPE K KEFF N NEFF. 
phsc.get.pairwise.relationships.keff.and.neff<- function(df, get.groups, w.slide=NA)
{
  stopifnot(c('ID1','ID2','W_FROM','W_TO','TYPE_BASIC')%in%colnames(df))
  #
  #	identify chunks of contiguous windows
  #	
  setkey(df, ID1, ID2, W_FROM)
  if(is.na(w.slide))
  {
    w.slide	<- df[, {
      ans	<- NA_integer_
      tmp	<- diff(W_FROM)
      if(length(tmp))
        ans	<- min(tmp)
      list(W_SLIDE=ans)
    }, by=c('ID1','ID2')]
    w.slide	<- subset(w.slide, !is.na(W_SLIDE))
    w.slide	<- ifelse(nrow(w.slide), w.slide[, min(W_SLIDE)], 1L)		
  }
  #	define chunks
  setkey(df, ID1, ID2, W_FROM)
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
phsc.get.pairwise.relationships.posterior<- function(df, n.type=2, n.obs=3, n.type.dir=n.type, n.obs.dir=n.obs, confidence.cut=0.5)
{
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
#' @title Generate a new column, or append to an existing column, in a \code{data.table} which determines whether each row belongs to one or more categories determiend by other columns in the table.
#' @param df A data.table
#' @param name The name of the column to be added or appended to
#' @param no.match.result A stringt that will identify rows that do not match any of the categories.
#' @param ... One or more lists. Each entry in each list, except \code{label} should take the name of a column in \code{df} and a value or vector of values which are the acceptable values for that column in this category. The \code{label} item gives the name for this category, which will be automatically generated in a long-winded fashion if \code{label} is absent. Categories must be mutually exclusive. 
#' @return If \code{name} is not an existing column in \code{df}, then \code{df} with a new column \code{name} which contains the category names. If \code{name} does exist, then its entries are appended with the category names following an underscore.
categorise <- function(df, name, no.match.result = NA_character_, ...){
  conditions <- list(...)
  adding.column <- !(name %in% colnames(df))
  
  if(adding.column){
    df[, c(name):= NA_character_]
  }
  
  leftovers <- 1:nrow(df)
  
  for(condition in conditions){
    
    if(!is.null(condition$label)){
      string <- condition$label
      if(!adding.column){
        string <- paste0("_", string)
      }
    } else {
      if(adding.column){
        string <- ""
      } else {
        string <- "_"
      }
    }
    
    first <- T
    cols <- 1:nrow(df)
    for(column.no in 1:length(condition)){
      
      if(is.null(condition$label)){
        if(first){
          first <- F
        } else {
          string <- paste0(string, "_")
        }
      }
      col.name <- labels(condition)[column.no]
      col.value <- condition[[column.no]]
      
      if(col.name != "label"){
        
        if(is.null(condition$label)){
          if(is.logical(col.value)){
            if(col.value){
              string <- paste0(string, col.name)
            } else {
              string <- paste0(string, "NOT", ".", col.name)
            }
          } else {
            string <- paste0(string, col.name, ".", paste(col.value, collapse="."))
          }
        }
        
        kept.cols <-  df[, which(get(col.name) %in% col.value)]
        
        cols <- intersect(cols, kept.cols)
      }
    }
    
    
    if(adding.column){
      set(df, cols, name, string)
    } else {
      set(df, cols, name, df[cols, paste0(get(name),string)])
    }
    leftovers <- setdiff(leftovers, cols)
    
  }
  
  if(length(leftovers)>0){
    if(adding.column){
      set(df, leftovers, name, no.match.result)
    } else {
      set(df, leftovers, name, df[leftovers, paste0(get(name),"_", no.match.result)])
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
phsc.get.prior.parameter.n0<- function(n.states, keff=2, neff=3, confidence.cut=0.66)
{
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

phsc.get.pairwise.relationships.numbers<- function()
{
  tmp	<- matrix(c('TYPE_PAIR_DI2','2',
                  'TYPE_PAIR_TO','2',
                  'TYPE_PAIR_TODI2x2','4',					
                  'TYPE_PAIR_TODI2','2',
                  'TYPE_DIR_TODI2','2',															
                  'TYPE_CHAIN_TODI','2',
                  'TYPE_ADJ_DIR_TODI2','2',
                  'TYPE_NETWORK_SCORES','3',
                  'TYPE_ADJ_NETWORK_SCORES','3',
                  'TYPE_BASIC','24'), ncol=2,byrow=TRUE)
  colnames(tmp)	<- c('GROUP','N_TYPE')
  tmp				<- as.data.table(tmp)
  set(tmp, NULL, 'N_TYPE', tmp[, as.integer(N_TYPE)])
  tmp
}

