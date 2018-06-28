

#' @title Calculate parameters of the posterior density for pairwise host relationships
#' @export
#' @param ptrees A list of class \code{phyloscanner.trees} produced by \code{phyloscanner.analyse.trees}.
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
  
  if(verbose) cat('Merging tip and read counts...\n')
  
  all.classifications <- all.classifications %>% 
    inner_join(host.tips.and.reads, by=c("host.1"="host.id", "tree.id")) %>% 
    rename(tips.1 = tips, reads.1=reads)
  all.classifications <- all.classifications %>% 
    inner_join(host.tips.and.reads, by=c("host.2"="host.id", "tree.id")) %>% 
    rename(tips.2 = tips, reads.2=reads)

  if(verbose) cat('Finding window coordinates...\n')
  
  all.classifications <- all.classifications %>% 
    mutate(window.start = map_int(tree.id, function(x) as.integer(gsub('[^0-9]*([0-9]+)_to_([0-9]+).*','\\1', x))))
  all.classifications <- all.classifications %>% 
    mutate(window.end = map_int(tree.id, function(x) as.integer(gsub('[^0-9]*([0-9]+)_to_([0-9]+).*','\\2', x))))
  
  if(verbose) cat('Assigning discrete proximity categories...\n')
  
  all.classifications <- all.classifications %>% 
    mutate(categorical.distance = map_chr(patristic.distance, function(x){
    if(x < close.threshold){
      "close"
    } else if(x >= distant.threshold){
      "distant"
    } else {
      "intermediate"
    }
  }))
  
  if(verbose) cat('Reducing transmission window stats to windows with at least',min.reads,'reads and at least',min.tips,'tips...\n')
  
  all.classifications	<- all.classifications %>% 
    filter(reads.1 >= min.reads & reads.2 >= min.reads & tips.1 >= min.tips & tips.2 >= min.tips)
  
  if(verbose) cat('Total number of windows with transmission assignments is ',nrow(all.classifications),'.\n', sep="")		
  
  if(verbose) cat('Calculating basic pairwise relationships for windows (n=',nrow(all.classifications),')...\n', sep="")
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
  
  if(verbose) cat('Calculating derived pairwise relationships for windows (n=',nrow(all.classifications),')...\n', sep="")
  all.classifications <- all.classifications %>% get.pairwise.relationships(relationship.types=relationship.types)
  
  if(verbose) cat('Calculating k.eff and n.eff for windows (n=',nrow(all.classifications),')...\n', sep="")
  category.parameters	<- all.classifications %>% 
    get.keff.and.neff(relationship.types)
  
  if(verbose) cat('Calculating posterior state probabilities for pairs and relationship groups (n=',nrow(category.parameters),')...\n', sep="")
  category.parameters	<- category.parameters %>% 
    get.posterior.scores()
  
  list(dwin=all.classifications, rplkl=category.parameters)
}

#' @title Calculate pairwise relationships
#' @description This function calculates pairwise relationships between each pair of two individuals in any window. Several different relationship groups can be calculated, for example just using pairwise distance, or using both pairwise distance and topology to define likely pairs.
#' @export    
#' @param df tibble to which new columns of relationship groups will be added. Must contain a column with name TYPE_BASIC. This column contains fundamental relationship states for every window, from which other relationships are derived. 
#' @param relationship.types names of relationship groups  
#' @keywords internal
#' @return input tibble with new columns. Each new column defines relationship states for a specific relationship group. 
 
get.pairwise.relationships <- function(df, relationship.types=c('proximity.3.way',
                                                        'any.ancestry',
                                                        'close.x.contiguous',
                                                        'close.and.contiguous',
                                                        'close.and.contiguous.and.directed',
                                                        'close.and.adjacent.and.directed',
                                                        'close.and.contiguous.and.ancestry.cat',
                                                        'close.and.adjacent.and.ancestry.cat',
                                                        'adjacent.and.proximity.cat'
                                                        )){
  
  if('proximity.3.way' %in% relationship.types) {
    df <- df %>% categorise("proximity.3.way", "intermediate",
                     list(categorical.distance = "close", label="close"),
                     list(categorical.distance = "distant", label="distant"))
  }	
  #	
  #	group to define likely pair just based on topology
  #	
  if('any.ancestry' %in% relationship.types) {
    df <- df %>% categorise("any.ancestry", "other",
                     list(contiguous = TRUE, ancestry = c("anc", "desc", "multiAnc", "multiDesc", "complex"), label = "anc.or.complex"))
  }
  
  #
  #	group for full 2x2 table of distance and topology
  #
  if('close.x.contiguous' %in% relationship.types) {		
    df <- df %>% categorise("close.x.contiguous", "not.close.other",
                     list(contiguous = TRUE, categorical.distance = "close", label = "contiguous_close"),
                     list(contiguous = TRUE, categorical.distance = c("intermediate", "distant"), label = "contiguous_not.close"),
                     list(contiguous = FALSE, categorical.distance = "close", label = "noncontiguous_close"))
  }
  #
  #	group to determine linkage - only 2 states, linked / unlinked
  #
  if('close.and.contiguous' %in% relationship.types) {
    df <- df %>% categorise("close.and.contiguous", "not.close.or.noncontiguous",
                     list(contiguous = TRUE, categorical.distance = "close", label="close_contiguous"))
  }	
  #
  #	group to determine direction of transmission in likely pairs
  #	based on contiguous linkage
  #
  if('close.and.contiguous.and.directed' %in% relationship.types) {
    df <- df %>% categorise("close.and.contiguous.and.directed", NA_character_,
                     list(contiguous = TRUE, ancestry=c("anc", "multiAnc"), categorical.distance = "close", label="12"),
                     list(contiguous = TRUE, ancestry=c("desc", "multiDesc"), categorical.distance = "close", label="21"))
  }
  #
  #	group to determine direction of transmission in likely pairs
  #	based on adjacent linkage
  #
  if('close.and.adjacent.and.directed' %in% relationship.types) {
    df	<- df %>% categorise("close.and.adjacent.and.directed", NA_character_,
                     list(adjacent = TRUE, ancestry=c("anc", "multiAnc"), categorical.distance = "close", label="12"),
                     list(adjacent = TRUE, ancestry=c("desc", "multiDesc"), categorical.distance = "close", label="21"))
  }
  #
  #	group to determine probabilities for transmission networks
  #	based on contiguous linkage
  #
  if('close.and.contiguous.and.ancestry.cat' %in% relationship.types) {				
    df	<- df %>% categorise("close.and.contiguous.and.ancestry.cat", "not.close.or.noncontiguous",
                     list(contiguous = TRUE, ancestry=c("anc", "multiAnc"), categorical.distance = "close", label = "12"),
                     list(contiguous = TRUE, ancestry=c("desc", "multiDesc"), categorical.distance = "close", label = "21"),
                     list(contiguous = TRUE, ancestry=c("none", "complex"), categorical.distance = "close", label = "complex.or.no.ancestry"))
  }
  #
  #	group to determine probabilities for transmission networks
  #	based on adjacent linkage
  #
  if('close.and.adjacent.and.ancestry.cat' %in% relationship.types) {				
    df <- df %>% categorise("close.and.adjacent.and.ancestry.cat", "not.close.or.nonadjacent",
                     list(adjacent = TRUE, ancestry=c("anc", "multiAnc"), categorical.distance = "close", label = "12"),
                     list(adjacent = TRUE, ancestry=c("desc", "multiDesc"), categorical.distance = "close", label = "21"),
                     list(adjacent = TRUE, ancestry=c("none", "complex"), categorical.distance = "close", label = "complex.or.no.ancestry"))
  }		
  #
  #	group to transmission networks
  #	
  if('adjacent.and.proximity.cat' %in% relationship.types) {
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
#' @param df tibble with basic relationship types for paired individuals across windows. Must contain columns 'ID1','ID2','W_FROM','W_TO','TYPE_BASIC'. 
#' @param relationship.types names of relationship groups
#' @import tidyverse  
#' @return new tibble with columns host.1 host.2 GROUP TYPE K KEFF N NEFF.
#'  
get.keff.and.neff <- function(df, relationship.types, w.slide=NA){
  
  stopifnot(c('host.1', 'host.2', 'window.start', 'window.end', 'basic.classification') %in% colnames(df))
  
  df <- df %>% arrange(host.1, host.2, window.start)
  
  #	identify chunks of contiguous windows
  
  if(is.na(w.slide)) {
    # There are warnings here but it's OK unless _every_ pair in the dataset occurs in only one window
    
    w.slide <- df %>% group_by(host.1, host.2) %>% summarise(slide.width = min(diff(window.start))) %>% ungroup()
    
    w.slide	<- w.slide %>% filter(!is.na(slide.width))
    w.slide	<- ifelse(nrow(w.slide), min(w.slide$slide.width), 1L)		
    
    if(w.slide == Inf){
      stop("Cannot calculate sliding window size from this data")
    }
  }
  #	define chunks

  df <- df %>% 
    group_by(host.1, host.2) %>% 
    mutate(new.block = (function(x){
    if_else(is.na(lag(x,1)), TRUE, lag(x, 1) != x-w.slide) 
    })(window.start)) %>% 
    mutate(chunk.no = cumsum(new.block)) %>% 
    select(-new.block) %>% 
    ungroup()
  
  #	define chunk length in terms of non-overlapping windows	& number of windows in chunk
  
  df <- df %>% 
    group_by(host.1, host.2, chunk.no) %>% 
    mutate(effective.n.windows = (max(window.end) + 1 - min(window.start))/(window.end + 1 - window.start), n.windows = n()) %>% 
    ungroup()
  
  categorisation <- df %>% 
    select("basic.classification", relationship.types) %>% 
    distinct()

  #	for each chunk, count: windows by type and effective length of chunk
  #	then sum chunks
  
  category.parameters <- df %>% 
    group_by(host.1, host.2, chunk.no, basic.classification) %>%
    summarise(k = n(), k.eff = n()/n.windows[1] * effective.n.windows[1]) %>%
    ungroup()
  
  category.parameters <- category.parameters %>%
    group_by(host.1, host.2, basic.classification) %>%
    gather(statistic, value, k:k.eff) %>%
    ungroup()
  
  category.parameters <- category.parameters %>%
    group_by(host.1, host.2, statistic, basic.classification) %>%
    summarise(value = sum(value)) %>%
    ungroup()

  #	add relationship types
  
  category.parameters <- category.parameters %>% 
    inner_join(categorisation, by="basic.classification")
  
  #	melt relationship groups
  
  category.parameters <- category.parameters %>% 
    gather(categorisation, type, c("basic.classification", relationship.types)) %>%
    filter(!is.na(type))
  
  #	sum k and k.eff of same relationship state
  
  category.parameters <- category.parameters %>% 
    group_by(host.1, host.2, categorisation, type, statistic) %>%
    summarise(value = sum(value)) %>%
    ungroup()
  
  #	add zero-count relationship states (change to wide table and set NA's to zero's)
  
  category.parameters <- category.parameters %>% 
    complete(nesting(categorisation, type), nesting(host.1, host.2, statistic), fill = list(value=0))

  #	expand KEFF and K columns now that everything is done
  
  category.parameters <- category.parameters %>% spread(statistic, value)
  
  #	calculate N and NEFF
  
  category.parameters <- category.parameters %>%
    group_by(host.1, host.2, categorisation) %>%
    mutate(n = sum(k), n.eff = sum(k.eff)) %>%
    ungroup()
  
  category.parameters
}

#' @title Calculate marginal posterior probability for two individuals being in a particular relationship state
#' @description This function calculates the parameters that specify the marginal posterior probability for two individuals being in a particular relationship state. The marginal posterior is Beta distributed and this function calculates the ALPHA and BETA parameters.
#' @export  
#' @param df Input tibble   
#' @param n.type Calibration parameter for the prior: minimum number of windows of state to select a pair of individuals with confidence of at least at least confidence.cut, if the total number of windows is n.obs
#' @param n.obs Calibration parameter for the prior: total number of windows. 
#' @param confidence.cut Calibration parameter for the prior: confidence cut off.  
#' @return Input tibble with two additional columns POSTERIOR_ALPHA and POSTERIOR_BETA
#' @keywords internal
get.posterior.scores <- function(df){
  
  stopifnot(c('categorisation') %in% colnames(df))
  category.counts <- get.pairwise.relationship.category.counts()	
  category.counts <- category.counts %>%
    mutate(par.prior = n.type)
  
  df <- df %>% inner_join(category.counts)
  
  df <- df %>%
    mutate(posterior.alpha = par.prior/n.type + k.eff, 
           posterior.beta = par.prior*(1-(1/n.type)) + n.eff - k.eff, 
           posterior.score = (posterior.alpha - 1)/(posterior.alpha + posterior.beta - n.type))

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

get.pairwise.relationship.category.counts <- function() {
  tmp	<- matrix(c('proximity.3.way', 3,
                  'any.ancestry',2 ,
                  'close.x.contiguous', 4,					
                  'close.and.contiguous', 2,
                  'close.and.contiguous.and.directed', 2,															
                  'adjacent.and.proximity.cat', 3,
                  'close.and.adjacent.and.directed', 2,
                  'close.and.contiguous.and.ancestry.cat', 4,
                  'close.and.adjacent.and.ancestry.cat', 4,
                  'basic.classification', 24 ), ncol=2,byrow=TRUE)
  colnames(tmp)	<- c('categorisation','n.type')
  tmp <- tmp %>% as_tibble() %>% type_convert()
  tmp
}

