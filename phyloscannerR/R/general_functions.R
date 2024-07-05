# list.files seems to behave differently on different systems

#' @keywords internal
#' @export list.files.mod

list.files.mod <- function(path = ".", pattern = NULL, all.files = FALSE,
                           full.names = FALSE, recursive = FALSE,
                           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE){
  original <- list.files(path, pattern, all.files, full.names, recursive, ignore.case, include.dirs, no..)
  sapply(original, function(x) if(substr(x, 1, 2)=="./") {substr(x, 3, nchar(x))} else {x})
}

#' @keywords internal
#' @export get.suffix

get.suffix <- function(file.name, prefix, extension){
  
  file.name <- basename(file.name)
  prefix <- basename(prefix)

  substr(file.name, nchar(prefix)+1, nchar(file.name)-nchar(extension)-1)
}

#' @keywords internal
#' @export get.window.coords

get.window.coords <- function(string, regex = "^.*(?:.*\\D)?([0-9]+)_to_([0-9]+).*$"){

  start <- if(length(grep(regex, string))>0) as.numeric(sub(regex, "\\1", string)) else NA
  end   <- if(length(grep(regex, string))>0) as.numeric(sub(regex, "\\2", string)) else NA

  if(any(is.na(start)) | any(is.na(end))) {
    stop(paste0("ERROR: cannot determine window coordinates"))
  }
  
  list(start = start, end = end)
}

#' @keywords internal
#' @export get.window.coords.string

get.window.coords.string <- function(string, regex = ".*(?:.*\\D)?([0-9]+)_to_([0-9]+).*$"){
  wc <- get.window.coords(string, regex)
  
  paste0(wc$start, "_to_", wc$end)
}

#' @keywords internal
#' @export all.hosts.from.trees

all.hosts.from.trees <- function(ptrees){
  
  hosts <- ptrees %>% map("hosts.for.tips")
  
  hosts <- unique(unlist(hosts))
  hosts <- hosts[!is.na(hosts)]
  hosts <- hosts[order(hosts)]
  hosts
}

#' @export
#' @import tibble
#' @title Find bam and corresponding reference files 
#' @description This function finds bam and corresponding reference files in a given directory,
#' and groups them by a common sample ID as well as by an individual ID. 
#' @param data.dir Full path of data directory 
#' @param regex.person Regular expression with one set of round brackets, which identifies the person ID in the file name of bams and references
#' @param regex.bam Regular expression that identifies bam files, with one set of round brackets that identifies the sample ID.
#' @param regex.ref Regular expression that identifies ref files, with one set of round brackets that identifies the sample ID.
#' @return tibble with rows 'IND' (individual identifier), 'SAMPLE' (sample identifier), 'BAM' (bam file), and 'REF' (reference file).
find.bam.and.references<- function(data.dir, regex.person='^([A-Z0-9]+-[A-Z0-9]+)-.*$', regex.bam='^(.*)\\.bam$', regex.ref='^(.*)_ref\\.fasta$', verbose=1) 
{	
  ptyd		<- tibble(BAM=list.files(data.dir, pattern=regex.bam, full.names=TRUE, recursive=TRUE))
  ptyd		<- ptyd %>% 
    mutate( IND:= gsub(regex.person,'\\1',basename(BAM)) ) %>%
    mutate( SAMPLE:= gsub(regex.bam,'\\1',basename(BAM)) )				
  tmp			<- tibble(REF=list.files(data.dir, pattern=regex.ref, full.names=TRUE, recursive=TRUE))
  tmp			<- tmp %>% 
    mutate( IND:= gsub(regex.person,'\\1',basename(REF)) ) %>%
    mutate( SAMPLE:= gsub(regex.ref,'\\1',basename(REF)) )
  ptyd		<- full_join(ptyd, tmp, by=c('IND','SAMPLE'))
  if(verbose)	cat('\nFound data for individuals, n=', length(na.omit(unique(ptyd$IND))), ', found bam files, n=', length(na.omit(unique(ptyd$BAM))), ', found reference files, n=', length(na.omit(unique(ptyd$REF))))
  if(any(is.na(ptyd$BAM)))
  {
    warning('\nCould not find bam file for all individuals, n=', length(which(is.na(ptyd$BAM))),'\nSamples with missing bam files are ignored. Please check.')
    print( ptyd %>% filter( is.na(BAM) ) )
  }				
  if(any(is.na(ptyd$REF)))
  {
    warning('\nCould not find reference file for all individuals, n=', length(which(is.na(ptyd$REF))),'\nSamples with missing reference files are ignored. Please check.')
    print( ptyd %>% filter( is.na(REF) ) )
  }
  ptyd		<- ptyd %>% filter( !is.na(BAM) & !is.na(REF) )
  if(verbose)	cat('\nFound bam *and* reference files for individuals, n=', length(unique(ptyd$IND)) )
  ptyd	
}

#' @export
#' @import tibble
#' @author Oliver Ratmann
#' @title Group individuals for batched phyloscanner analysis
#' @description This function groups individuals for phyloscanner analyses, so that phylogenetic linkage between every pair of individuals is assessed at least once.
#'   Specifically, individuals are grouped into batches of specified size, and then, all possible pairs of batches are formed. 
#'   Each of these pairs of batches defines a group of individuals between whom phylogenetic linkages are assessed in one phyloscanner run.
#'   The number of individuals in each group is twice the batch size.
#' @param x Character vector of individual identifiers. 
#' @param batch.size Batch size. Default is 50.
#' @return tibble with rows 'IND' (individual identifiers), 'PTY_RUN' group for phyloscanner analysis, and 'BATCH' batch of individuals (not used further, but there should be two batches of individuals in each phyloscanner analysis).
#' @examples
#' x <- c("15-01402","15-04719","16-00616","16-00801","16-01173","16-01191","16-01302","16-01408","16-01414","16-01465","16-01581","16-01663","16-03259","16-03499","16-03516","16-03644","16-03807")
#' pty.runs	<- phyloscannerR::assign.groups.for.batched.phyloscanner.analysis(x, batch.size=50)
assign.groups.for.batched.phyloscanner.analysis<- function(x, batch.size=50)	
{
  dind	<- tibble(IND=x)
  dind	<- dind %>% mutate( BATCH:= 1L+floor((seq_len(nrow(dind))-1)/batch.size) )
  batches	<- dind %>% distinct(BATCH) %>% unlist() %>% unname()
  pty.runs<- t(combn(batches,2))
  pty.runs<- pty.runs %>% 
    as_tibble() %>% 
    rename(BATCH=V1, BATCH2=V2) %>%
    mutate( PTY_RUN:= seq_along(BATCH) )
  tmp		<- full_join(pty.runs, dind, by='BATCH')
  dind	<- dind %>% rename(IND2=IND, BATCH2=BATCH)
  tmp2	<- full_join(pty.runs, dind, by='BATCH2')
  dind	<- dind %>% rename(IND=IND2, BATCH=BATCH2)
  pty.runs<- tmp2 %>% 
    rename(IND=IND2) %>% 
    rbind(tmp) %>% 
    arrange(PTY_RUN) %>%
    select(PTY_RUN, IND, everything()) 
  pty.runs
}
