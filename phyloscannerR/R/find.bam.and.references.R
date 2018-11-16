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